#!/usr/bin/env python3
"""
SeqMetrics orchestrator. Runs only the requested modules, checking only
the dependencies THOSE modules need before doing any work -- nobody has
to install the whole panel to run one module.

Two input modes -- pick whichever matches how your data is organized;
neither is more "correct," they're just different groupings:

  Single-file mode (this project's own usage -- one combined FASTA per
  trim type, many sequences from many species pooled together):
    python run_features.py --modules basic,tail_hydrophobicity \
        --nt combined.fna --aa combined.faa --out-dir outputs/
    -> outputs/<module>.tsv, one file per module

  Batch mode (e.g. one FASTA per species, matching how the MLA-chapter
  project organized its inputs -- SeqMetrics doesn't care what the
  grouping IS, it just runs each module once per file in the directory):
    python run_features.py --modules basic,codon_usage \
        --batch species_fastas/ --out-dir outputs/
    -> outputs/<module>/<input_filename_stem>.tsv, one per (module, file)
    Nucleotide modules look for *.fna in the batch dir; protein modules
    look for *.faa. A module is silently skipped for a batch dir with no
    matching files, not an error -- lets one dir serve both input types.

Either way, how you got your data into a FASTA with sane headers (e.g.
this project's stage5_build_inputs.py, or a per-species DATA dict like
the old project used) is YOUR project's own adapter step, not SeqMetrics'
job -- kept out of this tool on purpose.

--modules is comma-separated. See README.md for the full module list and
docs/install_<module>.md for what each one actually needs installed.
"""
import argparse
import functools
import importlib
import platform
import re
import shlex
import shutil
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).parent

# (module_name, required_input: "nt"|"aa", real external binaries/libs to
# check for on PATH before running -- empty list = pure Python, always ready)
MODULES = {
    "basic":               {"input": "nt", "requires": []},
    "tail_hydrophobicity": {"input": "aa", "requires": []},
    "composition":         {"input": "aa", "requires": ["pepstats"]},
    "aggregation":         {"input": "aa", "requires": ["hcatk", "tango"]},
    "disorder":            {"input": "aa", "requires": []},  # checked via python-import, not PATH binary -- see PYMODULE_REQUIRES
    "codon_usage":         {"input": "nt", "requires": ["codonw"]},
    "coding_potential":    {"input": "nt", "requires": ["cpat"]},
    "tm_domain":           {"input": "aa", "requires": []},  # checked via venv binary, not conda -- see VENV_REQUIRES
    "localization":        {"input": "aa", "requires": ["pepstats", "perl", "java"]},
}

# Modules whose real dependency is an importable python module (not a PATH
# binary) -- resolved via resolve_pymodule() instead of resolve_binary().
PYMODULE_REQUIRES = {
    "disorder": "iupred3_lib",
}

# Modules whose real dependency lives in a plain python venv (not a conda
# env, and not on bare PATH) -- resolved via resolve_venv() instead of
# resolve_binary(). DeepTMHMM2's dependency stack (torch, fair-esm, etc.)
# was deliberately isolated in its own venv this session after an earlier
# `pip install` without an active env silently upgraded numpy/pandas/
# biopython/matplotlib in the SYSTEM python -- see DECISIONS-equivalent
# note in docs/install_tm_domain.md.
VENV_REQUIRES = {
    "tm_domain": "dtm2",
}

# The documented default conda env name per binary -- see docs/install_<module>.md.
# Used only as a fallback when the binary isn't already on bare PATH and no
# --module-env override was given.
DEFAULT_CONDA_ENVS = {
    "pepstats": "em_boss",
    "hcatk": "hca_tango",
    "tango": "hca_tango",
    "codonw": "codon_w",
    "cpat": "cpat",
    "iupred3_lib": "iupred3",
    "perl": "em_boss",
    "java": "em_boss",
}

# The documented default venv PATH per binary -- see docs/install_<module>.md.
# Relative to this script's own directory, so it works regardless of the
# caller's cwd (same reasoning as every other absolute-path fix this session).
DEFAULT_VENVS = {
    "dtm2": str(HERE / ".venvs" / "deeptmhmm2"),
}

# iupred3_lib isn't pip-installable -- it's two files a user obtains directly
# from the authors (see docs/install_disorder.md) and drops somewhere on
# PYTHONPATH. This is the documented default location inside the WSL conda
# env this lab already uses; `~` is expanded by the WSL shell itself, not by
# this script, so it's correct for any username. Only meaningful for the
# "iupred3" env specifically -- there is no equivalent for a bare-PATH or
# native-conda resolution, since a normal `pip install` env wouldn't need it.
PYMODULE_PYTHONPATH_WSL = {
    "iupred3_lib": "~/miniconda3/envs/iupred3/apps/iupred3",
}


# Sourced before every WSL-bridged conda call -- a plain `wsl.exe conda ...`
# does NOT go through a login shell that sources conda's init script (the
# same gotcha hit earlier this session doing WSL environment checks by hand),
# so both the availability check and the real run need this, not just one.
_WSL_CONDA_INIT = (
    "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
    "source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null"
)


@functools.lru_cache(maxsize=None)
def _conda_available(via_wsl):
    """Cheap pre-check so we don't pay conda's own startup overhead (a real
    cost -- confirmed ~10s per WSL cold start) on every binary lookup on a
    system that simply doesn't have conda at all. Cached once per
    (native/WSL) -- always call this positionally, never by keyword, or
    lru_cache treats it as a different, uncached call (confirmed: this was
    the actual reason repeat calls weren't getting cached during testing)."""
    if via_wsl:
        if shutil.which("wsl.exe") is None:
            return False
        try:
            result = subprocess.run(["wsl.exe", "bash", "-lc", f"{_WSL_CONDA_INIT}; command -v conda"],
                                     capture_output=True, text=True, timeout=15)
            return result.returncode == 0 and bool(result.stdout.strip())
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False
    return shutil.which("conda") is not None


def _conda_has_binary(env_name, binary, via_wsl):
    """Check only -- does NOT return something reusable for the real run,
    since the WSL case needs the sourcing string rebuilt with the real args
    anyway (see build_command). Always calls _conda_available positionally."""
    if not _conda_available(via_wsl):
        return False
    if via_wsl:
        cmd = ["wsl.exe", "bash", "-lc", f"{_WSL_CONDA_INIT}; conda run -n {env_name} which {binary}"]
    else:
        cmd = ["conda", "run", "-n", env_name, "which", binary]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        return result.returncode == 0 and bool(result.stdout.strip())
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def resolve_binary(binary, overrides_by_module, module_name):
    """Returns (resolution, found). resolution is (kind, env_name) where
    kind is 'direct' | 'conda' | 'wsl_conda' -- pass to build_command() to
    get the actual argv for the real run, once it's time to call the tool
    for real rather than just checking it exists.

    Resolution order: --module-env override (conda env name, or a directory
    path for a non-conda install) -> bare PATH -> the documented default
    conda env, tried natively first, then (Windows only) bridged via WSL.
    """
    override = overrides_by_module.get(module_name)
    if override:
        looks_like_path = ("/" in override) or ("\\" in override)
        if not looks_like_path and _conda_has_binary(override, binary, False):
            return ("conda", override), True
        if not looks_like_path and platform.system() == "Windows" and _conda_has_binary(override, binary, True):
            return ("wsl_conda", override), True
        if not looks_like_path:
            print(f"[run_features] --module-env {module_name}={override} didn't resolve "
                  f"'{binary}' as a conda env; trying it as a PATH directory instead.",
                  file=sys.stderr)
        if shutil.which(binary, path=override):
            return ("direct", None), True
        return ("direct", None), False

    if shutil.which(binary) is not None:
        return ("direct", None), True

    default_env = DEFAULT_CONDA_ENVS.get(binary)
    if default_env:
        if _conda_has_binary(default_env, binary, False):
            return ("conda", default_env), True
        if platform.system() == "Windows" and _conda_has_binary(default_env, binary, True):
            return ("wsl_conda", default_env), True

    return ("direct", None), False


_WIN_PATH_RE = re.compile(r"^[A-Za-z]:[\\/]")


def _win_to_wsl_path(s):
    """C:\\Users\\x\\f.faa -> /mnt/c/Users/x/f.faa. Needed because a WSL-bridged
    call still receives Windows-style paths from the rest of this script (it
    runs on the Windows host) -- WSL's bash can reach them under /mnt/<drive>,
    but not under the raw Windows spelling. Only ever applied on Windows, and
    only to args that actually look like a Windows absolute path (leaves
    plain flags like "-f" or "--jobs" untouched)."""
    if not _WIN_PATH_RE.match(s):
        return s
    drive = s[0].lower()
    rest = s[2:].replace("\\", "/")
    return f"/mnt/{drive}{rest}"


def build_command(resolution, binary, args):
    """Turn a resolve_binary() result into the actual argv to subprocess.run
    for the real call, now that it's time to run the tool rather than just
    check it exists."""
    kind, env_name = resolution
    if kind == "direct":
        return [binary, *args]
    if kind == "conda":
        return ["conda", "run", "-n", env_name, binary, *args]
    if kind == "wsl_conda":
        wsl_binary = _win_to_wsl_path(binary)
        wsl_args = [_win_to_wsl_path(str(a)) for a in args]
        arg_str = " ".join(shlex.quote(a) for a in [wsl_binary, *wsl_args])
        return ["wsl.exe", "bash", "-lc", f"{_WSL_CONDA_INIT}; conda run -n {env_name} {arg_str}"]
    if kind == "venv":
        return [str(_venv_binary_path(env_name, binary)), *args]
    raise ValueError(f"unknown resolution kind: {kind}")


def _venv_binary_path(venv_dir, binary):
    """A venv's own binary, called directly by full path -- no `activate`
    needed for a single subprocess call, unlike conda run. Windows venvs put
    executables in Scripts\\ with a .exe suffix; everywhere else it's bin/
    with no suffix."""
    venv_dir = Path(venv_dir)
    if platform.system() == "Windows":
        return venv_dir / "Scripts" / f"{binary}.exe"
    return venv_dir / "bin" / binary


def resolve_venv(binary, overrides_by_module, module_name):
    """Same shape as resolve_binary()/resolve_pymodule(): --module-env
    override (a venv directory path this time, not a conda env name) ->
    the documented default venv location. No bare-PATH or conda fallback --
    a venv-isolated tool was deliberately kept out of any shared env."""
    override = overrides_by_module.get(module_name)
    if override:
        p = _venv_binary_path(override, binary)
        return ("venv", override), p.exists()

    default_venv = DEFAULT_VENVS.get(binary)
    if default_venv:
        p = _venv_binary_path(default_venv, binary)
        if p.exists():
            return ("venv", default_venv), True

    return ("venv", None), False


def _conda_has_pymodule(env_name, module, via_wsl):
    """Same idea as _conda_has_binary, but checks `python -c "import X"`
    instead of `which X` -- IUPred3 is two files a user drops on PYTHONPATH,
    not a binary. For the WSL case, if `module` needs a documented PYTHONPATH
    addition (see PYMODULE_PYTHONPATH_WSL), prepend it as an env-var-prefixed
    shell invocation so `conda run` inherits it."""
    if not _conda_available(via_wsl):
        return False
    pythonpath = PYMODULE_PYTHONPATH_WSL.get(module)
    if via_wsl:
        prefix = f"PYTHONPATH={pythonpath} " if pythonpath else ""
        cmd = ["wsl.exe", "bash", "-lc",
               f'{_WSL_CONDA_INIT}; {prefix}conda run -n {env_name} python -c "import {module}"']
    else:
        cmd = ["conda", "run", "-n", env_name, "python", "-c", f"import {module}"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def resolve_pymodule(module, overrides_by_module, module_name):
    """Same shape/order as resolve_binary(), for an importable python module
    instead of a PATH binary: --module-env override (conda env name only --
    a path override for a pymodule isn't supported yet, unlike resolve_binary)
    -> already importable in this interpreter -> the documented default conda
    env, native then (Windows only) WSL-bridged."""
    override = overrides_by_module.get(module_name)
    if override and "/" not in override and "\\" not in override:
        if _conda_has_pymodule(override, module, False):
            return ("conda", override), True
        if platform.system() == "Windows" and _conda_has_pymodule(override, module, True):
            return ("wsl_conda", override), True

    try:
        importlib.import_module(module)
        return ("direct", None), True
    except ImportError:
        pass

    default_env = DEFAULT_CONDA_ENVS.get(module)
    if default_env:
        if _conda_has_pymodule(default_env, module, False):
            return ("conda", default_env), True
        if platform.system() == "Windows" and _conda_has_pymodule(default_env, module, True):
            return ("wsl_conda", default_env), True

    return ("direct", None), False


def check_module_ready(name, spec, overrides_by_module):
    """Returns (missing: list[str], resolved: dict[binary, (kind, env_name)]) --
    pass an entry from `resolved` to build_command() to get the real argv."""
    missing = []
    resolved = {}
    for b in spec["requires"]:
        resolution, found = resolve_binary(b, overrides_by_module, name)
        if found:
            resolved[b] = resolution
        else:
            missing.append(b)
    pymodule = PYMODULE_REQUIRES.get(name)
    if pymodule:
        resolution, found = resolve_pymodule(pymodule, overrides_by_module, name)
        if found:
            resolved[pymodule] = resolution
        else:
            missing.append(f"{pymodule} (see docs/install_{name}.md)")
    venv_binary = VENV_REQUIRES.get(name)
    if venv_binary:
        resolution, found = resolve_venv(venv_binary, overrides_by_module, name)
        if found:
            resolved[venv_binary] = resolution
        else:
            missing.append(f"{venv_binary} (see docs/install_{name}.md)")
    return missing, resolved


# Non-interactive version flag per binary, verified by direct test against
# each real tool before adding it here -- not guessed. codonw is
# deliberately absent: it has no non-interactive version flag, and probing
# it with an unrecognized flag drops into its interactive menu
# ("Press return or enter to continue"), which would hang a real run
# waiting on stdin. hcatk and tango are also absent -- confirmed via a
# timeout-guarded test that neither hangs, but neither prints anything
# resembling a version string either (hcatk emits unrelated deprecation
# warnings; tango just reports a normal usage error). Anything not in this
# table returns "unknown" rather than a guessed or hung probe.
VERSION_PROBES = {
    "pepstats": ["-version"],
    "perl": ["-v"],
    "java": ["-version"],
    "cpat": ["--version"],
    "hmmsearch": ["-h"],
}


def get_tool_version(binary, resolution):
    """Best-effort version string for a resolved binary, for the audit log.
    15s timeout and a broad except -- a version probe must never be the
    thing that makes a real run fail or hang."""
    if binary not in VERSION_PROBES:
        return "unknown"
    try:
        cmd = build_command(resolution, binary, VERSION_PROBES[binary])
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
        combined = (result.stdout + result.stderr).strip()
        if not combined:
            return "unknown"
        if binary == "hmmsearch":
            # -h prints a multi-line banner; the version is on its own line,
            # not necessarily the first one.
            for line in combined.splitlines():
                if "HMMER" in line:
                    return line.strip("# ").strip()
            return "unknown"
        return combined.splitlines()[0].strip()
    except Exception:
        return "unknown"


def get_localizer_commit(ref_dir):
    """LOCALIZER has no version flag at all (it's a plain script, not a
    versioned release on PyPI/conda) -- the git commit hash of the cloned
    install is the actual reproducibility-relevant fact here, not a
    made-up version string. ref_dir is the path to LOCALIZER.py itself
    (see docs/install_localization.md); its repo root is two levels up
    (LOCALIZER/Scripts/LOCALIZER.py -> LOCALIZER/)."""
    repo_root = str(Path(ref_dir).parent.parent)
    try:
        if ref_dir.startswith("/"):
            cmd = ["wsl.exe", "bash", "-lc",
                   f"git -C {shlex.quote(repo_root)} rev-parse --short HEAD"]
        else:
            cmd = ["git", "-C", repo_root, "rev-parse", "--short", "HEAD"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        return result.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


class RunLog:
    """One append-only, plain-text log per run_features.py invocation --
    matches this repo's own wrapper scripts' logging style (timestamped
    lines to a .log file) rather than introducing a second, structured
    format. Thread-safe: --jobs > 1 means multiple (module, file) tasks
    write to this concurrently."""

    def __init__(self, out_dir):
        log_dir = out_dir / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.path = log_dir / f"run_{timestamp}.log"
        self._lock = threading.Lock()
        self._fh = open(self.path, "w", newline="\n")

    def write(self, msg):
        line = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n"
        with self._lock:
            self._fh.write(line)
            self._fh.flush()

    def close(self):
        self._fh.close()


def run_module(name, input_path, out_path, resolved, ref_dir=None):
    """Run one module on one input file, writing one output table. This is
    the single unit both single-file mode and batch mode call -- batch mode
    is just this function in a loop, not a separate code path. `resolved`
    is the {binary: invocation_prefix} dict from check_module_ready, for
    modules that need to wrap their real subprocess call in `conda run -n
    ENV ...` (or `wsl.exe conda run -n ENV ...`) instead of calling the
    binary bare -- not needed by basic/tail_hydrophobicity, which are pure
    Python with nothing external to resolve. `ref_dir` is only meaningful
    for modules that score against a pre-built per-species reference
    (codon_usage's cai.coa/fop.coa/cbi.coa) rather than running the tool
    cold -- see --module-ref."""
    if name == "basic":
        subprocess.run([sys.executable, str(HERE / "modules/basic/compute_basic_features.py"),
                         str(input_path), "--out", str(out_path)], check=True)
    elif name == "tail_hydrophobicity":
        subprocess.run([sys.executable, str(HERE / "modules/tail_hydrophobicity/compute_ctth.py"),
                         str(input_path), "--out", str(out_path)], check=True)

    elif name == "disorder":
        # batch_iupred_features_cysexcl.py's single-file mode is -i FILE -o
        # FILE, which is exactly this function's per-file contract already --
        # no scratch-dir relocation needed, unlike composition/aggregation.
        script = HERE / "modules/disorder/batch_iupred_features_cysexcl.py"
        kind, env_name = resolved.get("iupred3_lib", ("direct", None))
        if kind == "direct":
            cmd = [sys.executable, str(script), "-i", str(input_path), "-o", str(out_path), "-f"]
        elif kind == "conda":
            cmd = ["conda", "run", "-n", env_name, "python", str(script),
                   "-i", str(input_path), "-o", str(out_path), "-f"]
        elif kind == "wsl_conda":
            pythonpath = PYMODULE_PYTHONPATH_WSL.get("iupred3_lib", "")
            wsl_args = [_win_to_wsl_path(str(a)) for a in (script, "-i", input_path, "-o", out_path)]
            arg_str = " ".join(shlex.quote(a) for a in [wsl_args[0], wsl_args[1], wsl_args[2], wsl_args[3], wsl_args[4], "-f"])
            cmd = ["wsl.exe", "bash", "-lc",
                   f"{_WSL_CONDA_INIT}; PYTHONPATH={pythonpath} conda run -n {env_name} python {arg_str}"]
        else:
            raise ValueError(f"unknown resolution kind: {kind}")
        subprocess.run(cmd, check=True)

    elif name == "composition":
        # pepstats-flow writes to <outdir>/tables/<input_stem>.tsv, not a path
        # we control directly -- run it into a scratch dir, then relocate the
        # one file we asked for and discard the rest (pepstats_raw/, master/,
        # logs/) it always creates alongside it.
        script = HERE / "modules/composition/pepstats-flow"
        scratch_dir = out_path.parent / f"_composition_scratch_{out_path.stem}"
        scratch_dir.mkdir(parents=True, exist_ok=True)
        cmd = build_command(resolved["pepstats"], "bash",
                             [str(script), "-f", str(input_path), "-o", str(scratch_dir), "-t"])
        subprocess.run(cmd, check=True)
        produced = scratch_dir / "tables" / f"{input_path.stem}.tsv"
        if not produced.exists():
            raise FileNotFoundError(f"composition produced no output at {produced}")
        shutil.move(str(produced), str(out_path))
        shutil.rmtree(scratch_dir, ignore_errors=True)

    elif name == "aggregation":
        # run_hca_tango.py writes <input_stem>.HCA_TANGO.tsv under --work-dir
        # (archiving any pre-existing same-name output with a timestamp suffix
        # rather than overwriting) -- same relocate-from-scratch pattern as
        # composition, since we don't control its output filename directly.
        script = HERE / "modules/aggregation/run_hca_tango.py"
        resolution = resolved.get("hcatk") or resolved.get("tango") or ("direct", None)
        scratch_dir = out_path.parent / f"_aggregation_scratch_{out_path.stem}"
        scratch_dir.mkdir(parents=True, exist_ok=True)
        kind, _ = resolution
        if kind == "direct":
            cmd = [sys.executable, str(script), str(input_path), "--work-dir", str(scratch_dir), "--jobs", "4"]
        else:
            cmd = build_command(resolution, "python",
                                 [str(script), str(input_path), "--work-dir", str(scratch_dir), "--jobs", "4"])
        subprocess.run(cmd, check=True)
        produced = scratch_dir / f"{input_path.stem}.HCA_TANGO.tsv"
        if not produced.exists():
            candidates = sorted(scratch_dir.glob(f"{input_path.stem}.HCA_TANGO*.tsv"),
                                 key=lambda p: p.stat().st_mtime)
            if not candidates:
                raise FileNotFoundError(f"aggregation produced no output in {scratch_dir}")
            produced = candidates[-1]
        shutil.move(str(produced), str(out_path))
        shutil.rmtree(scratch_dir, ignore_errors=True)

    elif name == "codon_usage":
        # Unlike every other module, codonW cannot score against nothing --
        # obtain_indices.sh needs cai.coa/fop.coa/cbi.coa already sitting in
        # its working directory (codonW reads them as bare relative
        # filenames, not by path) before it will run at all. Building those
        # references is a separate, explicit upstream step (see
        # docs/install_codon_usage.md) -- this is scoring only, and refuses
        # to run rather than let codonW silently fall back to its own
        # generic built-in tables, which were judged wrong for cross-species
        # comparison. No fallback here mirrors obtain_indices.sh's own
        # "Zero-Default Failsafe" -- fail at the SeqMetrics level with a
        # clear message instead of two layers deep in a subprocess.
        if ref_dir is None:
            raise ValueError(
                "codon_usage needs a pre-built reference: pass "
                "--module-ref codon_usage=<dir containing cai.coa> "
                "(see docs/install_codon_usage.md -- Stage 1)")
        ref_dir = Path(ref_dir)
        cai = ref_dir / "cai.coa"
        if not cai.exists():
            raise FileNotFoundError(
                f"codon_usage: no cai.coa at {ref_dir} -- Stage 1 "
                f"(calculate_indices.sh) hasn't been run for this species yet")

        # species is derived from the input file's own name, same convention
        # orchestrate_codonw.sh itself uses -- this assumes the input file IS
        # already one species' sequences. For this project's own pooled,
        # multi-species composite FASTAs, that assumption does NOT hold yet
        # (a real, separately-tracked gap -- see the species-splitting
        # adapter discussion) and this module should not be pointed at those
        # files until that splitting step exists.
        species = input_path.stem
        script = HERE / "modules/codon_usage/obtain_indices.sh"
        scratch_dir = out_path.parent / f"_codon_usage_scratch_{out_path.stem}"
        scratch_dir.mkdir(parents=True, exist_ok=True)
        for coa_name in ("cai.coa", "fop.coa", "cbi.coa"):
            src = ref_dir / coa_name
            if src.exists():
                shutil.copy(str(src), str(scratch_dir / coa_name))
            else:
                (scratch_dir / coa_name).touch()  # matches calculate_indices.sh's own placeholder convention

        cmd = build_command(resolved["codonw"], "bash",
                             [str(script), str(input_path), species, str(scratch_dir)])
        subprocess.run(cmd, check=True)
        produced = scratch_dir / f"{species}_indices.out"
        if not produced.exists():
            raise FileNotFoundError(f"codon_usage produced no output at {produced}")
        shutil.move(str(produced), str(out_path))
        shutil.rmtree(scratch_dir, ignore_errors=True)

    elif name == "coding_potential":
        # Same shape as codon_usage: scoring only, against a pre-built
        # per-species reference (a hexamer table + logit model) that Stage 1
        # (make_hexamer_tab -> make_logitModel, upstream, not SeqMetrics' job
        # -- see docs/install_coding_potential.md) must already exist. No
        # fallback here either, for the same reason as codon_usage: CPAT
        # cannot meaningfully score without one.
        if ref_dir is None:
            raise ValueError(
                "coding_potential needs a pre-built reference: pass "
                "--module-ref coding_potential=<dir containing hexamer.tsv "
                "and logit.RData> (see docs/install_coding_potential.md)")
        ref_dir = Path(ref_dir)
        hexamer = ref_dir / "hexamer.tsv"
        logit_model = ref_dir / "logit.RData"
        if not hexamer.exists() or not logit_model.exists():
            raise FileNotFoundError(
                f"coding_potential: need both {hexamer} and {logit_model} -- "
                f"Stage 1 hasn't been run for this species yet")

        scratch_dir = out_path.parent / f"_coding_potential_scratch_{out_path.stem}"
        scratch_dir.mkdir(parents=True, exist_ok=True)
        prefix = scratch_dir / "cpat"
        cmd = build_command(resolved["cpat"], "cpat",
                             ["-g", str(input_path), "-d", str(logit_model),
                              "-x", str(hexamer), "-o", str(prefix)])
        # cwd=scratch_dir: confirmed by a real leftover file (CPAT_run_info.log
        # found sitting at this repo's own root, not in any scratch dir) that
        # CPAT writes at least one file to the process's actual cwd,
        # independent of the -o prefix that scopes everything else. Without
        # this, every real run pollutes whatever directory the user happened
        # to invoke run_features.py from. For the wsl_conda case, Python's
        # own subprocess cwd doesn't reach into the WSL bash session wsl.exe
        # starts, so cd into the (WSL-translated) scratch dir inside the
        # bash -lc string itself instead.
        kind = resolved["cpat"][0]
        if kind == "wsl_conda":
            wsl_scratch = _win_to_wsl_path(str(scratch_dir))
            cmd = [cmd[0], cmd[1], "-lc", f"cd {shlex.quote(wsl_scratch)} && {cmd[3]}"]
            subprocess.run(cmd, check=True)
        else:
            subprocess.run(cmd, check=True, cwd=str(scratch_dir))
        raw = Path(f"{prefix}.ORF_prob.best.tsv")
        if not raw.exists():
            raise FileNotFoundError(f"coding_potential produced no output at {raw}")

        # CPAT uppercases every ID it writes and truncates pipe-delimited
        # headers oddly (see repair_cpat_ids.py's module docstring for the
        # confirmed real bug reports this fixes) -- always restore original
        # case, using the species-specific extraction rule for the three
        # species that need it, "generic" (first whitespace token, verbatim)
        # for anything else, including this project's own composite headers.
        repair_script = HERE / "modules/coding_potential/repair_cpat_ids.py"
        subprocess.run([sys.executable, str(repair_script),
                         "--fasta", str(input_path), "--cpat-table", str(raw),
                         "--out", str(out_path), "--species", "generic"], check=True)
        shutil.rmtree(scratch_dir, ignore_errors=True)

    elif name == "tm_domain":
        # run_dtm2.py is fully self-contained (runs dtm2 into its own temp
        # dir, parses the four output files into one flat TSV, cleans up) --
        # unlike composition/aggregation, no scratch-dir relocation needed
        # here at this level.
        script = HERE / "modules/tm_domain/run_dtm2.py"
        kind, env_name = resolved.get("dtm2", ("venv", None))
        if kind != "venv" or env_name is None:
            raise ValueError(f"tm_domain: unexpected resolution {resolved.get('dtm2')}")
        dtm2_path = _venv_binary_path(env_name, "dtm2")
        subprocess.run([sys.executable, str(script), "-i", str(input_path),
                         "-o", str(out_path), "--dtm2", str(dtm2_path)], check=True)

    elif name == "localization":
        # ref_dir here is not a directory of pre-built references like
        # codon_usage/coding_potential -- it's the direct path to the
        # user's own LOCALIZER.py (never bundled, see
        # docs/install_localization.md). run_localizer.py itself already
        # does the run-then-parse two-step and writes
        # <input_stem>.localizer_binary.tsv into --output; relocate that
        # one file the same way composition/aggregation do.
        if ref_dir is None:
            raise ValueError(
                "localization needs --module-ref localization=<path to LOCALIZER.py> "
                "(see docs/install_localization.md)")
        script = HERE / "modules/localization/run_localizer.py"
        parser_script = HERE / "modules/localization/localizer_bin.py"
        scratch_dir = out_path.parent / f"_localization_scratch_{out_path.stem}"
        scratch_dir.mkdir(parents=True, exist_ok=True)
        cmd = build_command(resolved["pepstats"], "python",
                             [str(script), str(input_path),
                              "--localizer-script", str(ref_dir),
                              "--mode", "plant",
                              "--output", str(scratch_dir),
                              "--localizer-parser", str(parser_script),
                              "--quiet"])
        subprocess.run(cmd, check=True)
        produced = scratch_dir / f"{input_path.stem}.localizer_binary.tsv"
        if not produced.exists():
            raise FileNotFoundError(f"localization produced no output at {produced}")
        shutil.move(str(produced), str(out_path))
        shutil.rmtree(scratch_dir, ignore_errors=True)

    else:
        print(f"[run_features] '{name}' needs its own real invocation wired in "
              f"(see modules/{name}/ and docs/install_{name}.md) -- not yet "
              f"connected to this orchestrator, unlike basic/tail_hydrophobicity. "
              f"Resolved binaries for when it is: {resolved}", file=sys.stderr)
        return False
    return True


EXT_FOR_INPUT = {"nt": ".fna", "aa": ".faa"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--modules", required=True, help="Comma-separated module names")
    ap.add_argument("--nt", type=Path, help="Single-file mode: nucleotide FASTA")
    ap.add_argument("--aa", type=Path, help="Single-file mode: protein FASTA")
    ap.add_argument("--batch", type=Path, help="Batch mode: directory of FASTA files "
                     "(*.fna for nucleotide modules, *.faa for protein modules) -- "
                     "runs each requested module once per matching file")
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--module-env", action="append", default=[], metavar="MODULE=VALUE",
                     help="Override environment resolution for one module. VALUE is a "
                          "conda env name (tried via `conda run -n VALUE`), or a directory "
                          "path (added to PATH for a direct call, for non-conda installs). "
                          "Repeatable, e.g. --module-env composition=em_boss "
                          "--module-env codon_usage=/opt/codonw/bin")
    ap.add_argument("--module-ref", action="append", default=[], metavar="MODULE=PATH",
                     help="An external asset path a module needs beyond its resolved "
                          "binaries. For codon_usage/coding_potential: a pre-built "
                          "per-species reference DIRECTORY (cai.coa/fop.coa/cbi.coa, or "
                          "hexamer.tsv+logit.RData -- see their install docs). For "
                          "codon_usage specifically, this also assumes the input file IS "
                          "that one species' sequences, species taken from the filename. "
                          "For localization: the direct FILE path to LOCALIZER.py itself "
                          "(never bundled -- see docs/install_localization.md).")
    ap.add_argument("--jobs", type=int, default=1, metavar="N",
                     help="Run up to N (module, file) tasks concurrently -- e.g. "
                          "batch mode across many species files for one module, or "
                          "several modules on the same single-file input. Each task "
                          "already writes to its own scratch dir and output path, so "
                          "this is safe to raise. Default 1 (fully sequential, today's "
                          "existing behavior). WSL-bridged modules each spawn a "
                          "wsl.exe/conda session per task -- don't set this arbitrarily "
                          "high on a machine with limited WSL headroom.")
    args = ap.parse_args()

    if args.batch and (args.nt or args.aa):
        sys.exit("Use either --batch, or --nt/--aa, not both")
    if not args.batch and not args.nt and not args.aa:
        sys.exit("Need --batch DIR, or at least one of --nt/--aa")

    # Resolve to absolute before anything reaches run_module -- a relative
    # path is only valid relative to THIS process's cwd, but a WSL-bridged
    # module runs in WSL's own bash with a different cwd entirely (confirmed:
    # a relative --aa path silently produced "file does not exist" inside
    # WSL, since _win_to_wsl_path only recognizes an absolute Windows path
    # to translate -- a bare relative path passed through untouched and
    # meant nothing on the other side of the bridge).
    args.out_dir = args.out_dir.resolve()
    if args.nt:
        args.nt = args.nt.resolve()
    if args.aa:
        args.aa = args.aa.resolve()
    if args.batch:
        args.batch = args.batch.resolve()

    overrides_by_module = {}
    for entry in args.module_env:
        if "=" not in entry:
            sys.exit(f"--module-env expects MODULE=VALUE, got: {entry!r}")
        mod, val = entry.split("=", 1)
        overrides_by_module[mod] = val

    refs_by_module = {}
    for entry in args.module_ref:
        if "=" not in entry:
            sys.exit(f"--module-ref expects MODULE=DIR, got: {entry!r}")
        mod, val = entry.split("=", 1)
        # Resolve to absolute now, same reason as --nt/--aa/--batch/--out-dir
        # above -- confirmed via a real test that CPAT itself has no cwd
        # requirement (it embeds whatever path string it's given verbatim,
        # and works fine from an unrelated cwd when given absolute paths),
        # but a relative --module-ref value would still silently depend on
        # whatever cwd this script happened to be invoked from, which is
        # exactly the class of bug already found and fixed once this
        # session for the other path arguments.
        #
        # EXCEPT: a value that's already a POSIX-absolute path (starts with
        # "/") is left untouched. This script runs natively on Windows, so
        # Path(val).resolve() interprets a WSL-native path like
        # "/home/user/..." as relative to Windows' own filesystem root --
        # confirmed by a real failure, --module-ref localization=<a real
        # WSL path to LOCALIZER.py, which only exists inside WSL, not
        # mirrored on the Windows side> came out mangled into something
        # under Git's own install directory. _win_to_wsl_path() already
        # knows to leave a non-Windows-looking path alone; this needs the
        # same guard one step earlier, before Windows' own Path.resolve()
        # gets a chance to mangle it first.
        if val.startswith("/"):
            refs_by_module[mod] = val
        else:
            refs_by_module[mod] = str(Path(val).resolve())

    requested = args.modules.split(",")
    unknown = [m for m in requested if m not in MODULES]
    if unknown:
        sys.exit(f"Unknown module(s): {unknown}. Known: {sorted(MODULES)}")
    unknown_overrides = [m for m in overrides_by_module if m not in MODULES]
    if unknown_overrides:
        sys.exit(f"--module-env given for unknown module(s): {unknown_overrides}")
    unknown_refs = [m for m in refs_by_module if m not in MODULES]
    if unknown_refs:
        sys.exit(f"--module-ref given for unknown module(s): {unknown_refs}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    log = RunLog(args.out_dir)
    log.write(f"invoked: {' '.join(sys.argv)}")
    print(f"[run_features] logging to {log.path}", file=sys.stderr)

    # Flat list of independent (module, file) tasks across every requested,
    # ready module -- built up front so --jobs can run all of them through
    # one worker pool, not one pool per module. Each tuple is exactly
    # run_module()'s own argument shape plus a label for reporting.
    tasks = []

    for name in requested:
        spec = MODULES[name]

        if not args.batch:
            if spec["input"] == "nt" and not args.nt:
                sys.exit(f"Module '{name}' needs --nt (nucleotide FASTA)")
            if spec["input"] == "aa" and not args.aa:
                sys.exit(f"Module '{name}' needs --aa (protein FASTA)")

        missing, resolved = check_module_ready(name, spec, overrides_by_module)
        if missing:
            msg = (f"SKIPPING '{name}' -- missing: {missing}. "
                   f"See docs/install_{name}.md, or pass --module-env {name}=<env_or_path>.")
            print(f"[run_features] {msg}", file=sys.stderr)
            log.write(msg)
            continue

        if resolved:
            summary = {b: (kind if env is None else f"{kind}:{env}") for b, (kind, env) in resolved.items()}
            print(f"[run_features] '{name}' resolved: {summary}", file=sys.stderr)
            versions = {b: get_tool_version(b, res) for b, res in resolved.items()}
            log.write(f"'{name}' resolved: {summary} -- versions: {versions}")
        ref = refs_by_module.get(name)
        if ref:
            log.write(f"'{name}' --module-ref: {ref}")
            if name == "localization":
                log.write(f"'{name}' LOCALIZER commit: {get_localizer_commit(ref)}")

        if args.batch:
            ext = EXT_FOR_INPUT[spec["input"]]
            matches = sorted(args.batch.glob(f"*{ext}"))
            if not matches:
                msg = f"SKIPPING '{name}' -- no {ext} files in {args.batch}"
                print(f"[run_features] {msg}", file=sys.stderr)
                log.write(msg)
                continue
            module_out_dir = args.out_dir / name
            module_out_dir.mkdir(parents=True, exist_ok=True)
            for input_path in matches:
                out_path = module_out_dir / f"{input_path.stem}.tsv"
                tasks.append((name, input_path, out_path, resolved, ref))
        else:
            input_path = args.nt if spec["input"] == "nt" else args.aa
            out_path = args.out_dir / f"{name}.tsv"
            tasks.append((name, input_path, out_path, resolved, ref))

    def _run_one(task):
        # run_module() raises (CalledProcessError, FileNotFoundError,
        # ValueError, ...) for most real failures rather than returning
        # False -- fine when everything ran sequentially and one failure
        # was meant to stop the whole invocation, but with --jobs > 1 an
        # uncaught exception from one (module, file) task would blow up
        # the whole batch and lose visibility into every other task still
        # in flight, exactly when --jobs matters most (many species at
        # once). Catch here and report as a normal failed result instead.
        name, input_path, out_path, resolved, ref = task
        try:
            ok = run_module(name, input_path, out_path, resolved, ref)
        except Exception as e:
            return name, input_path, out_path, False, str(e)
        return name, input_path, out_path, ok, None

    if args.jobs <= 1:
        results = (_run_one(t) for t in tasks)
    else:
        pool = ThreadPoolExecutor(max_workers=args.jobs)
        futures = [pool.submit(_run_one, t) for t in tasks]
        results = (f.result() for f in as_completed(futures))

    any_failed = False
    for name, input_path, out_path, ok, err in results:
        if ok:
            msg = f"'{name}' on {input_path.name} -> {out_path}"
        else:
            any_failed = True
            msg = f"'{name}' on {input_path.name} FAILED: {err}" if err else \
                  f"'{name}' on {input_path.name} FAILED -- see stderr above"
        print(f"[run_features] {msg}", file=sys.stderr if not ok else sys.stdout)
        log.write(msg)

    log.close()
    if any_failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
