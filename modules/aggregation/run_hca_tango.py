#!/usr/bin/env python3
"""
HCA + TANGO production pipeline (auditable, parallel, timed)

Author:       Adekola Owoyemi
Maintained by: Casola Lab
Version:      0.1.0 (pre-release, unshared)

USAGE:
    ./run_hca_tango.py input.faa

OPTIONS:
    --jobs N        Parallel TANGO workers (default: 4)
    --work-dir DIR  Directory for all outputs (default: caller's cwd).
                     Pre-existing output for the same input is archived
                     with a timestamp suffix, never overwritten in place.
    -c, --check     Check required tools and exit
    --dry-run       Print commands without executing

OUTPUTS (under --work-dir):
    <input>.hca                 raw hcatk segmentation output (clean IDs)
    <input>.HCA_TANGO.tsv       final table, original sequence IDs restored
    <input>.tango_logs/         per-sequence TANGO cmd/stdout/stderr + failures.log
    <input>.run_summary.txt

--------------------------------------------------------------------------
Changelog vs. original (~/scripts/run_hca_tango.py) — see DECISIONS.md for
full rationale and validation:

  1. CRITICAL: fasta_iter() used to keep the *entire* raw header text after
     ">" as the sequence "name" (no whitespace split). hcatk normalizes its
     own output keys to the first whitespace-delimited token, so the join
     `if name in hca` silently failed for every header carrying extra
     fields (species descriptions, coordinates, etc.) — i.e. essentially
     all real headers in this project. This produced an empty or
     near-empty task list with no error raised anywhere.
     Confirmed by direct test: hcatk given ">gene1 | Symbols: FOO | desc"
     wrote back ">gene1 78 0.003 0.013" (first token only) — the original
     script's `name` would still have been the full ">gene1 | Symbols:
     FOO | desc" string, which is not "gene1" and never matches.

  2. CRITICAL: for pipe-delimited headers with NO whitespace before the
     first "|" (e.g. this project's O. sativa CDS FASTA convention,
     "ID|chromosome:...|..."), hcatk's own whitespace-only tokenization
     cannot separate the ID from the pipe-chain either — it would treat
     the whole chain as one "id". This is the identical bug already found
     and fixed in codonw_extracted/bin/remap_seq_ids.py this project. Two
     independent third-party binaries (hcatk, tango) each do their own
     uncontrolled ID parsing, so the robust fix is the same one used for
     CodonW: remap every sequence to a clean "seqN" ID before handing the
     FASTA to either tool, and restore original IDs only in the final
     output table. This guarantees hcatk and tango always agree on IDs,
     for any species' header convention, without per-species parsing.

  3. Unescaped `name` (and `seq`) were spliced into a `shell=True` command
     STRING for TANGO. Any header/ID containing spaces, "|", quotes, or
     other shell metacharacters (again: most real headers here) would be
     reinterpreted by the shell — e.g. "|" as an actual pipe — silently
     corrupting or breaking the invocation, independent of bug #1/#2.
     Fixed by building an argv LIST and calling subprocess without
     shell=True (verified empirically to produce identical TANGO output
     to the shell=True form, with none of the parsing risk). Remapping to
     clean seqN IDs (fix #2) removes the metacharacter risk at the source
     as well — this is defense in depth, not a replacement for it.

  4. A single TANGO failure raised RuntimeError inside a worker, which
     ProcessPoolExecutor would surface via fut.result() and crash the
     entire batch — one bad sequence kills a multi-hour run. Now caught
     per-sequence, logged to <prefix>.tango_logs/failures.log with the
     original ID, and skipped; the run continues and reports a final
     failure count instead of dying partway through.

  5. Only TANGO's AGG/AMYLO tokens were kept; TURN/HELIX/HELAGG/BETA are
     parsed from the exact same stdout at zero extra cost and were being
     discarded. Now all six are captured as columns, following the same
     "capture what the tool already gives you" pattern used in the
     IUPred rebuild (mean/median/frac/etc. all derived from output IUPred
     already produces).

  6. Sequences with a trailing "*" (stop-codon symbol, common in CDS-
     derived protein FASTA) or that are empty after stripping are now
     detected and skipped with a logged reason, rather than silently fed
     to HCA/TANGO as-is or crashing.

  7. --work-dir added (default: caller's cwd, not the script's own
     location), with timestamped archiving of any pre-existing output for
     the same input prefix — same convention as codonw_extracted and
     iupred_extracted in this project, instead of overwriting in place.
--------------------------------------------------------------------------
"""

import subprocess
import sys
import warnings
import argparse
import shutil
import time
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

warnings.filterwarnings("ignore")

# ----------------------------------
# Helpers
# ----------------------------------
def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(msg):
    print(f"[{now()}] {msg}", flush=True)


def require_cmd(cmd, help_text):
    if shutil.which(cmd) is None:
        sys.exit(f"\n[ERROR] Required program '{cmd}' not found.\n{help_text}\n")


def _check_hcatk_shebang():
    """Documented, recurring bug (not hypothetical -- confirmed on a real
    install): pip's shebang-rewrite step, run when installing pyHCA's
    `hcatk` console script, can truncate the installed file's own leading
    bytes -- the file ends up starting mid-word (e.g. 'TART LICENCE ###...'
    instead of '#!<python>\\n### START LICENCE ###...'). `which hcatk`
    still resolves it fine (the file exists, is executable, is on PATH),
    but running it hands bash a Python source file with no valid shebang,
    which then tries to interpret every line as a shell command and fails
    with cryptic 'command not found' / syntax errors that look nothing
    like a missing-shebang problem unless you already know to look for it.
    Checked here once, before every run, since `which` alone can't catch
    this and the failure otherwise reads as unrelated Python errors."""
    path = shutil.which("hcatk")
    if path is None:
        return  # require_cmd("hcatk", ...) already caught and reported this
    try:
        with open(path, "rb") as f:
            head = f.read(2)
    except OSError:
        return
    if head != b"#!":
        sys.exit(
            f"\n[ERROR] 'hcatk' resolves on PATH ({path}) but its file doesn't "
            f"start with a valid shebang line -- this is a known pip/setuptools "
            f"bug where installing pyHCA truncates the script's own leading "
            f"bytes during the shebang-rewrite step, not a SeqMetrics problem.\n"
            f"Fix: open {path} and restore its first line to a real shebang "
            f"(e.g. '#!{sys.executable}') followed by '### START LICENCE' (the "
            f"line was cut mid-word, into '### S' + 'TART LICENCE'), or try:\n"
            f"  pip install -e <path to your pyHCA clone> --force-reinstall --no-deps\n"
        )


def clean_original_id(header: str) -> str:
    """Extract the true sequence identifier from a raw FASTA header.

    Same logic as codonw_extracted/bin/remap_seq_ids.py's
    original_id_from_header(), reused here deliberately so every tool in
    this pipeline agrees on what "the ID" means for a given header. Split
    on whitespace first (drops description fields), then split on "|"
    (drops pipe-chains that immediately follow the ID with no whitespace,
    e.g. this project's O. sativa convention).
    """
    token = header.split()[0] if header.split() else ""
    if not token:
        raise ValueError("Encountered empty FASTA header identifier.")
    return token.split("|")[0]


def parse_fasta(path: Path):
    """Yield (clean_id, raw_header, sequence) for each record."""
    header, seq_lines = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield clean_original_id(header), header, "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            yield clean_original_id(header), header, "".join(seq_lines)


# ----------------------------------
# Argument parsing
# ----------------------------------
parser = argparse.ArgumentParser(
    description="Run pyHCA (foldability) and TANGO (aggregation) on a FASTA file."
)
parser.add_argument("fasta", nargs="?", help="Protein FASTA file")
parser.add_argument("--jobs", type=int, default=4, help="Parallel TANGO workers (default: 4)")
parser.add_argument("--work-dir", type=str, default=None,
                     help="Output directory (default: current working directory)")
parser.add_argument("-c", "--check", action="store_true", help="Check tools and exit")
parser.add_argument("--dry-run", action="store_true", help="Print commands only")
args = parser.parse_args()

require_cmd("hcatk", "Fix:\n  conda activate hca_tango\n  pip install -e pyHCA\n")
_check_hcatk_shebang()
require_cmd(
    "tango",
    "Fix:\n  Download from https://tango.crg.es/\n"
    "  chmod +x tango\n"
    "  ln -s /path/to/tango $CONDA_PREFIX/bin/tango\n"
)

if args.check:
    log("OK: hcatk found")
    log("OK: tango found")
    log("Environment looks good")
    sys.exit(0)

if not args.fasta:
    sys.exit("ERROR: FASTA file required")

FASTA = Path(args.fasta).resolve()
PREFIX = FASTA.stem

WORK_DIR = Path(args.work_dir).resolve() if args.work_dir else Path.cwd()
WORK_DIR.mkdir(parents=True, exist_ok=True)

REMAP_FASTA = WORK_DIR / f"{PREFIX}.remapped.faa"
HCA_FILE = WORK_DIR / f"{PREFIX}.hca"
OUT_TSV = WORK_DIR / f"{PREFIX}.HCA_TANGO.tsv"
LOG_DIR = WORK_DIR / f"{PREFIX}.tango_logs"
SUMMARY_FILE = WORK_DIR / f"{PREFIX}.run_summary.txt"
FAILURES_LOG = LOG_DIR / "failures.log"

# Archive any pre-existing output for this prefix instead of overwriting.
if not args.dry_run:
    existing = [p for p in (HCA_FILE, OUT_TSV, LOG_DIR, SUMMARY_FILE) if p.exists()]
    if existing:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        archive_dir = WORK_DIR / f"{PREFIX}_prior_run_{stamp}"
        archive_dir.mkdir(parents=True, exist_ok=True)
        for p in existing:
            p.rename(archive_dir / p.name)
        log(f"INFO Archived prior output for '{PREFIX}' -> {archive_dir}")

LOG_DIR.mkdir(parents=True, exist_ok=True)

# ----------------------------------
# Run summary header
# ----------------------------------
run_start = datetime.now()
with open(SUMMARY_FILE, "w") as s:
    s.write(f"Run started: {run_start}\n")
    s.write(f"Input FASTA: {FASTA}\n")
    s.write(f"Work dir: {WORK_DIR}\n")
    s.write(f"Parallel jobs: {args.jobs}\n")
    s.write(f"Dry run: {args.dry_run}\n")
    s.write("TANGO version: not queried (interactive binary)\n\n")

# ----------------------------------
# Step 0: remap to clean seqN IDs (fixes #1/#2 at the source — see
# DECISIONS.md). Every downstream tool (hcatk, tango) only ever sees
# clean, whitespace/pipe-free IDs, so their independent internal ID
# parsing can never disagree.
# ----------------------------------
id_map = {}          # seqN -> original_header (full, for reference)
id_to_original = {}  # seqN -> clean_original_id (used in final output)
records = []          # (seqN, seq)

for idx, (clean_id, raw_header, seq) in enumerate(parse_fasta(FASTA), start=1):
    seqn = f"seq{idx}"
    id_map[seqn] = raw_header
    id_to_original[seqn] = clean_id
    records.append((seqn, seq))

if not records:
    sys.exit(f"ERROR: no sequences found in {FASTA}")

log(f"INFO Loaded {len(records)} sequences, remapped to clean IDs")

if args.dry_run:
    log(f"DRY-RUN would write remapped FASTA -> {REMAP_FASTA}")
else:
    # newline="\n": this FASTA is fed straight to hcatk/tango as real sequence
    # input -- if this ever runs under Windows Python, plain text mode writes
    # CRLF, and a stray "\r" landing at the end of every sequence line is not
    # a safe no-op for a tool parsing raw sequence characters (same failure
    # class as the confirmed CRLF/HEG-ID bug in codon_usage, but here it can
    # corrupt the sequence data itself, not just an ID lookup).
    with open(REMAP_FASTA, "w", newline="\n") as f:
        for seqn, seq in records:
            f.write(f">{seqn}\n{seq}\n")

# ----------------------------------
# Run HCA
# ----------------------------------
log(f"INFO Starting HCA on {FASTA.name} (remapped IDs)")

hca_cmd = [
    "hcatk", "segment",
    "-i", str(REMAP_FASTA),
    "-o", str(HCA_FILE),
    "-m", "domain",
    "-t", "aminoacid",
]

if args.dry_run:
    log("DRY-RUN " + " ".join(hca_cmd))
else:
    subprocess.run(hca_cmd, check=True, stderr=subprocess.DEVNULL)

# ----------------------------------
# Parse Global_HCA_score.
# Header line format for `-m domain` is:
#   >protein_id protein_length pvalue score
# (confirmed against pyHCA source, pyHCA/core/annotateHCA.py:
#  outf.write(">{} {} {:.3f} {:.3f}\n".format(prot, len(sequence), pvalue, score))
# — score is the LAST field, i.e. parts[3]; the tool's own header comment
# undersells this as a 3-field format and doesn't mention pvalue, but the
# actual write call is 4 fields and this indexing is correct.)
# ----------------------------------
hca = {}
if not args.dry_run:
    with open(HCA_FILE) as f:
        for line in f:
            # Only ">"-prefixed lines are per-sequence score lines. Without
            # this check, `domain`/`cluster` annotation lines (also >=4
            # whitespace-split fields, in `-m domain` output) get miscounted
            # as pseudo-sequences keyed "domain"/"cluster" — confirmed via
            # smoke test: a 3-sequence input reported "5 sequences parsed"
            # before this fix. Harmless for the real join (those keys never
            # match a real seqN ID) but wrong bookkeeping worth closing.
            if not line.startswith(">"):
                continue
            parts = line.split()
            if len(parts) >= 4:
                hca[parts[0].lstrip(">")] = parts[3]

log(f"INFO Parsed HCA scores for {len(hca)} sequences")

# ----------------------------------
# TANGO worker
# ----------------------------------
TANGO_FIELDS = ["AGG", "AMYLO", "TURN", "HELIX", "HELAGG", "BETA"]


def run_tango(task):
    seqn, seq = task

    seq = seq.rstrip("*").strip()
    if not seq:
        return seqn, None, "empty sequence after stripping stop-codon symbol"

    cmd = ["tango", seqn, "ct=N", "nt=N", "ph=7.0", "te=298", "io=0.02", f"seq={seq}"]
    (LOG_DIR / f"{seqn}.tango.cmd").write_text(" ".join(cmd) + "\n")

    if args.dry_run:
        return seqn, {k: "DRY" for k in TANGO_FIELDS}, None

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=LOG_DIR,
            stdin=subprocess.DEVNULL,
            timeout=120,
        )
    except subprocess.TimeoutExpired:
        return seqn, None, "TANGO timed out after 120s"

    (LOG_DIR / f"{seqn}.tango.out").write_text(result.stdout)
    (LOG_DIR / f"{seqn}.tango.err").write_text(result.stderr)

    if result.returncode != 0:
        return seqn, None, f"TANGO exited {result.returncode} (see {seqn}.tango.err)"

    values = {k: "0.0" for k in TANGO_FIELDS}
    tokens = result.stdout.split()
    for i, tok in enumerate(tokens):
        if tok in values and i + 1 < len(tokens):
            values[tok] = tokens[i + 1]

    if all(v == "0.0" for v in values.values()) and result.stdout.strip():
        # Non-empty output but none of the expected labels were found —
        # TANGO's output format has drifted from what this script expects.
        # Report it rather than silently writing all-zero rows.
        return seqn, None, f"unrecognized TANGO output format: {result.stdout.strip()[:200]!r}"

    return seqn, values, None


# ----------------------------------
# Parallel TANGO execution
# ----------------------------------
log(f"INFO Starting TANGO with {args.jobs} workers")

tasks = [(seqn, seq) for seqn, seq in records if seqn in hca]
skipped_no_hca = len(records) - len(tasks)
if skipped_no_hca:
    log(f"WARNING {skipped_no_hca} sequences had no HCA score and will be skipped for TANGO")

total = len(tasks)
start_time = time.time()
results = {}
failures = []

with ProcessPoolExecutor(max_workers=args.jobs) as exe:
    futures = {exe.submit(run_tango, t): t[0] for t in tasks}

    for i, fut in enumerate(as_completed(futures), 1):
        seqn, values, err = fut.result()
        if err is not None:
            failures.append((seqn, err))
        else:
            results[seqn] = values

        if i % 1000 == 0 or i == total:
            elapsed = time.time() - start_time
            rate = i / elapsed if elapsed > 0 else 0
            remaining = (total - i) / rate if rate > 0 else 0
            log(
                f"PROGRESS {i}/{total} ({(i/total)*100:.1f}%) | "
                f"elapsed {elapsed/60:.1f} min | ETA {remaining/60:.1f} min | "
                f"failures so far: {len(failures)}"
            )

if failures:
    with open(FAILURES_LOG, "w") as f:
        for seqn, err in failures:
            f.write(f"{seqn}\t{id_to_original.get(seqn, seqn)}\t{err}\n")
    log(f"WARNING {len(failures)} sequences failed TANGO — see {FAILURES_LOG}")

# ----------------------------------
# Write final TSV — original sequence IDs restored here, not before.
# ----------------------------------
# newline="\n": this project routinely joins feature tables by exact ID
# match across stages (genes.tsv, composite headers, etc.) -- a trailing
# CRLF-corrupted ID would silently fail such a join with no error anywhere.
with open(OUT_TSV, "w", newline="\n") as out:
    out.write("Seq_ID\tGlobal_HCA_score\t" + "\t".join(f"TANGO_{k}" for k in TANGO_FIELDS) + "\n")
    for seqn in sorted(results, key=lambda s: int(s.replace("seq", ""))):
        values = results[seqn]
        orig_id = id_to_original[seqn]
        row = "\t".join(values[k] for k in TANGO_FIELDS)
        out.write(f"{orig_id}\t{hca[seqn]}\t{row}\n")

run_end = datetime.now()
with open(SUMMARY_FILE, "a") as s:
    s.write(f"\nTotal sequences: {len(records)}\n")
    s.write(f"Skipped (no HCA score): {skipped_no_hca}\n")
    s.write(f"TANGO attempted: {total}\n")
    s.write(f"TANGO succeeded: {len(results)}\n")
    s.write(f"TANGO failed: {len(failures)}\n")
    s.write(f"\nRun finished: {run_end}\n")
    s.write(f"Elapsed time: {run_end - run_start}\n")

log(f"DONE TSV written -> {OUT_TSV}")
log(f"DONE Logs -> {LOG_DIR}")
log(f"DONE Summary -> {SUMMARY_FILE}")
