# Install: aggregation (HCA + TANGO)

**Correction (2026-09-20): `pip install hcatk` does not work.** Neither
`hcatk` nor `pyhca` exists on public PyPI -- checked directly against
PyPI's own JSON API (`pypi.org/pypi/hcatk/json` -> 404,
`pypi.org/pypi/pyhca/json` -> 404), not just "couldn't get it installed."
The `hcatk` binary comes from **pyHCA** (Bitard-Feildel & Faure,
bornberglab.org, MIT-licensed), and this lab's own working install is an
*editable install from a locally cloned git repo*, not a pip package:

```
git clone https://github.com/DarkVador-HCA/pyHCA.git
cd pyHCA
pip install -e .
which hcatk   # confirms the console-script entry point pyHCA's setup.py defines
```

**Real, recurring bug, confirmed on a real install -- `which hcatk` succeeding
is not enough.** pip's shebang-rewrite step (it rewrites the installed
script's first line to point at the target interpreter) can truncate the
file's own leading bytes in the process -- the installed `hcatk` ends up
starting mid-word, e.g. `TART LICENCE ###...` instead of the real
`#!/usr/bin/env python3\n\n# START LICENCE ###...` (shebang, a blank line,
then a single-`#` comment -- verified against a real working install, not
guessed). `which` still finds it (the file exists, is executable, is on
PATH), but running it hands bash a Python source file with no valid
shebang, which tries to interpret every line as a shell command:
```
hcatk: line 1: TART: command not found
hcatk: line 32: syntax error near unexpected token `('
```
This reads like an unrelated Python/shell problem unless you already know
to check the shebang. `run_hca_tango.py` now checks this itself before
every run (`_check_hcatk_shebang()` reads the file's first 2 bytes) rather
than relying on `which` alone. Fix, verified against a real working
install:
```
HCATK=$(which hcatk)
{ printf '#!/usr/bin/env python3\n\n# S'; cat "$HCATK"; } > /tmp/hcatk_fixed
mv /tmp/hcatk_fixed "$HCATK" && chmod +x "$HCATK"
```
Or reinstall so pip regenerates it cleanly (not confirmed to actually avoid
the bug, since the bug is in pip's own rewrite step, not the source file --
the manual fix above is the one actually verified to work):
```
pip install -e <path to your pyHCA clone> --force-reinstall --no-deps
```

**Second known issue, from this lab's prior pyHCA integration (`htlcp`'s own
`BROAD_PIPELINE_SETUP.md`), not independently re-confirmed against the
current pyHCA source**: a Python `'rU'` file-open-mode compatibility
problem. `'rU'` (universal-newlines mode) was deprecated in Python 3.11 and
removed entirely in 3.12 -- if pyHCA's source opens any file with `'rU'`, it
will raise a `ValueError` on the `python=3.11` env this doc's own install
command creates. If `hcatk` fails with that error after the shebang fix
above, this is the likely cause -- either patch the offending `open()` call
to `'r'` in your pyHCA clone, or create the env with an older Python
(`python=3.8`, before this became an error) as a faster workaround.

pyHCA is MIT-licensed, so this part is genuinely freely redistributable
-- unlike TANGO below. **TANGO is not open** -- it's distributed under an
academic-use license by its original authors (Rousseau, Serrano,
Schymkowitz lab); you register and download it directly from them. This
repo ships the wrapper script only (`run_hca_tango.py`), never the TANGO
binary itself.

```
# obtain the tango binary yourself, place it on PATH
which hcatk tango   # confirm both resolve
```

This lab's existing WSL setup already has both working under the conda
env name `hca_tango`, confirmed via `which hcatk tango` and `pip show
pyHCA` (confirming the editable-install path above) on 2026-09-20 --
meaning a working TANGO binary already exists there from a prior,
properly-licensed install. Ask before assuming that install is
redistributable to a new machine; the license is per-installation for
TANGO specifically, not transferable by copying the binary (pyHCA itself
has no such restriction).

If your conda env for this module is named something other than
`hca_tango`, pass it to the orchestrator explicitly:
`--module-env aggregation=<your_env_name>`.

Run: `python run_hca_tango.py input.faa --jobs N --work-dir DIR`. See
`DECISIONS.md` in this module's directory for the real bugs already
found and fixed in this wrapper (a silent header-matching failure that
likely nulled out most past runs on real multi-species FASTA).
