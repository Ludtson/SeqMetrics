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
