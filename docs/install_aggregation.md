# Install: aggregation (HCA + TANGO)

**Correction (2026-09-20): `pip install hcatk` does not work.** Neither
`hcatk` nor `pyhca` exists on public PyPI -- checked directly against
PyPI's own JSON API (`pypi.org/pypi/hcatk/json` -> 404,
`pypi.org/pypi/pyhca/json` -> 404), not just "couldn't get it installed."
The `hcatk` binary comes from **pyHCA** (Bitard-Feildel & Faure,
bornberglab.org, MIT-licensed), and this lab's own working install is an
*editable install from a locally cloned git repo*, not a pip package:

**Clone pyHCA outside this repo, not inside it.** Confirmed the hard way:
running the clone command from inside a freshly-cloned SeqMetrics puts
`pyHCA/` inside SeqMetrics' own working tree -- it's not gitignored, so it
shows up as untracked clutter and a careless `git add -A` could commit a
whole separate third-party repo into this one. Clone it as a sibling, same
convention as this project's other third-party tools:
```
cd ~/Bioinformatics   # or wherever you keep pyHCA/LOCALIZER/etc. as siblings, not inside SeqMetrics/
git clone https://github.com/DarkVador-HCA/pyHCA.git
cd pyHCA
pip install -e .
which hcatk   # confirms the console-script entry point pyHCA's setup.py defines
```

**Real dependencies `pip install -e .` doesn't pull in, confirmed on a real
install -- install all three together, not one at a time.** pyHCA's source
imports `Bio.SubsMat.MatrixInfo` (removed in modern Biopython -- a bare
`pip install biopython` resolves too new and fails with
`ModuleNotFoundError: No module named 'Bio.SubsMat'`) and `pandas`
directly (`disorderHCA.py`), neither declared as an install-time
dependency by pyHCA's own `setup.py`. Matches this lab's own prior pin
for this exact env (`htlcp/BROAD_PIPELINE_SETUP.md`), which installs all
three in one line:
```
pip install pandas lightgbm "biopython==1.79"
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

**Second known issue -- confirmed on a real install (2026-09-21), not just
referenced from a prior lab doc.** `pyHCA/core/ioHCA.py`'s
`read_multifasta_it()` opens the input FASTA with `open(path, "rU")` --
`'U'` mode was removed entirely in Python 3.11, which is exactly what this
doc's own env creation command uses. Real error, confirmed:
```
ValueError: invalid mode: 'rU'
  File ".../pyHCA/core/ioHCA.py", line 78, in read_multifasta_it
    with open(path, "rU") as handle:
```
`hcatk` writes its output file's header comment *before* this call, so a
`.hca` output file with only a header and no real content is this exact
error, not a different problem. Fix -- confirmed safe, not a workaround:
`'rU'` is semantically identical to plain `'r'` in Python 3 (universal
newline handling is the text-mode default already), so this is a genuine
correction to code written for Python 2/early 3, not a functional change:
```
grep -rln '"rU"' <path to your pyHCA clone>/pyHCA/ | xargs sed -i 's/"rU"/"r"/g'
```
(Confirmed two occurrences, both in `ioHCA.py`, on the version cloned from
`github.com/DarkVador-HCA/pyHCA` -- check your own clone in case a
different fork or version has more.)

pyHCA is MIT-licensed, so this part is genuinely freely redistributable
-- unlike TANGO below. **TANGO is not open** -- it's distributed under an
academic-use license by its original authors (Rousseau, Serrano,
Schymkowitz lab); you register and download it directly from them. This
repo ships the wrapper script only (`run_hca_tango.py`), never the TANGO
binary itself.

Registering at tango.crg.es emails you a download link -- often wrapped
in a mail-security redirect (e.g. URL Defense) that a bare `wget`/`curl`
can't follow (confirmed: it silently returns the wrapper's HTML landing
page, not the file, with no error -- check `file` on whatever you
downloaded before trusting it). Open the link in a real browser instead,
pick the correct architecture (e.g. "Executable Tango - Linux 64bits"),
and `scp` the real download to your target machine. It'll be a small
`.zip` (confirmed ~78KB, not hundreds of MB -- this is a small compiled
binary, don't expect a large download):
```
unzip tango2_3_1.linux64.zip
chmod +x tango_x86_64_release
cp tango_x86_64_release "$CONDA_PREFIX/bin/tango"   # with hca_tango active
which hcatk tango   # confirm both resolve
tango                # no args -- should prompt interactively (Y/N, then a
                      # filename), not crash; Ctrl-C out once confirmed
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
