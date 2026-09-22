# Install: aggregation (HCA + TANGO)

`pip install hcatk` does not work — neither `hcatk` nor `pyhca` exists
on PyPI. The `hcatk` binary comes from **pyHCA** (Bitard-Feildel &
Faure, bornberglab.org, MIT-licensed), installed as an editable
install from a locally cloned git repo.

```
conda create -n hca_tango python=3.11 -y
conda activate hca_tango
```

Clone pyHCA outside this repo, not inside it — cloning it inside a
SeqMetrics checkout leaves it untracked and unguarded against an
accidental `git add -A`:
```
cd ~/Bioinformatics   # or wherever you keep third-party tools as siblings, not inside SeqMetrics/
git clone https://github.com/DarkVador-HCA/pyHCA.git
cd pyHCA
pip install -e .
which hcatk
```

**Real dependencies `pip install -e .` does not pull in** — install
all three together:
```
pip install pandas lightgbm "biopython==1.79"
```
pyHCA's source imports `Bio.SubsMat.MatrixInfo`, removed in modern
Biopython; a bare `pip install biopython` resolves too new and fails
with `ModuleNotFoundError: No module named 'Bio.SubsMat'`. `pandas` is
also a direct, undeclared dependency (`disorderHCA.py`).

**Known pip/setuptools bug**: `which hcatk` succeeding is not enough.
pip's shebang-rewrite step can truncate the installed script's own
leading bytes, leaving it starting mid-word (e.g. `TART LICENCE
###...` instead of `#!/usr/bin/env python3\n\n# START LICENCE
###...`). The file still resolves on PATH and is still executable, but
running it hands bash a Python source file with no valid shebang:
```
hcatk: line 1: TART: command not found
hcatk: line 32: syntax error near unexpected token `('
```
`run_hca_tango.py` checks this itself before every run
(`_check_hcatk_shebang()`). Fix:
```
HCATK=$(which hcatk)
{ printf '#!/usr/bin/env python3\n\n# S'; cat "$HCATK"; } > /tmp/hcatk_fixed
mv /tmp/hcatk_fixed "$HCATK" && chmod +x "$HCATK"
```

**Second known issue**: `pyHCA/core/ioHCA.py`'s `read_multifasta_it()`
opens the input FASTA with `open(path, "rU")` — `'U'` mode was removed
in Python 3.11. Error:
```
ValueError: invalid mode: 'rU'
```
`'rU'` is semantically identical to plain `'r'` in Python 3 (universal
newline handling is already the text-mode default), so this is a safe
fix, not a workaround:
```
grep -rln '"rU"' <path to your pyHCA clone>/pyHCA/ | xargs sed -i 's/"rU"/"r"/g'
```

pyHCA is MIT-licensed and freely redistributable. **TANGO is not** —
it is distributed under an academic-use license by its original
authors (Rousseau, Serrano, Schymkowitz lab); register and download it
directly from tango.crg.es. This repo ships the wrapper script only
(`run_hca_tango.py`), never the TANGO binary.

Registration emails a download link often wrapped in a mail-security
redirect that a bare `wget`/`curl` cannot follow (it silently returns
the wrapper's HTML landing page, not the file — check `file` on
whatever you downloaded). Open the link in a browser instead, pick the
correct architecture, and transfer the real download to your target
machine:
```
unzip tango2_3_1.linux64.zip
chmod +x tango_x86_64_release
cp tango_x86_64_release "$CONDA_PREFIX/bin/tango"   # with hca_tango active
which hcatk tango
tango   # no args -- should prompt interactively (Y/N, then a filename), not crash
```

If your conda env for this module is named something other than
`hca_tango`, pass it explicitly: `--module-env
aggregation=<your_env_name>`.

Run: `python run_hca_tango.py input.faa --jobs N --work-dir DIR`. See
`DECISIONS.md` in this module's directory for the header-matching fix
this wrapper depends on.
