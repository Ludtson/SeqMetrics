# Install: disorder (IUPred3 + ANCHOR2)

IUPred3 is academic-license, not open-source. Register at
iupred3.elte.hu and download the engine (`iupred3.py`, `iupred3_lib.py`,
and its `data/` directory — the energy matrices and histograms it needs
at runtime). This repo ships the wrapper script only
(`build_iupred_features.py`), never the engine.

**Primary method — drop the engine into the env's own site-packages,
no PYTHONPATH needed:**
```
conda create -n iupred3 python=3.11 scipy -y

SITE_PACKAGES=$(conda run -n iupred3 python -c 'import site; print(site.getsitepackages()[0])')

cp /path/to/iupred3/iupred3_lib.py /path/to/iupred3/iupred3.py "$SITE_PACKAGES/"
cp -r /path/to/iupred3/data "$SITE_PACKAGES/"

conda run -n iupred3 python -c "import iupred3_lib; print('OK')"
```

**Alternative**, if you don't want to touch the env's site-packages:
pass `--module-pythonpath disorder=/path/to/iupred3` to
`run_features.py`. Setting `PYTHONPATH` in your own shell is not
sufficient — `run_features.py` resolves this module through a
subprocess (`conda run -n ENV python ...`) that does not reliably
inherit your shell's `PYTHONPATH`.

If your env isn't named `iupred3`, pass `--module-env
disorder=<env_name>`.

Run: `python build_iupred_features.py proteins.faa output.tsv` (entry
point is `batch_iupred_features_cysexcl.py`). Uses the
cysteine-excluded variant, which corrects for disulfide-bond bias in
disorder scores; see `DECISIONS.md` for the other two variants shipped
in the source project and why they're not used here.
