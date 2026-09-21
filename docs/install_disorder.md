# Install: disorder (IUPred3 + ANCHOR2)

**IUPred3 is academic-license**, not open-source -- you register at the
authors' site (iupred3.elte.hu) and download the engine (`iupred3.py`,
`iupred3_lib.py`, and its `data/` directory -- the energy matrices and
histograms it needs at runtime, not optional) yourself. This repo ships
the wrapper script only (`build_iupred_features.py`), never the IUPred3
engine itself. The archive also includes `LICENSE` and a `README` --
worth checking the `README` names IUPred3's real authors (Erdos, Pajkos,
Dosztanyi) as a sanity check that you have the genuine release.

**Primary method, verified on a real install -- drop the engine straight
into your conda env's own site-packages, needs no PYTHONPATH at all:**
```
conda create -n iupred3 python=3.11 scipy -y

SITE_PACKAGES=$(conda run -n iupred3 python -c 'import site; print(site.getsitepackages()[0])')
echo "$SITE_PACKAGES"   # confirm it's a real path before copying anything

cp /path/to/iupred3/iupred3_lib.py /path/to/iupred3/iupred3.py "$SITE_PACKAGES/"
cp -r /path/to/iupred3/data "$SITE_PACKAGES/"

conda run -n iupred3 python -c "import iupred3_lib; print('OK')"
```

**Alternative method, if you'd rather not touch the env's site-packages**:
`export PYTHONPATH` in your own shell only guarantees the module is
importable in *that shell* -- `run_features.py` actually resolves this
module via a subprocess (`conda run -n ENV python ...`), and whether that
subprocess inherits your shell's PYTHONPATH isn't something to assume.
Use `--module-pythonpath disorder=/path/to/iupred3` instead (added
specifically for this) -- the orchestrator threads it into the resolution
check and the real run itself, not just a manual sanity check you run
once and hope matches:
```
python run_features.py --modules disorder ... --module-pythonpath disorder=/path/to/iupred3
```
This is a separate, purely opt-in mechanism from the WSL-only
`PYMODULE_PYTHONPATH_WSL` hardcoded path (this lab's own WSL setup, see
below) -- passing `--module-pythonpath` never affects that existing
WSL behavior.

This lab's existing WSL setup already has a working install under the
conda env name `iupred3`, with the engine files at
`~/miniconda3/envs/iupred3/apps/iupred3/{iupred3.py,iupred3_lib.py}` --
confirmed present on 2026-09-20, and `iupred3_lib.py` specifically
(the only file SeqMetrics actually imports) confirmed byte-identical
(md5sum) to a separately-obtained copy used on a different machine, if
you're ever unsure whether two copies are the same release. As with
TANGO, that install exists because of a prior license agreement; don't
assume it's freely copyable to a different machine without checking the
license terms yourself.

If your conda env for this module is named something other than
`iupred3`, pass it to the orchestrator explicitly:
`--module-env disorder=<your_env_name>`.

Run: `python build_iupred_features.py proteins.faa output.tsv` (this
repo's actual entry point is `batch_iupred_features_cysexcl.py`, same
CLI). Three variants exist in the source project (plain /
cysteine-excluded / cysteine-stripped) -- this repo ships the
cysteine-excluded (`batch_iupred_features_cysexcl.py`) variant as the one
actually used, per the documented manuscript methods (cysteine removal
corrects for disulfide-bond bias in disorder scores). See `DECISIONS.md`
for the full rationale.
