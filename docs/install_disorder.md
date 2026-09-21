# Install: disorder (IUPred3 + ANCHOR2)

**IUPred3 is academic-license**, not open-source -- you register at the
authors' site (iupred3.elte.hu) and download `iupred3.py` +
`iupred3_lib.py` yourself. This repo ships the wrapper script only
(`build_iupred_features.py`), never the IUPred3 engine itself.

```
# after obtaining iupred3.py + iupred3_lib.py from the authors:
export PYTHONPATH="$PYTHONPATH:/path/to/iupred3"
python -c "import iupred3_lib; print('OK')"
```

This lab's existing WSL setup already has a working install under the
conda env name `iupred3`, with the engine files at
`~/miniconda3/envs/iupred3/apps/iupred3/{iupred3.py,iupred3_lib.py}` --
confirmed present on 2026-09-20. As with TANGO, that install exists
because of a prior license agreement; don't assume it's freely copyable
to a different machine without checking the license terms yourself.

If your conda env for this module is named something other than
`iupred3`, pass it to the orchestrator explicitly:
`--module-env disorder=<your_env_name>` -- but note the PYTHONPATH
addition documented above is currently hardcoded in `run_features.py`
(`PYMODULE_PYTHONPATH_WSL`) to this env's conventional path; a
differently-named env still needs its engine files reachable on
PYTHONPATH some other way (e.g. installed straight into that env's own
site-packages) since the override doesn't carry a matching path override
yet.

Run: `python build_iupred_features.py proteins.faa output.tsv` (this
repo's actual entry point is `batch_iupred_features_cysexcl.py`, same
CLI). Three variants exist in the source project (plain /
cysteine-excluded / cysteine-stripped) -- this repo ships the
cysteine-excluded (`batch_iupred_features_cysexcl.py`) variant as the one
actually used, per the documented manuscript methods (cysteine removal
corrects for disulfide-bond bias in disorder scores). See `DECISIONS.md`
for the full rationale.
