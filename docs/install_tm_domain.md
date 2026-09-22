# Install: tm_domain (DeepTMHMM2)

The official DTU/BioLib DeepTMHMM requires registration for local
(non-Docker) use. This module uses `fteufel/DeepTMHMM2`
(github.com/fteufel/DeepTMHMM2) instead, an ungated reimplementation.
No license notice is present in that repo — verify redistribution
rights yourself before relying on it. CLI tool is `dtm2`, defaults to
CPU.

**Isolated venv, not conda, not the system Python:**
```
python -m venv .venvs/deeptmhmm2
.venvs/deeptmhmm2/Scripts/pip install git+https://github.com/fteufel/DeepTMHMM2.git   # Windows
.venvs/deeptmhmm2/bin/pip install git+https://github.com/fteufel/DeepTMHMM2.git        # Linux/macOS
```

Do not `pip install` this into an active/system Python without a venv
— its dependency stack (torch, fair-esm, lightning) can upgrade
numpy/pandas/biopython/matplotlib in that interpreter, breaking any
other project sharing it.

The default install pulls the full GPU-enabled `torch` build (several
GB of CUDA packages even without a GPU). For CPU-only:
```
.venvs/deeptmhmm2/bin/pip install torch --index-url https://download.pytorch.org/whl/cpu
.venvs/deeptmhmm2/bin/pip install git+https://github.com/fteufel/DeepTMHMM2.git
```

First run downloads model checkpoints (5x `topology_N.ckpt` + 5x
`memtype_N.ckpt`) and the ESM2 650M language model — expect a slow
first invocation, fast ones after.

**Resolution is by venv path, not conda env name** — `run_features.py`
looks for `SeqMetrics/.venvs/deeptmhmm2/` by default, or an explicit
`--module-env tm_domain=<venv_dir>` (a directory path, not an env
name, for this module specifically).

Run: `python modules/tm_domain/run_dtm2.py -i proteins.faa -o
output.tsv --dtm2 <path to dtm2>`.

Output columns: `Seq_ID, N_TMRs, Residues_<topology_label>` (one per
topology label actually seen in the run, e.g.
inside/outside/signal/TMhelix), `MembraneProb_<membrane_type>` (one
per membrane-type category). The `Residues_*`/`MembraneProb_*` columns
are not fixed — different batches can produce different column sets
depending on what topology/membrane labels appear. Use
`utils/combine_datasets.py` (see `utils/README.md`) to join files by
`Seq_ID` rather than concatenating them directly.
