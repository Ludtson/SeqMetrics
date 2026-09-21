# Install: tm_domain (DeepTMHMM2)

**The official DTU/BioLib DeepTMHMM requires registration** for any local
(non-Docker) use -- confirmed via a real report of someone else hitting
this (a Bowman Lab blog post), not just the DTU site's own claims. Docker
itself is also a real obstacle on RHEL specifically, independent of the
registration requirement.

This module uses **`fteufel/DeepTMHMM2`** instead
(github.com/fteufel/DeepTMHMM2) -- a separate, ungated, pip-installable
reimplementation by one of the original DeepTMHMM paper's authors. No
license notice found in that repo (not the same as confirmed-clear
license status -- worth checking yourself before redistributing anything
built from it, same caution as any other tool in this project). CLI tool
is `dtm2`, defaults to CPU, no GPU/Docker required.

## Isolated venv, not conda -- and not the system Python

```
python -m venv .venvs/deeptmhmm2
.venvs/deeptmhmm2/Scripts/pip install git+https://github.com/fteufel/DeepTMHMM2.git   # Windows
.venvs/deeptmhmm2/bin/pip install git+https://github.com/fteufel/DeepTMHMM2.git        # Linux/macOS
```

**Do not `pip install` this into an active/system Python without a venv.**
Confirmed the hard way this session: doing exactly that once silently
upgraded numpy (1.24->2.4), pandas (2.3->**3.0**, a major version jump),
biopython, and matplotlib in the system interpreter, with pip itself
warning about a resulting scipy/numpy incompatibility -- a real
regression risk for any other project sharing that interpreter. Full
remediation required reverting six packages to their pinned versions and
uninstalling the entire torch/lightning/fair-esm dependency stack by
hand. The isolated venv exists specifically so this never happens again.

First run downloads model checkpoints (5x `topology_N.ckpt` +
5x `memtype_N.ckpt` to `~/.cache/deeptmhmm2/`) and the ESM2 650M language
model (`~/.cache/torch/hub/checkpoints/`, a few hundred MB) -- expect a
slow first invocation, fast ones after.

## Resolution: "venv" kind, not conda

Unlike every other module in this repo, this one isn't found via
`--module-env <conda_env_name>` -- `run_features.py`'s `resolve_venv()`
looks for the venv itself at the documented default path
(`SeqMetrics/.venvs/deeptmhmm2/`, relative to this repo, so it works
regardless of the caller's cwd) or at an explicit
`--module-env tm_domain=<venv_dir>` override (a directory path, not an
env name, for this module specifically).

Run: `python modules/tm_domain/run_dtm2.py -i proteins.faa -o output.tsv
--dtm2 <path to dtm2>` -- this wrapper is self-contained: it runs `dtm2`
into its own temp directory, parses all four of its output files
(`TMRs.gff3`, `membrane_types.tsv`, `predicted_topologies.3line`,
`predictions.json`) into one flat table, and cleans up after itself.

Real output columns verified against real data: `Seq_ID, N_TMRs,
Residues_<topology_label> (one column per topology label DeepTMHMM2
actually emitted, e.g. inside/outside/signal/TMhelix),
MembraneProb_<membrane_type> (one per membrane-type category)`. Spot
check: the tracked locus `AT5G15843.1` (this project's own real DNG,
verified elsewhere against `genes.tsv`) predicts 0 TMRs here, matching
`genes.tsv`'s own `transmembrane_present=0` for that exact gene.

**The `Residues_*` columns are not fixed** -- they're built from whatever
topology labels actually appear in a given run's results, so two separate
runs on different sequence batches can produce tables with a different
number of `Residues_*` columns (e.g. a batch containing a real membrane
protein gets a `Residues_TMhelix` column that a batch of only soluble
proteins won't have). Concatenating such files naively will misalign
columns -- use `utils/combine_datasets.py` (see `utils/README.md`) to
join them by `Seq_ID` properly instead.
