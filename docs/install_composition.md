# Install: composition (EMBOSS PEPSTATS)

EMBOSS is open-source (GPL) -- freely installable, not something you
need to register for.

```
conda create -n em_boss -c bioconda -c conda-forge emboss
conda activate em_boss
which pepstats   # confirm it resolves
```

This lab's existing WSL setup already has this working under the conda
env name `em_boss`, confirmed via `which pepstats` on 2026-09-20.

Run: `modules/composition/pepstats-run <fasta> <output-dir>` (see
`modules/composition/README.md`, copied in from the source project, for
the full flag list). Only the plain `MolePct` output is needed for this
project -- the Dayhoff-stat variant is a per-amino-acid linear rescaling
of the same numbers and adds no new signal for a same-amino-acid,
different-time-point comparison (confirmed 2026-09-20, not run here).
