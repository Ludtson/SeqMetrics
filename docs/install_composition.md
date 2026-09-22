# Install: composition (EMBOSS PEPSTATS)

EMBOSS is open-source (GPL), no registration required.

```
conda create -n em_boss -c bioconda -c conda-forge emboss
conda activate em_boss
which pepstats
```

Run: `modules/composition/pepstats-run <fasta> <output-dir>` (see
`modules/composition/README.md` for the full flag list). Only the plain
`MolePct` output is used; the Dayhoff-stat variant is not (it is a
per-amino-acid linear rescaling of the same numbers and adds no
independent signal for this comparison).
