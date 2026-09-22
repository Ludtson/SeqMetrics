# Install: coding_potential (CPAT)

CPAT is open-source and pip-installable, no registration needed.

```
conda create -n cpat python=3.11 -y
conda activate cpat
pip install CPAT
conda install -c conda-forge r-base -y
which cpat make_hexamer_tab make_logitModel Rscript
```

Rscript (R) is needed on PATH even though you never invoke it directly
— CPAT's scoring step generates and runs a small R script internally.

## Two stages

### Stage 1 — reference-building (training), once per species

Builds a species-specific hexamer usage table and a trained logistic
model from known coding and non-coding sequences. CPAT ships no
default model.

Build the non-coding training set with `get_noncoding_data.py`
(`modules/coding_potential/`), from a genome + GFF3 + protein FASTA:
```
python get_noncoding_data.py \
    -p protein.faa -g genome.fa -a annotation.gff3 \
    -o out_dir -r <species_prefix>
# -> out_dir/<protein_basename>_nc.faa
```
Finds every non-CDS genome segment (intergenic/intronic), keeps only
ones ≥200nt (`--min-noncoding-length`), then rejects any segment
containing a ≥250nt ORF in any of 6 frames (`--min-orf-length`) or
whose 6-frame translation matches a 20-mer from the known protein set.
`-n/--subset-cds-count` optionally draws a random CDS subset for the
`coding_cds.fa` side in the same run.

Then:
```
make_hexamer_tab -c coding_cds.fa -n noncoding.fa > hexamer.tsv
make_logitModel -x hexamer.tsv -c coding_cds.fa -n noncoding.fa -o model
# -> model.logit.RData
```

**Before reusing any pre-built model across projects**, verify it was
trained against the same genome annotation version your own data
uses. A model trained on one annotation version, scored against
another for the "same" species, will not error — it will silently
score against mismatched gene/isoform calls.

### Stage 2 — scoring, against an already-built reference

This is SeqMetrics' own job:
```
run_features.py --modules coding_potential --module-ref coding_potential=<dir containing hexamer.tsv and logit.RData>
```
Runs `cpat -g target.fa -d logit.RData -x hexamer.tsv -o prefix`, then
restores original sequence IDs (CPAT uppercases every ID it writes and
mishandles some pipe-delimited headers — see `repair_cpat_ids.py`'s
module docstring). Unlike codonW, CPAT has no working-directory
requirement — every path handed to `cpat` must be absolute, which
`run_features.py` guarantees for `--module-ref`.

`hexamer.tsv` must be named exactly that. `logit.RData` is resolved
automatically — `make_logitModel -o model` produces `model.logit.RData`,
not a bare `logit.RData`, and `run_features.py` accepts that directly
as long as it is the only `*.logit.RData` file in the `--module-ref`
directory. If more than one is present, name the one you want to use
exactly `logit.RData`.

Output columns: `ID, mRNA, ORF_strand, ORF_frame, ORF_start, ORF_end,
ORF, Fickett, Hexamer, Coding_prob`, plus `Gene_ID` (restored from
CPAT's own `seq_ID`).
