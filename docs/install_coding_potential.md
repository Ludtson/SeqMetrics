# Install: coding_potential (CPAT)

CPAT is open-source and pip-installable, no registration needed. Verified
real: `CPAT` exists on PyPI (checked directly against
`pypi.org/pypi/CPAT/json`, not just assumed).

```
conda create -n cpat python=3.11 -y
conda activate cpat
pip install CPAT
conda install -c conda-forge r-base -y
which cpat make_hexamer_tab make_logitModel Rscript   # confirm all four resolve
```

Rscript (R) is needed on PATH even though you never invoke it directly --
CPAT's scoring step generates and runs a small R script internally
(`predict()` against the trained logistic model). `conda install -c
conda-forge r-base` is the confirmed-working way to get it into the same
env, verified on a real install (2026-09-21) -- all four binaries resolved
cleanly.

This lab's existing WSL setup already has this working under the conda
env name `cpat` (CPAT 3.0.5), confirmed on 2026-09-20.

## Two stages, same shape as codon_usage

**Stage 1 -- reference-building (training), once per species.** Builds
a species-specific hexamer usage table and a trained logistic model from
known coding and non-coding sequences -- CPAT ships no default model, it
cannot score anything without one:
```
make_hexamer_tab -c coding_cds.fa -n noncoding.fa > hexamer.tsv
make_logitModel -x hexamer.tsv -c coding_cds.fa -n noncoding.fa -o model
# -> model.logit.RData
```
**Correction (2026-09-21): `get_noncoding_data.py` was never actually
missing.** A prior version of this doc claimed it didn't exist anywhere
in this project's files, "confirmed via a full search" -- that claim was
simply wrong, not stale; the script (`modules/coding_potential/
get_noncoding_data.py`) is real, complete, and already in this repo. It
builds `noncoding.fa` straight from a genome + GFF3 + protein FASTA --
exactly the inputs `codon_usage`'s HMM-HEG reference-building already
needs, so no separate data-gathering step is required:
```
python get_noncoding_data.py \
    -p protein.faa -g genome.fa -a annotation.gff3 \
    -o out_dir -r <species_prefix>
# -> out_dir/<protein_basename>_nc.faa
```
Internally: finds every non-CDS genome segment (intergenic/intronic),
keeps only ones ≥200nt (`--min-noncoding-length`), then rejects any
segment containing a ≥250nt ORF in any of 6 frames
(`--min-orf-length`) or whose 6-frame translation matches a 20-mer from
the known protein set -- two independent checks against accidentally
including real coding sequence in the "noncoding" training set. Also
supports `-n/--subset-cds-count` to draw a random CDS subset for the
`coding_cds.fa` side in the same run.

Three species have trained models from the source (older) project's prior
run (Arabidopsis thaliana, Brassica rapa, Oryza sativa) --
`2026fall-mla-chapter/feature_outputs/cpat_output/<species>/*_cpat_output/{1_hexamer,2_model}/`.
Confirmed by direct test: SeqMetrics' scoring output using Athaliana's
existing model reproduces the historical run's `Coding_prob`/`Fickett`/
`Hexamer` values byte-for-byte across 25,278 overlapping genes.

**Before reusing any pre-built model across projects, verify it was
trained against the same genome annotation version your own data uses.**
A model trained on one annotation (e.g. an older TAIR/RefSeq/Ensembl
release) scored against a different one for the "same" species won't
error -- it'll just silently score against mismatched gene/isoform calls.
This is a real, easy-to-miss failure mode worth checking explicitly, not
something SeqMetrics itself can detect or warn about automatically.

**Stage 2 -- scoring, against an already-built reference.** This is
SeqMetrics' own job (`run_features.py --modules coding_potential
--module-ref coding_potential=<dir containing hexamer.tsv and
logit.RData>`), which runs `cpat -g target.fa -d logit.RData -x
hexamer.tsv -o prefix` and then restores original sequence IDs (CPAT
uppercases every ID it writes, and mishandles some pipe-delimited
headers -- see `repair_cpat_ids.py`'s module docstring). **Unlike
codonW, CPAT has no working-directory requirement** -- confirmed by
direct test, running it from an unrelated cwd with absolute paths for
every argument, and its generated internal R script embedded those
absolute paths correctly rather than assuming anything about cwd. As
long as every path handed to `cpat` is absolute (which `run_features.py`
now guarantees for `--module-ref`), there is nothing to isolate here the
way codonW's bare-relative-filename `.coa` convention requires.

Real output columns verified against real data: `ID, mRNA, ORF_strand,
ORF_frame, ORF_start, ORF_end, ORF, Fickett, Hexamer, Coding_prob`, plus
`Gene_ID` (renamed from CPAT's own `seq_ID` after ID restoration).
