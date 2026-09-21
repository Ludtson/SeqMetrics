# Install: coding_potential (CPAT)

CPAT is open-source and pip-installable, no registration needed. Verified
real: `CPAT` exists on PyPI (checked directly against
`pypi.org/pypi/CPAT/json`, not just assumed).

```
pip install CPAT
which cpat make_hexamer_tab make_logitModel   # confirm all three resolve
```

Also needs Rscript (R) on PATH -- CPAT's scoring step generates and runs
a small R script internally (`predict()` against the trained logistic
model), so `Rscript` must be installed even though you never invoke it
directly.

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
**Real, currently-unresolved gap**: building the `noncoding.fa` input
needs a `get_noncoding_data.py` script that the source project's own
`cpat_orchestrator.sh` depends on (`command -v get_noncoding_data.py` is
checked before it runs) -- but that script does not exist anywhere in
this project's files, only the reference to it. Confirmed via a full
search, not an oversight in copying. Until it's located or rewritten,
Stage 1 can only be run for a species where you already have a real
noncoding training set some other way.

Three species already have trained models from the source project's
prior run (Arabidopsis thaliana, Brassica rapa, Oryza sativa) --
`2026fall-mla-chapter/feature_outputs/cpat_output/<species>/*_cpat_output/{1_hexamer,2_model}/`.
Confirmed by direct test: SeqMetrics' scoring output using Athaliana's
existing model reproduces the historical run's `Coding_prob`/`Fickett`/
`Hexamer` values byte-for-byte across 25,278 overlapping genes.

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
