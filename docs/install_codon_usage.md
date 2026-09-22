# Install: codon_usage (codonW)

codonW is open-source and free, no registration needed.

```
conda create -n codon_w -c bioconda -c conda-forge codonw
conda activate codon_w
which codonw
```

## Two stages

codonW cannot score a sequence without a per-species reference built
first, in a separate stage, from data this repo does not ship.

### Stage 1 — reference-building (calibration), once per species

Builds `cai.coa`/`fop.coa`/`cbi.coa` from that species' own codon
usage. Three ways to seed it, in order of preference:

**1. HEG via GFF keyword annotation (preferred):**
```
extract_heg_ids.py --fasta-dir DIR --gff-dir DIR --out-dir DIR
```
Matches a fixed keyword set (ribosomal proteins, elongation factors,
Rubisco, chlorophyll a/b, GAPDH/enolase, HSP70/chaperones) against a
GFF3's functional annotations. Needs that species' own CDS FASTA and
matching GFF3.

**2. HEG via conserved-domain HMM search**, when GFF3 annotation is
sparse or missing usable keywords:

Needs HMMER (a separate dependency from codonW):
```
conda create -n hmmer -c bioconda hmmer
conda activate hmmer
which hmmsearch hmmfetch hmmpress
```

Build the marker set once (not per species) — Pfam-A must be
downloaded in full to extract the ~30 needed profiles (`hmmfetch` only
works against an already-downloaded, indexed database); delete the
full database afterward, keep only the small derived files:
```
cd /tmp/pfam_scratch
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
hmmfetch --index Pfam-A.hmm
cut -f1 heg_pfam_accessions.tsv | tail -n +2 > heg_accessions.txt
hmmfetch -f Pfam-A.hmm heg_accessions.txt > heg_markers.hmm
```
`hmmfetch` requires its own index (`hmmfetch --index`), separate from
`hmmpress`'s index — both are needed. Accessions in
`heg_pfam_accessions.tsv` are unversioned (e.g. `PF00163`); the live
database's keys are versioned (e.g. `PF00163.25`) and change between
Pfam releases, so resolve the current version per accession at fetch
time rather than hardcoding one.

Per species:
```
python extract_heg_ids_hmm.py --protein-fasta species.faa \
    --hmm-profile heg_markers.hmm --species Athaliana --out-dir OUT
```
Needs a protein FASTA whose headers match the species' CDS FASTA
headers exactly.

Not fully clade-agnostic: Rubisco and chlorophyll a/b-binding protein
are photosynthesis-specific. Pass `--exclude-photosynthesis` for a
non-plant genome.

**3. Statistical top-5%-Fop (automatic fallback, `coa` mode)** — needs
only the species' CDS FASTA. Runs `codonw -coa_cu` once, ranks genes
by the resulting Fop value, and rebuilds the reference from the top
5%. Self-referential (defines "high expression" from the species' own
codon bias, not real expression data) — prefer HEG mode when
available.

### Running Stage 1

```
calculate_indices.sh <cleaned_fasta_abs> <species> <workdir_abs> <heg_dir_or_NONE> <basis_mode: auto|heg|coa>
```

`auto` tries HEG first, falls back to top-5%-Fop only if no HEG file
is found, and logs which one ran.

`orchestrate_codonw.sh` is meant to drive `calculate_indices.sh` across
every species in a directory, but calls a `run_species.sh` that does
not exist in this repo — it will fail immediately. Call
`calculate_indices.sh` directly per species until that script exists.

**`<cleaned_fasta_abs>` requirement**: `calculate_indices.sh` matches
HEG IDs against FASTA headers by exact string comparison, not a
first-token split. A CDS FASTA straight from a genome annotation
(multi-field headers, e.g. `>ATCG00500.1 pacid=... locus=...`) will
match zero HEG genes against a bare-gene-ID HEG list. Use a CDS FASTA
with single-token headers (one sequence per gene) that matches the
protein FASTA the HEG list was built from 1:1 — verify with `grep -c
"^>"` on both before trusting the pairing. Alternatively, run
`remap_seq_ids.py prepare` on the raw FASTA and translate the HEG ID
list into the resulting `seqN` ID space before calling
`calculate_indices.sh`.

If Stage 1 fails to build a reference, `calculate_indices.sh` exits
non-zero with an error naming the cause.

### Stage 2 — scoring, against an already-built reference

This is SeqMetrics' own job:
```
run_features.py --modules codon_usage --module-ref codon_usage=<dir containing cai.coa>
```
`obtain_indices.sh` runs codonW with `-cai_file cai.coa -fop_file
fop.coa -cbi_file cbi.coa`. The `.coa` files must be in codonW's
working directory — codonW reads them as bare relative filenames, not
by path; `obtain_indices.sh` handles this by `cd`-ing into the workdir
first. Fails with a clear error if `cai.coa` is missing rather than
falling back to codonW's built-in generic tables.

Threading is per-species, via `orchestrate_codonw.sh -t N`, not
per-sequence.
