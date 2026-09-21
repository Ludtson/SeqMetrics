# Install: codon_usage (codonW)

codonW is open-source and free (originally by John Peden) -- freely
installable, no registration needed. Verified real: `codonw` exists on
bioconda (checked directly against `anaconda.org/bioconda/codonw`, not
just assumed from the conda-forge/bioconda channel convention).

```
conda create -n codon_w -c bioconda -c conda-forge codonw
conda activate codon_w
which codonw   # confirm it resolves
```

This lab's existing WSL setup already has this working under the conda
env name `codon_w`, confirmed via `which codonw` on 2026-09-20.

## This is not a one-shot "run the tool" module

Unlike composition/aggregation/disorder, codonW cannot score a sequence
against nothing -- it needs a **per-species reference** built first, in a
genuinely separate stage, from real data this repo does not ship or
generate on its own. There are two stages, and they are not the same
operation:

**Stage 1 -- reference-building (calibration), once per species.**
Builds `cai.coa`/`fop.coa`/`cbi.coa` from that species' own codon usage.
Two ways to seed it, in order of preference, with an automatic fallback
so nobody is ever fully blocked:

1. **HEG (biological), via keyword annotation** -- `extract_heg_ids.py
   --fasta-dir DIR --gff-dir DIR --out-dir DIR` matches a fixed,
   plant-appropriate keyword set (ribosomal proteins, elongation
   factors, Rubisco, chlorophyll a/b, GAPDH/enolase, HSP70/chaperones --
   the standard universal HEG markers) against a GFF3's functional
   annotations, propagated gene<->mRNA<->CDS. Works for **any** species
   as-is -- **you need that species' own CDS FASTA + matching GFF3**,
   which this repo does not provide. Preferred method.
2. **HEG (biological), via conserved-domain HMM search** --
   `extract_heg_ids_hmm.py`, matches a protein FASTA against a small
   curated Pfam-A profile subset (`heg_pfam_accessions.tsv` in this
   module -- ~30 accessions across ribosomal proteins, elongation
   factors, GAPDH, enolase, and HSP70/chaperonins, each accession
   verified live against the InterPro API, not recalled from memory).
   Annotation-independent -- works even when a species' GFF3 is sparse
   or full of "hypothetical protein" placeholders, since it finds the
   conserved domain itself rather than trusting someone else's label for
   it. **Mostly, not entirely, clade-agnostic**: the ribosomal/EF/GAPDH
   markers are near-universal across bacteria/archaea/eukaryotes, but
   two categories (Rubisco, chlorophyll a/b-binding protein) are
   photosynthesis-specific -- pass `--exclude-photosynthesis` for a
   non-plant genome, where those two can only ever produce false
   negatives, never real hits.

   Pfam is CC0-licensed (public domain -- verified directly against
   Pfam's own current documentation, not assumed), so redistributing
   even a full local copy would be fine, but there's no reason to ship
   the ~1.5GB+ full database when ~30 profiles are all this needs. Build
   the small subset yourself, once (not per-species):
   ```
   wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
   gunzip Pfam-A.hmm.gz
   hmmpress Pfam-A.hmm
   cut -f1 heg_pfam_accessions.tsv | tail -n +2 > heg_accessions.txt
   hmmfetch -f Pfam-A.hmm heg_accessions.txt > heg_markers.hmm
   ```
   Then per species:
   ```
   python extract_heg_ids_hmm.py --protein-fasta species.faa \
       --hmm-profile heg_markers.hmm --species Athaliana --out-dir OUT
   ```
   Needs `hmmsearch`/`hmmfetch`/`hmmpress` (HMMER, bioconda: `conda
   create -n hmmer -c bioconda hmmer`) and a protein FASTA whose headers
   match the species' CDS FASTA headers exactly (the standard case when
   both come from the same genome annotation).

   (An orthology-projection method via OrthoFinder's `N0.tsv` was
   considered and dropped -- newer OrthoFinder versions don't reliably
   produce `N0.tsv` at all, since it depends on successfully rooting a
   species tree during the hierarchical-orthogroup step. A
   BLASTP-to-single-reference-species shortcut was also considered and
   dropped -- accurate enough for one project's own species set, but not
   general enough to ship as SeqMetrics' answer for arbitrary future
   users on arbitrary taxa, since a single reference species' HEG set
   doesn't transfer across large evolutionary distances or non-plant
   marker categories.)
3. **Statistical top-5%-Fop (the automatic fallback, "coa" mode)** --
   needs nothing but the species' own CDS FASTA. Runs `codonw -coa_cu`
   once on the whole gene set, ranks every gene by the Fop value that
   first pass assigns it, and rebuilds the reference from just the
   top 5%. This is **self-referential by construction** -- it defines
   "high expression" as "whatever this species' own codon bias already
   looks most biased toward," not real expression data. It tracks real
   biology reasonably well in a genome with strong translational
   selection, and can drift from anything real in one that doesn't --
   there's no external check. Prefer HEG mode when you can.

Run via `calculate_indices.sh <cleaned_fasta_abs> <species> <workdir_abs>
<heg_dir_or_NONE> <basis_mode: auto|heg|coa>` for one species, or let
`orchestrate_codonw.sh -i <dir> -o <outdir> -e <heg_dir> -b <basis_mode>`
drive it across every species-named file in a directory (species is
parsed from each file's basename). `auto` tries HEG first and falls back
to top-5%-Fop only if no HEG file is found for that species -- and now
logs which one actually ran per species, rather than deciding it by
accident of which file happened to exist.

**Stage 2 -- scoring, against an already-built reference.** This is
SeqMetrics' own job (`run_features.py --modules codon_usage`) --
`obtain_indices.sh` runs codonW with `-cai_file cai.coa -fop_file
fop.coa -cbi_file cbi.coa` against whatever target sequences you give
it. Two things that will silently break this if missed:
- **The `.coa` files must be in codonW's working directory** -- codonW
  reads them as bare relative filenames, not by path. `obtain_indices.sh`
  already `cd`s into the workdir before invoking it; a custom invocation
  that doesn't will produce codonW's own generic built-in tables
  instead, with no error.
- **No silent fallback, by design** -- `obtain_indices.sh` hard-fails if
  `cai.coa` isn't already present rather than letting codonW quietly use
  its built-in defaults, which were judged wrong for this kind of
  cross-species comparison. If Stage 2 fails with "cai.coa not found,"
  that's Stage 1 not having been run for that species -- not a bug to
  route around.

**Threading is per-species**, via `orchestrate_codonw.sh -t N` (`xargs
-P N`, one worker per species file), not per-sequence.

Real output columns not yet verified against a sample file in this repo
-- test on real data before trusting the exact column set (see main
README).
