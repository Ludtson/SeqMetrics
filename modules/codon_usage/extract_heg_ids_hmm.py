#!/usr/bin/env python3
"""
extract_heg_ids_hmm.py

Identifies Highly Expressed Gene (HEG) seed IDs via conserved-domain HMM
search, as an alternative to extract_heg_ids.py's GFF3-keyword approach.
Annotation-independent: matches the protein sequence itself against a
small curated Pfam-A profile subset (see heg_pfam_accessions.tsv in this
directory), so it still finds real ribosomal/EF/GAPDH/etc. genes even in
a genome whose GFF3 is sparse or full of "hypothetical protein" entries.

Requires a protein FASTA whose headers are the SAME IDs as the species'
CDS FASTA (the standard convention when both come from the same genome
annotation) -- this script's output is meant to be handed straight to
codon_usage's calculate_indices.sh, which matches HEG IDs against CDS
FASTA headers by exact string, same as extract_heg_ids.py's output.

Setup (one-time, not per-species -- see docs/install_codon_usage.md):
  wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
  gunzip Pfam-A.hmm.gz
  hmmpress Pfam-A.hmm
  cut -f1 heg_pfam_accessions.tsv | tail -n +2 > heg_accessions.txt
  hmmfetch -f Pfam-A.hmm heg_accessions.txt > heg_markers.hmm

Usage:
  python extract_heg_ids_hmm.py --protein-fasta species.faa \
      --hmm-profile heg_markers.hmm --species Athaliana --out-dir OUT
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


def get_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--protein-fasta", required=True, type=Path,
                    help="Protein FASTA for one species (headers must match the CDS FASTA's headers)")
    p.add_argument("--hmm-profile", required=True, type=Path,
                    help="Curated HEG marker profile file (built via hmmfetch -- see module docstring)")
    p.add_argument("--species", required=True, help="Species name -- output is <species>_heg_ids.txt")
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--exclude-photosynthesis", action="store_true",
                    help="Skip the Rubisco/chlorophyll-a-b markers -- use for non-plant genomes "
                         "where those two categories can only ever produce false negatives, never hits")
    return p.parse_args()


def main():
    args = get_args()
    if not args.protein_fasta.exists():
        sys.exit(f"ERROR: protein FASTA not found: {args.protein_fasta}")
    if not args.hmm_profile.exists():
        sys.exit(f"ERROR: HMM profile not found: {args.hmm_profile} "
                 f"(see this script's module docstring for the one-time hmmfetch setup)")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.out_dir / f"{args.species}_heg_ids.txt"

    with tempfile.TemporaryDirectory() as tmp:
        tblout = Path(tmp) / "hits.tbl"
        # --cut_ga: use each Pfam family's own curator-set gathering threshold
        # rather than a single arbitrary e-value cutoff across very different
        # domain lengths/conservation levels -- this is Pfam's own documented
        # best practice for exactly this kind of family-membership call.
        cmd = ["hmmsearch", "--cut_ga", "--tblout", str(tblout),
               str(args.hmm_profile), str(args.protein_fasta)]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            sys.exit(f"ERROR: hmmsearch failed:\n{result.stderr}")

        excluded_accessions = set()
        if args.exclude_photosynthesis:
            acc_table = Path(__file__).parent / "heg_pfam_accessions.tsv"
            with open(acc_table) as f:
                next(f)
                for line in f:
                    fields = line.rstrip("\n").split("\t")
                    if fields[3].startswith("no"):
                        excluded_accessions.add(fields[0])

        hit_ids = set()
        with open(tblout) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.split()
                query_id, query_accession = fields[0], fields[3]
                if query_accession.split(".")[0] in excluded_accessions:
                    continue
                hit_ids.add(query_id)

    # newline="\n": force LF-only regardless of host OS -- calculate_indices.sh's
    # awk lookup (WSL/Linux) keys on an exact "ID" match against the CDS
    # FASTA header; a CRLF-written ID list (what plain Windows Python text
    # mode produces) keys on "ID\r" instead and silently matches nothing.
    # Same real bug confirmed and fixed in extract_heg_ids.py -- fixing here
    # too before it gets triggered by running this one under Windows Python.
    with open(out_path, "w", newline="\n") as f:
        for hit_id in sorted(hit_ids):
            f.write(hit_id + "\n")

    print(f"[{args.species}] {len(hit_ids)} HEG candidate(s) written to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
