#!/usr/bin/env python3
"""
Repair CPAT output IDs by restoring original header case (and, for pipe-
delimited headers, the clean gene ID).

Renames the first column header to 'Gene_ID'.

Author:   Adekola Owoyemi, Casola Lab, TAMU
Version:  0.1.0

Usage:
  repair_cpat_ids.py \
      --fasta transcripts.fa \
      --cpat-table species_cpat.tsv \
      --out repaired_best.tsv

Why this is needed:
  CPAT uppercases every ID it writes. For headers with no whitespace before
  the first "|" (e.g. "ID|chromosome:..."), it also keeps the whole pipe
  chain as the ID. Some headers (e.g. CNCB GWH) start with an accession and
  carry the gene ID in an "OriID=" field.

How IDs are matched:
  The lookup key on both sides is the first whitespace token, cut at the
  first "|" and lowercased. The value written back is the canonical ID from
  the FASTA header, extracted by the rule chosen with --species.
"""

import argparse
import re
import sys


# Header ID rules (duplicated rather than imported so the script runs standalone).
def fasta_id_athaliana(header): return header[1:].split()[0].split("|")[0].strip()
def fasta_id_oryza_sativa(header): return header[1:].split("|")[0].strip()
def fasta_id_brassica_rapa(header):
    m = re.search(r"OriID=([^\s\t]+)", header)
    if not m:
        return header[1:].split()[0].strip()
    oid = m.group(1)
    oid = re.sub(r"\.mRNA_OriTrasc.*$", "", oid)
    oid = re.sub(r"\.mRNA$", "", oid)
    return oid.strip()


def fasta_id_generic(header): return header[1:].split()[0].strip()


CANONICAL_ID_FN = {
    "athaliana": fasta_id_athaliana,
    "brassica_rapa": fasta_id_brassica_rapa,
    "oryza_sativa": fasta_id_oryza_sativa,
    # Default: first whitespace token, case restored, no further cleanup.
    "generic": fasta_id_generic,
}


def parse_args():
    p = argparse.ArgumentParser(
        description="Replace CPAT table IDs with original FASTA header IDs."
    )
    p.add_argument("--fasta", required=True,
                   help="Original transcript FASTA used as CPAT -t input.")
    p.add_argument("--cpat-table", required=True,
                   help="CPAT best ORF / classification table (TSV).")
    p.add_argument("--out", required=True,
                   help="Output path for repaired table.")
    p.add_argument("--species", required=True, choices=sorted(CANONICAL_ID_FN),
                   help="Selects the per-species canonical-ID extraction rule.")
    p.add_argument("--id-column", type=int, default=0,
                   help="0-based index of the ID column (default: 0, CPAT's 'seq_ID').")
    return p.parse_args()


def cpat_lookup_key(token: str) -> str:
    """Normalize a CPAT ID to the map key: cut at the first '|', lowercase."""
    return token.split("|")[0].lower()


def load_fasta_ids(fasta_path, canonical_fn):
    """Return dict: cpat_lookup_key(raw_first_token) -> canonical gene ID.

    The key mimics CPAT's own ID for the header (first token, cut at '|',
    lowercased); the value is the canonical gene ID from the chosen rule.
    """
    mapping = {}
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if not line[1:].strip():
                    continue
                raw_token = line[1:].split()[0] if line[1:].split() else ""
                if not raw_token:
                    continue
                key = cpat_lookup_key(raw_token)
                canonical = canonical_fn(line)
                # Only keep first occurrence.
                mapping.setdefault(key, canonical)
    return mapping


def repair_table(cpat_table, out_path, id_col, id_map):
    fixed = 0
    missing = 0
    header_processed = False

    # LF only: a trailing CR would break exact-ID joins with other tables.
    with open(cpat_table) as inp, open(out_path, "w", newline="\n") as out:
        for line in inp:
            if not line.strip():
                continue

            if line.startswith("#"):
                out.write(line)
                continue

            parts = line.rstrip("\n").split("\t")

            if not header_processed:
                if id_col < len(parts):
                    parts[id_col] = "Gene_ID"
                out.write("\t".join(parts) + "\n")
                header_processed = True
                continue

            if id_col >= len(parts):
                out.write(line)
                continue

            cur_id = parts[id_col]
            # Same normalization applied when the map was built, so both
            # sides of the lookup are on equal footing.
            key = cpat_lookup_key(cur_id)

            if key in id_map:
                parts[id_col] = id_map[key]
                fixed += 1
            else:
                missing += 1

            out.write("\t".join(parts) + "\n")

    sys.stderr.write(f"[repair_cpat_ids] Fixed IDs: {fixed}, Missing matches: {missing}\n")


def main():
    args = parse_args()
    id_map = load_fasta_ids(args.fasta, CANONICAL_ID_FN[args.species])
    if not id_map:
        sys.stderr.write(f"[repair_cpat_ids] ERROR: No IDs parsed from FASTA {args.fasta}\n")
        sys.exit(1)

    repair_table(args.cpat_table, args.out, args.id_column, id_map)


if __name__ == "__main__":
    main()
