#!/usr/bin/env python3
"""
Repair CPAT output IDs by restoring original header case (and, for pipe-
delimited headers, the clean gene ID).

Renames the first column header to 'Gene_ID'.

Author:   Adekola Owoyemi, Casola Lab, TAMU
Version:  0.1.0 (pre-release, unshared)

Usage:
  repair_cpat_ids.py \
      --fasta transcripts.fa \
      --cpat-table species_cpat.tsv \
      --out repaired_best.tsv

--------------------------------------------------------------------------
CHANGELOG vs. original (~/scripts/repair_cpat_ids.py)
--------------------------------------------------------------------------
CPAT uppercases every ID it writes to its own output tables, and — for
headers with no whitespace before the first "|" (this project's O. sativa
convention, e.g. "OSJAP01G00010.5|chromosome:IRGSP-1.0:...") — it also
keeps the entire pipe-chain as the ID, not just the gene ID before the
first "|". Confirmed by direct test:

    input header:  OSJAP01G00010.5|chromosome:IRGSP-1.0:1:1:1000:1|...
    CPAT seq_ID:   OSJAP01G00010.5|CHROMOSOME:IRGSP-1.0:1:1:1000:1|...

A second, DIFFERENT species-specific problem exists for B. rapa: its CDS
FASTA header's first whitespace token is an internal accession
("GWHTAAES000001"), not the canonical gene ID — the real ID only appears
later in the header as "OriID=BraA01g000010.3C.mRNA". CPAT (like every
other tool in this pipeline that reads a FASTA header naively) grabs that
first token as its own seq_ID, so CPAT's raw output is keyed on
"GWHTAAES000001", not "BraA01g000010.3C". This is the exact same
B. rapa-specific problem already solved once for pepstats in
pepstats_v2_extracted/remap_pepstats_ids.py — rather than re-derive it
here, this script now reuses that exact same per-species canonical-ID
extraction (fasta_id_athaliana / fasta_id_brassica_rapa /
fasta_id_oryza_sativa), selected via --species, instead of the generic
"first token, then split on |" rule the original script effectively used
(and which is only actually correct for A. thaliana and O. sativa, not
B. rapa).

The lookup key on the CPAT-table side is still built by lowercasing (to
undo CPAT's uppercasing) and, for O. sativa specifically, also stripping
everything from the first "|" onward — mirroring the same fix already
applied in codonw_extracted/bin/remap_seq_ids.py and
hca_tango_extracted/bin/run_hca_tango.py. Both sides of the match
(building the map from the FASTA, and looking up CPAT's own ID) go
through the same per-species extraction where it matters, so case
restoration, pipe-chain cleanup, and the B. rapa OriID= rewrite can't
silently drift out of sync with each other or with the rest of the
pipeline.
--------------------------------------------------------------------------
"""

import argparse
import re
import sys


# --- Reused verbatim from pepstats_v2_extracted/remap_pepstats_ids.py ---
# Same per-species canonical-ID extraction already validated there. Kept
# duplicated (not imported) so this script stays runnable standalone from
# any directory, same convention as the rest of this project's scripts.
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
    # Correct default for any header with no pipe-chain and no B. rapa-style
    # OriID= convention -- e.g. this project's own composite
    # "gene_id::role::node" headers, which have neither. Restores case only;
    # does no species-specific cleanup, since there's none to do for a
    # pipe-free, OriID=-free header. (fasta_id_athaliana would actually
    # produce the identical result for such headers -- split("|")[0] on a
    # string with no "|" is a no-op -- but naming this explicitly avoids
    # implying an Athaliana-specific assumption for genes from any species.)
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
    """Normalize CPAT's own (uppercased, possibly pipe-chain-preserving)
    ID into the same shape used as the map key: lowercase, and cut at the
    first '|' the way O. sativa's raw header token would be. A no-op for
    A. thaliana/B. rapa IDs, which never contain '|'."""
    return token.split("|")[0].lower()


def load_fasta_ids(fasta_path, canonical_fn):
    """Return dict: cpat_lookup_key(raw_first_token) -> canonical gene ID.

    The map KEY mimics what CPAT itself would have produced from this
    header (first whitespace token, cut at '|', lowercased) so it can be
    matched directly against CPAT's table output. The map VALUE is the
    real, canonical, pipeline-wide gene ID for this species — which for
    B. rapa is NOT the same string as the raw first token.
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

    # newline="\n": this table's Gene_ID column is exactly the kind of thing
    # this project joins by exact match against other feature tables (see
    # the same fix applied in codon_usage/aggregation) -- a CRLF-corrupted
    # ID would silently fail such a join with no error anywhere.
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
