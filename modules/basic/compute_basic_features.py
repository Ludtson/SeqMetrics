#!/usr/bin/env python3
"""
Basic nucleotide features -- length, GC% -- computed directly from the CDS
sequence, no external tool, no genomic/GFF context needed.

Originally also computed protein length + 20 amino-acid composition
percentages, but that duplicated PEPSTATS's own Residues and
{Ala,Cys,...}_MolePct columns exactly (same 20 percentages, different
naming convention) -- caught 2026-09-20, dropped rather than kept as a
redundant cross-check nobody asked for. Nucleotide composition is the one
thing no other module in this pipeline touches at all, which is the only
reason this module exists on its own.
"""
import argparse
import csv
import sys
from pathlib import Path


def fasta_iter(path):
    name, buf = None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line)
        if name is not None:
            yield name, "".join(buf)


def nt_features(seq):
    seq = seq.upper()
    length = len(seq)
    gc = sum(1 for c in seq if c in "GC")
    return {"Length_nt": length, "GC_pct": round(100 * gc / length, 3) if length else "NA"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("nt", type=Path, help="Nucleotide (CDS) FASTA")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    rows = {sid: nt_features(seq) for sid, seq in fasta_iter(args.nt)}

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=["Seq_ID", "Length_nt", "GC_pct"], delimiter="\t")
        w.writeheader()
        for sid in sorted(rows):
            row = {"Seq_ID": sid}
            row.update(rows[sid])
            w.writerow(row)

    print(f"[basic] {len(rows)} sequences -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
