#!/usr/bin/env python3
"""
compute_ctth.py

C-terminal tail hydrophobicity (CTTH) features from a protein FASTA.
Pure Python, no external dependencies.

Author:   Adekola Owoyemi, Casola Lab, TAMU
Version:  0.1.0

Features (Kyte-Doolittle hydropathy):
  CTTH_30       mean hydropathy of the last 30 residues (seq[-30:]).
  CTTH_25       mean hydropathy of seq[-30:-5]: the 25 residues before the
                terminal 5, not the last 25 residues.
  Tail5_charge  net charge of the last 5 residues: (K + R) - (D + E).

Proteins shorter than a window use what is available. Non-standard residues
are left out of the hydropathy mean.
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, Iterable, Tuple

KD = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6,
    "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}


def fasta_iter(path: Path) -> Iterable[Tuple[str, str]]:
    name = None
    seq = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.upper())
        if name is not None:
            yield name, "".join(seq)


def hydrophobicity(seq: str):
    if not seq:
        return "NA"
    vals = [KD[aa] for aa in seq if aa in KD]
    if not vals:
        return "NA"
    return round(sum(vals) / len(vals), 3)


def tail_charge(seq: str) -> int:
    pos = seq.count("K") + seq.count("R")
    neg = seq.count("D") + seq.count("E")
    return pos - neg


def compute_ctth(fasta: Path) -> Dict[str, Dict[str, str]]:
    out = {}
    for sid, seq in fasta_iter(fasta):
        last30 = seq[-30:]
        last30_minus5 = seq[-30:-5]
        out[sid] = {
            "Length": str(len(seq)),
            "CTTH_30": str(hydrophobicity(last30)),
            "CTTH_25": str(hydrophobicity(last30_minus5)),
            "Tail5_charge": str(tail_charge(seq[-5:])),
        }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", type=Path, help="Protein FASTA file")
    ap.add_argument("--out", required=True, type=Path, help="Output TSV path")
    args = ap.parse_args()

    print(f"[ctth] Reading {args.fasta}", file=sys.stderr)
    rows = compute_ctth(args.fasta)
    print(f"[ctth] Computed CTTH for {len(rows)} sequences", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as out:
        writer = csv.DictWriter(
            out, fieldnames=["Seq_ID", "Length", "CTTH_30", "CTTH_25", "Tail5_charge"], delimiter="\t",
            lineterminator="\n",  # csv defaults to \r\n on every OS
        )
        writer.writeheader()
        for sid in sorted(rows):
            row = {"Seq_ID": sid}
            row.update(rows[sid])
            writer.writerow(row)

    print(f"[ctth] Wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
