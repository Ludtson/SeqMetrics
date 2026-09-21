#!/usr/bin/env python3
"""
compute_ctth.py

C-terminal tail hydrophobicity (CTTH), unbundled from the HTLCP pipeline
(E:\\features-rebirth.zip -> htlcp/src/hcaTango-pipeline/bin/run_htlc_pipeline.py,
compute_ctth(), lines 26-95 of that file).

Author:   Adekola Owoyemi, Casola Lab, TAMU
Version:  0.1.0 (pre-release, unshared)

--------------------------------------------------------------------------
WHY THIS IS ITS OWN MODULE
--------------------------------------------------------------------------
The original HTLCP pipeline bundles CTTH alongside HCA/TANGO, LOCALIZER,
and a pepstats orchestrator behind one shared virtual environment — a
~6GB zip when moved between machines, even though CTTH itself is pure
Python with zero external dependencies (no hcatk, no tango binary, no
Perl/LOCALIZER, no EMBOSS). Lifted out here as its own module, matching
the same per-tool `pipeline/<tool>_extracted/` convention already used
for hca_tango, pepstats, codonw, and iupred this project. The logic below
is unchanged from the original (reviewed, not rewritten) — only the
surrounding CLI/I-O is new, so this can run standalone on any real FASTA
without the rest of the HTLCP bundle.

--------------------------------------------------------------------------
METHOD (matches manuscripts/Adekola_et_al_2025-6.v0.0_cc.docx methods text:
"C-terminal tail hydrophobicity scores were obtained using the GRAVY
calculator... in the C-terminal 30 amino acids... along the
Kyte-Doolittle scale.")
--------------------------------------------------------------------------
CTTH_30       : mean Kyte-Doolittle hydropathy over the last 30 residues
                (seq[-30:]). This is the one that matches the manuscript's
                stated definition exactly.
CTTH_25       : mean Kyte-Doolittle hydropathy over seq[-30:-5] — NOT the
                last 25 residues (that would be seq[-25:]). This is the
                25 residues immediately BEFORE the terminal 5, i.e. the
                "near-tail" window with the very end of the protein
                deliberately excluded (the terminal 5 are analyzed
                separately, below, for charge instead of hydrophobicity).
                Worth being explicit about this in any methods write-up
                that cites "CTTH_25" — the name reads as "last 25 aa,"
                but the actual window is not that.
Tail5_charge  : net charge (count of K+R minus count of D+E) over the
                last 5 residues only.

Proteins shorter than the window length are handled gracefully (Python
slicing clamps rather than erroring); residues outside the 20 standard
amino acids are silently excluded from the hydrophobicity average rather
than treated as zero.
--------------------------------------------------------------------------
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
            out, fieldnames=["Seq_ID", "Length", "CTTH_30", "CTTH_25", "Tail5_charge"], delimiter="\t"
        )
        writer.writeheader()
        for sid in sorted(rows):
            row = {"Seq_ID": sid}
            row.update(rows[sid])
            writer.writerow(row)

    print(f"[ctth] Wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
