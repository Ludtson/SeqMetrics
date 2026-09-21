#!/usr/bin/env python3
"""
run_dtm2.py

Wrapper around DeepTMHMM2 (fteufel/DeepTMHMM2 -- an ungated, pip-installable
reimplementation of DeepTMHMM by one of the original paper's authors; the
official DTU/BioLib DeepTMHMM requires registration for any local use, see
docs/install_tm_domain.md). Runs dtm2 into a scratch subdirectory and
collapses its four output files (TMRs.gff3, membrane_types.tsv,
predicted_topologies.3line, predictions.json) into one flat per-sequence
TSV, matching the single-table-per-module convention every other module
in this repo follows.

Usage:
  python run_dtm2.py -i proteins.faa -o output.tsv --dtm2 /path/to/dtm2
"""

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from pathlib import Path


def parse_gff3(gff3_path):
    """Returns {seq_id: n_tmrs} from the "## <id> Number of predicted TMRs: N"
    comment lines -- these are DeepTMHMM2's own count, not re-derived by
    counting "TMhelix" rows ourselves, so it can't silently drift from
    whatever counting rule the tool itself uses."""
    n_tmrs = {}
    pattern = re.compile(r"^## (\S+) Number of predicted TMRs: (\d+)")
    with open(gff3_path) as f:
        for line in f:
            m = pattern.match(line)
            if m:
                n_tmrs[m.group(1)] = int(m.group(2))
    return n_tmrs


def parse_topology_residue_counts(gff3_path):
    """Returns {seq_id: {topology_label: residue_count}} from the per-region
    data rows (id, source, topology_label, start, end, ...) -- e.g. how many
    residues are "inside" vs "outside" vs "TMhelix", not just the TMR count."""
    counts = {}
    with open(gff3_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip() or line.strip() == "//":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            seq_id, _, topology, start, end = fields[0], fields[1], fields[2], fields[3], fields[4]
            try:
                length = int(end) - int(start) + 1
            except ValueError:
                continue
            counts.setdefault(seq_id, {})[topology] = counts.setdefault(seq_id, {}).get(topology, 0) + length
    return counts


def parse_membrane_types(tsv_path):
    """Returns {seq_id: {membrane_type_column: probability}}."""
    result = {}
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            seq_id = row["id"]
            result[seq_id] = {k: v for k, v in row.items() if k != "id"}
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, type=Path, help="Protein FASTA")
    ap.add_argument("-o", "--out", required=True, type=Path, help="Output TSV path")
    ap.add_argument("--dtm2", required=True, help="Path to the dtm2 executable")
    ap.add_argument("--device", default="cpu", help="Device to run on (default: cpu -- "
                    "this wrapper exists specifically because the official DeepTMHMM "
                    "requires registration for CPU/local use; forcing GPU defeats that purpose "
                    "for anyone without one)")
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        cmd = [args.dtm2, str(args.input), tmp, "--device", args.device]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            sys.stderr.write(result.stderr)
            sys.exit(f"ERROR: dtm2 failed (exit {result.returncode})")

        gff3 = Path(tmp) / "TMRs.gff3"
        membrane_tsv = Path(tmp) / "membrane_types.tsv"
        if not gff3.exists():
            sys.exit(f"ERROR: dtm2 produced no TMRs.gff3 in {tmp}")

        n_tmrs = parse_gff3(gff3)
        topology_counts = parse_topology_residue_counts(gff3)
        membrane_probs = parse_membrane_types(membrane_tsv) if membrane_tsv.exists() else {}

        all_topology_labels = sorted({label for v in topology_counts.values() for label in v})
        all_membrane_cols = sorted({col for v in membrane_probs.values() for col in v})

        args.out.parent.mkdir(parents=True, exist_ok=True)
        # newline="\n": same CRLF-safety discipline as every other module in
        # this repo -- this table's Seq_ID gets joined by exact match against
        # other feature tables downstream.
        with open(args.out, "w", newline="\n") as f:
            writer = csv.writer(f, delimiter="\t")
            header = (["Seq_ID", "N_TMRs"]
                       + [f"Residues_{label}" for label in all_topology_labels]
                       + [f"MembraneProb_{col.replace(' ', '_')}" for col in all_membrane_cols])
            writer.writerow(header)
            for seq_id in sorted(n_tmrs):
                row = [seq_id, n_tmrs[seq_id]]
                row += [topology_counts.get(seq_id, {}).get(label, 0) for label in all_topology_labels]
                row += [membrane_probs.get(seq_id, {}).get(col, "NA") for col in all_membrane_cols]
                writer.writerow(row)


if __name__ == "__main__":
    main()
