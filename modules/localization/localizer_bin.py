#!/usr/bin/env python3
"""Convert LOCALIZER tabular output to a per-protein binary table.

Expected LOCALIZER input: tab-separated text with at least 4 columns:
ID, chloroplast signal, mitochondrial signal, NLS signal.

Any field that starts with 'Y' (case-insensitive) is interpreted as present.
"""

import argparse
import csv
import re
import sys
from pathlib import Path


def _is_present(value: str) -> int:
    return 1 if value.strip().upper().startswith("Y") else 0


def _extract_prob(value: str) -> str:
    # Examples: "Y (0.877 | 62-130)" or "-"
    m = re.search(r"\((\d*\.\d+|\d+)\s*\|", value)
    return m.group(1) if m else "NA"


def _extract_nls_motif(value: str) -> str:
    # Example: "Y (KRKR)"; LOCALIZER NLS does not expose probabilities.
    m = re.search(r"Y\s*\(([^)]+)\)", value.strip(), flags=re.IGNORECASE)
    return m.group(1).strip() if m else "NA"


def _parse_data_row(line: str):
    stripped = line.strip()
    if not stripped:
        return None
    if stripped.startswith("#"):
        return None
    if set(stripped) == {"-"}:
        return None

    # LOCALIZER outputs are usually tabular with either tabs or aligned spaces.
    fields = [x.strip() for x in re.split(r"\t+|\s{2,}", stripped) if x.strip()]
    if len(fields) < 4:
        return None

    # Skip header-like rows.
    if fields[0].lower() in {"identifier", "id"}:
        return None

    return fields[0], fields[1], fields[2], fields[3]


def parse_localizer(input_path: Path, output_path: Path):
    stats = {
        "CLS": 0,
        "MLS": 0,
        "NLS": 0,
        "C_M": 0,
        "C_N": 0,
        "M_N": 0,
        "C_M_N": 0,
    }

    with input_path.open("r", newline="") as inf, output_path.open("w", newline="") as outf:
        writer = csv.writer(outf, delimiter="\t")

        # Output is standardized for downstream merges.
        writer.writerow(
            [
                "ID",
                "CLS_binary_0no_1yes",
                "MLS_binary_0no_1yes",
                "NLS_binary_0no_1yes",
                "CLS_prob",
                "MLS_prob",
                "NLS_motif",
            ]
        )

        found_rows = False
        for line in inf:
            parsed = _parse_data_row(line)
            if parsed is None:
                continue

            found_rows = True
            gene_id, c_field, m_field, n_field = parsed

            cls = _is_present(c_field)
            mls = _is_present(m_field)
            nls = _is_present(n_field)
            cls_prob = _extract_prob(c_field)
            mls_prob = _extract_prob(m_field)
            nls_motif = _extract_nls_motif(n_field)

            if cls:
                stats["CLS"] += 1
            if mls:
                stats["MLS"] += 1
            if nls:
                stats["NLS"] += 1
            if cls and mls:
                stats["C_M"] += 1
            if cls and nls:
                stats["C_N"] += 1
            if mls and nls:
                stats["M_N"] += 1
            if cls and mls and nls:
                stats["C_M_N"] += 1

            writer.writerow([gene_id, cls, mls, nls, cls_prob, mls_prob, nls_motif])

        if not found_rows:
            raise ValueError("no LOCALIZER prediction rows were found in input")

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Convert LOCALIZER output to binary TSV and report summary stats."
    )
    parser.add_argument("-i", "--input", required=True, help="LOCALIZER input .txt/.tsv file")
    parser.add_argument("-o", "--output", required=True, help="Output binary TSV")

    args = parser.parse_args()

    input_path = Path(args.input).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()

    if not input_path.exists():
        print(f"ERROR: input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        stats = parse_localizer(input_path, output_path)
    except Exception as exc:  # pragma: no cover
        print(f"ERROR: failed to parse {input_path}: {exc}", file=sys.stderr)
        sys.exit(1)

    print("LOCALIZER binary conversion complete", file=sys.stderr)
    print(f"Input:  {input_path}", file=sys.stderr)
    print(f"Output: {output_path}", file=sys.stderr)
    print(
        "Counts: "
        f"CLS={stats['CLS']} MLS={stats['MLS']} NLS={stats['NLS']} "
        f"C_M={stats['C_M']} C_N={stats['C_N']} M_N={stats['M_N']} C_M_N={stats['C_M_N']}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
