#!/usr/bin/env python3
"""Create CodonW-safe FASTA IDs and restore original IDs in text outputs.

Commands:
  prepare: Build a 2-column mapping table (seq, id) and a renamed FASTA.
  restore: Replace first-column short IDs in output files with original IDs.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple


def parse_fasta_records(fasta_path: Path) -> List[Tuple[str, List[str]]]:
    records: List[Tuple[str, List[str]]] = []
    current_header: str | None = None
    current_seq_lines: List[str] = []

    with fasta_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if current_header is not None:
                    records.append((current_header, current_seq_lines))
                current_header = line[1:].strip()
                current_seq_lines = []
            else:
                if current_header is None:
                    if line.strip() == "":
                        continue
                    raise ValueError(
                        f"Invalid FASTA format in {fasta_path}: sequence line before first header."
                    )
                current_seq_lines.append(line)

    if current_header is not None:
        records.append((current_header, current_seq_lines))

    if not records:
        raise ValueError(f"No FASTA records found in {fasta_path}.")

    return records


_ORIID_RE = re.compile(r"OriID=([^\s\t]+)")


def original_id_from_header(header: str) -> str:
    # Some headers (e.g. CNCB GWH) start with an accession and carry the gene ID
    # in an "OriID=" field; use that when present.
    m = _ORIID_RE.search(header)
    if m:
        oid = m.group(1)
        oid = re.sub(r"\.mRNA_OriTrasc.*$", "", oid)
        oid = re.sub(r"\.mRNA$", "", oid)
        return oid.strip()

    # Keep the sequence identifier token (up to first whitespace) as the ID key.
    token = header.split()[0] if header.split() else ""
    if not token:
        raise ValueError("Encountered empty FASTA header identifier.")
    # Headers like "ID|chromosome:...|..." have no whitespace before the first
    # "|"; keep only the part before it.
    token = token.split("|")[0]
    return token


def command_prepare(input_fasta: Path, output_fasta: Path, output_map: Path) -> None:
    records = parse_fasta_records(input_fasta)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    output_map.parent.mkdir(parents=True, exist_ok=True)

    # LF only: CR characters would corrupt sequence lines and exact-ID matches.
    with output_fasta.open("w", encoding="utf-8", newline="\n") as fasta_out, output_map.open(
        "w", encoding="utf-8", newline="\n"
    ) as map_out:
        map_out.write("seq\tid\n")
        for idx, (header, seq_lines) in enumerate(records, start=1):
            short_id = f"seq{idx}"
            orig_id = original_id_from_header(header)
            map_out.write(f"{short_id}\t{orig_id}\n")
            fasta_out.write(f">{short_id}\n")
            for seq_line in seq_lines:
                fasta_out.write(f"{seq_line}\n")


def load_map(mapping_tsv: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with mapping_tsv.open("r", encoding="utf-8") as handle:
        header = handle.readline().strip()
        if header != "seq\tid":
            raise ValueError(
                f"Unexpected mapping header in {mapping_tsv}. Expected 'seq\\tid', got '{header}'."
            )
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                raise ValueError(f"Malformed mapping row in {mapping_tsv}: '{line}'")
            short_id, original_id = parts
            mapping[short_id] = original_id
    return mapping


_FIRST_TOKEN_RE = re.compile(r"^(\S+)(\s+)(.*)$")


def restore_file(mapping: Dict[str, str], file_path: Path) -> None:
    with file_path.open("r", encoding="utf-8") as handle:
        lines = handle.readlines()

    restored_lines: List[str] = []
    for line in lines:
        stripped = line.rstrip("\n")
        if stripped == "":
            restored_lines.append(line)
            continue

        # Keep the original delimiter after the first token (tab or spaces),
        # so tab-separated rows stay aligned.
        match = _FIRST_TOKEN_RE.match(stripped)
        if match:
            first, delimiter, rest = match.group(1), match.group(2), match.group(3)
        else:
            first, delimiter, rest = stripped, "", ""

        if first in mapping:
            replacement = mapping[first]
            if not rest:
                restored_lines.append(replacement + "\n")
            else:
                restored_lines.append(replacement + delimiter + rest + "\n")
        else:
            restored_lines.append(line)

    # LF only, as for the other outputs of this module.
    with file_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.writelines(restored_lines)


def command_restore(mapping_tsv: Path, files: List[Path]) -> None:
    mapping = load_map(mapping_tsv)
    for file_path in files:
        if not file_path.exists():
            continue
        restore_file(mapping, file_path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="CodonW sequence ID remapping helper")
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser(
        "prepare",
        help="Create a renamed FASTA and seq->original ID map",
    )
    prepare.add_argument("--input-fasta", required=True, type=Path)
    prepare.add_argument("--output-fasta", required=True, type=Path)
    prepare.add_argument("--output-map", required=True, type=Path)

    restore = subparsers.add_parser(
        "restore",
        help="Restore original IDs into output files using a seq->id map",
    )
    restore.add_argument("--mapping", required=True, type=Path)
    restore.add_argument("--files", nargs="+", required=True, type=Path)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "prepare":
        command_prepare(args.input_fasta, args.output_fasta, args.output_map)
    elif args.command == "restore":
        command_restore(args.mapping, args.files)


if __name__ == "__main__":
    main()