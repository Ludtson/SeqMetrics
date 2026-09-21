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
    # Bug fix (recurring): this project's B. rapa FASTA headers carry an
    # internal accession as the first whitespace token (e.g.
    # ">GWHTAAES000001  Protein=... OriID=BraA01g000010.3C.mRNA ...") --
    # the real, canonical gene ID used everywhere else in this pipeline
    # sits in an "OriID=" field further into the header, not in the first
    # token. Without this check, restored IDs come back as the internal
    # GWH accession instead of the real gene ID -- confirmed twice now:
    # once in the original full-scale B. rapa CodonW run (fixed inline at
    # the time, not saved back here), and again when recovering B. rapa's
    # missing genes, which silently repeated the exact same bug because
    # the inline fix never made it into this script. Same OriID=
    # extraction already used in pepstats_v2_extracted/remap_pepstats_ids.py's
    # fasta_id_brassica_rapa() -- kept consistent with that logic here so
    # the two don't drift apart again.
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
    # Bug fix: some headers (confirmed for this project's O. sativa CDS FASTA,
    # e.g. "OSJAP01G00010.5|chromosome:IRGSP-1.0:1:1:...|oge_maker_genes|...")
    # pack multiple pipe-delimited fields immediately after the ID with no
    # whitespace before the first "|". Splitting on whitespace alone captured
    # the entire pipe-chain as a single "token" instead of just the gene ID.
    # Confirmed by direct inspection of a real production run's output before
    # this fix: restored IDs were the full pipe-delimited string, not the
    # clean gene ID. Splitting on "|" here is a no-op for headers that don't
    # have this pattern (A. thaliana), since their first whitespace-token
    # never itself contains a "|".
    token = token.split("|")[0]
    return token


def command_prepare(input_fasta: Path, output_fasta: Path, output_map: Path) -> None:
    records = parse_fasta_records(input_fasta)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    output_map.parent.mkdir(parents=True, exist_ok=True)

    # newline="\n" on both: output_fasta is fed straight to codonw as real
    # sequence input (a stray CRLF from Windows Python text mode corrupts
    # sequence lines, not just an ID lookup), and output_map's IDs get
    # matched by exact string elsewhere -- same failure class as the
    # confirmed CRLF/HEG-ID bug already found and fixed in this module.
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

        # Bug fix: previously used stripped.split(None, 1), which finds the
        # first token correctly but discards what the ORIGINAL delimiter
        # after it was (tab, space, multiple spaces...). Reconstructing with
        # a hard-coded " " silently turned every tab-delimited data row into
        # a row where the ID and the first real column are space-joined
        # instead of tab-separated — every column after that stays correctly
        # tab-separated, so the corruption is invisible until something
        # parses the file by splitting on tabs and gets misaligned columns.
        # The header row was never affected (its first token is "title",
        # which never matches a seqN key, so it's never rewritten) — which is
        # exactly why this was easy to miss on a quick look.
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

    # newline="\n" for consistency with the other writes in this module --
    # this is the final restored results table, so the risk here is lower
    # (whatever reads it next is more likely Python/pandas, which tolerates
    # CRLF), but there's no reason to leave one unprotected write in a module
    # already confirmed to have this exact failure class.
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