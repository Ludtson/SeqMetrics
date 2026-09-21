#!/usr/bin/env python3
"""Run LOCALIZER and optionally convert predictions to binary TSV.

This script is designed for single-environment workflows where LOCALIZER,
pepstats (EMBOSS), and Perl are installed together.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def check_dependencies(localizer_script: Path):
    report = {
        "localizer_script": "OK" if localizer_script.exists() else "MISSING",
        "pepstats": "OK" if shutil.which("pepstats") else "MISSING",
        "perl": "OK" if shutil.which("perl") else "MISSING",
    }
    return report


def print_report(report):
    for key in sorted(report):
        print(f"CHECK {key}={report[key]}")


def log(msg: str, quiet: bool):
    if not quiet:
        print(msg)


def main():
    parser = argparse.ArgumentParser(description="Run LOCALIZER and convert output")
    parser.add_argument("fasta", help="Input FASTA")
    parser.add_argument(
        "--localizer-script",
        required=True,
        help="Path to LOCALIZER.py",
    )
    parser.add_argument(
        "--mode",
        choices=["plant", "effector"],
        default="plant",
        help="LOCALIZER mode",
    )
    parser.add_argument(
        "--output",
        default="results_localizer",
        help="Output directory",
    )
    parser.add_argument(
        "--localizer-parser",
        default="bin/localizer_bin.py",
        help="Path to localizer binary parser script",
    )
    parser.add_argument(
        "--no-binary",
        action="store_true",
        help="Do not generate binary TSV",
    )
    parser.add_argument(
        "--mature-sequences",
        action="store_true",
        help="Effector mode: pass -M to LOCALIZER",
    )
    parser.add_argument(
        "--signal-peptide-length",
        type=int,
        default=None,
        help="Effector mode: pass -S <x> to LOCALIZER",
    )
    parser.add_argument("--check", action="store_true", help="Dependency check only")
    parser.add_argument("--dry-run", action="store_true", help="Print commands only")
    parser.add_argument("--quiet", action="store_true", help="Minimal console output")

    args = parser.parse_args()

    fasta = Path(args.fasta).expanduser().resolve()
    if not fasta.exists():
        sys.exit(f"ERROR: FASTA not found: {fasta}")

    localizer_script = Path(args.localizer_script).expanduser().resolve()
    report = check_dependencies(localizer_script)
    if not args.quiet:
        print_report(report)

    missing = [k for k, v in report.items() if v != "OK"]
    if args.check:
        return
    if missing:
        sys.exit("ERROR: missing dependencies: " + ", ".join(missing))

    outdir = Path(args.output).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    prefix = fasta.stem
    raw_out = outdir / f"{prefix}.localizer.raw.txt"
    stderr_log = outdir / f"{prefix}.localizer.stderr.log"
    binary_out = outdir / f"{prefix}.localizer_binary.tsv"

    cmd = [sys.executable, str(localizer_script)]
    cmd.extend(["-p" if args.mode == "plant" else "-e", "-i", str(fasta)])

    if args.mode == "effector":
        if args.mature_sequences:
            cmd.append("-M")
        if args.signal_peptide_length is not None:
            cmd.extend(["-S", str(args.signal_peptide_length)])

    log("STEP LOCALIZER predict", args.quiet)
    log("CMD: " + " ".join(cmd), args.quiet)
    if not args.dry_run:
        with raw_out.open("w") as out, stderr_log.open("w") as err:
            subprocess.run(cmd, check=True, stdout=out, stderr=err)

    if not args.no_binary:
        parser_script = Path(args.localizer_parser).expanduser().resolve()
        conv_cmd = [
            sys.executable,
            str(parser_script),
            "-i",
            str(raw_out),
            "-o",
            str(binary_out),
        ]
        log("STEP LOCALIZER parse", args.quiet)
        log("CMD: " + " ".join(conv_cmd), args.quiet)
        if not args.dry_run:
            subprocess.run(conv_cmd, check=True)

    print(f"DONE raw: {raw_out}")
    print(f"DONE stderr: {stderr_log}")
    if not args.no_binary:
        print(f"DONE binary: {binary_out}")


if __name__ == "__main__":
    main()
