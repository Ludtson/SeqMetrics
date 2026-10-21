#!/usr/bin/env python3
"""Run LOCALIZER and optionally convert predictions to binary TSV.

This script is designed for single-environment workflows where LOCALIZER,
pepstats (EMBOSS), and Perl are installed together.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def split_fasta(fasta: Path, chunk_dir: Path, chunk_size: int):
    """Split FASTA into chunk_00001.faa, ... of chunk_size records each.
    Reuses existing chunks when the manifest says they came from the same
    input file (same path, size) and the same chunk size -- that is what
    lets a rerun resume instead of starting over."""
    chunk_dir.mkdir(parents=True, exist_ok=True)
    manifest = chunk_dir / "manifest.txt"
    signature = f"{fasta}\t{fasta.stat().st_size}\t{chunk_size}\n"
    if manifest.exists() and manifest.read_text() == signature:
        chunks = sorted(chunk_dir.glob("chunk_*.faa"))
        return [(c, sum(1 for line in c.open() if line.startswith(">"))) for c in chunks]
    for old in chunk_dir.glob("chunk_*"):
        old.unlink()
    chunks, n, out, idx = [], 0, None, 0
    with fasta.open() as f:
        for line in f:
            if line.startswith(">"):
                if out is None or n == chunk_size:
                    if out:
                        out.close(); chunks.append((path, n))
                    idx += 1; n = 0
                    path = chunk_dir / f"chunk_{idx:05d}.faa"
                    out = path.open("w", newline="\n")
                n += 1
            if out:
                out.write(line)
    if out:
        out.close(); chunks.append((path, n))
    manifest.write_text(signature)
    return chunks


def run_chunked(base_cmd, fasta: Path, outdir: Path, raw_out: Path, chunk_size: int, chunk_jobs: int, env=None):
    """Run LOCALIZER chunk by chunk (chunk_jobs at a time), print a PROGRESS
    line after each, skip chunks already finished (.done marker), then join
    the chunk outputs in order into raw_out. LOCALIZER itself only writes
    when its whole input is done, so without chunks a proteome gives no
    progress signal for hours and a crash loses everything."""
    chunks = split_fasta(fasta, outdir / "chunks", chunk_size)
    total_seqs = sum(n for _, n in chunks)
    lock, state, t0 = threading.Lock(), {"chunks": 0, "seqs": 0}, time.time()

    def report(n):
        with lock:
            state["chunks"] += 1; state["seqs"] += n
            el = (time.time() - t0) / 60
            eta = el / state["seqs"] * (total_seqs - state["seqs"]) if state["seqs"] else 0
            print(f"PROGRESS {state['chunks']}/{len(chunks)} chunks ({state['seqs']}/{total_seqs} seqs) "
                  f"| elapsed {el:.1f} min | ETA {eta:.1f} min", file=sys.stderr, flush=True)

    def one(item):
        chunk, n = item
        done, raw = chunk.with_suffix(".done"), chunk.with_suffix(".raw.txt")
        if done.exists():
            report(n)
            return
        tmp = chunk.with_suffix(".raw.part")
        with tmp.open("w", newline="\n") as out, chunk.with_suffix(".stderr.log").open("w", newline="\n") as err:
            subprocess.run(base_cmd + ["-i", str(chunk)], check=True, stdout=out, stderr=err, env=env)
        tmp.replace(raw)
        done.touch()
        report(n)

    print(f"LOCALIZER: {total_seqs} seqs in {len(chunks)} chunks of <= {chunk_size}, {chunk_jobs} at a time",
          file=sys.stderr, flush=True)
    with ThreadPoolExecutor(max_workers=max(1, chunk_jobs)) as pool:
        list(pool.map(one, chunks))  # list() re-raises the first chunk failure
    with raw_out.open("w", newline="\n") as out:
        for chunk, _ in chunks:
            out.write(chunk.with_suffix(".raw.txt").read_text())


def java_env(localizer_script: Path, java_cpus: int):
    """Environment that caps the CPU threads of every Java (WEKA) process LOCALIZER starts.
    Without a cap, each JVM sizes its garbage-collector and JIT thread pools to ALL cores of the
    machine, so N parallel chunks use far more than N cores (8 sets x 4 chunks reached ~170 cores
    on a 128-core workstation). -XX:ActiveProcessorCount makes the JVM act as if it had
    `java_cpus` cores; -XX:+UseSerialGC removes the parallel GC threads.
    The JVM announces JAVA_TOOL_OPTIONS with one line on stderr, which LOCALIZER mixes into
    WEKA's output; only the parser patch in docs/install_localization.md skips it. So the cap
    is applied only when localization.py carries that patch; otherwise (None, reason)."""
    if java_cpus <= 0:
        return None, "java CPU cap off (--java-cpus 0)"
    loc_py = localizer_script.parent / "localization.py"
    try:
        patched = "line.split()[0].isdigit()" in loc_py.read_text()
    except OSError:
        patched = False
    if not patched:
        return None, (f"WARNING: java CPU cap NOT applied: {loc_py} lacks the parser patch from "
                      "docs/install_localization.md; each JVM may use every core")
    opts = f"-XX:ActiveProcessorCount={java_cpus} -XX:+UseSerialGC"
    env = dict(os.environ)
    env["JAVA_TOOL_OPTIONS"] = (env.get("JAVA_TOOL_OPTIONS", "") + " " + opts).strip()
    return env, f"java CPU cap: JAVA_TOOL_OPTIONS={env['JAVA_TOOL_OPTIONS']}"


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
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=2000,
        help="Run LOCALIZER on chunks of this many sequences, with a PROGRESS line per chunk "
             "and resume from finished chunks (default 2000; 0 = one call on the whole FASTA)",
    )
    parser.add_argument(
        "--chunk-jobs",
        type=int,
        default=4,
        help="Chunks run at the same time (default 4). Needs the parallel-run patch "
             "from docs/install_localization.md when > 1",
    )
    parser.add_argument(
        "--java-cpus",
        type=int,
        default=1,
        help="CPU threads each Java (WEKA) process may use (default 1; 0 = no cap). Total load is "
             "about chunk-jobs x java-cpus per input file. Needs the parser patch from "
             "docs/install_localization.md; skipped with a warning without it",
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

    base_cmd = [sys.executable, str(localizer_script), "-p" if args.mode == "plant" else "-e"]
    if args.mode == "effector":
        if args.mature_sequences:
            base_cmd.append("-M")
        if args.signal_peptide_length is not None:
            base_cmd.extend(["-S", str(args.signal_peptide_length)])
    cmd = base_cmd + ["-i", str(fasta)]

    env, env_msg = java_env(localizer_script, args.java_cpus)
    print(env_msg, file=sys.stderr, flush=True)
    log("STEP LOCALIZER predict", args.quiet)
    log("CMD: " + " ".join(cmd) + (f"  [chunked: {args.chunk_size} seqs, {args.chunk_jobs} at a time]"
                                    if args.chunk_size > 0 else ""), args.quiet)
    if not args.dry_run:
        if args.chunk_size > 0:
            run_chunked(base_cmd, fasta, outdir, raw_out, args.chunk_size, args.chunk_jobs, env)
        else:
            with raw_out.open("w", newline="\n") as out, stderr_log.open("w", newline="\n") as err:
                subprocess.run(cmd, check=True, stdout=out, stderr=err, env=env)

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
