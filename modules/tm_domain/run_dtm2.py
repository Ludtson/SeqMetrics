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
import contextlib
import csv
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

try:
    import fcntl  # POSIX only; on Windows the GPU lock is skipped (see gpu_slot)
except ImportError:
    fcntl = None

DEFAULT_GPU_LOCK = os.environ.get("SEQMETRICS_GPU_LOCK", "/tmp/seqmetrics_gpu.lock")
STDERR_TAIL_LINES = 40
# IUPAC ambiguity / non-standard residues -> closest standard residue. DeepTMHMM2's
# encoder only knows the 20 standard amino acids (+X) and crashes on e.g. 'J'
# (KeyError: 'J'). J = Leu/Ile (both hydrophobic), B = Asn/Asp, Z = Gln/Glu,
# U = selenocysteine, O = pyrrolysine. Applied to the copy fed to dtm2 only;
# each protein's number of replaced residues is reported in the output.
AMBIGUOUS = {"J": "L", "B": "N", "Z": "Q", "U": "C", "O": "K"}


def write_standardized(src, dest):
    """Copy src to dest replacing AMBIGUOUS residues; return {seq_id: n_replaced}."""
    table = str.maketrans("".join(AMBIGUOUS), "".join(AMBIGUOUS.values()))
    counts, cur = {}, None
    with open(src) as inp, open(dest, "w", newline="\n") as out:
        for line in inp:
            if line.startswith(">"):
                tok = line[1:].split()
                cur = tok[0] if tok else None
                if cur is not None:
                    counts.setdefault(cur, 0)
                out.write(line.rstrip("\r\n") + "\n")
            else:
                seq = line.rstrip("\r\n")
                if cur is not None:
                    counts[cur] += sum(seq.upper().count(a) for a in AMBIGUOUS)
                out.write(seq.translate(table).translate(str.maketrans("jbzuo", "lnqck")) + "\n")
    return counts


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


def resolve_device(requested, dtm2_path):
    """'auto' asks the dtm2 venv's own torch whether a CUDA GPU is usable;
    any failure to answer falls back to cpu."""
    if requested != "auto":
        return requested
    bindir = Path(dtm2_path).parent
    for name in ("python", "python3", "python.exe"):
        py = bindir / name
        if py.exists():
            try:
                r = subprocess.run([str(py), "-c", "import torch; print(torch.cuda.is_available())"],
                                   capture_output=True, text=True, timeout=120)
                return "cuda" if r.stdout.strip() == "True" else "cpu"
            except (OSError, subprocess.SubprocessError):
                return "cpu"
    return "cpu"


def gpu_free_mb():
    """Smallest free memory (MiB) across visible GPUs per nvidia-smi, or None
    if nvidia-smi isn't available -- the caller then skips the memory wait."""
    try:
        r = subprocess.run(["nvidia-smi", "--query-gpu=memory.free", "--format=csv,noheader,nounits"],
                           capture_output=True, text=True, timeout=30)
        values = [int(v) for v in r.stdout.split() if v.strip().isdigit()]
        return min(values) if r.returncode == 0 and values else None
    except (OSError, subprocess.SubprocessError, ValueError):
        return None


@contextlib.contextmanager
def gpu_slot(lock_path, min_free_mb, poll_s=60):
    """Serialize GPU use across processes, then wait for enough free memory.

    One dtm2 run loads the ESM2-650M model (~4-8 GB on the GPU); several at
    once on one card ran out of memory and crashed. The exclusive flock makes
    every other SeqMetrics tm_domain run wait its turn -- released by the OS
    even if this process dies, so nothing stays stuck. The free-memory check
    then covers GPU jobs the lock can't see (other users, other tools)."""
    lock = None
    if fcntl is None:
        print("[run_dtm2] GPU lock unavailable on this OS; running without it", file=sys.stderr)
    else:
        try:
            fd = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o666)
            with contextlib.suppress(OSError):
                os.chmod(lock_path, 0o666)  # shared machine: other users' runs need the same lock
            lock = os.fdopen(fd, "r+")
        except OSError as e:
            print(f"[run_dtm2] cannot open GPU lock {lock_path} ({e}); running without it", file=sys.stderr)
    try:
        if lock is not None:
            print(f"[run_dtm2] waiting for GPU lock {lock_path}", file=sys.stderr)
            fcntl.flock(lock, fcntl.LOCK_EX)
            print("[run_dtm2] GPU lock acquired", file=sys.stderr)
        while True:
            free = gpu_free_mb()
            if free is None or free >= min_free_mb:
                break
            print(f"[run_dtm2] GPU free {free} MiB < {min_free_mb} MiB needed; checking again in {poll_s} s",
                  file=sys.stderr)
            time.sleep(poll_s)
        yield
    finally:
        if lock is not None:
            fcntl.flock(lock, fcntl.LOCK_UN)
            lock.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, type=Path, help="Protein FASTA")
    ap.add_argument("-o", "--out", required=True, type=Path, help="Output TSV path")
    ap.add_argument("--dtm2", required=True, help="Path to the dtm2 executable")
    ap.add_argument("--device", default="auto",
                    help="auto (default: GPU if the dtm2 venv's torch can see one, else cpu), cpu, or cuda")
    ap.add_argument("--gpu-lock", default=DEFAULT_GPU_LOCK,
                    help="GPU runs only: lock file that queues tm_domain runs one at a time "
                         "(default: $SEQMETRICS_GPU_LOCK or /tmp/seqmetrics_gpu.lock)")
    ap.add_argument("--gpu-min-free-mb", type=int, default=8000,
                    help="GPU runs only: wait until this much GPU memory is free before starting (default 8000)")
    args = ap.parse_args()
    device = resolve_device(args.device, args.dtm2)
    print(f"[run_dtm2] device: {device}", file=sys.stderr)

    with tempfile.TemporaryDirectory() as tmp:
        std_input = Path(tmp) / "input_standard.faa"
        replaced = write_standardized(args.input, std_input)
        n_prot = sum(1 for v in replaced.values() if v)
        if n_prot:
            print(f"[run_dtm2] {n_prot} proteins contained ambiguity codes "
                  f"({sum(replaced.values())} residues) -> replaced for dtm2 only: "
                  + ", ".join(f"{a}->{b}" for a, b in AMBIGUOUS.items()), file=sys.stderr)
        cmd = [args.dtm2, str(std_input), tmp, "--device", device]
        slot = gpu_slot(args.gpu_lock, args.gpu_min_free_mb) if device == "cuda" else contextlib.nullcontext()
        with slot:
            result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            # dtm2's stderr includes a tqdm bar (one line per update -- ~30k
            # lines for a proteome); only the tail holds the actual error.
            tail = result.stderr.splitlines()[-STDERR_TAIL_LINES:]
            sys.stderr.write("\n".join(tail) + "\n")
            if "OutOfMemoryError" in result.stderr:
                sys.stderr.write("[run_dtm2] hint: GPU out of memory -- another GPU job may be running; "
                                 "raise --gpu-min-free-mb or use --device cpu\n")
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
        # LF only: this table's Seq_ID gets joined by exact match against other
        # feature tables downstream. newline="\n" alone is not enough -- csv.writer
        # writes its own \r\n terminator on every OS unless lineterminator is set.
        with open(args.out, "w", newline="\n") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            header = (["Seq_ID", "N_TMRs"]
                       + [f"Residues_{label}" for label in all_topology_labels]
                       + [f"MembraneProb_{col.replace(' ', '_')}" for col in all_membrane_cols]
                       + ["Ambiguous_residues_replaced"])
            writer.writerow(header)
            for seq_id in sorted(n_tmrs):
                row = [seq_id, n_tmrs[seq_id]]
                row += [topology_counts.get(seq_id, {}).get(label, 0) for label in all_topology_labels]
                row += [membrane_probs.get(seq_id, {}).get(col, "NA") for col in all_membrane_cols]
                row += [replaced.get(seq_id, 0)]
                writer.writerow(row)


if __name__ == "__main__":
    main()
