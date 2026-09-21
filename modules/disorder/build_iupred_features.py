#!/usr/bin/env python

"""build_iupred_features.py

Usage:
  python build_iupred_features.py proteins.faa output.tsv

Reads a multi-FASTA of proteins, runs IUPred3 (long, short, glob)
and ANCHOR2 via iupred3_lib, and aggregates per-protein features.

Sequences with length < 30 amino acids are SKIPPED and do not
appear in the output.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from datetime import datetime

import pandas as pd
import iupred3_lib  # must be importable in this env (same as iupred3.py)


# ----------------------
# Logging helper
# ----------------------

def log(msg: str) -> None:
    """Print a timestamped message to stderr."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sys.stderr.write(f"[{ts}] {msg}\n")


# ----------------------
# FASTA parsing
# ----------------------

def read_fasta(path: str):
    """Simple FASTA parser: yields (id, seq) for each record."""
    sid = None
    seq_chunks: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                # flush previous
                if sid is not None:
                    yield sid, "".join(seq_chunks)
                sid = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if sid is not None:
        yield sid, "".join(seq_chunks)


# ----------------------
# Segment utilities
# ----------------------

def find_segments(mask: pd.Series, min_len: int) -> List[Tuple[int, int]]:
    """Return list of (start, end) 1-based inclusive ranges
    where mask is True for runs of length >= min_len.
    """
    segments: List[Tuple[int, int]] = []
    in_seg = False
    start = None

    for idx, val in mask.items():
        if val and not in_seg:
            in_seg = True
            start = idx
        elif not val and in_seg:
            end = idx - 1
            if end - start + 1 >= min_len:
                segments.append((start, end))
            in_seg = False
            start = None

    if in_seg and start is not None:
        end = mask.index[-1]
        if end - start + 1 >= min_len:
            segments.append((start, end))

    return segments


# ----------------------
# Feature computation
# ----------------------

def compute_features_for_protein(
    pid: str,
    seq: str,
    long_scores,
    short_scores,
    glob_scores,
    anchor_scores,
    disorder_threshold: float,
    binding_threshold: float,
    min_idr_len: int,
    min_bind_len: int,
) -> Dict:
    """Compute per-protein features from per-residue scores."""
    feat: Dict = {}
    feat["protein_id"] = pid

    length = len(seq)
    feat["length"] = length

    pos = pd.Series(range(1, length + 1), name="pos")
    aa = pd.Series(list(seq), name="aa")

    df = pd.DataFrame(
        {
            "pos": pos,
            "aa": aa,
            "long": pd.Series(long_scores[:length]),
            "short": pd.Series(short_scores[:length]),
            "glob": pd.Series(glob_scores[:length]),
            "anchor2": pd.Series(anchor_scores[:length]) if anchor_scores is not None else None,
        }
    )

    # Per-mode basic stats
    for mode in ("long", "short", "glob"):
        s = df[mode]
        base = mode
        feat[f"{base}_mean"] = float(s.mean())
        feat[f"{base}_median"] = float(s.median())
        feat[f"{base}_min"] = float(s.min())
        feat[f"{base}_max"] = float(s.max())
        feat[f"{base}_std"] = float(s.std(ddof=0))
        feat[f"{base}_frac_ge_{disorder_threshold:.2f}"] = float((s >= disorder_threshold).mean())

    # Disorder-related features, using long mode
    s = df["long"]
    dis_mask = s >= disorder_threshold
    feat["disorder_frac"] = float(dis_mask.mean())

    # N/C/central fractions
    nwindow = min(50, length)
    cwindow = min(50, length)

    n_mask = df["pos"] <= nwindow
    c_mask = df["pos"] > (length - cwindow)
    mid_mask = ~(n_mask | c_mask)

    feat["N_frac_disordered"] = float(dis_mask[n_mask].mean()) if n_mask.any() else 0.0
    feat["C_frac_disordered"] = float(dis_mask[c_mask].mean()) if c_mask.any() else 0.0
    feat["mid_frac_disordered"] = float(dis_mask[mid_mask].mean()) if mid_mask.any() else 0.0
    feat["NC_disorder_asymmetry"] = feat["N_frac_disordered"] - feat["C_frac_disordered"]

    # Disordered segments
    dis_mask_indexed = pd.Series(dis_mask.values, index=df["pos"])
    idr_segments = find_segments(dis_mask_indexed, min_idr_len)

    feat["n_idr_segments"] = len(idr_segments)
    if idr_segments:
        lengths = [e - b + 1 for b, e in idr_segments]
        feat["total_idr_len"] = int(sum(lengths))
        feat["longest_idr_len"] = int(max(lengths))
        feat["mean_idr_len"] = float(pd.Series(lengths).mean())
        feat["n_idr_ge30"] = int(sum(L >= 30 for L in lengths))
        feat["n_idr_ge50"] = int(sum(L >= 50 for L in lengths))
        feat["has_long_idr_ge30"] = int(any(L >= 30 for L in lengths))
        feat["has_long_idr_ge50"] = int(any(L >= 50 for L in lengths))
    else:
        feat["total_idr_len"] = 0
        feat["longest_idr_len"] = 0
        feat["mean_idr_len"] = 0.0
        feat["n_idr_ge30"] = 0
        feat["n_idr_ge50"] = 0
        feat["has_long_idr_ge30"] = 0
        feat["has_long_idr_ge50"] = 0

    # Binding-related features (ANCHOR2)
    if anchor_scores is not None:
        a = df["anchor2"]
        bind_mask = a >= binding_threshold

        feat["anchor2_frac_binding"] = float(bind_mask.mean())
        feat["anchor2_max"] = float(a.max())
        feat["anchor2_mean"] = float(a.mean())

        bind_mask_indexed = pd.Series(bind_mask.values, index=df["pos"])
        bind_segments = find_segments(bind_mask_indexed, min_bind_len)

        feat["n_bind_segments"] = len(bind_segments)
        if bind_segments:
            blengths = [e - b + 1 for b, e in bind_segments]
            feat["total_bind_len"] = int(sum(blengths))
            feat["longest_bind_len"] = int(max(blengths))
            feat["mean_bind_len"] = float(pd.Series(blengths).mean())

            # Overlap with IDRs
            dis_positions = set(df.loc[dis_mask, "pos"].tolist())
            bind_positions = set(df.loc[bind_mask, "pos"].tolist())
            overlap = dis_positions & bind_positions

            feat["n_disordered_binding_residues"] = int(len(overlap))
            feat["disordered_binding_frac_of_disordered"] = float(
                len(overlap) / max(len(dis_positions), 1)
            )
            feat["disordered_binding_frac_of_binding"] = float(
                len(overlap) / max(len(bind_positions), 1)
            )
        else:
            feat["total_bind_len"] = 0
            feat["longest_bind_len"] = 0
            feat["mean_bind_len"] = 0.0
            feat["n_disordered_binding_residues"] = 0
            feat["disordered_binding_frac_of_disordered"] = 0.0
            feat["disordered_binding_frac_of_binding"] = 0.0
    else:
        feat["anchor2_frac_binding"] = 0.0
        feat["anchor2_max"] = 0.0
        feat["anchor2_mean"] = 0.0
        feat["n_bind_segments"] = 0
        feat["total_bind_len"] = 0
        feat["longest_bind_len"] = 0
        feat["mean_bind_len"] = 0.0
        feat["n_disordered_binding_residues"] = 0
        feat["disordered_binding_frac_of_disordered"] = 0.0
        feat["disordered_binding_frac_of_binding"] = 0.0

    return feat


# ----------------------
# Main
# ----------------------

def parse_args():
    p = argparse.ArgumentParser(description="Build IUPred3/ANCHOR2 features per protein")
    p.add_argument("fasta", help="Input multi-FASTA file of proteins")
    p.add_argument("output", help="Output TSV/CSV file (extension decides separator)")
    p.add_argument(
        "--disorder-threshold",
        type=float,
        default=0.5,
        help="Threshold for disordered residues (default 0.5)",
    )
    p.add_argument(
        "--binding-threshold",
        type=float,
        default=0.5,
        help="Threshold for binding residues (ANCHOR2; default 0.5)",
    )
    p.add_argument(
        "--min-idr-len",
        type=int,
        default=10,
        help="Minimum length for a disordered segment (default 10)",
    )
    p.add_argument(
        "--min-bind-len",
        type=int,
        default=5,
        help="Minimum length for a binding segment (default 5)",
    )
    p.add_argument(
        "--smoothing",
        choices=["no", "medium", "strong"],
        default="medium",
        help="IUPred3 smoothing type (default medium)",
    )
    return p.parse_args()


def main():
    args = parse_args()

    fasta_path = Path(args.fasta)
    out_path = Path(args.output)

    log(f"Starting build_iupred_features.py")
    log(f"Input FASTA: {fasta_path}")
    log(f"Output file: {out_path}")

    if not fasta_path.is_file():
        log(f"ERROR: Input FASTA not found: {fasta_path}")
        sys.exit(1)

    # First pass: count total sequences and how many are long enough
    all_records: List[Tuple[str, str]] = list(read_fasta(str(fasta_path)))
    total_seqs = len(all_records)
    log(f"Total sequences in FASTA: {total_seqs}")

    # Filter by length >= 30
    records = [(pid, seq) for pid, seq in all_records if len(seq) >= 30]
    long_enough = len(records)
    skipped = total_seqs - long_enough

    if skipped > 0:
        log(f"Skipping {skipped} sequences with length < 30")

    if long_enough == 0:
        log("No sequences >= 30 aa to process; writing empty file.")
        pd.DataFrame().to_csv(
            out_path,
            sep="\t" if out_path.suffix.lower() in (".tsv", ".txt") else ",",
        )
        return

    log(f"Will process {long_enough} sequences (length >= 30)")

    feats: List[Dict] = []

    # Setup progress checkpoints at ~10%, 25%, 50%, 75%, 100%
    checkpoints = [0.10, 0.25, 0.50, 0.75, 1.0]
    next_cp_index = 0

    for idx, (pid, seq) in enumerate(records, start=1):
        # Optional: log a few sequence IDs at the start
        if idx == 1:
            log(f"Processing first sequence: {pid} (length {len(seq)})")

        smoothing = args.smoothing

        long_scores = iupred3_lib.iupred(seq, "long", smoothing=smoothing)[0]
        short_scores = iupred3_lib.iupred(seq, "short", smoothing=smoothing)[0]
        glob_scores = iupred3_lib.iupred(seq, "glob", smoothing=smoothing)[0]
        anchor_scores = iupred3_lib.anchor2(seq)

        feat = compute_features_for_protein(
            pid=pid,
            seq=seq,
            long_scores=long_scores,
            short_scores=short_scores,
            glob_scores=glob_scores,
            anchor_scores=anchor_scores,
            disorder_threshold=args.disorder_threshold,
            binding_threshold=args.binding_threshold,
            min_idr_len=args.min_idr_len,
            min_bind_len=args.min_bind_len,
        )
        feats.append(feat)

        # Progress logging
        frac = idx / long_enough
        if next_cp_index < len(checkpoints) and frac >= checkpoints[next_cp_index]:
            pct = int(checkpoints[next_cp_index] * 100)
            log(f"Progress: processed {idx}/{long_enough} sequences (~{pct}%)")
            next_cp_index += 1

    feat_df = pd.DataFrame(feats).set_index("protein_id")

    if out_path.suffix.lower() in (".tsv", ".txt"):
        feat_df.to_csv(out_path, sep="\t")
    else:
        feat_df.to_csv(out_path, sep=",")

    log(f"Wrote features to {out_path}")
    log("Finished build_iupred_features.py")


if __name__ == "__main__":
    main()
