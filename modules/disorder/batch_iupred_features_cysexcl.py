#!/usr/bin/env python3

"""build_iupred_features.py

Usage:
  # Single file mode
  python build_iupred_features.py -i proteins.faa -o output.tsv -f

  # Batch mode (directory)
  python build_iupred_features.py -i input_folder/ -o output_folder/ -f --concat --jobs 4

Description:
  Reads FASTA files (single or batch), runs IUPred3 (long, short, glob)
  and ANCHOR2 via iupred3_lib, and aggregates per-protein features.

  Removes dependencies on non-native libraries (like pandas) for easier portability.
"""

import argparse
import sys
import os
import csv
import glob
import time
import statistics
import concurrent.futures
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from datetime import datetime

# Try to import the library; handle failure gracefully if not in env
try:
    import iupred3_lib
except ImportError:
    sys.stderr.write("ERROR: iupred3_lib not found. Make sure you are in the correct environment.\n")
    sys.exit(1)


# ----------------------
# Logging & Utilities
# ----------------------

def log(msg: str) -> None:
    """Print a timestamped message to stderr."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    sys.stderr.write(f"[{ts}] {msg}\n")

def safe_mean(values: List[float]) -> float:
    return statistics.mean(values) if values else 0.0

def safe_median(values: List[float]) -> float:
    return statistics.median(values) if values else 0.0

def safe_stdev(values: List[float]) -> float:
    return statistics.stdev(values) if len(values) > 1 else 0.0

def safe_max(values: List[float]) -> float:
    return max(values) if values else 0.0

def safe_min(values: List[float]) -> float:
    return min(values) if values else 0.0

# ----------------------
# FASTA parsing
# ----------------------

def read_fasta(path: str):
    """Simple FASTA parser: yields (id, seq) for each record."""
    sid = None
    seq_chunks = []
    
    try:
        with open(path, 'r') as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if sid is not None:
                        yield sid, "".join(seq_chunks)
                    sid = line[1:].split()[0]
                    seq_chunks = []
                else:
                    seq_chunks.append(line)
        if sid is not None:
            yield sid, "".join(seq_chunks)
    except Exception as e:
        log(f"Error reading FASTA {path}: {e}")

# ----------------------
# Segment utilities (Native Python)
# ----------------------

def find_segments(bool_list: List[bool], min_len: int) -> List[Tuple[int, int]]:
    """
    Return list of (start, end) 1-based inclusive ranges
    where bool_list is True for runs of length >= min_len.
    """
    segments = []
    in_seg = False
    start_idx = -1

    # Iterate through the boolean mask
    for i, val in enumerate(bool_list):
        if val and not in_seg:
            in_seg = True
            start_idx = i
        elif not val and in_seg:
            # End of a segment
            length = i - start_idx
            if length >= min_len:
                # Store 1-based indices (start + 1, current i)
                segments.append((start_idx + 1, i))
            in_seg = False
    
    # Check if a segment continues to the very end
    if in_seg:
        length = len(bool_list) - start_idx
        if length >= min_len:
            segments.append((start_idx + 1, len(bool_list)))

    return segments

# ----------------------
# Feature computation
# ----------------------

def compute_features_for_protein(
    pid: str,
    seq: str,
    long_scores: List[float],
    short_scores: List[float],
    glob_scores: List[float],
    anchor_scores: Optional[List[float]],
    config: Dict[str, Any]
) -> Dict[str, Any]:
    """Compute per-protein features using standard library only."""
    
    # Unpack config
    disorder_thresh = config['disorder_threshold']
    binding_thresh = config['binding_threshold']
    min_idr = config['min_idr_len']
    min_bind = config['min_bind_len']
    do_advanced = config['compute_engineered']

    length = len(seq)

    # --- Cysteine exclusion (Option D) ---
    # IUPred runs on the REAL, unmodified sequence above (Cys is a standard
    # residue with well-defined behavior in the algorithm; nothing about the
    # sequence itself is altered). What changes here is which residues are
    # trusted for the AGGREGATE statistics: since IUPred cannot see disulfide
    # bonds, a Cys held rigid by a real disulfide can still be scored as
    # disorder-prone from local composition alone. Rather than guessing a
    # substitute identity (which risks its own composition bias — Ser, for
    # instance, is itself disorder-promoting, and DNG/AG measurably differ in
    # Cys content, so a substitution's bias would land unevenly between
    # classes), Cys positions are simply excluded from every statistic below.
    # No sequence length change, no invented residue identity.
    keep_idx = [i for i, ch in enumerate(seq) if ch != "C"]
    n_kept = len(keep_idx)

    # Initialize feature dict
    feat = {
        "protein_id": pid,
        "length": length,
        "n_cys_excluded": length - n_kept,
    }

    # 1. Basic Stats for all modes (over non-cysteine positions only)
    modes = {
        "long": long_scores,
        "short": short_scores,
        "glob": glob_scores
    }

    for mode_name, scores in modes.items():
        # Truncate to length just in case library returns extra, then drop Cys positions
        full = scores[:length]
        s = [full[i] for i in keep_idx]

        feat[f"{mode_name}_mean"] = safe_mean(s)
        feat[f"{mode_name}_median"] = safe_median(s)
        feat[f"{mode_name}_min"] = safe_min(s)
        feat[f"{mode_name}_max"] = safe_max(s)
        feat[f"{mode_name}_std"] = safe_stdev(s)

        # Fraction >= threshold, of non-cysteine residues only
        count_ge = sum(1 for x in s if x >= disorder_thresh)
        feat[f"{mode_name}_frac_ge_{disorder_thresh:.2f}"] = count_ge / n_kept if n_kept > 0 else 0.0

    # ---------------------------------------------------------
    # If engineered features are NOT requested, return early
    # ---------------------------------------------------------
    if not do_advanced:
        return feat

    # ---------------------------------------------------------
    # Advanced / Engineered Features
    # ---------------------------------------------------------
    
    # Focus on "long" scores for disorder segments.
    # dis_mask_full is computed over EVERY real position (including Cys —
    # the score there is real, just not trusted for aggregation), so that
    # positional windows (N/C/mid) still reflect true sequence architecture.
    l_scores = modes["long"][:length]
    dis_mask_full = [x >= disorder_thresh for x in l_scores]
    keep_set = set(keep_idx)

    feat["disorder_frac"] = (
        sum(1 for i in keep_idx if dis_mask_full[i]) / n_kept if n_kept > 0 else 0.0
    )

    # Region fractions (N-term, C-term, Mid) — restricted to non-cysteine
    # residues within each positional window.
    nwindow = min(50, length)
    cwindow = min(50, length)
    mid_start = nwindow
    mid_end = length - cwindow

    n_part = [dis_mask_full[i] for i in range(nwindow) if i in keep_set]
    c_part = [dis_mask_full[i] for i in range(length - cwindow, length) if i in keep_set]
    mid_part = ([dis_mask_full[i] for i in range(mid_start, mid_end) if i in keep_set]
                if mid_end > mid_start else [])

    feat["N_frac_disordered"] = safe_mean([float(x) for x in n_part])
    feat["C_frac_disordered"] = safe_mean([float(x) for x in c_part])
    feat["mid_frac_disordered"] = safe_mean([float(x) for x in mid_part])
    feat["NC_disorder_asymmetry"] = feat["N_frac_disordered"] - feat["C_frac_disordered"]

    # IDR Segments — computed on the COMPACTED mask (cysteine positions
    # omitted entirely, not treated as an "ordered" break), so a disordered
    # stretch interrupted only by a single Cys is counted as one continuous
    # run rather than artificially split.
    dis_mask = [dis_mask_full[i] for i in keep_idx]
    idr_segments = find_segments(dis_mask, min_idr)
    feat["n_idr_segments"] = len(idr_segments)
    
    idr_lengths = [(end - start + 1) for start, end in idr_segments]
    
    if idr_lengths:
        feat["total_idr_len"] = sum(idr_lengths)
        feat["longest_idr_len"] = max(idr_lengths)
        feat["mean_idr_len"] = safe_mean(idr_lengths)
        feat["n_idr_ge30"] = sum(1 for L in idr_lengths if L >= 30)
        feat["n_idr_ge50"] = sum(1 for L in idr_lengths if L >= 50)
        feat["has_long_idr_ge30"] = 1 if feat["n_idr_ge30"] > 0 else 0
        feat["has_long_idr_ge50"] = 1 if feat["n_idr_ge50"] > 0 else 0
    else:
        feat["total_idr_len"] = 0
        feat["longest_idr_len"] = 0
        feat["mean_idr_len"] = 0.0
        feat["n_idr_ge30"] = 0
        feat["n_idr_ge50"] = 0
        feat["has_long_idr_ge30"] = 0
        feat["has_long_idr_ge50"] = 0

    # ANCHOR2 (Binding) Features — same cysteine-exclusion treatment
    if anchor_scores:
        a_scores_full = anchor_scores[:length]
        a_scores = [a_scores_full[i] for i in keep_idx]
        bind_mask_full = [x >= binding_thresh for x in a_scores_full]
        bind_mask = [bind_mask_full[i] for i in keep_idx]  # compacted, for segment-finding

        feat["anchor2_frac_binding"] = sum(bind_mask) / n_kept if n_kept > 0 else 0.0
        feat["anchor2_max"] = safe_max(a_scores)
        feat["anchor2_mean"] = safe_mean(a_scores)

        bind_segments = find_segments(bind_mask, min_bind)
        feat["n_bind_segments"] = len(bind_segments)

        bind_lengths = [(end - start + 1) for start, end in bind_segments]

        if bind_lengths:
            feat["total_bind_len"] = sum(bind_lengths)
            feat["longest_bind_len"] = max(bind_lengths)
            feat["mean_bind_len"] = safe_mean(bind_lengths)

            # Overlap calculation on the compacted (non-cysteine) index space,
            # consistent with dis_mask above (also compacted).
            dis_indices = {i for i, val in enumerate(dis_mask) if val}
            bind_indices = {i for i, val in enumerate(bind_mask) if val}

            overlap_count = len(dis_indices.intersection(bind_indices))

            feat["n_disordered_binding_residues"] = overlap_count

            denom_dis = len(dis_indices) if len(dis_indices) > 0 else 1
            denom_bind = len(bind_indices) if len(bind_indices) > 0 else 1

            feat["disordered_binding_frac_of_disordered"] = overlap_count / denom_dis
            feat["disordered_binding_frac_of_binding"] = overlap_count / denom_bind
        else:
            feat["total_bind_len"] = 0
            feat["longest_bind_len"] = 0
            feat["mean_bind_len"] = 0.0
            feat["n_disordered_binding_residues"] = 0
            feat["disordered_binding_frac_of_disordered"] = 0.0
            feat["disordered_binding_frac_of_binding"] = 0.0
    else:
        # Default null values for anchor
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
# Worker Function
# ----------------------

def process_single_file(fasta_path: str, output_path: str, config: Dict[str, Any]) -> Tuple[str, List[Dict]]:
    """
    Worker function to process one FASTA file.
    Returns (output_path, list_of_features) for aggregation.
    """
    try:
        # Read seqs >= 30
        records = []
        for pid, seq in read_fasta(fasta_path):
            if len(seq) >= 30:
                records.append((pid, seq))
        
        if not records:
            # Create empty file
            with open(output_path, 'w', newline='') as f:
                pass
            return output_path, []

        all_feats = []
        smoothing = config.get("smoothing", "medium")

        for pid, seq in records:
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
                config=config
            )
            feat['source_file'] = Path(fasta_path).name # Add provenance
            all_feats.append(feat)

        # Write to individual TSV
        if all_feats:
            keys = all_feats[0].keys()
            with open(output_path, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=keys, delimiter='\t')
                writer.writeheader()
                writer.writerows(all_feats)
        
        return output_path, all_feats

    except Exception as e:
        log(f"CRITICAL ERROR processing {fasta_path}: {e}")
        return output_path, []

# ----------------------
# Main / Batch Logic
# ----------------------

def parse_args():
    p = argparse.ArgumentParser(description="Build IUPred3/ANCHOR2 features (Batch/Parallel)")
    
    # Input/Output
    p.add_argument("-i", "--input", required=True, help="Input FASTA file OR directory containing FASTA files")
    p.add_argument("-o", "--output", required=True, help="Output filename (if input is file) OR directory (if input is dir)")
    
    # Flags
    p.add_argument("-f", "--features", action="store_true", help="Enable calculation of engineered features (segments, IDRs, asymmetry)")
    p.add_argument("--concat", action="store_true", help="In batch mode, generate a master concatenated file")
    p.add_argument("-j", "--jobs", type=int, default=1, help="Number of parallel jobs (default 1)")
    
    # Parameters
    p.add_argument("--disorder-threshold", type=float, default=0.5, help="Threshold for disorder (0.5)")
    p.add_argument("--binding-threshold", type=float, default=0.5, help="Threshold for binding (0.5)")
    p.add_argument("--min-idr-len", type=int, default=10, help="Min length for IDR (10)")
    p.add_argument("--min-bind-len", type=int, default=5, help="Min length for binding segment (5)")
    p.add_argument("--smoothing", choices=["no", "medium", "strong"], default="medium", help="Smoothing type")

    return p.parse_args()

def main():
    args = parse_args()
    
    # Config dictionary to pass to workers
    config = {
        "disorder_threshold": args.disorder_threshold,
        "binding_threshold": args.binding_threshold,
        "min_idr_len": args.min_idr_len,
        "min_bind_len": args.min_bind_len,
        "smoothing": args.smoothing,
        "compute_engineered": args.features
    }

    input_path = Path(args.input)
    output_path = Path(args.output)
    
    # -----------------------
    # 1. Determine Mode
    # -----------------------
    is_batch = input_path.is_dir()
    
    if is_batch:
        log(f"Batch mode detected. Input dir: {input_path}")
        
        # Setup directories
        # Structure: 
        #   output_path/
        #       individual_files/
        #       master_combined.tsv (optional)
        
        individual_dir = output_path / "individual_files"
        individual_dir.mkdir(parents=True, exist_ok=True)
        
        # Find files
        extensions = ['*.fa', '*.faa', '*.fasta']
        input_files = []
        for ext in extensions:
            input_files.extend(input_path.glob(ext))
        
        if not input_files:
            log("No FASTA files found in input directory.")
            sys.exit(0)
            
        log(f"Found {len(input_files)} files to process.")
        
        # Prepare Tasks
        tasks = []
        for f in input_files:
            # Output name: input_name + .tsv
            out_name = f.stem + "_features.tsv"
            out_f = individual_dir / out_name
            tasks.append((str(f), str(out_f)))
            
        # Run Parallel
        all_results_flat = []
        
        start_time = time.time()
        
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as executor:
            # Map returns iterator in order
            futures = [executor.submit(process_single_file, inp, outp, config) for inp, outp in tasks]
            
            total = len(futures)
            completed = 0
            
            for future in concurrent.futures.as_completed(futures):
                out_p, file_feats = future.result()
                completed += 1
                if completed % 5 == 0 or completed == total:
                    log(f"Progress: {completed}/{total} files processed.")
                
                if args.concat:
                    all_results_flat.extend(file_feats)
        
        elapsed = time.time() - start_time
        log(f"Batch processing finished in {elapsed:.2f} seconds.")
        
        # Concatenate
        if args.concat and all_results_flat:
            master_file = output_path / "master_combined.tsv"
            log(f"Writing concatenated master file to {master_file}...")
            keys = all_results_flat[0].keys()
            with open(master_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=keys, delimiter='\t')
                writer.writeheader()
                writer.writerows(all_results_flat)
        elif args.concat:
            log("Warning: No results found to concatenate.")

    else:
        # -----------------------
        # 2. Single File Mode
        # -----------------------
        if not input_path.exists():
            log(f"Error: Input file {input_path} does not exist.")
            sys.exit(1)
            
        # Ensure output dir exists if path has parents
        if output_path.parent:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
        log(f"Single file mode. Processing {input_path} -> {output_path}")
        process_single_file(str(input_path), str(output_path), config)
        log("Done.")

if __name__ == "__main__":
    main()
    