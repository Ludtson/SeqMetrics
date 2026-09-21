#!/usr/bin/env python3
# ==============================================================================
# Script Name:  extract_heg_ids.py
# Description:  Identifies and extracts Highly Expressed Gene (HEG) seed IDs 
#               by parsing CDS FASTA headers and cross-referencing functional 
#               annotations inside corresponding GFF3 files via a BFS queue.
# Credits:      Adekola Owoyemi, Casola Lab, TAMU
# Version:      0.1.0 (pre-release, unshared)
# ==============================================================================

import re
import sys
import os
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

"""
extract_heg_ids.py (CDS-Driven Version)

Identifies Highly Expressed Genes (HEGs) by:
1. Reading IDs from a CDS FASTA file (longest isoforms).
2. Verifying their functional status in a corresponding GFF3 file.
3. Using hierarchical propagation (Gene <-> mRNA <-> CDS) to find keywords.

Keywords include Ribosomal proteins, Elongation factors, Rubisco, and Chaperones.

Usage:
  python extract_heg_ids.py --fasta-dir <dir> --gff-dir <dir> --out-dir <dir> --threads 25
"""

def get_args():
    parser = argparse.ArgumentParser(
        description="Extract HEG IDs from CDS FASTA by verifying status in GFF3.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--fasta-dir", required=True, help="Directory containing CDS FASTA files")
    parser.add_argument("--gff-dir", required=True, help="Directory containing GFF3 files")
    parser.add_argument("--out-dir", required=True, help="Output directory for ID lists")
    parser.add_argument("--threads", type=int, default=25, help="Number of parallel threads")
    parser.add_argument("--suffix", default="_heg_ids.txt", help="Suffix for output files")
    return parser.parse_args()

def parse_gff_heg_keys(gff_path):
    """
    Optimized GFF parser: 
    1. Builds hierarchy (ID->Parent) and identifies seeds.
    2. Propagates HEG status using a queue.
    3. Extracts keys ONLY for HEG features.
    """
    patterns = [
        r"ribosomal\s+protein", r"\bRPS\d+", r"\bRPL\d+",
        r"elongation\s+factor", r"\beEF\d", r"\btuf",
        r"rubisco", r"ribulose", r"rbcL", r"rbcS",
        r"chlorophyll\s+a/b", r"\bCAB\b", r"\bLHCB",
        r"GAPDH", r"glyceraldehyde", r"enolase", r"\beno\b",
        r"HSP70", r"HSP60", r"chaperone"
    ]
    heg_regex = re.compile("|".join(patterns), re.IGNORECASE)

    # Relationships and Seeds
    parent_map = {} # ChildID -> [ParentIDs]
    children_map = defaultdict(list) # ParentID -> [ChildIDs]
    heg_seeds = set()
    
    # Temporarily store attributes to avoid double parsing later
    # Only for IDs that end up being HEG
    all_features = {} # ID -> attr_string

    re_id = re.compile(r'ID=([^;]+)')
    re_parent = re.compile(r'Parent=([^;]+)')

    with open(gff_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) < 9: continue
            
            attr_str = parts[8]
            m_id = re_id.search(attr_str)
            m_parent = re_parent.search(attr_str)
            
            fid = m_id.group(1) if m_id else None
            # If no ID, we use a coordinate-based hash to track it
            if not fid:
                fid = f"{parts[0]}_{parts[3]}_{parts[4]}"

            all_features[fid] = attr_str

            if m_parent:
                parents = m_parent.group(1).split(',')
                parent_map[fid] = parents
                for p in parents:
                    children_map[p].append(fid)

            if heg_regex.search(attr_str):
                heg_seeds.add(fid)

    # Fast Propagation using a Queue (BFS)
    final_heg_ids = set()
    queue = list(heg_seeds)
    
    while queue:
        curr = queue.pop(0)
        if curr in final_heg_ids:
            continue
        final_heg_ids.add(curr)
        
        # Add Parents
        if curr in parent_map:
            for p in parent_map[curr]:
                if p not in final_heg_ids:
                    queue.append(p)
        
        # Add Children
        if curr in children_map:
            for c in children_map[curr]:
                if c not in final_heg_ids:
                    queue.append(c)

    # Final targeted harvesting of keys
    heg_keys = set()
    for fid in final_heg_ids:
        attr_str = all_features.get(fid, "")
        for item in attr_str.split(';'):
            if '=' in item:
                key_type, val = item.split('=', 1)
                if key_type.strip().lower() in ['id', 'name', 'parent', 'gene_id', 'transcript_id', 'protein_id']:
                    # Handle comma-separated attribute values
                    for sub_val in val.split(','):
                        clean_val = sub_val.strip()
                        if clean_val:
                            heg_keys.add(clean_val)    
    return heg_keys

def clean_id(identifier):
    """
    Normalizes IDs by removing common prefixes and suffixes.
    e.g. 'gene:Bo1g007080' -> 'Bo1g007080'
    e.g. 'Bo1g007080.1' -> 'Bo1g007080'
    """
    if not identifier: return ""
    # Remove common prefixes
    clean = re.sub(r'^(gene|rna|transcript|cds)[:\-]', '', identifier, flags=re.IGNORECASE)
    # Remove common version suffixes (.1, .t1, -mRNA-1)
    clean = re.sub(r'(\.\d+|\.t\d+|\-mRNA\-\d+)$', '', clean)
    return clean

def process_species(task_data):
    fasta_path, gff_path, out_path = task_data
    species = os.path.basename(fasta_path).replace('.fna', '')
    
    try:
        # 1. Parse GFF to get HEG Keys
        heg_keys = parse_gff_heg_keys(gff_path)
        # Create a set of cleaned keys for fuzzy matching
        cleaned_heg_keys = {clean_id(k) for k in heg_keys if k}
        
        # 2. Process FASTA
        heg_ids_found = []
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Get ID (first word after >)
                    fasta_id = line.strip().split()[0][1:]
                    
                    # Try exact match first
                    if fasta_id in heg_keys:
                        heg_ids_found.append(fasta_id)
                    else:
                        # Try cleaned/fuzzy match
                        if clean_id(fasta_id) in cleaned_heg_keys:
                            heg_ids_found.append(fasta_id)

        # 3. Write Output
        # newline='\n': force LF-only regardless of host OS. Without this,
        # running on Windows Python writes CRLF, and calculate_indices.sh's
        # awk lookup (running in WSL/Linux) then keys on "ID\r" -- which
        # never matches a CDS FASTA's clean "ID" header, so every sequence
        # extraction silently comes back empty and cai.coa never gets built.
        # Confirmed via a real run: 98 HEG IDs found here, zero sequences
        # extracted downstream, until this fix.
        with open(out_path, 'w', newline='\n') as out:
            for hid in sorted(list(set(heg_ids_found))):
                out.write(f"{hid}\n")
        
        return f"{species}: Extracted {len(heg_ids_found)} HEG IDs."
    except Exception as e:
        import traceback
        return f"{species}: ERROR - {str(e)}\n{traceback.format_exc()}"

def main():
    args = get_args()
    
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # Match FASTA and GFF files
    tasks = []
    fasta_files = [f for f in os.listdir(args.fasta_dir) if f.endswith('.fna')]
    print(f"Found {len(fasta_files)} FASTA files in {args.fasta_dir}")
    
    for ffile in fasta_files:
        species_base = ffile.replace('.fna', '')
        gff_file = species_base + ".gff"
        gff_path = os.path.join(args.gff_dir, gff_file)
        
        if os.path.exists(gff_path):
            fasta_path = os.path.join(args.fasta_dir, ffile)
            out_path = os.path.join(args.out_dir, species_base + args.suffix)
            tasks.append((fasta_path, gff_path, out_path))
        else:
            print(f"Warning: Missing GFF for {species_base}")

    print(f"Starting batch processing with {args.threads} threads...")
    
    total_files = len(tasks)
    files_with_hegs = 0
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(process_species, t) for t in tasks]
        for future in as_completed(futures):
            res = future.result()
            print(res)
            if "Extracted" in res:
                count = int(res.split("Extracted ")[1].split()[0])
                if count > 0:
                    files_with_hegs += 1

    print("\n" + "="*30)
    print("FINAL SUMMARY")
    print("="*30)
    print(f"Total Species processed: {total_files}")
    print(f"Species with HEGs found: {files_with_hegs}")
    print(f"Species with 0 HEGs:    {total_files - files_with_hegs}")
    print("="*30)

if __name__ == "__main__":
    main()
    