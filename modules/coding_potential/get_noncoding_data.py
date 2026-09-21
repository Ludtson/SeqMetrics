#!/usr/bin/env python3
"""
get_noncoding_data.py

Builds a high-confidence, conservative non-coding training set from genome + GFF.

Required inputs:
 -p / --protein-fasta      Protein FASTA file with known coding sequences.
 -g / --genome-fasta       Genome FASTA file (nucleotide).
 -a / --annotation-gff     Annotation file in GFF3 format containing CDS features.
 -o / --out-dir            Output directory where all result files will be saved.

Optional:
 -l / --min-noncoding-length Minimum length (nt) of non-coding segments to keep (default: 200).
 -m / --min-orf-length      Minimum ORF length (nt) that causes a segment to be rejected (default: 250).
 -r / --organism-prefix     Prefix for output FASTA headers (e.g., 'A_thaliana'). Default: 'Nc'.
 -n / --subset-cds-count    If set (>0), randomly selects N coding sequences from the input CDS for output.

Output Filename Convention (within -o directory):
 - Non-coding: <base_name>_nc.<ext>
 - CDS Subset: <base_name>_subset_<N>.<ext>
"""

import sys
import argparse
import os 
import random # Native Python library

MAX_SEG_LEN = 5000  # max segment length to output; longer are split

# --- Utility Functions ---

def read_fasta(path):
    """Reads a FASTA file and returns a dict of sequence names to sequences."""
    seqs = {}
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        # Handle the last sequence in the file
        if name is not None:
            seqs[name] = "".join(chunks).upper()
    return seqs


def parse_gff_cds(path):
    """Return dict: chrom -> list of (start, end) CDS intervals (1-based, inclusive), merged."""
    cds = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, _, ftype, start, end, _, _, _, _ = parts
            if ftype.lower() != "cds":
                continue
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue
            cds.setdefault(chrom, []).append((start, end))
            
    # sort and merge intervals
    for chrom in cds:
        intervals = sorted(cds[chrom])
        merged = []
        if not intervals:
            continue
        
        cur_start, cur_end = intervals[0]
        for s, e in intervals[1:]:
            # Check for overlap or immediate adjacency
            if s <= cur_end + 1:
                cur_end = max(cur_end, e)
            else:
                merged.append((cur_start, cur_end))
                cur_start, cur_end = s, e
        merged.append((cur_start, cur_end))
        cds[chrom] = merged
    return cds


def complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


genetic_code = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def translate_frame(seq, frame):
    """Translate one reading frame on a single strand."""
    aa = []
    n = len(seq)
    i = frame
    while i + 3 <= n:
        codon = seq[i:i+3]
        # Use 'X' for ambiguous bases (like N) or non-canonical codons
        aa.append(genetic_code.get(codon, 'X'))
        i += 3
    return "".join(aa)


def find_orf_longer_than(seq, min_nt):
    """Return True if there is an ATG-initiated ORF >= min_nt in any frame, on either strand."""
    stops = {"TAA", "TAG", "TGA"}

    def scan_one_strand(s):
        n = len(s)
        for frame in range(3):
            i = frame
            in_orf = False
            start_pos = None
            while i + 3 <= n:
                codon = s[i:i+3]
                
                # Treat 'N' or ambiguous codons as an immediate break in the ORF
                if 'N' in codon: 
                    in_orf = False
                    start_pos = None
                    i += 3
                    continue
                
                if not in_orf:
                    if codon == "ATG":
                        in_orf = True
                        start_pos = i
                else:
                    if codon in stops:
                        orf_len = i + 3 - start_pos
                        if orf_len >= min_nt:
                            return True
                        in_orf = False
                        start_pos = None
                i += 3
            # Check for ORF that runs off end
            if in_orf:
                orf_len = n - start_pos
                if orf_len >= min_nt:
                    return True
        return False

    if scan_one_strand(seq):
        return True
    rc = complement(seq)
    if scan_one_strand(rc):
        return True
    return False


def load_protein_kmers(protein_fasta, k=20):
    """Store simple protein k-mers (without '*' or 'X') for exact-match exclusion."""
    prots = read_fasta(protein_fasta)
    kmers = set()
    for p in prots.values():
        p = p.strip().upper()
        if len(p) < k:
            continue
        for i in range(0, len(p) - k + 1):
            kmer = p[i:i+k]
            # Only store clean k-mers from known proteins
            if 'X' not in kmer and '*' not in kmer and len(kmer) == k:
                 kmers.add(kmer)
    return kmers


def peptide_matches_known(seq, protein_kmers, k=20):
    """Crude filter: if any k-mer from translated segment matches known protein k-mers, reject."""
    if not protein_kmers:
        return False
        
    strands = [seq, complement(seq)]
    for s in strands:
        for frame in range(3):
            pep = translate_frame(s, frame)
            if len(pep) < k:
                continue
            for i in range(0, len(pep) - k + 1):
                kmer = pep[i:i+k]
                # Check only clean k-mers against the known set
                if 'X' not in kmer and '*' not in kmer and kmer in protein_kmers:
                    return True
    return False


def get_non_cds_segments(chrom_len, cds_intervals):
    """Return list of (start, end) segments not covered by CDS (1-based inclusive)."""
    segments = []
    if not cds_intervals:
        segments.append((1, chrom_len))
        return segments
    
    # 1. Segment before the first CDS
    first_start, _ = cds_intervals[0]
    if first_start > 1:
        segments.append((1, first_start - 1))
        
    # 2. Segments between CDS features (intergenic or intronic regions)
    for (_, e1), (s2, _) in zip(cds_intervals, cds_intervals[1:]):
        if s2 > e1 + 1:
            segments.append((e1 + 1, s2 - 1))
            
    # 3. Segment after the last CDS
    _, last_end = cds_intervals[-1]
    if last_end < chrom_len:
        segments.append((last_end + 1, chrom_len))
        
    return segments


def split_segment(start, end, max_len):
    """Split long segment into <= max_len chunks."""
    chunks = []
    s = start
    while s <= end:
        e = min(s + max_len - 1, end)
        chunks.append((s, e))
        s = e + 1
    return chunks


# --- Main Execution ---

def parse_args():
    parser = argparse.ArgumentParser(
        description="Builds a high-confidence, conservative non-coding training set from genome + GFF."
    )
    parser.add_argument(
        "-p", "--protein-fasta", required=True,
        help="Protein FASTA file with known coding sequences."
    )
    parser.add_argument(
        "-g", "--genome-fasta", required=True,
        help="Genome FASTA file (nucleotide)."
    )
    parser.add_argument(
        "-a", "--annotation-gff", required=True,
        help="Annotation file in GFF3 format containing CDS features."
    )
    parser.add_argument(
        "-o", "--out-dir", required=True,
        help="Output directory where all result files will be saved."
    )
    parser.add_argument(
        "-l", "--min-noncoding-length", type=int, default=200,
        help="Minimum length (nt) of non-coding segments to keep (default: 200)."
    )
    parser.add_argument(
        "-m", "--min-orf-length", type=int, default=250,
        help="Minimum ORF length (nt) that causes a segment to be rejected as coding-like (default: 250)."
    )
    parser.add_argument(
        "-r", "--organism-prefix", type=str, default="Nc",
        help="Prefix for output FASTA headers (e.g., 'Hsap', 'Bstricta'). Default: 'Nc'."
    )
    parser.add_argument(
        "-n", "--subset-cds-count", type=int, default=0,
        help="If set (>0), randomly selects N coding sequences from the input CDS for output."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Load parameters
    protein_fa = args.protein_fasta # Input CDS/Protein file
    out_dir = args.out_dir          # Output directory
    subset_count = args.subset_cds_count
    min_len = args.min_noncoding_length
    min_orf_nt = args.min_orf_length
    organism_prefix = args.organism_prefix
    
    # --- 1. Define Base Name and Output Paths ---
    # Get the basename from the input CDS file (e.g., 'Arabidopsis_thaliana.fna')
    input_basename = os.path.basename(protein_fa) 
    base_name, ext = os.path.splitext(input_basename)
    
    # Define output names
    NONCODING_OUT_FA = os.path.join(out_dir, f"{base_name}_nc{ext}") 
    
    # Define the CDS subset path (if needed)
    if subset_count > 0:
        SUBSET_OUT_FA = os.path.join(out_dir, f"{base_name}_subset_{subset_count}{ext}")
    
    # --- 2. Create the Output Directory ---
    sys.stderr.write(f"Creating output directory: {out_dir}\n")
    os.makedirs(out_dir, exist_ok=True) 

    # --- 3. Run CDS Subsetting (If requested) ---
    if subset_count > 0:
        sys.stderr.write(f"\n3. Generating random subset of {subset_count} CDS sequences...\n")
        sys.stderr.write(f"   Reading all CDS sequences from: {protein_fa}\n")
        
        all_cds = read_fasta(protein_fa)
        all_ids = list(all_cds.keys())
        total_cds = len(all_ids)
        
        if subset_count > total_cds:
            sys.stderr.write(f"   Warning: Requested {subset_count} > Total {total_cds}. Outputting all sequences.\n")
            selected_ids = all_ids
        else:
            selected_ids = random.sample(all_ids, subset_count)

        # newline="\n": real sequence data fed to make_hexamer_tab/
        # make_logitModel -- a Windows-Python CRLF write would land a stray
        # \r inside sequence content (same bug class fixed elsewhere in
        # SeqMetrics this session).
        with open(SUBSET_OUT_FA, "w", newline="\n") as out_subset:
            for seq_id in selected_ids:
                seq = all_cds[seq_id]
                out_subset.write(f">{seq_id}\n")
                # Write sequence on a single line
                out_subset.write(f"{seq}\n")

        sys.stderr.write(f"   --- Done. Wrote {len(selected_ids)} CDS sequences to {SUBSET_OUT_FA} ---\n")
    # -------------------------------------------------------------------

    sys.stderr.write("4. Reading genome...\n")
    genome = read_fasta(args.genome_fasta)

    sys.stderr.write("5. Parsing and merging GFF CDS intervals...\n")
    cds_intervals = parse_gff_cds(args.annotation_gff)

    sys.stderr.write("6. Loading known protein k-mers (k=20)...\n")
    protein_kmers = load_protein_kmers(args.protein_fasta, k=20)

    try:
        # newline="\n": same reason as SUBSET_OUT_FA above.
        out = open(NONCODING_OUT_FA, "w", newline="\n")
    except Exception as e:
        sys.stderr.write(f"Error opening output file {NONCODING_OUT_FA}: {e}\n")
        sys.exit(1)
        
    noncoding_id = 0

    sys.stderr.write("7. Scanning genome and applying filters...\n")
    for chrom, seq in genome.items():
        # Removed verbose "Processing..." messages as requested.
        
        chrom_len = len(seq)
        cds_ints = cds_intervals.get(chrom, [])
        
        # Step 7a: Get all non-CDS segments
        non_cds = get_non_cds_segments(chrom_len, cds_ints)

        for start, end in non_cds:
            # Step 7b: Split long segments and iterate over chunks
            for s, e in split_segment(start, end, MAX_SEG_LEN):
                segment_len = e - s + 1
                
                # Filter 1: Minimum Length Check
                if segment_len < min_len:
                    continue
                
                subseq = seq[s-1:e] 

                # Filter 2: Long ORF Check (6 frames)
                if find_orf_longer_than(subseq, min_orf_nt):
                    continue

                # Filter 3: Known Protein K-mer Check (6 frames)
                if peptide_matches_known(subseq, protein_kmers, k=20):
                    continue

                # Segment is accepted
                noncoding_id += 1
                
                # Output Header Format: >{Prefix}_nc_{ID}_chr{CHROM}_{START}_{END}
                header = f">{organism_prefix}_nc_{noncoding_id}_chr{chrom}_{s}_{e}"
                
                out.write(header + "\n")
                # Write sequence on a single line
                out.write(subseq + "\n") 

    out.close()
    sys.stderr.write(f"--- Done. Written {noncoding_id} noncoding sequences to {NONCODING_OUT_FA} ---\n")


if __name__ == "__main__":
    main()
    