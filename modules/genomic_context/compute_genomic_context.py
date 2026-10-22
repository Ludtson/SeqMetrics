#!/usr/bin/env python3
"""
Genomic-context features for each sequence of a CDS FASTA, from the
genome annotation (GFF3) and the genome FASTA it was made on.

Each FASTA ID (first whitespace token), or the [protein_id=...] tag of
NCBI cds_from_genomic headers, is looked up in the GFF as a transcript ID
or Name, or as a CDS protein_id / Protein_Accession / Name, each also
tried with or without a '.p' suffix and an 'rna-' or 'transcript:' prefix. Only CDS coordinates are used; UTR annotation is
too uneven between genomes to compare.

Columns:
  cds_exons          number of CDS segments
  cds_span           genomic span from the first to the last coding base
  cds_intron_total   total length of the introns between CDS segments
  kozak_score        sum of log10 base frequencies over the 12-nt window
                     -6..+6 around the start codon (6 nt upstream, ATG,
                     3 nt downstream). Frequencies per position come from
                     all sequences in the input with a valid window
                     (pseudocount 1), so scores are relative to that set.
                     Empty if the CDS does not start with ATG or the window
                     runs off the sequence.
  overlaps_gene      1 if the CDS span overlaps another annotated gene on
                     either strand, else 0
  nearest_gene_dist  distance (bp) from the CDS span to the nearest other
                     gene on the same sequence; 0 if overlapping, empty if
                     there is none
  ssr_count, ssr_bp, ssr_frac
                     perfect tandem repeats in the CDS, motif 1-6 nt, at
                     least 10/6/5/5/5/5 copies for motif length 1/2/3/4/5/6
                     (MISA defaults); ssr_bp merges overlapping repeats
Rows for IDs that cannot be placed have empty feature columns.

  kozak_logodds      only with --kozak-background: sum over the same window of
                     log10 f(base) - log10 f_background(base), where f is the
                     table behind kozak_score and f_background a saved table,
                     e.g. from the species' igORFs. Base composition cancels,
                     so scores compare across genomes.

Sequences without an annotation (e.g. IGOR igORFs) can be located from their IDs with --id-coords. Kozak
frequencies can be saved (--kozak-table-out) and reused (--kozak-table), e.g. to score igORFs against the
species' genes.

The GFF and genome must use the same sequence names. If they don't (for
example GenBank names in one file and RefSeq names in the other), pass
--seqid-map, a tab-separated file whose first two columns are the GFF name
and the genome name. Without it the run stops with an error when most
linked sequences name a sequence absent from the genome.
"""
import argparse
import bisect
import csv
import gzip
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
MIN_COPIES = {1: 10, 2: 6, 3: 5, 4: 5, 5: 5, 6: 5}
SSR_RE = [re.compile(r"([ACGT]{%d})\1{%d,}" % (k, MIN_COPIES[k] - 1)) for k in range(1, 7)]
COLUMNS = ["Seq_ID", "cds_exons", "cds_span", "cds_intron_total", "kozak_score", "overlaps_gene",
           "nearest_gene_dist", "ssr_count", "ssr_bp", "ssr_frac"]


def _open(path):
    path = str(path)
    return gzip.open(path, "rt", encoding="utf-8", errors="ignore") if path.endswith(".gz") \
        else open(path, encoding="utf-8", errors="ignore")


def fasta_iter(path, headers=False):
    """Yield (ID, sequence), or (ID, full header) with headers=True."""
    name, head, buf = None, "", []
    with _open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, head if headers else "".join(buf)
                head = line[1:].strip()
                name, buf = head.split()[0] if head else "", []
            elif not headers:
                buf.append(line.strip().upper())
        if name is not None:
            yield name, head if headers else "".join(buf)


def _attrs(field):
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def read_gff(path, seqid_map):
    """Return (key -> transcript ID, transcript ID -> CDS segments,
    transcript ID -> parent gene ID, sequence -> sorted genes)."""
    tx_keys, cds, parent, genes = {}, defaultdict(list), {}, defaultdict(list)
    with _open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            seq = seqid_map.get(c[0], c[0])
            a = _attrs(c[8])
            start, end = int(c[3]), int(c[4])
            if c[2] == "gene":
                genes[seq].append((start, end, a.get("ID", "")))
            elif c[2] in ("mRNA", "transcript"):
                tid = a.get("ID")
                if tid:
                    keys = (tid, a.get("Name"), re.sub(r"^(rna-|transcript:)", "", tid))
                    tx_keys.setdefault(tid, set()).update(k for k in keys if k)
                    if a.get("Parent"):
                        parent[tid] = a["Parent"].split(",")[0]
            elif c[2] == "CDS":
                pid = a.get("protein_id") or a.get("Protein_Accession") or a.get("Name")
                for par in a.get("Parent", "").split(","):
                    if not par:
                        continue
                    cds[par].append((seq, start, end, c[6]))
                    if pid:
                        tx_keys.setdefault(par, set()).add(pid)
    key2tx = {}
    for tid, keys in tx_keys.items():
        for k in keys:
            for v in (k, k + ".p", re.sub(r"\.p$", "", k)):
                key2tx.setdefault(v, tid)
    for s in genes:
        genes[s].sort()
    return key2tx, cds, parent, genes


def read_seqid_map(path):
    m = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 2 and not line.startswith("#"):
                m[c[0]] = c[1]
    return m


def read_kozak_table(path):
    """Counts per window position (12 rows of A, C, G, T) as written by --kozak-table-out."""
    with open(path, encoding="utf-8") as f:
        rows = [{b: float(r[b]) for b in "ACGT"} for r in csv.DictReader(f, delimiter="\t")]
    if len(rows) != 12:
        sys.exit(f"[genomic_context] {path}: expected 12 positions, found {len(rows)}")
    return rows


def primitive(motif):
    """False if the motif is a repeat of a shorter unit, e.g. 'ATAT'."""
    k = len(motif)
    return not any(k % d == 0 and motif == motif[:d] * (k // d) for d in range(1, k))


def ssr(seq):
    hits = sorted((m.start(), m.end()) for rx in SSR_RE for m in rx.finditer(seq) if primitive(m.group(1)))
    merged = []
    for s, e in hits:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    bp = sum(e - s for s, e in merged)
    return len(hits), bp, round(bp / len(seq), 6) if seq else 0.0


def neighbours(genes, lo, hi, own_id, envelope=True):
    """(overlaps another gene, distance to the nearest other gene) for the span lo..hi. The sequence's own gene is
    skipped by ID, or, without an ID and with envelope=True, as the gene that spans the whole CDS."""
    if not genes:
        return None, None
    ov, near = 0, None
    j = bisect.bisect_left(genes, (lo, 0, ""))
    for s, e, gid in genes[max(0, j - 50):j + 50]:
        if (gid == own_id) if own_id else (envelope and s <= lo and e >= hi):
            continue
        d = 0 if (s <= hi and e >= lo) else (s - hi if s > hi else lo - e)
        ov = ov or int(d == 0)
        near = d if near is None else min(near, d)
    return ov, near


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("nt", type=Path, help="CDS FASTA (its IDs select the genes)")
    ap.add_argument("--gff", required=True, type=Path, help="GFF3 annotation (.gz accepted)")
    ap.add_argument("--genome", required=True, type=Path, help="Genome FASTA (.gz accepted)")
    ap.add_argument("--seqid-map", type=Path, help="GFF name -> genome name, tab-separated")
    ap.add_argument("--id-coords", action="store_true",
                    help="take each sequence's location from its ID instead of the GFF: '|'-separated fields 2-5 = "
                         "sequence, start, end, strand, 0-based half-open (BED) coordinates including the stop codon, "
                         "as in IGOR igORF IDs (e.g. Athaliana|Chr3|5931873|5931990|-|24). The GFF is still used for "
                         "neighbouring genes.")
    ap.add_argument("--kozak-table", type=Path,
                    help="score Kozak windows against this frequency table (from --kozak-table-out) instead of one "
                         "estimated from the input, e.g. igORFs against the species' genes")
    ap.add_argument("--kozak-table-out", type=Path, help="write the Kozak base counts used for scoring (12 x ACGT)")
    ap.add_argument("--kozak-background", type=Path,
                    help="background Kozak table (from --kozak-table-out, e.g. of the species' igORFs); adds the "
                         "kozak_logodds column")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    seqid_map = read_seqid_map(args.seqid_map) if args.seqid_map else {}
    key2tx, cds, parent, genes = read_gff(args.gff, seqid_map)
    genome = dict(fasta_iter(args.genome))
    ids, placed_in = [], []  # placed_in: (Seq_ID, sequence, strand, [(start, end) 1-based], own gene ID or None)
    for sid, head in fasta_iter(args.nt, headers=True):
        ids.append(sid)
        if args.id_coords:
            f = sid.split("|")
            if len(f) >= 5 and f[2].isdigit() and f[3].isdigit() and f[4] in "+-":
                placed_in.append((sid, f[1], f[4], [(int(f[2]) + 1, int(f[3]))], None))
            continue
        pid = re.search(r"\[protein_id=([^\]\s]+)\]", head)  # NCBI cds_from_genomic headers
        tid = next((key2tx[k] for k in (sid, pid and pid.group(1)) if k in key2tx), None)
        if tid in cds:
            chrom, strand = cds[tid][0][0], cds[tid][0][3]
            segs = sorted((s[1], s[2]) for s in cds[tid] if s[0] == chrom)
            placed_in.append((sid, chrom, strand, segs, parent.get(tid)))
    if ids and not placed_in:
        sys.exit(f"[genomic_context] none of the {len(ids)} FASTA IDs (e.g. {ids[:3]}) "
                 + ("has sequence|start|end|strand in fields 2-5" if args.id_coords else
                    f"matches a transcript ID/Name or CDS protein ID in {args.gff}"))
    absent = [p for p in placed_in if p[1] not in genome]
    if placed_in and len(absent) > len(placed_in) / 2:
        names = sorted({p[1] for p in absent})[:3]
        sys.exit(f"[genomic_context] {len(absent)} of {len(placed_in)} sequences lie on sequences missing from the "
                 f"genome FASTA (e.g. {names}; genome names e.g. {list(genome)[:3]}). "
                 f"Use --seqid-map if the GFF and genome name sequences differently.")

    rows, windows, n_placed = {sid: {} for sid in ids}, {}, 0
    for sid, chrom, strand, segs, own in placed_in:
        if chrom not in genome:
            continue
        lo, hi = segs[0][0], segs[-1][1]
        g = genome[chrom]
        seq = "".join(g[s - 1:e] for s, e in segs)
        if strand == "-":
            seq = seq.translate(COMP)[::-1]
            # start codon at hi-2..hi (1-based); window hi-5..hi+6, reverse-complemented
            win = g[hi - 6:hi + 6].translate(COMP)[::-1] if hi + 6 <= len(g) else ""
        else:
            # start codon at lo..lo+2 (1-based); window lo-6..lo+5
            win = g[lo - 7:lo + 5] if lo >= 7 else ""
        if len(win) == 12 and win[6:9] == "ATG" and all(b in "ACGT" for b in win):
            windows[sid] = win
        ov, near = neighbours(genes.get(chrom, []), lo, hi, own, envelope=not args.id_coords)
        n_ssr, ssr_bp, ssr_frac = ssr(seq)
        rows[sid] = {"cds_exons": len(segs), "cds_span": hi - lo + 1,
                     "cds_intron_total": sum(segs[i + 1][0] - segs[i][1] - 1 for i in range(len(segs) - 1)),
                     "overlaps_gene": ov, "nearest_gene_dist": near,
                     "ssr_count": n_ssr, "ssr_bp": ssr_bp, "ssr_frac": ssr_frac}
        n_placed += 1

    if args.kozak_table:
        counts = read_kozak_table(args.kozak_table)
    else:
        counts = [defaultdict(lambda: 1) for _ in range(12)]
        for w in windows.values():
            for i, b in enumerate(w):
                counts[i][b] += 1
    if args.kozak_table_out:
        with open(args.kozak_table_out, "w", newline="") as o:
            o.write("position\tA\tC\tG\tT\n")
            for i, p in enumerate(list(range(-6, 0)) + list(range(1, 7))):
                o.write(f"{p}\t" + "\t".join(str(counts[i][b]) for b in "ACGT") + "\n")
    total = [sum(counts[i][b] for b in "ACGT") for i in range(12)]
    bg = read_kozak_table(args.kozak_background) if args.kozak_background else None
    bg_total = [sum(bg[i][b] for b in "ACGT") for i in range(12)] if bg else None
    for sid, w in windows.items():
        lp = [math.log10(counts[i][b] / total[i]) for i, b in enumerate(w)]
        rows[sid]["kozak_score"] = round(sum(lp), 4)
        if bg:
            rows[sid]["kozak_logodds"] = round(sum(lp) - sum(math.log10(bg[i][b] / bg_total[i])
                                                                for i, b in enumerate(w)), 4)

    columns = COLUMNS + (["kozak_logodds"] if bg else [])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for sid in ids:
            w.writerow({"Seq_ID": sid, **{k: ("" if v is None else v) for k, v in rows[sid].items()}})

    print(f"[genomic_context] {len(ids)} sequences: {len(placed_in)} located "
          f"({'from IDs' if args.id_coords else 'in the GFF'}), {n_placed} placed on the genome, "
          f"{len(windows)} with a Kozak window -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
