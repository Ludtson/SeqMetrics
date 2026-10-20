#!/usr/bin/env python3
"""
merge_modules.py

Merge one SeqMetrics run's module outputs into one table per input file,
joined on sequence ID. Knows each module's own ID column (they differ:
Seq_ID, Gene_ID, protein_id, ID, codonW's padded 'title'), cleans IDs
(whitespace, codonW padding, stray \\r), and prefixes every feature column
with its module name (composition.Isoelectric_Point, basic.GC_pct, ...) so
the several length/GC columns from different modules never collide and
every value stays traceable to the module that produced it.

Outer join: a sequence missing from one module keeps its row, with empty
cells for that module, and is counted in the coverage report -- nothing
silently disappears.

Layouts (same as run_features.py):
  single-file  OUT/<module>.tsv           -> OUT/merged.tsv (+ merged.coverage.tsv)
  batch        OUT/<module>/<stem>.tsv    -> OUT/merged/<stem>.tsv (+ <stem>.coverage.tsv)
--stack (batch only) also writes OUT/merged/ALL.tsv: every stem stacked, with a
'set' column; columns are the union, and tm_domain Residues_* columns absent
from one batch are filled with 0 (DeepTMHMM2 only reports labels it saw --
see docs/install_tm_domain.md), every other missing value stays empty.

Usage:
  python utils/merge_modules.py --run-dir outputs/
  python utils/merge_modules.py --run-dir outputs/ --modules composition,disorder --stack
"""

import argparse
import csv
import sys
from pathlib import Path

# each module's ID column, as written by that module (first column is the fallback)
MODULE_ID = {
    "basic": "Seq_ID", "tail_hydrophobicity": "Seq_ID", "composition": "Gene_ID",
    "aggregation": "Seq_ID", "disorder": "protein_id", "tm_domain": "Seq_ID",
    "localization": "ID", "codon_usage": "title", "coding_potential": "Gene_ID",
}
# secondary ID-like columns that are not features
DROP = {"coding_potential": {"ID"}}
KNOWN = list(MODULE_ID)
ZERO_FILL_PREFIX = "tm_domain.Residues_"


def read_module(path, module, warn):
    """Return (feature_columns, {id: {col: value}}, ids_in_order) for one module TSV."""
    with open(path, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows:
        return [], {}, []
    header = [h.strip().replace("\r", "") for h in rows[0]]
    key = MODULE_ID.get(module, header[0])
    if key not in header:
        warn(f"{path}: expected ID column '{key}' not found; using first column '{header[0]}'")
        key = header[0]
    k = header.index(key)
    keep = [i for i, h in enumerate(header) if i != k and h and h not in DROP.get(module, set())]
    cols = [f"{module}.{header[i]}" for i in keep]
    data, order, dups = {}, [], 0
    for r in rows[1:]:
        r = [x.strip().replace("\r", "") for x in r]
        if not r or not any(r) or len(r) <= k or not r[k]:
            continue
        sid = r[k]
        if sid in data:
            dups += 1
            continue
        data[sid] = {c: (r[i] if i < len(r) else "") for c, i in zip(cols, keep)}
        order.append(sid)
    if dups:
        warn(f"{path}: {dups} duplicate IDs, kept the first occurrence of each")
    return cols, data, order


def merge_set(paths, warn):
    """paths: {module: Path}. Returns (columns, rows as dicts in first-seen order, coverage rows)."""
    columns, tables, order, seen = [], {}, [], set()
    for module in [m for m in KNOWN if m in paths] + [m for m in paths if m not in KNOWN]:
        cols, data, ids = read_module(paths[module], module, warn)
        columns += cols
        tables[module] = (cols, data)
        for sid in ids:
            if sid not in seen:
                seen.add(sid); order.append(sid)
    merged = []
    for sid in order:
        row = {"Seq_ID": sid}
        for cols, data in tables.values():
            # explicit blanks for a sequence a module skipped, so --stack's 0-fill only
            # ever touches tm_domain columns that a whole batch never produced
            row.update(data.get(sid) or {c: "" for c in cols})
        merged.append(row)
    coverage = [{"module": m, "rows": len(d), "missing_vs_union": len(order) - len(d)} for m, (_, d) in tables.items()]
    return columns, merged, coverage


def write_tsv(path, columns, rows, fill=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="\n") as f:  # LF on every OS, like the module outputs
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        for r in rows:
            w.writerow([r[c] if c in r else (fill(c, r) if fill else "") for c in columns])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, type=Path, help="run_features.py --out-dir of one run")
    ap.add_argument("--modules", help="Comma-separated subset (default: every module output found)")
    ap.add_argument("--stack", action="store_true", help="Batch layout: also write merged/ALL.tsv with a 'set' column")
    args = ap.parse_args()
    run = args.run_dir
    wanted = set(args.modules.split(",")) if args.modules else None
    warn = lambda m: print(f"[merge_modules] WARN {m}", file=sys.stderr)

    single = {p.stem: p for p in run.glob("*.tsv") if p.stem in MODULE_ID or (wanted and p.stem in wanted)}
    batch_dirs = [d for d in run.iterdir() if d.is_dir() and d.name in MODULE_ID]
    if wanted:
        single = {m: p for m, p in single.items() if m in wanted}
        batch_dirs = [d for d in batch_dirs if d.name in wanted]

    if single and not batch_dirs:
        cols, rows, cov = merge_set(single, warn)
        write_tsv(run / "merged.tsv", ["Seq_ID"] + cols, rows)
        write_tsv(run / "merged.coverage.tsv", ["module", "rows", "missing_vs_union"], cov)
        print(f"[merge_modules] {len(rows)} sequences x {len(cols)} features from {len(single)} modules -> {run / 'merged.tsv'}")
        return
    if not batch_dirs:
        sys.exit(f"No module outputs found in {run} (expected <module>.tsv or <module>/<stem>.tsv)")

    stems = sorted({p.stem for d in batch_dirs for p in d.glob("*.tsv")})
    stacked, all_cols = [], []
    for stem in stems:
        paths = {d.name: d / f"{stem}.tsv" for d in batch_dirs if (d / f"{stem}.tsv").exists()}
        cols, rows, cov = merge_set(paths, warn)
        write_tsv(run / "merged" / f"{stem}.tsv", ["Seq_ID"] + cols, rows)
        write_tsv(run / "merged" / f"{stem}.coverage.tsv", ["module", "rows", "missing_vs_union"], cov)
        print(f"[merge_modules] {stem}: {len(rows)} sequences x {len(cols)} features from {len(paths)} modules")
        if args.stack:
            all_cols += [c for c in cols if c not in all_cols]
            stacked += [{"set": stem, **r} for r in rows]
    if args.stack:
        # 0 only for a topology label this row's batch never reported AND only if the row
        # has a tm_domain result at all; a sequence tm_domain skipped stays empty everywhere
        has_tm = lambda r: r.get("tm_domain.N_TMRs", "") != ""
        write_tsv(run / "merged" / "ALL.tsv", ["set", "Seq_ID"] + all_cols, stacked,
                  fill=lambda c, r: "0" if c.startswith(ZERO_FILL_PREFIX) and has_tm(r) else "")
        print(f"[merge_modules] stacked {len(stacked)} rows x {len(all_cols)} features -> {run / 'merged' / 'ALL.tsv'}")


if __name__ == "__main__":
    main()
