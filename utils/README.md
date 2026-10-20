# utils/

General-purpose helpers that aren't specific to any one module.

## merge_modules.py

Merges one run's module outputs into **one table per input file**, in a
single step. Use this rather than chaining `combine_datasets.py` once per
module.

```
python utils/merge_modules.py --run-dir outputs/            # single-file or batch layout
python utils/merge_modules.py --run-dir outputs/ --stack    # batch: also merged/ALL.tsv
```

- **Knows each module's ID column.** They differ: `Seq_ID`, `Gene_ID`
  (composition, coding_potential), `protein_id` (disorder), `ID`
  (localization), `title` (codon_usage). All are joined into one `Seq_ID`.
- **Prefixes every feature column with its module**, for example
  `composition.Isoelectric_Point` or `basic.GC_pct`. Several modules
  report their own length or GC, so nothing collides, and each value
  shows which module produced it.
- **Outer join.** A sequence missing from one module keeps its row with
  empty cells. `merged.coverage.tsv` (or `merged/<stem>.coverage.tsv`)
  lists rows per module and how many are missing. On a real
  *A. thaliana* run, disorder lacked 2 sequences and coding_potential 1.
- **`--stack`** (batch only) also writes `merged/ALL.tsv`, with every input
  file stacked and a `set` column. Columns are the union across files.
  tm_domain `Residues_*` columns that a batch never produced are filled
  with 0 (DeepTMHMM2 only reports labels it saw); all other gaps stay empty.
- IDs are stripped of whitespace. Output is TSV with LF line endings.
  Only the Python standard library is used.

## combine_datasets.py

Left-joins two tabular files (CSV or TSV, auto-detected) on a key column
-- the tool for combining separate module output tables (e.g.
`composition.tsv` + `aggregation.tsv`) into one enriched table, or for
merging a `tm_domain` output batch that has extra `Residues_TMhelix`-style
columns (see docs/install_tm_domain.md's note on why those columns are
dynamic) against one that doesn't.

```
python utils/combine_datasets.py \
    -b composition.tsv -e aggregation.tsv \
    -k Gene_ID:Seq_ID -o merged.tsv
```

**Two things worth knowing before using it, both confirmed by direct
test, not assumed:**

- **Output is always TSV, regardless of input delimiter** (fixed this
  session -- it previously silently defaulted to comma regardless of
  whether the inputs were TSV, which would have produced a comma file
  from two genuinely tab-delimited inputs; now hardcoded to tab-out,
  matching this project's own convention, confirmed by test even when
  the base file is CSV).
- **Output line endings are LF (`\n`) on every OS**, like every other
  output in this repo. Python's `csv` writer defaults to `\r\n` on every
  OS, even when the file is opened with `newline="\n"`, so the writer
  sets `lineterminator="\n"` explicitly. The same default had put `\r` on
  the basic, tail_hydrophobicity, disorder, tm_domain and localization
  outputs until it was fixed in each writer.
