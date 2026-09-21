# utils/

General-purpose helpers that aren't specific to any one module.

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
- **Output line endings are `\r\n`, on every OS, by design** -- this is
  Python's `csv` module's own default dialect (`lineterminator='\r\n'`),
  not a mistake, and not the same LF-only guarantee every other file in
  this project has after this session's CRLF-bug fixes elsewhere. Fine
  for anything read back by `csv`/pandas (which handle `\r\n`
  transparently), but if you ever feed a merged file's IDs into a
  non-Python exact-match tool (awk, grep -x, another script that reads
  lines naively), be aware the line terminator differs from this
  project's other output files.
