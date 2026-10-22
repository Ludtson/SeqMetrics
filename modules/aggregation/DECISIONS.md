# Aggregation module (HCA + TANGO): design notes

## Sequence IDs

`hcatk` keys its output on the first whitespace-separated token of each
header, and does not split headers such as `ID|chromosome:...` that have no
whitespace before the first `|`. TANGO receives the ID on its command line.
To keep both tools consistent for any header style, every sequence is
remapped to a short ID (`seqN`) before either tool runs. Original IDs
(first whitespace token, cut at the first `|`) are restored only in the
final TSV.

## Calling TANGO

TANGO is called with an argument list, without a shell, so IDs and
sequences are never interpreted as shell syntax. Output is identical to the
shell-string form.

## Failures

A sequence for which TANGO fails is logged to
`<prefix>.tango_logs/failures.log` (seqN, original ID, reason) and skipped.
The run continues, and the failure count is reported in
`<prefix>.run_summary.txt`.

## Parsing

- **HCA score:** for `hcatk segment -m domain`, each sequence line is
  `>id length pvalue score` (pyHCA `core/annotateHCA.py`), so the score is
  the fourth field. Only lines starting with `>` are read; `domain` and
  `cluster` annotation lines are skipped.
- **TANGO:** with `ct=N nt=N ph=7.0 te=298 io=0.02`, TANGO prints
  `AGG`, `AMYLO`, `TURN`, `HELIX`, `HELAGG` and `BETA` totals. All six are
  reported.

## Input cleaning

A trailing `*` (stop symbol) is removed. Sequences that are empty after
this are skipped and logged.

## Re-runs

Existing output for the same input prefix is moved to
`<prefix>_prior_run_<timestamp>/` rather than overwritten.

## Tests

- Three headers in different styles (`ID | description`, `ID|chromosome:...`,
  plain `ID`): all matched between HCA and TANGO, with original IDs in the
  output.
- A batch with one sequence that is empty after cleaning: that sequence is
  logged and skipped; the others are complete.
- Re-running archives the previous output.
