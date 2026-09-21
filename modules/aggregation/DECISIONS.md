# HCA + TANGO pipeline — decisions and validation

Source: `~/scripts/run_hca_tango.py` (WSL, user `ludtson`), pulled and reviewed
2026-09-03. Rewritten version: `bin/run_hca_tango.py` (0.1.0, pre-release,
unshared). Original preserved unmodified alongside as
`bin/run_hca_tango_original.py` for reference/diff.

## Why this needed rework, not just a light pass

The user asked for this to be made "more robust like the IUPred script."
Reviewing it surfaced one bug so severe it likely made every past run of
this pipeline on real multi-species FASTA produce empty or near-empty
output silently, with no error anywhere in the run. Everything below was
verified empirically before being asserted (same discipline as pepstats,
IUPred, CodonW this project) — including retracting two initial suspicions
that turned out, on direct testing, to be wrong.

## Confirmed bugs and fixes

### 1. Header/ID mismatch silently dropped nearly all sequences (critical)

`fasta_iter()` set `name = line[1:]` — the **entire** raw header text after
`>`, with no whitespace split. Meanwhile `hcatk`'s own segmentation output
normalizes each ID to its first whitespace-delimited token. The join
`if name in hca` therefore compared a full messy header against a clean
token and failed for any FASTA with descriptive headers.

Confirmed by direct test: given `>gene1 | Symbols: FOO | desc here`,
`hcatk segment -m domain` wrote back `>gene1 78 0.003 0.013` — first token
only. The original script's `name` for that record would still have been
the full `"gene1 | Symbols: FOO | desc here"` string, which never equals
`"gene1"`.

All three species' real protein FASTA files carry exactly this kind of
multi-field header (established repeatedly earlier in this project's
pepstats/IUPred/CodonW work), so this was not a hypothetical edge case.

### 2. Pipe-chain-immediately-after-ID headers break `hcatk`'s own parsing too

For headers with no whitespace before the first `|` (this project's O.
sativa convention, e.g. `ID|chromosome:IRGSP-1.0:...`), `hcatk`'s
whitespace-only tokenization can't separate the ID from the pipe-chain
either — its own output key would be the whole chain, not the gene ID.
Identical bug already found and fixed in
`codonw_extracted/bin/remap_seq_ids.py` (`original_id_from_header()`) this
project.

Because **two independent third-party binaries** (`hcatk`, `tango`) each do
their own uncontrolled header parsing, patching the *output* side per-tool
would mean re-deriving and maintaining two separate ID-recovery hacks
forever, one per species convention. Instead, the fix mirrors the CodonW
pattern: remap every sequence to a clean `seqN` ID **before** either tool
ever sees the FASTA, and restore original IDs only in the final TSV. Both
tools then agree on IDs by construction, for any species' header style,
with no per-species logic needed. Directly addresses the user's "should be
usable for other species" goal stated earlier in this project.

### 3. Unescaped ID/sequence spliced into a `shell=True` string

```python
cmd = f'tango {name} ct="N" nt="N" ph="7.0" te="298" io="0.02" seq="{seq}"'
subprocess.run(cmd, shell=True, ...)
```

Any header containing spaces, `|`, or quotes (again: most real headers
here) would be reinterpreted by the shell — e.g. `|` as an actual pipe —
independently of bug #1/#2, corrupting or breaking the TANGO invocation.

Fixed by building an argv list and calling `subprocess.run` without
`shell=True`. Verified empirically that list-form invocation produces
identical TANGO stdout to the shell-string form for the same inputs — no
behavior change, only the parsing risk removed. The seqN remap (fix #2)
also removes the metacharacter risk at the source; the argv-list change is
defense in depth on top of that, not a substitute for it.

### 4. One failed sequence killed the entire batch

```python
if result.returncode != 0:
    raise RuntimeError(f"TANGO failed for {name}")
```

Raised inside a `ProcessPoolExecutor` worker, this propagates out of
`fut.result()` and stops the whole run — a single bad sequence in a
multi-hour, thousands-of-sequences job loses all completed work. Now caught
per-sequence: logged to `<prefix>.tango_logs/failures.log` (seqN, original
ID, reason) and skipped; the run finishes and reports a failure count in
`<prefix>.run_summary.txt` instead of dying partway through. Same
"per-item error isolation, not batch-killing" principle used in the IUPred
rebuild.

### 5. HCA-file parsing counted `domain`/`cluster` lines as sequences

`-m domain` output has, per sequence, a `>id length pvalue score` line
followed by `domain ...` and `cluster ...` annotation lines — the latter
two also have >=4 whitespace-split fields. The original (and my first
draft) unconditionally applied `if len(parts) >= 4: hca[parts[0]...] = ...`
without checking for the `>` prefix, so those annotation lines got
miscounted as pseudo-sequences keyed `"domain"`/`"cluster"`. Caught by the
smoke test itself: a 3-sequence input reported "Parsed HCA scores for 5
sequences" before this fix. Harmless for the real join (those keys never
match a real `seqN` ID) but wrong bookkeeping, closed while already in this
code. Fix: skip any line not starting with `>` before parsing.

### 6. Empty / stop-codon-only sequences

CDS-derived protein FASTA commonly carries a trailing `*` (stop codon
symbol). Not handled before; now stripped, and a sequence that is empty
after stripping is skipped with a logged reason rather than fed to
TANGO/HCA as-is.

## Suspicions raised, then retracted after direct verification

- **"AGG"/"AMYLO" token parsing might not match real TANGO output.**
  Suspected from an earlier, differently-parameterized standalone test
  whose output looked like `... HELAGG 0 ...` with no isolated `AGG` token.
  Re-tested with this script's *exact* parameters
  (`ct=N nt=N ph=7.0 te=298 io=0.02`): real output is
  `AGG 424.614 AMYLO 5773.45 TURN 60.4393 HELIX 17.4114 HELAGG 0 BETA
  188.046` — `AGG` and `AMYLO` ARE separate, correctly-matched tokens.
  The original parsing logic for these two fields was correct all along.
  Retracted before ever asserting it to the user.

- **`parts[3]` for the HCA score might be the wrong column** (looked like
  it could be a p-value instead, since the score-line comment in `hcatk`'s
  own output only documents 3 fields but real output has 4). Checked
  against `pyHCA` source directly
  (`pyHCA/core/annotateHCA.py:986`):
  `outf.write(">{} {} {:.3f} {:.3f}\n".format(prot, len(sequence), pvalue, score))`
  — field order is `id, length, pvalue, score`, so `score` is the last
  (4th) field, i.e. `parts[3]`. The original indexing was correct; only the
  header comment in `hcatk`'s own output is stale/incomplete documentation.
  Retracted before asserting a bug here.

## Additional improvement (not a bug fix, an enrichment)

TANGO's stdout already contains `TURN`, `HELIX`, `HELAGG`, `BETA` alongside
`AGG`/`AMYLO`, at zero extra runtime cost — the original script parsed and
discarded them. All six are now captured as separate output columns. Same
"use everything the tool already gives you" principle applied to the
IUPred rebuild.

## Validation performed

Smoke-tested against a synthetic 3-record FASTA deliberately covering all
three real header conventions seen in this project:
- A. thaliana-style: `>AT1G51370.2 | Symbols:  | ... | chr1:...`
- O. sativa-style (pipe with no preceding space): `>OSJAP01G00010.5|chromosome:...`
- B. rapa-style (clean, no extra fields): `>BraA01g000010.3C`

Result: all 3 sequences correctly matched between HCA and TANGO and
produced full rows in the final TSV, with original (not seqN) IDs restored
— under the original script, the first two would have matched 0 times.

Also tested:
- Re-running on the same input archives the prior output
  (`<prefix>_prior_run_<timestamp>/`) instead of overwriting it.
- A batch containing one stop-codon-only (empty after stripping) sequence
  alongside two valid ones: the bad sequence is logged to
  `failures.log` and skipped; the two valid sequences still produce
  correct, complete rows — the run does not crash.
- List-form (no `shell=True`) TANGO invocation produces stdout identical
  in content/format to the original shell-string form for the same inputs.

Not yet done (deferred per the same "trust historical values for the
remaining ~14 untouched columns" decision the user made for the historical
dataset): a full-scale re-run of all three species' real FASTA files
through this corrected pipeline, and a comparison against whatever
historical HCA/TANGO values exist in the ~20-column historical dataset.
That full run + validation is the natural next step once prioritized.
