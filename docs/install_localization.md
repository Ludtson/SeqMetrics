# Install: localization (LOCALIZER)

LOCALIZER predicts chloroplast/mitochondrial transit peptides
(cTP/mTP) and nuclear localization signals (NLS) in plant protein
sequences. Plant-specific — no general-organism alternative (e.g.
TargetP) is built into this repo.

GPL-3.0 licensed, no registration needed. This repo ships wrapper
scripts only (`run_localizer.py`, `localizer_bin.py`), never LOCALIZER
itself, per its own redistribution policy.

## Dependencies

1. **LOCALIZER itself.** Clone it as a sibling of this repo, not
   somewhere ad hoc — every `--module-ref localization=...` example
   below, and every project script that calls `run_features.py`, uses
   this exact path:
   ```
   cd ..                      # next to SeqMetrics, not inside it
   git clone https://github.com/JanaSperschneider/LOCALIZER.git
   cd LOCALIZER/Scripts
   unzip weka-3-6-12.zip
   ```
   Result: `<SeqMetrics's parent dir>/LOCALIZER/Scripts/LOCALIZER.py`.
2. **Java** — required for `java -cp weka.jar
   weka.classifiers.functions.SMO`, called directly by
   `localization.py`.
3. **Perl** — for NLStradamus (`nlstradamus.pl`, bundled with
   LOCALIZER).
4. **Python 3** — runs `LOCALIZER.py`.
5. **EMBOSS's `pepstats`** — LOCALIZER does not look for `pepstats` on
   PATH. It hardcodes `SCRIPT_PATH + '/EMBOSS-6.5.7/emboss/pepstats'`
   and refuses to run unless that exact path exists. Fix, without a
   full EMBOSS rebuild:
   ```
   mkdir -p <LOCALIZER>/Scripts/EMBOSS-6.5.7/emboss
   ln -s $(which pepstats) <LOCALIZER>/Scripts/EMBOSS-6.5.7/emboss/pepstats
   ```
6. **Patch for parallel runs** — needed when running several
   species at once (`--jobs` > 1), and for chunked runs with
   `--chunk-jobs` > 1 (the default is 4, see "Chunked runs" below).
   Without it, use `--chunk-jobs 1`. Concurrent Java processes collide
   on Java's shared performance-data folder in `/tmp`, and the losing
   process prints a warning into WEKA's output. `parse_weka_output()`
   in `localization.py` skips exactly 5 header lines, so that extra
   line pushes the `inst#` header into the data and the run fails with
   `ValueError: invalid literal for int() with base 10: 'inst#'`. It
   is random, so a different species can fail each run. Fix:
   ```
   cd <LOCALIZER>/Scripts
   cp localization.py localization.py.bak
   sed -i "s/\['java', '-cp'/['java', '-XX:-UsePerfData', '-cp'/g" localization.py
   python3 - <<'EOF'
   p = "localization.py"
   s = open(p).read()
   old = "        content = content[5:]\n\n        for line in content:\n            if line.strip():"
   new = "\n        for line in content:\n            if line.strip() and line.split()[0].isdigit():"
   assert s.count(old) == 1, s.count(old)
   open(p, "w").write(s.replace(old, new))
   print("parser patched")
   EOF
   ```
   The first edit stops Java writing the shared file. The second makes
   the parser read only lines that start with an instance number.

## Environment

EMBOSS, Perl, Java, and Python all live in the `em_boss` env (shared
with `composition`), not a separate env:
```
conda create -n em_boss -c bioconda -c conda-forge emboss   # skip if already done for composition
conda install -n em_boss -c conda-forge perl openjdk python=3.11
```

## Usage

```
python run_localizer.py proteins.faa \
    --localizer-script ../LOCALIZER/Scripts/LOCALIZER.py \
    --mode plant \
    --output OUTDIR \
    --localizer-parser localizer_bin.py
```

Same path via `run_features.py`: `--module-ref localization=../LOCALIZER/Scripts/LOCALIZER.py`
(or the absolute equivalent, `<SeqMetrics's parent dir>/LOCALIZER/Scripts/LOCALIZER.py`).

### Chunked runs: progress, resume, parallel

LOCALIZER itself takes the whole FASTA in one call and writes nothing
until it finishes. On a full proteome (30–50k proteins) that meant
hours with no sign of progress, and a crash lost everything. So
`run_localizer.py` now splits the input into chunks and runs LOCALIZER
on each:

```
python run_localizer.py proteins.faa --localizer-script ... --output OUTDIR \
    --localizer-parser localizer_bin.py --chunk-size 2000 --chunk-jobs 4
```

- **Progress:** one line per finished chunk on stderr, e.g.
  `PROGRESS 5/16 chunks (10000/30401 seqs) | elapsed 22.0 min | ETA 48.0 min`.
- **Resume:** chunks live in `OUTDIR/chunks/`. A finished chunk gets a
  `.done` marker, and a rerun with the same input and chunk size skips
  it. Through `run_features.py` the scratch dir is kept after a failure,
  so rerunning the same command (or adding `--resume`) continues from
  the last finished chunk.
- **Parallel:** `--chunk-jobs N` chunks at once (default 4). This needs
  the parallel-run patch, step 6 above.
- **Same output:** chunk results are joined in input order into the
  usual `<stem>.localizer.raw.txt`, then parsed once. Output is identical
  to a single call. LOCALIZER scores each protein on its own, so the
  chunk boundaries don't affect results.
- **Failures:** a failing chunk stops the run with a non-zero exit, and
  no merged output is written. Its stderr is in
  `chunks/chunk_NNNNN.stderr.log`.
- `--chunk-size 0` restores the old single call.
- **CPU cap (`--java-cpus`, default 1):** each Java (WEKA) process
  LOCALIZER starts would otherwise size its garbage-collector and
  compiler threads to every core of the machine, so N chunks use far
  more than N cores. On a 128-core machine, 8 input files × 4 chunks
  reached about 170 cores. The wrapper now sets
  `JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=1 -XX:+UseSerialGC"` for
  every LOCALIZER call. Expect about `chunk-jobs × java-cpus` cores per
  input file, plus short pepstats/perl steps. `--java-cpus 0` turns the
  cap off. The JVM announces `JAVA_TOOL_OPTIONS` with one line on stderr,
  which LOCALIZER mixes into WEKA's output. The parser patch (step 6)
  skips that line, so the cap is applied only when `localization.py`
  carries the patch; otherwise the run prints a warning and goes ahead
  uncapped. To check that a JVM honours the flags:
  `JAVA_TOOL_OPTIONS="-XX:ActiveProcessorCount=1" java -XX:+PrintFlagsFinal -version | grep ActiveProcessorCount`.
- **Whole-machine budget:** the cap limits each process, not the machine.
  On a shared machine, also bound the whole run, e.g.
  `taskset -c 0-31 nice -n 10 python run_features.py ...`, and count every
  other job already running before choosing how many files run in
  parallel.

Output columns: `ID, CLS_binary_0no_1yes, MLS_binary_0no_1yes,
NLS_binary_0no_1yes, CLS_prob, MLS_prob, NLS_motif`. NLS is
motif-based and typically has no probability.

## Git Bash on Windows

If your LOCALIZER install lives inside WSL and you invoke
`run_features.py` from Git Bash specifically, Git Bash's MSYS layer
rewrites any argument that looks like an absolute POSIX path (e.g.
`--module-ref localization=/home/.../LOCALIZER.py`) into a Windows
path before Python sees it. PowerShell, cmd.exe, and native WSL bash
do not do this. If you hit it, prefix the command with
`MSYS_NO_PATHCONV=1`.
