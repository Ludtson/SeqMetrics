# Install: localization (LOCALIZER)

LOCALIZER predicts chloroplast/mitochondrial transit peptides (cTP/mTP)
and nuclear localization signals (NLS) in plant protein sequences.
**Plant-specific** -- see README's note on a general-organism alternative
(TargetP or similar) for non-plant use, not yet built.

GPL-3.0 licensed, no registration needed. This repo ships the wrapper
scripts only (`run_localizer.py`, `localizer_bin.py`) -- never LOCALIZER
itself, per its own upstream policy (a prior integration attempt in this
lab already documented this exact rule: "LOCALIZER must be installed by
users from its official GitHub repository... do not redistribute bundled
LOCALIZER installs").

## Real dependency list, verified by direct install and test -- not just
## LOCALIZER's own documentation, which misses one of these

1. **LOCALIZER itself**:
   ```
   git clone https://github.com/JanaSperschneider/LOCALIZER.git
   cd LOCALIZER/Scripts
   unzip weka-3-6-12.zip
   ```
2. **Java** -- genuinely required (`java -cp weka.jar
   weka.classifiers.functions.SMO`, called directly from
   `localization.py`). **Not mentioned in this lab's own earlier HTLCP
   integration doc's dependency list** -- only found by reading
   LOCALIZER's actual source, not by trusting the existing setup guide.
3. **Perl** -- for NLStradamus, bundled as `nlstradamus.pl` inside the
   distribution; only the interpreter needs installing.
4. **Python 3** -- runs `LOCALIZER.py` itself.
5. **EMBOSS's `pepstats`** -- real gotcha, confirmed by direct test:
   LOCALIZER does **not** look for `pepstats` on PATH. It hardcodes an
   expected path in `LOCALIZER.py` (`SCRIPT_PATH + '/EMBOSS-6.5.7/emboss/'`,
   then string-concatenates `'pepstats'` directly onto it) and refuses to
   run at all if that exact directory doesn't exist -- even if a perfectly
   good `pepstats` is already on PATH. Compiling EMBOSS from source there
   (LOCALIZER's own documented approach) works but is slow and duplicates
   an EMBOSS install this project already has via conda for the
   `composition` module. Cheaper fix, used here:
   ```
   mkdir -p <LOCALIZER>/Scripts/EMBOSS-6.5.7/emboss
   ln -s $(which pepstats) <LOCALIZER>/Scripts/EMBOSS-6.5.7/emboss/pepstats
   ```
   Zero extra disk, satisfies LOCALIZER's hardcoded check exactly.

## One conda env, not two

All of the above (EMBOSS, Perl, Java, Python) live in the **same
`em_boss` env** already used by `composition` -- not a separate
`localizer`-specific env. A dedicated env was created and tested first,
then deliberately merged back into `em_boss` once it became clear it was
just duplicating a ~200MB EMBOSS install for no reason:
```
conda install -n em_boss -c conda-forge perl openjdk python=3.11
```

## Redistribution policy

LOCALIZER itself, WEKA, and NLStradamus are all third-party -- this repo
never bundles any of them, only the two wrapper scripts. Consistent with
the exact same rule already established (independently) in a prior
integration attempt in this lab (`htlcp`'s own
`README.md`/`BROAD_PIPELINE_SETUP.md`): fetch major dependencies from
official sources, never re-host them.

## Usage

```
python run_localizer.py proteins.faa \
    --localizer-script <path>/LOCALIZER/Scripts/LOCALIZER.py \
    --mode plant \
    --output OUTDIR \
    --localizer-parser localizer_bin.py
```

Real output columns, verified against real LOCALIZER predictions on its
own bundled Arabidopsis test set: `ID, CLS_binary_0no_1yes,
MLS_binary_0no_1yes, NLS_binary_0no_1yes, CLS_prob, MLS_prob, NLS_motif`
(chloroplast/mitochondrial/nuclear signal presence, probabilities where
LOCALIZER reports them -- NLS is motif-based and typically has none).

## Git-Bash-for-Windows gotcha with `--module-ref`

If your LOCALIZER install lives inside WSL (e.g.
`/home/user/.../LOCALIZER.py`, since Java/WEKA/Perl are more naturally
installed there) and you invoke `run_features.py` from **Git Bash on
Windows** specifically, Git Bash's own MSYS layer auto-converts any
argument that looks like an absolute POSIX path into a Windows path
*before Python ever sees it* -- confirmed by a real failure, the same
`--module-ref localization=/home/.../LOCALIZER.py` value arrived at
`run_features.py` mangled into a path under Git's own install directory.
This is not a SeqMetrics bug (PowerShell, cmd.exe, and native WSL bash
don't do this), but if you hit it from Git Bash, prefix the command with
`MSYS_NO_PATHCONV=1`.
