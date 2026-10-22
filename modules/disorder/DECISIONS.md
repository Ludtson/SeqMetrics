# Disorder module: cysteine handling

## Problem

IUPred predicts disorder from sequence alone and cannot tell whether a
cysteine forms a disulfide bond. A region held rigid by a disulfide can
therefore be scored as disorder-prone. Removing cysteines before estimating
disorder is common practice (Uversky and Dunker 2010); IUPred's own
documentation does not address it.

## Approach used: exclude cysteines from the statistics

`batch_iupred_features_cysexcl.py` runs IUPred on the unmodified sequence,
then leaves cysteine positions out of every summary statistic:

- per-mode mean, median, min, max, SD and fraction above threshold;
- `disorder_frac` and the N-terminal, C-terminal and middle fractions;
- IDR and ANCHOR2 segments, computed on a mask with cysteine positions
  removed, so a disordered stretch interrupted by one cysteine counts as one
  segment.

Sequence length and residue identities are unchanged.

## Alternatives considered

- **Deleting cysteines** shortens sequences. Groups that differ in cysteine
  content are shortened by different amounts, adding a length bias.
- **Substituting serine for cysteine** keeps the length but replaces an
  order-promoting residue with a disorder-promoting one, which raises
  disorder scores.
- **No correction** (`batch_iupred_features.py`) remains available as a
  robustness check.

## Validation

Against uncorrected scores on three plant proteomes, the cysteine-excluded
scores correlate at r = 0.9999 (mean absolute difference about 0.002, with
roughly 7 cysteines excluded per protein). Serine substitution gave
r = 0.94–0.96 and a mean shift of about 0.05, with scores rising
systematically.
