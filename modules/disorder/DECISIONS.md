# IUPred Cysteine Handling — Decisions & Rationale

Covers how disulfide-bond bias in IUPred disorder scores was handled, and
why. Companion to `pipeline/gff_core_features/DECISIONS.md` and
`pipeline/pepstats/DECISIONS.md`.

## The problem

Your manuscript's Methods already state: "Protein structural disorder was
calculated using IUPred2... after removing cysteines from the protein
sequences in order to account for the possible presence of the disulfide
bonds, which can strongly affect ISD estimates (Uversky and Dunker 2010)."
IUPred is sequence-only and structure-blind — it has no way to know whether
a given cysteine forms a real disulfide bond, so a region held rigid by a
disulfide can still be scored as disorder-prone from local composition
alone. The correction exists to stop that specific miscall from entering
the disorder statistics.

**Checked first: does IUPred itself say anything about this?** No. Pulled
the actual README shipped with the IUPred3 install (written by Erdős,
Pajkos, and Dosztányi, the tool's authors) — it documents modes, smoothing,
and the method's statistical basis, with zero mention of cysteine or
disulfide handling. This is not an IUPred recommendation; it's a practice
from a different paper (Uversky & Dunker 2010) that researchers apply on
top of IUPred, which your manuscript already cites correctly. There is no
tool-authority answer to defer to — the choice has to be justified on the
mechanism itself.

## Options considered, and why three of four were rejected or deprecated

**Deletion** (literally removing `C` from the sequence, shortening it) —
rejected before implementation. DNG and AG measurably differ in cysteine
content (2.29% vs. 2.05% in A. thaliana, ~11% relative). Deletion would
shrink the two classes by different amounts as a side effect, reintroducing
a class-correlated *length* confound — precisely the kind of artifact this
entire pipeline rebuild exists to eliminate.

**Substitution** (`batch_iupred_features_cysstrip.py`, C→S before running
IUPred) — implemented, run on all three species, then deprecated after
inspection. Serine keeps sequence length constant (avoiding deletion's
problem), but Serine is a classically disorder-*promoting* residue
(Uversky/Radivojac classification), while cysteine is order-promoting. The
substitution doesn't neutralize the signal, it likely inflates it — and
since DNG/AG differ in Cys content, that inflation lands unevenly between
classes. Confirmed empirically: this version showed r=0.94–0.96 against the
no-strip baseline and a mean shift of ~0.05, both classes' means rising
under substitution — a real, non-trivial, composition-correlated effect,
not a neutral correction. **Output retained in `iupred_output_cysstrip/`
for the record, but not used for any analysis.**

**Exclusion from statistics** (`batch_iupred_features_cysexcl.py`) — the
version used going forward. IUPred runs on the real, unmodified sequence
(cysteine is a standard residue with well-defined behavior in the
algorithm — there's no computational reason to alter it). Cysteine
positions are then excluded from every aggregate statistic: per-mode
mean/median/min/max/std/frac≥threshold, `disorder_frac`, the N/C/mid
positional fractions, and IDR/ANCHOR2 segment-finding (segments computed on
a compacted mask with cysteine positions omitted entirely, not treated as
an ordered break — so a disordered stretch interrupted by a single
cysteine counts as one continuous run rather than being artificially
split). No invented residue identity, no sequence-length change.

**No correction at all** (`batch_iupred_features.py`, unmodified) — kept as
the supplementary robustness check, per plan (Option D primary, Option A
supplementary).

## Validation: exclusion vs. no-strip

| Species | r | mean\|diff\| | Cys excluded/protein (mean) |
|---|---|---|---|
| A. thaliana | 0.9999 | 0.0016 | 7.71 |
| B. rapa | 0.9999 | 0.0017 | 6.78 |
| O. sativa | 0.9999 | 0.0019 | 7.39 |

For comparison, the deprecated substitution version showed r=0.94–0.96 and
mean|diff|≈0.05 against the same baseline — roughly 30x larger. Exclusion
produces exactly the small, honest shift expected from dropping a handful
of untrusted data points out of an average over hundreds of residues;
substitution produced a shift an order of magnitude larger, which is itself
evidence the substitution was adding a real (unwanted) signal rather than
correcting one.

## Robustness result: DNG-vs-AG gap under both surviving versions

| Species | Gap, no-strip | Gap, cys-excluded |
|---|---|---|
| A. thaliana | −0.0184 | −0.0191 |
| B. rapa | +0.0329 | +0.0326 |
| O. sativa | +0.0987 | +0.0990 |

Sign and magnitude essentially unchanged in all three species. Two
conclusions: (1) the cysteine-handling choice does not materially affect
this dataset's disorder signal once done correctly — worth stating plainly
in the methods rather than treating as a live uncertainty; (2) the A.
thaliana anomaly (DNG showing *lower* mean disorder than AG, opposite to B.
rapa/O. sativa and to the literature's general direction) is not a
cysteine-handling artifact — it survives identically under both the flawed
and the corrected treatment, so whatever explains it, it isn't this. See
the separate note on A. thaliana's deeper DNG-curation filtering as the
more likely lead.

## Wording implication for the Methods text

"After removing cysteines from the protein sequences" (current text)
implies sequence-level editing. The implemented approach doesn't edit the
sequence — it excludes cysteine positions from the summary statistics. The
methods sentence should be updated to something like: "cysteine residues
were excluded from ISD summary statistics, since IUPred cannot account for
disulfide-bond-mediated rigidity (Uversky and Dunker 2010)." Same citation,
more accurate description of what was actually run.

## Output locations

- `iupred_output/` — no-strip (Option A), supplementary robustness check
- `iupred_output_cysexcl/` — cysteine-excluded-from-statistics (Option D), **primary**
- `iupred_output_cysstrip/` — Ser-substitution (Option C), **deprecated, do not use for analysis** — retained only so the reasoning trail above is checkable against real output
