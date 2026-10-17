> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# FGFR2 active/inactive structural-proxy pilot

## Fixed question

Does adding a mutation-sensitive contrast between prepared active and inactive
FGFR2 structural environments improve continuous activation prediction over the
exact prior biology-plus-ESM ridge model? The task is retrospective preclinical
assay prediction on previously reused development folds. No external assay
outcomes, original final evaluation, or drug-resistance labels are used.

## Structure research and decision

The original active/inactive crystal-pair assumption is too strong. The 2017
mechanistic study explains that 2PSQ adopts a partially active conformation
influenced by crystal packing and sulfate. A 2023 PNAS study states that an
authentically autoinhibited FGFR2 crystal structure was unavailable and prepares
inactive FGFR2 from the FGFR1 3KY2 template, with active FGFR2 derived from 2PVF.
An inactive label in a structural catalog is insufficient evidence by itself.

Sources: [eLife 2017](https://elifesciences.org/articles/21137),
[PNAS 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC9974468/),
[RCSB 2PVF](https://www.rcsb.org/structure/2PVF),
[RCSB 3KY2](https://www.rcsb.org/structure/3KY2).

Use the prepared starting PDBs from the authors' two WT-labeled archives, linked
at the [original laboratory index](https://yzhang.hpc.nyu.edu/MD_data_FGFR/).
Only the first PDB member of each archive was extracted. No trajectories or
mutant archives are used. Coordinate sequences reveal that both preparations
retain C491A, so they are not exact assay-WT structures despite the archive
names. Both lack phosphotyrosine and have ATP. The inactive conformation remains
a homology model; these are matched structural proxies, not an experimental
WT state pair or an equilibrium ensemble.

Keep the common residue/atom set, mask engineered position 491 and uncertain
regions 509–511 and 581–596 in both query and contact environments. Do not drop
assay rows lacking geometry: use the same paired-availability flag and zero
geometry convention for every structural arm. Exact file hashes, sequence
mapping, state checks and exclusions are in `state-spec.json` and
`structure-audit.json`. Both states are DFG-in; their DFG phenylalanine side-chain
orientations and activation-loop geometries differ. Generic DFG-in/out labels
would misdescribe this comparison.

## Fixed features and models

The feature block contains 31 columns: six smooth local CA contact contexts,
eight mutation-chemistry interactions, twelve motif-proximity terms and their
interactions, four common-heavy-atom contact terms, and one common availability
flag. CA kernels use 8 Å and heavy-atom kernels use 5 Å, both with sixth-power
decay. Immediate sequence neighbors within two residues are excluded from
contact contexts. Mutation size uses side-chain heavy-atom-count change, not an
unverified volume or energy table. No known-activator labels or residue-specific
measured activities enter these descriptors.

The three structural arms have the same feature names, count and missingness
policy; effective matrix rank and variance may still differ. Contrasts are
active minus inactive before training-only preprocessing. The availability
flag is shared, not subtracted. Structural displacements or contact differences
are not free energies, catalytic rates or relaxed mutant structures. Paired-state
kinase research motivates testing the hypothesis, but its Potts/FEP calculations
are substantially more sophisticated than these static descriptors.
[Gizzio et al., Nature Communications 2024](https://www.nature.com/articles/s41467-024-50812-0).

Four fixed trained arms:

1. **baseline:** unchanged previous biology + ESM ridge.
2. **inactive:** baseline plus inactive-context descriptors.
3. **active:** baseline plus active-context descriptors.
4. **contrast, primary:** baseline plus active-minus-inactive descriptors.

Every arm uses ridge alpha=100, SVD solver, the same training-only biology
vectorization/imputation/scaling, and separate training-only PCA16 branches for
WT and signed mutant-minus-WT ESM embeddings. No hyperparameter selection,
class weighting, clipping, fine-tuning, or post-score revision. The baseline
must reproduce all previous cached ridge predictions to 1e-12 before new
structural arms are fitted. Keep the historical biology HGB predictions as an
unchanged reference-only row: beating the weaker matched ridge does not mean
beating the best historical continuous predictor.

## Cohort, target and evaluation

Retain exactly 668 complete eligible FGFR2 variants across ten occupied aligned
blocks and the prior global fold map: 275, 180 and 213 scored variants. Never
regenerate fold assignments after restricting to FGFR2. Continue all original
final and reserved external-panel group exclusions.

Continuous target is the unchanged log2(1 + median of four protein-level
replicate summaries), with synonymous codons collapsed within each replicate.
The primary metrics are equal-weight mean-fold Spearman correlation and MAE on
this scale. Pooled metrics are descriptive. Frozen secondary endpoints remain
high activation (>5 median and >=1.5 in three of four summaries) and strong
activation (>10 with the same consistency rule). Each fold nominates its top
ceil(10% of N), totaling 68 picks; ties resolve by candidate ID.

Incremental development gate, fixed before fitting:

- Mean-fold Spearman gain at least 0.03 over the matched baseline.
- Spearman improves in at least two of three folds.
- Mean-fold MAE no more than 5% worse than baseline.
- No reduction in the high-activation recovery count.

Separate state-contrast check: primary mean-fold Spearman at least as high as
both single-state arms; MAE no more than 5% worse than the better single-state
arm. Passing only the incremental gate does not support the need to subtract
states. Neither check establishes a molecular mechanism or external validity.

Use 1,000 paired aligned-block resamples within the original folds, seed 17,
for descriptive Spearman-gain and MAE-reduction percentiles versus baseline and
each single-state control. Count undefined draws. No refitting or independent
validation claims. Report every arm and failed criterion without promoting a
secondary result after scoring.

## Verification boundary

Before fitting, test structural identity/common-mask handling, rigid-transform
invariance, sign reversal on state swap, zero contrast for identical states,
feature-cache integrity, train-only preprocessing and exact matched-baseline
reproduction. Freeze source, feature, research and prior-result hashes in a new
run protocol. Preserve failed attempts and verify all protected hashes again
after scoring. Do not modify previously frozen modules or outcome datasets.
