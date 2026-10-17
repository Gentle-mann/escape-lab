> Historical research record. See [current results](../../RESULTS.md) for the release interpretation.

# FGFR revision v2: research and locked analysis plan

## What failed and what this revision tests

The first Nemotron revision changed logistic C from 1 to 10. It improved its
development and validation average precision (AP), then decreased held-out macro
AP from 0.10923 to 0.09481. That final set is spent and will not be reused.
Development-only diagnostics show coefficient norms increased from roughly 8–10
to 20–23 and probability error (Brier) worsened in every fold for both drugs.
Prevalence varied substantially between mutation regions. This is evidence of
unstable extrapolation and a weak selection rule, not proof of a single cause.

The old representation combined amino-acid categories, local sequence, and
absolute position kernels without scaling or drug-specific geometry. It could
memorize regional associations but did not encode whether a residue contacts the
particular drug. A larger C weakened shrinkage. The previous planner also
misdescribed kernel width as bins; the new prompt specifies exact semantics.

## Research decisions

1. Add experimentally grounded geometry from FGFR1–pemigatinib [7WCL](https://www.rcsb.org/structure/7WCL)
   and FGFR1–futibatinib [6MZW](https://www.rcsb.org/structure/6MZW).
   Align actual structure sequences to assay references; do not assume equal
   numbering. These are homolog-template proxies for FGFR2–4, not mutant poses.
   Both structures have engineered C584S. Mask side-chain features at mismatches;
   preserve missing-geometry flags. Futibatinib's covalent C488 contact differs
   materially from pemigatinib. Exclude the C488A 6MZQ template.
2. Use compact chemistry, ligand proximity, solvent exposure, local packing, and
   kinase-anchor distances, with at most twelve specified chemistry/geometry
   interactions and four explicit gene indicators. Family conservation was
   considered but is not implemented; no large evolutionary alignment is claimed.
   [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/methodology.html)
   and [Bio.PDB SASA](https://biopython.org/docs/latest/api/Bio.PDB.SASA.html)
   informed mapping and exposure computation.
3. Fit imputation and scaling inside each training fold. Compare regularized
   logistic models with shallow histogram boosting. Generic protein-language
   model pathogenicity is deferred: loss of protein function is not equivalent
   to drug resistance. [ESM variant prediction](https://github.com/facebookresearch/esm)
   is a plausible later feature, not ground truth or a substitute for this assay.
4. Separate experiment planning from numerical prediction. Nemotron selects a
   bounded candidate with a falsifiable hypothesis, sees only inner-fold
   aggregates, and retains experiment memory within that search. Reset memory
   between outer folds. This follows evaluated-search ideas in
   [AIDE](https://arxiv.org/abs/2502.13138) and
   [AIRA](https://arxiv.org/abs/2507.02554), without allowing executable agent edits.
5. Selection itself can overfit noisy cross-validation
   ([Cawley and Talbot](https://www.jmlr.org/papers/v11/cawley10a.html)). Evaluate
   the complete selection procedure in nested grouped folds
   ([Varma and Simon](https://doi.org/10.1186/1471-2105-7-91)). Compare with
   [random search](https://jmlr.org/papers/v13/bergstra12a.html) and exhaustive
   numerical search. Brier measures total probability error, not just calibration;
   avoid unstable post-hoc isotonic calibration on so few groups
   ([scikit-learn](https://scikit-learn.org/stable/modules/calibration.html)).

## Implementation and verification boundary

New modules implement v2; preserve the old frozen predictor, evaluator, receipts,
and failed attempts. Tests cover leakage boundaries, train-only preprocessing,
missing coordinates, strict proposals, selection regressions, and the one-shot
external evaluator. Existing tests remain required. A source/data/feature hash
manifest is written before running v2 scoring. A separate committed prediction
hash is written before opening new external numerical outcomes.

## Predeclared comparison (before v2 scoring)

Training pool: original development plus validation, excluding the old final and
all nine aligned blocks containing the new external single-mutant panel. This
leaves 3,730 protein variants in 15 blocks. Blocks group aligned positions across
all four genes. No outcome-based split construction. Deterministic hashed block
assignment gives three outer folds; each outer training pool has three inner folds.
These old development data were already used in v1, so nested results are exploratory.

Catalog: chemistry, structure, and interaction representations each with logistic
C = 0.01, 0.1, 1; structure boosting depth 2 or 3; interaction boosting depth 2.
Twelve candidates. The reference is chemistry-logistic-c0.1. Fixed biology is
structure-logistic-c0.1. The v1 comparator is the original fixed C1 model retrained
on exactly the same eligible examples, not the old fitted model.

All selectors use the same acceptance rule. A candidate must improve pooled
macro AP over the incumbent; versus the fixed reference, neither drug's AP may
fall by more than 0.01, each drug's Brier must be at most reference × 1.05 + 0.001,
at least two of three fold macro AP values must be no worse, and no fold macro AP
may fall by more than 0.01. Among eligible candidates retain the greatest pooled
macro AP; ties retain the incumbent. These thresholds are engineering guards,
not statistical significance tests. If none qualifies, retain the reference.

Controls: fixed v1, fixed biology, random search with three new candidates, full
twelve-candidate grid, and Nemotron with three proposal opportunities. Repeat
selection independently inside each outer fold and once on all eligible data
for the external prediction. The grid has a larger search budget and is an upper
control for this finite catalog, not a cost-matched comparison. Cache numerical
evaluations, but never show hidden control results to the planner. Maximum 12
proposal opportunities, each with at most one format repair; one shared client
token budget, no automatic HTTP retries, 60-second HTTP timeout. Record actual
calls, failures, token usage, and wall time.
The shared budget is 180,000 tokens (a stop threshold, not a billing ceiling)
and 1,200 seconds, checked between bounded model fits and HTTP calls. Search
status is completed with three evaluated proposals, degraded with one or two,
and unavailable with none. External values remain unopened if any search is
unavailable. Degraded comparisons are labeled explicitly. Shared numerical
caching makes wall time unsuitable as a direct per-method efficiency comparison.

## Independent challenge, not another cellular classification test

The source is [Saifudeen et al., Nature Biotechnology (2026)](https://doi.org/10.1038/s41587-026-03090-8),
Supplementary Table S13: purified-protein Reaction Biology HotSpot measurements.
Metadata inspection identified 25 WT-matched single mutants across FGFR1–4.
Exclude double mutants, a fusion, and FGFR2 M420I, whose WT does not match our
reference. Freeze this exclusion before reading values. Preserve the separate
RLY-4008 IC50 panel unopened for later work.

All five methods commit predictions before the S13 numerical cells are opened.
Primary external measure: per-drug Spearman correlation between predicted cellular
resistance probability and measured residual biochemical kinase activity, plus
their mean. At least five finite measurements and nonconstant scores/outcomes
are required per drug; otherwise the result is inconclusive. No external AP,
Brier, clinical accuracy, or invented binary resistance threshold. WT subtraction
is disabled because matched ATP and concentration conditions are not established.
Missing/censored cells are explicit, never zeros. Reserve a one-shot receipt
before label parsing; failures after reservation also spend the evaluation.

This is independent measurement transfer, not a new biological discovery.
Published hotspot knowledge and model-pretraining contamination remain possible;
some discovery snippets were already visible. New workbook numerical cells were
unopened during design. Twenty-five selected mutants in nine regions are a small,
biased sample. Better internal AP without better external rank correlation is
not generalization success. An agent that merely chooses the same numerical
configuration as random/grid has not demonstrated an independent AI advantage.

## Success, failure, and stopping

Report feature gain versus retrained v1, and agent gain versus fixed biology,
matched random, and grid separately. Any improvement is descriptive unless its
uncertainty and replication support more. Negative or mixed results remain in
the report. Freeze and evaluate once, then stop; do not revise against opened
external outcomes. The system remains a preclinical research prototype, with no
patients or treatment recommendations.

For descriptive uncertainty, use 2,000 paired bootstrap resamples of entire
aligned blocks (seed 6100), recalculating the primary metric difference for each
method against v1 and for the agent against the other controls. Report percentile
2.5–97.5% intervals and valid replicate counts. Undefined correlations remain
undefined. These intervals do not correct for prior development-data reuse,
the selected external panel, or broader assay/domain shifts; they are not a
claim of independent clinical evidence. Bootstrap results never affect selection.
