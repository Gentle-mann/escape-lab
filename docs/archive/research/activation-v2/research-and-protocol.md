> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# Activation predictor diagnosis and second pilot

This work uses only previously authorized development data. RLY, MANO and original final aligned blocks stay excluded. Research sources inform feature design; no published mutation outcome list is used as an input feature. The original and first-pilot code/results remain unchanged.

## What failed locally

The saved first-pilot out-of-fold predictions establish three concrete weaknesses:

- Predictions compress the activation range: standard deviation 0.209 versus measured 0.458; pooled R² 0.070.
- Top-decile recovery within FGFR1/2/4 is 7.25%, 8.96%, 7.35%. Among the 136 development variants resistant to either drug, mean predicted activation is 0.920 versus measured 1.470 on the log2(1+enrichment) scale.
- The model has no position or local sequence context. Some 384 variants participate in identical-feature collisions despite different positions. Cα geometry is missing for 19.9% of variants; correlation in that subset is 0.019. Geometry-present variants also have tail failures, so this is only a partial explanation.

FGFR2 contributes 72.3% of squared prediction error. Pooled correlation can hide this weakness, although removing gene means still leaves correlation around 0.369. Gene identity alone therefore does not explain the entire signal.

See `diagnostics.json` and `diagnostics.md` for the complete aggregate audit. These diagnostics are post-hoc development analysis, not independent evidence.

## Research conclusions

1. **Activation direction must be learned explicitly.** The source FGFR study found that AlphaMissense, EVE and REVEL correlated better with depleted/loss-of-function variants than enriched/gain-of-function variants. It defines activating nucleotide variants by enrichment of at least 1.5 in at least three of four replicates. It also reports large variance for strongly activating variants. These findings support a separate activation-class head and argue against blindly downweighting noisy extremes. [Tangermann et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12807871/)

2. **Evolutionary scores can add information, but are not activation labels.** EVE learns evolutionary constraint from sequence alignments; AlphaMissense learns pathogenicity using structural/evolutionary representations and population-frequency supervision. Neither score should be equated with increased FGFR activity. In this pilot their relationship to the assay is learned within each training fold. [EVE](https://www.nature.com/articles/s41586-021-04043-8), [AlphaMissense](https://doi.org/10.1126/science.adg7492)

3. **Mechanistic context should describe regulatory regions, not only ligand proximity.** The FGFR molecular-brake study supports representing the conserved brake and activation-loop context. Our inexpensive implementation uses sequence distances and mutation-property interactions. It does not calculate mutation-induced energetics or an active/inactive structural transition. [Chen et al.](https://doi.org/10.1016/j.molcel.2007.06.028), [FGFR1 conformational dynamics](https://www.nature.com/articles/ncomms8877)

4. **Larger protein models are plausible future representations, not an automatic solution.** PreMode separates effect magnitude from gain/loss direction and highlights difficult gene-specific cases, including FGFR2. Its supervised labels and those of kinase-specific Activark require an overlap audit before importing pretrained predictions. Activark includes known-variant neighborhood features that could contaminate this evaluation. Borrowing outcome-independent feature ideas is safer for this pilot. [PreMode](https://www.nature.com/articles/s41467-025-62318-4), [Activark](https://link.springer.com/article/10.1186/s13073-025-01564-z)

5. **Simple augmented models are worth testing before large fine-tuning jobs.** Protein-fitness work shows the utility of combining evolutionary scores and assay observations with compact supervised models. ProteinGym distinguishes random and position-based evaluation; performance under easier random variant splits should not be assumed to transfer here. [Hsu et al.](https://www.nature.com/articles/s41587-021-01146-5), [ProteinGym](https://github.com/OATML-Markslab/ProteinGym)

6. **Replicate uncertainty is real but not identifiable from a single pooled target.** Enrich2 and DiMSum model experimental uncertainty using richer count information. Four summary replicates support reproducibility diagnostics here, but not their full measurement-error models or a rigorous prediction ceiling. [Enrich2](https://link.springer.com/article/10.1186/s13059-017-1272-5), [DiMSum](https://link.springer.com/article/10.1186/s13059-020-02091-3)

## Fixed small experiment

The cohort, exclusion groups, three outer folds, inner grouped cross-fitting and continuous target are exactly those of the first pilot. All 2,719 complete variants must be recovered; otherwise the run fails before scoring. No hyperparameter search or winner selection is performed.

Four arms:

1. `original`: exact first-pilot activation model.
2. `capacity`: original features, HGB with 200 iterations, depth 4, at most 15 leaves, minimum leaf size 10 and L2 regularization 1. Original settings were 100 iterations, depth 3, seven leaves, minimum leaf size 30 and L2 10.
3. `biology`: same higher-capacity regressor plus exact-reference local sequence context, smooth aligned-position features, sequence proximity to validated FGFR1 brake N546/E562/K638 and activation-loop Y653/Y654, mutation-property interactions, and audited AlphaMissense/EVE scores with missingness indicators.
4. `biology_tail`: the same continuous regressor plus a classifier for consistently high activation, using the same features and fixed tree settings. This is the predeclared primary challenger. Its continuous predictions are identical to `biology`; they are not two independent findings.

The binary target is a protein-level adaptation: at least three of four **per-replicate codon medians** reach 1.5 enrichment. It is not asserted to reproduce the paper's nucleotide-level class labels. It does not use resistance labels. A missing training class gives a deterministic training-prevalence fallback and an explicit diagnostic.

Both heads produce inner out-of-fold features for resistance training. The resistance predictor remains the fixed structural logistic model with C=0.1. The new head contributes its score and score × ligand-proximity interaction. Comparing `biology_tail` with `biology` isolates this additional head's contribution. No held-out activation measurements enter sequence-only predictions.

The primary component gate requires: macro within-gene Spearman gain ≥0.05, positive gains in at least two folds, no gene Spearman loss >0.05, MAE increase ≤5%, macro activation-class AP gain ≥0.02, and no gene activation-class AP loss >0.02. Undefined gene/fold class metrics or tail-training fallbacks make the gate inconclusive. These thresholds are engineering triage rules, not statistical significance.

The downstream gate is separate: resistance macro AP gain ≥0.005 over the structural baseline, gains in at least two folds, no drug AP loss >0.01, and drug Brier ≤ baseline ×1.05 +0.001. We also report whether the additional head helps versus the enriched regression alone.

All arms are reported, including failures. Metrics include per-gene and per-fold rank correlation, activation-class AP and top-10% recall, continuous error/spread, gene-prior references, split-half replicate agreement, and descriptive paired block-bootstrap intervals. Repeated use of these eleven related development blocks cannot establish external generalization.

## Deferred alternatives

- Frozen protein-language-model embeddings plus ridge: promising, but requires model download, isoform alignment and training-overlap checks. The existing evolutionary predictions are the smaller first test.
- Full pretrained gain/loss models: defer because curated FGFR training labels may overlap our reserved experiments.
- Molecular dynamics or active/inactive free-energy models: require validated structures, sampling and experimental calibration; this pilot cannot identify those parameters.
- Inverse-variance weighting: defer because it might suppress the high-activation variants we need to recover.

This document is prepared before the new numerical run. The run's `protocol.json` records exact settings, source hashes, exclusions and the code state before outcome ingestion/scoring.
