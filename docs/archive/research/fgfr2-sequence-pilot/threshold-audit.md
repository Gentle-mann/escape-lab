> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# FGFR2 high-activation endpoint and evaluation audit

This is a read-only audit of existing, authorized development targets, not a new model fit. No external numerical outcomes, original final outcomes, or source workbook cells were opened. The saved activation-revision-v2 target records were selected only after reconstructing the RLY/MANO metadata exclusions and requiring original development/validation membership. `threshold-audit.json` records counts and the target-source hash.

## Fixed biological definition

The original study requires enrichment at least 1.5 in three of four biological replicates, then categorizes median enrichment as potentially activating (1.5–2.5), weakly activating (>2.5–5), activating (>5–10), and strongly activating (>10). We use the paper's fixed >5 boundary for the primary high-activation target and >10 as the prespecified strong-activation secondary target. These choices were communicated before the count audit, not selected by model performance. [Tangermann et al., Nature Genetics](https://www.nature.com/articles/s41588-025-02431-8)

Our protein-level version collapses synonymous codon representations within each replicate before applying the replicate rule and median threshold. It is an adaptation, not an assertion that we exactly reproduce the published nucleotide-level classifications. Here all >5 and >10 variants also satisfy the consistency rule. That observation is a cohort property, not a reason to omit the rule from code.

| Existing fold | Complete variants | Replicate-consistent >=1.5 | Consistent and median >2.5 | Primary >5 | Secondary >10 |
|---|---:|---:|---:|---:|---:|
| pilot-0 | 275 | 23 | 17 | 7 | 2 |
| pilot-1 | 180 | 36 | 25 | 7 | 0 |
| pilot-2 | 213 | 37 | 18 | 12 | 5 |
| Total | 668 | 96 | 60 | 26 | 7 |

There are ten occupied FGFR2 aligned blocks, while the prior all-gene cohort occupied eleven. **Preserve the saved global fold assignment.** Calling the hash-and-round-robin fold function on FGFR2 alone can change fold membership because one block disappears. Every variant from one aligned block must stay together, and any cross-gene training arm must use the same global exclusion.

## Recommended narrow experiment

Freeze one primary challenger and its numerical settings before fitting. Compare it with both the saved enriched-regression ranking and an FGFR2-only version of the same downstream head using the existing biological features. The latter isolates new sequence information from a change in cohort, objective, class weighting, or model class. A sequence-only arm is an interpretable ablation; do not choose the best of several arms after scoring and call it the primary.

Extract frozen protein-language-model features without activation labels. Freeze the checkpoint, exact assay sequence, residue numbering, feature definition and checkpoint hash. Train scaling, imputation, PCA and any feature selection only on each training fold. Frozen unsupervised protein pretraining may contain the wild-type sequence, which is different from supervised access to this assay's targets; no curated mutation outcomes or pretrained FGFR activity annotations should enter the features.

A small regularized supervised head is more appropriate for this first test than training or fine-tuning a large sequence model on 26 positive examples. Protein transfer-learning experiments find useful pretrained features but do not show uniform benefit from increasing model size. This supports testing a compact representation before assuming a larger one solves directionality. [Li et al., ICML 2024](https://proceedings.mlr.press/v235/li24a.html)

## Metrics and advancement recommendation

The primary question is whether a fixed assay budget retrieves more high-activation mutations in unseen sequence blocks. Report the count of positives retrieved from the top 10% **within each fold**, totaling 28 + 18 + 22 = 68 nominated variants. Ties need a prespecified deterministic rule. This avoids conflating fold-specific score scales when ranking the pooled out-of-fold predictions.

Also report equal-weight mean fold average precision (AP), each fold's AP and prevalence, and pooled AP as a secondary descriptive measure. AP exposes the false-positive burden relevant to this rare-positive retrieval task; it is not percentage accuracy. ROC AUC or partial AUC at FPR <=0.1 can be secondary but should not substitute for the actual selection-budget metric. PR and ROC summaries answer related but different questions. [Davis and Goadrich, ICML 2006](https://ftp.cs.wisc.edu/machine-learning/shavlik-group/davis.icml06.pdf)

Suggested engineering gate, to be frozen by the evaluator before fitting: at least three additional high-activation recoveries among the fixed 68 nominees; mean-fold AP gain at least 0.03 against the existing enriched ranking; AP gains in at least two of three folds; and no decline in recovery against the matched existing-feature classifier. Require no single-class training fallback. These are triage thresholds, not statistical significance, and may be adjusted before fits with the reason recorded.

Strong activation (>10) has only seven positives and none in pilot-1. Report its retrieval count and the undefined fold AP explicitly. Do not require three defined strong-activation fold metrics or silently replace the zero-positive fold with an easier split. Do not let this secondary endpoint become the success criterion after seeing the primary result.

Use paired block-bootstrap intervals on the fixed out-of-fold predictions as descriptive uncertainty only, with invalid draws counted. With ten related development blocks reused across several model revisions, a narrow interval would still not constitute independent validation. Improvements qualify the next developmental step, not a claim that the model generalizes externally or predicts drug resistance better. This experiment need not train any resistance model or read resistance labels.

## Adversarial failure checks

- If the new model merely recognizes each residue position, grouping should expose that failure; random variant splits would hide it.
- If improvement appears only after pooling scores from different fitted folds, require the within-fold retrieval and AP evidence before advancing.
- If classifier probabilities are distorted by class balancing, call outputs ranking scores unless calibrated using training-only splits. Calibration is not required for this retrieval pilot.
- If the model finds generic damaging variants but not activating variants, zero-shot evolutionary plausibility is not a sufficient surrogate for direction. Report the negative result.
- If the >10 endpoint improves while >5 fails, report an exploratory secondary signal, not a passed primary test.
- If a pretrained model includes supervised mutation labels from this study or the reserved panels, reject those predictions for this pilot even if their apparent performance is excellent.
