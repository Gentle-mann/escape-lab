> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# FGFR2 sequence-representation pilot

## Hypothesis and scope

The previous activation predictor learned broad trends but retrieved few activating FGFR2 mutations. We will test whether a frozen protein-language-model representation, coupled to a small supervised head trained specifically on FGFR2, improves retrieval of high-activation mutations in held-out sequence blocks. This is a developmental component test, not a cell simulation, external validation, or drug-resistance experiment.

All models and analysis rules below are fixed before fitting. No hyperparameter search, label-dependent feature extraction, post-result threshold selection, external numerical outcomes, original final outcomes, or resistance-label ingestion is allowed. The previously used development outcomes and folds are reused, so a positive result advances development rather than proving generalization.

## Research rationale

1. **Protein context can carry useful information beyond hand-designed features.** ESM-2 learns amino-acid representations from protein sequences; its research establishes that representations encode structural information. That is a reason to test a richer representation, not evidence that ESM-2 predicts FGFR activation directly. We use embeddings, not ESMFold coordinates or simulated molecular motion. [Lin et al., Science 2023](https://doi.org/10.1126/science.ade2574); [official ESM implementation](https://github.com/facebookresearch/esm).

2. **Mutation effects depend on sequence context.** Protein language models have demonstrated zero-shot variant-effect prediction on experimental assays. Our supervised test uses both the wild-type representation and the change induced by the mutant sequence. A language-model score measures sequence compatibility; its direction is not automatically gain versus loss of kinase activity. [Meier et al., NeurIPS 2021](https://papers.neurips.cc/paper_files/paper/2021/hash/f51338d736f95dd42427296047067694-Abstract.html).

3. **A larger model is not automatically a better small-data model.** A systematic transfer-learning study found useful pretrained features but no uniform scaling benefit across downstream protein tasks. A 35-million-parameter frozen checkpoint with a strongly regularized head is therefore a reasonable inexpensive experiment before larger models or fine-tuning. [Li et al., ICML 2024](https://proceedings.mlr.press/v235/li24a.html).

4. **The relevant biological endpoint is activation direction and strength.** The FGFR study uses replicate consistency and fixed enrichment categories. It also shows why generic deleteriousness predictions are insufficient proxies for activating mutations. We train on activation outcomes directly inside training folds and do not import known activating-mutation annotations. [Tangermann et al., Nature Genetics](https://www.nature.com/articles/s41588-025-02431-8).

5. **Retrieval needs a budget-based evaluation.** Precision-recall analysis characterizes the false-positive burden that matters for a rare-positive selection task. We report actual discoveries among a fixed nomination budget as well as AP, rather than optimizing broad correlation alone. [Davis and Goadrich, ICML 2006](https://ftp.cs.wisc.edu/machine-learning/shavlik-group/davis.icml06.pdf).

## Cohort and targets

Use the same complete FGFR2 development cohort from activation-revision-v2: **668 protein variants across ten occupied aligned blocks**. Reconstruct the original development/validation membership and reserved RLY/MANO block exclusions from metadata. Require identity with the previous complete cohort. Reuse the **saved global fold mapping** rather than rerunning block assignment after restricting to FGFR2. The original three folds contain 275, 180 and 213 variants.

The protein-level replicate summaries are medians across nucleotide representations within each replicate. A variant is replicate-consistently active if at least three of four summaries are >=1.5. This is a protein-level adaptation of the publication's nucleotide-level rule, not an exact reconstruction of its published classifications.

- **Primary high activation:** replicate-consistent active and median enrichment **>5**. This combines the paper's activating and strongly activating categories. The 26 positives are distributed 7, 7 and 12 across the existing folds.
- **Secondary strong activation:** the same consistency rule and median enrichment **>10**, matching the paper's strongly activating category. There are seven positives, distributed 2, 0 and 5. The zero-positive fold's AP is undefined and must be reported explicitly.
- **Continuous target for the ridge ablation:** the existing `log2(1 + median enrichment)` target, unchanged from the prior pilot.

These cutoffs were selected from the paper before the threshold-count audit and without model-performance inspection. See `threshold-audit.json` for provenance and counts. They must not be adjusted after fitting.

## Frozen representations and four arms

The checkpoint is `facebook/esm2_t12_35M_UR50D`, pinned to revision `6fbf070e65b0b7291e7bbcd451118c216cff79d8`. Use the exact assay reference sequence with residue identity and numbering validation. The model is frozen and extracts wild-type and mutant-sequence representations without outcomes. Represent each variant with two branches: wild-type context and the mutant-minus-wild-type representation. Fit a separate scaler and **16 principal components per branch**, giving **32 sequence components total**, using training-fold data only. Record the exact extraction implementation, checkpoint-file hashes and cache manifest before scoring.

The new supervised models use FGFR2 training variants only. Existing biological features are the outcome-independent activation-v2 bundle; their feature definitions remain unchanged. Training-fold-only imputation, scaling and PCA are mandatory. No pooled-label feature selection or full-cohort PCA is permitted.

| Arm | Representation and model | Role |
|---|---|---|
| `reference` | Saved enriched activation-v2 regressor's FGFR2 out-of-fold scores | Historical matched-cohort reference; trained in the original all-gene experiment |
| `biology_logistic` | Existing biological features; logistic regression, C=0.1, balanced class weights | Matched FGFR2 classifier, isolating the effect of new sequence features |
| `esm_biology_logistic` | Existing biological features plus 16 PCA components from each sequence branch, 32 total; same C=0.1 balanced logistic regression | **Primary challenger** |
| `esm_biology_ridge` | Same combined representation; ridge regression, alpha=100, continuous target | Fixed objective ablation, reported regardless of result |

The classifier target is primary high activation. Class weights are computed from training labels only. A single-class training fold fails closed before scoring; it cannot yield a passed pilot. Class-balanced probabilities are not claimed to be calibrated risks. The continuous ridge scores are rankings, not class probabilities. The primary model is not replaced by the best observed arm after results arrive.

## Metrics, gate and uncertainty

Nominate the top `ceil(0.10 * fold_n)` variants **within each fold**. The three budgets are 28, 18 and 22, totaling **68 nominations**. Resolve ties deterministically by candidate ID. Report recovered high-activation positives and recall, precision, AP and prevalence per fold. Mean-fold AP gives equal weight to the three fold scores; pooled out-of-fold AP is secondary because independently fitted models can have different score scales. Report the >10 endpoint separately with undefined metrics preserved.

The primary challenger advances only if all conditions hold:

1. Against `reference`, recover at least **three additional high-activation positives** within the 68 nominations and improve mean-fold AP by **at least 0.03**.
2. Against `biology_logistic`, do not lose recovered positives and improve mean-fold AP by **at least 0.03**.
3. Improve AP against `reference` in at least **two of three folds**.
4. All required metrics are defined, predictions are complete, boundaries pass, and no training fallback occurred.

These are engineering triage rules, not significance thresholds. Failure of the primary gate remains a failure even if the ridge ablation or secondary >10 endpoint looks better. Report every arm, runtime and diagnostic, including unsuccessful revisions.

Paired block resampling of fixed out-of-fold predictions provides descriptive uncertainty only. Count invalid draws and avoid treating the 668 variants as 668 independent biological experiments. Ten related sequence blocks, reused across development rounds, are insufficient for a claim of external generalization. A passed component gate would support another development step; it would not establish improved resistance prediction or justify revealing the protected external tests.

## Verification before execution

Tests must establish exclusion enforcement, exact saved fold membership, train-only preprocessing, held-out-label invariance, deterministic nominations, strict >5 and >10 boundaries, correct handling of an empty secondary-positive fold, and gate behavior at its thresholds. Preserve original sources and results. The run writes a new protocol and immutable result directory, records hashes, and verifies protected artifacts are unchanged after execution.
