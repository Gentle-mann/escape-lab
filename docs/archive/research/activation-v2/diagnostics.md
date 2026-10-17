> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# Why the first activation predictor was weak

This is a post-hoc diagnosis of the existing activation-transfer-v1 out-of-fold predictions. It uses only the same 2,719 authorized development variants, with no new model fitting, no external numerical outcomes and no original final outcomes. All numerical outputs below are aggregates. Full results are in `diagnostics.json`.

## Main finding: the model misses the biologically interesting high tail

The pooled Spearman correlation of 0.376 hides a narrow prediction range and poor recovery of activation extremes. Measured activation has standard deviation 0.458, versus only 0.209 for predictions. Pooled R² is 0.070. The true top decile has mean activation 1.798, but the model predicts only 1.035 for those variants, in the pilot's log2(1 + enrichment) units.

| Gene | Spearman | Model RMSE | Training-gene-mean RMSE | R² | True SD | Predicted SD | Top-decile overlap |
|---|---:|---:|---:|---:|---:|---:|---:|
| FGFR1 | 0.240 | 0.278 | 0.271 | -0.088 | 0.267 | 0.124 | 7.25% |
| FGFR2 | 0.195 | 0.758 | 0.747 | -0.056 | 0.738 | 0.229 | 8.96% |
| FGFR3 | 0.662 | 0.273 | 0.364 | 0.409 | 0.356 | 0.213 | 15.94% |
| FGFR4 | 0.259 | 0.250 | 0.245 | -0.068 | 0.242 | 0.119 | 7.35% |

The mean reference predicts each gene's mean activation from the complementary outer training folds. It uses no held-out outcomes to form predictions and is computed without fitting a new model. The complex model improves substantially over this reference only for FGFR3. Global model RMSE is 0.442 versus 0.452 for this reference. Top-decile overlap ranks all out-of-fold predictions within a gene and compares the top `ceil(0.1*n)` sets; random expected overlap is approximately 10%. This is a descriptive retrospective statistic, not a preregistered test.

FGFR3's useful rank signal is not good evidence of gain-of-function recovery: even its true top-decile mean is 1.259 while the corresponding predicted mean is only 0.886. The model appears better at broad variation, including low activation, than at the high-activation events needed downstream.

## The lost signal is relevant to resistance

Among the 136 variants resistant to at least one drug, measured activation averages 1.470 but predicted activation averages 0.920. Mean underprediction is 0.550, and activation prediction Spearman within that subset is -0.134. Among the other 2,583 variants, measured/predicted means are 0.901/0.921 and Spearman is 0.403. These are descriptive outcome-defined subgroups, not independent validation.

| Direct ranking score, no fitted resistance model | Pem AP | Fut AP |
|---|---:|---:|
| Measured activation | 0.2363 | 0.1081 |
| Predicted activation | 0.0412 | 0.0183 |
| Resistance prevalence | 0.0440 | 0.0197 |

This diagnostic shows that the first predictor discards activation signal useful for resistance ranking. Measured activation at inference is an assay-assisted task, not sequence-only prediction. These direct rankings are diagnostic comparisons and do not replace the original frozen resistance evaluation.

## The representation has a concrete blind spot

The activation feature builder includes substitution chemistry, intrinsic distances from one inhibitor-bound FGFR1 structure and gene indicators. It excludes ligand and covalent-cysteine features, but contains **no aligned-position or local-sequence-context features**.

- C-alpha geometry is unavailable for 542/2,719 variants (19.9%). Within this subset, activation Spearman is 0.019, compared with 0.388 where geometry is available.
- Sidechain geometry is unavailable for 900/2,719 variants (33.1%).
- 384 variants occur in 144 repeated feature-vector groups containing distinct sequence positions. The model literally cannot distinguish those positions from the supplied inputs. There are 2,479 unique feature vectors for 2,719 variants.
- The mean within-collision activation standard deviation is 0.125. This quantifies information lost by the representation, but does not establish how much performance a new feature will recover.

Missing geometry cannot explain the whole failure: the geometry-present subset still has R² only 0.054 and substantially underestimates the high tail.

## Error is concentrated and pooled statistics can mislead

FGFR2 accounts for 72.3% of squared activation error despite being 24.6% of rows. Aligned blocks 056 and 060 together account for 47.9% of squared error. The pooled score therefore responds strongly to a few high-variance sequence regions. The training-gene-mean baseline also has weak ranking, and removing gene means leaves pooled Spearman at 0.369, so gene offsets alone do not explain the model's moderate correlation.

Fold Spearman values are 0.400, 0.366 and 0.344; R² values are 0.047, 0.182 and 0.001. The third fold has little magnitude-prediction benefit even though rank correlation looks moderate.

## Bounded next comparison

The evidence supports a small, fixed follow-up with the same cohort and grouped folds:

1. Keep the existing model as the anchor and add a fold-trained gene-mean reference.
2. Add exact-reference sequence context, aligned-position or motif-relative encodings, and interpretable mutation-by-intrinsic-geometry interactions. Preserve missingness indicators and train-only preprocessing. No external phenotype labels or hand-entered known activating mutations are needed.
3. Evaluate a modest change in model structure, such as a regularized nonlinear model with the new features and a separate high-activation classification output. Tail thresholds must be computed from each training partition, within gene if used; training features for the resistance model must remain nested out-of-fold predictions.
4. Report per-gene rank, RMSE, top-tail recovery and downstream resistance ranking together. A higher global correlation without better high-tail recovery is not sufficient progress.

Do not infer that a stronger language model can fix this numerical representation, that a new simulation is identified, or that lower global error establishes better resistance prediction. Replicate noise and the transfer between activation and resistance cell backgrounds remain unresolved by this OOF-only diagnosis.
