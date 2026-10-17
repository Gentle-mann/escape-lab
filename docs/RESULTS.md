# Results and claim boundaries

Escape Lab has a working autonomous experiment-selection pipeline and a reproducible improvement to one FGFR2 ranking implementation. **We have not demonstrated that an LLM improves resistance prediction, that the latest activation predictor generalizes to a new study, or that this system simulates or cures cancer.**

The complete [experiment ledger](EXPERIMENTS.md) preserves negative results, operational failures, controls, and experiments still in progress. These results were produced on 17 October 2026.

## Latest result: recover ranking information lost to numerical saturation

We changed the ranking score from a logistic classifier's sigmoid probability to that same fitted classifier's raw decision margin. All features, coefficients, folds, labels, and training settings remained the same. In the first fold, 136 of 275 probabilities had rounded to exactly 1.0. The raw margin recovered their ordering.

| FGFR2 high-activation ranking | Mean fold average precision | High variants retrieved /26 | Strong variants retrieved /7 |
|---|---:|---:|---:|
| Sequence classifier, sigmoid score | 0.09732 | 1 | 0 |
| **Same classifier, raw margin** | **0.11135** | **4** | **0** |
| Historical biology tree reference | 0.07186 | 3 | 0 |

The retrieval budget was 68 nominations, the top 10% rounded up separately within each of three fixed folds. The cohort contained 668 FGFR2 variants across ten occupied aligned-position blocks. “High” required at least three of four protein-level replicate summaries to be at least 1.5 and the four-replicate median to exceed 5. “Strong” additionally required the median to exceed 10. These are operational adaptations of the published assay categories.

**AP is a ranking metric, not percentage accuracy.** The new AP is 0.01402 higher. This is a technical correction implemented during development, not a new biological representation or an autonomous Nemotron/GLM discovery. The improvement occurred in one fold; the other two were unchanged. The descriptive paired-block interval for AP gain was approximately [0, 0.03261], including zero. The folds have been reused during development.

This is not the highest AP of every historical arm: the unbounded sequence ridge had AP 0.11708 and an inactive-structure ridge had AP 0.12227, but their continuous predictions suffered extreme extrapolation error. We preserve those results rather than selecting a favorable leaderboard. The margin correction does not resolve the classifier's underlying feature extrapolation or establish calibrated probabilities.

The tiny rerun used three fits: 0.31 seconds fitting/scoring and 2.92 seconds total with cached features. It reproduced the original probabilities exactly; all 594 protected artifacts remained unchanged. No new external assay outcomes were accessed.

Source: [ranking implementation](../src/escape_lab/fgfr2_sequence_ap.py). Evidence: [AP repair report](archive/reports/mechanistic-pilot/fgfr2-sequence-ap-v1/report.md), [decision and limitations](archive/reports/mechanistic-pilot/fgfr2-sequence-ap-v1/decision.md).

## Autonomous agent comparisons

These endpoints concern **drug resistance**, a different task from FGFR2 high activation. Their AP values are not interchangeable with the table above.

| Experiment | Comparator | Agent result | Conclusion |
|---|---:|---:|---|
| Original one-shot final resistance test | Fixed macro AP 0.10923 | Nemotron 0.09481 | Agent's development gain did not generalize |
| Revised v2 internal grouped evaluation | Fixed structural AP 0.03463; random/grid 0.03494 | Nemotron 0.03052; GLM-5.3 0.03052 | Stronger language model did not improve internal predictions |
| S13 biochemical panel, mean Spearman | Fixed structure 0.40568 | Nemotron 0.20892; GLM-5.3 0.29695 | GLM point estimate increased, but no established advantage; panel reused for GLM |
| You.com evidence pilot, reused S13 | GLM control 0.29695 | GLM with reviewed evidence 0.29695 | Identical selected model and predictions |

The tested swap was **GLM-5.3**, not GLM-5.2. Under the matched opportunity budget, GLM produced 6/12 valid proposals versus Nemotron's 12/12. Recorded wall times were 393.63 and 243.55 seconds; these are run timings, not a controlled serving-speed benchmark. Model-specific reasoning and sampling defaults were not controlled.

A repeated control/relevant-evidence/placebo study was started in the parallel research chat. No completed aggregate result was available at this publication snapshot; no benefit is claimed for it.

## Reported external-score references

The table reports average precision for **drug-resistance ranking** under
pemigatinib and futibatinib.

| Method | Pemigatinib AP | Futibatinib AP | Average AP |
|---|---:|---:|---:|
| EVE | 0.0894 | 0.0648 | 0.0771 |
| REVEL | 0.0671 | 0.0499 | 0.0585 |
| AlphaMissense | 0.0627 | 0.0500 | 0.0564 |

The exact scoring cohort, split, label definitions, score direction and
missing-score handling for this separate calculation are not archived in this
release. Our **0.11135 activation AP cannot be placed in this resistance table
as evidence of superiority**. Even our older resistance results require matching
those evaluation details before a head-to-head claim is justified. AP depends
on the evaluated cases and positive prevalence, as well as the model.

[EVE](https://www.nature.com/articles/s41586-021-04043-8),
[REVEL](https://pmc.ncbi.nlm.nih.gov/articles/PMC5065685/) and
[AlphaMissense](https://deepmind.google/research/publications/21083/) are
pretrained variant-effect prediction methods. Their scores serve here as
reference benchmarks for drug-resistance ranking.

The demonstrated six-hour achievement is the implemented research workflow and
the matched development ranking repair. Private specialized systems may perform
better; they were not evaluated. More time would allow stronger representations
and independent validation to be tested, without guaranteeing a performance gain.

## Evaluation status

| Evidence source | Status |
|---|---|
| Original resistance final holdout | **Consumed** in the initial Nemotron experiment; excluded from subsequent tuning |
| S13 biochemical panel | **Consumed** in v2; explicitly retrospective when reused for model swap and You.com |
| Development folds used by activation/sequence/structure pilots | **Reused development evidence**, not fresh confirmation |
| Reserved RLY-4008 S5 and MANO Data 1 outcome panels | **Not evaluated by these pilots**; remain reserved for a separately frozen external protocol |

The next credible scientific claim requires a frozen model and a suitable untouched study, with assay and endpoint compatibility established in advance. Repeatedly testing the same development folds cannot supply that confirmation. This is a preclinical research tool, with no patients, treatment recommendations, or claims of validated clinical performance.
