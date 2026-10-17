> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# Second activation pilot: what improved and what did not

The revised predictor learns broader activation trends better, especially in FGFR2. It still does not recover the rare activating variants reliably enough to demonstrate improved resistance prediction. Keep the reserved external tests closed.

One fixed experiment took **75.00 seconds** on the same **2,719 variants and 11 aligned blocks** as pilot v1. The original activation predictions reproduced exactly, and both original resistance comparisons reproduced within 1e-12. No new LLM calls or paid inference were used.

## Controlled comparisons

| Model | Macro within-gene activation Spearman | FGFR2 Spearman | Activation MAE | Resistance macro AP |
|---|---:|---:|---:|---:|
| Original activation model | 0.33902 | 0.19510 | 0.24706 | 0.05961 |
| More capacity, same features | 0.36730 | 0.24814 | 0.26181 | 0.05817 |
| More capacity + biological features | 0.38101 | 0.34645 | 0.26329 | 0.06012 |
| Biological model + activation classifier | 0.38101 | 0.34645 | 0.26329 | 0.06067 |

The structural-only resistance baseline, without an activation predictor, has AP **0.06016**. Higher Spearman and AP are better; lower MAE is better. AP is not percentage accuracy. The two enriched arms share their continuous regressor, so the repeated correlation improvement is one result.

Most of the observed overall correlation improvement came from extra model capacity: +0.02828, followed by +0.01371 from the biological feature bundle. We did not isolate the contribution of individual features within that bundle. FGFR2 improved substantially, but FGFR3 correlation declined from 0.66227 to 0.62545.

The new feature bundle includes exact-reference flanking sequence, smooth aligned-position coordinates, sequence proximity to the molecular brake and activation-loop residues, mutation-property interactions, and AlphaMissense/EVE predictions with missingness flags. These are contextual proxies, not mutant structures or simulated free energies. AlphaMissense covers all 2,719 eligible variants; EVE covers 2,001, with 718 explicitly missing. Source score identities were reconciled against the study's protein labels and our assay sequences; the upstream predictor providers' full isoform-mapping procedures were not independently reconstructed.

## Why this is not a successful activation-tail model

There are only **122 activating protein variants** under the predeclared replicate-consistent definition:

| Gene | Activating variants | Complete variants |
|---|---:|---:|
| FGFR1 | 5 | 690 |
| FGFR2 | 96 | 668 |
| FGFR3 | 2 | 686 |
| FGFR4 | 19 | 675 |

This is a protein-level adaptation of the study's nucleotide-level rule, not an assertion that we reproduced its published class labels. One FGFR3 evaluation fold contains no activating examples. The predeclared component gate therefore returns **inconclusive**; no undefined gene was silently removed from a macro average.

Missing examples are not the only obstacle. The observed correlation gain of **0.04199** misses the **0.05** threshold, the activation-class AP gain of about **0.00796** misses **0.02**, and continuous MAE worsens **6.57%**, exceeding the **5%** guard. Across genes, activation recall among the top 10% of predictions falls from **9.45% to 5.78%** in the primary two-head arm.

The classifier is particularly unhelpful for FGFR2: the enriched continuous regressor retrieves **12 of 96** activating variants in its top 67 predictions, while the classifier retrieves **2 of 96**. The extra head must not be promoted as an improvement merely because its pooled downstream AP is marginally higher.

The resistance gate returns **do_not_advance**. The primary model improves AP over the structural baseline by only **0.00052**, far below the **0.005** threshold. Two folds improve, but the third loses **0.02729 AP**. The descriptive paired-block difference interval is **[-0.00553, 0.00684]**. This is not established resistance improvement.

## What we learned about the weakness

1. **The original representation had real gaps.** It omitted sequence position/context, and missing structure caused different variants to share identical features. The new model addresses those gaps, but doing so is insufficient for high-activation prioritization.
2. **The original capacity restricted useful trends.** The capacity-only ablation improves rank correlation. More capacity also worsens numerical error and does not help resistance, so simply enlarging the model is not the answer.
3. **Activation direction remains difficult.** Evolutionary/pathogenicity information can identify disrupted proteins without distinguishing increased from decreased activity. The source study explicitly finds stronger association with loss-of-function than gain-of-function. Our result is consistent with that limitation, not proof that any one feature caused it. [Original FGFR study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12807871/)
4. **The rare positive class is poorly supported outside FGFR2.** Two FGFR3 positives and five FGFR1 positives cannot support reliable three-fold per-gene classification. The 11 aligned blocks also provide little independent regional diversity.
5. **Assay variability contributes uncertainty.** Split-half replicate correlations are 0.399, 0.489, 0.855 and 0.684 for FGFR1–4. Median replicate spread is greater among activating examples. These summaries do not establish measurement noise as the primary cause or a hard performance ceiling.

## Decision

Retain the enriched regressor as a developmental candidate and the experiment as a reproducible ablation. Do not promote the extra classifier or claim a validated biological simulator. Do not consume external outcomes to select the next model.

The next defensible development study would focus on activation direction and extreme-variant ranking, with FGFR2 reported separately because it supplies most positive examples. Richer frozen sequence representations and a compact supervised model are candidates, but require provenance checks and the same position-grouped evaluation. Stronger evidence will ultimately require more independent activating examples; tuning the same eleven blocks cannot substitute for that.

## Verification and artifacts

- **348 tests passed:** the established frozen suite plus 53 new tests. Scoped Ruff passed. One existing Starlette/httpx deprecation warning remains.
- **493 protected files unchanged.** All resistance arms have identical candidate IDs, labels, groups and folds. Original evaluation receipts are preserved.
- RLY, MANO and original final groups remained excluded from all fits. No reserved external numerical outcomes were opened.
- The primary challenger, feature definitions, parameters and gates were frozen before the new run. No post-score model tuning occurred.
- The bootstrap uses fixed predictions, not refitting, across 11 related blocks. It is descriptive, not independent confirmation. There were 108/1,000 undefined activation-classification draws; all are counted. No resistance draws lacked positives in this run. The evaluator convention assigns zero AP difference to such resistance draws if they occur.
- A separate chat was concurrently developing You.com integration. Its new files and tests were outside this experiment's protected dependency inventory; this work did not modify them.

Files: frozen protocol (`protocol.json`; local artifact), [complete results](results.json), [full metrics](report.md), integrity verification (`verification.json`; local artifact), and [research and design rationale](../../../research/activation-v2/research-and-protocol.md).

Reproduce in a fresh directory from the project root:

```bash
uv run python -m escape_lab.activation_revision_run --output reports/mechanistic-pilot/activation-revision-reproduction-01
```

Repeating the run reproduces the same development experiment; it does not create new evidence.
