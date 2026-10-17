> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# What the small pilot establishes

The activation-informed approach has a promising measured-input diagnostic, but the current sequence-only implementation did not pass its advancement rule. Do not spend the reserved external tests on this implementation yet.

The run took 14.26 seconds, including ingestion, fitting, scoring and descriptive uncertainty. Model fitting took 9.60 seconds. It used no LLM calls or paid inference. Of 2,749 eligible development variants, 2,719 had all four activation replicate summaries; every comparison used this same cohort across 11 aligned blocks and three outer folds.

| Inputs available to resistance predictor | Macro average precision |
|---|---:|
| Existing structural features | 0.06016 |
| Structural features + predicted activation | 0.05961 |
| Structural features + measured activation | 0.11718 |
| Structural features + predicted shuffled activation | 0.05763 |

Average precision measures ranking of the rare resistant variants. These are not percentages of correct predictions, and the measured-input gain does not mean nearly twice the accuracy.

## Interpretation

**The measured activation assay adds useful information in this developmental comparison.** It improves ranking for both drugs, improves macro average precision in all three folds, and lowers macro Brier error from 0.04949 to 0.04164. Its descriptive block-bootstrap difference in macro average precision is positive, with 2.5th–97.5th percentiles of 0.01308 to 0.13811. This remains one study with only 11 related blocks, not proof of independent biological generalization.

**Our current predicted activation does not deliver that benefit.** Macro average precision falls by 0.00055, with improvement in only one of three folds. Its descriptive difference interval spans zero, from -0.00513 to 0.00320. The predeclared decision is `do_not_advance`.

**Activation prediction is a plausible immediate bottleneck, not a proven unique cause.** Predicted and measured activation correlate moderately overall (Spearman 0.376), but only 0.195 within FGFR2 and 0.240 within FGFR1. The predictor could miss the rare, extreme activation changes that matter for resistance even while learning some broad trends. Cross-fitting noise, the downstream model and transfer between cell backgrounds also remain possible contributors. This pilot does not distinguish them causally.

The measured condition assumes a new laboratory measurement is already available for each held-out variant. It is a valid assay-assisted task, but it cannot be presented as a sequence-only prediction. Both augmented conditions add an activation feature and its fixed interaction with ligand proximity, so the measured gain does not isolate the main activation feature from that interaction.

**Probability quality remains weak.** A post-hoc sanity check using only the existing aggregate counts shows that predicting each outer-training set's resistance prevalence would give macro Brier error 0.03083, better than even the measured-activation arm's 0.04164. The improved ranking therefore does not establish accurate resistance probabilities. This additional descriptive check involved no refitting and did not change the predeclared gate; its calculation is saved in calibration-sanity.json (`calibration-sanity.json`; local artifact).

## Next decision

Keep the full external evaluation closed. The next bounded development experiment should target activation-prediction quality, with explicit checks within each gene and for activation extremes, while retaining the same leakage controls and fixed resistance comparator. Probability calibration is also unresolved. Another defensible product scope is assay-assisted resistance prediction, which would require access to activation measurements at deployment.

This pilot provides no evidence about mutant binding-energy simulation, kinase kinetic parameters, cell-population evolution, drug combinations or cures. It tests only the first activation-informed part of the proposed architecture.

## Reproducibility and integrity

- Frozen protocol (`protocol.json`; local artifact), [aggregate results](results.json), and [full report](report.md) contain the exact settings and results.
- Eleven new tests and the full suite passed: 295 tests total; Ruff passed. One existing Starlette/httpx deprecation warning remains.
- All 475 protected source, test, data and report artifacts stayed hash-identical, including the original evaluation receipts. The original v2 protocol's 28 protected hashes also match.
- Original final variants and all prospective RLY/MANO aligned blocks were excluded. No external numerical outcomes were opened.
- This was one preregistered developmental run. No post-score model changes or tuning were performed.

To reproduce in a **new** output directory from the project root:

```bash
uv run python -m escape_lab.mechanistic_pilot_run --output reports/mechanistic-pilot/reproduction-01
```

A repeat is a reproduction of the same developmental test, not new evidence or an untouched final test.
