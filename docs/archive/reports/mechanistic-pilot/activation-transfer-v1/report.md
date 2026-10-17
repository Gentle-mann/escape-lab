> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# Activation-transfer viability pilot

Decision: **do_not_advance**. This is development triage, not a fresh test.

Candidates: 2719; aligned blocks: 11; total run: 14.26 seconds.

Four fixed arms use the same three outer group folds and resistance labels. Predicted/shuffled activation features are generated with nested group cross-fitting. Measured activation is an assay-assisted diagnostic, unavailable for sequence-only inference.

| Arm | Pem AP | Fut AP | Macro AP | Macro Brier |
|---|---:|---:|---:|---:|
| structure | 0.06276 | 0.05755 | 0.06016 | 0.04949 |
| predicted_activation | 0.06219 | 0.05703 | 0.05961 | 0.05048 |
| measured_activation | 0.15444 | 0.07992 | 0.11718 | 0.04164 |
| shuffled_activation | 0.05799 | 0.05728 | 0.05763 | 0.04878 |

AP is average precision, a ranking metric; it is not percentage accuracy. Higher AP and lower Brier are better.

## Sample size and fold results

- pem: 2702 conclusive variants, 119 resistant, 17 inconclusive labels; prevalence 0.04404.
- fut: 2690 conclusive variants, 53 resistant, 29 inconclusive labels; prevalence 0.01970.

| Fold | Baseline macro AP | Predicted activation | Measured activation | Shuffled activation | Pem positives | Fut positives |
|---|---:|---:|---:|---:|---:|---:|
| pilot-0 | 0.08905 | 0.08902 | 0.14998 | 0.08876 | 37 | 29 |
| pilot-1 | 0.04669 | 0.04747 | 0.11907 | 0.04751 | 55 | 14 |
| pilot-2 | 0.06846 | 0.06018 | 0.17777 | 0.03556 | 27 | 10 |

## Activation prediction

Target: log2(1 + median activation enrichment). The four assay replicates are summarized after collapsing nucleotide representations within each replicate.

| Arm / gene | Spearman | Mean absolute error | N |
|---|---:|---:|---:|
| predicted_activation / pooled | 0.3758 | 0.2471 | 2719 |
| predicted_activation / FGFR1 | 0.2401 | 0.1594 | 690 |
| predicted_activation / FGFR2 | 0.1951 | 0.4881 | 668 |
| predicted_activation / FGFR3 | 0.6623 | 0.1898 | 686 |
| predicted_activation / FGFR4 | 0.2586 | 0.1563 | 675 |
| shuffled_activation / pooled | 0.0352 | 0.2869 | 2719 |
| shuffled_activation / FGFR1 | -0.0939 | 0.1662 | 690 |
| shuffled_activation / FGFR2 | -0.1367 | 0.5046 | 668 |
| shuffled_activation / FGFR3 | -0.2380 | 0.3168 | 686 |
| shuffled_activation / FGFR4 | -0.0562 | 0.1644 | 675 |

## Predeclared decision

```json
{
  "status": "do_not_advance",
  "failed_guards": [
    "macro_ap_gain_below_0.005",
    "gain_in_fewer_than_two_folds"
  ],
  "inconclusive_reasons": [],
  "fold_macro_ap_deltas": {
    "pilot-0": -2.923075207697523e-05,
    "pilot-1": 0.0007773586175472275,
    "pilot-2": -0.008281120009413633
  },
  "interpretation": "Development triage only; thresholds are not significance tests."
}
```

## Descriptive uncertainty

- predicted_activation minus structure, block-bootstrap 2.5th–97.5th percentiles: [-0.00513, 0.00320].
- measured_activation minus structure, block-bootstrap 2.5th–97.5th percentiles: [0.01308, 0.13811].
- shuffled_activation minus structure, block-bootstrap 2.5th–97.5th percentiles: [-0.00883, 0.00091].

Descriptive only: 11 related blocks, no refitting, no biological replicate resampling; not a generalization confidence guarantee.

## Scope and integrity

- Original development and validation data have already informed prior work. This run is exploratory, despite grouped cross-validation.
- The 11 blocks are related sequence regions, not 11 independent biological experiments.
- Activation and resistance come from different cell backgrounds. Enrichment is not an isolated kinase activity or absolute growth-rate measurement.
- Structural distances are proxies, not simulated mutant binding affinities. This pilot does not test kinetics, drug combinations, cell evolution or patient outcomes.
- The predicted branch is nonlinear auxiliary supervision. Its features include predicted activation and its interaction with fixed ligand proximity; no LLM is used.
- A single within-gene shuffle is a capacity control, not a permutation significance test.
- Complete activation coverage defines the common cohort. Missing assay values are not imputed as zero.
- Original final and prospective RLY/MANO aligned blocks are excluded from activation and resistance fits. External numerical tables are not read.
- The trusted evaluator parses the monolithic original outcome JSON and source workbook but exposes only authorized development rows/columns to this pilot.
- 475 pre-existing source, test, data and report artifacts remained hash-identical. Original evaluation receipts are unchanged.

Protocol was written and hashed before numerical ingestion. No hyperparameter search, LLM calls, external-test opening or post-score model revision occurred.

Run directory: `reports/mechanistic-pilot/activation-transfer-v1`.
