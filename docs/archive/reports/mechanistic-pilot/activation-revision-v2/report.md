> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# Second activation pilot

Primary challenger: `biology_tail`. Component: **inconclusive**. Resistance: **do_not_advance**.

Same 2719 variants and 11 aligned blocks as pilot v1. Fitting: 69.08s; total: 75.00s.

| Arm | Macro within-gene rho | MAE | Activation-class macro AP | Activation-class macro recall at 10% | Resistance macro AP | Macro Brier |
|---|---:|---:|---:|---:|---:|---:|
| original | 0.33902 | 0.24706 | 0.04632 | 0.09454 | 0.05961 | 0.05048 |
| capacity | 0.36730 | 0.26181 | 0.04817 | 0.04715 | 0.05817 | 0.05154 |
| biology | 0.38101 | 0.26329 | 0.05527 | 0.04441 | 0.06012 | 0.04790 |
| biology_tail | 0.38101 | 0.26329 | 0.05428 | 0.05784 | 0.06067 | 0.04859 |

Structural-only resistance baseline: AP 0.06016; Brier 0.04949.

Higher correlation, AP and recall are better; lower MAE and Brier are better. AP is not percentage accuracy. The two enriched arms share the same continuous regressor. Their difference is the additional activation classifier.

## Per-gene activation results

| Gene | Arm | Rho | Activation-class AP | Recall at 10% | Active / N |
|---|---|---:|---:|---:|---:|
| FGFR1 | original | 0.24011 | 0.00969 | 0.20000 | 5 / 690 |
| FGFR1 | capacity | 0.27024 | 0.00590 | 0.00000 | 5 / 690 |
| FGFR1 | biology | 0.24263 | 0.00618 | 0.00000 | 5 / 690 |
| FGFR1 | biology_tail | 0.24263 | 0.00958 | 0.00000 | 5 / 690 |
| FGFR2 | original | 0.19510 | 0.13469 | 0.07292 | 96 / 668 |
| FGFR2 | capacity | 0.24814 | 0.14793 | 0.08333 | 96 / 668 |
| FGFR2 | biology | 0.34645 | 0.18550 | 0.12500 | 96 / 668 |
| FGFR2 | biology_tail | 0.34645 | 0.15514 | 0.02083 | 96 / 668 |
| FGFR3 | original | 0.66227 | 0.00440 | 0.00000 | 2 / 686 |
| FGFR3 | capacity | 0.65082 | 0.00509 | 0.00000 | 2 / 686 |
| FGFR3 | biology | 0.62545 | 0.00373 | 0.00000 | 2 / 686 |
| FGFR3 | biology_tail | 0.62545 | 0.00462 | 0.00000 | 2 / 686 |
| FGFR4 | original | 0.25860 | 0.03651 | 0.10526 | 19 / 675 |
| FGFR4 | capacity | 0.30000 | 0.03377 | 0.10526 | 19 / 675 |
| FGFR4 | biology | 0.30952 | 0.02568 | 0.05263 | 19 / 675 |
| FGFR4 | biology_tail | 0.30952 | 0.04777 | 0.21053 | 19 / 675 |

## Decisions and uncertainty

```json
{
  "component_gate": {
    "status": "inconclusive",
    "reasons": [
      "undefined_gene_or_fold_metric"
    ]
  },
  "resistance_gate": {
    "status": "do_not_advance",
    "reasons": [
      "macro_ap_gain_below_0.005"
    ],
    "fold_macro_ap_deltas": [
      0.0014356938271851971,
      0.005315286001061667,
      -0.027290036586460745
    ]
  },
  "uncertainty": {
    "draws": 1000,
    "seed": 17,
    "groups": 11,
    "contrasts": {
      "macro_gene_spearman": {
        "percentile_2_5": -0.018395962726819967,
        "percentile_97_5": 0.11113778517503069,
        "valid_draws": 1000,
        "invalid_draws": 0
      },
      "macro_gene_active_ap": {
        "percentile_2_5": -0.004321277111154985,
        "percentile_97_5": 0.029268915730391015,
        "valid_draws": 892,
        "invalid_draws": 108
      },
      "resistance_macro_ap": {
        "percentile_2_5": -0.005528276996093139,
        "percentile_97_5": 0.006841497195013687,
        "valid_draws": 1000,
        "invalid_draws": 0
      }
    },
    "interpretation": "Conditional descriptive percentiles on 11 related development blocks, not independent biological confirmation. No refitting. Activation contrasts versus original; resistance versus structure."
  }
}
```

## Interpretation limits

- Original development outcomes and folds have been reused. This is a hypothesis-directed development comparison, not a new final test.
- The binary activation target is our protein-level adaptation of the study's replicate rule; it is not asserted to reproduce its nucleotide-level annotations.
- Added motif features are sequence distances, not mutant binding/free-energy predictions. AlphaMissense/EVE are auxiliary predictions, not activation labels. No curated outcome annotations or pretrained FGFR-specific supervised models were imported.
- Missing EVE scores are explicitly flagged and imputed only in training folds. Both activation heads are generated out of fold for downstream training.
- No LLM calls, hyperparameter search, external-test reveal or post-score model tuning occurred.
- 493 protected files stayed hash-identical; original activation and resistance predictions reproduce to 1e-12.

See results.json for all per-drug, fold, prior-reference, replicate and fitting diagnostics.
