> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# Tiny sequence-classifier AP test

Technical ranking screen: **pass**.

Only change: rank by the same fitted classifier's decision margin instead of its sigmoid probability. Three original classifier fits, no parameter search.

| Ranking | Mean-fold high-activation AP ↑ | High hits /26 at 68 picks | Strong hits /7 |
|---|---:|---:|---:|
| probability | 0.09732 | 1 | 0 |
| margin | 0.11135 | 4 | 0 |
| historical_tree_reference | 0.07186 | 3 | 0 |

| Fold | Probability AP | Margin AP | Probability hits | Margin hits |
|---|---:|---:|---:|---:|
| pilot-0 | 0.02463 | 0.06671 | 0 | 3 |
| pilot-1 | 0.21570 | 0.21570 | 1 | 1 |
| pilot-2 | 0.05163 | 0.05163 | 0 | 0 |

## Scope

- Same 668 FGFR2 variants, original folds, target, embeddings, PCA, scaler, head and class weights.
- Raw margins resolve sigmoid rounding ties; they do not fix unsafe feature extrapolation, improve calibration or add biological information.
- Sequence AP superiority was concentrated in pilot-1; it was not consistent across folds.
- Development outcomes are reused. External assay outcomes and final evaluation remain reserved.
- Exact original probability reproduction: {'n': 668, 'max_abs_difference': 0.0, 'tolerance': 1e-12}.
- Three fits/scoring: 0.306s; total: 2.923s.
- Protected files unchanged: 594.

```json
{
  "gate": {
    "status": "pass",
    "checks": {
      "mean_fold_ap_gain_at_least_001": true,
      "no_fold_ap_loss": true,
      "no_high_hit_loss": true
    },
    "mean_fold_ap_gain": 0.014024831900723678,
    "high_hits_gain": 3,
    "interpretation": "Technical ranking-recovery screen, not independent model validation"
  },
  "uncertainty": {
    "draws": 1000,
    "valid_draws": 1000,
    "invalid_draws": 0,
    "seed": 17,
    "percentile_2_5": -1.1564823173178713e-18,
    "percentile_97_5": 0.03261486237766165,
    "interpretation": "Descriptive paired-block mean-fold AP gains on reused development data"
  },
  "diagnostics": [
    {
      "fold": "pilot-0",
      "n": 275,
      "probability_unique": 140,
      "margin_unique": 275,
      "probability_exact_one": 136,
      "margin_min": -31.762902310215782,
      "margin_max": 1084487.1128465498,
      "model": {
        "arm": "esm_biology_logistic",
        "training_n": 393,
        "high_count": 19,
        "training_groups": [
          "aligned-block-047",
          "aligned-block-056",
          "aligned-block-057",
          "aligned-block-060",
          "aligned-block-063",
          "aligned-block-064"
        ],
        "feature_count": 298,
        "biology_feature_count": 266,
        "sequence_components": {
          "wt": 16,
          "delta": 16
        },
        "coefficient_norm": 1.7206047661107624,
        "iterations": 88,
        "convergence_warning": false
      }
    },
    {
      "fold": "pilot-1",
      "n": 180,
      "probability_unique": 180,
      "margin_unique": 180,
      "probability_exact_one": 0,
      "margin_min": -35.658455487468174,
      "margin_max": 32.46199810133347,
      "model": {
        "arm": "esm_biology_logistic",
        "training_n": 488,
        "high_count": 19,
        "training_groups": [
          "aligned-block-038",
          "aligned-block-042",
          "aligned-block-052",
          "aligned-block-056",
          "aligned-block-060",
          "aligned-block-061",
          "aligned-block-063"
        ],
        "feature_count": 296,
        "biology_feature_count": 264,
        "sequence_components": {
          "wt": 16,
          "delta": 16
        },
        "coefficient_norm": 1.62002555432966,
        "iterations": 91,
        "convergence_warning": false
      }
    },
    {
      "fold": "pilot-2",
      "n": 213,
      "probability_unique": 213,
      "margin_unique": 213,
      "probability_exact_one": 0,
      "margin_min": -16.66027484126748,
      "margin_max": 5.306192043845723,
      "model": {
        "arm": "esm_biology_logistic",
        "training_n": 455,
        "high_count": 14,
        "training_groups": [
          "aligned-block-038",
          "aligned-block-042",
          "aligned-block-047",
          "aligned-block-052",
          "aligned-block-057",
          "aligned-block-061",
          "aligned-block-064"
        ],
        "feature_count": 300,
        "biology_feature_count": 268,
        "sequence_components": {
          "wt": 16,
          "delta": 16
        },
        "coefficient_norm": 1.588817095345273,
        "iterations": 71,
        "convergence_warning": false
      }
    }
  ]
}
```
