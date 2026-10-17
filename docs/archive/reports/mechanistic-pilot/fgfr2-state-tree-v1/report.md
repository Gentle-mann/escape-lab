> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# Tiny FGFR2 structure-contrast tree pilot

Fixed development screening gate: **fail**.

Six fits only: pooled FGFR1–4 biology HGB baseline and the same HGB with cached FGFR2 active-minus-inactive descriptors, each on the three original folds. No tuning, agent loop, nested crossfits, resistance evaluation or external assay reads.

| Model | Mean-fold Spearman ↑ | Mean-fold MAE ↓ | High hits /26 | Strong hits /7 |
|---|---:|---:|---:|---:|
| baseline | 0.38350 | 0.49053 | 3 | 0 |
| contrast | 0.32900 | 0.54613 | 1 | 0 |

668 FGFR2 predictions; 68 total nominations. MAE uses the unchanged log2(1 + activation) target. Spearman is not percentage accuracy.

| Fold | Baseline rho | Contrast rho | Baseline MAE | Contrast MAE |
|---|---:|---:|---:|---:|
| pilot-0 | 0.42334 | 0.36410 | 0.49386 | 0.62039 |
| pilot-1 | 0.34259 | 0.29945 | 0.47867 | 0.49819 |
| pilot-2 | 0.38455 | 0.32347 | 0.49907 | 0.51982 |

## Gate and descriptive uncertainty

```json
{
  "gate": {
    "status": "fail",
    "checks": {
      "rho_gain_at_least_003": false,
      "rho_improves_two_folds": false,
      "mae_within_5pct": false,
      "no_high_hit_loss": false
    },
    "rho_gain": -0.05449109544727693,
    "rho_fold_wins": 0,
    "mae_reduction": -0.05560153125753592,
    "high_hits_gain": -2
  },
  "uncertainty": {
    "draws": 1000,
    "seed": 17,
    "unit": "aligned block within original fold",
    "contrasts": {
      "rho_gain": {
        "valid_draws": 1000,
        "invalid_draws": 0,
        "percentile_2_5": -0.16660979774503562,
        "percentile_97_5": 0.04254189288940501
      },
      "mae_reduction": {
        "valid_draws": 1000,
        "invalid_draws": 0,
        "percentile_2_5": -0.11337348287832626,
        "percentile_97_5": -0.0014219834244397149
      }
    },
    "interpretation": "Descriptive fixed-prediction resampling of reused development data; no refitting or external validation."
  }
}
```

## Verification and limits

- Baseline reproduction: {'n': 668, 'max_abs_difference': 0.0, 'tolerance': 1e-12}.
- 583 protected files unchanged.
- Six fits/scoring: 5.278s; total: 6.671s.
- Cached development labels only; no workbook ingestion or new outcome datasets.
- Same training rows/order, HGB hyperparameters and train-only imputation in both arms.
- Additional features are available only for FGFR2; other genes have zero geometry and availability=0. No new cross-gene structural mapping was introduced.
- No ESM, standardization, PCA, range clamp, output clipping or new biological features.
- Reused development folds cannot establish generalization; all external assays remain reserved.
- Structural features remain static prepared-coordinate proxies, not mutant energy or activity simulations.
