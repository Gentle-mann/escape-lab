> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# FGFR2 bounded-feature sensitivity

**Post-result development sensitivity. The original failed run remains preserved.**

Biology features are limited to their training ranges at inference. Training matrices, coefficients, sequence PCA and all four feature arms are unchanged.

Incremental gate: **pass**; state-contrast check: **pass**.

| Arm | Mean-fold Spearman | Mean-fold MAE | High hits /26 at 68 picks | Strong hits /7 |
|---|---:|---:|---:|---:|
| baseline | 0.23373 | 0.53482 | 2 | 0 |
| inactive | 0.23298 | 0.53251 | 3 | 0 |
| active | 0.22697 | 0.54118 | 2 | 0 |
| contrast | 0.26784 | 0.53265 | 2 | 1 |
| historical_biology_hgb | 0.38350 | 0.49053 | 3 | 0 |

Spearman measures ranking, not percentage accuracy. MAE is on the unchanged log2(1 + activation) target. Historical HGB is an unchanged cached reference.

## Fixed gates

```json
{
  "incremental": {
    "status": "pass",
    "checks": {
      "mean_fold_rho_gain_003": true,
      "rho_improves_two_folds": true,
      "mae_within_5pct": true,
      "no_high_hit_loss": true
    },
    "rho_gain": 0.03410690578681416,
    "rho_fold_wins": 2,
    "high_hits_gain": 0
  },
  "state_contrast": {
    "status": "pass",
    "checks": {
      "rho_at_least_inactive": true,
      "rho_at_least_active": true,
      "mae_within_5pct_best_single": true
    },
    "interpretation": "Whether contrast beats the single-state controls; not a causal mechanism test."
  }
}
```

## Per-fold results

| Fold | Arm | Spearman | MAE | High hits |
|---|---|---:|---:|---:|
| pilot-0 | baseline | 0.14819 | 0.55482 | 1 |
| pilot-1 | baseline | 0.36340 | 0.49328 | 1 |
| pilot-2 | baseline | 0.18961 | 0.55635 | 0 |
| pilot-0 | inactive | 0.12860 | 0.54478 | 1 |
| pilot-1 | inactive | 0.41317 | 0.48212 | 2 |
| pilot-2 | inactive | 0.15716 | 0.57063 | 0 |
| pilot-0 | active | 0.14240 | 0.56163 | 1 |
| pilot-1 | active | 0.37682 | 0.49159 | 1 |
| pilot-2 | active | 0.16168 | 0.57031 | 0 |
| pilot-0 | contrast | 0.25082 | 0.53268 | 2 |
| pilot-1 | contrast | 0.46689 | 0.47088 | 0 |
| pilot-2 | contrast | 0.08580 | 0.59439 | 0 |
| pilot-0 | historical_biology_hgb | 0.42334 | 0.49386 | 2 |
| pilot-1 | historical_biology_hgb | 0.34259 | 0.47867 | 1 |
| pilot-2 | historical_biology_hgb | 0.38455 | 0.49907 | 0 |

## Verification and limits

- Unbounded score reproduction for every arm: {'baseline': {'n': 668, 'max_abs_difference': 0.0, 'tolerance': 1e-12}, 'inactive': {'n': 668, 'max_abs_difference': 0.0, 'tolerance': 1e-12}, 'active': {'n': 668, 'max_abs_difference': 0.0, 'tolerance': 1e-12}, 'contrast': {'n': 668, 'max_abs_difference': 0.0, 'tolerance': 1e-12}}.
- Protected files unchanged: 570.
- Fitting/scoring: 1.568s; total: 4.542s.
- Exactly the same 668 variants and three original folds; no new assay outcomes.
- Training-range bounds are fitted without labels and never updated using test rows.
- No parameter search, output clipping, external evaluation or selection of a new primary arm.
- Structural inputs remain prepared C491A-background proxies; inactive is a homology model.
- Static contacts are not free energies or mutant relaxation. Masking C491A cannot remove possible long-range effects of that shared engineered background.
- Block-resampling intervals are descriptive on reused development predictions.
