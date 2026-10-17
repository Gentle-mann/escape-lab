> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# FGFR2 conformation-descriptor pilot

Primary incremental gate: **fail**. State-contrast check: **fail**.

Continuous target: log2(1 + median activation), unchanged four-replicate protein summaries. Same 668 variants, ten occupied position blocks, and three original folds.

| Arm | Mean-fold Spearman | Mean-fold MAE | High hits /26 at 68 picks | Strong hits /7 |
|---|---:|---:|---:|---:|
| baseline | 0.19585 | 3933.44035 | 4 | 0 |
| inactive | 0.20696 | 1102.51852 | 5 | 0 |
| active | 0.19078 | 1777.89393 | 5 | 0 |
| contrast | 0.22328 | 21494.56594 | 5 | 0 |
| historical_biology_hgb | 0.38350 | 0.49053 | 3 | 0 |

Higher Spearman and retrieval, lower MAE are better. These are not percentage accuracy. Historical HGB is an unchanged cached reference, not a newly trained arm.

## Fixed decisions

```json
{
  "incremental": {
    "status": "fail",
    "checks": {
      "mean_fold_rho_gain_003": false,
      "rho_improves_two_folds": true,
      "mae_within_5pct": false,
      "no_high_hit_loss": true
    },
    "rho_gain": 0.027433464629071397,
    "rho_fold_wins": 2,
    "high_hits_gain": 1
  },
  "state_contrast": {
    "status": "fail",
    "checks": {
      "rho_at_least_inactive": true,
      "rho_at_least_active": true,
      "mae_within_5pct_best_single": false
    },
    "interpretation": "Whether contrast beats the single-state controls; not a causal mechanism test."
  }
}
```

## Fold results

| Fold | Arm | Spearman | MAE | High hits |
|---|---|---:|---:|---:|
| pilot-0 | baseline | 0.17358 | 11798.00972 | 3 |
| pilot-1 | baseline | 0.21538 | 1.75335 | 1 |
| pilot-2 | baseline | 0.19858 | 0.55797 | 0 |
| pilot-0 | inactive | 0.17661 | 3304.99596 | 3 |
| pilot-1 | inactive | 0.28583 | 1.98713 | 2 |
| pilot-2 | inactive | 0.15845 | 0.57247 | 0 |
| pilot-0 | active | 0.16647 | 5331.21367 | 3 |
| pilot-1 | active | 0.23829 | 1.89610 | 2 |
| pilot-2 | active | 0.16758 | 0.57201 | 0 |
| pilot-0 | contrast | 0.26320 | 64480.82726 | 3 |
| pilot-1 | contrast | 0.32493 | 2.18308 | 1 |
| pilot-2 | contrast | 0.08172 | 0.68748 | 1 |
| pilot-0 | historical_biology_hgb | 0.42334 | 0.49386 | 2 |
| pilot-1 | historical_biology_hgb | 0.34259 | 0.47867 | 1 |
| pilot-2 | historical_biology_hgb | 0.38455 | 0.49907 | 0 |

## Scope and verification

- All four fitted arms use the same continuous ridge head and preprocessing. Structural arms have equal feature count and a common geometry mask.
- Inactive coordinates are a published FGFR1-templated FGFR2 homology model; active coordinates derive from prepared 2PVF. Both retain C491A and lack phosphotyrosine.
- Engineered position 491 and declared uncertain loops are excluded from queries and contact environments. Only common residue/atom identities contribute.
- Static contact/packing differences are mechanistic proxies, not energy, state population, relaxed mutant geometry or simulated catalytic activity.
- Reused development data only; external assay outcomes and original final evaluation remain excluded. No post-score parameter tuning.
- Exact previous ridge reproduction: {'n': 668, 'max_abs_difference': 0.0, 'tolerance': 1e-12}.
- Protected artifacts unchanged: 558.
- Fitting/scoring: 1.133s; evaluator total: 3.893s.
