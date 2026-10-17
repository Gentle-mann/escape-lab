> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# FGFR2 sequence representation pilot

Primary development gate: **fail**.

668 FGFR2 variants; 26 high-activation positives; ten occupied position blocks.
High means >=3/4 replicate medians >=1.5 AND median >5. This is a protein-level adaptation of the study categories. Strong (>10) is a secondary endpoint.

| Model | Mean fold AP | High hits / 68 picks | Recall of 26 | Strong hits / 7 |
|---|---:|---:|---:|---:|
| cached_v2_biology | 0.0719 | 3 | 11.5% | 0 |
| biology_logistic | 0.0471 | 0 | 0.0% | 0 |
| esm_biology_logistic | 0.0973 | 1 | 3.8% | 0 |
| esm_biology_ridge | 0.1171 | 4 | 15.4% | 0 |

AP is ranking quality, not accuracy. Selections are top 10% within each fold (28+18+22=68), with candidate-ID tie breaking. Pooled scores are secondary.

```json
{
  "status": "fail",
  "checks": {
    "three_extra_hits_vs_v2": false,
    "ap_gain_003_vs_v2": false,
    "ap_gain_003_vs_biology": true,
    "no_hit_loss_vs_biology": true,
    "ap_improves_two_folds_vs_v2": false
  },
  "high_hits_gain_vs_v2": -2,
  "mean_fold_ap_gain_vs_v2": 0.02545734538724327,
  "mean_fold_ap_gain_vs_biology": 0.05022878500373031,
  "fold_wins": 1
}
```

## Per-fold results

| Fold | Model | AP | High hits | Positive / N |
|---|---|---:|---:|---:|
| pilot-0 | cached_v2_biology | 0.0453 | 2 | 7 / 275 |
| pilot-1 | cached_v2_biology | 0.1133 | 1 | 7 / 180 |
| pilot-2 | cached_v2_biology | 0.0571 | 0 | 12 / 213 |
| pilot-0 | biology_logistic | 0.0246 | 0 | 7 / 275 |
| pilot-1 | biology_logistic | 0.0579 | 0 | 7 / 180 |
| pilot-2 | biology_logistic | 0.0587 | 0 | 12 / 213 |
| pilot-0 | esm_biology_logistic | 0.0246 | 0 | 7 / 275 |
| pilot-1 | esm_biology_logistic | 0.2157 | 1 | 7 / 180 |
| pilot-2 | esm_biology_logistic | 0.0516 | 0 | 12 / 213 |
| pilot-0 | esm_biology_ridge | 0.0718 | 3 | 7 / 275 |
| pilot-1 | esm_biology_ridge | 0.2234 | 1 | 7 / 180 |
| pilot-2 | esm_biology_ridge | 0.0560 | 0 | 12 / 213 |

Head fitting/scoring: 0.64s; feature extraction: 49.62s.

- ESM2-35M remains frozen. No LLM calls or supervised PLM fine-tuning.
- Scalers, imputers and separate PCA16 branches fit only each training fold.
- Class-balanced logistic outputs are ranking scores, not calibrated probabilities.
- Cached V2 reference was trained across genes. The matched biology-only classifier isolates the representation change in FGFR2-only training.
- No external assay outcomes, original final evaluation, resistance labels or post-score hyperparameter tuning used. These development folds have been reused.
- Strong activation has no positives in pilot-1; its AP there is undefined.
- All 518 protected artifacts hash-identical.

See protocol.json, results.json and feature manifest for full provenance.
