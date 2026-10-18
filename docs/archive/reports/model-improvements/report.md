# FGFR2 model improvements

A fixed 12-fit experiment tests numerical stability and identifies which sequence representations contribute to activation ranking.

The same 668 FGFR2 variants, ten position blocks and three saved development folds are used throughout. Each fit uses only its training partition.

| Arm | Mean-fold AP | High hits / 68 picks | Maximum absolute margin |
|---|---:|---:|---:|
| legacy_both | 0.111346 | 4 | 1.084e+06 |
| floor_both | 0.059510 | 2 | 20.13 |
| legacy_wt | 0.114019 | 4 | 9.094e+05 |
| legacy_delta | 0.113766 | 5 | 2.198e+06 |

Stability advancement gate: **fail**.

The historical probability and raw-margin predictions are reproduced before any challenger is fitted. Ranking uses raw margins in every arm.

The variance floor changes only biology feature scaling: fitted standard deviations below 0.01 use a denominator of 0.01. Sequence scaling and PCA stay unchanged. The WT-only and mutation-delta-only arms keep legacy scaling to isolate sequence content.

Representation arms are diagnostic ablations. They do not change the default model automatically. Reused development folds measure this experiment; an external biological validation remains the next test.

Fits: 12; fitting/scoring: 0.919s; total: 2.133s.

Protected artifacts: 529; changed: 0. No final or external evaluation was run.

## Predeclared decision

```json
{
  "stability": {
    "status": "fail",
    "checks": {
      "mean_fold_ap_gain_at_least_001": false,
      "no_fold_ap_loss": false,
      "no_high_hit_loss": false,
      "maximum_absolute_margin_reduced_tenfold": true
    },
    "mean_fold_ap_gain": -0.05183592287423172,
    "high_hits_gain": -2,
    "margin_max_before": 1084487.1128465498,
    "margin_max_after": 20.129137691660976
  },
  "representation": {
    "status": "diagnostic_only",
    "default_changed": false,
    "arms": [
      "legacy_both",
      "legacy_wt",
      "legacy_delta"
    ]
  }
}
```

## What this establishes

- Both historical probability scores and raw margins reproduce exactly across all 668 variants. The shared ranking API preserves the earlier AP improvement without changing the fitted head.
- The single predeclared variance floor removes all 136 saturated probabilities in the affected fold and reduces the largest absolute margin from 1,084,487 to 20.13. Its AP falls from 0.111346 to 0.059510 and retrieval falls from four to two hits. It remains an optional rejected challenger.
- WT-only AP is 0.114019 with four hits; mutation-delta-only AP is 0.113766 with five hits. Both differences from the two-branch model are small. These arms retain extreme extrapolation, so their numbers do not establish improved biological generalization.
- Combining the two sequence branches does not improve this fixed development comparison. The next research question is how to preserve mutation information while learning a stable relationship across positions. No new threshold or configuration was selected after scoring.

## Reproduce this experiment

With the existing approved local assay inputs and frozen ESM2 cache:

```sh
uv run --frozen python -m escape_lab.fgfr2_sequence_revision --output reports/model-improvements/a-fresh-run
```

The runner writes its protocol before loading approved development outcomes, reproduces the historical baseline before fitting challengers, keeps row-level predictions in ignored local reports, and checks the protected input hashes after scoring. This archive contains aggregate results and a sanitized protocol only.
