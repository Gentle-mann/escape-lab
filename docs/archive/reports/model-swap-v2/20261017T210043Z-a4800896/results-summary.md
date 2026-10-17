> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# GLM-5.3 model-swap retrospective comparison

This is a single matched model-swap run on reused v2 measurements, not a fresh blind validation. The primary internal endpoint is pooled macro AP; mean fold AP is a separate descriptive endpoint. Higher AP/Spearman and lower Brier are better.

Nemotron reference: `reports/revision-v2/20261017T201321Z-a29b17a8`

GLM run: `reports/model-swap-v2/20261017T210043Z-a4800896`

| Endpoint | v1_fixed | biology_fixed | random | grid | nemotron | glm |
|---|---:|---:|---:|---:|---:|---:|
| Pooled macro AP | 0.034413 | 0.034628 | 0.034938 | 0.034938 | 0.030516 | 0.030516 |
| Macro Brier | 0.035374 | 0.035781 | 0.034593 | 0.034593 | 0.035174 | 0.035174 |
| Mean fold macro AP | 0.040507 | 0.041569 | 0.044647 | 0.044647 | 0.040944 | 0.040944 |
| Fold AP std | 0.011882 | 0.016528 | 0.014032 | 0.014032 | 0.016113 | 0.016113 |
| External macro Spearman | -0.018063 | 0.405681 | 0.283740 | 0.182598 | 0.208915 | 0.296952 |
| Internal pem average_precision | 0.048517 | 0.049223 | 0.047570 | 0.047570 | 0.041003 | 0.041003 |
| Internal pem brier | 0.049245 | 0.050717 | 0.048630 | 0.048630 | 0.049389 | 0.049389 |
| Internal pem n | 3707 | 3707 | 3707 | 3707 | 3707 | 3707 |
| Internal pem positive_count | 175 | 175 | 175 | 175 | 175 | 175 |
| Internal pem inconclusive_count | 23 | 23 | 23 | 23 | 23 | 23 |
| External pem spearman | 0.254280 | 0.298904 | 0.473917 | 0.412320 | 0.103481 | 0.508559 |
| External pem n | 25 | 25 | 25 | 25 | 25 | 25 |
| External pem aligned_groups | 9 | 9 | 9 | 9 | 9 | 9 |
| Internal fut average_precision | 0.020310 | 0.020032 | 0.022307 | 0.022307 | 0.020029 | 0.020029 |
| Internal fut brier | 0.021503 | 0.020845 | 0.020556 | 0.020556 | 0.020959 | 0.020959 |
| Internal fut n | 3695 | 3695 | 3695 | 3695 | 3695 | 3695 |
| Internal fut positive_count | 76 | 76 | 76 | 76 | 76 | 76 |
| Internal fut inconclusive_count | 35 | 35 | 35 | 35 | 35 | 35 |
| External fut spearman | -0.290406 | 0.512458 | 0.093563 | -0.047123 | 0.314349 | 0.085345 |
| External fut n | 25 | 25 | 25 | 25 | 25 | 25 |
| External fut aligned_groups | 9 | 9 | 9 | 9 | 9 | 9 |
| outer-0 macro_ap | 0.054927 | 0.055910 | 0.054186 | 0.054186 | 0.054186 | 0.054186 |
| outer-0 macro_brier | 0.051760 | 0.054225 | 0.051343 | 0.051343 | 0.051343 | 0.051343 |
| outer-0 pem average_precision | 0.073484 | 0.079489 | 0.072738 | 0.072738 | 0.072738 | 0.072738 |
| outer-0 pem brier | 0.070230 | 0.075215 | 0.069827 | 0.069827 | 0.069827 | 0.069827 |
| outer-0 fut average_precision | 0.036370 | 0.032331 | 0.035634 | 0.035634 | 0.035634 | 0.035634 |
| outer-0 fut brier | 0.033291 | 0.033236 | 0.032860 | 0.032860 | 0.032860 | 0.032860 |
| outer-1 macro_ap | 0.040770 | 0.050383 | 0.054947 | 0.054947 | 0.050383 | 0.050383 |
| outer-1 macro_brier | 0.027595 | 0.025541 | 0.024697 | 0.024697 | 0.025541 | 0.025541 |
| outer-1 pem average_precision | 0.061685 | 0.047438 | 0.065502 | 0.065502 | 0.047438 | 0.047438 |
| outer-1 pem brier | 0.044138 | 0.041385 | 0.039936 | 0.039936 | 0.041385 | 0.041385 |
| outer-1 fut average_precision | 0.019855 | 0.053329 | 0.044392 | 0.044392 | 0.053329 | 0.053329 |
| outer-1 fut brier | 0.011052 | 0.009697 | 0.009457 | 0.009457 | 0.009697 | 0.009697 |
| outer-2 macro_ap | 0.025825 | 0.018413 | 0.024807 | 0.024807 | 0.018264 | 0.018264 |
| outer-2 macro_brier | 0.024230 | 0.024227 | 0.024506 | 0.024506 | 0.025493 | 0.025493 |
| outer-2 pem average_precision | 0.034650 | 0.023861 | 0.029644 | 0.029644 | 0.021717 | 0.021717 |
| outer-2 pem brier | 0.031640 | 0.032423 | 0.033217 | 0.033217 | 0.034274 | 0.034274 |
| outer-2 fut average_precision | 0.017000 | 0.012964 | 0.019971 | 0.019971 | 0.014810 | 0.014810 |
| outer-2 fut brier | 0.016819 | 0.016032 | 0.015796 | 0.015796 | 0.016713 | 0.016713 |

Selected configurations:

| Method | outer-0 | outer-1 | outer-2 | Full search |
|---|---|---|---|---|
| v1_fixed | v1_fixed | v1_fixed | v1_fixed | v1_fixed |
| biology_fixed | structure-logistic-c0.1 | structure-logistic-c0.1 | structure-logistic-c0.1 | structure-logistic-c0.1 |
| random | chemistry-logistic-c0.1 | structure-hgb-depth3 | structure-hgb-depth3 | interactions-hgb-depth2 |
| grid | chemistry-logistic-c0.1 | structure-hgb-depth3 | structure-hgb-depth3 | structure-hgb-depth3 |
| nemotron | chemistry-logistic-c0.1 | structure-logistic-c0.1 | chemistry-logistic-c0.1 | structure-logistic-c1 |
| glm | chemistry-logistic-c0.1 | structure-logistic-c0.1 | chemistry-logistic-c0.1 | interactions-logistic-c0.1 |

Observed agent budgets and usage:

| Run | Valid / opportunities | Requests / repairs | Input / output tokens | Budget tokens | Seconds | Estimated USD |
|---|---:|---:|---:|---:|---:|---:|
| nemotron | 12 / 12 | 14 / 2 | 95604 / 53354 | 148958 | 243.551 | 0.030902 |
| glm | 6 / 12 | 17 / 5 | 84687 / 80813 | 165500 | 393.632 | 0.474139 |

nemotron search statuses (outer-0, outer-1, outer-2, full): `["completed", "completed", "completed", "completed"]`. Proposal statuses: `{"proposed": 12}`. Failures: `{}`. Failed or truncated responses consumed their original opportunity; the recorded format-repair policy was retained.

glm search statuses (outer-0, outer-1, outer-2, full): `["degraded", "degraded", "degraded", "degraded"]`. Proposal statuses: `{"proposed": 6, "truncated_response": 4, "invalid_response": 2}`. Failures: `{"The response reached the completion token limit.": 4, "No valid catalog experiment was proposed after one format repair.": 2}`. Failed or truncated responses consumed their original opportunity; the recorded format-repair policy was retained.

Costs use uncached Crusoe list prices retrieved 2026-10-17: Nemotron $0.10 input / $0.40 output and GLM $1.40 input / $4.40 output per million tokens. These are estimates, not invoices; caching discounts are not applied. [Pricing source](https://www.crusoe.ai/cloud/pricing).

Paired descriptive 95% block-bootstrap percentile intervals (2,000 identical group draws, seed 6100), conditional on committed predictions:

| Endpoint | Contrast | Delta | 2.5% | 97.5% | Valid / invalid draws |
|---|---|---:|---:|---:|---:|
| internal | glm_minus_v1_fixed | -0.003897 | -0.012646 | 0.003915 | 2000 / 0 |
| internal | glm_minus_biology_fixed | -0.004112 | -0.010859 | 0.001459 | 2000 / 0 |
| internal | glm_minus_random | -0.004422 | -0.011090 | 0.000055 | 2000 / 0 |
| internal | glm_minus_grid | -0.004422 | -0.011090 | 0.000055 | 2000 / 0 |
| internal | glm_minus_nemotron | 0.000000 | 0.000000 | 0.000000 | 2000 / 0 |
| external | glm_minus_v1_fixed | 0.315015 | -0.183845 | 0.601315 | 2000 / 0 |
| external | glm_minus_biology_fixed | -0.108729 | -0.519605 | 0.082645 | 2000 / 0 |
| external | glm_minus_random | 0.013212 | -0.483642 | 0.387008 | 2000 / 0 |
| external | glm_minus_grid | 0.114354 | -0.302247 | 0.460094 | 2000 / 0 |
| external | glm_minus_nemotron | 0.088037 | -0.358545 | 0.368979 | 2000 / 0 |

The four numerical baselines, their search choices, partitions, committed predictions and metrics were checked for agreement. All internal metrics and external per-drug Spearman/N were recomputed from released predictions. No refitting, fresh holdout access, outcome-based selection of runs, or statistical significance claim is made. Intervals omit fitting and LLM selection uncertainty; related blocks are not independent biological replicates. One run cannot establish a general model-quality effect.

The original protocol fields changed only in model ID and creation time. The comparison adds its provenance and retrospective authorization. Agent instructions, tools, numerical evaluator, folds, data, features, seeds, selection guards, opportunity/repair budgets and scientific limits retain the reference protocol.

Explicit request settings: `{"tool_choice": "required", "max_tokens": 8192, "stream": false}`. Unspecified in both: temperature, top_p, seed, reasoning_effort, thinking. Model-specific sampling and reasoning defaults are not controlled; remaining_seconds and subsequent feedback histories are dynamic.

Runtime is measured local wall time, including instrumentation and any concurrent local work; it is not a controlled serving-speed benchmark.


Reproduce analysis after scientific completion:

```sh
uv run python -m escape_lab.model_swap_analysis 'reports/revision-v2/20261017T201321Z-a29b17a8' 'reports/model-swap-v2/20261017T210043Z-a4800896'
```

Exact inputs, hashes, counts, all contrasts and full metrics are in `comparison.json`. Existing scientific artifacts are not rewritten.
