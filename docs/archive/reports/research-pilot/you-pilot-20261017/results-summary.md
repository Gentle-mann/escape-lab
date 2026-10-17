> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# Reviewed You.com research pilot

Status: **completed**. One paired GLM run, three opportunities per arm.

Both arms used the original fixed v2 numerical catalog and features. Only the research arm received reviewed general biological evidence cards.

External research-minus-control correlation difference: **0.00000**; paired descriptive 95% interval [0.00000, 0.00000].

| Arm / method | Selected configuration | Inner selection AP | Pem correlation | Fut correlation | Mean correlation |
|---|---|---:|---:|---:|---:|
| control / v1_fixed | v1_fixed | n/a | 0.25428 | -0.29041 | -0.01806 |
| control / biology_fixed | structure-logistic-c0.1 | 0.04178 | 0.29890 | 0.51246 | 0.40568 |
| control / random | interactions-hgb-depth2 | 0.06027 | 0.47392 | 0.09356 | 0.28374 |
| control / grid | structure-hgb-depth3 | 0.06053 | 0.41232 | -0.04712 | 0.18260 |
| control / agent | interactions-logistic-c0.1 | 0.04217 | 0.50856 | 0.08535 | 0.29695 |
| research / agent | interactions-logistic-c0.1 | 0.04217 | 0.50856 | 0.08535 | 0.29695 |

Inner AP is a selection diagnostic, not an independent internal test. The four numerical controls reproduced original selections, folds, inner scores and committed external predictions to absolute tolerance 1e-12.

| Arm | Valid proposals / opportunities | Requests | Charged tokens | Seconds | Estimated inference USD |
|---|---:|---:|---:|---:|---:|
| control | 2 / 3 | 3 | 33899 | 86.26 | 0.09844 |
| research | 2 / 3 | 3 | 33960 | 80.12 | 0.09783 |

Inference cost uses previously recorded Crusoe GLM prices of $1.40/M input and $4.40/M output tokens, without cache discounts. It is an estimate, not an invoice; retrieval cost is accounted separately in research artifacts.

Both selected models and predictions were frozen before released S13 scoring.

This is retrospective reuse of an already spent small biochemical panel. One stochastic pair cannot establish a causal or average benefit from research. Bootstrap intervals condition on these fitted predictions and exclude model, literature-selection and biological-replicate uncertainty. No fresh external holdout was consumed.

Run directory: `reports/research-pilot/you-pilot-20261017`.

Reproduction:

```sh
cd '.'
uv run python -m escape_lab.research_pilot --root . --evidence 'research/you-evidence/you-pilot-20261017/evidence.json'
```

Each invocation creates a new recorded scientific attempt; do not rerun to select a favorable outcome.
