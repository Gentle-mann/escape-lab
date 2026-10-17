> Historical research record. See [current results](../../../../RESULTS.md) for the release interpretation.

# Decision from the tiny matched tree pilot

Do not promote the current active-minus-inactive descriptor block. It worsened
FGFR2 activation prediction in this test, including ranking and absolute error
on all three existing folds. Keep the historical biology tree as the best
measured baseline for now.

| Measure | Exact historical tree baseline | Same tree + frozen contrast |
|---|---:|---:|
| Mean-fold Spearman | 0.38350 | 0.32900 |
| Mean-fold MAE | 0.49053 | 0.54613 |
| High activators found, 68 nominations | 3 of 26 | 1 of 26 |
| Strong activators found | 0 of 7 | 0 of 7 |

Six small fits took 5.28 seconds; the complete evaluator took 6.67 seconds.
Both models used identical training data/order, fixed folds, continuous targets,
HGB settings and preprocessing. Only the 31 cached FGFR2 contrast features were
added. All 668 baseline predictions exactly reproduced the historical results.
No parameter search, additional structures, model calls, assay workbook reads,
external data outcomes or final evaluation were used.

The unchanged training setup pools FGFR1–4; the new structural coordinates cover
FGFR2 only, with explicit zeros and an availability flag for unsupported rows.
This is a test of that specific feature implementation and shared model, not a
refutation of active/inactive-state biology or thermodynamic simulation.

The correlation-gain block-resampling interval spans -0.167 to +0.043, so this
small reused-development sample is uncertain about the wider population effect.
Nevertheless, every measured fold deteriorated, all preset screening criteria
failed, and the observed MAE increased by about 11.3 percent. There is no positive
result to carry forward from this comparison.

With the hackathon deadline close, stop further tuning of these static features
against the same reused outcomes. Preserve the stronger baseline and the visible
rejection of this attempted improvement. The external datasets remain reserved.

See report.md, protocol.json, results.json and independent-verification.json for
reproducible details. The relevant regression checks passed 52 tests and scoped
Ruff; all 583 protected files stayed unchanged at completion.
