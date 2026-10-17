> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# Tiny matched tree-model pilot

Run only six small HGB fits: baseline and one structural challenger on each of
the three saved outer folds. This isolates the activation-prediction component;
omit agent inference, resistance models, nested crossfits and hyperparameter
search. Preserve the difficult third fold rather than choosing a favorable
slice of already-reused development data.

The exact historical biology baseline trained on 2,719 complete pooled FGFR1–4
variants. Preserve its original candidate order, global aligned-block folds,
HGB parameters, original biological features, median imputation and target.
Training fold sizes are 1,819 / 1,745 / 1,874. Score only the same 668 FGFR2
variants, split 275 / 180 / 213. Reproduce all historical FGFR2 baseline scores
to 1e-12 before fitting the challenger. Use cached authorized development
targets; no workbook ingestion is necessary.

The challenger adds exactly the 31 existing frozen active-minus-inactive
structural descriptors for FGFR2. Other genes have zero descriptors and
availability=0; do not invent cross-gene geometry. Keep HGB capacity identical.
No ESM, PCA, scaling, range clamp, output clipping, new structures or new
biological features. Trees do not perform the unsafe standardized linear
extrapolation diagnosed in the previous ridge run.

Keep the previous development screening requirements: mean-fold Spearman gain
at least 0.03, improvement in at least two folds, mean-fold MAE no more than 5%
worse, and no loss of high-activation hits at the frozen 68-nomination budget.
High/strong thresholds and continuous target are unchanged. Use 1,000 paired
block resamples within folds, seed 17, as descriptive uncertainty only.

Report all six fits, every fold and failed criterion. Do not tune after scoring.
This small test can identify a developmental lead; it cannot establish external
generalization or a molecular mechanism. External RLY2023/MANO2021 assay data
and the original final evaluation remain unopened/excluded.
