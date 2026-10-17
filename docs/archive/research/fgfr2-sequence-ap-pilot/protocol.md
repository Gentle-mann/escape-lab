> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# One targeted AP test: recover classifier ranking before sigmoid saturation

The frozen sequence-plus-biology classifier's mean-fold AP is 0.09732 versus
0.04709 for the matched biology-only classifier. This is dominated by pilot-1:
0.21570 versus 0.05792. Pilot-0 is effectively unchanged (0.02463), and pilot-2
is worse (0.05163 versus 0.05872). Therefore the existing ablation supports a
sequence-representation effect mainly on one fold; it does not demonstrate
consistent generalization. Comparison with the pooled tree additionally changes
training objective and cohort, so that comparison does not isolate sequence.

Probability diagnostics reveal 136 of 275 sequence-classifier outputs equal
exactly 1.0 in pilot-0. This destroys distinctions between large finite logits.
The matched biology-only classifier has a similar 135-row saturation tie. The
previous standardized-feature extrapolation diagnosis explains how very large
classifier margins can arise. The current targeted test asks whether preserving
raw decision margins recovers useful ranking within those numerical ties.

Run exactly three unchanged sequence-classifier fits on the original FGFR2
cohort/folds. Retain the same ESM cache, training-only PCA/scalers, class balance,
logistic settings and labels. Each fitted model produces original probabilities
and raw margins. Require all 668 original probabilities to reproduce to 1e-12
before accepting margin metrics. Do not change parameters, features, targets,
training cohort, output calibration, or clip extrapolated inputs.

This is a technical ranking sensitivity, not a new biological model. Its fixed
screen requires mean-fold AP gain at least 0.01, no fold AP degradation, and no
loss of high-activation hits at the existing 68-pick budget. A two-fold gain is
not expected because only the saturated fold should materially change. Preserve
failure if breaking the tie worsens ranking. Use 1,000 paired block resamples,
seed17, as descriptive uncertainty. Report every fold and the unchanged tree
reference; do not select a different endpoint after scoring.

Only previously authorized development activation values may be read from the
original assay. Reserved external assays and final outcomes remain excluded.
All prior files and these source/test/protocol files are hash-frozen before fits.
