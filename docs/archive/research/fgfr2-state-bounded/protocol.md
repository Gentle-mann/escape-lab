> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# Numerical-stability sensitivity after the preserved FGFR2 state pilot

This is a post-result development sensitivity check, not an independent test.
The original `fgfr2-state-v1` failed both gates and is preserved without changes.
Its unbounded ridge baseline reproduced all previous scores exactly. Continuous
MAE exposed a failure that ranking metrics alone had obscured.

An outcome-free feature-matrix diagnosis traced the largest extrapolation to
`activation_position_rbf_450`: training standard deviation 1.171e-8 and held-out
maximum 0.923 produced a standardized value near 78.8 million. Related proximity
features also extrapolated far outside their training ranges. This is not a
floating-point crash or evidence about molecular activation; it is unsupported
linear extrapolation under the held-out-position evaluation.

## One fixed correction

Keep each model's training matrix, ridge coefficients, alpha, feature set,
sequence PCA, target, cohort, folds and evaluation gates unchanged. After the
training-only biology scaler, record every biology column's training minimum
and maximum. At inference, clamp those standardized biology columns to that
training range. This applies to the same old biology features in every arm and
to the additional structural features in each structural arm. Sequence PCA
columns remain unchanged. No output clipping or label-dependent thresholds.

This deliberately conservative policy avoids extrapolating localized features
into regions absent from training. It can suppress useful out-of-range signal;
it is not a calibrated uncertainty estimate or proof that interpolation is
biologically correct. There is no search over clipping values or other models.

For each of all four arms, require the same fitted model's **unbounded** scores
to reproduce all 668 saved original-run predictions to 1e-12. Only then accept
the separately computed bounded scores. Report clipping counts, all original
metrics/gates and the unchanged historical HGB reference. Use the same fixed
1,000 paired block resamples for descriptive uncertainty.

## Interpretation boundary

Any gain over the previous unbounded run may come from this preprocessing
correction. Evidence for added structural information must compare the bounded
contrast arm against the equally bounded baseline and single-state controls.
Both attempts reuse development outcomes. External RLY2023/MANO2021 assays and
the original final evaluation remain unopened/excluded. Preserve source hashes
and both results; do not relabel a favorable sensitivity result as a successful
preregistered experiment.
