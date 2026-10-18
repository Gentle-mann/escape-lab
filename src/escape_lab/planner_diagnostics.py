"""Evaluator-produced, identity-free diagnostics allowed into diagnostic planners.

Only a trusted evaluator can establish that aggregates come from its permitted
training partition. Sanitization restricts the shape; it cannot prove provenance.
"""

import math

import numpy as np

REJECTION_CODES = frozenset({
    'undefined_metric', 'no_pooled_ap_improvement', 'pem_ap_regression',
    'fut_ap_regression', 'pem_brier_regression', 'fut_brier_regression',
    'invalid_inner_folds', 'inconsistent_fold_gain', 'worst_fold_regression',
})
FRACTIONS = frozenset({
    'probability_saturated_fraction', 'probability_unique_fraction',
    'feature_shift_fraction', 'nonfinite_score_fraction', 'convergence_warning_fraction',
    'fallback_fraction', 'probability_exact_one_fraction',
})
COUNTS = frozenset({'near_constant_feature_count', 'score_count', 'fit_count', 'below_floor_feature_count'})
MAGNITUDES = frozenset({
    'scaled_feature_abs_max', 'ranking_score_abs_max', 'coefficient_norm_max',
    'coefficient_norm_mean', 'fit_seconds', 'minimum_training_scale',
})


def sanitize_diagnostics(value):
    if not isinstance(value, dict):
        return {}
    result = {}
    for key in sorted(FRACTIONS | COUNTS | MAGNITUDES):
        number = value.get(key)
        if type(number) not in (int, float) or not math.isfinite(number):
            continue
        if key in FRACTIONS and 0 <= number <= 1 or key in COUNTS and type(number) is int and 0 <= number <= 1_000_000_000 or key in MAGNITUDES and 0 <= number <= 1e12:
            result[key] = number
    return result


def rejection_codes(value):
    if not isinstance(value, list):
        return []
    return sorted({code for code in value if isinstance(code, str) and code in REJECTION_CODES})


def numerical_diagnostics(probabilities, *, scaled_features=None, training_variances=None,
                          ranking_scores=None):
    """No rows, names, labels or identities survive this aggregate calculation."""
    scores = np.asarray(probabilities, dtype=float).reshape(-1)
    finite = scores[np.isfinite(scores)]
    result = {'score_count': len(scores),
              'nonfinite_score_fraction': float(1 - len(finite) / len(scores)) if len(scores) else 0.,
              'probability_saturated_fraction': float(np.mean((finite == 0) | (finite == 1)))
              if len(finite) else 0.,
              'probability_unique_fraction': float(len(np.unique(finite)) / len(finite))
              if len(finite) else 0.}
    if scaled_features is not None:
        features = np.asarray(scaled_features, dtype=float)
        known = features[np.isfinite(features)]
        if len(known):
            result['scaled_feature_abs_max'] = float(min(np.max(np.abs(known)), 1e12))
            result['feature_shift_fraction'] = float(np.mean(np.abs(known) > 8.))
    if training_variances is not None:
        variances = np.asarray(training_variances, dtype=float)
        result['near_constant_feature_count'] = int(np.count_nonzero(
            np.isfinite(variances) & (variances >= 0) & (variances < 1e-12)))
    if ranking_scores is not None:
        ranking = np.asarray(ranking_scores, dtype=float)
        known = ranking[np.isfinite(ranking)]
        if len(known):
            result['ranking_score_abs_max'] = float(min(np.max(np.abs(known)), 1e12))
    return sanitize_diagnostics(result)


def aggregate_fit_diagnostics(records):
    """Summarize a fixed number of fits without carrying fold or drug identities."""
    result = {'fit_count': len(records)}
    for key in FRACTIONS:
        values = [sanitize_diagnostics(record).get(key) for record in records]
        known = [value for value in values if value is not None]
        if known:
            result[key] = float(np.mean(known))
    for key in COUNTS - {'fit_count'}:
        values = [sanitize_diagnostics(record).get(key) for record in records]
        known = [value for value in values if value is not None]
        if known:
            # Count of near-constant features is a maximum, not a count of unique columns.
            result[key] = sum(known) if key == 'score_count' else max(known)
    for key in MAGNITUDES:
        values = [sanitize_diagnostics(record).get(key) for record in records]
        known = [value for value in values if value is not None]
        if known:
            result[key] = min(known) if key == 'minimum_training_scale' else max(known)
    norms = [record['coefficient_norm'] for record in records
             if type(record.get('coefficient_norm')) in (int, float)
             and math.isfinite(record['coefficient_norm']) and 0 <= record['coefficient_norm'] <= 1e12]
    if norms:
        result.update(coefficient_norm_mean=float(np.mean(norms)), coefficient_norm_max=max(norms))
    convergence = [record['convergence_warning'] for record in records
                   if type(record.get('convergence_warning')) is bool]
    if convergence:
        result['convergence_warning_fraction'] = float(np.mean(convergence))
    fallback = [record['fallback_reason'] is not None for record in records
                if 'fallback_reason' in record]
    if fallback:
        result['fallback_fraction'] = float(np.mean(fallback))
    return sanitize_diagnostics(result)
