import copy

import pytest

from escape_lab.fgfr2_state_run import (
    compare_gates,
    continuous_summary,
    paired_uncertainty,
    verify_cached_predictions,
)


def rows(scores=(0., 1., 2.), truth=(0., 1., 2.)):
    return [{"candidate_id": str(i), "group_id": f"g{i}", "fold": "pilot-0",
             "score": p, "truth": y, "high": i == 2, "strong": False}
            for i, (p, y) in enumerate(zip(scores, truth, strict=True))]


def test_continuous_metric_direction_and_no_clipping():
    result = continuous_summary(rows((-1., 1., 4.)))
    assert result["mean_fold_spearman"] == 1.
    assert result["mean_fold_mae"] == 1.
    assert result["high_hits"] == 1
    assert result["pooled_mae"] == 1.


def test_mean_fold_not_pooled_metric():
    data = rows()
    other = [{**r, "candidate_id": "x" + r["candidate_id"], "fold": "pilot-1",
              "score": 2. - r["score"]} for r in rows()]
    result = continuous_summary(data + other)
    assert result["mean_fold_spearman"] == 0.
    assert result["mean_fold_mae"] == pytest.approx(2 / 3)


def metric(rho, mae=1., hits=4, folds=(.2, .2, .2)):
    return {"mean_fold_spearman": rho, "mean_fold_mae": mae, "high_hits": hits,
            "by_fold": {str(i): {"spearman": r} for i, r in enumerate(folds)}}


def test_two_gates_distinguish_general_gain_from_state_contrast():
    values = {"baseline": metric(.20), "inactive": metric(.24),
              "active": metric(.22), "contrast": metric(.23, folds=(.23, .23, .23))}
    result = compare_gates(values)
    assert result["incremental"]["status"] == "pass"
    assert result["state_contrast"]["status"] == "fail"
    values["contrast"] = metric(.25, 1.051, folds=(.25, .25, .25))
    assert compare_gates(values)["incremental"]["status"] == "fail"
    values["contrast"] = metric(.25, 1.0, 3, folds=(.25, .25, .25))
    assert compare_gates(values)["incremental"]["status"] == "fail"
    values["contrast"]["mean_fold_spearman"] = None
    assert compare_gates(values)["incremental"]["status"] == "inconclusive"


def test_cached_reproduction_checks_partition_truth_and_prediction():
    original = rows()
    assert verify_cached_predictions(original, list(reversed(original)))["max_abs_difference"] == 0.
    for key, value in (("score", 999.), ("truth", 999.), ("fold", "wrong"), ("group_id", "wrong")):
        bad = copy.deepcopy(original)
        bad[0][key] = value
        with pytest.raises(ValueError):
            verify_cached_predictions(bad, original)


def test_duplicate_or_nonfinite_predictions_rejected():
    with pytest.raises(ValueError):
        verify_cached_predictions(rows() + rows(), rows())
    with pytest.raises(ValueError):
        continuous_summary(rows((float("nan"), 1., 2.)))


def test_constant_predictions_do_not_silently_pass_correlation_gate():
    assert continuous_summary(rows((1., 1., 1.)))["mean_fold_spearman"] is None


def test_paired_bootstrap_is_invariant_to_row_order_and_identical_predictions():
    data = [{**r, "group_id": "g"} for r in rows()]
    arms = {arm: copy.deepcopy(data) for arm in ("baseline", "inactive", "active", "contrast")}
    arms["active"].reverse()
    result = paired_uncertainty(arms, draws=10)
    for contrast in result["contrasts"].values():
        for metric in contrast.values():
            assert metric["valid_draws"] == 10
            assert metric["invalid_draws"] == 0
            assert metric["percentile_2_5"] == 0.
            assert metric["percentile_97_5"] == 0.
    arms["active"][0]["truth"] += 1.
    with pytest.raises(ValueError, match="Unpaired"):
        paired_uncertainty(arms, draws=10)


def test_paired_bootstrap_counts_undefined_correlations():
    arms = {arm: rows((1., 1., 1.)) for arm in ("baseline", "inactive", "active", "contrast")}
    result = paired_uncertainty(arms, draws=10)
    for contrast in result["contrasts"].values():
        assert contrast["rho_gain"]["valid_draws"] == 0
        assert contrast["rho_gain"]["invalid_draws"] == 10
        assert contrast["rho_gain"]["percentile_2_5"] is None
        assert contrast["mae_reduction"]["valid_draws"] == 10
