import copy

import numpy as np
import pytest
from sklearn.metrics import average_precision_score

from escape_lab.research_repeat_analysis import METHODS, analyze
from escape_lab.revision_evaluation import summarize


def records(directions=None, *, missing_drug=None):
    directions = directions or {}
    rows = []
    for repeat in range(10):
        for method in METHODS:
            direction = directions.get(method, 0)
            for group in range(4):
                for variant in range(4):
                    label = variant >= 2
                    # Every complete group contains both classes; a small repeat effect
                    # changes ranking in the neutral condition without changing truth.
                    neutral = [.15, .70, .40, .80][(variant + repeat % 2) % 4]
                    score = (.9 if label else .1) if direction == 1 else (
                        (.1 if label else .9) if direction == -1 else neutral)
                    for drug in ("pem", "fut"):
                        rows.append({
                            "repeat": repeat, "method": method,
                            "candidate_id": f"v{group}-{variant}", "group_id": f"g{group}",
                            "gene": "FGFR2", "drug": drug, "fold": f"outer-{group % 3}",
                            "probability": score,
                            "label": None if drug == missing_drug else label,
                        })
    return rows


def test_identical_predictions_have_zero_deltas_and_intervals():
    result = analyze(records(), replicates=25, seed=81)
    for contrast in result["comparisons"].values():
        assert contrast["delta"] == 0
        assert contrast["interval_95"] == [0, 0]
        assert contrast["interval_99_bonferroni"] == [0, 0]
        assert contrast["repeat_only"]["ties"] == 10
        assert contrast["repeat_only"]["delta_sd"] == 0
    assert result["experiments"]["research"]["status"] == "no_observed_improvement"
    assert result["experiments"]["agent_selection"]["status"] == "no_observed_improvement"
    assert result == analyze(list(reversed(records())), replicates=25, seed=81)


def test_paired_group_and_repeat_draws_match_independent_brute_calculation():
    rows = records({"relevant": 1})
    result = analyze(rows, replicates=20, seed=81)
    rng = np.random.default_rng(81)
    by_cell = {(repeat, method): [r for r in rows if r["repeat"] == repeat
                                and r["method"] == method]
               for repeat in range(10) for method in ("relevant", "control")}
    differences = []
    for _ in range(20):
        group_draw = rng.integers(4, size=4)
        repeat_draw = rng.integers(10, size=10)
        means = {}
        for method in ("relevant", "control"):
            scores = []
            for repeat in repeat_draw:
                sampled = [r for group in group_draw for r in by_cell[(repeat, method)]
                           if r["group_id"] == f"g{group}"]
                scores.append(summarize(sampled)["macro_ap"])
            means[method] = np.mean(scores)
        differences.append(means["relevant"] - means["control"])
    comparison = result["comparisons"]["relevant_minus_control"]
    assert comparison["interval_95"] == pytest.approx(np.percentile(differences, [2.5, 97.5]))
    assert comparison["interval_99_bonferroni"] == pytest.approx(
        np.percentile(differences, [.5, 99.5]))
    assert result["bootstrap"]["repeat_draws_shared_across_methods"] is True
    assert result["bootstrap"]["group_draws_shared_across_methods"] is True


@pytest.mark.parametrize("field,value", [("label", True), ("group_id", "wrong"),
                                       ("fold", "outer-99"), ("gene", "wrong")])
def test_paired_truth_and_metadata_must_match(field, value):
    rows = records()
    rows[0][field] = value
    with pytest.raises(ValueError, match="paired|metadata"):
        analyze(rows, replicates=2)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -.1, 1.1, True])
def test_invalid_probabilities_rejected(value):
    rows = records()
    rows[0]["probability"] = value
    with pytest.raises(ValueError, match="probability"):
        analyze(rows, replicates=2)


def test_coverage_missing_repeat_duplicate_and_wrong_label_rejected():
    source = records()
    for rows in (source[1:], source + [source[0]], [r for r in source if r["repeat"] != 9]):
        with pytest.raises(ValueError, match="coverage|Duplicate|repeat"):
            analyze(rows, replicates=2)
    rows = copy.deepcopy(source)
    rows[0]["label"] = 1
    with pytest.raises(ValueError, match="label"):
        analyze(rows, replicates=2)


def test_macro_ap_is_mean_of_repeat_pooled_oof_scores_with_ties():
    rows = records({"relevant": 1})
    result = analyze(rows, replicates=5)
    expected = []
    for repeat in range(10):
        subset = [r for r in rows if r["repeat"] == repeat and r["method"] == "control"]
        metrics = summarize(subset)
        expected.append(metrics["macro_ap"])
        assert result["methods"]["control"]["per_repeat"][str(repeat)] == metrics
    assert result["methods"]["control"]["mean"]["macro_ap"] == pytest.approx(np.mean(expected))
    assert result["methods"]["relevant"]["mean"]["macro_ap"] == 1
    assert result["methods"]["control"]["mean"]["by_fold"]


def test_weighted_ap_tie_handling_matches_sklearn():
    from escape_lab.research_repeat_analysis import _weighted_ap_draws

    truth = np.array([0, 1, 1, 0, 1, 0], dtype=float)
    scores = np.array([.5, .5, .7, .7, .8, .1])
    groups = np.array([0, 1, 0, 1, 2, 2])
    counts = np.array([[1, 1, 1], [2, 0, 1], [0, 3, 0]])
    observed = _weighted_ap_draws(truth, scores, groups, counts)
    expected = [average_precision_score(truth, scores, sample_weight=c[groups]) for c in counts]
    assert observed == pytest.approx(expected)


def test_positive_and_negative_experiment_gates():
    research = analyze(records({"relevant": 1}), replicates=100)
    assert research["experiments"]["research"]["status"] == "supported_development_improvement"
    agent = analyze(records({"control": 1}), replicates=100)
    assert agent["experiments"]["agent_selection"]["status"] == (
        "supported_development_improvement")
    worse = analyze(records({"relevant": -1, "control": -1}), replicates=100)
    assert worse["experiments"]["research"]["status"] == "no_observed_improvement"
    assert worse["experiments"]["agent_selection"]["status"] == "no_observed_improvement"


def test_positive_point_without_both_research_contrasts_is_inconclusive():
    result = analyze(records({"relevant": 1, "placebo": 1}), replicates=100)
    assert result["comparisons"]["relevant_minus_control"]["passes_improvement_rule"]
    assert not result["comparisons"]["relevant_minus_placebo"]["passes_improvement_rule"]
    assert result["experiments"]["research"]["status"] == "inconclusive"


def test_undefined_metrics_are_inconclusive_without_dropping_them():
    result = analyze(records(missing_drug="fut"), replicates=10)
    assert result["methods"]["control"]["mean"]["macro_ap"] is None
    contrast = result["comparisons"]["relevant_minus_control"]
    assert contrast["delta"] is None
    assert contrast["interval_99_bonferroni"] is None
    assert contrast["invalid_replicates"] == 10
    assert result["experiments"]["research"]["status"] == "inconclusive"


def test_grid_is_optional_but_any_supplied_grid_requires_full_coverage():
    rows = [r for r in records() if r["method"] != "grid"]
    assert set(analyze(rows, replicates=2)["methods"]) == set(METHODS) - {"grid"}
    rows.append(next(r for r in records() if r["method"] == "grid"))
    with pytest.raises(ValueError, match="coverage"):
        analyze(rows, replicates=2)


def test_positive_ap_gain_does_not_override_calibration_guard():
    rows = records({"relevant": 1})
    for row in rows:
        if row["method"] == "relevant":
            row["probability"] = .99 if row["label"] else .98
    result = analyze(rows, replicates=100)
    comparison = result["comparisons"]["relevant_minus_control"]
    assert comparison["delta"] > .002
    assert comparison["interval_99_bonferroni"][0] > 0
    assert "pem_brier_regression" in comparison["failed_rules"]
    assert "fut_brier_regression" in comparison["failed_rules"]
    assert result["experiments"]["research"]["status"] == "inconclusive"


def test_positive_macro_gain_does_not_override_per_drug_regression():
    rows = records({"relevant": 1})
    for row in rows:
        if row["method"] == "relevant" and row["drug"] == "fut":
            variant = int(row["candidate_id"].split("-")[-1])
            row["probability"] = [.8, .6, .9, .5][variant]
    result = analyze(rows, replicates=100)
    comparison = result["comparisons"]["relevant_minus_control"]
    assert comparison["delta"] > .002
    assert "fut_ap_regression" in comparison["failed_rules"]
    assert result["experiments"]["research"]["status"] == "inconclusive"


def test_some_undefined_bootstrap_draws_do_not_produce_valid_only_interval():
    rows = records({"relevant": 1})
    for row in rows:
        if row["drug"] == "fut" and row["group_id"] != "g0":
            row["label"] = None
    result = analyze(rows, replicates=100)
    comparison = result["comparisons"]["relevant_minus_control"]
    assert comparison["delta"] is not None
    assert 0 < comparison["invalid_replicates"] < 100
    assert comparison["interval_95"] is None
    assert comparison["interval_99_bonferroni"] is None
    assert result["experiments"]["research"]["status"] == "inconclusive"
