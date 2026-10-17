import pytest

from escape_lab.revision_evaluation import (
    accepts,
    grouped_folds,
    retained_candidates,
    summarize,
)


def metric(ap=.2, brier=.1, folds=(.2, .2, .2)):
    return {"macro_ap": ap, "per_drug": {
        d: {"average_precision": ap, "brier": brier} for d in ("pem", "fut")},
        "by_fold": {f"inner-{i}": {"macro_ap": value} for i, value in enumerate(folds)}}


def test_ap_gain_with_probability_regression_rejected():
    assert not accepts(metric(.25, .15, (.25, .25, .25)), metric())
    assert accepts(metric(.25, .1, (.25, .25, .25)), metric())


def test_fold_instability_and_incumbent_regression_rejected():
    assert not accepts(metric(.25, .1, (.4, .18, .18)), metric())
    assert not accepts(metric(.25, .1, (.4, .21, .18)), metric())
    assert not accepts(metric(.25, .1, (.25, .25, .25)), metric(), metric(.3))
    assert not accepts(metric(.25, .1, (.25, None, .25)), metric())


def test_external_whole_groups_and_old_final_excluded():
    rows = [{"candidate_id": str(i), "split": split, "group_id": group}
            for i, (split, group) in enumerate([
                ("development", "g1"), ("validation", "g1"),
                ("development", "g2"), ("validation", "g3"), ("final", "g4")])]
    kept = retained_candidates(rows, {"training_exclusion": {"external_aligned_groups": ["g1"]}})
    assert [r["candidate_id"] for r in kept] == ["2", "3"]


def test_group_folds_are_disjoint_complete_and_deterministic():
    rows = [{"candidate_id": f"{g}-{i}", "group_id": f"g{g}"}
            for g in range(9) for i in range(2)]
    folds = grouped_folds(rows, "test")
    assert folds == grouped_folds(list(reversed(rows)), "test")
    seen = []
    for fold in folds:
        train_groups = {r["group_id"] for r in fold["training"]}
        pool_groups = {r["group_id"] for r in fold["pool"]}
        assert not train_groups & pool_groups
        seen.extend(r["candidate_id"] for r in fold["pool"])
    assert sorted(seen) == sorted(r["candidate_id"] for r in rows)
    with pytest.raises(ValueError):
        grouped_folds(rows[:2], "too-small")


def test_metrics_do_not_treat_unknown_as_negative():
    rows = [{"drug": drug, "label": label, "probability": p, "fold": f"inner-{i}"}
            for drug in ("pem", "fut") for i in range(3)
            for label, p in ((True, .8), (False, .2), (None, .9))]
    result = summarize(rows)
    assert result["macro_ap"] == 1
    assert result["macro_brier"] == pytest.approx(.04)
    assert result["per_drug"]["pem"]["n"] == 6
    assert result["per_drug"]["pem"]["inconclusive_count"] == 3


def test_search_rejects_outer_labels_before_evaluation(tmp_path):
    from escape_lab.revision_evaluation import search_candidates
    rows = [{"candidate_id": str(i), "group_id": f"g{i}"} for i in range(6)]
    labels = {str(i): {"pem": False, "fut": False} for i in range(7)}
    with pytest.raises(ValueError, match="outside"):
        search_candidates(rows, labels, None, None, 4, 0, tmp_path)


def test_nested_search_gets_fresh_planner_and_no_outer_labels(monkeypatch, tmp_path):
    import escape_lab.revision_evaluation as evaluation

    rows = [{"candidate_id": str(i), "group_id": f"g{i}", "gene": "FGFR1"}
            for i in range(9)]
    labels = {str(i): {"pem": bool(i % 2), "fut": bool(i % 2)} for i in range(9)}
    planners, searches = [], []

    def factory():
        planner = object()
        planners.append(planner)
        return planner

    def search(training, inner_labels, store, planner, seed, deadline, output_dir):
        assert set(inner_labels) == {c["candidate_id"] for c in training}
        assert len(inner_labels) < len(labels)
        searches.append(planner)
        return {"selections": {m: "fixed" for m in evaluation.METHODS},
                "agent_status": "completed"}

    def predict(selections, training, inner_labels, pool, store, seed, deadline=None):
        assert not set(inner_labels) & {c["candidate_id"] for c in pool}
        return {"predictions": [
            {"method": method, "candidate_id": c["candidate_id"], "drug": drug,
             "probability": .5}
            for method in selections for c in pool for drug in ("pem", "fut")]}

    monkeypatch.setattr(evaluation, "search_candidates", search)
    monkeypatch.setattr(evaluation, "predict_selected", predict)
    result = evaluation.nested_compare(rows, labels, None, factory, tmp_path, float("inf"))
    assert len(set(planners)) == 3
    assert planners == searches
    assert len(result["predictions"]) == 9 * 2 * 5


def test_expired_budget_stops_before_any_numerical_evaluation(tmp_path):
    from escape_lab.revision_evaluation import search_candidates
    rows = [{"candidate_id": str(i), "group_id": f"g{i}"} for i in range(6)]
    labels = {str(i): {"pem": False, "fut": False} for i in range(6)}
    with pytest.raises(TimeoutError):
        search_candidates(rows, labels, None, None, 4, 0, tmp_path)
