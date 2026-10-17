import copy

import pytest

from escape_lab.fgfr2_sequence_run import (
    advancement_gate,
    fold_partitions,
    paired_uncertainty,
    ranking_metrics,
    validate_reference,
)


def candidate(cid, group, split="development"):
    return {"candidate_id": cid, "group_id": group, "gene": "FGFR2", "split": split,
            "position": int(cid), "aligned_position": int(cid), "wt": "A", "mut": "V"}


def test_saved_global_fold_map_keeps_empty_fgfr2_group():
    rows = [candidate("1", "a"), candidate("2", "c"), candidate("3", "d")]
    mapping = {"pilot-0": ["a", "b"], "pilot-1": ["c"], "pilot-2": ["d"]}
    folds = fold_partitions(rows, mapping)
    assert [r["candidate_id"] for r in folds[0][2]] == ["1"]
    assert [r["candidate_id"] for r in folds[1][2]] == ["2"]
    assert {folds[0][1][i]["group_id"] for i in range(2)} == {"c", "d"}


@pytest.mark.parametrize("bad", ["final", "FGFR1", "missing", "duplicate"])
def test_cohort_boundaries(bad):
    rows = [candidate("1", "a"), candidate("2", "c"), candidate("3", "d")]
    mapping = {"pilot-0": ["a"], "pilot-1": ["c"], "pilot-2": ["d"]}
    if bad == "final":
        rows[0]["split"] = "final"
    elif bad == "FGFR1":
        rows[0]["gene"] = "FGFR1"
    elif bad == "missing":
        rows[0]["group_id"] = "unknown"
    else:
        mapping["pilot-2"].append("a")
    with pytest.raises(ValueError):
        fold_partitions(rows, mapping)


def test_ranking_ties_and_undefined_positive_class():
    rows = [{"candidate_id": str(i).zfill(2), "score": 0.5, "high": i == 0,
             "strong": False, "truth": i} for i in range(11)]
    result = ranking_metrics(rows)
    assert result["budget"] == 2
    assert result["high"]["hits"] == 1
    assert result["high"]["ap"] == pytest.approx(1 / 11)
    assert result["strong"]["ap"] is None
    assert result["strong"]["recall"] is None
    assert ranking_metrics(list(reversed(rows))) == result


def test_reference_identity_and_truth_are_checked():
    c = candidate("1", "a")
    labels = {"1": {"value": 2.0, "active": True}}
    saved = [{**c, "fold": "pilot-0", "truth": 2.0, "active": True, "prediction": 0.7}]
    assert validate_reference([c], labels, saved, {"a": "pilot-0"})["1"] == 0.7
    for field, value in (("truth", 3.0), ("group_id", "x"), ("fold", "pilot-1")):
        bad = copy.deepcopy(saved)
        bad[0][field] = value
        with pytest.raises(ValueError):
            validate_reference([c], labels, bad, {"a": "pilot-0"})


def summary(ap, hits, folds):
    return {"mean_fold_high_ap": ap, "high_hits": hits,
            "by_fold": {str(i): {"high": {"ap": v}} for i, v in enumerate(folds)}}


def test_gate_requires_all_declared_conditions():
    metrics = {"cached_v2_biology": summary(.1, 2, [.1, .1, .1]),
               "biology_logistic": summary(.12, 4, [.12, .12, .12]),
               "esm_biology_logistic": summary(.2, 5, [.2, .2, .2])}
    assert advancement_gate(metrics)["status"] == "pass"
    metrics["esm_biology_logistic"]["high_hits"] = 4
    assert advancement_gate(metrics)["status"] == "fail"
    metrics["esm_biology_logistic"]["mean_fold_high_ap"] = None
    assert advancement_gate(metrics)["status"] == "inconclusive"


def test_exact_ap_threshold_is_inclusive():
    metrics = {"cached_v2_biology": summary(.1, 2, [.1, .1, .1]),
               "biology_logistic": summary(.1, 4, [.1, .1, .1]),
               "esm_biology_logistic": summary(.13, 5, [.13, .13, .13])}
    assert advancement_gate(metrics)["status"] == "pass"


def test_bootstrap_rejects_unpaired_targets():
    r = {"candidate_id": "1", "fold": "pilot-0", "group_id": "a", "high": True,
         "strong": False, "truth": 3.0, "score": .4}
    rows = {arm: [copy.deepcopy(r)] for arm in
            ("cached_v2_biology", "biology_logistic", "esm_biology_logistic")}
    rows["esm_biology_logistic"][0]["high"] = False
    with pytest.raises(ValueError, match="target or partition"):
        paired_uncertainty(rows, draws=1)
