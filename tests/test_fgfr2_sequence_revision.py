import copy
import json

import pytest

from escape_lab.fgfr2_sequence_revision import (
    ARM_CONFIGS,
    approved_cohort,
    public_result,
    stability_gate,
)


def candidates():
    return [{"candidate_id": f"variant-{i}", "gene": "FGFR2", "group_id": f"g{i}",
             "split": "development"} for i in range(3)]


def protocol():
    return {"eligible_ids": [c["candidate_id"] for c in candidates()],
            "fold_groups": {f"fold-{i}": [f"g{i}"] for i in range(3)},
            "reserved_groups": ["reserved"]}


def test_cohort_is_exactly_prior_allowlist_and_does_not_read_other_outcomes():
    rows = candidates() + [{"candidate_id": "sealed", "split": "final", "outcome": object()}]
    assert approved_cohort(rows, protocol(), expected_n=3) == candidates()
    rows[0]["split"] = "final"
    with pytest.raises(ValueError, match="development"):
        approved_cohort(rows, protocol(), expected_n=3)


def test_missing_duplicate_and_reserved_candidates_fail_closed():
    with pytest.raises(ValueError, match="cohort"):
        approved_cohort(candidates()[:-1], protocol(), expected_n=3)
    with pytest.raises(ValueError, match="cohort"):
        approved_cohort(candidates() + [candidates()[0]], protocol(), expected_n=3)
    rows = candidates()
    rows[0]["group_id"] = "reserved"
    with pytest.raises(ValueError, match="reserved"):
        approved_cohort(rows, protocol(), expected_n=3)


def test_predeclared_design_changes_one_factor_and_has_twelve_fits():
    assert ARM_CONFIGS == {
        "legacy_both": {"stability": "legacy", "sequence_branches": "both"},
        "floor_both": {"stability": "variance_floor", "sequence_branches": "both"},
        "legacy_wt": {"stability": "legacy", "sequence_branches": "wt"},
        "legacy_delta": {"stability": "legacy", "sequence_branches": "delta"},
    }
    assert len(ARM_CONFIGS) * 3 == 12


def metric(ap, hits=4):
    return {"mean_fold_high_ap": ap, "high_hits": hits,
            "by_fold": {"fold": {"high": {"ap": ap}}}}


def test_stability_gate_needs_ranking_gain_and_stability_and_no_retrieval_loss():
    before = metric(.11)
    assert stability_gate(before, metric(.13), 100000., 10.)["status"] == "pass"
    for after, maximum in ((metric(.11), 10.), (metric(.13, hits=3), 10.),
                           (metric(.13), 50000.)):
        assert stability_gate(before, after, 100000., maximum)["status"] == "fail"


def test_public_result_uses_allowlist_and_drops_private_records_and_groups():
    raw = {"version": "test", "metrics": {}, "gates": {}, "fit_count": 12,
           "protocol_sha256": "abc", "baseline_reproduction": {}, "integrity": {},
           "fitting_seconds": 1., "elapsed_seconds": 2., "failures": [],
           "private_predictions": [{"truth": 99}], "eligible_ids": ["secret"],
           "diagnostics": [{"arm": "legacy_both", "fold": "f", "margin_abs_max": 3.,
                            "model": {"training_groups": ["secret"], "coefficient_norm": 1.}}]}
    original = copy.deepcopy(raw)
    result = public_result(raw)
    assert "secret" not in json.dumps(result)
    assert "truth" not in json.dumps(result)
    assert raw == original
