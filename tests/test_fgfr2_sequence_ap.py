import numpy as np
import pytest
from scipy.special import expit

from escape_lab.fgfr2_sequence_ap import paired_ap_interval, ranking_outputs, screening_gate


class Head:
    def decision_function(self, matrix):
        return matrix[:, 0]


class Model:
    arm = "esm_biology_logistic"
    model = Head()

    def predict(self, candidates):
        if candidates[0].get("split") == "final":
            raise ValueError("final forbidden")
        return expit([c["margin"] for c in candidates])

    def _predict_matrix(self, candidates, ids):
        return np.array([[c["margin"]] for c in candidates])


def test_raw_scores_resolve_probability_saturation_without_changing_model():
    candidates = [{"candidate_id": str(i), "margin": margin} for i, margin in enumerate((40., 60., 80.))]
    out = ranking_outputs(Model(), candidates)
    np.testing.assert_array_equal(out["probability"], [1., 1., 1.])
    np.testing.assert_array_equal(out["margin"], [40., 60., 80.])


def test_unsaturated_ranking_and_scope_guards_preserved():
    candidates = [{"candidate_id": str(i), "margin": margin} for i, margin in enumerate((-3., 1., 2.))]
    out = ranking_outputs(Model(), candidates)
    np.testing.assert_array_equal(np.argsort(out["probability"]), np.argsort(out["margin"]))
    candidates[0]["split"] = "final"
    with pytest.raises(ValueError, match="final"):
        ranking_outputs(Model(), candidates)


def test_nonfinite_margin_rejected():
    with pytest.raises(ValueError):
        ranking_outputs(Model(), [{"candidate_id": "a", "margin": float("inf")}])


def metric(ap, hits=1, folds=(.03, .21, .05)):
    return {"mean_fold_high_ap": ap, "high_hits": hits,
            "by_fold": {str(i): {"high": {"ap": value}} for i, value in enumerate(folds)}}


def test_gate_requires_ap_gain_and_no_retrieval_or_fold_loss():
    base = metric(.097)
    assert screening_gate(base, metric(.107, folds=(.06, .21, .05)))["status"] == "pass"
    for proposed in (metric(.10), metric(.12, hits=0), metric(.12, folds=(.2, .20, .05))):
        assert screening_gate(base, proposed)["status"] == "fail"


def test_paired_bootstrap_identity_and_target_alignment():
    rows = [{"candidate_id": str(i), "group_id": "g", "fold": "f", "truth": float(i),
             "high": i == 2, "strong": False, "score": float(i)} for i in range(3)]
    arms = {"probability": rows, "margin": [dict(r) for r in reversed(rows)]}
    result = paired_ap_interval(arms, draws=10)
    assert result["valid_draws"] == 10
    assert result["percentile_2_5"] == 0.
    assert result["percentile_97_5"] == 0.
    arms["margin"][0]["group_id"] = "bad"
    with pytest.raises(ValueError):
        paired_ap_interval(arms, draws=10)
