import json

import numpy as np
import pytest
from sklearn.metrics import average_precision_score

from escape_lab.data import sha256
from escape_lab.revision_uncertainty import (
    METHODS,
    _metric,
    _paired_arrays,
    paired_group_bootstrap,
    run_uncertainty,
)


def records(endpoint, same_scores=False):
    rows = []
    for method in METHODS:
        for i in range(12):
            for drug in ("pem", "fut"):
                score = (i + 1) / 13
                if method != "agent" and not same_scores:
                    score = 1 - score
                row = {"method": method, "drug": drug, "group_id": f"g{i // 2}"}
                if endpoint == "internal":
                    row.update(candidate_id=f"v{i}", label=i >= 6, probability=score)
                else:
                    row.update(variant_id=f"v{i}", residual_activity=float(i),
                               measurement_status="observed", prediction_score=score)
                rows.append(row)
    return rows


@pytest.mark.parametrize("endpoint", ["internal", "external"])
def test_identical_methods_have_exact_zero_paired_intervals(endpoint):
    result = paired_group_bootstrap(records(endpoint, same_scores=True), endpoint,
                                    replicates=20, seed=1)
    assert result["aligned_group_count"] == 6
    assert set(result["comparisons"]) == {
        "biology_fixed_minus_v1_fixed", "random_minus_v1_fixed", "grid_minus_v1_fixed",
        "agent_minus_v1_fixed", "agent_minus_biology_fixed", "agent_minus_random",
        "agent_minus_grid",
    }
    for comparison in result["comparisons"].values():
        assert comparison["delta"] == 0
        assert comparison["percentile_2_5"] == 0
        assert comparison["percentile_97_5"] == 0
        assert comparison["valid_replicates"] == 20


def test_external_direction_determinism_and_pairing():
    rows = records("external")
    first = paired_group_bootstrap(rows, "external", replicates=20, seed=7)
    second = paired_group_bootstrap(list(reversed(rows)), "external", replicates=20, seed=7)
    assert first == second
    assert first["point_estimates"]["agent"] == pytest.approx(1)
    assert first["comparisons"]["agent_minus_v1_fixed"]["delta"] == pytest.approx(2)


def test_biology_baseline_contrast_is_independent_of_agent_scores():
    rows = records("external")
    for row in rows:
        if row["method"] == "biology_fixed":
            row["prediction_score"] = 1 - row["prediction_score"]
        elif row["method"] == "agent":
            row["prediction_score"] = .5
    result = paired_group_bootstrap(rows, "external", replicates=20, seed=7)
    biology = result["comparisons"]["biology_fixed_minus_v1_fixed"]
    assert biology["method"] == "biology_fixed"
    assert biology["comparator_method"] == "v1_fixed"
    assert biology["delta"] == pytest.approx(2)
    assert biology["valid_replicates"] == 20
    assert biology["percentile_2_5"] == pytest.approx(2)
    assert result["comparisons"]["agent_minus_v1_fixed"]["valid_replicates"] == 0


def test_block_weights_equal_repeating_whole_blocks():
    rows = records("internal")
    groups, arrays = _paired_arrays(rows, "internal")
    counts = np.array([0, 2, 1, 0, 1, 2])
    assert len(groups) == len(counts)
    actual = _metric(arrays, "v1_fixed", counts, "internal")
    truth, scores = [], []
    for row in rows:
        if row["method"] == "v1_fixed" and row["drug"] == "pem":
            for _ in range(counts[int(row["group_id"][1:])]):
                truth.append(row["label"])
                scores.append(row["probability"])
    assert actual == pytest.approx(average_precision_score(truth, scores))


def test_unknown_and_censored_are_excluded_without_changing_pairing():
    internal, external = records("internal"), records("external")
    for row in internal:
        if row["candidate_id"] == "v0":
            row["label"] = None
    for row in external:
        if row["variant_id"] == "v0":
            row.update(residual_activity=None, measurement_status="censored")
    assert paired_group_bootstrap(internal, "internal", 5)["observed_pairs_per_drug"] == {
        "pem": 11, "fut": 11}
    assert paired_group_bootstrap(external, "external", 5)["observed_pairs_per_drug"] == {
        "pem": 11, "fut": 11}


def test_constant_external_draws_are_counted_as_invalid():
    rows = records("external")
    for row in rows:
        row["prediction_score"] = .5
    result = paired_group_bootstrap(rows, "external", 20)
    for comparison in result["comparisons"].values():
        assert comparison["delta"] is None
        assert comparison["valid_replicates"] == 0
        assert comparison["invalid_replicates"] == 20
        assert comparison["percentile_2_5"] is None


@pytest.mark.parametrize("change", ["missing", "duplicate", "label", "group", "nan"])
def test_unpaired_or_invalid_inputs_rejected(change):
    rows = records("internal")
    if change == "missing":
        rows.pop()
    elif change == "duplicate":
        rows.append(dict(rows[0]))
    elif change == "label":
        rows[0]["label"] = True
    elif change == "group":
        rows[0]["group_id"] = "different"
    else:
        rows[0]["probability"] = float("nan")
    with pytest.raises(ValueError):
        paired_group_bootstrap(rows, "internal", 5)


def test_incomplete_run_rejects_before_reading_outcome_artifacts(tmp_path):
    (tmp_path / "run-status.json").write_text(json.dumps({"status": "started"}))
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps({"status": "completed"}))
    with pytest.raises(ValueError, match="must be completed"):
        run_uncertainty(tmp_path, receipt)


def test_released_record_hash_tampering_rejected(tmp_path):
    external = tmp_path / "external"
    external.mkdir()
    record_path = external / "external-records.json"
    record_path.write_text(json.dumps({"records": records("external")}))
    digest = sha256(record_path)
    report_path = external / "external-report.json"
    report_path.write_text(json.dumps({"records_path": str(record_path), "records_sha256": digest}))
    results = tmp_path / "results.json"
    results.write_text("{}")
    (tmp_path / "run-status.json").write_text(json.dumps({
        "status": "completed", "results_sha256": sha256(results)}))
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps({"status": "completed", "report_path": str(report_path),
                                  "report_sha256": sha256(report_path), "records_sha256": digest}))
    record_path.write_text("{}")
    with pytest.raises(ValueError, match="checksum mismatch"):
        run_uncertainty(tmp_path, receipt)
