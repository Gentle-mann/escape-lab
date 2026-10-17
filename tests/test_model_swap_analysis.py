import copy
import json

import pytest

from escape_lab import model_swap_analysis
from escape_lab.data import sha256
from escape_lab.improvement import write_json
from escape_lab.model_swap_analysis import (
    BASELINES,
    METHODS,
    assert_equivalent,
    compare_predictions,
    estimate_cost,
    load_completed_run,
    paired_model_swap_bootstrap,
    render_summary,
    run_comparison,
)
from escape_lab.revision_evaluation import summarize
from escape_lab.revision_uncertainty import paired_group_bootstrap


def records(endpoint, agent_direction=1):
    rows = []
    for method in (*BASELINES, "agent"):
        for i in range(12):
            for drug in ("pem", "fut"):
                score = (i + 1) / 13
                if method == "agent" and agent_direction == -1:
                    score = 1 - score
                row = {"method": method, "drug": drug, "group_id": f"g{i // 2}"}
                if endpoint == "internal":
                    row.update(candidate_id=f"v{i}", label=i >= 6, probability=score,
                               fold=f"outer-{i // 4}")
                else:
                    row.update(variant_id=f"v{i}", residual_activity=float(i),
                               measurement_status="observed", prediction_score=score)
                rows.append(row)
    return rows


@pytest.mark.parametrize("endpoint", ["internal", "external"])
def test_same_draws_and_old_contrasts_match_frozen_framework(endpoint):
    reference = records(endpoint, -1)
    glm = records(endpoint)
    expected = paired_group_bootstrap(reference, endpoint, replicates=30, seed=6100)
    actual = paired_model_swap_bootstrap(reference, glm, endpoint, replicates=30, seed=6100)
    assert set(actual["point_estimates"]) == set(METHODS)
    for key, value in expected["comparisons"].items():
        translated = key.replace("agent", "nemotron")
        renamed = {k: v.replace("agent", "nemotron") if isinstance(v, str) else v
                   for k, v in value.items()}
        assert actual["comparisons"][translated] == renamed
    delta = actual["comparisons"]["glm_minus_nemotron"]
    assert delta["delta"] > 0
    if endpoint == "external":
        assert delta["delta"] == pytest.approx(2)
        assert delta["percentile_2_5"] == pytest.approx(2)
    assert actual == paired_model_swap_bootstrap(
        list(reversed(reference)), list(reversed(glm)), endpoint, 30, 6100)


@pytest.mark.parametrize("endpoint", ["internal", "external"])
def test_identical_agent_predictions_have_exact_zero_intervals(endpoint):
    actual = paired_model_swap_bootstrap(records(endpoint), records(endpoint), endpoint, 20)
    contrast = actual["comparisons"]["glm_minus_nemotron"]
    assert contrast["delta"] == 0
    assert contrast["percentile_2_5"] == contrast["percentile_97_5"] == 0


@pytest.mark.parametrize("field", ["probability", "fold", "label", "group_id"])
def test_changed_baseline_predictions_rejected(field):
    reference = records("internal")
    glm = copy.deepcopy(reference)
    glm[0][field] = {"probability": .9, "fold": "outer-1", "label": True,
                     "group_id": "wrong"}[field]
    with pytest.raises(ValueError, match="baseline|paired"):
        compare_predictions(reference, glm, "internal")


def test_cross_run_agent_outcome_mismatch_rejected():
    reference = records("external")
    glm = copy.deepcopy(reference)
    for row in glm:
        if row["variant_id"] == "v0":
            row["residual_activity"] = 99
    with pytest.raises(ValueError):
        paired_model_swap_bootstrap(reference, glm, "external", 2)


def test_metric_comparison_uses_strict_absolute_tolerance():
    assert_equivalent({"x": .5}, {"x": .5 + 1e-13}, "metrics")
    with pytest.raises(ValueError, match="metrics"):
        assert_equivalent({"x": .5}, {"x": .5 + 1e-9}, "metrics")


def test_incomplete_run_rejected_before_other_reads(tmp_path):
    (tmp_path / "run-status.json").write_text('{"status":"started"}')
    with pytest.raises(ValueError, match="completed"):
        load_completed_run(tmp_path)


def test_tampered_result_rejected_before_outcome_reads(tmp_path):
    result = tmp_path / "results.json"
    result.write_text("{}")
    digest = sha256(result)
    (tmp_path / "run-status.json").write_text(json.dumps({
        "status": "completed", "results_sha256": digest}))
    result.write_text('{"changed":true}')
    with pytest.raises(ValueError, match="checksum"):
        load_completed_run(tmp_path)


def test_uncached_cost_basis_is_explicit():
    accounting = {"usage_records": [
        {"prompt_tokens": 95604, "completion_tokens": 53354, "usage_estimated": False}]}
    result = estimate_cost(accounting, "nemotron")
    assert result["usd"] == pytest.approx(.030902)
    assert result["invoice_verified"] is False
    assert result["cached_discount_applied"] is False
    assert result["source_url"] == "https://www.crusoe.ai/cloud/pricing"
    with pytest.raises(ValueError):
        estimate_cost(accounting, "unrecognized")


def test_summary_identifies_reused_measurements_and_primary_metric():
    # Empty endpoint tables are sufficient to test mandatory interpretation text.
    payload = {"reference_run": "/ref", "glm_run": "/glm", "methods": {},
               "runs": {}, "uncertainty": {}, "audit": {"limitations": []},
               "reproduce_command": "uv run python -m escape_lab.model_swap_analysis /ref /glm"}
    result = render_summary(payload)
    assert "reused" in result.lower()
    assert "single" in result.lower()
    assert "pooled macro AP" in result
    assert "fresh blind" in result


def completed_fixture(path, model, agent_direction):
    path.mkdir()
    internal, external = records("internal", agent_direction), records("external", agent_direction)
    metrics = {method: summarize([row for row in internal if row["method"] == method])
               for method in (*BASELINES, "agent")}
    methods = {}
    for method in (*BASELINES, "agent"):
        score = -1.0 if method == "agent" and agent_direction == -1 else 1.0
        methods[method] = {"macro_spearman": score, "per_drug": {
            drug: {"spearman": score, "n": 12, "aligned_groups": 6, "reason": None}
            for drug in ("pem", "fut")}}
    selections = {method: "baseline" for method in (*BASELINES, "agent")}
    selections["agent"] = f"selected-{model}"
    search = {"selections": selections, "random_candidates": ["a", "b", "c"],
              "inner_groups": {"inner-0": ["g0", "g1"]}, "inner_metrics": {"baseline": {}},
              "agent_proposals": [{"status": "proposed", "metrics": {}, "accepted": False}] * 3,
              "agent_attempts": [{"attempt": 1}] * 3,
              "successful_proposals": 3, "failed_opportunities": 0}
    for fold in range(3):
        write_json(path / f"nested/outer-{fold}/search.json", search)
    write_json(path / "full-search/search.json", search)
    write_json(path / "nested/nested-results.json", {
        "metrics": metrics, "predictions": internal,
        "selections": {f"outer-{fold}": selections for fold in range(3)},
        "outer_groups": {f"outer-{fold}": [f"g{fold * 2}", f"g{fold * 2 + 1}"]
                         for fold in range(3)}})
    write_json(path / "external/external-records.json", {"records": external})
    report = {"status": "completed", "methods": methods,
              "records_path": str(path / "external/external-records.json"),
              "records_sha256": sha256(path / "external/external-records.json")}
    write_json(path / "external/external-report.json", report)
    write_json(path / "external-predictions.json", {"predictions": [
        {"method": r["method"], "drug": r["drug"], "variant_id": r["variant_id"],
         "score": r["prediction_score"]} for r in external]})
    model_id = "zai-org/GLM-5.3" if model == "glm" else "nvidia/NVIDIA-Nemotron-3.5-Super-VL"
    write_json(path / "protocol.json", {"model": model_id, "created_at": model,
                                       "proposal_opportunities_per_search": 3})
    commitment = {str(path / name): sha256(path / name)
                  for name in ("protocol.json", "external-predictions.json")}
    write_json(path / "external-freeze.json", {"retrospective_reuse_authorized": True,
               "selected_configurations": selections, "hashes": commitment})
    write_json(path / "results.json", {
        "nested_metrics": metrics, "external": report, "full_selections": selections,
        "agent_statuses": ["completed"] * 4, "seconds": 2,
        "token_accounting": {"budget_charge_tokens": 36,
                             "usage_records": [{"prompt_tokens": 2, "completion_tokens": 1,
                                                "total_tokens": 3}] * 12}})
    reseal_fixture(path)


def reseal_fixture(path):
    write_json(path / "run-status.json", {
        "status": "completed", "results_sha256": sha256(path / "results.json"),
        "protocol_sha256": sha256(path / "protocol.json"),
        "artifact_sha256": {str(p): sha256(p) for p in path.rglob("*.json")
                            if p.name != "run-status.json"}})


def test_completed_analysis_end_to_end_and_no_input_changes(tmp_path, monkeypatch):
    reference, glm = tmp_path / "reference", tmp_path / "glm"
    completed_fixture(reference, "nemotron", -1)
    completed_fixture(glm, "glm", 1)
    original = {p: sha256(p) for p in tmp_path.rglob("*.json")}
    bootstrap = paired_model_swap_bootstrap

    def short_bootstrap(reference_records, glm_records, endpoint):
        return bootstrap(reference_records, glm_records, endpoint, replicates=20)

    monkeypatch.setattr(model_swap_analysis, "paired_model_swap_bootstrap", short_bootstrap)
    result = run_comparison(reference, glm)
    assert result["uncertainty"]["external"]["comparisons"]["glm_minus_nemotron"]["delta"] == 2
    assert result["runs"]["glm"]["valid_evaluated_proposals"] == 12
    assert (glm / "comparison.json").exists()
    assert (glm / "results-summary.md").exists()
    assert all(sha256(p) == digest for p, digest in original.items())
    with pytest.raises(FileExistsError):
        run_comparison(reference, glm)


def test_completed_status_cannot_hide_prediction_tampering(tmp_path):
    completed_fixture(tmp_path / "run", "glm", 1)
    path = tmp_path / "run/external/external-records.json"
    value = json.loads(path.read_text())
    value["records"][0]["prediction_score"] = .999
    write_json(path, value)
    with pytest.raises(ValueError, match="checksum"):
        load_completed_run(tmp_path / "run")


def test_unbound_nested_predictions_are_rejected(tmp_path):
    run = tmp_path / "run"
    completed_fixture(run, "glm", 1)
    status = json.loads((run / "run-status.json").read_text())
    del status["artifact_sha256"][str(run / "nested/nested-results.json")]
    write_json(run / "run-status.json", status)
    with pytest.raises(ValueError, match="Missing completed artifact checksum"):
        load_completed_run(run)
