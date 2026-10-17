"""Synthetic regressions for review findings; no real data or inference access."""

import json
import time
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest

from escape_lab import revision_cli, revision_evaluation
from escape_lab.experiment_planner import ExperimentProposal
from escape_lab.revision_models import catalog


def metrics(ap=0.2):
    return {
        "macro_ap": ap,
        "per_drug": {
            drug: {"average_precision": ap, "brier": 0.1, "n": 10, "positive_count": 2}
            for drug in ("pem", "fut")
        },
        "by_fold": {f"inner-{i}": {"macro_ap": ap} for i in range(3)},
    }


class SyntheticPlanner:
    def __init__(self, successful):
        self.successful = successful
        self.calls = 0
        self.attempts = []

    def propose(self, configurations, feedback):
        choices = [key for key in configurations if key != revision_evaluation.REFERENCE]
        candidate = choices[self.calls] if self.calls < self.successful else None
        self.calls += 1
        return ExperimentProposal(
            candidate_id=candidate,
            parent_id=revision_evaluation.REFERENCE if candidate else None,
            hypothesis="Synthetic diagnostic hypothesis." if candidate else "",
            expected_effect="AP improvement." if candidate else "",
            rejection_condition="No improvement." if candidate else "",
            status="proposed" if candidate else "failed",
            error=None if candidate else "Synthetic failure",
            calls=1,
            model_id="synthetic-no-network",
        )


@pytest.mark.parametrize(
    ("successes", "status"),
    [(0, "unavailable"), (1, "degraded"), (2, "degraded"), (3, "completed")],
)
def test_failed_opportunities_are_reported_separately_from_baseline_selection(
    monkeypatch, tmp_path, successes, status
):
    rows = [{"candidate_id": str(i), "group_id": f"g{i}", "gene": "FGFR1"} for i in range(9)]
    labels = {row["candidate_id"]: {"pem": True, "fut": False} for row in rows}

    def evaluate(config_id, tasks, supplied_labels, store, seed, deadline):
        assert set(supplied_labels) == set(labels)
        ap = 0.2 if config_id == revision_evaluation.REFERENCE else 0.25
        return {"candidate_id": config_id, "metrics": metrics(ap), "seconds": 0.01}

    monkeypatch.setattr(revision_evaluation, "evaluate_cv", evaluate)
    result = revision_evaluation.search_candidates(
        rows, labels, None, SyntheticPlanner(successes), 1, time.monotonic() + 60, tmp_path
    )
    assert result["agent_status"] == status
    assert result["successful_proposals"] == successes
    assert result["failed_opportunities"] == 3 - successes
    assert len(result["agent_proposals"]) == 3
    if successes == 0:
        assert result["selections"]["agent"] == revision_evaluation.REFERENCE
        assert result["agent_status"] != "completed"


def test_expired_search_does_not_start_numerical_fits(monkeypatch, tmp_path):
    rows = [{"candidate_id": str(i), "group_id": f"g{i}"} for i in range(9)]
    labels = {row["candidate_id"]: {"pem": True, "fut": False} for row in rows}

    def should_not_evaluate(*args, **kwargs):
        pytest.fail("Numerical evaluation must not begin after expiration")

    monkeypatch.setattr(revision_evaluation, "evaluate_cv", should_not_evaluate)
    with pytest.raises(TimeoutError):
        revision_evaluation.search_candidates(
            rows, labels, None, SyntheticPlanner(3), 1, 0, tmp_path
        )


@pytest.mark.parametrize("undefined", [None, float("nan"), float("inf")])
def test_undefined_fold_metric_rejects_without_crashing(undefined):
    candidate = metrics(0.25)
    candidate["by_fold"]["inner-0"]["macro_ap"] = undefined
    assert revision_evaluation.acceptance_reasons(candidate, metrics()) == ["undefined_metric"]


@pytest.mark.parametrize(
    ("outer_status", "full_status", "permitted"),
    [
        ("completed", "completed", True),
        ("degraded", "degraded", True),
        ("unavailable", "completed", False),
        ("completed", "unavailable", False),
    ],
)
def test_external_gate_requires_a_success_in_every_search_and_commits_before_reveal(
    monkeypatch, tmp_path, outer_status, full_status, permitted
):
    root = tmp_path / "project"
    output = tmp_path / "run"
    panel_dir = root / "data/external_discovery"
    processed = root / "data/processed"
    panel_dir.mkdir(parents=True)
    processed.mkdir(parents=True)
    evaluator_source = root / "src/escape_lab/external_panel.py"
    evaluator_source.parent.mkdir(parents=True)
    evaluator_source.write_text("# Synthetic evaluator fixture; never executed.\n")
    rows = [
        {"candidate_id": str(i), "group_id": f"g{i % 15}", "split": "development"}
        for i in range(3730)
    ]
    (processed / "candidates.json").write_text(json.dumps(rows))
    panel = {
        "training_exclusion": {"external_aligned_groups": ["external-group"]},
        "document_sha256": "synthetic-document-hash",
        "variants": [{"variant_id": "FGFR1:A1V", "group_id": "external-group"}],
    }
    (panel_dir / "s13-panel.json").write_text(json.dumps(panel))
    calls = []
    selections = {method: revision_evaluation.REFERENCE for method in revision_evaluation.METHODS}
    nested_metrics = {method: metrics() for method in revision_evaluation.METHODS}

    def inspect_panel(path):
        calls.append("inspect_metadata")
        return {}

    def feature_store(path):
        return object()

    def freeze_protocol(project_root, output_dir, candidates, panel_path):
        assert len(candidates) == 3730
        protocol = {"hashes": {}}
        revision_cli.write_json(output_dir / "protocol.json", protocol)
        calls.append("freeze")
        return protocol

    def load_labels(data_dir, candidates):
        assert "freeze" in calls
        return {row["candidate_id"]: {"pem": True, "fut": False} for row in candidates}

    def nested_compare(*args):
        calls.append("nested")
        return {
            "agent_status": {
                "outer-0": outer_status,
                "outer-1": "completed",
                "outer-2": "completed",
            },
            "metrics": deepcopy(nested_metrics),
        }

    def search(*args):
        calls.append("full_search")
        return {"agent_status": full_status, "selections": selections}

    def predict(selected, training, labels, pool, store, seed, model_dir, *, deadline):
        calls.append("predict")
        assert set(selected) == set(revision_evaluation.METHODS)
        assert len(labels) == 3730
        assert not {row["candidate_id"] for row in pool} & set(labels)
        return {
            "predictions": [
                {"method": method, "candidate_id": "FGFR1:A1V", "drug": drug, "probability": 0.5}
                for method in revision_evaluation.METHODS
                for drug in ("pem", "fut")
            ],
            "model_hashes": {},
            "diagnostics": [],
        }

    def verify_hashes(hashes):
        for filename, expected in hashes.items():
            assert Path(filename).is_file()
            assert revision_cli.sha256(filename) == expected

    def reveal(panel_path, freeze_path, prediction_path, *, receipt_path, output_dir):
        calls.append("reveal")
        assert calls.index("predict") < calls.index("reveal")
        freeze = json.loads(Path(freeze_path).read_text())
        committed = json.loads(Path(prediction_path).read_text())
        assert freeze["agent_statuses"] == [outer_status, "completed", "completed", full_status]
        assert len(committed["predictions"]) == 10
        assert (Path(prediction_path).parent / "fitted-models.json").is_file()
        assert freeze["external_evaluation"]["predictions_sha256"] == revision_cli.sha256(
            prediction_path
        )
        receipt_path.parent.mkdir(parents=True, exist_ok=True)
        receipt_path.write_text(json.dumps({"status": "completed", "synthetic": True}))
        return {"synthetic": True}

    monkeypatch.setattr(revision_cli, "inspect_external_panel", inspect_panel)
    monkeypatch.setattr(revision_cli, "BiologicalFeatureStore", feature_store)
    monkeypatch.setattr(revision_cli, "freeze_protocol", freeze_protocol)
    monkeypatch.setattr(revision_cli, "load_labels", load_labels)
    monkeypatch.setattr(revision_cli, "nested_compare", nested_compare)
    monkeypatch.setattr(revision_cli, "search_candidates", search)
    monkeypatch.setattr(revision_cli, "predict_selected", predict)
    monkeypatch.setattr(revision_cli, "verify_hashes", verify_hashes)
    monkeypatch.setattr(revision_cli, "evaluate_external_panel", reveal)
    client = SimpleNamespace(total_tokens=31, usage_records=[])
    if permitted:
        result = revision_cli.run(root, client, output)
        assert "reveal" in calls
        assert result["agent_statuses"] == [outer_status, "completed", "completed", full_status]
    else:
        with pytest.raises(RuntimeError, match="no evaluated proposal"):
            revision_cli.run(root, client, output)
        assert "predict" not in calls and "reveal" not in calls
        assert not (root / "data/revision-v2/external-receipt.json").exists()
    status = json.loads((output / "run-status.json").read_text())
    assert status["status"] == ("completed" if permitted else "failed")


def test_planner_catalog_contract_matches_all_twelve_numerical_configs():
    from escape_lab.experiment_planner import _safe_catalog

    configurations = {key: value.model_dump() for key, value in catalog().items()}
    assert _safe_catalog(configurations) == configurations
