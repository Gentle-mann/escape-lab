"""One matched v2 GLM substitution, with explicitly retrospective S13 scoring.

Original experiment modules, receipts and artifacts are read-only. The unchanged
planner gets only fresh inner-CV feedback. Released S13 outcomes are read by the
evaluator only after predictions and model files have been hash-committed.
"""

import argparse
import hashlib
import json
import time
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from uuid import uuid4

from escape_lab.biological_features import BiologicalFeatureStore
from escape_lab.cli import credentials
from escape_lab.crusoe import CrusoeClient
from escape_lab.data import sha256
from escape_lab.experiment_planner import ExperimentPlanner
from escape_lab.external_panel import _predictions, _summaries
from escape_lab.improvement import load_labels, write_json
from escape_lab.revision_cli import MODEL as REFERENCE_MODEL
from escape_lab.revision_cli import RUN_SECONDS, TOKEN_BUDGET, verify_hashes
from escape_lab.revision_evaluation import (
    METHODS,
    check_deadline,
    nested_compare,
    predict_selected,
    retained_candidates,
    search_candidates,
)
from escape_lab.revision_models import catalog

MODEL = "zai-org/GLM-5.3"
REFERENCE_RUN = "20261017T201321Z-a29b17a8"


class AuditedClient(CrusoeClient):
    """Same wire payload and failure policy; extra safe request/usage evidence."""

    def __init__(self, *args, audit_path, **kwargs):
        super().__init__(*args, **kwargs)
        self.audit_path = Path(audit_path)

    def _request(self, method, path, payload=None):
        if method != "POST":
            return super()._request(method, path, payload)
        start = time.monotonic()
        record = {"started_at": datetime.now(UTC).isoformat(), "request": payload,
                  "request_sha256": hashlib.sha256(
                      json.dumps(payload, sort_keys=True).encode()).hexdigest()}
        try:
            result = super()._request(method, path, payload)
            record["response"] = {key: result.get(key) for key in
                                  ("id", "model", "created", "system_fingerprint", "usage")}
            record["finish_reasons"] = [c.get("finish_reason")
                                        for c in result.get("choices", [])]
            return result
        except BaseException as error:
            record["error_type"] = type(error).__name__
            raise
        finally:
            record["seconds"] = time.monotonic() - start
            with self.audit_path.open("a") as handle:
                handle.write(json.dumps(record, allow_nan=False) + "\n")


def read_json(path):
    return json.loads(Path(path).read_text())


def verify_commitment(freeze_path):
    freeze = read_json(freeze_path)
    if freeze.get("retrospective_reuse_authorized") is not True:
        raise ValueError("Explicit retrospective reuse authorization required")
    if not freeze.get("hashes"):
        raise ValueError("Prediction commitment is empty")
    verify_hashes(freeze["hashes"])
    return freeze


def released_measurements(rows, manifest):
    """Require all original methods to agree on each released measurement."""
    expected = {(m, v["variant_id"], d) for m in METHODS
                for v in manifest["variants"] for d in ("pem", "fut")}
    indexed = {(r["method"], r["variant_id"], r["drug"]): r for r in rows}
    if len(indexed) != len(rows) or set(indexed) != expected:
        raise ValueError("Released measurement coverage differs")
    measurements = {}
    fields = ("residual_activity", "measurement_status", "source_cell", "group_id",
              "bound", "inclusive_bound", "outside_nominal_range")
    for variant in manifest["variants"]:
        for drug in ("pem", "fut"):
            values = [indexed[m, variant["variant_id"], drug] for m in METHODS]
            baseline = values[0]
            if any(any(row.get(f) != baseline.get(f) for f in fields) for row in values):
                raise ValueError("Released methods disagree on measurements")
            if baseline["group_id"] != variant["group_id"]:
                raise ValueError("Released group differs from manifest")
            value = {"value": baseline["residual_activity"],
                     "status": baseline["measurement_status"],
                     "source_cell": baseline["source_cell"]}
            value.update({f: baseline[f] for f in fields[4:] if f in baseline})
            measurements[variant["variant_id"], drug] = value
    return measurements


def score_released(root, reference, output, manifest, freeze_path):
    freeze = verify_commitment(freeze_path)
    prediction_path = output / "external-predictions.json"
    if str(prediction_path) not in freeze["hashes"]:
        raise ValueError("External prediction path is absent from commitment")
    # No workbook is opened. The original completed receipt authenticates released cells.
    receipt_path = root / "data/revision-v2/external-receipt.json"
    verify_hashes(freeze["released_artifact_sha256"])
    receipt = read_json(receipt_path)
    report_path = reference / "external/external-report.json"
    if receipt.get("status") != "completed" or Path(receipt["report_path"]) != report_path:
        raise ValueError("Reference external receipt must be completed for this run")
    verify_hashes({str(report_path): receipt["report_sha256"]})
    report = read_json(report_path)
    records_path = reference / "external/external-records.json"
    if (Path(report["records_path"]) != records_path
            or report["records_sha256"] != receipt["records_sha256"]):
        raise ValueError("Released external record binding differs")
    verify_hashes({str(records_path): receipt["records_sha256"]})
    predictions = _predictions(prediction_path, METHODS, manifest)
    measurements = released_measurements(read_json(records_path)["records"], manifest)
    coverage, metrics, records = _summaries(METHODS, predictions, measurements, manifest)
    destination = output / "external/external-records.json"
    write_json(destination, {"panel_id": manifest["panel_id"], "records": records})
    result = {"status": "completed", "retrospective_reuse": True,
              "panel_id": manifest["panel_id"], "primary_metric": "mean_per_drug_spearman",
              "minimum_n_per_drug": 5, "measurement_coverage": coverage, "methods": metrics,
              "primary_evaluable": all(m["macro_spearman"] is not None for m in metrics.values()),
              "wt_adjustment_applied": False, "classification_metrics_computed": False,
              "records_path": str(destination), "records_sha256": sha256(destination),
              "interpretation": "Retrospective reuse of spent S13 biochemical measurements; "
                                "not fresh blind validation or clinical accuracy."}
    write_json(output / "external/external-report.json", result)
    write_json(output / "retrospective-evaluation.json", {
        "status": "completed", "completed_at": datetime.now(UTC).isoformat(),
        "fresh_blind_validation": False, "original_receipt_unchanged": True,
        "commitment_sha256": sha256(freeze_path),
        "released_artifact_sha256": freeze["released_artifact_sha256"],
        "report_sha256": sha256(output / "external/external-report.json")})
    return result


def prepare_protocol(root, reference, output, preservation_path):
    preserved = read_json(preservation_path)["hashes"]
    verify_hashes(preserved)
    original = read_json(reference / "protocol.json")
    verify_hashes(original["hashes"])
    if original["model"] != REFERENCE_MODEL:
        raise ValueError("Unexpected reference inference model")
    if {k: version(k) for k in original["dependencies"]} != original["dependencies"]:
        raise ValueError("Locked numerical dependency versions changed")
    if original["catalog"] != {k: v.model_dump() for k, v in catalog().items()}:
        raise ValueError("Candidate catalog changed")
    protocol = dict(original)
    protocol.update(created_at=datetime.now(UTC).isoformat(), model=MODEL)
    additions = [Path(__file__), root / "src/escape_lab/revision_uncertainty.py",
                 root / "docs/model-swap-v2-protocol.md"]
    protocol["model_swap"] = {
        "reference_run": str(reference), "reference_protocol_sha256": sha256(reference / "protocol.json"),
        "retrospective_reuse_authorized": True, "fresh_blind_validation": False,
        "changed_original_fields": ["created_at", "model"],
        "preservation_manifest": str(preservation_path),
        "preservation_manifest_sha256": sha256(preservation_path),
        "additional_source_sha256": {str(p): sha256(p) for p in additions},
        "request_settings": {"tool_choice": "required", "max_tokens": 8192, "stream": False},
        "unspecified_in_both": ["temperature", "top_p", "seed", "reasoning_effort", "thinking"],
        "serving_caveat": "Model-specific sampling and reasoning defaults are not controlled; "
                          "remaining_seconds and subsequent feedback histories are dynamic.",
        "scientific_attempts": 1,
    }
    protocol["reference_artifact_sha256"] = {
        p: h for p, h in preserved.items() if Path(p).is_relative_to(reference)}
    write_json(output / "protocol.json", protocol)
    return protocol, preserved


def run(root, reference, output, preservation_path):
    output.mkdir(parents=True, exist_ok=False)
    protocol, preserved = prepare_protocol(root, reference, output, preservation_path)
    panel_path = root / "data/external_discovery/s13-panel.json"
    panel = read_json(panel_path)
    candidates = retained_candidates(read_json(root / "data/processed/candidates.json"), panel)
    if ([c["candidate_id"] for c in candidates] != protocol["training_ids"]
            or sorted({c["group_id"] for c in candidates}) != protocol["training_groups"]):
        raise ValueError("Training pool differs from exact frozen protocol")
    store = BiologicalFeatureStore(root / "data/revision-v2/features")
    labels = load_labels(root / "data/processed", candidates)
    key, _ = credentials(root)
    status = {"status": "started", "started_at": datetime.now(UTC).isoformat(),
              "protocol_sha256": sha256(output / "protocol.json"),
              "retrospective_reuse": True,
              "reference_artifact_sha256": protocol["reference_artifact_sha256"]}
    write_json(output / "run-status.json", status)
    with AuditedClient(key, model_id=MODEL, max_tokens=8192, max_retries=0,
                       max_total_tokens=TOKEN_BUDGET, timeout=60,
                       audit_path=output / "requests.jsonl") as client:
        started = None
        try:
            client.discover_models(nemotron_only=False)
            selected = client.configure_model(MODEL)
            if selected.tools_supported is not True:
                raise ValueError("GLM tool support must be advertised before scientific run")
            raw = client._request("GET", "models")
            write_json(output / "discovery.json", {
                "retrieved_at": datetime.now(UTC).isoformat(),
                "models": [{k: row.get(k) for k in
                            ("id", "pricing", "supported_parameters", "context_length", "top_provider")}
                           for row in raw["data"] if row["id"] in (MODEL, REFERENCE_MODEL)]})
            started = time.monotonic()
            deadline = started + RUN_SECONDS

            def planner_factory():
                return ExperimentPlanner(client, max_proposals=3)

            nested = nested_compare(candidates, labels, store, planner_factory,
                                    output / "nested", deadline)
            check_deadline(deadline)
            search = search_candidates(candidates, labels, store, planner_factory(), 5100,
                                       deadline, output / "full-search")
            statuses = list(nested["agent_status"].values()) + [search["agent_status"]]
            if "unavailable" in statuses:
                raise RuntimeError("At least one search unavailable; external not scored")
            pool = [{**c, "candidate_id": c["variant_id"]} for c in panel["variants"]]
            committed = predict_selected(search["selections"], candidates, labels, pool, store,
                                         5100, output / "models", deadline=deadline)
            prediction_path = output / "external-predictions.json"
            write_json(prediction_path, {"predictions": [
                {"method": r["method"], "variant_id": r["candidate_id"], "drug": r["drug"],
                 "score": r["probability"]} for r in committed["predictions"]]})
            write_json(output / "fitted-models.json", {
                "selections": search["selections"], "model_hashes": committed["model_hashes"],
                "diagnostics": committed["diagnostics"]})
            verify_hashes(preserved)
            freeze_path = output / "external-freeze.json"
            released_paths = [root / "data/revision-v2/external-receipt.json",
                              reference / "external/external-report.json",
                              reference / "external/external-records.json"]
            commitment = {str(prediction_path): sha256(prediction_path),
                          str(output / "protocol.json"): status["protocol_sha256"],
                          **committed["model_hashes"],
                          **protocol["model_swap"]["additional_source_sha256"]}
            selected_paths = [output / "nested/nested-results.json",
                              output / "full-search/search.json", output / "fitted-models.json"]
            selected_paths += list((output / "nested").glob("outer-*/search.json"))
            commitment.update({str(p): sha256(p) for p in selected_paths})
            write_json(freeze_path, {"retrospective_reuse_authorized": True,
                                    "committed_at": datetime.now(UTC).isoformat(),
                                    "selected_configurations": search["selections"],
                                    "hashes": commitment,
                                    "released_artifact_sha256": {str(p): preserved[str(p)]
                                                                 for p in released_paths}})
            print("Predictions committed; scoring already released S13 measurements.", flush=True)
            check_deadline(deadline)
            external = score_released(root, reference, output, panel, freeze_path)
            verify_hashes(preserved)
            verify_commitment(freeze_path)
            result = {"nested_metrics": nested["metrics"], "external": external,
                      "full_selections": search["selections"], "agent_statuses": statuses,
                      "seconds": time.monotonic() - started,
                      "token_accounting": {"budget_charge_tokens": client.total_tokens,
                                           "usage_records": client.usage_records},
                      "interpretation": "Single model substitution; exploratory internal CV and "
                                        "retrospective reuse of spent S13 measurements."}
            write_json(output / "results.json", result)
            status.update(status="completed", results_sha256=sha256(output / "results.json"))
            print(json.dumps({"run_dir": str(output), "seconds": result["seconds"],
                              "tokens": client.total_tokens, "agent_statuses": statuses}), flush=True)
            return result
        except BaseException as error:
            status.update(status="failed", error_type=type(error).__name__)
            raise
        finally:
            status.update(seconds=time.monotonic() - started if started else 0,
                          budget_charge_tokens=client.total_tokens, usage_records=client.usage_records)
            status["artifact_sha256"] = {str(p): sha256(p) for p in output.rglob("*")
                                         if p.is_file() and p.name != "run-status.json"}
            write_json(output / "run-status.json", status)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    root = args.root.resolve()
    run_id = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    output = root / "reports/model-swap-v2" / run_id
    run(root, root / "reports/revision-v2" / REFERENCE_RUN, output,
        root / "reports/model-swap-v2/preservation-before.json")


if __name__ == "__main__":
    main()
