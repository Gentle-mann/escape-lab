"""Run the predeclared v2 comparison and consume the external panel exactly once.

Usage: uv run python -m escape_lab.revision_cli --root .
No option changes the scientific protocol or resets a receipt.
"""

import argparse
import json
import time
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from uuid import uuid4

import httpx

from escape_lab.biological_features import BiologicalFeatureStore
from escape_lab.cli import configured_client
from escape_lab.data import sha256
from escape_lab.experiment_planner import ExperimentPlanner
from escape_lab.external_panel import evaluate_external_panel, inspect_external_panel
from escape_lab.improvement import load_labels, write_json
from escape_lab.revision_evaluation import (
    METHODS,
    check_deadline,
    nested_compare,
    predict_selected,
    retained_candidates,
    search_candidates,
)
from escape_lab.revision_models import catalog

MODEL = "nvidia/NVIDIA-Nemotron-3.5-Super-VL"
RUN_SECONDS = 1200
TOKEN_BUDGET = 180_000


def verify_hashes(hashes):
    for filename, digest in hashes.items():
        path = Path(filename)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"Frozen artifact changed: {path.name}")


def freeze_protocol(root, output_dir, candidates, panel_path):
    # The optional bootstrap module is post-run analysis, never imported by this
    # runner or used for selection. It records its own hash with its report.
    paths = [p for p in (root / "src/escape_lab").glob("*.py")
             if p.name != "revision_uncertainty.py"]
    paths += [root / "pyproject.toml", root / "uv.lock", root / "docs/revision-v2-research.md",
              root / "data/processed/candidates.json", root / "data/processed/outcomes.json",
              root / "data/processed/final-evaluation.json", panel_path]
    paths += list((root / "data/revision-v2/features").glob("*.json"))
    hashes = {str(p.resolve()): sha256(p) for p in sorted(paths)}
    result = {"schema_version": "fgfr-revision-v2", "created_at": datetime.now(UTC).isoformat(),
              "hashes": hashes, "catalog": {key: value.model_dump() for key, value in catalog().items()},
              "training_ids": [c["candidate_id"] for c in candidates],
              "training_groups": sorted({c["group_id"] for c in candidates}),
              "outer_salt": "outer-20261017", "outer_seeds": [4100, 4101, 4102],
              "full_search_seed": 5100, "proposal_opportunities_per_search": 3,
              "max_searches": 4, "max_format_repairs_per_proposal": 1,
              "model": MODEL, "seconds_budget": RUN_SECONDS, "token_stop_threshold": TOKEN_BUDGET,
              "http_timeout_seconds": 60, "http_retries": 0,
              "deadline_semantics": "Checked between bounded fits and HTTP calls, not an OS kill.",
              "external_reveal_requires": "At least one evaluated proposal in each of four searches; failures remain reported.",
              "dependencies": {name: version(name) for name in
                               ("numpy", "scipy", "scikit-learn", "biopython", "pydantic", "openpyxl")}}
    write_json(output_dir / "protocol.json", result)
    return result


def run(root, client, output_dir):
    root, output_dir = Path(root).resolve(), Path(output_dir).resolve()
    receipt = root / "data/revision-v2/external-receipt.json"
    if receipt.exists():
        raise FileExistsError("External evaluation is already reserved/spent; cannot run again")
    panel_path = root / "data/external_discovery/s13-panel.json"
    panel = json.loads(panel_path.read_text())
    # Metadata-only validation; the external evaluator exclusively owns outcome reads.
    inspect_external_panel(panel_path)
    candidates = retained_candidates(json.loads((root / "data/processed/candidates.json").read_text()),
                                     panel)
    if len(candidates) != 3730 or len({c["group_id"] for c in candidates}) != 15:
        raise ValueError("Dataset differs from the predeclared 3730-variant / 15-block pool")
    store = BiologicalFeatureStore(root / "data/revision-v2/features")
    output_dir.mkdir(parents=True, exist_ok=False)
    protocol = freeze_protocol(root, output_dir, candidates, panel_path)
    labels = load_labels(root / "data/processed", candidates)
    started = time.monotonic()
    deadline = started + RUN_SECONDS
    run_status = {"status": "started", "started_at": datetime.now(UTC).isoformat(),
                  "protocol_sha256": sha256(output_dir / "protocol.json")}
    write_json(output_dir / "run-status.json", run_status)
    try:
        def planner_factory():
            return ExperimentPlanner(client, max_proposals=3)

        nested = nested_compare(candidates, labels, store, planner_factory,
                                output_dir / "nested", deadline)
        check_deadline(deadline)
        search = search_candidates(candidates, labels, store, planner_factory(), 5100, deadline,
                                   output_dir / "full-search")
        statuses = list(nested["agent_status"].values()) + [search["agent_status"]]
        if "unavailable" in statuses:
            raise RuntimeError("At least one agent search had no evaluated proposal; external remains unopened")
        pool = [{**c, "candidate_id": c["variant_id"]} for c in panel["variants"]]
        committed = predict_selected(search["selections"], candidates, labels, pool, store, 5100,
                                     output_dir / "models", deadline=deadline)
        predictions = {"predictions": [
            {"method": row["method"], "variant_id": row["candidate_id"],
             "drug": row["drug"], "score": row["probability"]}
            for row in committed["predictions"]]}
        prediction_path = output_dir / "external-predictions.json"
        write_json(prediction_path, predictions)
        write_json(output_dir / "fitted-models.json", {
            "selections": search["selections"], "model_hashes": committed["model_hashes"],
            "diagnostics": committed["diagnostics"]})
        verify_hashes(protocol["hashes"])
        verify_hashes(committed["model_hashes"])
        verify_hashes({str(output_dir / "protocol.json"): run_status["protocol_sha256"]})
        external_freeze = {"protocol_sha256": run_status["protocol_sha256"],
                          "model_hashes": committed["model_hashes"],
                          "selected_configurations": search["selections"],
                          "agent_statuses": statuses,
                          "external_evaluation": {
                              "authorized": True,
                              "panel_manifest_sha256": sha256(panel_path),
                              "document_sha256": panel["document_sha256"],
                              "predictions_sha256": sha256(prediction_path),
                              "evaluator_sha256": sha256(root / "src/escape_lab/external_panel.py"),
                              "methods": list(METHODS),
                              "primary_metric": "mean_per_drug_spearman",
                              "minimum_n_per_drug": 5}}
        freeze_path = output_dir / "external-freeze.json"
        write_json(freeze_path, external_freeze)
        print("All five methods committed predictions; opening independent measurements once.", flush=True)
        check_deadline(deadline)
        external = evaluate_external_panel(panel_path, freeze_path, prediction_path,
                                           receipt_path=receipt, output_dir=output_dir / "external")
        verify_hashes(protocol["hashes"])
        verify_hashes(committed["model_hashes"])
        verify_hashes({str(output_dir / "protocol.json"): run_status["protocol_sha256"]})
        report = {"nested_metrics": nested["metrics"], "external": external,
                  "full_selections": search["selections"], "agent_statuses": statuses,
                  "seconds": time.monotonic() - started,
                  "token_accounting": {"budget_charge_tokens": client.total_tokens,
                                       "usage_records": client.usage_records},
                  "interpretation": "Exploratory internal CV and one-shot biochemical assay transfer; no clinical accuracy claim."}
        write_json(output_dir / "results.json", report)
        run_status.update(status="completed", results_sha256=sha256(output_dir / "results.json"))
        print(json.dumps({"nested_macro_ap": {m: nested["metrics"][m]["macro_ap"] for m in METHODS},
                          "external": external, "full_selections": search["selections"]}, indent=2),
              flush=True)
        return report
    except BaseException as error:
        run_status.update(status="failed", error_type=type(error).__name__)
        raise
    finally:
        run_status.update(seconds=time.monotonic() - started,
                          budget_charge_tokens=client.total_tokens,
                          usage_records=client.usage_records)
        write_json(output_dir / "run-status.json", run_status)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    root = args.root.resolve()
    run_id = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid4().hex[:8]
    with configured_client(root, max_total_tokens=TOKEN_BUDGET, model_id=MODEL) as client:
        client.max_tokens = 8192
        client.max_retries = 0
        client._http.timeout = httpx.Timeout(60)
        run(root, client, root / "reports/revision-v2" / run_id)


if __name__ == "__main__":
    main()
