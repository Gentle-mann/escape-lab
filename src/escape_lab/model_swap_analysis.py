"""Post-run retrospective GLM/Nemotron comparison of already released v2 outcomes.

No model selection, model fitting, source workbook reads, or inference calls occur
here. The existing five-method validator and metric are reused without modifying
their globals; a sixth score vector is added only after both runs are validated.
"""

import argparse
import json
import math
import shlex
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from escape_lab import revision_uncertainty
from escape_lab.data import sha256
from escape_lab.improvement import write_json
from escape_lab.revision_evaluation import summarize
from escape_lab.revision_uncertainty import DRUGS, REPLICATES, SEED, _metric, _paired_arrays

BASELINES = ("v1_fixed", "biology_fixed", "random", "grid")
METHODS = (*BASELINES, "nemotron", "glm")
CONTRASTS = tuple((a.replace("agent", "nemotron"), b.replace("agent", "nemotron"))
                  for a, b in revision_uncertainty.CONTRASTS) + tuple(
                      ("glm", method) for method in (*BASELINES, "nemotron"))
SEARCHES = ("outer-0", "outer-1", "outer-2", "full")
PRICES = {"nemotron": {"input": .10, "output": .40},
          "glm": {"input": 1.40, "output": 4.40}}
ABS_TOL = 1e-12


def assert_equivalent(left, right, context):
    """Recursively compare reports, allowing only tiny absolute float noise."""
    if isinstance(left, dict) and isinstance(right, dict):
        if set(left) != set(right):
            raise ValueError(f"{context}: keys differ")
        for key in left:
            assert_equivalent(left[key], right[key], f"{context}.{key}")
    elif isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            raise ValueError(f"{context}: lengths differ")
        for i, (a, b) in enumerate(zip(left, right, strict=True)):
            assert_equivalent(a, b, f"{context}[{i}]")
    elif type(left) in (int, float) and type(right) in (int, float):
        if not math.isclose(left, right, rel_tol=0, abs_tol=ABS_TOL):
            raise ValueError(f"{context}: numerical mismatch ({left!r}, {right!r})")
    elif left != right or type(left) is not type(right):
        raise ValueError(f"{context}: values differ")


def _indexed(records, endpoint, method):
    identity = "candidate_id" if endpoint == "internal" else "variant_id"
    result = {}
    for row in records:
        if row["method"] != method:
            continue
        key = (row[identity], row["drug"])
        if key in result:
            raise ValueError("Duplicate prediction record")
        result[key] = row
    return result


def compare_predictions(reference, glm, endpoint):
    """Validate within-run pairing, then require identical numerical baselines."""
    ref_groups, ref_arrays = _paired_arrays(reference, endpoint)
    glm_groups, glm_arrays = _paired_arrays(glm, endpoint)
    if ref_groups != glm_groups:
        raise ValueError("Cross-run paired groups differ")
    for method in BASELINES:
        left, right = _indexed(reference, endpoint, method), _indexed(glm, endpoint, method)
        assert_equivalent(left, right, f"{endpoint} baseline {method}")
    for drug in DRUGS:
        for field in ("truth", "groups"):
            if not np.array_equal(ref_arrays[drug][field], glm_arrays[drug][field]):
                raise ValueError("Cross-run paired outcomes or group assignments differ")
    # Internal fold identity must also match for the two LLM-selected predictors.
    if endpoint == "internal":
        left = _indexed(reference, endpoint, "agent")
        right = _indexed(glm, endpoint, "agent")
        for key in left:
            if left[key].get("fold") != right[key].get("fold"):
                raise ValueError("Cross-run agent outer folds differ")
    return ref_groups, ref_arrays, glm_arrays


def paired_model_swap_bootstrap(reference, glm, endpoint, replicates=REPLICATES, seed=SEED):
    """Same v2 paired group draws and metric, with direct GLM contrasts added."""
    if type(replicates) is not int or replicates < 1:
        raise ValueError("Positive replicate count required")
    groups, arrays, glm_arrays = compare_predictions(reference, glm, endpoint)
    for drug in DRUGS:
        scores = arrays[drug]["scores"]
        scores["nemotron"] = scores.pop("agent")
        scores["glm"] = glm_arrays[drug]["scores"]["agent"]
    points = {method: _metric(arrays, method, np.ones(len(groups), dtype=int), endpoint)
              for method in METHODS}
    draws = {contrast: [] for contrast in CONTRASTS}
    rng = np.random.default_rng(seed)
    for _ in range(replicates):
        counts = np.bincount(rng.integers(len(groups), size=len(groups)), minlength=len(groups))
        scores = {method: _metric(arrays, method, counts, endpoint) for method in METHODS}
        for (method, comparator), differences in draws.items():
            if scores[method] is not None and scores[comparator] is not None:
                differences.append(scores[method] - scores[comparator])
    comparisons = {}
    for (method, comparator), differences in draws.items():
        valid = points[method] is not None and points[comparator] is not None
        comparisons[f"{method}_minus_{comparator}"] = {
            "method": method, "comparator_method": comparator,
            "estimate": points[method], "comparator_estimate": points[comparator],
            "delta": points[method] - points[comparator] if valid else None,
            "percentile_2_5": float(np.percentile(differences, 2.5)) if differences else None,
            "percentile_97_5": float(np.percentile(differences, 97.5)) if differences else None,
            "valid_replicates": len(differences),
            "invalid_replicates": replicates - len(differences),
        }
    return {
        "endpoint": endpoint,
        "metric": "macro_average_precision" if endpoint == "internal"
                  else "mean_per_drug_spearman",
        "aligned_group_count": len(groups), "aligned_groups": groups,
        "observed_pairs_per_drug": {drug: len(arrays[drug]["truth"]) for drug in DRUGS},
        "replicates": replicates, "seed": seed, "point_estimates": points,
        "comparisons": comparisons,
        "interpretation": "Descriptive paired block-bootstrap percentiles conditional on committed "
                          "predictions. Related sequence blocks are not independent biological "
                          "replicates. No refitting, model selection, p-values, or significance claim.",
    }


def _read(path):
    return json.loads(Path(path).read_text())


def _verified(path, digest, inputs):
    path = Path(path).resolve()
    if not isinstance(digest, str) or sha256(path) != digest:
        raise ValueError(f"Released artifact checksum mismatch: {path.name}")
    inputs[str(path)] = digest
    return _read(path)


def _manifest_digest(manifest, path, run_dir):
    return manifest.get(str(path)) or manifest.get(str(path.relative_to(run_dir)))


def load_completed_run(run_dir, reference_artifact_sha256=None):
    """Verify completed status and released-artifact chains before scoring.

    Old v2 nested records are bound by its existing uncertainty report; new runs
    bind every analysis input in run-status.artifact_sha256. Original searches
    can additionally be bound by the model-swap's pre-run reference snapshot.
    Only released run artifacts and their completed receipt are read.
    """
    run_dir = Path(run_dir).resolve()
    status_path = run_dir / "run-status.json"
    status = _read(status_path)
    if status.get("status") != "completed":
        raise ValueError("Run must be completed before analysis")
    inputs = {str(status_path): sha256(status_path)}
    if reference_artifact_sha256:
        expected = _manifest_digest(reference_artifact_sha256, status_path, run_dir)
        if expected:
            _verified(status_path, expected, inputs)
    results = _verified(run_dir / "results.json", status.get("results_sha256"), inputs)
    protocol = _verified(run_dir / "protocol.json", status.get("protocol_sha256"), inputs)
    manifest = dict(status.get("artifact_sha256", {}))
    manifest.update(reference_artifact_sha256 or {})
    limitations = []
    uncertainty_path = run_dir / "uncertainty.json"
    uncertainty = None
    if uncertainty_path.exists():
        digest = _manifest_digest(manifest, uncertainty_path, run_dir)
        uncertainty = _verified(uncertainty_path, digest or sha256(uncertainty_path), inputs)
        if uncertainty.get("source_sha256") != sha256(Path(revision_uncertainty.__file__)):
            raise ValueError("Original uncertainty implementation checksum differs")
        for name, digest in uncertainty["input_sha256"].items():
            path = Path(name).resolve()
            if path.is_relative_to(run_dir):
                _verified(path, digest, inputs)
                manifest[str(path)] = digest
            elif path.name == "external-receipt.json":
                receipt = _verified(path, digest, inputs)
                if receipt.get("status") != "completed":
                    raise ValueError("Original external receipt must be completed")
                report_path = run_dir / "external/external-report.json"
                if Path(receipt["report_path"]).resolve() != report_path:
                    raise ValueError("External receipt belongs to another run")
                _verified(report_path, receipt["report_sha256"], inputs)
                manifest[str(report_path)] = receipt["report_sha256"]
                manifest[str(run_dir / "external/external-records.json")] = (
                    receipt["records_sha256"])
                for field in ("predictions", "frozen_protocol"):
                    binding = receipt["bound_files"][field]
                    bound = Path(binding["path"]).resolve()
                    if not bound.is_relative_to(run_dir):
                        raise ValueError("Receipt binds an artifact outside its run")
                    _verified(bound, binding["sha256"], inputs)
                    manifest[str(bound)] = binding["sha256"]
    loaded = {}
    required = ["nested/nested-results.json", "external/external-report.json",
                "external/external-records.json", "external-predictions.json",
                "external-freeze.json"]
    search_paths = {search: ("full-search/search.json" if search == "full"
                            else f"nested/{search}/search.json") for search in SEARCHES}
    for relative in required + list(search_paths.values()):
        path = run_dir / relative
        digest = _manifest_digest(manifest, path, run_dir)
        if not digest:
            if uncertainty is not None and relative in search_paths.values():
                digest = sha256(path)
                limitations.append(f"Original {relative} lacked a historical completion hash; "
                                   "current bytes and cross-artifact consistency were checked.")
            else:
                raise ValueError(f"Missing completed artifact checksum: {relative}")
        loaded[relative] = _verified(path, digest, inputs)
    report = loaded["external/external-report.json"]
    records_path = run_dir / "external/external-records.json"
    if report.get("status") != "completed":
        raise ValueError("External report must be completed")
    if Path(report["records_path"]).resolve() != records_path:
        raise ValueError("External report records belong to another run")
    _verified(records_path, report["records_sha256"], inputs)
    assert_equivalent(results["external"], report, "result/report external agreement")
    nested = loaded["nested/nested-results.json"]
    assert_equivalent(results["nested_metrics"], nested["metrics"], "nested/result metrics")
    frozen = loaded["external-freeze.json"]
    assert_equivalent(frozen["selected_configurations"], results["full_selections"],
                      "frozen/result selections")
    predictions_path = run_dir / "external-predictions.json"
    if frozen.get("retrospective_reuse_authorized") is True:
        commitments = frozen["hashes"]
        for path in (run_dir / "protocol.json", predictions_path):
            if commitments.get(str(path)) != sha256(path):
                raise ValueError("Retrospective commitment checksum differs")
        for name, digest in commitments.items():
            if sha256(Path(name)) != digest:
                raise ValueError("Retrospective commitment artifact changed")
            inputs[str(Path(name).resolve())] = digest
    else:
        if frozen["protocol_sha256"] != status["protocol_sha256"]:
            raise ValueError("External freeze protocol checksum differs")
        if frozen["external_evaluation"]["predictions_sha256"] != sha256(predictions_path):
            raise ValueError("Committed external predictions checksum differs")
    committed = {(r["method"], r["variant_id"], r["drug"]): r["score"]
                 for r in loaded["external-predictions.json"]["predictions"]}
    records = loaded["external/external-records.json"]["records"]
    scored = {(r["method"], r["variant_id"], r["drug"]): r["prediction_score"] for r in records}
    if len(committed) != len(loaded["external-predictions.json"]["predictions"]):
        raise ValueError("Duplicate committed external predictions")
    assert_equivalent(committed, scored, "committed/scored external predictions")
    searches = {search: loaded[path] for search, path in search_paths.items()}
    for search, saved in searches.items():
        selected = results["full_selections"] if search == "full" else nested["selections"][search]
        assert_equivalent(saved["selections"], selected, f"{search} search/result selections")
    return {"path": str(run_dir), "status": status, "protocol": protocol, "results": results,
            "nested": nested, "external_records": records, "searches": searches,
            "uncertainty": uncertainty, "input_sha256": inputs, "limitations": limitations}


def _verify_metrics(run):
    for method in (*BASELINES, "agent"):
        internal = [row for row in run["nested"]["predictions"] if row["method"] == method]
        assert_equivalent(summarize(internal), run["results"]["nested_metrics"][method],
                          f"{method} recomputed internal metrics")
        correlations = []
        for drug in DRUGS:
            rows = [row for row in run["external_records"] if row["method"] == method
                    and row["drug"] == drug and row["measurement_status"] == "observed"]
            expected = run["results"]["external"]["methods"][method]["per_drug"][drug]
            truth = [row["residual_activity"] for row in rows]
            scores = [row["prediction_score"] for row in rows]
            value = None
            if len(rows) >= 5 and len(set(truth)) > 1 and len(set(scores)) > 1:
                value = float(spearmanr(scores, truth).statistic)
                if not math.isfinite(value):
                    value = None
            assert_equivalent(value, expected["spearman"], f"{method}/{drug} external Spearman")
            correlations.append(value)
            assert_equivalent(len(rows), expected["n"], f"{method}/{drug} external N")
            assert_equivalent(len({r["group_id"] for r in rows}), expected["aligned_groups"],
                              f"{method}/{drug} external groups")
        macro = sum(correlations) / len(DRUGS) if all(v is not None for v in correlations) else None
        assert_equivalent(macro, run["results"]["external"]["methods"][method]["macro_spearman"],
                          f"{method} external macro Spearman")


def estimate_cost(accounting, model):
    rates = PRICES.get(model, None)
    if rates is None:
        raise ValueError("Unknown model price")
    records = accounting["usage_records"]
    complete = all(type(r.get(k)) is int and r[k] >= 0
                   for r in records for k in ("prompt_tokens", "completion_tokens"))
    prompt = sum(r["prompt_tokens"] for r in records) if complete else None
    completion = sum(r["completion_tokens"] for r in records) if complete else None
    return {"usd": ((prompt * rates["input"] + completion * rates["output"]) / 1_000_000
                    if complete else None),
            "prompt_tokens": prompt, "completion_tokens": completion,
            "input_usd_per_million": rates["input"], "output_usd_per_million": rates["output"],
            "cached_discount_applied": False, "invoice_verified": False,
            "usage_estimated": any(r.get("usage_estimated", False) for r in records),
            "source_url": "https://www.crusoe.ai/cloud/pricing", "retrieved_date": "2026-10-17",
            "basis": "Uncached list-price estimate from recorded request input/output tokens; "
                     "not an invoice. Provider billing, caching and failed-request charging may differ."}


def _run_summary(run, model):
    searches, results = run["searches"], run["results"]
    attempts = [r for s in searches.values() for r in s["agent_attempts"]]
    proposals = [p for s in searches.values() for p in s["agent_proposals"]]
    accounting = results["token_accounting"]
    return {"model_id": run["protocol"]["model"], "seconds": results["seconds"],
            "status_seconds": run["status"].get("seconds"),
            "proposal_opportunities": len(proposals),
            "valid_evaluated_proposals": sum(p.get("status") == "proposed" and "metrics" in p
                                             for p in proposals),
            "accepted_proposals": sum(p.get("accepted", False) for p in proposals),
            "failed_opportunities": sum(s["failed_opportunities"] for s in searches.values()),
            "planner_attempts": len(attempts),
            "format_repair_requests": sum(r.get("attempt", 1) > 1 for r in attempts),
            "recorded_requests": len(accounting["usage_records"]),
            "budget_charge_tokens": accounting["budget_charge_tokens"],
            "reported_total_tokens": sum(r.get("total_tokens", 0) or 0
                                         for r in accounting["usage_records"]),
            "request_statuses": [r.get("request_status") for r in accounting["usage_records"]],
            "proposal_status_counts": dict(Counter(p.get("status", "unknown") for p in proposals)),
            "proposal_error_counts": dict(Counter(str(p["error"]) for p in proposals
                                                  if p.get("error") is not None)),
            "attempt_status_counts": dict(Counter(r.get("status", "unknown") for r in attempts)),
            "agent_statuses": results["agent_statuses"],
            "searches": {search: {"agent_status": value.get("agent_status"),
                                  "valid_proposals": value["successful_proposals"],
                                  "failed_opportunities": value["failed_opportunities"],
                                  "proposal_statuses": [p.get("status")
                                                        for p in value["agent_proposals"]],
                                  "proposal_errors": [p.get("error")
                                                      for p in value["agent_proposals"]]}
                         for search, value in searches.items()},
            "selections": {s: value["selections"]["agent"] for s, value in searches.items()},
            "cost": estimate_cost(accounting, model)}


def _value(value):
    if value is None:
        return "undefined"
    if type(value) is int:
        return str(value)
    return f"{value:.6f}"


def render_summary(payload):
    lines = ["# GLM-5.3 model-swap retrospective comparison", "",
             ("This is a single matched model-swap run on reused v2 measurements, not a fresh blind "
             "validation. The primary internal endpoint is pooled macro AP; mean fold AP is a "
             "separate descriptive endpoint. Higher AP/Spearman and lower Brier are better."), "",
             f"Nemotron reference: `{payload['reference_run']}`", "",
             f"GLM run: `{payload['glm_run']}`", ""]
    methods = payload["methods"]
    if methods:
        columns = list(METHODS)
        header = "| Endpoint | " + " | ".join(columns) + " |"
        rule = "|---|" + "---:|" * len(columns)
        lines.extend([header, rule])
        specs = [("Pooled macro AP", ("internal", "macro_ap")),
                 ("Macro Brier", ("internal", "macro_brier")),
                 ("Mean fold macro AP", ("internal", "mean_fold_macro_ap")),
                 ("Fold AP std", ("internal", "fold_ap_std")),
                 ("External macro Spearman", ("external", "macro_spearman"))]
        for drug in DRUGS:
            for metric in ("average_precision", "brier", "n", "positive_count", "inconclusive_count"):
                specs.append((f"Internal {drug} {metric}", ("internal", "per_drug", drug, metric)))
            for metric in ("spearman", "n", "aligned_groups"):
                specs.append((f"External {drug} {metric}", ("external", "per_drug", drug, metric)))
        for fold in SEARCHES[:3]:
            for metric in ("macro_ap", "macro_brier"):
                specs.append((f"{fold} {metric}", ("internal", "by_fold", fold, metric)))
            for drug in DRUGS:
                for metric in ("average_precision", "brier"):
                    specs.append((f"{fold} {drug} {metric}",
                                  ("internal", "by_fold", fold, "per_drug", drug, metric)))
        for label, path in specs:
            values = []
            for method in columns:
                value = methods[method]
                for part in path:
                    value = value[part]
                values.append(_value(value))
            lines.append(f"| {label} | " + " | ".join(values) + " |")
        lines.extend(["", "Selected configurations:", "",
                      "| Method | outer-0 | outer-1 | outer-2 | Full search |",
                      "|---|---|---|---|---|"])
        for method in METHODS:
            lines.append(f"| {method} | " + " | ".join(
                methods[method]["selections"][search] for search in SEARCHES) + " |")
    lines.extend(["", "Observed agent budgets and usage:", "",
                  ("| Run | Valid / opportunities | Requests / repairs | Input / output tokens | "
                  "Budget tokens | Seconds | Estimated USD |"), "|---|---:|---:|---:|---:|---:|---:|"])
    for name, run in payload["runs"].items():
        cost = run["cost"]
        lines.append(f"| {name} | {run['valid_evaluated_proposals']} / "
                     f"{run['proposal_opportunities']} | {run['recorded_requests']} / "
                     f"{run['format_repair_requests']} | {cost['prompt_tokens']} / "
                     f"{cost['completion_tokens']} | {run['budget_charge_tokens']} | "
                     f"{run['seconds']:.3f} | {_value(cost['usd'])} |")
    for name, run in payload["runs"].items():
        lines.extend(["", (f"{name} search statuses (outer-0, outer-1, outer-2, full): "
                      f"`{json.dumps(run['agent_statuses'])}`. Proposal statuses: "
                      f"`{json.dumps(run['proposal_status_counts'])}`. Failures: "
                      f"`{json.dumps(run['proposal_error_counts'])}`. "
                      "Failed or truncated responses consumed their original opportunity; "
                      "the recorded format-repair policy was retained.")])
    lines.extend(["", ("Costs use uncached Crusoe list prices retrieved 2026-10-17: Nemotron "
                  "$0.10 input / $0.40 output and GLM $1.40 input / $4.40 output per million "
                  "tokens. These are estimates, not invoices; caching discounts are not applied. "
                  "[Pricing source](https://www.crusoe.ai/cloud/pricing)."), "",
                  ("Paired descriptive 95% block-bootstrap percentile intervals (2,000 identical "
                  "group draws, seed 6100), conditional on committed predictions:"), "",
                  "| Endpoint | Contrast | Delta | 2.5% | 97.5% | Valid / invalid draws |",
                  "|---|---|---:|---:|---:|---:|"])
    for endpoint, result in payload["uncertainty"].items():
        for contrast, values in result["comparisons"].items():
            if not contrast.startswith("glm_minus_"):
                continue
            lines.append(f"| {endpoint} | {contrast} | {_value(values['delta'])} | "
                         f"{_value(values['percentile_2_5'])} | "
                         f"{_value(values['percentile_97_5'])} | "
                         f"{values['valid_replicates']} / {values['invalid_replicates']} |")
    lines.extend(["", ("The four numerical baselines, their search choices, partitions, committed "
                  "predictions and metrics were checked for agreement. All internal metrics and "
                  "external per-drug Spearman/N were recomputed from released predictions. "
                  "No refitting, fresh holdout access, outcome-based selection of runs, or "
                  "statistical significance claim is made. Intervals omit fitting and LLM "
                  "selection uncertainty; related blocks are not independent biological "
                  "replicates. One run cannot establish a general model-quality effect."), ""])
    serving = payload["audit"].get("model_swap_protocol", {})
    if serving:
        lines.extend([("The original protocol fields changed only in model ID and creation time. "
                      "The comparison adds its provenance and retrospective authorization. "
                      "Agent instructions, tools, numerical evaluator, folds, data, features, "
                      "seeds, selection guards, opportunity/repair budgets and scientific limits "
                      "retain the reference protocol."), "",
                      (f"Explicit request settings: `{json.dumps(serving.get('request_settings'))}`. "
                      f"Unspecified in both: {', '.join(serving.get('unspecified_in_both', []))}. "
                      f"{serving.get('serving_caveat', '')}"), "",
                      ("Runtime is measured local wall time, including instrumentation and "
                      "any concurrent local work; it is not a controlled serving-speed benchmark."), ""])
    for limitation in payload["audit"]["limitations"]:
        lines.append(f"- {limitation}")
    lines.extend(["", "Reproduce analysis after scientific completion:", "", "```sh",
                  payload["reproduce_command"], "```", "",
                  ("Exact inputs, hashes, counts, all contrasts and full metrics are in "
                  "`comparison.json`. Existing scientific artifacts are not rewritten."), ""])
    return "\n".join(lines)


def run_comparison(reference_dir, glm_dir):
    reference_dir, glm_dir = Path(reference_dir).resolve(), Path(glm_dir).resolve()
    if reference_dir == glm_dir:
        raise ValueError("Reference and comparison runs must differ")
    targets = [glm_dir / "comparison.json", glm_dir / "results-summary.md"]
    if any(path.exists() for path in targets):
        raise FileExistsError("Comparison outputs already exist")
    glm = load_completed_run(glm_dir)
    reference_hashes = glm["status"].get("reference_artifact_sha256")
    if reference_hashes is None:
        reference_hashes = glm["protocol"].get("reference_artifact_sha256")
    reference = load_completed_run(reference_dir, reference_hashes)
    if reference["protocol"]["model"] != "nvidia/NVIDIA-Nemotron-3.5-Super-VL":
        raise ValueError("Expected exact Nemotron reference model")
    if glm["protocol"]["model"] != "zai-org/GLM-5.3":
        raise ValueError("Expected exact GLM-5.3 comparison model")
    for key, value in reference["protocol"].items():
        if key not in ("model", "created_at"):
            assert_equivalent(value, glm["protocol"].get(key), f"matched protocol {key}")
    assert_equivalent(reference["nested"]["outer_groups"], glm["nested"]["outer_groups"],
                      "outer group partitions")
    for search in SEARCHES:
        left, right = reference["searches"][search], glm["searches"][search]
        for key in ("random_candidates", "inner_groups", "inner_metrics"):
            assert_equivalent(left[key], right[key], f"{search} baseline {key}")
        for method in BASELINES:
            assert_equivalent(left["selections"][method], right["selections"][method],
                              f"{search} baseline {method} selection")
        for run in (reference, glm):
            opportunities = len(run["searches"][search]["agent_proposals"])
            if opportunities != run["protocol"]["proposal_opportunities_per_search"]:
                raise ValueError("Proposal opportunities differ from frozen protocol")
    for run in (reference, glm):
        _verify_metrics(run)
    uncertainty = {}
    for endpoint in ("internal", "external"):
        ref_records = (reference["nested"]["predictions"] if endpoint == "internal"
                       else reference["external_records"])
        glm_records = (glm["nested"]["predictions"] if endpoint == "internal"
                       else glm["external_records"])
        uncertainty[endpoint] = paired_model_swap_bootstrap(ref_records, glm_records, endpoint)
    methods = {}
    for method in METHODS:
        run = glm if method == "glm" else reference
        source_method = "agent" if method in ("nemotron", "glm") else method
        methods[method] = {"internal": run["results"]["nested_metrics"][source_method],
                           "external": run["results"]["external"]["methods"][source_method],
                           "selections": {s: v["selections"][source_method]
                                          for s, v in run["searches"].items()}}
        for endpoint, key in (("internal", "macro_ap"), ("external", "macro_spearman")):
            assert_equivalent(uncertainty[endpoint]["point_estimates"][method],
                              methods[method][endpoint][key], f"{method} bootstrap point estimate")
    # Reproducing old contrasts also verifies the same draw sequence and metric code.
    if reference["uncertainty"]:
        for endpoint in ("internal", "external"):
            for name, contrast in reference["uncertainty"][endpoint]["comparisons"].items():
                translated = {key: value.replace("agent", "nemotron")
                              if isinstance(value, str) else value for key, value in contrast.items()}
                assert_equivalent(uncertainty[endpoint]["comparisons"][name.replace(
                    "agent", "nemotron")], translated, f"original {endpoint} bootstrap contrast")
    inputs = {**reference["input_sha256"], **glm["input_sha256"]}
    for name, digest in inputs.items():
        if sha256(Path(name)) != digest:
            raise ValueError("Input artifact changed during analysis")
    payload = {"analysis_only": True, "retrospective_reused_test": True,
               "reference_run": str(reference_dir), "glm_run": str(glm_dir),
               "methods": methods, "runs": {"nemotron": _run_summary(reference, "nemotron"),
                                               "glm": _run_summary(glm, "glm")},
               "uncertainty": uncertainty, "input_sha256": inputs,
               "source_sha256": {str(Path(__file__).resolve()): sha256(Path(__file__)),
                                 str(Path(revision_uncertainty.__file__).resolve()):
                                     sha256(Path(revision_uncertainty.__file__))},
               "audit": {"baseline_absolute_tolerance": ABS_TOL,
                         "baseline_predictions_and_metrics_verified": True,
                         "original_bootstrap_reproduced": reference["uncertainty"] is not None,
                         "original_protocol_changed_fields": ["model", "created_at"],
                         "model_swap_protocol": glm["protocol"].get("model_swap", {}),
                         "limitations": reference["limitations"] + glm["limitations"]},
               "reproduce_command": "uv run python -m escape_lab.model_swap_analysis "
                                    f"{shlex.quote(str(reference_dir))} {shlex.quote(str(glm_dir))}"}
    summary = render_summary(payload)
    write_json(targets[0], payload)
    targets[1].write_text(summary)
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_run", type=Path)
    parser.add_argument("glm_run", type=Path)
    args = parser.parse_args()
    payload = run_comparison(args.reference_run, args.glm_run)
    print(json.dumps({"comparison": str(args.glm_run / "comparison.json"),
                      "runs": payload["runs"]}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
