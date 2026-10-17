"""Publish aggregate artifacts without changing frozen research originals.

Never add predictions, requests, outcomes, models, workbooks, or credentials.
"""

import hashlib
import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DESTINATION = ROOT / "docs/archive"
RESULT_KEYS = {
    "metrics", "activation_metrics", "resistance_metrics", "gene_prior_metrics",
    "component_gate", "resistance_gate", "gate", "gates", "uncertainty",
    "baseline_reproduction", "v1_reproduction", "unbounded_reproduction",
    "fit_count", "fitting_seconds", "elapsed_seconds", "feature_extraction_seconds",
    "protocol_sha256", "integrity", "post_result_sensitivity", "coverage", "n",
    "status", "comparisons", "prototype_margin_passed", "interpretation",
    "nested_metrics", "external", "full_selections", "agent_statuses", "seconds",
    "token_accounting", "arms", "proposal_opportunities", "retrospective_reuse",
    "fresh_blind_validation", "limits", "numerical_controls_reproduced",
    "joint_commitment_sha256", "total_wall_seconds", "scientific_arm_seconds",
    "preserved_artifacts", "total_requests", "budget_charge_tokens",
}
# Recursive schemas reject unknown nested fields and nonaggregate value types.
# A dict describes allowed fields, a one-element list describes list elements,
# and a tuple of Python types describes a scalar. Missing fields stay omitted.
NUMBER = (int, float, type(None))
INTEGER = (int,)
TEXT = (str,)
BOOLEAN = (bool,)
METHODS = ("control", "relevant", "placebo", "fixed", "random", "numerical", "grid")
CONTRASTS = ("relevant_minus_control", "relevant_minus_placebo", "control_minus_fixed",
             "control_minus_random", "control_minus_numerical")
DRUG_METRIC = {k: NUMBER for k in ("average_precision", "brier", "prevalence")}
DRUG_METRIC.update({k: INTEGER for k in ("n", "positive_count", "inconclusive_count")})
FOLD_METRIC = {"macro_ap": NUMBER, "macro_brier": NUMBER,
               "per_drug": {drug: DRUG_METRIC for drug in ("pem", "fut")}}
SUMMARY = {**FOLD_METRIC, "mean_fold_macro_ap": NUMBER, "fold_ap_std": NUMBER,
           "by_fold": {f"outer-{i}": FOLD_METRIC for i in range(3)}}
EXPERIMENT = {"status": TEXT, "required_contrasts": [TEXT], "interpretation": TEXT,
              "contrast_directions": {name: TEXT for name in CONTRASTS}}
EXPERIMENTS = {name: EXPERIMENT for name in ("research", "agent_selection")}
GUARD = {"ap_delta": NUMBER, "ap_noninferiority": (bool, type(None)), "brier": NUMBER,
         "brier_limit": NUMBER, "brier_noninferiority": (bool, type(None))}
COMPARISON = {
    **{k: TEXT for k in ("method", "comparator_method", "direction")},
    **{k: NUMBER for k in ("estimate", "comparator_estimate", "delta")},
    "interval_95": [NUMBER], "interval_99_bonferroni": [NUMBER],
    "valid_replicates": INTEGER, "invalid_replicates": INTEGER,
    "repeat_only": {"delta_sd": NUMBER, "wins": NUMBER, "ties": NUMBER, "losses": NUMBER,
                    "deltas": [NUMBER], "interpretation": TEXT},
    "guards": {drug: GUARD for drug in ("pem", "fut")},
    "passes_improvement_rule": BOOLEAN, "failed_rules": [TEXT],
}
USAGE = {k: INTEGER for k in (
    "searches", "successful_proposals", "failed_opportunities", "requests", "tokens",
    "prompt_tokens", "completion_tokens", "unknown_usage_requests")}
USAGE.update(seconds_sum=NUMBER, estimated_uncached_usd=NUMBER,
             statuses={k: INTEGER for k in ("proposed", "truncated_response", "failed")})
RESULT_SCHEMA = {
    **{k: TEXT for k in ("schema_version", "metric", "run_status", "joint_freeze_sha256")},
    **{k: INTEGER for k in ("repeat_count", "aligned_group_count", "candidate_count",
                            "protected_artifacts")},
    "methods": {m: {"mean": SUMMARY, "per_repeat": {str(i): SUMMARY for i in range(10)}}
                for m in METHODS},
    "comparisons": {name: COMPARISON for name in CONTRASTS}, "experiments": EXPERIMENTS,
    "bootstrap": {
        **{k: INTEGER for k in ("replicates", "seed", "unique_drug_score_vectors",
                                "unique_method_repeat_patterns")},
        **{k: BOOLEAN for k in ("group_draws_shared_across_methods",
                                "repeat_draws_shared_across_methods")},
        **{k: TEXT for k in ("resampling", "draw_sha256", "intervals")},
    },
    "improvement_rule": {"minimum_absolute_macro_ap_gain": NUMBER,
                         "99_percent_lower_bound": TEXT,
                         "maximum_per_drug_ap_decrease": NUMBER, "maximum_per_drug_brier": TEXT},
    "limitations": [TEXT], "usage": {arm: USAGE for arm in METHODS[:3]},
    "total_wall_seconds": NUMBER, "fresh_external_validation": BOOLEAN,
}
CONFIGS = {f"{rep}-logistic-c{penalty}" for rep in ("chemistry", "structure", "interactions")
           for penalty in ("0.01", "0.1", "1")}
CONFIGS.update(("structure-hgb-depth2", "structure-hgb-depth3", "interactions-hgb-depth2"))
AUDIT_SCHEMA = {
    **{k: TEXT for k in ("schema_version", "verified_at", "run_status", "report_sha256")},
    **{k: INTEGER for k in ("scientific_attempts", "valid_proposals", "failed_opportunities",
                            "known_provider_usage_tokens", "unknown_usage_requests")},
    "experiments": EXPERIMENTS,
    "verification": {
        **{k: INTEGER for k in ("run_artifact_hashes", "joint_commitment_hashes",
                                "protocol_input_hashes", "preserved_prior_hashes",
                                "credential_scan_files", "credential_scan_hits")},
        **{k: BOOLEAN for k in ("all_hashes_valid", "protocol_hash_matches",
                                "results_hash_matches", "joint_freeze_hash_matches")},
    },
    "independent_request_audit": {
        "source": TEXT, "tools": [TEXT],
        **{k: INTEGER for k in ("searches", "requests", "per_arm_requests", "own_history_checks",
                                "adapter_fingerprint_pairs", "repairs")},
        **{k: BOOLEAN for k in ("passed", "outcomes_in_proposal_inputs", "cross_arm_history",
                                "raw_source_text", "all_requests_finished_before_commitment",
                                "scored_outcomes_inspected_by_isolation_auditor")},
    },
    "independent_results_review": {k: TEXT for k in (
        "source", "metric_and_selection_consistency", "both_experiments")},
    "choices_by_method": {m: {c: INTEGER for c in CONFIGS} for m in METHODS},
    "relevant_choice_disagreements": {m: INTEGER for m in ("control", "placebo")},
    "failure_types": {k: INTEGER for k in ("truncated_response", "connection_error")},
    "known_usage_estimated_inference_usd": NUMBER, "new_you_retrieval_estimated_usd": NUMBER,
    "complete_invoice": BOOLEAN, "fresh_external_validation": BOOLEAN,
}
REPEATED_SCHEMAS = {
    "reports/repeated-research/20261017-experiments-1-2/results.json": RESULT_SCHEMA,
    "reports/repeated-research/20261017-experiments-1-2-audit.json": AUDIT_SCHEMA,
}
REPEATED_JSON_FIELDS = {name: set(schema) for name, schema in REPEATED_SCHEMAS.items()}


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def clean_text(value):
    return value.replace(str(ROOT) + "/", "").replace(str(ROOT), ".")


def clean(value):
    if isinstance(value, dict):
        return {clean_text(k): clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean(v) for v in value]
    if isinstance(value, str):
        return clean_text(value)
    return value


def require_aggregate_content(value, schema):
    """Accept only the reviewed aggregate shape; unknown nested content fails closed."""
    error = "Publication does not match the reviewed aggregate schema"
    if isinstance(schema, dict):
        if not isinstance(value, dict) or set(value) - set(schema):
            raise ValueError(error)
        for key, child in value.items():
            require_aggregate_content(child, schema[key])
    elif isinstance(schema, list):
        if not isinstance(value, list):
            raise TypeError(error)
        for child in value:
            require_aggregate_content(child, schema[0])
    elif type(value) not in schema:
        raise ValueError(error)


def sources():
    paths = set(ROOT.glob("reports/**/*.md"))
    paths.update(ROOT.glob("docs/*.md"))
    paths.difference_update(ROOT / "docs" / name for name in (
        "EXPERIMENTS.md", "RESULTS.md", "REPRODUCING.md", "HACKATHON.md", "PUBLICATION.md"
    ))
    for folder in (
        "activation-v2", "fgfr2-sequence-pilot", "fgfr2-state-pilot",
        "fgfr2-state-bounded", "fgfr2-state-tree-pilot", "fgfr2-sequence-ap-pilot"
    ):
        paths.update((ROOT / "research" / folder).glob("*.md"))
    paths.update(ROOT.glob("reports/mechanistic-pilot/*/results.json"))
    for name in (
        "reports/improvement/20261017T194348Z-5e66cbab/final-results.json",
        "reports/revision-v2/20261017T201321Z-a29b17a8/results.json",
        "reports/model-swap-v2/20261017T210043Z-a4800896/results.json",
        "reports/research-pilot/you-pilot-20261017/results.json",
    ):
        paths.add(ROOT / name)
    paths.update(ROOT / name for name in REPEATED_JSON_FIELDS)
    return sorted(p for p in paths if p.is_file())


def main():
    selected = sources()
    included = set(selected)
    manifest = []
    for path in selected:
        original = path.read_bytes()
        relative = path.relative_to(ROOT)
        if path.suffix == ".json":
            data = json.loads(original)
            allowed = REPEATED_JSON_FIELDS.get(relative.as_posix(), RESULT_KEYS)
            public = {k: v for k, v in data.items() if k in allowed}
            if relative.as_posix() in REPEATED_JSON_FIELDS:
                require_aggregate_content(public, REPEATED_SCHEMAS[relative.as_posix()])
            text = json.dumps(clean(public), indent=2, allow_nan=False) + "\n"
            transformation = "allowlisted aggregate fields; local root removed"
        else:
            text = clean_text(original.decode())

            def link(match, source_path=path):
                label, target = match.group(1), match.group(2).strip("<>")
                if target.startswith(("https://", "http://", "#", "mailto:")):
                    return match.group(0)
                bare = target.split("#", 1)[0]
                candidate = (source_path.parent / bare).resolve()
                root_candidate = (ROOT / bare).resolve()
                if candidate in included:
                    return match.group(0)
                if root_candidate in included:
                    return f"[{label}]({os.path.relpath(root_candidate, source_path.parent)})"
                return f"{label} (`{target}`; local artifact)"

            text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link, text)
            results_link = os.path.relpath(ROOT / "docs/RESULTS.md",
                                           (DESTINATION / relative).parent)
            text = (f"> Historical research record. See [current results]({results_link}) "
                    "for the release interpretation.\n\n" + text)
            text = "\n".join(line.rstrip() for line in text.splitlines()) + "\n"
            transformation = "local root removed; private artifact links annotated; banner added"
        destination = DESTINATION / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(text)
        manifest.append({"source": relative.as_posix(), "source_sha256": digest(original),
                         "public_path": destination.relative_to(ROOT).as_posix(),
                         "public_sha256": digest(text.encode()), "transformation": transformation})
    DESTINATION.mkdir(parents=True, exist_ok=True)
    (DESTINATION / "manifest.json").write_text(json.dumps({
        "scope": "Curated aggregate-only publication; no row outcomes, traces, or credentials",
        "artifacts": manifest,
    }, indent=2) + "\n")
    print(f"Exported {len(manifest)} reviewed summary artifacts.")


if __name__ == "__main__":
    main()
