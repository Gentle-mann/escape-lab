"""Trusted grouped replay evaluation; agents only receive authorized observations.

This module is an evaluator, not an agent tool. It may read hidden measurements
to score predictions, but never supplies those scores or labels to investigators.
The final split is deliberately unsupported until a separate freeze procedure exists.
"""

import hashlib
import json
import random
import statistics
import time
from collections import Counter
from copy import copy
from datetime import UTC, datetime
from operator import itemgetter
from pathlib import Path
from uuid import uuid4

import numpy as np
from sklearn.metrics import average_precision_score, brier_score_loss

from escape_lab.agent import PROMPT_VERSION, SYSTEM_PROMPT
from escape_lab.baseline import NumericalSearch
from escape_lab.crusoe import CrusoeError
from escape_lab.data import sha256
from escape_lab.vault import Vault

POLICIES = ("random", "gene_prior", "logistic", "extra_trees", "extra_explore", "extra_diverse")
ALLOWED_POLICIES = {*POLICIES, "agent"}
FOLD_SALT = "escape-lab-development-folds-v1"
REFIT_INTERVAL = 5


def _group_hash(group):
    return hashlib.sha256(f"{FOLD_SALT}:{group}".encode()).hexdigest()


def validate_partition(training, pool):
    train_ids = {candidate["candidate_id"] for candidate in training}
    pool_ids = {candidate["candidate_id"] for candidate in pool}
    train_groups = {candidate["group_id"] for candidate in training}
    pool_groups = {candidate["group_id"] for candidate in pool}
    if train_ids & pool_ids or train_groups & pool_groups:
        raise ValueError("Training and query pool must have disjoint biological groups")
    if len(train_ids) != len(training) or len(pool_ids) != len(pool):
        raise ValueError("Duplicate protein candidate identity in a partition")


def build_tasks(candidates, stage):
    """Build folds from public group identities, before inspecting any outcomes."""
    if stage == "final":
        raise ValueError("Final test remains sealed")
    if stage not in {"development", "validation"}:
        raise ValueError("Unknown evaluation stage")
    group_splits = {}
    for candidate in candidates:
        group, split = candidate["group_id"], candidate["split"]
        if group in group_splits and group_splits[group] != split:
            raise ValueError("A biological group crosses the declared split manifest")
        group_splits[group] = split
    development = [candidate for candidate in candidates if candidate["split"] == "development"]
    if stage == "validation":
        tasks = [
            {
                "task_id": "validation-all-groups",
                "training": development,
                "pool": [candidate for candidate in candidates if candidate["split"] == stage],
            }
        ]
    else:
        groups = sorted({candidate["group_id"] for candidate in development}, key=_group_hash)
        if len(groups) < 3:
            raise ValueError("Development evaluation requires at least three biological groups")
        folds = {group: index % 3 for index, group in enumerate(groups)}
        tasks = []
        for fold in range(3):
            tasks.append(
                {
                    "task_id": f"development-fold-{fold}",
                    "training": [
                        candidate
                        for candidate in development
                        if folds[candidate["group_id"]] != fold
                    ],
                    "pool": [
                        candidate
                        for candidate in development
                        if folds[candidate["group_id"]] == fold
                    ],
                }
            )
    for task in tasks:
        validate_partition(task["training"], task["pool"])
        if not task["training"] or not task["pool"]:
            raise ValueError("Every task requires nonempty training and query pools")
    return tasks


def _prediction(learner, candidates, policy):
    if policy == "random":
        return np.full(len(candidates), learner.prior)
    if policy == "gene_prior":
        return np.asarray(
            [learner.gene_rates.get(candidate["gene"], learner.prior) for candidate in candidates]
        )
    return learner.predict(candidates)[0]


def _select(policy, shortlists, rng, explored_groups):
    if policy == "random":
        return rng.choice(shortlists["exploit"])
    if policy == "gene_prior":
        rows, key = shortlists["gene_prior"], "gene_prior"
    elif policy == "extra_explore":
        rows, key = shortlists["explore"], "exploration_score"
    else:
        rows, key = shortlists["exploit"], "numerical_probability"
    if policy == "extra_diverse":
        unexplored = [row for row in rows if row["group_id"] not in explored_groups]
        if unexplored:
            rows = unexplored
    best = max(row[key] for row in rows)
    return rng.choice([row for row in rows if row[key] == best])


def _agent_observations(observations, candidates_by_id):
    """Flatten only revealed data; one compact observation per protein variant."""
    result = []
    for observation in observations:
        outcome = observation["outcome"]
        candidate = candidates_by_id[observation["candidate_id"]]
        measurements = []
        for record in outcome["records"]:
            item = {"nucleotide_variant": record.get("nucleotide_variant")}
            for drug in ("pem", "fut"):
                values = record.get(drug, {})
                item[f"{drug}_replicates"] = values.get("replicates")
                item[f"{drug}_median"] = values.get("median")
                item[f"{drug}_source_sha256"] = values.get("source_sha256")
            measurements.append(item)
        result.append(
            {
                "candidate_id": candidate["candidate_id"],
                "gene": candidate["gene"],
                "mutation": candidate["protein_variant"],
                "evidence_id": outcome["evidence_ids"][0] if outcome["evidence_ids"] else None,
                "dual_resistant": outcome["dual_resistant"],
                "measurements": measurements,
            }
        )
    return result


def _usage(investigator):
    client = getattr(investigator, "client", None)
    return int(getattr(client, "total_tokens", 0)), len(getattr(client, "usage_records", []))


def _trial(vault, task, outcomes, policy, seed, budget, stage, investigator, label_rule):
    start = time.monotonic()
    training, pool = task["training"], task["pool"]
    training = [
        candidate
        for candidate in training
        if type(outcomes[candidate["candidate_id"]]["dual_resistant"]) is bool
    ]
    validate_partition(training, pool)
    if not pool:
        raise ValueError("No query candidates remain")
    training = sorted(training, key=itemgetter("candidate_id"))
    pool = sorted(pool, key=itemgetter("candidate_id"))
    pool_by_id = {candidate["candidate_id"]: candidate for candidate in pool}
    training_labels = [
        outcomes[candidate["candidate_id"]]["dual_resistant"] for candidate in training
    ]
    kind = "logistic" if policy in {"random", "gene_prior", "logistic"} else "extra_trees"
    learner = NumericalSearch(kind=kind, seed=seed, trees=100, leaf=3)
    learner.fit(training, training_labels)
    # These scores stay in this trusted evaluator, never in agent state or tools.
    audit_pool = [
        candidate
        for candidate in pool
        if type(outcomes[candidate["candidate_id"]]["dual_resistant"]) is bool
    ]
    initial_prediction = _prediction(learner, audit_pool, policy)
    audit_labels = [
        outcomes[candidate["candidate_id"]]["dual_resistant"] for candidate in audit_pool
    ]
    initial_ap, initial_brier = None, None
    if audit_labels:
        initial_ap = (
            float(average_precision_score(audit_labels, initial_prediction))
            if any(audit_labels)
            else 0.0
        )
        initial_brier = float(brier_score_loss(audit_labels, initial_prediction))
    effective_budget = min(budget, len(pool))
    policy_version = f"{policy}-v1-trees100-leaf3-refit{REFIT_INTERVAL}"
    if policy == "agent":
        policy_version += f"-{PROMPT_VERSION}"
    run = vault.create_run(list(pool_by_id), effective_budget, policy_version, stage, seed)
    auth = (run["run_id"], run["run_token"])
    rng = random.Random(seed)
    selected, observations, trace, decision_attempts = [], [], [], []
    explored_groups, positive_groups = set(), set()
    gene_hits, gene_queries = Counter(), Counter()
    fallback_count, failure = 0, None
    live_investigator = copy(investigator) if policy == "agent" else None
    if live_investigator is not None and hasattr(live_investigator, "decisions_made"):
        live_investigator.decisions_made = 0
    initial_tokens, initial_usage_count = _usage(live_investigator)
    first_hit = None
    for step in range(effective_budget):
        shortlist_limit = 8 if policy == "agent" else len(pool)
        shortlists = learner.shortlists(
            pool, set(selected), limit=shortlist_limit, explored_groups=explored_groups
        )
        if policy == "agent":
            state = {
                "task_id": task["task_id"],
                "objective": "coverage_audit",
                "step": step,
                "budget_remaining": effective_budget - step,
                "revealed_ids": selected.copy(),
                "observations": _agent_observations(observations, pool_by_id),
                "contract": {
                    "objective": "coverage_audit",
                    "budget": effective_budget,
                    "query_unit": "One unique protein variant's paired archived drug measurements",
                    "assay": "Retrospective separate-drug FGFR cell-line resistance",
                    "resistance_criterion": label_rule,
                    "allowed_conclusions": ["measured counterexample", "unresolved coverage"],
                },
            }
            try:
                decision = live_investigator.choose_action(state, shortlists)
            except (CrusoeError, ValueError, TypeError, KeyError, RuntimeError) as error:
                failure = f"Agent decision failed ({type(error).__name__}); no silent replacement."
                decision_attempts.append(
                    {"step": step + 1, "action": "failed", "failure_type": type(error).__name__}
                )
                break
            decision_attempts.append(
                {
                    "step": step + 1,
                    "action": decision.action,
                    "candidate_id": getattr(decision, "candidate_id", None),
                    "reason": getattr(decision, "reason", None),
                    "fallback_used": bool(getattr(decision, "fallback_used", False)),
                    "fallback_reason": getattr(decision, "fallback_reason", None),
                    "prompt_version": getattr(decision, "prompt_version", PROMPT_VERSION),
                }
            )
            fallback_count += int(getattr(decision, "fallback_used", False))
            if decision.action != "select_candidate":
                failure = "Agent stopped before the coverage audit budget was exhausted."
                break
            permitted = {row["candidate_id"]: row for rows in shortlists.values() for row in rows}
            if decision.candidate_id not in permitted or decision.candidate_id in selected:
                failure = "Agent selected an ineligible or already measured candidate."
                break
            chosen = permitted[decision.candidate_id]
            reason = decision.reason
            used_fallback = bool(getattr(decision, "fallback_used", False))
            fallback_reason = getattr(decision, "fallback_reason", None)
            prompt_version = getattr(decision, "prompt_version", None)
        else:
            chosen = _select(policy, shortlists, rng, explored_groups)
            reason = f"Predefined {policy} acquisition with randomized exact-score ties."
            used_fallback = False
            fallback_reason, prompt_version = None, None
        candidate_id = chosen["candidate_id"]
        probability = float(_prediction(learner, [pool_by_id[candidate_id]], policy)[0])
        receipt = vault.commit(
            *auth, candidate_id, probability, reason, policy_version, f"commit-{step}"
        )
        measurement = vault.reveal(*auth, receipt["commit_id"], f"reveal-{step}")
        label = measurement["outcome"]["dual_resistant"]
        candidate = pool_by_id[candidate_id]
        selected.append(candidate_id)
        observations.append(measurement)
        explored_groups.add(candidate["group_id"])
        gene_queries[candidate["gene"]] += 1
        if label:
            positive_groups.add(candidate["group_id"])
            gene_hits[candidate["gene"]] += 1
            if first_hit is None:
                first_hit = step + 1
        trace.append(
            {
                "step": step + 1,
                "prediction_before_reveal": probability,
                "decision_reason": reason,
                "fallback_used": used_fallback,
                "fallback_reason": fallback_reason,
                "prompt_version": prompt_version,
                "commit_id": receipt["commit_id"],
                "committed_at": receipt["committed_at"],
                **measurement,
            }
        )
        if (step + 1) % REFIT_INTERVAL == 0 and step + 1 < effective_budget:
            disclosed = [
                item for item in observations if type(item["outcome"]["dual_resistant"]) is bool
            ]
            learner.fit(
                training + [pool_by_id[item["candidate_id"]] for item in disclosed],
                training_labels + [item["outcome"]["dual_resistant"] for item in disclosed],
            )
    state = vault.state(*auth)
    hits = sum(item["outcome"]["dual_resistant"] is True for item in observations)
    evidence_ids = sorted(
        {
            evidence
            for item in observations
            if item["outcome"]["dual_resistant"]
            for evidence in item["outcome"]["evidence_ids"]
        }
    )
    verdict = vault.verdict(*auth, "refuted" if hits else "unresolved", evidence_ids)
    final_tokens, final_usage_count = _usage(live_investigator)
    usage_fields = {
        "model_id",
        "prompt_version",
        "prompt_tokens",
        "completion_tokens",
        "total_tokens",
        "budget_charge_tokens",
        "usage_estimated",
    }
    usage_records = getattr(getattr(live_investigator, "client", None), "usage_records", [])
    usage = [
        {key: value for key, value in record.items() if key in usage_fields}
        for record in usage_records[initial_usage_count:final_usage_count]
    ]
    valid = failure is None and len(selected) == effective_budget
    return {
        "task_id": task["task_id"],
        "policy": policy,
        "policy_version": policy_version,
        "seed": seed,
        "run_id": run["run_id"],
        "requested_budget": budget,
        "query_budget": effective_budget,
        "queries": len(selected),
        "hits": hits,
        "precision": hits / len(selected) if selected else 0.0,
        "delivered_hits_per_budget": hits / effective_budget,
        "query_completion_fraction": len(selected) / effective_budget,
        "inconclusive_reveals": sum(
            item["outcome"]["dual_resistant"] is None for item in observations
        ),
        "first_hit_query": first_hit,
        "positive_region_coverage": len(positive_groups),
        "queried_region_coverage": len(explored_groups),
        "gene_hits": dict(gene_hits),
        "gene_queries": dict(gene_queries),
        "initial_average_precision": initial_ap,
        "initial_brier": initial_brier,
        "initial_audit_count": len(audit_pool),
        "pool_count": len(pool),
        "initial_audit_undefined_labels": len(pool) - len(audit_pool),
        "initial_training_count": len(training),
        "training_groups": sorted({candidate["group_id"] for candidate in training}),
        "pool_groups": sorted({candidate["group_id"] for candidate in pool}),
        "excluded_inconclusive_training": len(task["training"]) - len(training),
        "fit_count": learner.fit_count,
        "selected_candidate_ids": selected,
        "committed_count": state["committed_count"],
        "revealed_count": state["revealed_count"],
        "fallback_count": fallback_count,
        "agent_unassisted": policy == "agent" and fallback_count == 0 and valid,
        "execution_mode": (
            "numerical"
            if policy != "agent"
            else "agent_with_fallback"
            if fallback_count
            else "agent_without_fallback"
        ),
        "valid": valid,
        "failure": failure,
        "verdict": verdict,
        "trace": trace,
        "decision_attempts": decision_attempts,
        "elapsed_seconds": time.monotonic() - start,
        "model_id": getattr(getattr(live_investigator, "client", None), "model_id", None),
        "model_budget_charge_tokens": final_tokens - initial_tokens,
        "model_completion_calls": final_usage_count - initial_usage_count,
        "model_usage": usage,
    }


def _summarize(trials, policies):
    summary = {}
    for policy in policies:
        attempted = [trial for trial in trials if trial["policy"] == policy]
        valid = [trial for trial in attempted if trial["valid"]]
        unassisted = [
            trial for trial in attempted if policy == "agent" and not trial["fallback_count"]
        ]
        assisted = [trial for trial in attempted if trial["fallback_count"]]
        row = {
            "attempted_trials": len(attempted),
            "valid_trials": len(valid),
            "completion_rate": len(valid) / len(attempted),
            "fallback_assisted_trials": len(assisted),
            "unassisted_agent_attempts": len(unassisted),
            "unassisted_agent_completed": sum(trial["valid"] for trial in unassisted),
            "unassisted_mean_hits": statistics.mean(trial["hits"] for trial in unassisted)
            if unassisted
            else None,
            "fallback_assisted_mean_hits": statistics.mean(trial["hits"] for trial in assisted)
            if assisted
            else None,
            "completed_only_mean_hits": statistics.mean(trial["hits"] for trial in valid)
            if valid
            else None,
            "primary_mean_scope": "All attempted trials, including partial failures and marked fallback assistance",
        }
        for metric in (
            "hits",
            "precision",
            "positive_region_coverage",
            "initial_average_precision",
            "initial_brier",
            "elapsed_seconds",
            "model_budget_charge_tokens",
            "delivered_hits_per_budget",
            "query_completion_fraction",
            "inconclusive_reveals",
        ):
            values = [trial[metric] for trial in attempted if trial[metric] is not None]
            row[f"mean_{metric}"] = statistics.mean(values) if values else None
        summary[policy] = row
    paired = []
    for policy_index, left in enumerate(policies):
        for right in policies[policy_index + 1 :]:
            left_trials = {
                (trial["task_id"], trial["seed"]): trial
                for trial in trials
                if trial["policy"] == left
            }
            right_trials = {
                (trial["task_id"], trial["seed"]): trial
                for trial in trials
                if trial["policy"] == right
            }
            keys = sorted(left_trials.keys() & right_trials.keys())
            if not keys:
                continue
            differences = [right_trials[key]["hits"] - left_trials[key]["hits"] for key in keys]
            paired.append(
                {
                    "left": left,
                    "right": right,
                    "paired_trials": len(keys),
                    "mean_hits_difference_right_minus_left": statistics.mean(differences),
                    "differences": differences,
                    "scope": "Delivered system utility across all attempted matched trials; partial failures retained",
                    "both_completed_pairs": sum(
                        left_trials[key]["valid"] and right_trials[key]["valid"] for key in keys
                    ),
                    "includes_agent_fallback": any(
                        left_trials[key]["fallback_count"] or right_trials[key]["fallback_count"]
                        for key in keys
                    ),
                    "no_fallback_pairs": sum(
                        not left_trials[key]["fallback_count"]
                        and not right_trials[key]["fallback_count"]
                        for key in keys
                    ),
                    "mean_hits_difference_without_fallback": statistics.mean(
                        right_trials[key]["hits"] - left_trials[key]["hits"]
                        for key in keys
                        if not left_trials[key]["fallback_count"]
                        and not right_trials[key]["fallback_count"]
                    )
                    if any(
                        not left_trials[key]["fallback_count"]
                        and not right_trials[key]["fallback_count"]
                        for key in keys
                    )
                    else None,
                }
            )
    return summary, paired


def _markdown(report):
    lines = [
        f"# Escape Lab {report['stage']} benchmark",
        "",
        "Retrospective cell-line replay. Final evaluation remains sealed.",
        "",
        f"Requested budget: {report['configuration']['budget']} paired-record reveals per trial.",
        "All methods receive the same disclosed initial training data; those observations are not charged to this query budget.",
        "",
        "Mean utility includes every attempted trial, including partial failures. Fallback-assisted agent utility is hybrid system performance.",
        "| Policy | Completed trials | Mean delivered hits | Precision | Initial AP | Initial Brier | Fallback trials |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for policy, row in report["summary"].items():
        values = []
        for metric in (
            "mean_hits",
            "mean_precision",
            "mean_initial_average_precision",
            "mean_initial_brier",
        ):
            values.append("n/a" if row[metric] is None else f"{row[metric]:.3f}")
        lines.append(
            f"| {policy} | {row['valid_trials']}/{row['attempted_trials']} | "
            + " | ".join(values)
            + f" | {row['fallback_assisted_trials']} |"
        )
    lines += [
        "",
        "## Interpretation limits",
        "",
        "- Hits count unique measured protein variants, not independent mechanisms or patients.",
        "- Initial AP and Brier use the fixed blinded query pool before any result is revealed.",
        "- Inconclusive variants remain queryable and consume budget; undefined labels are omitted only from fitting and prediction scoring.",
        "- Repeated seeds measure computational variability, not additional biological replication.",
        "- The exploratory acquisition adds tree disagreement; it is not nonmyopic ENS.",
        "- Fallback-assisted agent runs are explicitly marked and do not establish unassisted Nemotron performance.",
        "- Paired differences are descriptive; no patient-level confidence interval is claimed.",
        "",
    ]
    return "\n".join(lines)


def run_benchmark(
    data_dir,
    output_dir,
    stage="development",
    policies=POLICIES,
    seeds=(0, 1, 2),
    budget=20,
    investigator=None,
):
    if stage == "final":
        raise ValueError("Final test remains sealed")
    if stage not in {"development", "validation"}:
        raise ValueError("Unknown evaluation stage")
    policies, seeds = tuple(policies), tuple(seeds)
    if not policies or len(set(policies)) != len(policies) or not set(policies) <= ALLOWED_POLICIES:
        raise ValueError("Choose unique supported policies")
    if "agent" in policies and investigator is None:
        raise ValueError("The agent policy requires a configured investigator")
    if (
        not seeds
        or len(set(seeds)) != len(seeds)
        or any(type(seed) is not int or seed < 0 for seed in seeds)
    ):
        raise ValueError("Choose unique nonnegative integer seeds")
    if type(budget) is not int or budget < 1:
        raise ValueError("Budget must be a positive integer")
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    code_paths = [
        Path(__file__).with_name(filename)
        for filename in (
            "benchmark.py",
            "baseline.py",
            "vault.py",
            "data.py",
            "agent.py",
            "crusoe.py",
        )
    ]
    source_code_fingerprints = {path.name: sha256(path) for path in code_paths}
    client = getattr(investigator, "client", None)
    agent_configuration = {
        key: getattr(client, key, None)
        for key in ("model_id", "max_tokens", "max_total_tokens", "max_retries")
    }
    agent_configuration.update(
        {
            "prompt_version": PROMPT_VERSION,
            "system_prompt_sha256": hashlib.sha256(SYSTEM_PROMPT.encode()).hexdigest(),
            "temperature": "Provider default; not explicitly set",
        }
    )
    manifest = json.loads((data_dir / "manifest.json").read_text())
    for filename in ("candidates", "outcomes"):
        if sha256(data_dir / f"{filename}.json") != manifest[f"{filename}_sha256"]:
            raise ValueError(f"Prepared {filename} checksum mismatch")
    candidates = json.loads((data_dir / "candidates.json").read_text())
    outcomes = json.loads((data_dir / "outcomes.json").read_text())
    tasks = build_tasks(candidates, stage)
    report_id = str(uuid4())
    private = output_dir / "private" / "vaults"
    private.mkdir(parents=True, exist_ok=True)
    private.chmod(0o700)
    vault = Vault(private / f"{stage}-{report_id}.sqlite", candidates, outcomes)
    trials = []
    for task in tasks:
        for seed in seeds:
            for policy in policies:
                trials.append(
                    _trial(
                        vault,
                        task,
                        outcomes,
                        policy,
                        seed,
                        budget,
                        stage,
                        investigator,
                        manifest["label_rule"],
                    )
                )
    summary, paired = _summarize(trials, policies)
    report = {
        "report_id": report_id,
        "stage": stage,
        "created_at": datetime.now(UTC).isoformat(),
        "final_test_status": "sealed",
        "dataset_fingerprints": {
            key: manifest[key] for key in ("candidates_sha256", "outcomes_sha256")
        },
        "source_code_fingerprints": source_code_fingerprints,
        "source_changed_during_run": source_code_fingerprints
        != {path.name: sha256(path) for path in code_paths},
        "configuration": {
            "policies": list(policies),
            "seeds": list(seeds),
            "budget": budget,
            "refit_interval": REFIT_INTERVAL,
            "trees": 100,
            "min_samples_leaf": 3,
            "development_fold_salt": FOLD_SALT,
            "label_rule": manifest["label_rule"],
            "agent_generation": agent_configuration,
            "seed_controls": "Numerical model and tie-breaking; provider decoding seed is not set",
            "inconclusive_policy": "Remain queryable and consume budget; omit undefined labels only from fitting and prediction scoring",
            "initial_training": "Identical disclosed training set per task; not charged as query reveals",
        },
        "summary": summary,
        "paired_comparisons": paired,
        "trials": trials,
        "interpretation": "Descriptive grouped retrospective replay; seeds are not independent biological samples.",
    }
    encoded = json.dumps(report, indent=2, allow_nan=False) + "\n"
    history = output_dir / "history"
    history.mkdir(parents=True, exist_ok=True)
    (history / f"{stage}-{report_id}.json").write_text(encoded)
    (output_dir / f"{stage}.json").write_text(encoded)
    (output_dir / f"{stage}.md").write_text(_markdown(report))
    return report
