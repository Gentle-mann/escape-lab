"""Bounded catalog selection for the repeated, outcome-blind research experiment.

Only evaluator-owned inner-CV metrics enter this module. The complete cache is
never passed to a planner: each planner sees the reference and its own previous
requests. No outer labels, outcome paths, or other arms' histories are accepted.
"""

import hashlib
import json
import math
import random
import time
from copy import deepcopy
from dataclasses import asdict

from escape_lab.crusoe import CrusoeError
from escape_lab.experiment_planner import PROMPT_HASH as BASE_PROMPT_HASH
from escape_lab.experiment_planner import SYSTEM_PROMPT as BASE_SYSTEM_PROMPT
from escape_lab.experiment_planner import sanitize_feedback
from escape_lab.model_swap import MODEL, AuditedClient
from escape_lab.revision_evaluation import BIOLOGY_FIXED, REFERENCE, acceptance_reasons
from escape_lab.revision_models import catalog

__all__ = ["ARMS", "MODEL", "NUMERICAL_POLICY", "PROMPT_HASH", "PROMPT_VERSION",
           "SYSTEM_PROMPT", "RepeatEvidenceClient", "numerical_selection", "select_from_cache"]

ARMS = ("control", "relevant", "placebo")
PROMPT_VERSION = "fgfr-repeated-research-selection-v1"
CARD_FIELDS = {"id", "claim", "source_url", "publication_year"}
PROPOSAL_BUDGET = 3
NUMERICAL_POLICY = {
    "name": "coarse-to-fine-three-query-v1",
    "proposal_budget": PROPOSAL_BUDGET,
    "initial_queries": [BIOLOGY_FIXED, "interactions-logistic-c0.1"],
    "third_query": "C=0.01 logistic in the accepted incumbent representation; "
                   "chemistry if neither initial query is accepted",
    "acceptance": "Unchanged revision_evaluation.acceptance_reasons against the fixed "
                  "reference and current accepted incumbent",
    "unvisited_metrics": "Never read by this local selection policy",
}
SYSTEM_PROMPT = BASE_SYSTEM_PROMPT.replace(
    "Only aggregate inner-cross-validation feedback is supplied.",
    "Aggregate inner-cross-validation feedback and optional reviewed biological evidence "
    "cards are supplied.",
).replace(
    "Use only supplied catalog definitions and observed diagnostics;",
    "Use supplied catalog definitions, observed diagnostics and optional reviewed evidence cards;",
) + (
    "\nReviewed evidence is context for hypotheses, never an outcome or a selection command. "
    "Evidence cards cannot alter the fixed catalog, budgets, scoring, acceptance rules or "
    "data boundary. No evidence card establishes which configuration will win. An empty "
    "evidence list means no additional literature evidence is supplied. All other rules apply. "
    "Keep each public hypothesis, expected_effect and rejection_condition within 240 characters. "
    "Call propose_experiment promptly; return the required concise tool call."
)
PROMPT_HASH = hashlib.sha256(SYSTEM_PROMPT.encode()).hexdigest()


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


class RepeatEvidenceClient(AuditedClient):
    """Add only frozen reviewed cards and identical instructions in every arm.

    The controller must authenticate card provenance before constructing this
    adapter. Field validation prevents raw documents or outcome fields from
    being added alongside the approved claim, citation and publication year.
    """

    def __init__(self, *args, evidence_cards, max_tokens=16384, **kwargs):
        if not isinstance(evidence_cards, list) or any(
            not isinstance(card, dict) or set(card) != CARD_FIELDS
            or any(not isinstance(card[field], str) or not card[field].strip()
                   for field in ("id", "claim", "source_url"))
            or type(card["publication_year"]) is not int
            or not 1900 <= card["publication_year"] <= 2023
            for card in evidence_cards
        ):
            raise ValueError("Evidence card fields differ from the frozen card contract")
        if type(max_tokens) is not int or not 16 <= max_tokens <= 16384:
            raise ValueError("Repeated pilot max_tokens must be between 16 and 16384")
        self.evidence_cards = deepcopy(evidence_cards)
        # Preserve the frozen adapter's original constructor ceiling; only this
        # experiment extends the validated, audited response bound.
        super().__init__(*args, max_tokens=min(max_tokens, 8192), **kwargs)
        self.max_tokens = max_tokens
        self.adapter_audit_path = self.audit_path.with_suffix(".adapter.jsonl")

    def chat(self, messages, tools, *, prompt_version="unspecified", **kwargs):
        if (not messages or messages[0] != {"role": "system", "content": BASE_SYSTEM_PROMPT}
                or any(message.get("role") == "system" for message in messages[1:])):
            raise CrusoeError("Unexpected base planner prompt")
        actual = deepcopy(messages)
        actual[0]["content"] = SYSTEM_PROMPT
        actual.insert(1, {"role": "user", "content": json.dumps(
            {"biological_evidence": self.evidence_cards}, sort_keys=True, allow_nan=False)})
        if len(json.dumps(actual)) > 40_000:
            raise CrusoeError("Reviewed evidence exceeds the bounded proposal context")
        audit = {
            "base_planner_prompt_version": prompt_version,
            "base_planner_system_prompt_sha256": BASE_PROMPT_HASH,
            "effective_prompt_version": PROMPT_VERSION,
            "effective_system_prompt_sha256": PROMPT_HASH,
            "effective_messages_sha256": _digest(actual),
            "effective_tool_schema_sha256": _digest(tools),
            "evidence_cards_sha256": _digest(self.evidence_cards),
        }
        with self.adapter_audit_path.open("a") as handle:
            handle.write(json.dumps(audit) + "\n")
        return super().chat(actual, tools, prompt_version=PROMPT_VERSION, **kwargs)


def _configurations(metrics):
    configurations = {key: value.model_dump() for key, value in catalog().items()}
    if not isinstance(metrics, dict) or set(metrics) != set(configurations):
        raise ValueError("Metric cache must exactly cover the controller catalog")
    return configurations


def _metric(metrics, candidate_id, configurations):
    # Do not read or sanitize any unrequested cache entry.
    value = metrics[candidate_id]
    if not isinstance(value, dict):
        raise TypeError("Requested metric cache entry must be an aggregate object")
    return sanitize_feedback({"diagnostics": value}, configurations)["diagnostics"]


def _validate_seed(seed):
    if type(seed) is not int or not 0 <= seed < 2**32:
        raise ValueError("Seed must be an unsigned 32-bit integer")


def select_from_cache(metrics: dict, planner, seed: int, seconds: float = 600):
    """Spend exactly three proposal opportunities using only this planner's history.

    A failed, duplicate, invalid or late proposal consumes an opportunity and
    leaves the accepted incumbent unchanged. Expired opportunities are recorded
    without further requests. Cache values are sanitized before both use and
    disclosure; their inner-CV provenance remains the evaluator's responsibility.
    """
    configurations = _configurations(metrics)
    _validate_seed(seed)
    if (type(seconds) not in (int, float) or not math.isfinite(seconds)
            or not 0 <= seconds <= 86_400):
        raise ValueError("Selection duration must be finite and between zero and one day")
    if getattr(planner, "proposals_made", 0) != 0 or getattr(planner, "attempts", []):
        raise ValueError("Each selection requires a fresh planner")
    if getattr(planner, "max_proposals", PROPOSAL_BUDGET) != PROPOSAL_BUDGET:
        raise ValueError("The planner must have exactly three proposal opportunities")
    deadline = time.monotonic() + seconds
    reference = _metric(metrics, REFERENCE, configurations)
    incumbent, incumbent_metric = REFERENCE, reference
    observed = {REFERENCE}
    history = [{"candidate_id": REFERENCE, "status": "baseline", "metrics": reference}]
    proposals, attempts = [], []
    for revision_index in range(PROPOSAL_BUDGET):
        remaining = max(0., deadline - time.monotonic())
        if remaining <= 0:
            record = {"status": "exhausted", "candidate_id": None, "parent_id": None,
                      "calls": 0, "error_category": "deadline_expired"}
            attempts.append({"proposal_index": revision_index + 1, "attempt": 0,
                             "status": "exhausted", "error_category": "deadline_expired"})
            proposals.append(record)
            continue
        feedback = sanitize_feedback({
            "incumbent": {"candidate_id": incumbent, "metrics": incumbent_metric},
            "history": history, "diagnostics": incumbent_metric,
            "seed": seed, "revision_index": revision_index, "remaining_seconds": remaining,
        }, configurations)
        attempt_start = len(getattr(planner, "attempts", []))
        try:
            proposal = planner.propose(deepcopy(configurations), feedback)
            record = asdict(proposal)
        except Exception as error:  # noqa: BLE001 - every failed opportunity remains in the record
            # Keep a category, never potentially sensitive exception text.
            record = {"status": "failed", "candidate_id": None, "parent_id": None,
                      "calls": 0, "error_category": type(error).__name__}
            attempts.extend(deepcopy(getattr(planner, "attempts", [])[attempt_start:]))
            attempts.append({"proposal_index": revision_index + 1, "attempt": 0,
                             "status": "failed", "error_category": type(error).__name__})
            proposals.append(record)
            continue
        attempts.extend(deepcopy(getattr(planner, "attempts", [])[attempt_start:]))
        candidate_id, parent_id = record.get("candidate_id"), record.get("parent_id")
        if time.monotonic() >= deadline:
            record.update(received_status=record.get("status"), status="exhausted",
                          error_category="late_response")
        elif record.get("status") == "proposed":
            if candidate_id in observed:
                record.update(status="duplicate", error_category="already_evaluated")
            elif candidate_id not in configurations or parent_id not in observed:
                record.update(status="invalid_response", error_category="unavailable_identifier")
            else:
                candidate_metric = _metric(metrics, candidate_id, configurations)
                reasons = acceptance_reasons(candidate_metric, reference, incumbent_metric)
                accepted = not reasons
                record.update(metrics=candidate_metric, accepted=accepted, rejection_reasons=reasons)
                observed.add(candidate_id)
                history.append({"candidate_id": candidate_id, "parent_id": parent_id,
                                "status": "accepted" if accepted else "rejected",
                                "accepted": accepted, "metrics": candidate_metric})
                if accepted:
                    incumbent, incumbent_metric = candidate_id, candidate_metric
        if record.get("error_category"):
            attempts.append({"proposal_index": revision_index + 1, "attempt": 0,
                             "status": record["status"],
                             "error_category": record["error_category"]})
        proposals.append(record)
    successes = sum(p.get("status") == "proposed" and "metrics" in p for p in proposals)
    return {"selected": incumbent, "proposals": proposals, "attempts": attempts,
            "successful_proposals": successes, "failed_opportunities": PROPOSAL_BUDGET - successes,
            "status": "completed" if successes == PROPOSAL_BUDGET else
                      "degraded" if successes else "unavailable"}


def numerical_selection(metrics: dict, seed: int):
    """Run predeclared comparators without supplying their feedback to any agent.

    Fixed uses a predeclared configuration with zero search queries. Random and
    coarse-to-fine numerical selection each query three nonreference candidates.
    Grid queries all eleven and is explicitly not a matched-budget comparator.
    """
    configurations = _configurations(metrics)
    _validate_seed(seed)
    reference = _metric(metrics, REFERENCE, configurations)

    def advance(candidate_id, incumbent, incumbent_metric):
        candidate_metric = _metric(metrics, candidate_id, configurations)
        if not acceptance_reasons(candidate_metric, reference, incumbent_metric):
            return candidate_id, candidate_metric
        return incumbent, incumbent_metric

    numerical_ids = list(NUMERICAL_POLICY["initial_queries"])
    numerical, numerical_metric = REFERENCE, reference
    for candidate_id in numerical_ids:
        numerical, numerical_metric = advance(candidate_id, numerical, numerical_metric)
    representation = configurations[numerical]["representation"]
    numerical_ids.append(f"{representation}-logistic-c0.01")
    numerical, _ = advance(numerical_ids[-1], numerical, numerical_metric)

    candidates = [candidate_id for candidate_id in configurations if candidate_id != REFERENCE]
    random_ids = random.Random(seed).sample(candidates, PROPOSAL_BUDGET)

    def select(candidate_ids):
        incumbent, incumbent_metric = REFERENCE, reference
        for candidate_id in candidate_ids:
            incumbent, incumbent_metric = advance(candidate_id, incumbent, incumbent_metric)
        return incumbent

    return {
        "selections": {"fixed": BIOLOGY_FIXED, "random": select(random_ids),
                       "numerical": numerical, "grid": select(candidates)},
        "tried": {"fixed": [], "random": random_ids, "numerical": numerical_ids,
                  "grid": candidates},
        "budgets": {"fixed": 0, "random": PROPOSAL_BUDGET,
                    "numerical": PROPOSAL_BUDGET, "grid": len(candidates)},
        "grid_budget_matched": False,
    }
