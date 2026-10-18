"""Runnable synthetic policy transfer benchmark with a frozen independent referee.

This tests persistent diagnostic policies on generated numerical failure families.
The default provider is scripted;
--provider crusoe explicitly enables billed calls on generated training cases only.
"""

import argparse
import json
import os
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict, Field
from scipy.special import expit

from escape_lab.activation_planner import ActivationDiagnostics, ActivationPlanner
from escape_lab.crusoe import CrusoeClient
from escape_lab.research_policy import (
    ResearchPolicy,
    content_hash,
    read_policy,
    write_content_addressed,
    write_policy,
)

BENCHMARK_VERSION = 'synthetic-diagnostic-transfer-v1'


class FrozenReferee(BaseModel):
    model_config = ConfigDict(frozen=True, strict=True, extra='forbid')
    version: Literal['synthetic-referee-v1'] = 'synthetic-referee-v1'
    minimum_heldout_cases: Literal[24] = 24
    minimum_success_rate_gain: Literal[0.1] = .1
    maximum_false_positive_rate_increase: Literal[0.0] = 0.
    case_metric: Literal['exact_injected_failure_diagnosis'] = 'exact_injected_failure_diagnosis'


REFEREE = FrozenReferee()


class PromotionReceipt(BaseModel):
    model_config = ConfigDict(frozen=True, strict=True, extra='forbid')
    parent_hash: str = Field(pattern=r'^[0-9a-f]{64}$')
    challenger_hash: str = Field(pattern=r'^[0-9a-f]{64}$')
    referee_hash: str = Field(pattern=r'^[0-9a-f]{64}$')
    suite_hash: str = Field(pattern=r'^[0-9a-f]{64}$')
    scope: Literal['synthetic_mechanism_test'] = 'synthetic_mechanism_test'
    heldout_cases: int = Field(ge=0)
    parent_success_rate: float = Field(ge=0, le=1)
    challenger_success_rate: float = Field(ge=0, le=1)
    parent_false_positive_rate: float = Field(ge=0, le=1)
    challenger_false_positive_rate: float = Field(ge=0, le=1)
    accepted: bool


@dataclass(frozen=True)
class SyntheticCase:
    seed: int
    kind: str
    probabilities: np.ndarray
    ranking_scores: np.ndarray
    labels: np.ndarray
    scaled_features: np.ndarray
    variances: np.ndarray
    expected: frozenset

    @property
    def digest(self):
        return content_hash({'seed': self.seed, 'kind': self.kind,
                             'probabilities': self.probabilities.tolist(),
                             'ranking_scores': self.ranking_scores.tolist(),
                             'labels': self.labels.tolist(),
                             'scaled_features': self.scaled_features.tolist(),
                             'variances': self.variances.tolist(),
                             'expected': sorted(self.expected)})

    def executor(self):
        return ActivationDiagnostics(scope='synthetic', probabilities=self.probabilities,
            ranking_scores=self.ranking_scores, labels=self.labels,
            scaled_features=self.scaled_features, training_variances=self.variances)


def generate_case(seed, kind):
    if kind not in ('clean', 'saturation', 'feature_shift', 'combined'):
        raise ValueError('Unknown synthetic failure family')
    rng = np.random.default_rng(seed)
    n = int(rng.integers(40, 120))
    score = rng.normal(size=n)
    labels = score > np.quantile(score, .85)
    features = rng.uniform(-2., 2., (n, 4))
    variances = rng.uniform(.4, 1.4, 4)
    expected = set()
    if kind in ('saturation', 'combined'):
        # Every probability rounds to one, but order is still present in the raw score.
        score = score + float(rng.uniform(100., 400.))
        expected.add('saturation')
    if kind in ('feature_shift', 'combined'):
        features[:max(2, n // 8), 0] = rng.uniform(20., 100., max(2, n // 8))
        variances[0] = 1e-16
        expected.add('feature_shift')
    return SyntheticCase(seed, kind, expit(score), score, labels, features, variances,
                         frozenset(expected))


class ScriptedDiagnosticClient:
    """Deterministic protocol fixture. Never presented as a language model run."""

    model_id = 'scripted/offline-diagnostic-fixture'
    total_tokens = 0

    def chat(self, messages, tools, *, prompt_version):
        context = json.loads(messages[1]['content'])
        used = {record['action'] for record in context['history']}
        arguments = {}
        if 'inspect_score_saturation' not in used:
            action = 'inspect_score_saturation'
        elif 'compare_rankings' not in used:
            action = 'compare_rankings'
        elif 'inspect_feature_shift' not in used:
            action = 'inspect_feature_shift'
        elif context['earned_policy_rules'] and 'propose_policy_revision' not in used:
            action = 'propose_policy_revision'
            arguments = {'rules': context['earned_policy_rules']}
        else:
            action, arguments = 'stop', {'reason': 'diagnostics_complete'}
        if action not in ('stop', 'propose_policy_revision'):
            arguments = {'target': 'bound_partition'}
        return SimpleNamespace(finish_reason='tool_calls', message={'tool_calls': [{
            'id': 'synthetic-call', 'function': {'name': action,
                                               'arguments': json.dumps(arguments)}}]})


def execute_policy(policy, case):
    """Deterministic interpreter for the same typed rules the LLM receives.

    The frozen benchmark evaluates retained policy behavior independently of model
    sampling. It never calls an LLM or changes rules while seeing heldout outcomes.
    """
    executor = case.executor()
    findings = set()
    actions = []
    if 'inspect_saturation_before_ranking' in policy.rules:
        diagnostic = executor.execute('inspect_score_saturation')
        actions.append('inspect_score_saturation')
        if diagnostic.get('probability_saturated_fraction', 0) > 0:
            comparison = executor.execute('compare_rankings')
            actions.append('compare_rankings')
            if comparison.get('delta_ap', 0) > 0:
                findings.add('saturation')
    if 'inspect_feature_shift_before_promotion' in policy.rules:
        diagnostic = executor.execute('inspect_feature_shift')
        actions.append('inspect_feature_shift')
        if diagnostic.get('feature_shift_fraction', 0) > 0:
            findings.add('feature_shift')
    actions.append('stop')
    return {'findings': frozenset(findings), 'actions': actions}


def score_policy(policy, cases):
    correct = false_positives = 0
    for case in cases:
        result = execute_policy(policy, case)
        correct += result['findings'] == case.expected
        false_positives += bool(result['findings'] - case.expected)
    return {'cases': len(cases), 'correct': correct,
            'success_rate': float(correct / len(cases)),
            'false_positive_rate': float(false_positives / len(cases))}


def judge_policies(parent, challenger, cases):
    if not cases or len({case.digest for case in cases}) != len(cases):
        raise ValueError('Referee requires distinct nonempty synthetic cases')
    if challenger.parent_hash != parent.digest or challenger.version != parent.version + 1:
        raise ValueError('Challenger lineage does not match the frozen parent')
    parent_score, challenger_score = score_policy(parent, cases), score_policy(challenger, cases)
    accepted = (
        len(cases) >= REFEREE.minimum_heldout_cases and
        challenger_score['success_rate'] - parent_score['success_rate'] >=
        REFEREE.minimum_success_rate_gain and
        challenger_score['false_positive_rate'] - parent_score['false_positive_rate'] <=
        REFEREE.maximum_false_positive_rate_increase
    )
    receipt = PromotionReceipt(parent_hash=parent.digest, challenger_hash=challenger.digest,
        referee_hash=content_hash(REFEREE.model_dump()),
        suite_hash=content_hash([case.digest for case in cases]), heldout_cases=len(cases),
        parent_success_rate=parent_score['success_rate'],
        challenger_success_rate=challenger_score['success_rate'],
        parent_false_positive_rate=parent_score['false_positive_rate'],
        challenger_false_positive_rate=challenger_score['false_positive_rate'], accepted=accepted)
    return {'parent': parent_score, 'challenger': challenger_score}, receipt


def run_policy_smoke(output_dir, *, client=None, max_actions=6):
    output_dir = Path(output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError('Choose a new output directory; previous attempts are immutable')
    output_dir.mkdir(parents=True, exist_ok=True)
    provider = 'crusoe' if client is not None else 'scripted_offline'
    actual_client = client if client is not None else ScriptedDiagnosticClient()
    parent = ResearchPolicy.initial()
    write_policy(output_dir / 'policies', parent)
    # Training cases are made first. Only these cases produce a rule proposal.
    training = [generate_case(111, 'combined'), generate_case(112, 'clean')]
    training_run = ActivationPlanner(actual_client, max_actions=max_actions,
        proposal_source="llm_proposed" if client is not None else "scripted").run(
        training[0].executor(), policy=parent)
    proposal = training_run['policy_proposal']
    if proposal is None:
        # Explicitly retain parent rather than inventing a successful model proposal.
        report = {'benchmark_version': BENCHMARK_VERSION,
                  'evidence_kind': 'synthetic_mechanism_test', 'provider': provider,
                  'status': 'no_policy_proposed', 'training_run': training_run,
                  'live_model_calls': training_run['calls'] if client is not None else 0}
        (output_dir / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
        return report
    # Fixture provenance is changed explicitly; an offline fake must not become llm_proposed.
    challenger = parent.revise(tuple(proposal['rules']),
                               source='llm_proposed' if client is not None else 'scripted')
    training_metrics = {'parent': score_policy(parent, training),
                        'challenger': score_policy(challenger, training)}
    write_policy(output_dir / 'policies', challenger)
    referee_path = write_content_addressed(output_dir / 'referees', REFEREE.model_dump())
    frozen_manifest = {'parent_hash': parent.digest, 'challenger_hash': challenger.digest,
                       'referee_hash': referee_path.stem, 'generator_version': BENCHMARK_VERSION,
                       'training_case_hashes': [case.digest for case in training],
                       'heldout_seeds': list(range(5000, 5048)),
                       'heldout_families': ['clean', 'saturation', 'feature_shift', 'combined']}
    freeze_path = write_content_addressed(output_dir / 'freezes', frozen_manifest)
    # No policy or threshold edits occur after this frozen manifest is written.
    families = ('clean', 'saturation', 'feature_shift', 'combined')
    heldout = [generate_case(seed, families[index % 4])
               for index, seed in enumerate(range(5000, 5048))]
    metrics, receipt = judge_policies(parent, challenger, heldout)
    write_content_addressed(output_dir / 'receipts', receipt.model_dump())
    report = {'benchmark_version': BENCHMARK_VERSION, 'status': 'completed',
              'evidence_kind': 'synthetic_mechanism_test', 'provider': provider,
              'live_model_calls': training_run['calls'] if client is not None else 0,
              'training_case_hashes': [case.digest for case in training],
              'heldout_case_hashes': [case.digest for case in heldout],
              'policy_frozen_before_heldout_generation': True,
              'freeze_manifest_hash': freeze_path.stem, 'frozen_referee_hash': referee_path.stem,
              'training': training_metrics, 'heldout': metrics,
              'training_run': training_run, 'promotion_receipt': receipt.model_dump(),
              'interpretation': 'Generated numerical failure cases test persistent diagnostic rules. '
              'They do not estimate biological prediction quality or clinical benefit. '
              'Heldout cases vary seeds, sizes and magnitudes within four predefined families. '
              + ('The offline policy proposal uses a scripted fixture.'
                 if client is None else 'The proposal used the recorded Crusoe model; '
                 'heldout execution uses the fixed rule interpreter.')}
    (output_dir / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def load_promoted_policy(policy_directory, receipt_path, parent):
    """Reload a policy only through an accepted, hash-verified referee receipt.

    Content addressing detects accidental mutation. It is not a digital signature
    and does not authenticate an adversary who can replace the entire artifact set.
    """
    receipt_path = Path(receipt_path)
    raw = json.loads(receipt_path.read_text())
    if content_hash(raw) != receipt_path.stem:
        raise ValueError('Promotion receipt hash mismatch')
    receipt = PromotionReceipt.model_validate_json(json.dumps(raw))
    gates_pass = (receipt.heldout_cases >= REFEREE.minimum_heldout_cases and
                  receipt.challenger_success_rate - receipt.parent_success_rate >=
                  REFEREE.minimum_success_rate_gain and
                  receipt.challenger_false_positive_rate - receipt.parent_false_positive_rate <=
                  REFEREE.maximum_false_positive_rate_increase)
    if (not receipt.accepted or not gates_pass or receipt.parent_hash != parent.digest or
            receipt.referee_hash != content_hash(REFEREE.model_dump())):
        raise ValueError('Promotion receipt is not accepted by the frozen referee')
    policy = read_policy(Path(policy_directory) / (receipt.challenger_hash + '.json'))
    if policy.parent_hash != parent.digest or policy.version != parent.version + 1:
        raise ValueError('Policy lineage does not match the receipt')
    return policy


def add_provider_arguments(parser):
    parser.add_argument('--provider', choices=('offline', 'crusoe'), default='offline')
    parser.add_argument('--max-actions', type=int, choices=range(1, 9), default=6)
    parser.add_argument('--model-id', default='nvidia/NVIDIA-Nemotron-3.5-Super-VL')
    parser.add_argument('--max-completion-tokens', type=int, default=8192)
    parser.add_argument('--max-total-tokens', type=int, default=60000)
    parser.add_argument('--timeout-seconds', type=float, default=60.)


def configured_live_client(args, parser):
    key = os.environ.get('CRUSOE_API_KEY', '')
    if not key:
        parser.error('Set CRUSOE_API_KEY explicitly for --provider crusoe; .env is not loaded')
    if not 1000 <= args.max_total_tokens <= 100000 or not 1 <= args.timeout_seconds <= 120:
        parser.error('Total tokens must be 1000..100000 and timeout 1..120 seconds')
    client = CrusoeClient(key, max_tokens=args.max_completion_tokens,
                          max_total_tokens=args.max_total_tokens, max_retries=0,
                          timeout=args.timeout_seconds, parallel_tool_calls=False)
    try:
        client.discover_models(nemotron_only=False)
        model = client.configure_model(args.model_id)
        if model.tools_supported is None and not client.run_smoke().passed:
            raise ValueError('Model failed the function-tool capability check')
        return client
    except BaseException:
        client.close()
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('reports/agent-policy-smoke'))
    add_provider_arguments(parser)
    args = parser.parse_args(argv)
    if args.provider == 'offline':
        result = run_policy_smoke(args.output, max_actions=args.max_actions)
    else:
        with configured_live_client(args, parser) as client:
            result = run_policy_smoke(args.output, client=client, max_actions=args.max_actions)
    print(json.dumps({'status': result['status'], 'provider': result['provider'],
                      'report': str(args.output / 'report.json'),
                      'promotion': result.get('promotion_receipt', {}).get('accepted', False)}))
    return 0 if result['status'] == 'completed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
