"""Run bounded activation planning on synthetic cases or published aggregate replay.

Default is offline and synthetic. --mode replay reads only two public report JSON
files; it performs no biological fits or outcome reads. --provider crusoe enables
explicit billed calls. Generated policies are proposals until separately promoted.
"""

import argparse
import json
from pathlib import Path

from escape_lab.activation_planner import ActivationPlanner
from escape_lab.activation_replay import load_development_replay
from escape_lab.policy_benchmark import (
    ScriptedDiagnosticClient,
    add_provider_arguments,
    configured_live_client,
    load_promoted_policy,
    run_policy_smoke,
)
from escape_lab.research_policy import ResearchPolicy, read_policy, write_policy


class ReplayDiagnosticClient(ScriptedDiagnosticClient):
    def chat(self, messages, tools, *, prompt_version):
        context = json.loads(messages[1]['content'])
        used = {record['action'] for record in context['history']}
        if {'inspect_score_saturation', 'compare_rankings', 'inspect_feature_shift'} <= used \
                and 'request_sequence_ablation' not in used:
            from types import SimpleNamespace
            return SimpleNamespace(finish_reason='tool_calls', message={'tool_calls': [{
                'id': 'offline-replay', 'function': {'name': 'request_sequence_ablation',
                                                    'arguments': '{"target": "bound_partition"}'}}]})
        return super().chat(messages, tools, prompt_version=prompt_version)


def run_replay(root, output, *, client=None, policy=None, max_actions=6, retention_verified=False):
    if policy is not None and policy.version > 1 and not retention_verified:
        raise ValueError('A retained policy requires a verified promotion receipt')
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError('Choose a new output directory; previous attempts are immutable')
    output.mkdir(parents=True, exist_ok=True)
    executor = load_development_replay(root)
    provider = 'crusoe' if client is not None else 'scripted_offline'
    actual_client = client if client is not None else ReplayDiagnosticClient()
    policy = policy or ResearchPolicy.initial()
    write_policy(output / 'policies', policy)
    result = ActivationPlanner(actual_client, max_actions=max_actions,
        proposal_source='llm_proposed' if client is not None else 'scripted').run(executor, policy=policy)
    proposal = result['policy_proposal']
    if proposal is not None:
        proposed = ResearchPolicy.model_validate_json(json.dumps(proposal))
        write_policy(output / 'proposals', proposed)
    report = {'evidence_kind': 'fixed_development_aggregate_replay', 'provider': provider,
              'scope': executor.scope, 'metric': executor.metric,
              'source_hashes': executor.source_hashes, 'new_fits': 0,
              'external_datasets_accessed': False, 'new_scientific_result': False,
              'policy_promoted': False,
              'input_policy_status': 'synthetic_referee_promoted' if retention_verified else 'initial',
              'run': result}
    (output / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('synthetic', 'replay'), default='synthetic')
    parser.add_argument('--root', type=Path, default=Path.cwd())
    parser.add_argument('--output', type=Path, default=Path('reports/activation-diagnostic-smoke'))
    parser.add_argument('--policy', type=Path, help='Hash-named retained policy artifact for replay')
    parser.add_argument('--promotion-receipt', type=Path, help='Accepted hash-named referee receipt')
    parser.add_argument('--parent-policy', type=Path, help='Parent artifact; omitted means initial policy')
    add_provider_arguments(parser)
    args = parser.parse_args(argv)
    if args.mode == 'synthetic' and args.policy is not None:
        parser.error('--policy is for replay; synthetic comparison freezes its initial parent')
    if (args.promotion_receipt or args.parent_policy) and args.policy is None:
        parser.error('--promotion-receipt and --parent-policy require --policy')
    policy = read_policy(args.policy) if args.policy is not None else None
    retention_verified = False
    if policy is not None and policy.version > 1:
        if args.promotion_receipt is None:
            parser.error('Reusing a retained policy requires --promotion-receipt')
        parent = read_policy(args.parent_policy) if args.parent_policy else ResearchPolicy.initial()
        verified = load_promoted_policy(args.policy.parent, args.promotion_receipt, parent)
        if verified.digest != policy.digest:
            parser.error('Promotion receipt belongs to a different policy')
        retention_verified = True

    def execute(client=None):
        if args.mode == 'synthetic':
            return run_policy_smoke(args.output, client=client, max_actions=args.max_actions)
        return run_replay(args.root, args.output, client=client, policy=policy,
                          max_actions=args.max_actions, retention_verified=retention_verified)

    if args.provider == 'crusoe':
        with configured_live_client(args, parser) as client:
            result = execute(client)
    else:
        result = execute()
    status = result.get('status', result.get('run', {}).get('status'))
    print(json.dumps({'status': status, 'provider': result['provider'],
                      'report': str(args.output / 'report.json')}))
    return 0 if status in ('completed', 'stopped') else 1


if __name__ == '__main__':
    raise SystemExit(main())
