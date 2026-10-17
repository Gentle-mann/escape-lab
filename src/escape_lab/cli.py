"""Explicit CLI: local preparation never loads credentials or calls an LLM."""

import argparse
import json
import os
import secrets
from dataclasses import asdict
from pathlib import Path

from dotenv import dotenv_values


def build_parser():
    parser = argparse.ArgumentParser(description='Escape Lab preclinical experimental replay')
    parser.add_argument('--root', type=Path, default=Path.cwd())
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('prepare', help='Download, verify and group the public dataset')
    commands.add_parser('status', help='Print nonsecret setup status')
    commands.add_parser('crusoe-models', help='Discover accessible Nemotron model IDs')
    commands.add_parser('crusoe-check', help='Explicitly run one tool-calling smoke test')
    improve = commands.add_parser('improve', help='Autonomously revise per-drug resistance predictors')
    improve.add_argument('--output', type=Path)
    improve.add_argument('--max-seconds', type=int, default=600)
    improve.add_argument('--model', help='Exact served model ID; overrides CRUSOE_MODEL')
    improve.add_argument('--request-timeout', type=int, default=60,
                         help='HTTP request timeout in seconds, from 1 to 90 (default: 60)')
    improve.add_argument('--max-completion-tokens', type=int, choices=[2048, 4096, 8192],
                         default=8192, help='Includes reasoning tokens used before the tool call')
    improve.add_argument('--development-only', action='store_true',
                         help='Run development and validation; never open the final test')
    for name in ['baseline', 'agent-eval']:
        command = commands.add_parser(name)
        command.add_argument('--stage', choices=['development', 'validation'], default='development')
        command.add_argument('--seeds', type=int, nargs='+', default=[0, 1, 2])
        command.add_argument('--budget', type=int, default=20)
        command.add_argument('--output', type=Path)
        command.add_argument('--policies', nargs='+', default=(['extra_trees', 'agent']
                             if name == 'agent-eval' else
                             ['random', 'gene_prior', 'logistic', 'extra_trees',
                              'extra_explore', 'extra_diverse']))
        if name == 'agent-eval':
            command.add_argument('--max-total-tokens', type=int, default=250_000)
    serve = commands.add_parser('serve', help='Start loopback-only experiment vault API')
    serve.add_argument('--port', type=int, default=8017)
    return parser


def credentials(root):
    values = dotenv_values(root / '.env')
    key = os.environ.get('CRUSOE_API_KEY') or values.get('CRUSOE_API_KEY')
    model = os.environ.get('CRUSOE_MODEL') or values.get('CRUSOE_MODEL')
    return key, model


def project_status(root):
    key, model = credentials(root)
    manifest_path = root / 'data/processed/manifest.json'
    result = {'project_root': str(root), 'crusoe_key_configured': bool(key and key.strip()),
              'configured_model': model or None, 'dataset_prepared': manifest_path.exists(),
              'final_test': 'sealed'}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        for name in ['paired_nucleotide_rows', 'unique_protein_variants', 'split_counts',
                     'split_group_counts', 'candidates_sha256']:
            result[name] = manifest[name]
    receipt = root / 'data/processed/final-evaluation.json'
    if receipt.exists():
        try:
            saved = json.loads(receipt.read_text())
            result['final_test'] = saved['status']
            result['final_evaluation_run'] = saved.get('run_id')
        except (ValueError, KeyError):
            result['final_test'] = 'reserved_unknown'
    return result


def configured_client(root, max_total_tokens=250_000, model_id=None):
    from escape_lab.crusoe import CrusoeClient, CrusoeError
    key, requested = credentials(root)
    requested = model_id or requested
    if not key or not key.strip():
        raise CrusoeError('CRUSOE_API_KEY is not configured. Add it to the local .env file.')
    client = CrusoeClient(key, max_tokens=4096, max_total_tokens=max_total_tokens)
    models = client.discover_models()
    allowed = [m for m in models if m.tools_supported is not False]
    if requested:
        client.configure_model(requested)
    elif allowed:
        preferred = next((m for m in allowed if 'super' in m.id.lower()), allowed[0])
        client.configure_model(preferred.id)
    else:
        client.close()
        raise CrusoeError('No compatible Nemotron model was returned for this account.')
    return client


def main(argv=None):
    args = build_parser().parse_args(argv)
    root = args.root.resolve()
    if args.command == 'status':
        print(json.dumps(project_status(root), indent=2))
    elif args.command == 'prepare':
        from escape_lab.data import download_sources, prepare_dataset
        download_sources(root / 'data/raw')
        print(json.dumps(prepare_dataset(root / 'data/raw', root / 'data/processed'), indent=2))
    elif args.command in {'crusoe-models', 'crusoe-check'}:
        from escape_lab.crusoe import CrusoeError
        try:
            with configured_client(root) as client:
                if args.command == 'crusoe-models':
                    print(json.dumps([m.model_dump() for m in client.models.values()
                                      if 'nemotron' in m.id.lower()], indent=2))
                else:
                    result = client.run_smoke()
                    print(json.dumps({'smoke': asdict(result), 'usage': client.usage_records}, indent=2))
                    if not result.passed:
                        raise SystemExit(1)
        except CrusoeError as error:
            raise SystemExit(str(error)) from None
    elif args.command == 'improve':
        import httpx

        from escape_lab.crusoe import CrusoeError
        from escape_lab.improvement import run_improvement
        from escape_lab.improvement_agent import PipelineProposer
        try:
            if not 1 <= args.request_timeout <= 90:
                raise ValueError('Request timeout must be between 1 and 90 seconds.')
            with configured_client(root, max_total_tokens=50_000, model_id=args.model) as client:
                client.max_retries = 0
                client.max_tokens = args.max_completion_tokens
                client._http.timeout = httpx.Timeout(args.request_timeout)
                client.request_timeout_seconds = args.request_timeout
                if (client.models[client.model_id].tools_supported is None
                        and not client.run_smoke().passed):
                    raise CrusoeError('Tool capability test failed.')
                result = run_improvement(root / 'data/processed',
                    args.output or root / 'reports/improvement', PipelineProposer(client),
                    max_seconds=args.max_seconds, evaluate_final=not args.development_only)
                print(json.dumps({key: result[key] for key in
                    ['run_id', 'run_dir', 'outcome', 'elapsed_seconds', 'within_time_budget',
                     'selection', 'final'] if key in result}, indent=2))
                if result['outcome'] in {'incomplete', 'time_budget_exceeded'}:
                    raise SystemExit(1)
        except (CrusoeError, ValueError) as error:
            raise SystemExit(str(error)) from None
    elif args.command in {'baseline', 'agent-eval'}:
        from escape_lab.benchmark import run_benchmark
        output = args.output or root / 'reports' / args.command
        kwargs = {'data_dir': root / 'data/processed', 'output_dir': output,
                  'stage': args.stage, 'policies': tuple(args.policies),
                  'seeds': tuple(args.seeds), 'budget': args.budget}
        if args.command == 'baseline':
            result = run_benchmark(**kwargs)
        else:
            from escape_lab.agent import NemotronInvestigator
            from escape_lab.crusoe import CrusoeError
            try:
                with configured_client(root, args.max_total_tokens) as client:
                    smoke = client.run_smoke()
                    if not smoke.passed:
                        raise CrusoeError('Tool-calling smoke test failed; evaluation not started.')
                    kwargs['investigator'] = NemotronInvestigator(client, max_decisions=1000)
                    result = run_benchmark(**kwargs)
            except CrusoeError as error:
                raise SystemExit(str(error)) from None
        print(json.dumps({k: v for k, v in result.items() if k != 'trials'}, indent=2))
        print(f'Report directory: {output.resolve()}')
    elif args.command == 'serve':
        import uvicorn

        from escape_lab.api import create_app
        from escape_lab.vault import Vault
        data = root / 'data/processed'
        private = root / 'reports/private'
        private.mkdir(parents=True, exist_ok=True)
        token_file = private / 'admin-token'
        if not token_file.exists():
            token_file.write_text(secrets.token_urlsafe(32))
            token_file.chmod(0o600)
        vault = Vault(private / 'vault.sqlite', json.loads((data/'candidates.json').read_text()),
                      json.loads((data/'outcomes.json').read_text()))
        app = create_app(vault, token_file.read_text().strip())
        print(f'Vault bound to http://127.0.0.1:{args.port}; admin token stored privately.')
        uvicorn.run(app, host='127.0.0.1', port=args.port)


if __name__ == '__main__':
    main()
