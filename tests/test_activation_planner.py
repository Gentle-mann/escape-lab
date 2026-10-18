import json
from types import SimpleNamespace

import numpy as np
import pytest

from escape_lab.activation_planner import ActivationDiagnostics, ActivationPlanner
from escape_lab.research_policy import ResearchPolicy


class FakeClient:
    model_id = 'offline/test'
    total_tokens = 0

    def __init__(self, actions):
        self.actions = iter(actions)
        self.requests = []

    def chat(self, messages, tools, *, prompt_version):
        self.requests.append(json.loads(json.dumps({'messages': messages, 'tools': tools})))
        action, args = next(self.actions)
        if action in ('inspect_score_saturation', 'inspect_feature_shift', 'compare_rankings',
                       'request_sequence_ablation') and args == {}:
            args = {'target': 'bound_partition'}
        return SimpleNamespace(finish_reason='tool_calls', message={'content': 'DO NOT LOG',
            'tool_calls': [{'id': 'call-1', 'function': {'name': action,
                                                        'arguments': json.dumps(args)}}]})


def fixture():
    def ablation():
        return {'combined_ap': .9, 'wt_only_ap': .6, 'delta_only_ap': .8, 'private_variant': 'HIDDEN'}

    return ActivationDiagnostics(scope='synthetic', probabilities=np.ones(6),
        ranking_scores=np.array([50., 45., 40., 39., 38., 37.]),
        labels=np.array([True, True, False, False, False, False]),
        scaled_features=np.array([[20., 1.], [1., 0.]]),
        training_variances=np.array([1e-20, 1.]),
        ablation_runner=ablation)


def test_planner_executes_tools_reuses_aggregates_and_proposes_a_rule():
    client = FakeClient([
        ('inspect_score_saturation', {}), ('compare_rankings', {}),
        ('propose_policy_revision', {'rules': ['inspect_saturation_before_ranking']}),
        ('stop', {'reason': 'diagnostics_complete'}),
    ])
    result = ActivationPlanner(client).run(fixture(), policy=ResearchPolicy.initial())
    assert result['status'] == 'stopped'
    assert result['policy_proposal']['rules'] == ['inspect_saturation_before_ranking']
    assert result['history'][1]['result']['ranking_ap'] == 1.
    assert result['history'][1]['result']['probability_ap'] == pytest.approx(1 / 3)
    sent = json.dumps(client.requests)
    assert 'DO NOT LOG' not in sent and 'HIDDEN' not in sent
    contexts = json.dumps([json.loads(request['messages'][1]['content']) for request in client.requests])
    assert 'ranking_scores' not in contexts and 'labels' not in contexts
    last = json.loads(client.requests[-1]['messages'][1]['content'])
    assert last['history'][1]['result']['delta_ap'] > 0
    assert last['policy']['version'] == 1


def test_stop_does_not_force_extra_tools_or_live_calls():
    client = FakeClient([('stop', {'reason': 'no_justified_experiment'})])
    result = ActivationPlanner(client).run(fixture())
    assert result['status'] == 'stopped'
    assert result['calls'] == 1 and len(client.requests) == 1


def test_invalid_unknown_tool_and_unearned_policy_fail_closed():
    for action, args in [('read_file', {'path': '/secret'}),
                         ('propose_policy_revision', {'rules': ['inspect_saturation_before_ranking']})]:
        result = ActivationPlanner(FakeClient([(action, args)])).run(fixture())
        assert result['status'] == 'invalid_response'
        assert result['policy_proposal'] is None


def test_executor_scope_and_invalid_probabilities_are_rejected():
    with pytest.raises(ValueError, match='scope'):
        ActivationDiagnostics(scope='external', probabilities=[.2], ranking_scores=[.3], labels=[True])
    with pytest.raises(ValueError, match='probabilities'):
        ActivationDiagnostics(scope='synthetic', probabilities=[8.], ranking_scores=[.3], labels=[True])


def test_ablation_executes_only_once_and_drops_arbitrary_output():
    calls = []

    def ablate():
        calls.append(True)
        return {'combined_ap': .8, 'wt_only_ap': .5, 'delta_only_ap': .7,
                'raw_labels': [1], 'notes': 'PRIVATE'}

    executor = ActivationDiagnostics(scope='synthetic', probabilities=[.2, .8],
                                     ranking_scores=[-1., 1.], labels=[False, True],
                                     ablation_runner=ablate)
    client = FakeClient([('request_sequence_ablation', {}),
                         ('request_sequence_ablation', {})])
    result = ActivationPlanner(client).run(executor)
    assert result['status'] == 'invalid_response'
    assert calls == [True]
    assert 'PRIVATE' not in json.dumps(client.requests)


def test_actual_crusoe_adapter_sends_bounded_function_tools_without_arrays():
    import httpx

    from escape_lab.crusoe import CrusoeClient

    requests = []

    def handle(request):
        if request.url.path.endswith('/models'):
            return httpx.Response(200, json={'data': [
                {'id': 'nvidia/test-nemotron', 'supported_parameters': ['tools']}]})
        payload = json.loads(request.content)
        requests.append(payload)
        return httpx.Response(200, json={'choices': [{'finish_reason': 'tool_calls', 'message': {
            'tool_calls': [{'id': 'c', 'function': {'name': 'stop', 'arguments': json.dumps({
                'reason': 'diagnostics_complete'})}}]}}], 'usage': {'total_tokens': 37}})

    with CrusoeClient('test-key', transport=httpx.MockTransport(handle)) as client:
        client.discover_models()
        client.configure_model('nvidia/test-nemotron')
        result = ActivationPlanner(client).run(fixture())
    assert result['status'] == 'stopped' and result['budget_charge_tokens'] == 37
    assert len(requests) == 1 and requests[0]['tool_choice'] == 'required'
    context = json.loads(requests[0]['messages'][1]['content'])
    assert not {'labels', 'ranking_scores', 'probabilities'} & set(context)


def test_real_report_replay_uses_matching_mean_fold_metrics_without_fits():
    from pathlib import Path

    from escape_lab.activation_replay import load_development_replay

    root = Path(__file__).resolve().parents[1]
    if not (root / 'docs/archive/reports/model-improvements/results.json').exists():
        pytest.skip('Optional published follow-up report not present')
    executor = load_development_replay(root)
    assert executor.scope == 'fixed-development-replay'
    comparison = executor.execute('compare_rankings')
    assert comparison['ranking_ap'] == pytest.approx(.11134574838153616)
    assert comparison['probability_ap'] == pytest.approx(.09732091648081248)
    ablation = executor.execute('request_sequence_ablation')
    assert set(ablation) == {'combined_ap', 'wt_only_ap', 'delta_only_ap'}
    assert executor.new_fits == 0


def test_executor_snapshots_inputs_and_does_not_trust_supplied_score_diagnostics():
    labels = np.array([False, True])
    features = np.array([[1., 2.], [0., 1.]])
    executor = ActivationDiagnostics(scope='synthetic', probabilities=[.2, .8],
        ranking_scores=[-1., 1.], labels=labels, scaled_features=features,
        prediction_diagnostics={'probability_saturated_fraction': 1., 'score_count': 900})
    labels[:] = False
    features[:] = 100.
    score = executor.execute('inspect_score_saturation')
    assert score['score_count'] == 2 and score['probability_saturated_fraction'] == 0
    assert executor.execute('compare_rankings')['ranking_ap'] == 1.
    assert executor.execute('inspect_feature_shift')['feature_shift_fraction'] == 0.
    assert not executor.labels.flags.writeable


def test_replay_rejects_mismatched_scoring_populations(tmp_path):
    from pathlib import Path

    from escape_lab.activation_replay import DevelopmentReplay
    root = Path(__file__).resolve().parents[1]
    source = root / 'docs/archive/reports/model-improvements/results.json'
    if not source.exists():
        pytest.skip('Optional report absent')
    model = json.loads(source.read_text())
    ranking = json.loads((root / 'docs/archive/reports/mechanistic-pilot/'
                          'fgfr2-sequence-ap-v1/results.json').read_text())
    model['metrics']['legacy_wt']['by_fold']['pilot-0']['n'] -= 1
    with pytest.raises(ValueError, match='case counts'):
        DevelopmentReplay(model, ranking, source_hashes={})


def test_parallel_and_unknown_tool_failures_are_audited_without_raw_response():
    class Client:
        model_id = 'offline/test'
        total_tokens = 0

        def chat(self, messages, tools, *, prompt_version):
            return SimpleNamespace(finish_reason='tool_calls', message={'content': 'SECRET',
                'tool_calls': [{'function': {'name': 'inspect_score_saturation', 'arguments': '{}'}},
                               {'function': {'name': '/PRIVATE/VARIANT', 'arguments': 'SECRET'}}]})

    result = ActivationPlanner(Client()).run(fixture())
    assert result['status'] == 'invalid_response'
    assert result['error_category'] == 'multiple_tool_calls'
    assert result['response_audit'][0]['tool_call_count'] == 2
    assert result['response_audit'][0]['known_actions'] == ['inspect_score_saturation', 'unknown']
    assert 'SECRET' not in json.dumps(result) and 'PRIVATE' not in json.dumps(result)


def test_crusoe_parallel_tool_setting_is_optional_and_strict():
    import httpx

    from escape_lab.crusoe import CrusoeClient
    payloads = []

    def handle(request):
        if request.url.path.endswith('/models'):
            return httpx.Response(200, json={'data': [
                {'id': 'test', 'supported_parameters': ['tools', 'parallel_tool_calls']}]})
        payloads.append(json.loads(request.content))
        return httpx.Response(200, json={'choices': [{'message': {'content': None},
                                                     'finish_reason': 'stop'}]})

    for setting in (None, False):
        with CrusoeClient('test', transport=httpx.MockTransport(handle),
                          parallel_tool_calls=setting) as client:
            client.discover_models(nemotron_only=False)
            client.configure_model('test')
            client.chat([], [])
    assert 'parallel_tool_calls' not in payloads[0]
    assert payloads[1]['parallel_tool_calls'] is False
    with pytest.raises(ValueError):
        CrusoeClient('test', parallel_tool_calls=0)


def test_numeric_tools_require_explicit_bound_partition_and_reject_wrappers():
    for bad in ({'arguments': '{}'}, {'target': 'external'}, {'target': 'bound_partition', 'rows': []}):
        result = ActivationPlanner(FakeClient([('inspect_score_saturation', bad)])).run(fixture())
        assert result['status'] == 'invalid_response'
        assert result['error_category'] == 'invalid_tool_arguments'
        assert result['history'] == []
    client = FakeClient([('inspect_score_saturation', {'target': 'bound_partition'}),
                         ('stop', {'reason': 'diagnostics_complete'})])
    result = ActivationPlanner(client).run(fixture())
    assert result['status'] == 'stopped'
    assert result['prompt_version'].endswith('v2')
    schema = next(tool['function']['parameters'] for tool in client.requests[0]['tools']
                  if tool['function']['name'] == 'inspect_score_saturation')
    assert schema['required'] == ['target']
    assert schema['properties']['target']['const'] == 'bound_partition'
    assert schema['additionalProperties'] is False
