import json
from types import SimpleNamespace

import numpy as np

from escape_lab.experiment_planner import ExperimentPlanner, sanitize_feedback
from escape_lab.planner_diagnostics import aggregate_fit_diagnostics, numerical_diagnostics

BASE = 'chemistry-logistic-c0.1'
NEXT = 'structure-logistic-c0.1'
CATALOG = {
    BASE: {'representation': 'chemistry', 'model': 'logistic', 'C': .1, 'max_depth': 2},
    NEXT: {'representation': 'structure', 'model': 'logistic', 'C': .1, 'max_depth': 2},
}


class FakeClient:
    model_id = 'nvidia/test-nemotron'
    total_tokens = 0

    def __init__(self, actions):
        self.actions = iter(actions)
        self.requests = []

    def chat(self, messages, tools, *, prompt_version):
        self.requests.append({'messages': json.loads(json.dumps(messages)), 'tools': tools,
                              'prompt_version': prompt_version})
        name, arguments = next(self.actions)
        return SimpleNamespace(finish_reason='tool_calls', message={'tool_calls': [{
            'id': 'call-1', 'function': {'name': name, 'arguments': json.dumps(arguments)}}]})


def test_bounded_diagnostics_and_reasons_survive_without_identities_or_text():
    marker = 'PRIVATE-ROW-OUTER-LABEL'
    state = {'history': [{'candidate_id': BASE, 'status': 'rejected', 'metrics': {'macro_ap': .2},
                         'rejection_reasons': ['pem_brier_regression', marker],
                         'fit_diagnostics': {'probability_saturated_fraction': .8,
                                             'scaled_feature_abs_max': 1e30,
                                             'feature_shift_fraction': float('nan'),
                                             'feature_names': marker, 'labels': [True]}}]}
    safe = sanitize_feedback(state, CATALOG, include_diagnostics=True)
    assert marker not in json.dumps(safe)
    assert safe['history'][0]['rejection_reasons'] == ['pem_brier_regression']
    assert safe['history'][0]['fit_diagnostics'] == {'probability_saturated_fraction': .8}


def test_new_feedback_reaches_next_call_and_stop_is_explicit():
    client = FakeClient([('stop_search', {'reason': 'diagnostics_complete'})])
    planner = ExperimentPlanner(client, protocol='diagnostics-v3')
    state = {'incumbent': {'candidate_id': BASE, 'metrics': {'macro_ap': .2}},
             'history': [{'candidate_id': BASE, 'metrics': {'macro_ap': .2},
                          'rejection_reasons': ['no_pooled_ap_improvement'],
                          'fit_diagnostics': {'probability_saturated_fraction': .5}}]}
    result = planner.propose(CATALOG, state)
    sent = json.loads(client.requests[0]['messages'][1]['content'])
    assert sent['history'][0]['rejection_reasons'] == ['no_pooled_ap_improvement']
    assert sent['history'][0]['fit_diagnostics']['probability_saturated_fraction'] == .5
    assert result.status == 'stopped'
    assert result.prompt_version.endswith('v3')
    assert planner.propose(CATALOG, state).status == 'stopped'
    assert len(client.requests) == 1


def test_legacy_prompt_and_feedback_remain_without_new_diagnostic_fields():
    args = {'candidate_id': NEXT, 'parent_id': BASE, 'hypothesis': 'A',
            'expected_effect': 'B', 'rejection_condition': 'C'}
    client = FakeClient([('propose_experiment', args)])
    planner = ExperimentPlanner(client)
    state = {'incumbent': {'candidate_id': BASE, 'metrics': {'macro_ap': .2},
                          'rejection_reasons': ['undefined_metric'],
                          'fit_diagnostics': {'probability_saturated_fraction': .8}}}
    planner.propose(CATALOG, state)
    assert 'fit_diagnostics' not in json.dumps(client.requests)
    assert client.requests[0]['prompt_version'].endswith('v2')
    assert len(client.requests[0]['tools']) == 1


def test_model_diagnostics_are_aggregate_and_finite():
    result = numerical_diagnostics(np.array([0., 1., .3, .3]),
                                   scaled_features=np.array([[20., 0.], [1., 1.]]),
                                   training_variances=np.array([0., 1e-16, 1.]))
    assert result['probability_saturated_fraction'] == .5
    assert result['probability_unique_fraction'] == .75
    assert result['feature_shift_fraction'] == .25
    assert result['near_constant_feature_count'] == 2
    result = aggregate_fit_diagnostics([result, {'coefficient_norm': 4., 'raw': 'secret'}])
    assert result['coefficient_norm_max'] == 4.
    assert 'secret' not in json.dumps(result)


def test_controller_passes_rejections_and_fit_feedback_then_honors_stop(monkeypatch, tmp_path):
    import escape_lab.revision_evaluation as evaluation
    from escape_lab.experiment_planner import ExperimentProposal

    states = []

    class Planner:
        def __init__(self):
            self.attempts = []

        def propose(self, catalog, state):
            states.append(state)
            status = 'proposed' if len(states) == 1 else 'stopped'
            return ExperimentProposal(NEXT if status == 'proposed' else None, BASE,
                                      '', '', '', status, None, 1, 'offline')

    def evaluate(candidate, tasks, labels, store, seed, deadline):
        metrics = {'macro_ap': .2, 'per_drug': {
            drug: {'average_precision': .2, 'brier': .1} for drug in ('pem', 'fut')},
            'by_fold': {f'inner-{i}': {'macro_ap': .2} for i in range(3)}}
        return {'metrics': metrics, 'diagnostics': [{'probability_saturated_fraction': .7}]}

    monkeypatch.setattr(evaluation, 'evaluate_cv', evaluate)
    rows = [{'candidate_id': str(i), 'group_id': f'g{i}'} for i in range(6)]
    labels = {str(i): {'pem': False, 'fut': False} for i in range(6)}
    result = evaluation.search_candidates(rows, labels, None, Planner(), 2, float('inf'), tmp_path)
    assert len(states) == 2
    rejected = states[1]['history'][-1]
    assert rejected['rejection_reasons'] == ['no_pooled_ap_improvement']
    assert rejected['fit_diagnostics']['probability_saturated_fraction'] == .7
    assert result['agent_status'] == 'stopped'
    assert result['unused_opportunities'] == 1
    assert result['failed_opportunities'] == 0


def test_missing_diagnostics_are_not_reported_as_clean_and_scales_use_minimum():
    result = aggregate_fit_diagnostics([
        {'minimum_training_scale': 1e-9, 'below_floor_feature_count': 10},
        {'minimum_training_scale': .1, 'below_floor_feature_count': 7},
    ])
    assert result['minimum_training_scale'] == 1e-9
    assert result['below_floor_feature_count'] == 10
    assert 'convergence_warning_fraction' not in result
    assert 'fallback_fraction' not in result
