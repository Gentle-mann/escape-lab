import json

import pytest
from pydantic import ValidationError

from escape_lab.policy_benchmark import run_policy_smoke
from escape_lab.research_policy import ResearchPolicy, read_policy, write_policy


def test_policy_is_frozen_content_addressed_and_rejects_unknown_rules(tmp_path):
    parent = ResearchPolicy.initial()
    challenger = parent.revise(('inspect_saturation_before_ranking',), source='scripted')
    assert challenger.version == 2 and challenger.parent_hash == parent.digest
    path = write_policy(tmp_path, challenger)
    assert path.name == challenger.digest + '.json'
    assert read_policy(path) == challenger
    assert write_policy(tmp_path, challenger) == path
    with pytest.raises(ValidationError):
        challenger.version = 9
    payload = json.loads(path.read_text())
    payload['version'] = 12
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match='hash'):
        read_policy(path)
    with pytest.raises(ValidationError):
        parent.revise(('modify_final_labels',), source='scripted')


def test_policy_cannot_hold_labels_free_text_or_biological_identity():
    with pytest.raises(ValidationError):
        ResearchPolicy(version=1, parent_hash=None, rules=(), source='scripted', labels=[1])
    with pytest.raises(ValidationError):
        ResearchPolicy(version=1, parent_hash=None, rules=('FGFR2-V564F',), source='scripted')
    with pytest.raises(ValidationError):
        ResearchPolicy(version=2, parent_hash=None, rules=(), source='scripted')


def test_offline_smoke_freezes_before_distinct_heldout_cases_and_writes_receipt(tmp_path):
    report = run_policy_smoke(tmp_path)
    assert report['evidence_kind'] == 'synthetic_mechanism_test'
    assert report['provider'] == 'scripted_offline'
    assert report['live_model_calls'] == 0
    assert report['training_case_hashes']
    assert not set(report['training_case_hashes']) & set(report['heldout_case_hashes'])
    assert report['policy_frozen_before_heldout_generation'] is True
    assert report['heldout']['challenger']['success_rate'] > report['heldout']['parent']['success_rate']
    receipt = report['promotion_receipt']
    assert receipt['accepted'] is True
    assert receipt['referee_hash'] == report['frozen_referee_hash']
    assert read_policy(tmp_path / 'policies' / (receipt['challenger_hash'] + '.json')).version == 2
    assert (tmp_path / 'report.json').is_file()
    # The archived artifact is a pure rule object, not a model/label container.
    payload = json.loads((tmp_path / 'policies' / (receipt['challenger_hash'] + '.json')).read_text())
    assert set(payload) == {'version', 'parent_hash', 'rules', 'source'}


def test_promotion_receipt_reloads_policy_and_rejects_tampering(tmp_path):
    from escape_lab.policy_benchmark import load_promoted_policy
    report = run_policy_smoke(tmp_path)
    receipt = next((tmp_path / 'receipts').glob('*.json'))
    loaded = load_promoted_policy(tmp_path / 'policies', receipt, ResearchPolicy.initial())
    assert loaded.digest == report['promotion_receipt']['challenger_hash']
    value = json.loads(receipt.read_text())
    value['challenger_success_rate'] = .2
    receipt.write_text(json.dumps(value))
    with pytest.raises(ValueError, match='hash'):
        load_promoted_policy(tmp_path / 'policies', receipt, ResearchPolicy.initial())


def test_frozen_referee_rejects_duplicate_cases_and_empty_improvement():
    from escape_lab.policy_benchmark import REFEREE, generate_case, judge_policies
    parent = ResearchPolicy.initial()
    with pytest.raises(ValidationError):
        REFEREE.minimum_success_rate_gain = 0
    unchanged = parent.revise((), source='scripted')
    cases = [generate_case(seed, 'combined') for seed in range(24)]
    _, receipt = judge_policies(parent, unchanged, cases)
    assert not receipt.accepted
    with pytest.raises(ValueError, match='distinct'):
        judge_policies(parent, unchanged, [cases[0]] * 24)


def test_hash_valid_receipt_with_false_acceptance_still_fails_referee(tmp_path):
    from escape_lab.policy_benchmark import load_promoted_policy
    from escape_lab.research_policy import write_content_addressed
    report = run_policy_smoke(tmp_path)
    receipt = dict(report['promotion_receipt'])
    receipt['accepted'] = True
    receipt['heldout_cases'] = 0
    path = write_content_addressed(tmp_path / 'forged', receipt)
    with pytest.raises(ValueError, match='frozen referee'):
        load_promoted_policy(tmp_path / 'policies', path, ResearchPolicy.initial())


def test_smoke_never_overwrites_an_existing_attempt(tmp_path):
    run_policy_smoke(tmp_path)
    before = (tmp_path / 'report.json').read_bytes()
    with pytest.raises(FileExistsError, match='previous attempts'):
        run_policy_smoke(tmp_path)
    assert (tmp_path / 'report.json').read_bytes() == before


def test_replay_cli_refuses_unpromoted_policy_before_reading_reports(tmp_path):
    from escape_lab.activation_diagnostic_cli import main
    policy = ResearchPolicy.initial().revise(('inspect_saturation_before_ranking',), source='scripted')
    path = write_policy(tmp_path / 'policies', policy)
    with pytest.raises(SystemExit) as error:
        main(['--mode', 'replay', '--root', str(tmp_path / 'missing'),
              '--policy', str(path), '--output', str(tmp_path / 'run')])
    assert error.value.code == 2
    assert not (tmp_path / 'run').exists()
