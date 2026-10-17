import json

import pytest

from escape_lab.cli import build_parser, project_status


def test_final_evaluation_is_not_a_development_cli_option():
    with pytest.raises(SystemExit):
        build_parser().parse_args(['baseline', '--stage', 'final'])


def test_status_never_returns_credential_values(tmp_path):
    secret = 'secret-should-never-appear'
    (tmp_path / '.env').write_text(f'CRUSOE_API_KEY={secret}\n')
    result = project_status(tmp_path)
    assert result['crusoe_key_configured'] is True
    assert secret not in json.dumps(result)


def test_status_does_not_call_a_truncated_final_receipt_sealed(tmp_path):
    folder = tmp_path / 'data/processed'
    folder.mkdir(parents=True)
    (folder / 'final-evaluation.json').write_text('{')
    assert project_status(tmp_path)['final_test'] == 'reserved_unknown'


def test_improvement_has_a_bounded_default_and_explicit_sealed_mode():
    args = build_parser().parse_args(['improve', '--development-only'])
    assert args.development_only is True
    assert args.max_seconds == 600
    assert args.max_completion_tokens == 8192


def test_improvement_can_pin_model_and_request_timeout():
    args = build_parser().parse_args([
        'improve', '--model', 'nvidia/exact-nemotron-id', '--request-timeout', '60'])
    assert args.model == 'nvidia/exact-nemotron-id'
    assert args.request_timeout == 60


def test_explicit_model_selection_overrides_environment_and_catalog_order(monkeypatch, tmp_path):
    from types import SimpleNamespace

    from escape_lab.cli import configured_client

    (tmp_path / '.env').write_text('CRUSOE_API_KEY=fake-test-key\nCRUSOE_MODEL=other-model\n')

    class FakeClient:
        def __init__(self, api_key, **kwargs):
            self.model_id = None

        def discover_models(self):
            return [SimpleNamespace(id='other-model', tools_supported=True),
                    SimpleNamespace(id='exact-model', tools_supported=True)]

        def configure_model(self, requested):
            self.model_id = requested

    monkeypatch.setattr('escape_lab.crusoe.CrusoeClient', FakeClient)
    client = configured_client(tmp_path, model_id='exact-model')
    assert client.model_id == 'exact-model'
