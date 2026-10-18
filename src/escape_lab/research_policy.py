"""Content-addressed research policies: typed rules, never examples or outcomes."""

import hashlib
import json
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

Rule = Literal['inspect_saturation_before_ranking', 'inspect_feature_shift_before_promotion']
RULES = ('inspect_saturation_before_ranking', 'inspect_feature_shift_before_promotion')


def content_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


class ResearchPolicy(BaseModel):
    model_config = ConfigDict(frozen=True, strict=True, extra='forbid')
    version: int = Field(ge=1, le=1_000_000)
    parent_hash: str | None = Field(pattern=r'^[0-9a-f]{64}$')
    rules: tuple[Rule, ...] = Field(max_length=2)
    source: Literal['scripted', 'llm_proposed']

    @model_validator(mode='after')
    def valid_lineage(self):
        if (self.version == 1) != (self.parent_hash is None):
            raise ValueError('Only the initial policy has no parent')
        if len(set(self.rules)) != len(self.rules) or self.rules != tuple(sorted(self.rules)):
            raise ValueError('Rules must be unique and sorted')
        return self

    @property
    def digest(self):
        return content_hash(self.model_dump(mode='json'))

    @classmethod
    def initial(cls):
        return cls(version=1, parent_hash=None, rules=(), source='scripted')

    def revise(self, rules, *, source):
        return ResearchPolicy(version=self.version + 1, parent_hash=self.digest,
                              rules=tuple(sorted(rules)), source=source)


def write_content_addressed(directory, payload):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    digest = content_hash(payload)
    path = directory / (digest + '.json')
    data = json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + '\n'
    try:
        with path.open('x') as handle:
            handle.write(data)
    except FileExistsError:
        if json.loads(path.read_text()) != payload:
            raise ValueError('Refusing to overwrite a corrupt content-addressed artifact') from None
    return path


def write_policy(directory, policy):
    if not isinstance(policy, ResearchPolicy):
        raise TypeError('A validated frozen policy is required')
    return write_content_addressed(directory, policy.model_dump(mode='json'))


def read_policy(path):
    path = Path(path)
    payload = json.loads(path.read_text())
    if path.stem != content_hash(payload):
        raise ValueError('Policy content hash mismatch')
    return ResearchPolicy.model_validate_json(json.dumps(payload))
