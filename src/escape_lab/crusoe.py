"""Small, explicit Crusoe adapter. Importing this module never performs inference.

The caller supplies credentials; this module never reads .env or logs request bodies.
Transport injection permits deterministic tests without network access or payment.
"""

from __future__ import annotations

import json
import math
import time
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, Self

import httpx
from pydantic import BaseModel, ConfigDict, ValidationError

DEFAULT_BASE_URL = "https://api.inference.crusoecloud.com/v1"


class CrusoeError(RuntimeError):
    """Sanitized operational failure, without provider bodies or credentials."""


class ModelInfo(BaseModel):
    model_config = ConfigDict(extra="ignore")
    id: str
    supported_parameters: list[str] | dict[str, Any] | None = None

    @property
    def tools_supported(self) -> bool | None:
        if self.supported_parameters is None:
            return None
        return "tools" in self.supported_parameters


class ReadyArguments(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)
    ready: bool


@dataclass(frozen=True)
class SmokeResult:
    passed: bool
    model_id: str
    detail: str


@dataclass(frozen=True)
class ChatResult:
    message: dict[str, Any]
    finish_reason: str | None
    usage: dict[str, int | None]


class CrusoeClient:
    """Bounded sync HTTP adapter with explicit discovery and smoke verification.

    max_total_tokens is a stop threshold checked before a new completion. A final
    request may exceed that threshold because prompt tokens are provider-counted.
    max_tokens separately bounds each response. Missing usage is charged using a
    conservative estimate and marked as estimated; it is never called zero cost.
    """

    def __init__(
        self,
        api_key: str,
        *,
        model_id: str | None = None,
        base_url: str = DEFAULT_BASE_URL,
        transport: httpx.BaseTransport | None = None,
        max_tokens: int = 1024,
        max_retries: int = 2,
        max_total_tokens: int = 100_000,
        timeout: float = 60,
        sleep: Callable[[float], None] = time.sleep,
        parallel_tool_calls: bool | None = None,
    ):
        if not api_key.strip():
            raise ValueError("A nonempty Crusoe API key is required.")
        if not 16 <= max_tokens <= 8192:
            raise ValueError("max_tokens must be between 16 and 8192.")
        if not 0 <= max_retries <= 3 or max_total_tokens <= 0:
            raise ValueError("Invalid retry or token budget.")
        if parallel_tool_calls is not None and type(parallel_tool_calls) is not bool:
            raise ValueError("parallel_tool_calls must be boolean or omitted")
        self.parallel_tool_calls = parallel_tool_calls
        self.model_id = model_id
        self.max_tokens = max_tokens
        self.max_retries = max_retries
        self.max_total_tokens = max_total_tokens
        self.total_tokens = 0
        self.usage_records: list[dict[str, Any]] = []
        self.models: dict[str, ModelInfo] = {}
        self._smoke_verified: set[str] = set()
        self._sleep = sleep
        self._http = httpx.Client(
            base_url=base_url.rstrip("/") + "/",
            headers={"Authorization": f"Bearer {api_key}", "Content-Type": "application/json"},
            transport=transport,
            timeout=timeout,
            follow_redirects=False,
        )

    def close(self) -> None:
        self._http.close()

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def _request(self, method: str, path: str, payload: dict | None = None) -> dict:
        for attempt in range(self.max_retries + 1):
            try:
                response = self._http.request(method, path, json=payload)
            except httpx.HTTPError as error:
                # Ambiguous transport failure may follow a charged request. Do not retry it.
                raise CrusoeError(f"Crusoe transport failed ({type(error).__name__}).") from None
            if response.status_code in (429, 503) and attempt < self.max_retries:
                try:
                    delay = float(response.headers.get("retry-after", 2 ** attempt))
                except ValueError:
                    delay = 2 ** attempt
                if not math.isfinite(delay):
                    delay = 2 ** attempt
                self._sleep(min(max(delay, 0), 10))
                continue
            if response.status_code >= 300:
                raise CrusoeError(f"Crusoe returned HTTP {response.status_code}.")
            try:
                result = response.json()
            except (ValueError, json.JSONDecodeError):
                raise CrusoeError("Crusoe returned invalid JSON.") from None
            if not isinstance(result, dict):
                raise CrusoeError("Crusoe returned an invalid response object.")
            return result
        raise CrusoeError("Crusoe retry budget exhausted.")

    def discover_models(self, *, nemotron_only: bool = True) -> list[ModelInfo]:
        payload = self._request("GET", "models")
        try:
            models = [ModelInfo.model_validate(row) for row in payload["data"]]
        except (KeyError, TypeError, ValidationError):
            raise CrusoeError("Crusoe model discovery returned an invalid model list.") from None
        self.models = {model.id: model for model in models}
        if nemotron_only:
            return [model for model in models if "nemotron" in model.id.lower()]
        return models

    def configure_model(self, exact_id: str) -> ModelInfo:
        if exact_id not in self.models:
            raise CrusoeError("Requested model ID was not present in model discovery.")
        self.model_id = exact_id
        return self.models[exact_id]

    def _require_model(self, *, smoke: bool = False) -> ModelInfo:
        if not self.model_id or self.model_id not in self.models:
            raise CrusoeError("Discover models and configure an exact served model ID first.")
        model = self.models[self.model_id]
        if model.tools_supported is False:
            raise CrusoeError("Selected model explicitly does not advertise tools support.")
        if model.tools_supported is None and model.id not in self._smoke_verified and not smoke:
            raise CrusoeError("Tool capability is unknown; call run_smoke() explicitly first.")
        if self.total_tokens >= self.max_total_tokens:
            raise CrusoeError("Crusoe total token budget exhausted.")
        return model

    def chat(
        self,
        messages: list[dict],
        tools: list[dict],
        *,
        prompt_version: str = "unspecified",
        _smoke: bool = False,
    ) -> ChatResult:
        model = self._require_model(smoke=_smoke)
        payload = {
            "model": model.id,
            "messages": messages,
            "tools": tools,
            "tool_choice": "required",
            "max_tokens": self.max_tokens,
            "stream": False,
        }
        if self.parallel_tool_calls is not None:
            payload["parallel_tool_calls"] = self.parallel_tool_calls
        try:
            result = self._request("POST", "chat/completions", payload)
        except CrusoeError:
            # A transport/parser failure can follow charged inference. Conservatively
            # reserve budget while keeping unknown usage distinct from provider usage.
            self._record_usage(payload, None, prompt_version, "failed_or_unknown")
            raise
        usage = self._record_usage(payload, result.get("usage"), prompt_version, "completed")
        try:
            choice = result["choices"][0]
            raw_message = choice["message"]
            if not isinstance(raw_message, dict):
                raise TypeError
            # Retain only public final content and tool calls, never private reasoning fields.
            message = {"role": "assistant", "content": raw_message.get("content")}
            if isinstance(raw_message.get("tool_calls"), list):
                message["tool_calls"] = raw_message["tool_calls"]
            return ChatResult(message, choice.get("finish_reason"), usage)
        except (KeyError, IndexError, TypeError):
            raise CrusoeError("Crusoe completion did not contain a valid message.") from None

    def _record_usage(self, payload, raw_usage, prompt_version, request_status):
        if not isinstance(raw_usage, dict):
            raw_usage = {}
        usage = {}
        for field in ("prompt_tokens", "completion_tokens", "total_tokens"):
            value = raw_usage.get(field)
            usage[field] = value if type(value) is int and value >= 0 else None
        estimated = usage["total_tokens"] is None
        charge = usage["total_tokens"]
        if charge is None:
            charge = len(json.dumps(payload)) + self.max_tokens
        self.total_tokens += charge
        self.usage_records.append({
            "model_id": payload["model"],
            "prompt_version": prompt_version,
            **usage,
            "budget_charge_tokens": charge,
            "usage_estimated": estimated,
            "request_status": request_status,
        })
        return usage

    def run_smoke(self) -> SmokeResult:
        """Explicitly perform one paid, bounded tool-call capability check."""
        tool = {"type": "function", "function": {
            "name": "report_ready", "description": "Report readiness using ready=true.",
            "parameters": ReadyArguments.model_json_schema(),
        }}
        result = self.chat(
            [{"role": "user", "content": "Call report_ready with ready=true."}],
            [tool], prompt_version="smoke-v1", _smoke=True,
        )
        calls = result.message.get("tool_calls", [])
        try:
            if len(calls) != 1 or calls[0]["function"]["name"] != "report_ready":
                raise ValueError
            arguments = ReadyArguments.model_validate_json(calls[0]["function"]["arguments"])
            if not arguments.ready:
                raise ValueError
        except (KeyError, TypeError, ValueError, ValidationError):
            return SmokeResult(False, self.model_id or "", "A valid readiness tool call was absent.")
        self._smoke_verified.add(self.model_id)
        return SmokeResult(True, self.model_id, "Structured tool call validated locally.")
