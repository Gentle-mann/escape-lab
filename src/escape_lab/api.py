"""Small HTTP boundary exposing run capabilities, never the hidden dataset."""

import hmac
from typing import Annotated, Literal

from fastapi import FastAPI, Header, HTTPException
from fastapi.responses import JSONResponse
from pydantic import BaseModel, ConfigDict, Field

from escape_lab.vault import Vault, VaultError


class Request(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)


class RunRequest(Request):
    candidate_ids: list[str]
    budget: int = Field(ge=1)
    policy_version: str = Field(min_length=1, max_length=200)
    split: Literal["development", "validation", "final"]
    seed: int = 0


class CommitRequest(Request):
    candidate_id: str
    predicted_probability: float = Field(ge=0, le=1, allow_inf_nan=False)
    reason: str = Field(min_length=1, max_length=8000)
    policy_version: str
    idempotency_key: str = Field(min_length=1, max_length=200)


class RevealRequest(Request):
    commit_id: str
    idempotency_key: str = Field(min_length=1, max_length=200)


class VerdictRequest(Request):
    status: Literal["refuted", "unresolved"]
    evidence_ids: list[str]


RunToken = Annotated[str | None, Header(alias="X-Run-Token")]
AdminToken = Annotated[str | None, Header(alias="X-Admin-Token")]


def create_app(vault: Vault, admin_token: str) -> FastAPI:
    if not isinstance(admin_token, str) or not admin_token.strip():
        raise ValueError("A nonempty admin token is required")
    app = FastAPI(title="Escape Lab Experiment Vault", docs_url=None, redoc_url=None)

    @app.exception_handler(VaultError)
    async def handle_vault_error(request, error):
        return JSONResponse(status_code=error.status_code, content={"detail": str(error)})

    @app.post("/runs")
    def create_run(body: RunRequest, token: AdminToken = None):
        if token is None or not hmac.compare_digest(token.encode(), admin_token.encode()):
            raise HTTPException(status_code=401, detail="Unauthorized administration")
        return vault.create_run(**body.model_dump())

    @app.get("/runs/{run_id}/candidates")
    def candidates(run_id: str, token: RunToken = None):
        return vault.candidates(run_id, token)

    @app.get("/runs/{run_id}/state")
    def state(run_id: str, token: RunToken = None):
        return vault.state(run_id, token)

    @app.get("/runs/{run_id}/events")
    def events(run_id: str, token: RunToken = None):
        return vault.events(run_id, token)

    @app.post("/runs/{run_id}/commits")
    def commit(run_id: str, body: CommitRequest, token: RunToken = None):
        return vault.commit(run_id, token, **body.model_dump())

    @app.post("/runs/{run_id}/reveals")
    def reveal(run_id: str, body: RevealRequest, token: RunToken = None):
        return vault.reveal(run_id, token, **body.model_dump())

    @app.post("/runs/{run_id}/verdict")
    def verdict(run_id: str, body: VerdictRequest, token: RunToken = None):
        return vault.verdict(run_id, token, **body.model_dump())

    return app
