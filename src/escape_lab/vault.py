"""Trusted, run-scoped experiment replay. Never give agents access to this database."""

import hashlib
import hmac
import json
import math
import secrets
import sqlite3
from contextlib import closing, contextmanager
from datetime import UTC, datetime
from pathlib import Path
from uuid import uuid4

PUBLIC_FIELDS = {
    "candidate_id",
    "gene",
    "protein_variant",
    "position",
    "wt",
    "mut",
    "aligned_position",
    "group_id",
    "split",
}
OPTIONAL_PUBLIC_FIELDS = {"context", "reference_sequence_sha256"}
SPLITS = {"development", "validation", "final"}


class VaultError(ValueError):
    def __init__(self, message, status_code=400):
        super().__init__(message)
        self.status_code = status_code


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _hash(value):
    return hashlib.sha256(value.encode()).hexdigest()


def _now():
    return datetime.now(UTC).isoformat()


def _text(value, name, limit=2000):
    if not isinstance(value, str) or not value.strip() or len(value) > limit:
        raise VaultError(f"Invalid {name}", 422)


class Vault:
    """Provision only from trusted code; HTTP clients get run-scoped capabilities.

    The database contains hidden outcomes. Its path, process credentials, and the
    admin token must remain outside the agent's filesystem and tool environment.
    SQLite transactions enforce reservation and replay accounting atomically.
    """

    def __init__(self, db_path, candidates=None, outcomes=None, *, allow_final=False):
        self.db_path = Path(db_path)
        self.allow_final = allow_final
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        if (candidates is None) != (outcomes is None):
            raise VaultError("Supply both candidates and outcomes")
        with closing(self._connect()) as connection, connection:
            connection.executescript("""
                PRAGMA journal_mode=WAL;
                CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS candidates (
                    candidate_id TEXT PRIMARY KEY, public_json TEXT NOT NULL,
                    outcome_json TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS runs (
                    run_id TEXT PRIMARY KEY, token_hash TEXT NOT NULL,
                    candidate_ids_json TEXT NOT NULL, budget INTEGER NOT NULL,
                    policy_version TEXT NOT NULL, split TEXT NOT NULL, seed INTEGER NOT NULL,
                    created_at TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS commits (
                    commit_id TEXT PRIMARY KEY, run_id TEXT NOT NULL REFERENCES runs(run_id),
                    candidate_id TEXT NOT NULL REFERENCES candidates(candidate_id),
                    predicted_probability REAL NOT NULL, reason TEXT NOT NULL,
                    policy_version TEXT NOT NULL, committed_at TEXT NOT NULL,
                    UNIQUE(run_id, candidate_id));
                CREATE TABLE IF NOT EXISTS reveals (
                    commit_id TEXT PRIMARY KEY REFERENCES commits(commit_id),
                    run_id TEXT NOT NULL REFERENCES runs(run_id), response_json TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS requests (
                    run_id TEXT NOT NULL REFERENCES runs(run_id), operation TEXT NOT NULL,
                    idempotency_key TEXT NOT NULL, payload_hash TEXT NOT NULL,
                    response_json TEXT NOT NULL,
                    PRIMARY KEY(run_id, operation, idempotency_key));
                CREATE TABLE IF NOT EXISTS events (
                    sequence INTEGER PRIMARY KEY AUTOINCREMENT,
                    run_id TEXT NOT NULL REFERENCES runs(run_id), event_type TEXT NOT NULL,
                    event_json TEXT NOT NULL, created_at TEXT NOT NULL);
            """)
            for table in (
                "metadata",
                "candidates",
                "runs",
                "commits",
                "reveals",
                "requests",
                "events",
            ):
                for action in ("UPDATE", "DELETE"):
                    connection.execute(
                        f"CREATE TRIGGER IF NOT EXISTS {table}_no_{action.lower()} "
                        f"BEFORE {action} ON {table} BEGIN "
                        "SELECT RAISE(ABORT, 'append-only vault'); END;"
                    )
        self.db_path.chmod(0o600)
        if candidates is not None:
            self._import(candidates, outcomes)
        else:
            with self._transaction() as connection:
                if (
                    connection.execute("SELECT 1 FROM metadata WHERE key='dataset_hash'").fetchone()
                    is None
                ):
                    raise VaultError("New vault requires an initial dataset")

    def _connect(self):
        connection = sqlite3.connect(self.db_path, timeout=30)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys=ON")
        return connection

    @contextmanager
    def _transaction(self):
        connection = self._connect()
        try:
            connection.execute("BEGIN IMMEDIATE")
            yield connection
            connection.commit()
        except BaseException:
            connection.rollback()
            raise
        finally:
            connection.close()

    def _import(self, candidates, outcomes):
        if not isinstance(candidates, list) or not isinstance(outcomes, dict):
            raise VaultError("Dataset must contain a candidate list and outcome mapping")
        ids = []
        for candidate in candidates:
            if (
                not isinstance(candidate, dict)
                or not PUBLIC_FIELDS <= set(candidate)
                or not set(candidate) <= PUBLIC_FIELDS | OPTIONAL_PUBLIC_FIELDS
            ):
                raise VaultError("Candidate public fields do not match the allowed schema")
            _text(candidate["candidate_id"], "candidate_id")
            if candidate["split"] not in SPLITS:
                raise VaultError("Invalid candidate split")
            ids.append(candidate["candidate_id"])
        if len(set(ids)) != len(ids) or set(ids) != set(outcomes):
            raise VaultError("Candidate and outcome identifiers must match uniquely")
        for outcome in outcomes.values():
            if not isinstance(outcome, dict):
                raise VaultError("Invalid outcome")
            label = outcome.get("dual_resistant")
            evidence = outcome.get("evidence_ids")
            if (
                label is not None
                and type(label) is not bool
                or not isinstance(outcome.get("records"), list)
                or not isinstance(evidence, list)
                or any(not isinstance(item, str) or not item for item in evidence)
                or label is True
                and not evidence
            ):
                raise VaultError("Outcome requires a boolean/null label, records, and evidence IDs")
        ordered = sorted(candidates, key=lambda candidate: candidate["candidate_id"])
        try:
            fingerprint = _hash(_json({"candidates": ordered, "outcomes": outcomes}))
        except (TypeError, ValueError) as error:
            raise VaultError("Dataset must contain finite JSON-compatible values") from error
        with self._transaction() as connection:
            old = connection.execute(
                "SELECT value FROM metadata WHERE key='dataset_hash'"
            ).fetchone()
            if old is not None:
                if old["value"] != fingerprint:
                    raise VaultError(
                        "Dataset fingerprint mismatch; existing vault is immutable", 409
                    )
                return
            connection.executemany(
                "INSERT INTO candidates VALUES (?, ?, ?)",
                [
                    (item["candidate_id"], _json(item), _json(outcomes[item["candidate_id"]]))
                    for item in ordered
                ],
            )
            connection.execute("INSERT INTO metadata VALUES ('dataset_hash', ?)", (fingerprint,))

    @staticmethod
    def _authorize(connection, run_id, run_token):
        row = connection.execute("SELECT * FROM runs WHERE run_id=?", (run_id,)).fetchone()
        if (
            row is None
            or not isinstance(run_token, str)
            or not hmac.compare_digest(row["token_hash"], _hash(run_token))
        ):
            raise VaultError("Unauthorized run access", 401)
        return row

    @staticmethod
    def _event(connection, run_id, event_type, payload):
        connection.execute(
            "INSERT INTO events(run_id,event_type,event_json,created_at) VALUES(?,?,?,?)",
            (run_id, event_type, _json(payload), _now()),
        )

    @staticmethod
    def _retry(connection, run_id, operation, key, payload):
        _text(key, "idempotency key", 200)
        row = connection.execute(
            "SELECT payload_hash,response_json FROM requests "
            "WHERE run_id=? AND operation=? AND idempotency_key=?",
            (run_id, operation, key),
        ).fetchone()
        if row is None:
            return None
        if row["payload_hash"] != _hash(_json(payload)):
            raise VaultError("Conflicting intent for idempotency key", 409)
        return json.loads(row["response_json"])

    @staticmethod
    def _remember(connection, run_id, operation, key, payload, response):
        connection.execute(
            "INSERT INTO requests VALUES(?,?,?,?,?)",
            (run_id, operation, key, _hash(_json(payload)), _json(response)),
        )

    def create_run(self, candidate_ids, budget, policy_version, split, seed=0):
        """Trusted administration only. All candidates must belong to the requested split."""
        if (
            not isinstance(candidate_ids, list)
            or not candidate_ids
            or any(not isinstance(item, str) for item in candidate_ids)
            or len(set(candidate_ids)) != len(candidate_ids)
        ):
            raise VaultError("Candidate pool must be a nonempty unique identifier list")
        if type(budget) is not int or budget < 1 or budget > len(candidate_ids):
            raise VaultError("Invalid budget")
        if split not in SPLITS or type(seed) is not int:
            raise VaultError("Invalid split or seed")
        if split == "final" and not self.allow_final:
            raise VaultError("Final split is sealed", 403)
        _text(policy_version, "policy version", 200)
        run_id, token = str(uuid4()), secrets.token_urlsafe(32)
        with self._transaction() as connection:
            for candidate_id in candidate_ids:
                row = connection.execute(
                    "SELECT public_json FROM candidates WHERE candidate_id=?", (candidate_id,)
                ).fetchone()
                if row is None or json.loads(row["public_json"])["split"] != split:
                    raise VaultError("Unknown candidate or candidate belongs to a different split")
            connection.execute(
                "INSERT INTO runs VALUES(?,?,?,?,?,?,?,?)",
                (
                    run_id,
                    _hash(token),
                    _json(candidate_ids),
                    budget,
                    policy_version,
                    split,
                    seed,
                    _now(),
                ),
            )
            self._event(
                connection,
                run_id,
                "run_created",
                {"budget": budget, "policy_version": policy_version, "split": split, "seed": seed},
            )
        return {"run_id": run_id, "run_token": token}

    def candidates(self, run_id, run_token):
        with self._transaction() as connection:
            run = self._authorize(connection, run_id, run_token)
            return [
                json.loads(
                    connection.execute(
                        "SELECT public_json FROM candidates WHERE candidate_id=?", (candidate_id,)
                    ).fetchone()["public_json"]
                )
                for candidate_id in json.loads(run["candidate_ids_json"])
            ]

    def commit(
        self,
        run_id,
        run_token,
        candidate_id,
        predicted_probability,
        reason,
        policy_version,
        idempotency_key,
    ):
        if (
            type(predicted_probability) not in (int, float)
            or not math.isfinite(predicted_probability)
            or not 0 <= predicted_probability <= 1
        ):
            raise VaultError(
                "Invalid predicted probability; expected a finite number in [0,1]", 422
            )
        _text(reason, "reason", 8000)
        payload = {
            "candidate_id": candidate_id,
            "predicted_probability": float(predicted_probability),
            "reason": reason,
            "policy_version": policy_version,
        }
        with self._transaction() as connection:
            run = self._authorize(connection, run_id, run_token)
            retry = self._retry(connection, run_id, "commit", idempotency_key, payload)
            if retry is not None:
                return retry
            if candidate_id not in json.loads(run["candidate_ids_json"]):
                raise VaultError("Candidate is outside this run's pool", 404)
            if policy_version != run["policy_version"]:
                raise VaultError("Commit policy does not match the frozen run policy", 409)
            if connection.execute(
                "SELECT 1 FROM commits WHERE run_id=? AND candidate_id=?", (run_id, candidate_id)
            ).fetchone():
                raise VaultError("Candidate already committed", 409)
            count = connection.execute(
                "SELECT COUNT(*) FROM commits WHERE run_id=?", (run_id,)
            ).fetchone()[0]
            if count >= run["budget"]:
                raise VaultError("Measurement budget exhausted or reserved by pending commits", 409)
            response = {"commit_id": str(uuid4()), **payload, "committed_at": _now()}
            connection.execute(
                "INSERT INTO commits VALUES(?,?,?,?,?,?,?)",
                (
                    response["commit_id"],
                    run_id,
                    candidate_id,
                    payload["predicted_probability"],
                    reason,
                    policy_version,
                    response["committed_at"],
                ),
            )
            self._remember(connection, run_id, "commit", idempotency_key, payload, response)
            self._event(connection, run_id, "committed", response)
            return response

    def reveal(self, run_id, run_token, commit_id, idempotency_key):
        payload = {"commit_id": commit_id}
        with self._transaction() as connection:
            self._authorize(connection, run_id, run_token)
            retry = self._retry(connection, run_id, "reveal", idempotency_key, payload)
            if retry is not None:
                return retry
            commit = connection.execute(
                "SELECT * FROM commits WHERE run_id=? AND commit_id=?", (run_id, commit_id)
            ).fetchone()
            if commit is None:
                raise VaultError("Commit not found in this run", 404)
            existing = connection.execute(
                "SELECT response_json FROM reveals WHERE commit_id=?", (commit_id,)
            ).fetchone()
            if existing:
                response = json.loads(existing["response_json"])
            else:
                outcome = connection.execute(
                    "SELECT outcome_json FROM candidates WHERE candidate_id=?",
                    (commit["candidate_id"],),
                ).fetchone()
                response = {
                    "commit_id": commit_id,
                    "candidate_id": commit["candidate_id"],
                    "outcome": json.loads(outcome["outcome_json"]),
                    "revealed_at": _now(),
                }
                connection.execute(
                    "INSERT INTO reveals VALUES(?,?,?)", (commit_id, run_id, _json(response))
                )
                self._event(connection, run_id, "revealed", response)
            self._remember(connection, run_id, "reveal", idempotency_key, payload, response)
            return response

    def state(self, run_id, run_token):
        with self._transaction() as connection:
            run = self._authorize(connection, run_id, run_token)
            commits = connection.execute(
                "SELECT commit_id FROM commits WHERE run_id=? ORDER BY rowid", (run_id,)
            ).fetchall()
            reveals = connection.execute(
                "SELECT response_json FROM reveals WHERE run_id=? ORDER BY rowid", (run_id,)
            ).fetchall()
            observations = [json.loads(row["response_json"]) for row in reveals]
            revealed = {item["commit_id"] for item in observations}
            return {
                "run_id": run_id,
                "budget": run["budget"],
                "split": run["split"],
                "seed": run["seed"],
                "policy_version": run["policy_version"],
                "candidate_count": len(json.loads(run["candidate_ids_json"])),
                "committed_count": len(commits),
                "revealed_count": len(reveals),
                "remaining_budget": run["budget"] - len(commits),
                "pending_commit_ids": [
                    row["commit_id"] for row in commits if row["commit_id"] not in revealed
                ],
                "observations": observations,
            }

    def events(self, run_id, run_token):
        with self._transaction() as connection:
            self._authorize(connection, run_id, run_token)
            rows = connection.execute(
                "SELECT * FROM events WHERE run_id=? ORDER BY sequence", (run_id,)
            ).fetchall()
            return [
                {
                    "sequence": row["sequence"],
                    "event_type": row["event_type"],
                    "payload": json.loads(row["event_json"]),
                    "created_at": row["created_at"],
                }
                for row in rows
            ]

    def verdict(self, run_id, run_token, status, evidence_ids):
        if status not in {"refuted", "unresolved"}:
            raise VaultError("Only refuted or unresolved verdicts are permitted", 422)
        if not isinstance(evidence_ids, list) or any(
            not isinstance(item, str) for item in evidence_ids
        ):
            raise VaultError("Evidence IDs must be a list of strings", 422)
        with self._transaction() as connection:
            self._authorize(connection, run_id, run_token)
            observations = [
                json.loads(row["response_json"])
                for row in connection.execute(
                    "SELECT response_json FROM reveals WHERE run_id=?", (run_id,)
                )
            ]
            positives = [item for item in observations if item["outcome"]["dual_resistant"] is True]
            allowed = {
                evidence for item in positives for evidence in item["outcome"]["evidence_ids"]
            }
            if status == "refuted" and (not evidence_ids or not set(evidence_ids) <= allowed):
                raise VaultError(
                    "Refutation requires evidence from revealed positive measurements", 422
                )
            if status == "unresolved" and (evidence_ids or positives):
                raise VaultError(
                    "Unresolved requires no positive counterexamples or claimed evidence", 422
                )
            result = {
                "status": status,
                "evidence_ids": sorted(set(evidence_ids)),
                "verified_candidate_ids": sorted(item["candidate_id"] for item in positives),
                "scope": "Retrospective separate-drug cell-line assay; no patient, combination, "
                "sequence, or universal coverage conclusion.",
            }
            self._event(connection, run_id, "verdict", result)
            return result
