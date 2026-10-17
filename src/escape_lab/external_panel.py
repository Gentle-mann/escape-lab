"""One-shot evaluation of a separately generated biochemical FGFR panel.

Metadata inspection reads header and identifier cells only. Outcome cells are
read only after a frozen protocol, complete committed predictions, and an
exclusive receipt have been verified. This is a trusted evaluator boundary,
not a filesystem sandbox. The original MAVE final receipt is never accessed.
"""

import hashlib
import json
import math
import os
import re
import uuid
from collections import Counter
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import openpyxl
from scipy.stats import spearmanr

DEFAULT_RECEIPT = Path(__file__).resolve().parents[2] / "data/revision-v2/external-receipt.json"
DRUGS = ("pem", "fut")
_NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
_METHOD = re.compile(r"[a-z][a-z0-9_-]{0,63}\Z")


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _read_json(path):
    return json.loads(Path(path).read_text())


def _write_json(path, payload):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp-" + uuid.uuid4().hex)
    try:
        temporary.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _now():
    return datetime.now(UTC).isoformat()


def _manifest(manifest_path):
    manifest = _read_json(manifest_path)
    if manifest.get("schema_version") != "s13-external-panel-v1":
        raise ValueError("Unsupported external panel manifest")
    document = Path(manifest["document_path"])
    if not document.is_absolute():
        raise ValueError("The panel document path must be absolute")
    variants = manifest["variants"]
    if not variants or manifest.get("variant_count") != len(variants):
        raise ValueError("Panel variant count does not match its metadata")
    ids, columns = set(), set()
    for variant in variants:
        identity, column = variant["variant_id"], variant["column"]
        if not re.fullmatch(r"FGFR[1-4]:[A-Z][1-9]\d*[A-Z]", identity):
            raise ValueError("Invalid external variant identity")
        if not re.fullmatch(r"[A-Z]{1,3}", column) or column == "A":
            raise ValueError("Invalid external variant column")
        if identity in ids or column in columns:
            raise ValueError("Duplicate external variant or column")
        if not isinstance(variant["source_label"], str) or not variant["source_label"]:
            raise ValueError("Missing source variant label")
        if not isinstance(variant["group_id"], str) or not variant["group_id"]:
            raise ValueError("Missing aligned group identity")
        ids.add(identity)
        columns.add(column)
    header = manifest["header_row"]
    if type(header) is not int or header < 1:
        raise ValueError("Invalid panel header row")
    if set(manifest["drug_rows"]) != set(DRUGS):
        raise ValueError("Both frozen drug rows are required")
    rows = list(manifest["drug_rows"].values())
    if any(type(row) is not int or row <= header for row in rows) or len(set(rows)) != 2:
        raise ValueError("Invalid drug row metadata")
    if manifest["drug_labels"] != {"pem": "Pemigatinib", "fut": "Futibatinib"}:
        raise ValueError("Drug labels do not identify the expected compounds")
    return manifest


def inspect_external_panel(manifest_path):
    """Validate source headers and identifiers without accessing outcome cells."""
    manifest = _manifest(manifest_path)
    if _sha256(manifest["document_path"]) != manifest["document_sha256"]:
        raise ValueError("External document checksum mismatch")
    book = openpyxl.load_workbook(manifest["document_path"], read_only=True, data_only=False)
    try:
        sheet = book[manifest["sheet"]]
        for drug in DRUGS:
            actual = sheet[f"A{manifest['drug_rows'][drug]}"].value
            if actual != manifest["drug_labels"][drug]:
                raise ValueError("External drug header mismatch")
        for variant in manifest["variants"]:
            actual = sheet[f"{variant['column']}{manifest['header_row']}"].value
            if actual != variant["source_label"]:
                raise ValueError("External variant header mismatch")
    finally:
        book.close()
    return {"panel_id": manifest["panel_id"], "variant_count": len(manifest["variants"]),
            "variant_ids": [v["variant_id"] for v in manifest["variants"]],
            "group_ids": sorted({v["group_id"] for v in manifest["variants"]}),
            "drug_rows": manifest["drug_rows"], "numeric_outcomes_inspected": False}


def parse_residual_activity(value):
    """Parse a measurement without substituting missing or censored values.

Finite measurements outside the nominal 0–100 range are retained and flagged:
normalization noise must not silently be clipped or converted into missing data.
"""
    if value is None:
        return {"status": "missing", "value": None}
    if isinstance(value, bool):
        return {"status": "invalid", "value": None}
    if isinstance(value, (int, float)):
        number = float(value)
        if math.isnan(number):
            return {"status": "missing", "value": None}
        if not math.isfinite(number):
            return {"status": "invalid", "value": None}
    elif isinstance(value, str):
        text = value.strip().replace("−", "-")
        if text.casefold() in {"", "na", "n/a", "nan", "nd", "n.d.", "not determined", "—", "-"}:
            return {"status": "missing", "value": None}
        text = re.sub(r"\s*[×x]\s*10\s*\^\s*([+-]?\d+)\s*", r"e\1", text)
        censored = re.fullmatch(rf"(<=|>=|[<>≤≥])\s*({_NUMBER})\s*%?", text)
        if censored:
            bound = float(censored.group(2))
            if not math.isfinite(bound):
                return {"status": "invalid", "value": None}
            status = "censored_lower_bound" if censored.group(1) in {">", "≥", ">="} else "censored_upper_bound"
            return {"status": status, "value": None, "bound": bound,
                    "inclusive_bound": censored.group(1) in {"≤", "≥", "<=", ">="}}
        numeric = re.fullmatch(rf"({_NUMBER})\s*%?", text)
        if numeric is None:
            return {"status": "invalid", "value": None}
        number = float(numeric.group(1))
        if not math.isfinite(number):
            return {"status": "invalid", "value": None}
    else:
        return {"status": "invalid", "value": None}
    return {"status": "observed", "value": number,
            "outside_nominal_range": not 0 <= number <= 100}


def _predictions(path, methods, manifest):
    payload = _read_json(path)
    if set(payload) != {"predictions"} or not isinstance(payload["predictions"], list):
        raise ValueError("Invalid committed prediction document")
    variants = {v["variant_id"] for v in manifest["variants"]}
    expected = {(method, variant, drug) for method in methods for variant in variants for drug in DRUGS}
    parsed = {}
    for row in payload["predictions"]:
        if not isinstance(row, dict) or set(row) != {"method", "variant_id", "drug", "score"}:
            raise ValueError("Invalid committed prediction fields")
        key = (row["method"], row["variant_id"], row["drug"])
        if key not in expected:
            raise ValueError("Unknown prediction method, variant, or drug")
        if key in parsed:
            raise ValueError("Duplicate committed prediction")
        score = row["score"]
        if isinstance(score, bool) or not isinstance(score, (int, float)):
            raise TypeError("Prediction score must be numeric")
        if not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError("Prediction score must be a finite probability")
        parsed[key] = float(score)
    if set(parsed) != expected:
        raise ValueError("Committed prediction coverage is incomplete")
    return parsed


@dataclass(frozen=True)
class _Authorization:
    receipt_path: Path
    identifier: str
    files: dict


def _verify_authorization(authorization):
    receipt = _read_json(authorization.receipt_path)
    if receipt.get("status") != "reserved" or receipt.get("authorization_id") != authorization.identifier:
        raise ValueError("External outcome read is not authorized by an active receipt")
    if receipt.get("bound_files") != authorization.files:
        raise ValueError("External receipt commitment mismatch")
    for name, record in authorization.files.items():
        if _sha256(record["path"]) != record["sha256"]:
            raise ValueError(f"Frozen {name} changed before external outcome access")


def _read_numeric_cells(manifest, authorization):
    _verify_authorization(authorization)
    book = openpyxl.load_workbook(manifest["document_path"], read_only=True, data_only=False)
    measurements = {}
    try:
        sheet = book[manifest["sheet"]]
        for variant in manifest["variants"]:
            for drug in DRUGS:
                coordinate = f"{variant['column']}{manifest['drug_rows'][drug]}"
                cell = sheet[coordinate]
                if cell.data_type == "f":
                    parsed = {"status": "formula_not_evaluated", "value": None}
                elif cell.data_type == "e":
                    parsed = {"status": "spreadsheet_error", "value": None}
                else:
                    parsed = parse_residual_activity(cell.value)
                measurements[(variant["variant_id"], drug)] = {**parsed, "source_cell": coordinate}
    finally:
        book.close()
    _verify_authorization(authorization)
    return measurements


def _summaries(methods, predictions, measurements, manifest):
    coverage, results, records = {}, {}, []
    for drug in DRUGS:
        selected = [measurements[(v["variant_id"], drug)] for v in manifest["variants"]]
        counts = Counter(row["status"] for row in selected)
        coverage[drug] = {"expected": len(selected), "observed": counts.get("observed", 0),
                          **dict(counts), "outside_nominal_range": sum(
                              bool(row.get("outside_nominal_range")) for row in selected)}
    for method in methods:
        per_drug = {}
        for drug in DRUGS:
            truth, scores, groups = [], [], set()
            for variant in manifest["variants"]:
                value = measurements[(variant["variant_id"], drug)]
                score = predictions[(method, variant["variant_id"], drug)]
                record = {"method": method, "variant_id": variant["variant_id"], "drug": drug,
                          "group_id": variant["group_id"], "prediction_score": score,
                          "residual_activity": value["value"], "measurement_status": value["status"],
                          "source_cell": value["source_cell"]}
                for field in ("bound", "inclusive_bound", "outside_nominal_range"):
                    if field in value:
                        record[field] = value[field]
                records.append(record)
                if value["status"] == "observed":
                    truth.append(value["value"])
                    scores.append(score)
                    groups.add(variant["group_id"])
            correlation, reason = None, None
            if len(truth) < 5:
                reason = "fewer_than_five_observed_pairs"
            elif len(set(truth)) < 2:
                reason = "constant_measurements"
            elif len(set(scores)) < 2:
                reason = "constant_predictions"
            else:
                correlation = float(spearmanr(scores, truth).statistic)
                if not math.isfinite(correlation):
                    correlation, reason = None, "undefined_spearman"
            per_drug[drug] = {"spearman": correlation, "n": len(truth),
                              "aligned_groups": len(groups), "reason": reason}
        values = [per_drug[drug]["spearman"] for drug in DRUGS]
        macro = sum(values) / 2 if all(v is not None for v in values) else None
        results[method] = {"macro_spearman": macro, "per_drug": per_drug}
    return coverage, results, records


def evaluate_external_panel(manifest_path, frozen_protocol_path, committed_predictions_path,
                            receipt_path=DEFAULT_RECEIPT, output_dir=None):
    """Reveal once under an explicit hash-bound contract; return aggregate results.

The protocol's ``external_evaluation`` object must contain ``authorized: true``,
``panel_manifest_sha256``, ``document_sha256``, ``predictions_sha256``,
``evaluator_sha256``, a unique list of safe ``methods``,
``primary_metric: mean_per_drug_spearman``, and ``minimum_n_per_drug: 5``.
Any exception after exclusive reservation permanently spends the receipt.
"""
    receipt_path = Path(receipt_path).resolve()
    if receipt_path.exists():
        raise FileExistsError("External panel receipt already exists; this test is spent")
    manifest_path = Path(manifest_path).resolve()
    frozen_protocol_path = Path(frozen_protocol_path).resolve()
    committed_predictions_path = Path(committed_predictions_path).resolve()
    manifest = _manifest(manifest_path)
    contract = _read_json(frozen_protocol_path).get("external_evaluation", {})
    if contract.get("authorized") is not True:
        raise ValueError("Frozen protocol does not authorize external evaluation")
    if (contract.get("primary_metric") != "mean_per_drug_spearman"
            or type(contract.get("minimum_n_per_drug")) is not int
            or contract["minimum_n_per_drug"] != 5):
        raise ValueError("Unsupported external evaluation metric contract")
    methods = contract.get("methods")
    if (not isinstance(methods, list) or not methods
            or any(not isinstance(m, str) or not _METHOD.fullmatch(m) for m in methods)
            or len(set(methods)) != len(methods)):
        raise ValueError("Invalid frozen external method list")
    bindings = {"panel_manifest": manifest_path, "document": Path(manifest["document_path"]),
                "predictions": committed_predictions_path, "evaluator": Path(__file__).resolve()}
    files = {}
    for name, path in bindings.items():
        actual = _sha256(path)
        if contract.get(name + "_sha256") != actual:
            raise ValueError(f"Frozen {name} checksum mismatch")
        files[name] = {"path": str(path), "sha256": actual}
    if manifest["document_sha256"] != files["document"]["sha256"]:
        raise ValueError("Manifest document checksum mismatch")
    files["frozen_protocol"] = {"path": str(frozen_protocol_path),
                                "sha256": _sha256(frozen_protocol_path)}
    predictions = _predictions(committed_predictions_path, methods, manifest)
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt = {"status": "reserved", "authorization_id": uuid.uuid4().hex,
               "reserved_at": _now(), "panel_id": manifest["panel_id"], "bound_files": files,
               "predictions_sha256": files["predictions"]["sha256"]}
    with receipt_path.open("x") as handle:
        json.dump(receipt, handle, indent=2)
    authorization = _Authorization(receipt_path, receipt["authorization_id"], files)
    try:
        inspect_external_panel(manifest_path)
        _verify_authorization(authorization)
        output = Path(output_dir).resolve() if output_dir else receipt_path.parent / "external-results"
        output.mkdir(parents=True, exist_ok=True)
        records_path, report_path = output / "external-records.json", output / "external-report.json"
        if records_path.exists() or report_path.exists():
            raise FileExistsError("External output artifacts already exist")
        measurements = _read_numeric_cells(manifest, authorization)
        coverage, methods_result, records = _summaries(methods, predictions, measurements, manifest)
        _write_json(records_path, {"panel_id": manifest["panel_id"], "records": records})
        result = {"status": "completed", "panel_id": manifest["panel_id"],
                  "primary_metric": "mean_per_drug_spearman", "minimum_n_per_drug": 5,
                  "measurement_coverage": coverage, "methods": methods_result,
                  "primary_evaluable": all(m["macro_spearman"] is not None for m in methods_result.values()),
                  "wt_adjustment_applied": False, "classification_metrics_computed": False,
                  "records_path": str(records_path), "records_sha256": _sha256(records_path),
                  "interpretation": "Descriptive cross-assay biochemical transfer on related FGFR variants; "
                                    "not cellular-resistance accuracy or independent clinical validation."}
        _write_json(report_path, result)
        receipt.update(status="completed", completed_at=_now(), report_path=str(report_path),
                       report_sha256=_sha256(report_path), records_sha256=result["records_sha256"])
        _write_json(receipt_path, receipt)
        return result
    except BaseException as error:
        receipt.update(status="failed_spent", failed_at=_now(), error_type=type(error).__name__)
        _write_json(receipt_path, receipt)
        raise
