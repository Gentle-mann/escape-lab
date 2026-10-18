"""Rebuild the frozen development experiment from public sources, without credentials."""

import argparse
import csv
import hashlib
import importlib.metadata
import io
import json
import math
import platform
import time
import zipfile
from datetime import UTC, datetime
from pathlib import Path

import httpx
import numpy as np
import openpyxl
from scipy.special import expit
from sklearn.metrics import average_precision_score

from escape_lab.activation_data import load_development_activation
from escape_lab.activation_features import ActivationFeatureAugmenter
from escape_lab.biological_features import (
    SCORE_SETS,
    BiologicalFeatureStore,
    build_feature_store,
)
from escape_lab.data import sha256
from escape_lab.fgfr2_sequence_features import extract_features, load_features
from escape_lab.fgfr2_sequence_models import FGFR2SequenceModel, activation_label
from escape_lab.fgfr2_sequence_run import fold_partitions, summarize

RECIPE_PATH = Path(__file__).parent / "reproduction_assets/fgfr2-headline-v1.json"
CANDIDATE_FIELDS = (
    "candidate_id", "gene", "protein_variant", "position", "wt", "mut", "aligned_position",
    "group_id", "split", "context", "reference_sequence_sha256",
)
PINNED_PACKAGES = {
    "torch": "2.8.0", "transformers": "4.57.1", "huggingface-hub": "0.36.2",
    "safetensors": "0.8.0",
}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _write_json(path, value):
    with Path(path).open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def validate_recipe(recipe):
    candidates = recipe["candidates"]
    if _digest(candidates) != recipe["candidate_identity_sha256"]:
        raise ValueError("Frozen candidate identity digest mismatch")
    if (recipe["schema_version"] != "fgfr2-headline-reproduction-v1"
            or len(candidates) != 668
            or any(set(candidate) != set(CANDIDATE_FIELDS) for candidate in candidates)
            or any(candidate["gene"] != "FGFR2"
                   or candidate["split"] not in {"development", "validation"}
                   for candidate in candidates)):
        raise ValueError("Recipe must contain exactly the approved development identities")
    groups = {candidate["group_id"] for candidate in candidates}
    if len(groups) != 10 or groups & set(recipe["reserved_groups"]):
        raise ValueError("Frozen cohort overlaps reserved groups or has changed")
    fold_partitions(candidates, recipe["fold_groups"])
    return recipe


def load_recipe():
    return validate_recipe(json.loads(RECIPE_PATH.read_text()))


def write_verified(path, payload, expected_sha256):
    """Validate before writing; no partial or silently replaced source inputs."""
    if hashlib.sha256(payload).hexdigest() != expected_sha256:
        raise ValueError(f"Source checksum mismatch: {Path(path).name}")
    with Path(path).open("xb") as handle:
        handle.write(payload)


def fetch_bytes(url, *, max_bytes=32_000_000):
    if not url.startswith("https://"):
        raise ValueError("Public input URLs must use HTTPS")
    chunks, size = [], 0
    with httpx.stream("GET", url, timeout=120, follow_redirects=True) as response:
        response.raise_for_status()
        for chunk in response.iter_bytes():
            size += len(chunk)
            if size > max_bytes:
                raise ValueError("Public input exceeds bounded download size")
            chunks.append(chunk)
    return b"".join(chunks)


def supplement_bytes(payload, member, expected_sha256):
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        matches = [name for name in archive.namelist() if Path(name).name == member]
        if len(matches) != 1 or archive.getinfo(matches[0]).file_size > 8_000_000:
            raise ValueError("Expected exactly one bounded source workbook in supplement archive")
        workbook = archive.read(matches[0])
    if hashlib.sha256(workbook).hexdigest() != expected_sha256:
        raise ValueError("Supplement workbook checksum mismatch")
    return workbook


def predictor_rows(workbook_path, candidates):
    """Read identities and pretrained scores only; never use assay target columns here."""
    if not candidates or any(c["split"] not in {"development", "validation"} for c in candidates):
        raise ValueError("Predictor export permits development candidates only")
    allowed = {(c["gene"], f"{c['wt']}{c['position']}{c['mut']}"): c for c in candidates}
    if len(allowed) != len(candidates):
        raise ValueError("Duplicate predictor identity")
    values = {c["candidate_id"]: {"alphamissense": set(), "eve": set()} for c in candidates}
    workbook = openpyxl.load_workbook(workbook_path, read_only=True, data_only=True)
    try:
        sheet = workbook["S10_Screen data table"]
        header = next(sheet.iter_rows(min_row=10, max_row=10, max_col=43, values_only=True))
        if (header[0], header[4], header[41], header[42]) != (
                "FGFR", "AA MUT", "AlphaMissense score", "EVE score"):
            raise ValueError("Unexpected predictor source columns")
        identities = sheet.iter_rows(min_row=11, max_col=5, values_only=True)
        scores = sheet.iter_rows(min_row=11, min_col=42, max_col=43, values_only=True)
        for identity, scores_row in zip(identities, scores, strict=True):
            candidate = allowed.get((identity[0], str(identity[4])))
            if candidate is None:
                continue
            for name, raw in zip(("alphamissense", "eve"), scores_row, strict=True):
                if raw is None or str(raw).strip().lower() in {"", ".", "na", "n/a", "nan"}:
                    continue
                value = float(raw)
                if not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError("Predictor scores require finite values between zero and one")
                values[candidate["candidate_id"]][name].add(value)
    finally:
        workbook.close()
    result = []
    for candidate in candidates:
        row = {"candidate_id_name_match": candidate["candidate_id"], "gene": candidate["gene"],
               "protein_variant": f"{candidate['wt']}{candidate['position']}{candidate['mut']}"}
        for name, distinct in values[candidate["candidate_id"]].items():
            row[name] = next(iter(distinct)) if len(distinct) == 1 else ""
            row[f"{name}_status"] = ("available" if len(distinct) == 1 else
                                    "missing" if not distinct else "discordant_not_collapsed")
        result.append(row)
    return result


def prepare_inputs(output, *, workbook_path=None):
    """Fetch public references and reconstruct local development inputs without fitting."""
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    recipe, started = load_recipe(), time.monotonic()
    source_records = []
    try:
        raw, structures = output / "metadata", output / "structures"
        raw.mkdir()
        structures.mkdir()
        for gene, drugs in SCORE_SETS.items():
            for drug, accession in drugs.items():
                url = f"https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00001256-a-{accession}"
                payload = fetch_bytes(url, max_bytes=2_000_000)
                path = raw / f"{gene}-{drug}.json"
                # API metadata can acquire annotations. Exact assay sequences are verified below.
                write_verified(path, payload, hashlib.sha256(payload).hexdigest())
                source_records.append({"url": url, "sha256": sha256(path),
                                       "license": "CC BY-SA 4.0"})
        for filename, expected in recipe["biology"]["structure_sha256"].items():
            url = f"https://files.rcsb.org/download/{filename}"
            write_verified(structures / filename, fetch_bytes(url, max_bytes=4_000_000), expected)
            source_records.append({"url": url, "sha256": expected, "license": "CC0 1.0"})
        supplement = recipe["supplement"]
        if workbook_path is None:
            payload = supplement_bytes(fetch_bytes(supplement["url"]), supplement["archive_member"],
                                       supplement["sha256"])
            source_records.append({**supplement, "obtained": "downloaded this run"})
        else:
            payload = Path(workbook_path).read_bytes()
            source_records.append({**supplement, "obtained": "checksum-verified local workbook"})
        workbook = output / "supplement.xlsx"
        write_verified(workbook, payload, supplement["sha256"])
        workbook.chmod(0o600)
        manifest = build_feature_store(raw, structures, output / "biology")
        for key in ("feature_cache_sha256", "references_sha256"):
            if manifest[key] != recipe["biology"][key]:
                raise ValueError(f"Rebuilt biological features differ from frozen input: {key}")
        references = json.loads((output / "biology/references.json").read_text())
        if references["sequence_sha256"] != recipe["sequence_sha256"]:
            raise ValueError("Public assay sequences differ from the frozen reference")
        rows = predictor_rows(workbook, recipe["candidates"])
        with (output / "predictor-scores.csv").open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        files = {str(path.relative_to(output)): sha256(path) for path in sorted(output.rglob("*"))
                 if path.is_file()}
        status = {"status": "prepared", "recipe_sha256": sha256(RECIPE_PATH),
                  "sources": source_records, "input_sha256": files,
                  "prepared_at": datetime.now(UTC).isoformat(),
                  "elapsed_seconds": time.monotonic() - started,
                  "scope": "No activation labels loaded or models fitted during preparation"}
        _write_json(output / "prepared.json", status)
        return status
    except Exception as error:
        _write_json(output / "failed.json", {"status": "failed", "type": type(error).__name__,
                                             "reason": str(error), "sources": source_records})
        raise


def synthetic_saturation_demo():
    margins, target = np.array([45., 44., 43., 42.]), np.array([1, 1, 0, 0])
    probabilities = expit(margins)
    return {"data_kind": "synthetic", "probability_ap": float(average_precision_score(
                target, probabilities)), "margin_ap": float(average_precision_score(target, margins)),
            "probability_unique": len(np.unique(probabilities)), "margin_unique": len(margins),
            "fits": 0, "api_calls": 0,
            "meaning": "A deterministic illustration of sigmoid rounding, not a biological result"}


def doctor(inputs=None):
    recipe = load_recipe()
    packages = {}
    for name, expected in PINNED_PACKAGES.items():
        try:
            actual = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            actual = None
        packages[name] = {"expected": expected, "installed": actual, "matches": actual == expected}
    preparation_error = None
    prepared_complete = False
    if inputs:
        try:
            _verified_input_receipt(Path(inputs).resolve())
            prepared_complete = True
        except (ValueError, OSError, KeyError) as error:
            preparation_error = str(error)
    return {"recipe": {"variants": len(recipe["candidates"]), "folds": len(recipe["fold_groups"]),
                       "sha256": sha256(RECIPE_PATH)}, "python": platform.python_version(),
            "sequence_packages": packages,
            "sequence_runtime_ready": all(p["matches"] for p in packages.values()),
            "prepared_inputs_complete": prepared_complete,
            "prepared_inputs_error": preparation_error,
            "next": "uv run --frozen --extra sequence escape-lab-reproduce prepare --output PATH",
            "network_requests": 0, "outcomes_read": 0}


def _verified_input_receipt(inputs):
    if not (inputs / "prepared.json").is_file():
        raise ValueError("Run prepare first; a complete input receipt is required")
    receipt = json.loads((inputs / "prepared.json").read_text())
    if receipt.get("status") != "prepared" or receipt.get("recipe_sha256") != sha256(RECIPE_PATH):
        raise ValueError("Prepared inputs belong to a different frozen recipe")
    required = {"supplement.xlsx", "predictor-scores.csv", "biology/manifest.json",
                "biology/references.json", "biology/feature-cache.json"}
    hashes = receipt.get("input_sha256", {})
    if not required.issubset(hashes):
        raise ValueError("Prepared receipt is missing required inputs")
    for relative, expected in hashes.items():
        path = (inputs / relative).resolve()
        if not path.is_relative_to(inputs) or not path.is_file() or sha256(path) != expected:
            raise ValueError("Prepared input checksum mismatch")
    return receipt


def headline_comparison(metrics, recipe):
    checks = {}
    for arm, expected in recipe["expected_headline"].items():
        measured = metrics[arm]
        checks[arm] = {
            "ap_matches_at_1e_8": abs(measured["mean_fold_high_ap"]
                                     - expected["mean_fold_high_ap"]) <= 1e-8,
            "high_hits_match": measured["high_hits"] == expected["high_hits"],
            "budget_matches": measured["budget"] == expected["budget"],
        }
    return {"matches_published_aggregate": all(all(c.values()) for c in checks.values()),
            "checks": checks, "expected": recipe["expected_headline"],
            "comparison_scope": "Published aggregate development metrics; no historical predictions read"}


def reproduce_headline(inputs, output, *, device="cpu", feature_cache=None):
    inputs, output = Path(inputs).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError("Use a fresh output directory; previous attempts remain intact")
    recipe = load_recipe()
    receipt = _verified_input_receipt(inputs)
    output.mkdir(parents=True)
    started, candidates = time.monotonic(), recipe["candidates"]
    protocol = {"recipe_sha256": sha256(RECIPE_PATH), "prepared_sha256": sha256(inputs / "prepared.json"),
                "model": "esm_biology_logistic", "stability": "legacy", "sequence_branches": "both",
                "fit_budget": 3, "device": device, "feature_mode": "verified cache" if feature_cache
                else "fresh extraction", "created_at": datetime.now(UTC).isoformat()}
    _write_json(output / "protocol.json", protocol)
    try:
        references = inputs / "biology/references.json"
        if feature_cache:
            sequence = load_features(candidates, references, Path(feature_cache))
        else:
            if not doctor()["sequence_runtime_ready"]:
                raise ValueError("Install the locked sequence extra before fresh feature extraction")
            sequence = extract_features(candidates, references, output / "sequence-features",
                                        device=device, batch_size=4)
        checkpoint = sequence.manifest["checkpoint"]
        if {name: value["sha256"] for name, value in checkpoint["files"].items()} != (
                recipe["sequence"]["checkpoint_sha256"]):
            raise ValueError("Checkpoint bytes differ from pinned historical checkpoint")
        records, coverage = load_development_activation(inputs / "supplement.xlsx", candidates)
        if set(records) != {candidate["candidate_id"] for candidate in candidates}:
            raise ValueError("Activation source does not cover the frozen complete cohort")
        store = BiologicalFeatureStore(inputs / "biology")
        augmenter = ActivationFeatureAugmenter(store, candidates, references, inputs / "predictor-scores.csv")
        rows, diagnostics = {"probability": [], "margin": []}, []
        fitting_start = time.monotonic()
        for fold, train, test in fold_partitions(candidates, recipe["fold_groups"]):
            model = FGFR2SequenceModel("esm_biology_logistic", store, augmenter, sequence,
                                       stability="legacy", sequence_branches="both")
            model.fit(train, {c["candidate_id"]: records[c["candidate_id"]] for c in train})
            predictions = model.prediction_outputs(test)
            diagnostics.append({"fold": fold, "model": model.diagnostics,
                                "coefficients_sha256": hashlib.sha256(model.model.coef_.tobytes()).hexdigest()})
            for arm, scores in predictions.items():
                for candidate, score in zip(test, scores, strict=True):
                    record = records[candidate["candidate_id"]]
                    rows[arm].append({"candidate_id": candidate["candidate_id"], "fold": fold,
                                      "group_id": candidate["group_id"], "truth": record["value"],
                                      "high": activation_label(record), "strong": activation_label(record, 10),
                                      "score": float(score)})
        fitting_seconds = time.monotonic() - fitting_start
        metrics = {arm: summarize(values) for arm, values in rows.items()}
        for relative, expected in receipt["input_sha256"].items():
            if sha256(inputs / relative) != expected:
                raise ValueError("Prepared inputs changed while the experiment was running")
        comparison = headline_comparison(metrics, recipe)
        result = {"status": "completed", "metrics": metrics, "comparison": comparison,
                  "fit_count": 3, "fitting_seconds": fitting_seconds,
                  "elapsed_seconds": time.monotonic() - started, "coverage": coverage,
                  "feature_mode": protocol["feature_mode"],
                  "feature_sha256": sequence.manifest["features_sha256"],
                  "historical_feature_bytes_match": sequence.manifest["features_sha256"]
                  == recipe["sequence"]["historical_feature_sha256"],
                  "checkpoint_device": checkpoint["device"], "checkpoint_packages": checkpoint["packages"],
                  "head_packages": {name: importlib.metadata.version(name) for name in
                                    ("numpy", "scipy", "scikit-learn", "biopython")},
                  "diagnostics": diagnostics,
                  "scope": "Fixed development reproduction; original final and external assays were not scored"}
        _write_json(output / "private-predictions.json", rows)
        (output / "private-predictions.json").chmod(0o600)
        _write_json(output / "results.json", result)
        _write_json(output / "status.json", {"status": "completed", "results_sha256": sha256(output / "results.json")})
        return result
    except Exception as error:
        _write_json(output / "status.json", {"status": "failed", "type": type(error).__name__,
                                             "reason": str(error)})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    diagnostics = commands.add_parser("doctor", help="Check local requirements without network or outcomes")
    diagnostics.add_argument("--inputs", type=Path)
    commands.add_parser("saturation-demo", help="Offline synthetic demonstration, no fitting")
    prepare = commands.add_parser("prepare", help="Download and rebuild frozen public inputs, no fitting")
    prepare.add_argument("--output", type=Path, required=True)
    prepare.add_argument("--workbook", type=Path, help="Optional original workbook, verified by checksum")
    reproduce = commands.add_parser("headline", help="Fit three fixed heads on development data only")
    reproduce.add_argument("--inputs", type=Path, required=True)
    reproduce.add_argument("--output", type=Path, required=True)
    reproduce.add_argument("--device", choices=("cpu", "mps"), default="cpu")
    reproduce.add_argument("--feature-cache", type=Path, help="Optional verified ESM cache; reported as cached")
    args = parser.parse_args()
    if args.command == "doctor":
        result = doctor(args.inputs)
    elif args.command == "saturation-demo":
        result = synthetic_saturation_demo()
    elif args.command == "prepare":
        result = prepare_inputs(args.output, workbook_path=args.workbook)
    else:
        result = reproduce_headline(args.inputs, args.output, device=args.device,
                                    feature_cache=args.feature_cache)
        result = {key: result[key] for key in ("status", "comparison", "fit_count", "fitting_seconds",
                                               "elapsed_seconds", "feature_mode")}
    print(json.dumps(result, indent=2, allow_nan=False))
    if args.command == "headline" and not result["comparison"]["matches_published_aggregate"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
