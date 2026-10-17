"""Publish aggregate artifacts without changing frozen research originals.

Never add predictions, requests, outcomes, models, workbooks, or credentials.
"""

import hashlib
import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DESTINATION = ROOT / "docs/archive"
RESULT_KEYS = {
    "metrics", "activation_metrics", "resistance_metrics", "gene_prior_metrics",
    "component_gate", "resistance_gate", "gate", "gates", "uncertainty",
    "baseline_reproduction", "v1_reproduction", "unbounded_reproduction",
    "fit_count", "fitting_seconds", "elapsed_seconds", "feature_extraction_seconds",
    "protocol_sha256", "integrity", "post_result_sensitivity", "coverage", "n",
    "status", "comparisons", "prototype_margin_passed", "interpretation",
    "nested_metrics", "external", "full_selections", "agent_statuses", "seconds",
    "token_accounting", "arms", "proposal_opportunities", "retrospective_reuse",
    "fresh_blind_validation", "limits", "numerical_controls_reproduced",
    "joint_commitment_sha256", "total_wall_seconds", "scientific_arm_seconds",
    "preserved_artifacts", "total_requests", "budget_charge_tokens",
}


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def clean_text(value):
    return value.replace(str(ROOT) + "/", "").replace(str(ROOT), ".")


def clean(value):
    if isinstance(value, dict):
        return {clean_text(k): clean(v) for k, v in value.items()}
    if isinstance(value, list):
        return [clean(v) for v in value]
    if isinstance(value, str):
        return clean_text(value)
    return value


def sources():
    paths = set(ROOT.glob("reports/**/*.md"))
    paths.update(ROOT.glob("docs/*.md"))
    paths.difference_update(ROOT / "docs" / name for name in (
        "EXPERIMENTS.md", "RESULTS.md", "REPRODUCING.md", "HACKATHON.md", "PUBLICATION.md"
    ))
    for folder in (
        "activation-v2", "fgfr2-sequence-pilot", "fgfr2-state-pilot",
        "fgfr2-state-bounded", "fgfr2-state-tree-pilot", "fgfr2-sequence-ap-pilot"
    ):
        paths.update((ROOT / "research" / folder).glob("*.md"))
    paths.update(ROOT.glob("reports/mechanistic-pilot/*/results.json"))
    for name in (
        "reports/improvement/20261017T194348Z-5e66cbab/final-results.json",
        "reports/revision-v2/20261017T201321Z-a29b17a8/results.json",
        "reports/model-swap-v2/20261017T210043Z-a4800896/results.json",
        "reports/research-pilot/you-pilot-20261017/results.json",
    ):
        paths.add(ROOT / name)
    return sorted(p for p in paths if p.is_file())


def main():
    selected = sources()
    included = set(selected)
    manifest = []
    for path in selected:
        original = path.read_bytes()
        relative = path.relative_to(ROOT)
        if path.suffix == ".json":
            data = json.loads(original)
            public = {k: v for k, v in data.items() if k in RESULT_KEYS}
            text = json.dumps(clean(public), indent=2, allow_nan=False) + "\n"
            transformation = "allowlisted aggregate fields; local root removed"
        else:
            text = clean_text(original.decode())

            def link(match, source_path=path):
                label, target = match.group(1), match.group(2).strip("<>")
                if target.startswith(("https://", "http://", "#", "mailto:")):
                    return match.group(0)
                bare = target.split("#", 1)[0]
                candidate = (source_path.parent / bare).resolve()
                root_candidate = (ROOT / bare).resolve()
                if candidate in included:
                    return match.group(0)
                if root_candidate in included:
                    return f"[{label}]({os.path.relpath(root_candidate, source_path.parent)})"
                return f"{label} (`{target}`; local artifact)"

            text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link, text)
            results_link = os.path.relpath(ROOT / "docs/RESULTS.md",
                                           (DESTINATION / relative).parent)
            text = (f"> Historical research record. See [current results]({results_link}) "
                    "for the release interpretation.\n\n" + text)
            text = "\n".join(line.rstrip() for line in text.splitlines()) + "\n"
            transformation = "local root removed; private artifact links annotated; banner added"
        destination = DESTINATION / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(text)
        manifest.append({"source": relative.as_posix(), "source_sha256": digest(original),
                         "public_path": destination.relative_to(ROOT).as_posix(),
                         "public_sha256": digest(text.encode()), "transformation": transformation})
    DESTINATION.mkdir(parents=True, exist_ok=True)
    (DESTINATION / "manifest.json").write_text(json.dumps({
        "scope": "Curated aggregate-only publication; no row outcomes, traces, or credentials",
        "artifacts": manifest,
    }, indent=2) + "\n")
    print(f"Exported {len(manifest)} reviewed summary artifacts.")


if __name__ == "__main__":
    main()
