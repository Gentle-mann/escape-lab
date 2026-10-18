"""Verify public archive hashes, documentation links, and the historical AP card.

Reads only explicitly published documentation and aggregate reports. No scientific
evaluation, private outcomes, environment files, or network requests are involved.
"""

import argparse
import hashlib
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlsplit


class PublicationError(ValueError):
    pass


PUBLIC_DOCS = ("RESULTS.md", "EXPERIMENTS.md", "REPRODUCING.md", "HACKATHON.md", "PUBLICATION.md")
AP_RESULT = "docs/archive/reports/mechanistic-pilot/fgfr2-sequence-ap-v1/results.json"


def _inside(path, parent):
    return path.resolve().is_relative_to(parent.resolve())


def verify_headline(root):
    root = Path(root)
    metrics = json.loads((root / AP_RESULT).read_text())["metrics"]
    before, after = metrics["probability"], metrics["margin"]
    gain = 100 * (after["mean_fold_high_ap"] / before["mean_fold_high_ap"] - 1)
    for relative in ("README.md", "docs/RESULTS.md"):
        text = (root / relative).read_text()
        rows = [line.replace("**", "").split("|")[1:-1]
                for line in text.splitlines() if line.startswith("|")]
        rows = [[cell.strip() for cell in row] for row in rows]
        for label, data in (("probabilit", before), ("raw decision margin", after)):
            matches = [row for row in rows if row and label in row[0].lower()
                       and (label != "probabilit" or "original" in row[0].lower())]
            if len(matches) != 1:
                raise PublicationError(f"Missing or ambiguous headline row in {relative}: {label}")
            row = matches[0]
            expected_ap = f"{data['mean_fold_high_ap']:.5f}"
            if len(row) < 4 or row[1] != expected_ap:
                raise PublicationError(f"AP headline mismatch in {relative}: {label}")
            for column, endpoint in ((2, "high"), (3, "strong")):
                positives = sum(fold[endpoint]["positives"]
                                for fold in data["by_fold"].values())
                expected = f"{data[endpoint + '_hits']}/{positives}"
                if row[column].replace(" ", "") != expected:
                    raise PublicationError(f"Retrieval headline mismatch in {relative}: {label}")
            if not re.search(rf"\b{data['budget']}\b", text):
                raise PublicationError(f"Missing headline budget in {relative}")
        if f"{gain:.1f}%" not in text:
            raise PublicationError(f"Relative-gain headline mismatch in {relative}")
    return {"before_ap": before["mean_fold_high_ap"], "after_ap": after["mean_fold_high_ap"],
            "relative_gain_percent": gain}


def check_publication(root, check_headline=True):
    root = Path(root).resolve()
    archive = root / "docs/archive"
    manifest = archive / "manifest.json"
    entries = json.loads(manifest.read_text())["artifacts"]
    registered = set()
    markdown = {root / "README.md"}
    markdown.update(root / "docs" / name for name in PUBLIC_DOCS
                    if (root / "docs" / name).exists())
    markdown.update((root / "docs/reproduction").glob("*.md"))
    for entry in entries:
        relative = entry["public_path"]
        path = root / relative
        if not _inside(path, archive) or path.resolve() == manifest.resolve():
            raise PublicationError(f"Manifest path outside archive: {relative}")
        if path.resolve() in registered:
            raise PublicationError(f"Duplicate archive entry: {relative}")
        if not path.is_file():
            raise PublicationError(f"Missing archive artifact: {relative}")
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["public_sha256"]:
            raise PublicationError(f"Archive hash mismatch: {relative}")
        registered.add(path.resolve())
        if path.suffix == ".md":
            markdown.add(path)
    actual = {p.resolve() for p in archive.rglob("*") if p.is_file() and p != manifest}
    if actual - registered:
        raise PublicationError("Archive contains unregistered artifacts: " +
                               ", ".join(sorted(str(p.relative_to(root))
                                                for p in actual - registered)))
    links = 0
    for path in sorted(markdown):
        text = re.sub(r"```.*?```", "", path.read_text(), flags=re.DOTALL)
        text = re.sub(r"`[^`]*`", "", text)
        for match in re.finditer(r"\[[^\]]*\]\(([^)]+)\)", text):
            target = match.group(1).strip().strip("<>")
            parsed = urlsplit(target)
            if parsed.scheme or parsed.netloc or not parsed.path:
                continue
            destination = path.parent / unquote(parsed.path)
            if not _inside(destination, root) or not destination.exists():
                raise PublicationError(f"Broken public link in {path.relative_to(root)}: {target}")
            links += 1
    result = {"status": "pass", "artifacts": len(registered),
              "markdown_files": len(markdown), "relative_links": links}
    if check_headline:
        result["headline"] = verify_headline(root)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    print(json.dumps(check_publication(args.root), indent=2))


if __name__ == "__main__":
    main()
