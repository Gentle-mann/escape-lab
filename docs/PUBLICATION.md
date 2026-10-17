# Publication scope

This release contains project source, synthetic-fixture tests, the complete experiment ledger, and sanitized aggregate reports. Source data, per-variant outcomes/predictions, evaluation databases, trained-model pickles, full inference transcripts, raw papers, API credentials and pretrained weights stay outside Git.

The original local artifacts remain unchanged. `scripts/export_public_results.py` creates copies under `docs/archive/`, removes the local project root from paths, and turns references to unpublished local artifacts into annotated text. Aggregate JSON uses an explicit top-level field allowlist. `docs/archive/manifest.json` records SHA-256 hashes for both source and published copies. Hashes establish correspondence, not independent scientific correctness.

Some historical plans mention future capabilities and superseded designs. The current README and `docs/RESULTS.md` govern the release's claims. In particular, do not present an implemented agent loop as demonstrated agent superiority, a margin repair as new biological knowledge, or a repeated evaluation as a fresh holdout.

The larger repeated research/placebo study was in progress when the release was assembled. Its implemented source and protocol are included; it has no completed result in this snapshot unless the ledger explicitly identifies a completed result artifact. Later results should be published in a separate reviewed update.

## Licenses and attribution

Project source is MIT licensed under `LICENSE`. Documentation and derived aggregate measurement summaries are offered under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). Preserve attribution to the FGFR study authors, original datasets and this project's analysis when redistributing these summaries. Third-party weights, structures and dependencies retain their source licenses; no rights to those artifacts are implied by the code license.

The primary FGFR paper is Tangermann and colleagues, https://doi.org/10.1038/s41588-025-02431-8. Original source licenses and checksums are retained by the ingestion pipeline. Structure and sequence sources are cited in the corresponding archived protocols.

## Release verification

- Original local workspace: 583 tests passed; Ruff passed.
- Exact staged public snapshot in an empty temporary checkout: 567 tests passed, 16 optional local-reference tests skipped; Ruff passed. One existing Starlette/httpx deprecation warning remains.
- Reviewed 156 staged text files, approximately 1.62 MB, before the verification paragraph: no detected credential patterns, forbidden raw-data paths, oversized files, personal filesystem roots in docs, or row-outcome keys in archived aggregate JSON.
- All 54 archive entries match their recorded public SHA-256 hashes; public Markdown links resolve.
- Scoped pre-publication review covered vault reservations, capability checks, bounded proposal validation and research-evidence isolation. This is not a comprehensive security audit.
- No scientific model rerun or external holdout evaluation was performed for publication.
