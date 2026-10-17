# Publication scope

This release contains project source, synthetic-fixture tests, the complete experiment ledger, and sanitized aggregate reports. Source data, per-variant outcomes/predictions, evaluation databases, trained-model pickles, full inference transcripts, raw papers, API credentials and pretrained weights stay outside Git.

The original local artifacts remain unchanged. `scripts/export_public_results.py` creates copies under `docs/archive/`, removes the local project root from paths, and turns references to unpublished local artifacts into annotated text. Aggregate JSON uses an explicit top-level field allowlist. `docs/archive/manifest.json` records SHA-256 hashes for both source and published copies. Hashes establish correspondence, not independent scientific correctness.

Some historical plans mention future capabilities and superseded designs. The current README and `docs/RESULTS.md` govern the release's claims. In particular, do not present an implemented agent loop as demonstrated agent superiority, a margin repair as new biological knowledge, or a repeated evaluation as a fresh holdout.

The repeated research/placebo and selector experiments completed on 17 October 2026. This update publishes their [aggregate results](archive/reports/repeated-research/20261017-experiments-1-2/results.json), [completed review](archive/reports/repeated-research/20261017-experiments-1-2-review.md), and [verification and choice summary](archive/reports/repeated-research/20261017-experiments-1-2-audit.json). Both experiments are inconclusive under the predeclared improvement gates. The 90 searches reused development data and scored no final or external panel. Raw outcomes, per-variant predictions, fitted models and request logs remain local.

## Licenses and attribution

Project source is MIT licensed under `LICENSE`. Documentation and derived aggregate measurement summaries are offered under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). Preserve attribution to the FGFR study authors, original datasets and this project's analysis when redistributing these summaries. Third-party weights, structures and dependencies retain their source licenses; no rights to those artifacts are implied by the code license.

The primary FGFR paper is Tangermann and colleagues, https://doi.org/10.1038/s41588-025-02431-8. Original source licenses and checksums are retained by the ingestion pipeline. Structure and sequence sources are cited in the corresponding archived protocols.

## Original publication snapshot verification

These checks describe the original publication snapshot, before the completed repeated-study results were added. They are historical checks, not verification counts for this update.

- Original local workspace: 583 tests passed; Ruff passed.
- Exact staged public snapshot in an empty temporary checkout: 567 tests passed, 16 optional local-reference tests skipped; Ruff passed. One existing Starlette/httpx deprecation warning remains.
- Reviewed 156 staged text files, approximately 1.62 MB, before the verification paragraph: no detected credential patterns, forbidden raw-data paths, oversized files, personal filesystem roots in docs, or row-outcome keys in archived aggregate JSON.
- All 54 archive entries match their recorded public SHA-256 hashes; public Markdown links resolve.
- Scoped pre-publication review covered vault reservations, capability checks, bounded proposal validation and research-evidence isolation. This is not a comprehensive security audit.
- No scientific model rerun or external holdout evaluation was performed for publication.

## Completed repeated-study update verification

- Full local suite: 591 tests passed; Ruff passed. Eight synthetic exporter tests cover aggregate publication and rejection of nested outcome records, request bodies, credential fields and scalar-to-array injection. The existing Starlette/httpx deprecation warning remains.
- All 57 archive entries match their source and publication hashes. The 54 previously published entries are unchanged; three new files contain the completed review, aggregate metrics and audit summary.
- Checked 110 relative links across 64 public documentation files. No personal filesystem roots, configured credentials or API-key patterns were found. The two new JSON summaries pass explicit recursive field and type schemas; an independent review found no raw outcomes or agent payloads in the exports.
- Reverified the original run's 314 artifact hashes, 388 commitment entries and 292 preserved prior-artifact hashes. Raw outcomes, row-level predictions, fitted models and request logs remain outside Git.
- This publication update reran software tests only. It did not rerun a scientific experiment or score another holdout.
