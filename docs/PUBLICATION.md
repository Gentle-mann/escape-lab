# Publication scope

This release contains project source, synthetic-fixture tests, the complete experiment ledger, and sanitized aggregate reports. Source data, per-variant outcomes/predictions, evaluation databases, trained-model pickles, full inference transcripts, raw papers, API credentials and pretrained weights stay outside Git.

The original local artifacts remain unchanged. `scripts/export_public_results.py` creates copies under `docs/archive/`, removes the local project root from paths, and turns references to unpublished local artifacts into annotated text. Aggregate JSON uses an explicit top-level field allowlist. `docs/archive/manifest.json` records SHA-256 hashes for both source and published copies. Hashes establish correspondence, not independent scientific correctness.

The current README and `docs/RESULTS.md` distinguish implemented capabilities, measured outcomes and the research agenda. Historical plans preserve the design decisions that led here. Each result retains its endpoint, evaluation setting and evidence record.

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

## Development after submission

This increment adds canonical ranking/probability APIs, one fixed stability experiment, sequence ablations, diagnostic agent tools, retained research policies, a public-source reproduction command, a locked sequence environment and publication checks in CI. It also updates the project framing around the ambition to anticipate cancer's next move. The submitted six-hour result remains identifiable.

- Full local software suite: **659 tests passed**; Ruff passed. The local demo suite also passed **86 tests**, and the revised ambition page and presentation opening were checked in the in-app browser. The existing Starlette/httpx deprecation warning remains.
- The exact staged public snapshot in an empty temporary directory passed **643 tests**, with 16 optional local-reference tests skipped. Ruff, archive/metric/link checks, the no-network environment doctor and the synthetic saturation command all passed there without private inputs or credentials.
- Fresh public downloads, fresh ESM features and three fixed head fits reproduced the published AP and hit counts exactly in a source-only temporary workspace. Preparation plus extraction/evaluation took **100.55 seconds**, excluding dependency installation. The [verification receipt](archive/reports/reproduction/verification.json) records packages, inputs, source hashes and the retained preparation failure.
- The fixed 12-fit development comparison rejected the stability challenger and reported both representation ablations. Historical scores reproduced exactly; 529 protected artifact hashes stayed unchanged.
- Live Nemotron completed a five-call synthetic diagnostic/policy run and a five-call replay using the retained policy. The synthetic referee evaluated 48 newly generated cases. These are policy mechanism and integration checks, separately reported from biological prediction results. Operational failures and the corrective schema probe are preserved in the [agent records](archive/reports/agent-improvements/README.md).
- The archive now contains **98 hash-registered artifacts**. All 57 entries from the preceding publication snapshot remain unchanged. CI checks the manifest, public relative links, and headline AP/retrieval figures against the archived measurements.
- No original final or external biological assay was scored. Local presentation files remain separate from the GitHub source publication; this increment does not replace the submitted Loom video or deploy a new version of the hackathon demo.
