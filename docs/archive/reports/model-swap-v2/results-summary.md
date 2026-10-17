> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# GLM-5.3 substitution: no demonstrated overall prediction improvement

One matched run of Crusoe `zai-org/GLM-5.3` completed against the most recent
Nemotron v2 experiment on October 17, 2026. Internal predictions were identical.
External mean correlation increased descriptively, but the paired uncertainty
interval includes substantial worsening and improvement, and GLM remained below
the fixed structural baseline. GLM also completed fewer useful proposals under
the unchanged budget. This does not establish a stronger-model advantage.

| Method | Internal pooled macro AP ↑ | Internal macro Brier ↓ | External mean Spearman ↑ | External pem | External fut |
|---|---:|---:|---:|---:|---:|
| Original fixed numerical model | 0.034413 | 0.035374 | −0.018063 | 0.254280 | −0.290406 |
| Fixed structural numerical model | 0.034628 | 0.035781 | **0.405681** | 0.298904 | **0.512458** |
| Matched random search | **0.034938** | **0.034593** | 0.283740 | 0.473917 | 0.093563 |
| Exhaustive numerical grid | **0.034938** | **0.034593** | 0.182598 | 0.412320 | −0.047123 |
| Nemotron 3.5 Super VL | 0.030516 | 0.035174 | 0.208915 | 0.103481 | 0.314349 |
| GLM-5.3 | 0.030516 | 0.035174 | 0.296952 | **0.508559** | 0.085345 |

The external change was uneven: pemigatinib correlation increased while
futibatinib correlation decreased. Both models selected chemistry C0.1, structure
C0.1, and chemistry C0.1 for outer folds 0–2. For the external prediction GLM
selected `interactions-logistic-c0.1`; Nemotron selected `structure-logistic-c1`.

The original paired whole-block bootstrap (2,000 draws, seed 6100) gives:

- GLM minus Nemotron internal AP: **0.000000**, interval **[0, 0]**, because the
  committed predictions are identical.
- GLM minus Nemotron external correlation: **+0.088037**, 95% descriptive
  percentile interval **[−0.358545, +0.368979]**.
- GLM minus fixed structure external correlation: **−0.108729**, interval
  **[−0.519605, +0.082645]**.

All 2,000 draws were defined for each comparison. Internal evaluation covers 3,730
variants in 15 aligned blocks, with 3,707 observed pem and 3,695 observed fut
labels. External evaluation reuses 25 measurements per drug in nine blocks.
These intervals condition on committed predictions and exclude model-fitting,
LLM-selection and independent biological replication uncertainty. The external
panel was already spent: this is an explicitly authorized retrospective comparison,
not fresh blind validation. It is a preclinical cross-assay result, not clinical
accuracy or a treatment claim.

| Operational measure | Nemotron | GLM-5.3 |
|---|---:|---:|
| Allocated proposal opportunities | 12 | 12 |
| Evaluated proposals | 12 | 6 |
| Requests, including format repairs | 14 | 17 |
| Format repair requests | 2 | 5 |
| Opportunities lost to response truncation | 0 | 4 |
| Opportunities still invalid after repair | 0 | 2 |
| Input tokens | 95,604 | 84,687 |
| Output tokens | 53,354 | 80,813 |
| Total reported tokens | 148,958 | 165,500 |
| Scientific wall time | 243.55 s | 393.63 s |
| Uncached list-price estimate | $0.030902 | $0.474139 |
| Estimate using reported cached tokens | unavailable | $0.429925 |

GLM's four searches were all degraded, with 1, 1, 2 and 2 evaluated proposals.
Each still met the original requirement of at least one evaluated proposal before
external scoring. There were no HTTP failures, operational retries, scientific
restarts or alternative scientific attempts. Four responses reached the original
8,192-token cap; no extra opportunities or repairs were added. Neither the shared
180,000-token stop threshold nor the 1,200-second run deadline was reached.

Costs use [Crusoe's official prices](https://www.crusoe.ai/cloud/pricing), also
confirmed by account model discovery: GLM $1.40 input / $4.40 output / $0.26 cached
input, and Nemotron $0.10 / $0.40 / $0.08 per million tokens. GLM reported 38,784
cached input tokens. Its cache-adjusted estimate is
`((84687 - 38784) * 1.40 + 38784 * 0.26 + 80813 * 4.40) / 1000000`.
Original Nemotron records lack cached-token subcounts. These are estimates, not
verified invoices. GLM reported 73,141 reasoning tokens as a subset of its 80,813
completion tokens; they are not charged twice. Cached-token semantics follow
[Crusoe's metrics definitions](https://docs.crusoecloud.com/serverless-inference/inference-metrics/index.html).

The audit verified all 28 original protocol hashes, 334 preserved original files,
unchanged dependency versions and original receipts. All four numerical controls'
predictions, metrics, folds, selected configurations and random candidate lists
reproduced within 1e-12. The original seven bootstrap contrasts also reproduced.
The full suite passed **284 tests** and Ruff passed. The sole test warning is an
existing Starlette/httpx deprecation.

Implementation added two separate modules (`model_swap.py`, `model_swap_analysis.py`),
two test files and a comparison protocol document. No original scientific source,
dataset, feature, model artifact or receipt was modified. The wrapper uses the
original planner, numerical fitting, selection and external metric functions.
Only inference model and creation time differ in the original protocol fields;
additional fields record provenance and retrospective authorization.

Agent prompts and tool schemas are unchanged; the planner saw only its own
allowed inner-fold aggregates. Original Nemotron winning configurations, test
predictions and external outcomes never entered its requests. External scoring
read authenticated released S13 records only after new predictions, model files
and selections were committed. It did not open workbook cells, evaluate the old
v1 final set, or access reserved RLY-4008 numerical outcomes.

Wire settings remained `tool_choice=required`, `max_tokens=8192`, `stream=false`,
60-second timeout and zero retries. Temperature, inference seed, top-p and
reasoning controls were unspecified in both runs. Thus model-specific serving
defaults are not controlled, and requests are not byte-identical because history
and remaining time change. GLM's discovery metadata listed a completion limit of
zero; no limit was inferred from that value, and all actual 8,192-token requests
were accepted. Every response identified `zai-org/GLM-5.3`. Timing includes shared
numerical work and audit overhead, with a brief test collection and read-only
status checks overlapping the run; it is not a controlled serving-speed benchmark.

Recorded execution commands, from the Escape Lab project directory:

```sh
uv run python -m escape_lab.model_swap --root .
uv run python -m escape_lab.model_swap_analysis reports/revision-v2/20261017T201321Z-a29b17a8 reports/model-swap-v2/20261017T210043Z-a4800896
uv run pytest -q
uv run ruff check src tests
```

The analysis refuses to overwrite existing outputs. The scientific command
creates a new paid attempt and was executed only once for this comparison.

Full evidence: [all original endpoints and fold/configuration tables](20261017T210043Z-a4800896/results-summary.md),
machine-readable comparison and uncertainty (`20261017T210043Z-a4800896/comparison.json`; local artifact),
verification and cost audit (`20261017T210043Z-a4800896/verification-audit.json`; local artifact),
saved protocol (`20261017T210043Z-a4800896/protocol.json`; local artifact), and
actual request audit (`20261017T210043Z-a4800896/requests.jsonl`; local artifact).
