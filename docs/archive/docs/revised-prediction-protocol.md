> Historical research record. See [current results](../../RESULTS.md) for the release interpretation.

# Revised FGFR prediction protocol

This is the current protocol for `uv run escape-lab improve`. It supersedes the
active-search architecture and evaluation portions of the original consolidated
project plan. The earlier Word and Markdown plan files remain historical records.
The completed run's [results and audit record](../reports/improvement/results-summary.md)
show a functioning autonomous loop but worse final prediction performance. The
final evaluation is consumed. This document records the protocol used, not a
claim of successful improvement.

## Question and scope

Can a Nemotron agent improve two numerical resistance predictors by proposing
bounded configurations from aggregate development errors, compared with fixed
models and a matched two-trial random configuration search?

The two prediction targets are pemigatinib and futibatinib resistance measured
separately in engineered cancer-cell assays. Neither target measures patient
response, mutation emergence, or the efficacy of combining or sequencing drugs.
The earlier experiment vault and active-replay commands remain separate tools.

## Data and features

Eight public MAVEdb score sets provide 8,407 paired nucleotide-variant records,
normalized to 7,522 protein variants. Distinct nucleotide representations remain
in the evidence archive. An operational positive requires three finite,
nonnegative replicate measurements, median enrichment above 1.5, and enrichment
above 1 in at least two replicates. Missing, invalid, or discordant representations
are inconclusive. Inconclusive labels are excluded separately for each drug and
their counts are reported; they are not relabeled negative.

Features use only wild-type and mutant amino acids, optional gene identity,
sequence context, biochemical substitution attributes, and aligned position.
Assay scores, outcome-derived categories, candidate IDs, and split assignments
are not predictor features. Feature vocabularies are fitted on training data.

Twelve-residue aligned blocks are kept together across FGFR1–4. The fixed split
contains 17 development blocks, five validation blocks, and seven final blocks.
Published illustration blocks are restricted to development. All configurations
use the same three grouped development folds and the same training information.

## Bounded search and development selection

The runner first evaluates three fixed configurations with full features and
gene identity enabled:

| Family | Fixed settings |
|---|---|
| Logistic regression | Inverse regularization `C=1`, at most 1,000 iterations |
| Extra Trees | 100 trees, minimum leaf size 3, `max_features=0.8` |
| Histogram gradient boosting | At most 100 iterations, maximum depth 3, learning rate 0.1, L2=1 |

The fixed reference is the configuration with the highest development macro AP.
For each drug, AP is computed across all pooled out-of-fold predictions; macro AP
is the arithmetic mean of those two scores. It is not a fold-average score.
Brier scores, per-gene and per-fold summaries, inconclusive counts, timing, and
any numerical prior fallback are retained as diagnostics.

The ordinary control samples two distinct effective configurations from the
same allowed space. This is a small random-search comparator, not exhaustive
AutoML. Its evaluations run before agent proposals so inference latency cannot
reduce the control's opportunities.

Nemotron receives two proposal opportunities. Each must return a complete,
schema-validated configuration and a brief rationale. A single format repair is
allowed per proposal. The agent can choose model family, feature set, gene
inclusion, and enumerated numerical settings; it cannot write or execute code,
change labels or splits, access measurements, or alter scoring. Effective-config
hashes ignore inactive settings, so changing an irrelevant parameter is not a
new experiment. Duplicate or failed proposals consume an opportunity and remain
explicit failures. They are never silently replaced by baseline configurations.

Feedback contains only development aggregate metrics, permitted configuration
records, and gene/fold error summaries. Validation and final results never enter
agent feedback. Random-search results are not agent-visible; an independently
matching proposal reuses the same deterministic evaluation without receiving an
extra trial.

For both search methods, a development revision is retained only when:

1. Its macro AP strictly exceeds its current incumbent's macro AP.
2. Neither drug's AP is more than 0.02 below the fixed reference's AP.

Both searches start from the same fixed reference. All attempts, including
rejections, timeouts, and invalid proposals, remain in the report. Two matched
configuration trials do not imply identical hardware time; actual timings are
reported.

## One validation pass and the final gate

Validation begins only after every required control trial and both agent
proposal evaluations complete. Failed, skipped, or duplicate agent proposals
make the run `incomplete`; validation and final evaluation are both skipped.
The CLI exits with status 1 for an incomplete run or a time-budget overrun.

After this completeness check, the selected fixed, random, and agent
configurations are each evaluated on validation once. Identical configurations
share a cached evaluation. These models are fitted on development data; there is
no subsequent agent iteration from validation results.

The same promotion rule applies to the random and agent selections: validation
macro AP must strictly beat the fixed reference, with no per-drug drop greater
than 0.02. A selection that fails promotion reverts to the fixed reference for
the final comparison. The report preserves the unsuccessful candidate's metrics
and explicitly records that it was not promoted.

The final test is eligible only when all control trials and both agent proposal
evaluations completed, the agent was promoted on validation, final evaluation
was enabled, and more than 30 seconds remain. Otherwise it stays sealed. A failed
proposal or skipped opportunity cannot be presented as a complete two-trial run.

Before final labels are requested, the runner:

1. Freezes source, protocol, data, and selected model checksums.
2. Exclusively reserves `data/processed/final-evaluation.json`, a receipt shared
   across report directories.
3. Uses the already fitted, checksum-verified models to predict every final case.
4. Saves all methods' final predictions and records their checksum in the receipt.
5. Re-verifies the freeze, then requests final labels for scoring.

The archive is a monolithic JSON file. The trusted evaluator logically restricts
the labels it exposes; development loading may parse the archive's bytes. This
is a capability boundary around the agent, not an OS sandbox or a claim that
the evaluator never parsed final bytes.

The preset prototype success criterion is an absolute macro-AP gain of at least
0.05 over **both** the fixed and random-search controls, with neither per-drug AP
more than 0.02 below either comparator. A positive result is a retrospective
benchmark result, not evidence of a cure or general superiority in biology.

Once the final receipt exists, fresh tuning is blocked. An interrupted or failed
reserved final run remains spent; do not delete the receipt or choose another
output directory to obtain a replacement test. A changed scientific protocol
requires new untouched evidence, not repeated optimization against these labels.

## Running and inspecting it

From the project root, with the dataset prepared and the existing local `.env`
configured with `CRUSOE_API_KEY`:

```sh
uv run escape-lab status
uv run escape-lab improve --model nvidia/NVIDIA-Nemotron-3.5-Super-VL --request-timeout 60 --max-completion-tokens 8192 --max-seconds 600
```

`--model` must be an exact accessible model ID and overrides `CRUSOE_MODEL`.
`--request-timeout` defaults to 60 seconds, accepts integers from 1 through 90,
and is persisted as `generation.request_timeout_seconds` in `protocol.json`.
`--max-completion-tokens` defaults to 8,192 and accepts 2,048, 4,096, or 8,192.
It includes reasoning tokens used before a structured tool call and is recorded
as `generation.max_tokens`. An explicit `length` finish reason immediately
returns a failed `truncated_response`; it does not receive an identical-cap
format repair. Attempt diagnostics retain allowlisted finish-reason metadata,
without raw assistant content or private reasoning.
Use `--development-only` to keep final evaluation disabled; validation still
requires all development proposal evaluations to complete.

The runner's allowed wall-clock budget is 60–600 seconds, starting after client
setup and model discovery. Numerical jobs are run in terminable worker processes.
The CLI uses the selected per-call completion allowance, no transient HTTP retry,
the selected request timeout, and a 50,000-token stop threshold. The two proposal
opportunities and maximum one format repair per proposal are unchanged. Token and HTTP limits
are not a hard billing ceiling. A POSIX process alarm interrupts the overall run
at its wall-clock deadline; worker cleanup and writing the failure report can
take slightly longer. Elapsed time and overruns are recorded.

Three operational failures and the declared fourth attempt are preserved
in [the inference retry note](../reports/improvement/inference-retry-note.md).
The retries change inference configuration, error handling, and tool schema
representation, not the admissible configuration space or other scientific
choices. Validation reuse is not independent replication. No improvement
should be inferred from the earlier run's inherited-baseline validation row.

Artifacts are stored under `reports/improvement/<RUN_ID>/`:

- `protocol.json`: predeclared criteria, configurations, runtime, and fingerprints.
- `proposal-1.json` and `proposal-2.json`: attempted agent proposals when started.
- `development/` and `validation/`: metrics, predictions, diagnostics, and models.
- `freeze.json`: selected configurations and source/data/model checksums.
- `final-predictions-committed.json`, `final-results.json`, and
  `final-predictions.csv`: created only if final evaluation reaches those stages.
- `report.md` and `report.json`: run outcome, trial history, selection, timing,
  usage, and whether the final test stayed sealed.

Only load the runner's own trusted local `.pkl` artifacts. Pickle can execute
code; a checksum does not make an arbitrary uploaded model safe.

## Limits and data rights

The final set has only seven aligned sequence blocks across four related genes.
The paired block bootstrap is descriptive sensitivity analysis. Neither blocks,
random seeds, nor thousands of variants constitute independent patient or
laboratory replication. Broad clinical generalization and freedom from public
data memorization remain unverified. No comparison with unavailable historical
CAGI submissions is claimed.

MAVEdb source measurements and derived measurement data use CC BY-SA 4.0.
Preserve attribution, accessions, and applicable share-alike terms when sharing
data or adaptations. New project code does not copy unlicensed study analysis
code. Exact sources and checksums are retained in `data/raw/source_manifest.json`.

- [Original FGFR study](https://doi.org/10.1038/s41588-025-02431-8)
- [Example MAVEdb record](https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00001256-a-9)
- [Experimental challenge criteria](https://genomeinterpretation.org/cagi7-fgfr.html)
- [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
