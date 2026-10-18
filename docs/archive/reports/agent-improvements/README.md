# Research policies that improve and persist

Escape Labs is building a research system that learns how to investigate failures,
retains useful procedures, and checks them on new tasks before reuse. This increment
connects numerical diagnostics to the planner and makes research policies versioned,
verifiable artifacts.

## Implemented workflow

1. An independent numerical executor measures score saturation, training-scale
   issues, permitted development feature shift, and probability versus raw-score AP.
2. Nemotron on Crusoe chooses one bounded tool at a time. It can inspect those
   diagnostics, request the fixed WT/mutation-delta ablation, propose an earned
   policy rule, or stop. Failed and rejected calls remain in the run record.
3. Policies contain only a parent hash, version, typed diagnostic rules, and proposal
   provenance. Their filenames are their content hashes. They cannot contain
   examples, labels, protein identities, prompts, or executable code.
4. The policy and referee are frozen before the synthetic heldout cases are
   generated. The referee evaluates the parent and challenger and writes an immutable
   promotion receipt. Reusing a retained policy through the CLI requires that receipt.

The original resistance planner also has an opt-in `diagnostics-v3` protocol. It
receives allowlisted rejection codes and aggregate model diagnostics in subsequent
calls and can explicitly stop. Its acceptance guards are unchanged. The historical
`historical-v2` protocol remains the default so archived comparisons keep their
original prompt and tool contract. The old external-panel runner is unchanged.

## Measured mechanism test

The [live Nemotron verification](nemotron-live/report.json) completed five calls on
Crusoe: inspect saturation, compare rankings, inspect feature shift, propose two
typed diagnostic rules, and stop. The frozen referee accepted the policy after
**48/48** generated cases were diagnosed correctly, compared with **12/48** for
the initial empty policy, with zero false positives. The run used 7,641 reported
tokens. These are predefined numerical failure families, evaluated by a deterministic
rule interpreter; the result verifies policy retention, not new biological accuracy.

The [operational attempt record](live-attempts.json) preserves two rejected live
attempts and the synthetic schema probe that identified the argument wrapper. Replacing the empty
parameter schema with an explicit bound-partition target fixed the observed tool
argument wrapper. No failed attempt was silently retried or removed.

The [live retained-policy replay](nemotron-retained-replay/report.json) then loaded
that policy through its promotion receipt. Nemotron used five calls to inspect the
published FGFR saturation and scaling summaries, compare the recorded rankings,
retrieve the fixed sequence ablation and stop. It used 10,024 reported tokens and
ran zero fits. Its input policy status is `synthetic_referee_promoted`; it made no
new policy proposal or scientific prediction. Fractions in the replay diagnostics
are means over the three fitted folds, while AP is the specified mean-fold metric.

The current [offline synthetic report](synthetic-transfer-v3/report.json) uses the
same executable function-tool protocol with a clearly labeled scripted client.
It proposes both diagnostic rules from generated training failures, freezes them,
and evaluates 48 different generated cases:

| Policy | Correct complete diagnoses | Success rate | Cases with false positives |
| --- | ---: | ---: | ---: |
| Initial policy | 12 / 48 | 25% | 0 |
| Revised policy | 48 / 48 | 100% | 0 |

The cases cover four predefined families: clean scores, saturated probabilities,
feature shift, and both failures together. Heldout seeds, sample sizes, and magnitudes
differ from training cases. This proves that the implemented policy can persist and
transfer within these generated failure families. It is a software mechanism test;
it does not measure biological accuracy or demonstrate an LLM discovering a new
scientific strategy. The baseline policy initially has no diagnostic rules.

The [current development replay](development-replay-v3/report.json) runs those
bounded tools on published aggregate results, without model fitting or row-level
outcome access. It reports mean-fold activation AP **0.09732 to 0.11135** for the
original ranking repair, and the already-measured sequence ablations. It creates
an unpromoted rule proposal. It does not turn those reused folds into unseen data.

The original `synthetic-transfer/` and `development-replay/` records preserve the
first offline runs. The `*-v2/` records add allowlisted response-shape audit fields and
fix missing convergence/fallback diagnostics being incorrectly summarized as
zero in the first replay; absent values are now omitted. The AP results and policy
transfer counts are unchanged. The current `*-v3/` records use activation planner protocol v2: every numerical tool requires `target="bound_partition"`. This avoids the empty-object schema wrapper seen in a live synthetic provider probe while rejecting arbitrary targets and extra fields. All offline calls have zero inference charges.

## Commands

These commands work without credentials or private datasets. Use a fresh output
directory for every attempt; an existing nonempty directory is never overwritten.

```sh
uv run --frozen python -m escape_lab.activation_diagnostic_cli \
  --mode synthetic --output reports/policy-smoke-01

uv run --frozen python -m escape_lab.activation_diagnostic_cli \
  --mode replay --root . --output reports/development-replay-01
```

For a billed Nemotron run on generated synthetic cases, explicitly provide
`CRUSOE_API_KEY` in the process environment and use:

```sh
uv run --frozen python -m escape_lab.activation_diagnostic_cli \
  --mode synthetic --provider crusoe \
  --model-id nvidia/NVIDIA-Nemotron-3.5-Super-VL \
  --max-actions 6 --max-completion-tokens 8192 \
  --max-total-tokens 60000 --timeout-seconds 60 \
  --output reports/nemotron-synthetic-01
```

Each numerical tool requires `target="bound_partition"`; the attached evaluator determines that partition. The CLI does not load `.env`. Model discovery and, when capability metadata is
missing, one explicit function-tool capability probe precede planning. There are
zero HTTP retries. `parallel_tool_calls=false` enforces one tool request per step.
The total-token setting is a pre-request stop threshold; an in-flight completion
can exceed it. The per-completion token allowance and request timeout remain fixed.

Use `--mode replay` for a real model to investigate the recorded public aggregates.
This is replay rather than a new scientific evaluation. An invalid model response
ends the attempt; it is not silently retried until successful. Reports retain only
allowlisted failure categories, finish reasons, call counts, and known tool names,
never raw provider messages or private argument bodies.

To reuse a promoted policy, pass its content-addressed policy file and promotion
receipt to `--mode replay`:

```sh
uv run --frozen python -m escape_lab.activation_diagnostic_cli \
  --mode replay --root . \
  --policy reports/policy-smoke-01/policies/CHALLENGER_HASH.json \
  --promotion-receipt reports/policy-smoke-01/receipts/RECEIPT_HASH.json \
  --output reports/retained-policy-replay-01
```

The initial parent is inferred for version 2. Later versions require
`--parent-policy PARENT_HASH.json`. The loader checks hashes, lineage, the exact
frozen referee, minimum heldout cases, improvement, false positives, and acceptance.
Content hashes detect mutation; they are not signatures authenticating a party that
can replace an entire artifact set.

## Integration contract

`ActivationDiagnostics` binds numerical capabilities to a trusted evaluator's
current `inner-development` partition or generated `synthetic` cases. It copies and
freezes labels and arrays locally. The planner receives sanitized aggregate outputs.
Use `FGFR2SequenceModel.prediction_outputs()` once and pass its probability and
ranking outputs, plus its aggregate `prediction_diagnostics`. Biology scaling
statistics are training-derived; no outcome is used to fit preprocessing. A bounded
ablation callback may return only combined, WT-only, and mutation-delta AP.

`load_development_replay()` reads exactly the two published activation report files
and checks consistent fold/count summaries and baseline AP. Aggregate summaries do
not establish row identities; both source hashes are retained. Its explicit scope is
`fixed-development-replay`, and its metric is mean-fold high-activation AP. Diagnostic fractions in the aggregate replay are arithmetic means across fits; counts sum across scored rows and extreme magnitudes use the worst observed value. Missing diagnostics remain missing.

`ExperimentPlanner(client, protocol="diagnostics-v3")` opts a new resistance
controller into actionable rejection feedback and explicit stopping. It must remain
inside the caller's permitted development partition. This increment adds no command
that consumes or resets an external or final evaluation.

Focused tests cover sanitizer boundaries, immutable input snapshots, actual Crusoe
HTTP payloads via a mock transport, stop behavior, parser failure audits, policy
hash corruption, forged acceptance flags, duplicate heldout cases, frozen referee
thresholds, receipt-required reuse, and immutable run directories.
