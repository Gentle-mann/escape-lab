> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# Inference-only repair before the second run

Prior run: `20261017T193042Z-7992cb52`.

Both Nemotron requests failed before producing a predictor configuration. The
94.79-second run therefore supplies no evidence for or against agent improvement.
The original report is preserved. Its `no_agent_validation_improvement` outcome
and agent validation row should be interpreted as an operational failure with an
inherited baseline, not a scored agent proposal.

The second run explicitly selects `nvidia/NVIDIA-Nemotron-3.5-Super-VL`, which
previously passed the tool-calling smoke test, and increases the HTTP timeout
from 30 to 60 seconds. The first run used `nvidia/NVIDIA-Nemotron-3-Super-120B-A12B`
because automatic selection depended on catalog order. Failure diagnostics now
retain only safe categories. An incomplete proposal sequence skips validation.

These repairs do not alter scientific choices: dataset, labels, splits, feature
space, fixed predictors, two random configurations, model seed 0, search seed
20261017, two agent proposals, prompt, tool schema, metrics, promotion rules,
final success threshold, and 600-second run cap are unchanged. No first-run
validation scores are supplied to Nemotron. Validation reuse is not independent
replication. The final test is still sealed before this retry.

Declared retry command:

```sh
uv run escape-lab improve --model nvidia/NVIDIA-Nemotron-3.5-Super-VL --request-timeout 60 --max-seconds 600
```

## Completion-budget repair before the third run

Second run: `20261017T193621Z-d6e64c8c`, 71.68 seconds. The endpoint responded
successfully four times, but every response consumed exactly the 2,048-token
completion allowance and no complete valid configuration was returned. This is
consistent with the reasoning model being truncated before its tool call. No
agent configuration was scored; validation and final evaluation were skipped.
Measured provider usage was 22,548 tokens. The second run is preserved as another
operational failure, not a negative prediction result.

The third attempt keeps the same pinned model and 60-second HTTP timeout. It
raises the completion allowance to 8,192 tokens and the total-token stop threshold
to 50,000 so reasoning can finish before emitting a structured tool call. Token
limits are not a hard billing ceiling. The attempt log now records only safe
finish-reason metadata and treats explicit truncation as a failure without an
identical format-repair call. No raw private reasoning is retained.

The two-configuration search budget and all scientific choices listed above are
unchanged. The extra inference budget and two failed attempts will be disclosed;
no unsuccessful run is discarded and no final outcome has been used for tuning.

Declared third command:

```sh
uv run escape-lab improve --model nvidia/NVIDIA-Nemotron-3.5-Super-VL --request-timeout 60 --max-completion-tokens 8192 --max-seconds 600
```

## Schema compatibility repair before the fourth run

Third run: `20261017T193915Z-ee008058`, 85.97 seconds, 27,095 measured provider
tokens. Increasing the completion allowance allowed complete tool calls: all
four responses finished with `tool_calls` and used 1,519–4,328 output tokens.
However, no proposal passed type validation, so validation and final evaluation
were skipped again. The original artifacts remain unchanged.

The synthetic-only probe in `../timing/synthetic-schema-compatibility.json`
reproduced the transport issue. With the original `$ref` schema, `config` arrived
as a string containing JSON. Inlining the exact referenced definition caused
the endpoint to return a configuration object that passed strict validation.
The probe used fabricated aggregate scores and did not train or select a model
on the real data. Its candidate choices are not adopted by the real runner.

The fourth run uses this inlined schema, retaining every original field, enum,
required key, and prohibition of extra fields. All fields remain explicitly
required, matching the existing strict local completeness requirement. The tool schema hash changes;
the prompt and scientific search space do not. No permissive string conversion
is added. The model, token/time allowances, two proposed configurations, data,
seeds, labels, metrics, controls, and acceptance thresholds remain unchanged
from the third attempt. No validation score is sent to the model.

The same declared third command above is used after this repair. Three failed
real runs and three synthetic calls are preserved separately from the subsequent
scientific comparison. They are not independent replications.

## Completed comparison

Fourth run `20261017T194348Z-5e66cbab` completed both agent proposals and the
one-shot final evaluation in 100.74 seconds. Its final macro AP was 0.09481,
below the fixed baseline's 0.10923. The predictive improvement goal failed.
See the [results summary](results-summary.md) for the complete comparison,
runtime/usage accounting, and frozen evidence. The final test is now spent.
