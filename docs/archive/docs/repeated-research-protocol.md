> Historical research record. See [current results](../../RESULTS.md) for the release interpretation.

# Repeated research and selector experiments

This is one predeclared batch of two experiments on reused development data. No
old final, S13, RLY or MANO outcome is scored. All original frozen modules and
previous run artifacts remain unchanged.

Experiment 1 compares control (no cards), relevant (the two previously reviewed
You.com cards), and placebo (two truthful, unrelated biological cards). Placebo
describes bee communication and leaf gas exchange, with dates and support checked
against retrieved sources. Raw pages, excerpts, arm names and review commentary
are never sent to proposal agents. Relevant/placebo packets use the same schema,
two citations and approximately matching serialized character counts. Exact model
token lengths are not controlled. Additional source-acquisition attempts are
preserved; an incomplete source extraction is never treated as dated evidence.

Use the original v2 3,730-variant, 15-aligned-block development cohort and original
three outer folds (salt outer-20261017). Numerical fitting seeds are 4100, 4101,
4102; inner folds use the original inner-seed rule. Ten fresh planner repeats per
fold and arm produce 90 isolated searches. Repeats do not change numerical folds
or fitted results and are not new biological replications. Schedule matched
triplets with seed 7300 and rotate arm order by repeat plus fold index. Up to six
searches run concurrently; every search owns its client, planner and audit logs.
Provider LLM sampling defaults remain unspecified, so matching numerical seeds
does not imply identical LLM random draws.

Each search has the same twelve-configuration catalog, common chemistry-C0.1
reference, three new-proposal opportunities and original acceptance rules. Only
its reference and its own queried inner-CV aggregates can enter feedback.
Evaluator-owned numerical caches avoid refitting identical models, but neither
their contents nor other selectors' scores are supplied to agents. Public
summary fields are requested to be at most 240 characters. The output ceiling
increases equally from 8192 to 16384 tokens; this mitigates, but does not guarantee
elimination of, truncation. A search has a 120000-token stop threshold, 600-second
deadline, 120-second request timeout, zero HTTP retries and at most one format
repair per opportunity. Failures consume their original opportunity. Never rerun
a search because of its scientific result. Every search retains its incumbent,
including the original reference when all proposals fail, and enters analysis.

Experiment 2's primary agent is the no-evidence control, chosen in advance. It
compares against:

- Fixed structural C0.1: no adaptive search and no LLM cost.
- Random3: sample three distinct nonreference configurations uniformly without
  replacement using seed 7100 + 100*repeat + fold index.
- Numerical3: evaluate structure C0.1, then interactions C0.1, then C0.01 on the
  currently accepted incumbent's representation. It observes only these three
  evaluations and the shared reference. This tests this specific coarse-to-fine
  rule, not every possible numerical optimizer.
- Full grid: eleven nonreference evaluations, reported only as a larger-budget
  reference, never described as budget matched.

Agent, random and numerical search have the same three evaluation opportunities;
agents additionally incur inference cost. Fixed structural selection has a cost
advantage. Cache-amortized wall time is not a fair stand-alone optimizer cost.
Record actual requests, tokens, failures and selected configurations.

Finish every search, fit selected configurations using each outer training
partition only, and commit all unlabeled outer predictions, model files, search
logs, source/data hashes, protocol, jobs, cards and provenance in one manifest.
The scorer verifies all 90 jobs, all 10x3x7 selections, full two-drug prediction
coverage and every fitted-model hash before requesting outer labels. The existing
trusted load_labels reads a monolithic archive internally but returns only the
requested partition. This is logical separation; it does not claim filesystem
isolation of the trusted evaluator. Remote agents have no filesystem or web tools.

Primary outcome is the mean of ten repeat-level pooled out-of-fold macro average
precisions, each averaging the two drugs. Report per-drug AP, Brier, per-repeat
metrics, selected models and failure rates. Inner scores are selection feedback.
The five primary contrasts are relevant-control, relevant-placebo, control-fixed,
control-random and control-numerical. Evidence-agent versus other baselines is not
a substitute for these predeclared contrasts.

Use 10000 paired crossed bootstrap draws (seed 8100): independently resample whole
aligned groups and repeat IDs, sharing both draws across every method and drug.
Recompute repeat-level AP then average; do not pool duplicate repeat rows as
independent observations. Report 95% descriptive and 99% five-comparison-adjusted
percentile intervals, plus repeat-only dispersion. Intervals condition on fitted
predictions and do not include refitting uncertainty or independent-study variation.

A required contrast passes only when mean absolute AP gain is at least 0.002,
its 99% lower bound exceeds zero, neither drug's AP decreases by more than 0.01,
and each drug's mean Brier is no worse than comparator*1.05+0.001. These are predeclared
engineering criteria with approximate conditional intervals, not proof of
generalization. Research must pass both relevant contrasts; the control agent
must pass all three selector contrasts. Report positive but insufficient gains
as inconclusive, and report every directional difference. Undefined metrics or
intervals cannot pass. Repeated positive values on old development data are not
an untouched external confirmation.

Run once after synthetic checks and independent gate review. Freeze the exact
implementation and evidence hashes before the first scientific model request.
Preserve failures and report both experiments even if neither improves. The old
development cohort is not a substitute for the separately specified training
exclusions required by a future RLY/MANO transfer evaluation.
