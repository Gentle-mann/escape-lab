> Historical research record. See [current results](../../RESULTS.md) for the release interpretation.

# You.com evidence pilot

This pilot tests whether reviewed literature context changes GLM's selection from
the existing numerical-model catalog. It does not add biological features or
claim a new biological model. Continue from the original v2 3,730-variant pool,
not the separate mechanistic pilot's altered cohort.

You.com Search runs the four frozen general-kinase queries in
research/you-evidence-plan.json with a primary-source domain allowlist and a
2023 publication cutoff. Contents retrieves at most four selected URLs. Raw
results remain quarantined on the research/evaluator side. Provider date filters
are not trusted as proof: a curator checks source dates and exact supporting
passages. The original cellular study, S13 study, RLY panel and mutation/drug
response claims are excluded from the released context.

Only four predeclared general claims can be released: kinase conformational
plasticity, local ligand-contact chemistry, solvent environment and the limits
of interpreting static structures as kinetics. A retrieved paper must support
the claim or the card is omitted. Claims cannot contain measured responses,
variant-level effects, previous winning configurations or selection instructions.
Short supporting excerpts and raw research are kept in a separate provenance
ledger; the proposing agent receives only canonical claims and reviewed citations.
This tests curated, literature-supported context, not unrestricted web browsing.

Control and research arms use fresh GLM-5.3 clients and planners, with three
proposal opportunities each and one format repair permitted per opportunity.
The same revised prompt explicitly permits optional reviewed evidence; control
receives an empty card list. All original catalog, feature, preprocessing, seed
5100, inner folds, numerical models, promotion rules and metrics stay fixed.
Each arm gets a 90,000-token stop threshold, 600 seconds, max_tokens8192, a
60-second request timeout and zero HTTP retries. These are equal small-pilot
budgets, not a claim to reproduce all four searches in the earlier experiment.
Control runs first, then research. Inference sampling/reasoning defaults remain
unspecified as before. Numerical seeds do not fix LLM randomness.

Both arms' selected configurations, model files, predictions, evidence and
protocol must be jointly hash-committed before either arm is scored. Each arm
must evaluate at least one proposal; otherwise report an unavailable comparison.
Use only already released S13 measurements through the existing retrospective
scorer. Original one-shot receipts remain untouched. Do not evaluate the old
v1 final set or open any fresh/reserved panel. No source retrieval occurs during
selection, and neither arm receives the other's proposals or control/grid scores.

Primary pilot endpoint is research minus control external mean per-drug Spearman;
also report each drug, the four deterministic baselines, all selections, failures,
requests, tokens and costs. Inner AP is selection feedback, not held-out accuracy.
Use 2,000 paired whole-block bootstrap draws, seed6100, conditional on committed
predictions. One stochastic pair on a small reused panel cannot establish causal
or general improvement. A length/format failure consumes its original opportunity.
Preserve failures and stop after one paired attempt; no outcome-driven retries.

Scope all changes to new modules, tests and reports. Verify original frozen hashes,
baseline prediction reproduction, prompt isolation, unsupported-claim rejection,
joint commitments and immutable receipts. Record actual retrieval and inference
usage, including failed requests. Missing credentials block live retrieval and
the real pilot; they must never be replaced by simulated evidence.
