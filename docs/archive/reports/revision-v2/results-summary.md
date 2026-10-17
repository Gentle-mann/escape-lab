> Historical research record. See [current results](../../../RESULTS.md) for the release interpretation.

# FGFR v2: biological features improved one external result; agent advantage remains unproven

Run: `20261017T201321Z-a29b17a8`, completed 17 October 2026.

We built and tested the revision. The fixed structural model had a substantially
better point estimate on an independent biochemical panel, but essentially no
gain on internal cellular-assay holdouts. Nemotron did not beat the fixed
structural model or the matched random search. All uncertainty intervals for
the predeclared comparisons include zero. This is a mixed research result,
not a demonstrated autonomous improvement in generalization.

## What we learned about the original failure

The original agent changed logistic regression's C from 1 to 10, weakening
regularization. Its development coefficient norms increased from approximately
8–10 to 20–23, and Brier probability error worsened in every development fold for
both drugs. The AP-only selector accepted it anyway. The original unseen final
macro AP then fell from 0.10923 to 0.09481.

Those development diagnostics support unstable extrapolation and an inadequate
selection rule. They do not establish a unique biological cause. The original
features also omitted drug-binding geometry, used unscaled mixed features, and
relied heavily on position/sequence categories. Regional prevalence differences
and the very small number of independent mutation blocks made selection noisy.

Evidence: development diagnosis (`v1-development-diagnosis.json`; local artifact),
[original result](../improvement/results-summary.md).

## What changed

- Drug-specific geometry from the experimental FGFR1 structures
  [7WCL, pemigatinib](https://www.rcsb.org/structure/7WCL) and
  [6MZW, futibatinib](https://www.rcsb.org/structure/6MZW): ligand distances,
  contact counts, solvent exposure, packing, and distances to kinase motifs.
- Explicit substitution chemistry, four gene indicators, and twelve defined
  chemistry–geometry interactions. The three model representations have 56, 87,
  and 99 numerical columns. No conservation embedding or protein-language model
  feature is implemented.
- Structure sequences aligned to the exact assay references. Missing coordinates
  are flagged and imputed from training partitions only. Engineered substitutions
  and mismatched template side chains are masked. FGFR2–4 geometry remains a
  homologous FGFR1 proxy, not an experimentally observed mutant pose.
- Training-only imputation and scaling; stronger regularization options and
  shallow boosting; no class reweighting or post-hoc calibration.
- A catalog-based Nemotron planner with public, falsifiable experiment hypotheses,
  numerical feedback, retained rejected proposals, and fresh memory in every
  outer fold. It cannot access observations, change metrics, or execute code.
- Selection guards for per-drug AP, Brier probability error, consistency across
  inner folds, and worst-fold regressions. These are engineering safeguards,
  not significance tests.

The research rationale and frozen thresholds are in
[the research and protocol document](../../docs/revision-v2-research.md).
Relevant methodological evidence includes
[selection overfitting](https://www.jmlr.org/papers/v11/cawley10a.html),
[nested evaluation](https://doi.org/10.1186/1471-2105-7-91),
[AIDE](https://arxiv.org/abs/2502.13138), and
[AIRA](https://arxiv.org/abs/2507.02554).

## Evaluation boundary

The original final test remains spent and was not reused. V2 uses original
development and validation examples only, then excludes all nine aligned blocks
containing the external variants across all four FGFR genes. This leaves 3,730
variants in 15 groups. The v1 comparator is retrained on exactly the same eligible
data; its old 0.10923 final score is not directly comparable to this new task.

Three outer grouped folds each contain a fresh three-fold inner search. Controls
are fixed v1, fixed structural logistic C0.1, random search with three new trials,
and all twelve numerical candidates. Nemotron gets three proposal opportunities
per search. The exhaustive grid is a larger-budget control. Shared caching
prevents interpreting total wall time as a per-method speed comparison.

The external source is [Saifudeen et al., Nature Biotechnology 2026](https://doi.org/10.1038/s41587-026-03090-8),
Supplementary Table S13. Twenty-five single mutants have WT identities matching
the assay references. Doubles, a fusion, and the WT-mismatched FGFR2 M420I were
excluded before reading numerical outcomes. All 25 values per drug were observed.

The new assay measures residual biochemical kinase activity, not cellular
enrichment, patient response, or mutation emergence. Its primary metric is
per-drug Spearman rank correlation, averaged across the two drugs. Raw residual
activity includes intrinsic gene-specific potency as well as mutation effects.
WT subtraction was disabled because matched ATP and drug concentration conditions
were not established. No external classification threshold, AP, or Brier was
invented. Geometry covered 25/25 pemigatinib and 21/25 futibatinib examples;
missing geometry remained flagged.

## Results

Higher AP and correlation are better; lower Brier is better. AP and correlation
are different metrics and must not be compared to one another as percentages.

| Method | Internal macro AP | Internal Brier | External mean rank correlation |
|---|---:|---:|---:|
| Original fixed model, retrained | 0.03441 | 0.03537 | −0.01806 |
| Fixed structural model, C0.1 | 0.03463 | 0.03578 | **0.40568** |
| Three-trial random search | **0.03494** | **0.03459** | 0.28374 |
| Twelve-candidate grid | **0.03494** | **0.03459** | 0.18260 |
| Nemotron selection | 0.03052 | 0.03517 | 0.20892 |

Internal defined labels: pemigatinib 3,707 examples /175 positives; futibatinib
3,695 /76 positives. Their prevalence-based constant-score AP averages 0.03389,
so all internal methods are close to an uninformative ranking reference. The
internal results are exploratory because these development data were already
used in v1.

| Method | External pemigatinib correlation | External futibatinib correlation |
|---|---:|---:|
| Original fixed model | 0.25428 | −0.29041 |
| Fixed structural model | 0.29890 | **0.51246** |
| Random search | **0.47392** | 0.09356 |
| Grid search | 0.41232 | −0.04712 |
| Nemotron selection | 0.10348 | 0.31435 |

The fixed structural model's advantage is largely in futibatinib. This comparison
changes representation and preprocessing together; it does not isolate a causal
effect of any one feature. It is also based on a small selected cross-assay panel.

## What the agent actually did

For outer folds 0, 1, and 2, the agent retained chemistry C0.1, selected structure
C0.1, and retained chemistry C0.1 respectively. Eight of nine outer-search
proposals were rejected. The full-data search accepted structure C0.1, then
interactions C0.1, then structure C1.

The last step illustrates the remaining problem. Inner pooled AP increased from
0.04178 for fixed structure C0.1 to 0.04274 for structure C1. Yet the independently
measured external correlation fell from 0.40568 to 0.20892. The selected models
differ only in regularization here. Better inner AP still failed to identify the
better transferable model, even with the new guards.

The numerical grid selected structure boosting depth 3 on the full pool; random
selected interaction boosting depth 2. Neither beat the fixed structural model
externally. This implicates the relationship between selection evidence and the
target assay, in addition to the agent's limited search choices. It is not an
inference-service failure: all twelve proposal opportunities yielded evaluated
configurations, including two successful format repairs.

## Uncertainty and claim limits

Predeclared analysis used 2,000 paired resamples of entire aligned blocks, seed
6100. Every replicate was defined. These intervals condition on the committed
predictions; they do not refit models, correct for repeated development use, or
represent independent biological replicates.

| Contrast | Internal AP delta [95% percentile interval] | External correlation delta [interval] |
|---|---|---|
| Fixed structure − original | +0.00021 [−0.00927, +0.00837] | +0.42374 [−0.08953, +0.69897] |
| Agent − original | −0.00390 [−0.01265, +0.00392] | +0.22698 [−0.31510, +0.61792] |
| Agent − fixed structure | −0.00411 [−0.01086, +0.00146] | −0.19677 [−0.34629, +0.02050] |
| Agent − random | −0.00442 [−0.01109, +0.00005] | −0.07483 [−0.35861, +0.28651] |

All intervals include zero. We cannot claim robust generalization improvement
or an agent advantage. Public hotspot knowledge and model pretraining exposure
cannot be excluded. Publication discovery exposed some qualitative snippets;
the numerical workbook cells used here were withheld until predictions were
committed. This is archival validation infrastructure, not a newly discovered
cancer treatment.

## Runtime, verification, and artifacts

- Live run: **243.55 seconds**, approximately four minutes, excluding research,
  feature preparation, tests, and uncertainty analysis.
- **12/12** evaluated proposals; **14** completion requests including two format
  repairs; **148,958 provider-reported tokens**, none estimated.
- Full suite: **262 tests passed**; Ruff passed. Tests include malformed/failed
  proposals, missing measurements, train-only preprocessing, grouped separation,
  prediction commitments, source tampering, and one-shot receipts.
- Frozen v1 source/data/models remained hash-identical. V2 binds 28 runtime/data
  artifacts and ten fitted model files. All five prediction sets were committed
  before the exclusive external receipt was reserved and measurements were read.
- No runtime source or scientific policy was changed after scoring. New test
  fixture and reporting changes do not alter the frozen experiment.

Artifacts under the completed run directory (`20261017T201321Z-a29b17a8/`; local artifact):
protocol (`20261017T201321Z-a29b17a8/protocol.json`; local artifact),
[results](20261017T201321Z-a29b17a8/results.json),
uncertainty (`20261017T201321Z-a29b17a8/uncertainty.json`; local artifact),
agent experiment history (`20261017T201321Z-a29b17a8/full-search/search.json`; local artifact), and
external commitment (`20261017T201321Z-a29b17a8/external-freeze.json`; local artifact).

## What to carry forward

Keep the fixed structural model as the most promising observed candidate, with
fresh validation still required. Keep the agent as an experimental planner;
do not market it as improving predictions. The next scientific priority is a
larger collection of matched cellular and biochemical datasets, with entire
studies held out. That would let the selection objective reward actual transfer
instead of repeated small regional gains in one screen.

The separate RLY-4008 IC50 panel remains numerically unopened. It could support
another predeclared replication, but is also small and selected. Neither the old
final receipt nor the newly consumed external receipt should be reset. The v2
runner rejects another execution against this spent panel.
