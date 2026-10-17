> Historical research record. See [current results](../../RESULTS.md) for the release interpretation.

# Escape Lab Project Plan

## An autonomous research system for cancer drug resistance

Prepared for the Escape Lab team • 17 October 2026

We will build a preclinical research system that tests whether a proposed cancer-drug coverage claim survives real experimental evidence. The first version will investigate published FGFR mutation measurements through a blinded experiment vault, establish a strong numerical search baseline, and then test whether a Nemotron investigator running on Crusoe improves that search.

The long-term ambition is to discover and close cancer's escape routes before a treatment reaches patients. The hackathon contribution is a small, falsifiable part of that ambition: an agent commits to an investigation, receives actual laboratory measurements, recognizes a counterexample, and tests whether a revised investigation procedure works better on held-out cases.

This document records the agreed product, scientific scope, architecture, evaluation protocol, implementation sequence, and demonstration. It is an implementation plan, not a benchmark report. The biological examples below are published measurements; no agent performance numbers are asserted.

### What cancer resistance means

A treatment does not necessarily create resistant cancer cells. A tumor can already contain rare resistant subpopulations that survive treatment and expand. Cells can also acquire genetic changes or enter other biological states that reduce sensitivity. Resistance can involve the drug target, alternative pathways, cell state, the surrounding tissue, or several mechanisms together. One mutation screen captures only part of that biology.

Our first question is narrower: does a variant that resists one drug also resist another in the measured assay? A counterexample can refute a blanket drug-coverage hypothesis. It cannot establish that a patient will develop that mutation, that a combination fails, or that a particular sequence of treatments works.

### The first product

A researcher submits a claim such as “Drug B covers the mutations that allow resistance to Drug A.” Escape Lab produces an evidence-backed resistance audit: verified counterexamples, measured supporting cases, unresolved coverage, and the experiment that would most reduce a consequential uncertainty. The system serves drug-development research. It does not diagnose a patient or recommend a drug, dose, or care plan.

### Decisions already made

- Build data integrity, the experiment vault, and numerical baselines before the language-model loop.
- Let statistical tools rank candidates; let the investigator choose strategies and interpret evidence; let independent code enforce experimental truth and budgets.
- Use real archived measurements for the demo, clearly labeled as retrospective replay.
- Improve on development data, freeze a selected procedure, and preserve a final test that has never informed revisions.
- Keep useful conventional components even if agents fail to improve them. A credible negative result is preferable to a staged win.

<!-- PAGEBREAK -->

## The FGFR experimental foundation

The primary source is “Saturation mutagenesis identifies activating and resistance-inducing FGFR kinase domain mutations.” Its public supplementary workbook contains the sheet `S10_Screen data table`, with 8,407 missense nucleotide-variant rows across FGFR1–4. These are rows, not 8,407 independent tumors or necessarily 8,407 unique protein substitutions. The full screen includes 11,520 single-nucleotide substitutions, including other consequence classes. [1, 2]

The experiment used engineered NCI-H1581 cells. Activation has four biological replicates; the two resistance assays have three each. Pemigatinib and futibatinib were tested separately. The table includes replicate measurements, medians, classifications, significance fields, and existing predictors. Outcome-derived columns must never enter public candidate features. [1, 2]

### Two verified examples

| FGFR2 variant | Pemigatinib median | Futibatinib median | Published interpretation |
| --- | ---: | ---: | --- |
| V564G | 12.2924 | 0.9688 | Strong resistance to pemigatinib; neutral for futibatinib |
| N549K | 44.9882 | 14.1473 | Strong resistance to both drugs |

The values are relative enrichment under selection. They are not percent survival, mutation probabilities, or patient outcomes. The second example refutes the universal claim that the second drug covers every escape variant of the first in this assay. It does not test a drug combination or sequence. “Neutral” is an assay classification, not proof of clinical efficacy. [2]

### Import and provenance requirements

Download the source once, record its URL and checksum, retain the original file privately, and produce a schema-validated manifest. Preserve gene, nucleotide and protein identities, replicate values, source sheet and row, assay, and applicable license. Resolve duplicate protein substitutions deterministically without silently treating duplicate nucleotide encodings as independent discoveries.

Use a fixed, documented resistance definition. The associated CAGI challenge describes median enrichment greater than 1.5 with supporting results in at least two of three replicates. Reconcile that rule with the actual workbook columns and paper before freezing the benchmark. A statistical significance flag alone is insufficient because significant depletion can also be flagged. [3]

The article and its covered supplementary material carry CC BY 4.0 terms, subject to any separate notices. Inspected MAVEdb score sets use CC BY-SA 4.0. Preserve attribution, record the exact source license, and honor share-alike obligations for adapted MAVEdb data. Write our own loader; publication of research code does not by itself grant a reuse license. [1, 4]

### Scope of the evidence

This dataset supports resistance-variant discovery in a limited engineered assay. It cannot validate tumor evolution forecasts, patient digital twins, combination efficacy, or a cure. CAGI already used the measurements for a prediction challenge, so a new heatmap or prediction wrapper is not the project's novelty. [3]

<!-- PAGEBREAK -->

## Existing approaches and the contribution to test

Current resistance research includes molecular profiling, development of drugs for specific resistance mechanisms, repeated measurements during treatment, organoid and cell-line testing, evolution experiments, and adaptive treatment research. Computational models can prioritize experiments, estimate responses, and generate hypotheses. Escape Lab should connect to this work rather than suggest that resistance anticipation is a new field.

ARPA-H's ADAPT program is already pursuing ways to outpace cancer through repeated biological measurement and intervention. PERCEPTION uses single-cell expression information for response prediction and resistance-related analysis. These are important precedents for the larger ambition; they do not establish our proposed system's performance. [5, 6]

| Existing system | Established contribution | Consequence for Escape Lab |
| --- | --- | --- |
| BATCHIE | Bayesian active learning for efficient cancer drug-combination experiments | Compare against conventional statistical acquisition; longitudinal outcomes alone would not establish novelty. |
| BioDiscoveryAgent | Language-model agents select biological experiments and update from observations | An agent choosing the next experiment is already prior art. |
| Co-Scientist | Hypothesis generation, criticism, ranking, and revision with experimental follow-up | Agent debate and scientific role names are insufficient differentiation. |
| STELLA | Persistent scientific workflows and evolving research capabilities | Store reusable procedures, but independently verify whether revisions help. |
| CAGI FGFR challenge | Prediction tasks on the same FGFR measurements | Do not present the dataset, resistance atlas, or basic predictor as our invention. |

Sources: BATCHIE [7], BioDiscoveryAgent [8], Co-Scientist [9], STELLA [10], and CAGI [3].

### The concrete performance bet

Escape Lab tests whether an adversarial, auditable investigation finds verified resistance counterexamples more effectively than a strong numerical active-search policy. It also tests whether persistent revisions to the investigation procedure improve later runs while preserving evidence quality.

The system must distinguish three outputs: an unmeasured prediction, a verified counterexample, and an unresolved question. Every confident verdict must be supported by disclosed measurements and the frozen assay criteria. No amount of agreement among agents can turn a prediction into an observation.

### Architecture lessons from the research

More agents and tools do not guarantee better results. BioDiscoveryAgent reports model-dependent tool effects, including worse performance from the full tool package for a stronger model on several datasets. Multi-agent research architectures can benefit from parallel exploration, while the next experiment in this task depends strongly on the previous measurement. [8, 11, 12]

Use parallel work for independent source gathering. Keep the experimental loop under one controller. Add a skeptic only at decisions where it can prevent a real error, then measure whether its benefit exceeds its cost. Scientific novelty and commercial advantage remain hypotheses until supported by comparative results.

<!-- PAGEBREAK -->

## Hackathon eligibility and demonstration constraints

The participant guide's medical prohibition is explicit: “Any project using AI to generate and give medical advice.” A preclinical claim-testing tool appears compatible because its inputs and outputs concern laboratory evidence. The guide does not explicitly guarantee a research exemption, so this is a defensible interpretation, not organizer approval. [13]

A disclaimer alone is insufficient. The actual product must avoid patient records, personalized treatment choices, diagnosis, dosing, and care instructions. Its interface should ask for an experimental hypothesis and return assay-bounded evidence.

### The full prohibited categories

The live guide lists eleven exclusions. The categories below are paraphrased except for the quoted medical rule above. [13]

1. AI advice for mental health.
2. Basic retrieval-augmented generation applications.
3. Applications built with Streamlit.
4. Image-analysis projects.
5. Education chatbots.
6. AI screening of job applicants.
7. AI nutrition coaching.
8. Personality analysis.
9. AI-generated medical advice.
10. Projects whose central feature is a dashboard.
11. Sports analysis or coaching.

Escape Lab will use numerical data rather than tumor images, a lightweight web interface rather than Streamlit, and a live investigation rather than a dashboard as its main experience. Search supports provenance; the central function is execution and independent verification. The remaining coaching, screening, and personality categories do not describe this product.

### Track alignment

**Done Means Done is the primary track.** The team runs on Crusoe Managed Inference, defines a testable objective, uses an independent verifier, detects invalid conclusions and loops, and decides whether to continue, replan, or stop.

**Version N+1 is earned through measured revision.** The editor may change a playbook, acquisition configuration, or permitted tool. It retains a change only when separate validation supports improvement. Updating a prediction during a run is ordinary learning and does not alone demonstrate persistent self-improvement.

**Shared Brain is optional.** The typed evidence ledger prevents conflicting or stale claims, but a track claim should show a measurable quality benefit from shared memory rather than merely storing messages.

### Submission obligations

The guide requires a public repository, no more than four people, appropriate rights to code and data, and clear identification of the functionality built during the hackathon. Existing work cannot be presented as our own event-built contribution. The schedule allocates roughly six hours of hacking, a one-minute submission video, and a live demonstration of about three minutes. Sponsor examples are not exceptions to the main exclusions. [13, 14]

<!-- PAGEBREAK -->

## Agent architecture and control flow

Escape Lab is a small scientific control system. Language models manage the investigation, numerical tools rank experiments, and separately enforced code determines what the measurements establish. The system has a preparation phase, a sequential experiment phase, and a development-only revision phase.

### Preparation

Use You.com Search and Contents to gather primary sources, assay definitions, and permitted biological features. Store a frozen dossier with exact citations and a manifest. Check it for hidden outcome information before evaluation. Define the candidate universe, query unit, assay criteria, starting observations, budgets, and permitted verdicts in an immutable research contract.

### The experiment loop

`PLAN → RANK → COMMIT → REVEAL → UPDATE → CHECK → CONTINUE or STOP`

The investigator chooses a strategy and a candidate from a numerical shortlist. Before revelation, it records its prediction and rationale. The vault returns real measurements and a stable evidence identifier. The controller updates the evidence ledger and checks the budget, task objective, and termination rules. Tool failures produce explicit recoverable errors rather than fabricated results.

| Component | Responsibility | Authority boundary |
| --- | --- | --- |
| Investigator | Formulate tests, choose strategies, interpret observations, propose a verdict | Cannot read hidden outcomes or edit scoring |
| Numerical engine | Fit models and rank candidates reproducibly | Uses permitted features and disclosed observations only |
| Skeptic | Challenge assay interpretation and proposed stopping decisions | Criticism is not experimental evidence |
| Experiment vault | Reveal authorized archived records and enforce query accounting | Sole source of hidden measurements |
| Referee and controller | Validate actions, criteria, evidence references, and transitions | Independent of agent persuasion |
| Procedure editor | Propose persistent changes from development failures | Cannot alter data boundaries, labels, or evaluator |

### Stopping conditions

A falsification task stops after a qualifying counterexample. A coverage audit continues to discover additional escape variants until its budget or another predefined termination rule is reached. These objectives must be selected before a run and never exchanged midway to improve the apparent result.

The verdict can be “counterexample verified,” “no counterexample found within this search,” or “unresolved because evidence is insufficient.” Exhausting a finite search budget does not prove safety or complete coverage. An always-abstaining system also cannot win: evaluation considers both incorrect confidence and useful conclusions.

The skeptic runs at proposed stopping, scope-changing conclusions, major prediction failures, and procedure revisions. It receives the contract and measurements before persuasive investigator prose where practical, reducing anchoring. Its contribution must pass an ablation against the same system without the skeptic.

<!-- PAGEBREAK -->

## Experiment vault and evidence memory

The vault is the core trust boundary. Candidate metadata belongs on the agent side; hidden replicate outcomes and final evaluation labels remain on the evaluator side. A prompt telling the agent not to read a nearby answer file is not a security boundary. Agent-accessible tools, file mounts, APIs, and code execution must enforce the separation.

### Proposed tool contract

| Tool | Required behavior |
| --- | --- |
| `get_research_state()` | Return the frozen contract, public candidates, disclosed evidence, remaining budget, and procedure version. |
| `rank_candidates(policy_id, limit)` | Return a compact reproducible shortlist, feature provenance, and model estimates. |
| `inspect_evidence(evidence_ids)` | Retrieve existing disclosed measurements or permitted cited passages. |
| `commit_experiment(candidate_id, prediction, reason)` | Validate eligibility, persist the pre-result commitment, and return an immutable commit ID. |
| `reveal_experiment(commit_id)` | Reveal only the committed record, attach source provenance, and apply budget accounting exactly once. |
| `submit_verdict(status, evidence_ids)` | Validate that cited evidence supports an allowed assay-bounded verdict. |

Exact endpoint names may change during implementation; these behavioral contracts must remain. One benchmark query reveals a variant's paired records for both separately tested drugs, with available replicates. It is a paired-record reveal, not a claim of one physical laboratory experiment saved.

### Correctness and recovery

Revelation must be idempotent. A network retry for the same committed query returns the same response without charging again. Reject unknown candidates, absent commitments, duplicate discoveries, over-budget actions, and evidence references outside the run. A single global query budget covers investigator and skeptic actions. Persist ordered events so a run can resume without silently losing commitments or counting a result twice.

### Typed memory

| Record type | Stored facts | Who establishes it |
| --- | --- | --- |
| Observation | Assay, variant, replicate values, source file hash and row, reveal event | Experiment service |
| Literature prior | Claim, exact source, retrieval date, scope, and limitations | Preparation pipeline |
| Prediction or hypothesis | Probability or qualitative claim, model version, timestamp, evidence available then | Numerical model or investigator |
| Procedure version | Parent version, configuration hash, motivating failures, evaluation and acceptance record | Editor proposes; evaluator accepts |

SQLite is sufficient initially. Keep the durable record outside the model context. Give the model a concise current-state view and retrieval identifiers, not an unbounded transcript. Use Pydantic schemas for tool inputs and outputs. Never log credentials or send the private workbook through the inference API. Public provenance and actual agent-visible inputs should be inspectable in the run report.

<!-- PAGEBREAK -->

## Numerical search and baseline design

The main optimization problem is active search: find valuable counterexamples under a fixed measurement budget. Learning a globally accurate model is useful but is not identical to maximizing discoveries in the next twenty queries. [15]

### Candidate features

Start with features available without the hidden outcomes: gene, amino-acid substitution, aligned residue position, predefined kinase-region annotations, and substitution properties. Record the provenance of every external annotation. Do not include medians, significance calls, resistance categories, study-derived rankings, or published predictors whose inputs or training expose the hidden labels.

Represent the joint dual-resistance outcome directly. Multiplying separate probabilities for resistance to each drug assumes independence and is not justified by default. Preserve replicate information for evaluation and uncertainty handling. Define how incomplete, contradictory, or unclassifiable records are treated before scoring.

### Baselines to implement

| Method | Purpose |
| --- | --- |
| Random selection | Establish expected discovery from prevalence and budget alone |
| Gene or known-hotspot ranking | Detect easy gains caused by gene prevalence or familiar mutations |
| Regularized statistical classifier | Provide a transparent supervised baseline |
| Small tree ensemble with active acquisition | Capture nonlinear features and update after observations |
| Frozen tuned active-search policy | Main comparator for the investigator and procedure editor |

Select the strongest statistical configuration through development evaluation. Do not declare a specific estimator best before measuring it. Hyperparameter choices, seed sets, training examples, feature transformations, and candidate restrictions belong in versioned configuration.

### Acquisition strategies

Compare exploitation of estimated hit probability with bounded exploration, sequence-region diversity, and model disagreement. Tune these choices on development data. An LLM can choose among supported strategies, but numerical tools must turn that choice into an auditable ranking. The model's verbal confidence is not a calibrated biological probability.

All methods receive the same initial observations and candidate universe. Account for the cost of any common starting measurements explicitly. Baseline methods may not use a different gene distribution, more measurements, or a cleaner subset than the agent. When aggregating outcomes, count a unique protein substitution once rather than rewarding duplicate nucleotide encodings.

### What would count as useful agent work

The investigator may improve strategy selection after surprising evidence, identify an invalid assay inference, or recognize that the current search is too concentrated. Test each contribution. If the statistical baseline remains better at selecting experiments, retain it and report any separately measured auditing benefit of the agent. A persuasive explanation is not evidence of better acquisition.

<!-- PAGEBREAK -->

## Accuracy evaluation and protection against leakage

Verification has four layers: source-data correctness, pre-disclosure prediction quality, search effectiveness, and scientific correctness of conclusions. A system can perform well at one and poorly at another. We will evaluate each explicitly.

### Metrics

| Metric | Meaning and limitation |
| --- | --- |
| Hits at 20 | Unique verified dual-resistance protein variants found within twenty paired-record queries |
| Queries to first counterexample | How quickly a falsification task exposes a failure |
| Sequence-region coverage | Breadth of discoveries; not automatically diversity of biological mechanisms |
| Precision and recall based metrics | Quality of resistance predictions, especially under class imbalance |
| Calibration | Whether predicted probabilities agree with observed frequencies on a blinded audit set |
| Unsupported conclusions and useful verdicts | Penalize false confidence without rewarding universal abstention |
| Cost and latency | Tokens, inference cost, invalid actions, elapsed time, and recovery overhead |

Ordinary percent accuracy is insufficient when negative outcomes dominate. Prediction quality must also be assessed on a fixed blinded audit set, not only selected queries, because the selection policy changes which outcomes become visible.

### Data partitions

Use development data for modeling and debugging, separate selection-validation tasks for a limited recorded number of revisions, and an untouched final test for confirmation. Group duplicate protein substitutions and their replicate records. Group aligned residue positions across related FGFR genes to reduce leakage through close relatives. Freeze split manifests and document the resulting support and class counts.

During final evaluation, disable web access, inspect the frozen dossier for leaked outcomes, clear cross-run agent memory, and prevent final traces from reaching the editor. Agents may learn from legitimate measurements revealed within their own run. They may not transfer hidden final outcomes to subsequent runs or revisions.

The data are public, so pretraining exposure cannot be ruled out. The defensible claim is runtime-blinded retrospective evaluation. Results from four related genes in one engineered assay do not establish generalization to all cancer or to patients.

### Paired comparison and acceptance

Run every method on matched task instances and seeds. Report paired differences, variation by region and gene, and uncertainty intervals where justified. Repeated seeds measure computational variability, not new independent biological samples. Compare the agent editor with ordinary configuration tuning as well as the frozen predecessor.

After measuring baseline performance, define a minimum worthwhile gain and acceptable latency before tuning the agent or opening the final test. A revision must improve the predefined objective without breaking the evidence and budget invariants. If the final result disappoints, report it. Further iteration requires fresh untouched confirmation data; the same final test cannot be repeatedly optimized into success.

<!-- PAGEBREAK -->

## Sponsor integration and implementation stack

The first implementation uses Python, FastAPI, Pydantic, SQLite, and a numerical modeling library. A lightweight web execution view comes after the experiment loop works. Plain Python is sufficient to own the state transitions. LangGraph is optional if the team already knows it; checkpointing does not remove the need for idempotent tools. [16]

| Technology | Role in this project | Initial decision |
| --- | --- | --- |
| Nemotron | Investigator, selective skeptic, and development procedure editor | Benchmark tool use and judgment before selecting a model |
| Crusoe Managed Inference | Hosted model execution | Required inference path for the agent loop |
| You.com Search and Contents | Source discovery and extraction | Preparation only during scored replay |
| You.com Research | Separate hosted research agent | Optional background support, outside the experimental controller |
| OpenShell | Restricted execution for generated code | Add when agents can execute their own tools |
| NemoClaw | Deployment and lifecycle stack over supported runtimes | Defer unless already working on the team's machine |
| NeMo Agent Toolkit | Workflow evaluation and profiling | Consider after the core loop is reliable |

### Nemotron on Crusoe

Use the OpenAI-compatible Chat Completions interface at `https://api.inference.crusoecloud.com/v1`. Discover exact available model IDs and parameters through `/v1/models` rather than hard-coding names from marketing material. Nemotron 3 Super and a lower-latency alternative are candidates from the researched catalog, not proven winners for this task. [14, 17]

Before scientific evaluation, test schema-valid tool calls, budget adherence, contradictory evidence, error recovery, and honest unresolved verdicts. Validate arguments locally. Keep `CRUSOE_API_KEY` in the project's ignored `.env` file; the user has agreed to add it there. Never print the key, place it in run traces, or commit it. Confirm live access without exposing credentials.

### Execution and evaluation boundaries

If generated Python becomes necessary, run it with permitted features and disclosed observations only, no outbound network, and no access to evaluator files or credentials. OpenShell can enforce these boundaries. Use a compatible pinned NemoClaw and OpenShell stack; mixing unrelated latest versions adds avoidable integration risk. [18, 19]

NeMo Agent Toolkit can wrap functions and profile costs, latency, and tool use. Its custom evaluator must score the vault's records directly. Keep evaluator-only fields out of workflow inputs; settings that pass entire dataset rows can expose labels. The editor must not rewrite the full workflow configuration if that configuration controls splits or scoring. [20]

The system should remain runnable without optional orchestration frameworks. Sponsor integrations should support the measured scientific loop rather than determine its architecture.

<!-- PAGEBREAK -->

## Persistent revision and the build sequence

Ordinary within-run updating changes predictions after a measurement. Persistent revision changes how later investigations operate. The latter is the Version N+1 capability we intend to test.

### Bounded procedure revision

The editor receives development traces and measured failures. It may propose changes to acquisition settings, shortlist diversity, allowed feature selections, or the investigator playbook. Initially use validated configuration patches rather than arbitrary code. Each patch records its parent procedure, motivating failure IDs, changed fields, and configuration hash.

The editor cannot change assay thresholds, labels, splits, query accounting, access permissions, or scoring. Validate the patch, run regression checks, compare on separate selection-validation tasks, and accept or reject under the frozen rule. Keep rejected revisions in the record. A small revision budget helps prevent untracked overfitting. This adapts the execution-feedback idea used by GEPA without requiring its entire optimizer. [21]

### Seven implementation stages

| Stage | Deliverable | Completion evidence |
| --- | --- | --- |
| 1 Data and contract | Verified source import, labels, public/private schema, grouped splits | Source reconciliation, duplicate and missing-data audit, leakage checks |
| 2 Experiment vault | Commit, reveal, budget, evidence, and replay service | Isolation tests, idempotent retries, correct global accounting |
| 3 Numerical baseline | Random, informed ranking, tuned statistical active search | Reproducible grouped development report and saved configurations |
| 4 Nemotron loop | Crusoe client and investigator tool loop | Live model access, valid tools, recovery, correct use of evidence |
| 5 Agent comparison | Matched numerical and agent-assisted runs | Paired performance, cost, latency, and failure analysis |
| 6 Skeptic and editor | Selective criticism and bounded procedure revisions | Ablations and independently accepted improvements |
| 7 Final test and demo | Frozen procedure, untouched evaluation, live execution view | Reproducible report and honest demonstration of measured behavior |

### Testing priorities

Write focused tests for source mapping, duplicate identity, label construction, hidden-data access, pre-result commitments, budget limits, retry behavior, verdict evidence, and final-test isolation. Run the real pipeline end to end. Passing unit tests does not establish scientific effectiveness; the benchmark must provide that evidence.

Reserve the first milestone for the verified FGFR dataset, experiment vault, and numerical benchmark report. Add the agent only after that foundation is reliable. If time is short, cut optional framework integration, generated code, or extra agents before cutting measurement integrity. Implementation status and actual benchmark results should be maintained separately from this plan.

<!-- PAGEBREAK -->

## Demonstration and next research milestones

The audience should see evidence change the investigation. The principal interaction is a running experiment with a hypothesis, committed query, revealed measurement, and checked verdict. Graphs support that sequence; a static dashboard is not the product.

### Three minute live demonstration

| Time | Action | What it demonstrates |
| --- | --- | --- |
| 0 to 30 seconds | Explain the drug-coverage claim and the two published FGFR examples | One encouraging example cannot establish universal coverage |
| 30 to 90 seconds | Start a blinded run; show commitments and actual measurement reveals | Agent decisions are judged by laboratory evidence |
| 90 to 135 seconds | Show the counterexample or an honest unresolved result and evidence trail | Scientific conclusions are constrained by the assay |
| 135 to 165 seconds | Compare the frozen baseline and a validated revised procedure | Persistent improvement is measured, or a failed revision is rejected |
| 165 to 180 seconds | Show the next laboratory integration step and the actual benchmark summary | The prototype is a concrete first step toward the larger research platform |

Keep the known examples separate from blind evaluation. Do not claim those published discoveries, stage a weak initial agent, or select only flattering runs. The one-minute video can compress the same evidence sequence.

### Risks that change implementation choices

| Risk | Response |
| --- | --- |
| Hidden answers leak through files, source passages, or derived columns | Enforce a separate vault and inspect all public inputs |
| Easy gene prevalence produces impressive-looking scores | Include informed gene and hotspot baselines |
| Technical replication is mistaken for biological generalization | Report the actual unit of independence and assay scope |
| Agent prose overstates a measurement | Validate verdicts and keep combination, sequence, and patient claims outside scope |
| Repeated validation or final-test use overfits the procedure | Record revision budgets and preserve untouched confirmation data |
| Public data appeared in model pretraining | State the limitation and pursue new experimental validation later |

### Beyond the hackathon

Next, propose a previously unmeasured experiment and evaluate its real laboratory result under a frozen protocol. Expand to distinct cell systems, longitudinal experiments, and measured combinations or sequences. Replace the archival vault with a laboratory job queue while preserving commitments, provenance, and independent evaluation.

The pitch is: “Cancer evolves around our drugs. We are building an automated research lab that searches for escape routes before a treatment reaches patients. Today we test drug-coverage claims against real measurements and measure whether the system can improve its own investigation.”

<!-- PAGEBREAK -->

## Scientific and event sources

References support the facts and architecture discussed in this plan. Event rules and service capabilities can change; verify them again before submission or integration. The publication examples are not validation results for Escape Lab.

1. **FGFR saturation mutagenesis paper.** Experimental design, assay interpretation, and article license. https://www.nature.com/articles/s41588-025-02431-8

2. **FGFR supplementary workbook.** Primary measured examples and source table. https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41588-025-02431-8/MediaObjects/41588_2025_2431_MOESM3_ESM.xlsx

3. **CAGI 7 FGFR challenge.** Existing prediction task and resistance criteria. https://genomeinterpretation.org/cagi7-fgfr.html

4. **MAVEdb experiment and example score set.** Public resistance records and license metadata. https://www.mavedb.org/experiments/urn:mavedb:00001256-a/ and https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00001256-a-7/scores

5. **ARPA-H ADAPT clinical research program.** Existing efforts to outpace cancer through adaptive measurement and intervention. https://arpa-h.gov/news-and-events/arpa-h-opens-first-interventional-clinical-trials-outpace-cancer-real-time

6. **PERCEPTION.** Single-cell expression models of treatment response and resistance-related analysis. https://www.nature.com/articles/s43018-024-00756-7

7. **BATCHIE.** Bayesian active learning for cancer drug combinations. https://www.nature.com/articles/s41467-024-55287-7

8. **BioDiscoveryAgent.** Agent-based biological experiment selection and tool ablations. https://arxiv.org/html/2405.17631v3 and https://github.com/snap-stanford/BioDiscoveryAgent

9. **Co-Scientist.** Scientific hypothesis generation and revision with experimental follow-up. https://www.nature.com/articles/s41586-026-10644-y

10. **STELLA.** Persistent and evolving scientific workflows. https://github.com/zaixizhang/STELLA

11. **Anthropic multi-agent research system.** Parallel research architecture and coordination. https://www.anthropic.com/engineering/multi-agent-research-system

12. **Agent-scaling study.** Task-dependent benefits and costs of multi-agent coordination. https://arxiv.org/html/2512.08296v2

13. **Hackathon participant guide.** Tracks, exclusions, schedule, team and submission requirements. https://cerebralvalley.ai/e/crusoe-recursive-agents-hackathon/details

14. **Crusoe event developer reference.** Managed inference integration and event setup. https://github.com/crusoecloud/crusoe-developer-hub/tree/main/examples/hackathons/crusoe-cv-recursive-agents-hackathon

<!-- PAGEBREAK -->

## Engineering and evaluation sources

15. **Efficient Nonmyopic Active Search.** Formal search objective under a labeling budget. https://proceedings.mlr.press/v70/jiang17d.html

16. **LangGraph persistence.** Checkpoints and execution recovery. https://docs.langchain.com/oss/python/langgraph/persistence

17. **Crusoe model catalog.** Available managed inference models. https://docs.crusoecloud.com/serverless-inference/available-models/index.html

18. **OpenShell policy documentation.** Filesystem, network, and execution boundaries. https://docs.nvidia.com/openshell/latest/how-it-works/policies/overview

19. **NemoClaw ecosystem and blueprint.** Runtime relationships, lifecycle, and compatible versions. https://docs.nvidia.com/nemoclaw/latest/user-guide/openclaw/about/ecosystem and https://github.com/NVIDIA/NemoClaw/blob/main/nemoclaw-blueprint/blueprint.yaml

20. **NeMo Agent Toolkit evaluation and profiling.** Workflow evaluation, custom evaluators, and data contracts. https://docs.nvidia.com/nemo/agent-toolkit/latest/improve-workflows/evaluate.html and https://docs.nvidia.com/nemo/agent-toolkit/latest/extend/custom-components/custom-evaluator.html

21. **GEPA.** Reflective optimization using execution feedback. https://github.com/gepa-ai/gepa

22. **You.com Search and Contents.** Source retrieval and extraction. https://you.com/docs/guides/search and https://you.com/docs/guides/contents

23. **You.com Research.** Separate hosted research API. https://you.com/docs/api-reference/research/v1-research

24. **Tool design for agents.** Clear schemas, concise outputs, and evaluation of tool usefulness. https://www.anthropic.com/engineering/writing-tools-for-agents

25. **Managed-agent architecture.** Durable session, agent loop, and execution environment. https://www.anthropic.com/engineering/managed-agents

26. **Agent evaluation guidance.** Outcome-based grading and trace analysis. https://www.anthropic.com/engineering/demystifying-evals-for-ai-agents

### Project location

The project directory is ``. This plan is stored in its `docs` directory as an editable Word document and a Markdown source file. Credentials belong only in the ignored project `.env` file, not in either document.

### Reproducibility record to produce during implementation

The implementation should generate a dataset manifest with checksums and licenses; split manifests and label definitions; versioned baseline and agent configurations; immutable experiment commitments and reveals; model and tool versions; paired evaluation reports; costs and latency; and a record of accepted and rejected revisions. These artifacts will establish what was actually built and measured independently of the ambition described here.
