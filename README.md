# Escape Lab

**An autonomous research workflow that proposes experiments, checks its predictions, and keeps a record of what failed.** Built for the Recursive Agents Hackathon on Crusoe Managed Inference.

The long-term ambition is to anticipate cancer's escape routes before resistance emerges. This prototype tests a smaller, falsifiable question using published FGFR assays: can research agents improve numerical predictors under an independent evaluation contract?

**Preclinical research only. No patients, medical advice, treatment recommendations, or demonstrated cancer cure.**

## Latest measured result

On 668 FGFR2 variants across three fixed development folds, ranking the same sequence classifier by its raw decision score instead of a rounded probability improved mean-fold high-activation average precision:

| Ranking | Mean-fold AP ↑ | High-activation hits in 68 picks | Strong-activation hits |
|---|---:|---:|---:|
| Original sequence classifier probabilities | 0.09732 | 1 / 26 | 0 / 7 |
| Same classifier, raw decision margin | **0.11135** | **4 / 26** | 0 / 7 |
| Historical biological tree reference | 0.07186 | 3 / 26 | 0 / 7 |

AP measures ranking quality; **0.11135 does not mean 11.1% classification accuracy**. This was a numerical ranking repair on reused development folds. It resolved probability-rounding ties without changing fitted coefficients or adding biology. The descriptive paired interval includes zero, feature extrapolation remains a limitation, and this is not an autonomous LLM-discovered improvement. Some earlier ridge models had higher AP but severe continuous-value extrapolation problems; all are reported.

The drug-resistance experiments are a **different endpoint**. Neither Nemotron nor GLM-5.3 established superiority over the numerical baselines. A You.com evidence pilot produced identical predictions with and without research context. Negative results are part of the release.

- [Current results and interpretation](docs/RESULTS.md)
- [Complete experiment ledger, including failed attempts](docs/EXPERIMENTS.md)
- [Archived protocols, reports, aggregate JSON and source checksums](docs/archive/manifest.json)
- [Reproduction and data prerequisites](docs/REPRODUCING.md)
- [Hackathon demo and submission checklist](docs/HACKATHON.md)

## How it works

```mermaid
flowchart LR
    R[Reviewed literature via You.com] --> A[Nemotron or GLM on Crusoe]
    A -->|bounded proposal| N[Numerical training engine]
    N --> E[Independent evaluator]
    E -->|permitted development feedback| A
    E --> V[Accept or reject and record evidence]
    D[Private assay outcomes] --> E
```

The agent proposes validated model configurations. Numerical tools perform fitting and prediction. The evaluator owns labels, split grouping, metrics, budgets and promotion rules. Commitments precede held-out scoring; rejected proposals and operational failures remain recorded. Raw outcomes do not enter proposal prompts. The separate replay vault enforces authenticated, budgeted measurement reveals.

These are application-level controls. The prototype does not deploy NemoClaw/OpenShell or claim operating-system isolation from an adversary with filesystem access.

## Run the code checks

```sh
git clone https://github.com/Gentle-mann/escape-lab.git
cd escape-lab
uv sync --frozen --python 3.12
uv run --frozen pytest -q
uv run --frozen ruff check src tests scripts
uv run --frozen escape-lab --help
```

Tests use synthetic fixtures; they do not require API keys or consume scientific holdouts. Data-dependent local tests may skip in a clean clone.

For a new local data installation:

```sh
uv run --frozen escape-lab prepare
uv run --frozen escape-lab status
cp .env.example .env
```

Set `CRUSOE_API_KEY` and optional `YDC_API_KEY` locally. Never commit `.env` or private run traces. `prepare` downloads the original eight public MAVEdb score sets; it does **not** rebuild all later activation, sequence and structure experiments. See [reproduction limits and prerequisites](docs/REPRODUCING.md).

Run the local vault service after preparing data:

```sh
uv run --frozen escape-lab serve --port 8017
```

Agent and research runs use paid/credited inference. Do not run final-evaluation commands as a demo shortcut. The original final set and S13 panel are already spent; RLY and MANO outcomes remain reserved. A new checkout does not turn a public, previously used dataset into an independent validation.

## Source map

- `src/escape_lab/`: ingestion, vault/API, numerical baselines, Crusoe agents, research gates, model comparison and biological pilots.
- `tests/`: evaluator, isolation, numerical and agent-contract checks.
- `docs/`: result interpretation, full experiment ledger and curated historical archive.
- `scripts/export_public_results.py`: aggregate-only publication export from the local research workspace.

Local data, workbooks, model weights, row-level predictions, databases and request logs are excluded. The publication manifest records hashes of the original records and sanitized copies. This repo exposes the implementation and aggregate evidence; it is not a complete redistribution of all experimental inputs.

## Attribution and licensing

The primary FGFR source is [Tangermann and colleagues' study](https://doi.org/10.1038/s41588-025-02431-8), distributed through [MAVEdb](https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00001256-a-9). See the [CAGI FGFR challenge](https://genomeinterpretation.org/cagi7-fgfr.html) for experimental criteria. Protein representations use frozen ESM2-35M; structure inputs and limitations are cited in the archived protocols.

Project code is MIT licensed. Documentation and derived aggregate measurement summaries are CC BY-SA 4.0; underlying data, pretrained weights, publications and dependencies retain their own licenses. See [LICENSE](LICENSE), [NOTICE](NOTICE), and [publication scope](docs/PUBLICATION.md).
