> Historical research record. See [current results](../../RESULTS.md) for the release interpretation.

# You.com research phase

This integration uses You.com Search to discover papers and Contents to retrieve
selected sources. Research text stays in an evaluator-owned quarantine. The GLM
proposal agent receives only source-linked, reviewed general biology cards.

Set `YDC_API_KEY` in this project's ignored `.env` file. `YOU_API_KEY` is also
accepted. Never put a key in a command, source file, report or chat message.
The key needs Search and Contents access. Existing Crusoe configuration is used
for inference. Importing these modules does not make paid calls.

From the Escape Lab directory, run one frozen retrieval plan into a new directory:

```sh
uv run python -m escape_lab.you_research search \
  --plan research/you-evidence-plan.json \
  --output research/you-evidence/RUN_ID --env .env
```

Inspect `search-summary.json` as an evaluator. Select at most four eligible URLs
that identify general scientific sources published by 2023. Avoid assay-source,
reserved-panel and supplemental outcome pages. Save the exact observed URLs as
a JSON array, then retrieve them once:

```sh
uv run python -m escape_lab.you_research contents \
  --research-dir research/you-evidence/RUN_ID \
  --urls research/you-evidence/RUN_ID/selected-source-urls.json --env .env
```

The retrieval files contain untrusted text. Do not pass the directory, snippets,
titles or full papers to the proposal agent. API filtering is not a proof that a
paper contains no test-relevant facts.

A trusted curator creates `reviewed-claims.json` with schema
`reviewed-research-claims-v1`, the absolute `research_dir`, `plan_sha256`,
`content_sha256`, and one to four `cards`. Every reviewed card has exactly:

- `id` and `claim`, matching `research_evidence.CLAIMS` exactly.
- `source_url`, a canonical PMID, PMCID or PDB URL actually returned by Search
  and successfully retrieved by Contents.
- `publication_year` and `publication_year_support`, checked against the source.
- `supporting_quote`, an exact supporting excerpt of at most 25 words.

Use a different source for each card. Exact text matching verifies provenance;
the curator must separately verify that the paper supports the general claim.
Unsupported cards must be omitted. Then compile the immutable dossier:

```sh
uv run python -m escape_lab.research_evidence \
  research/you-evidence/RUN_ID/reviewed-claims.json \
  research/you-evidence/RUN_ID/evidence.json
```

Only after successful real retrieval and review, execute the paired pilot once:

```sh
uv run python -m escape_lab.research_pilot --root . \
  --evidence research/you-evidence/RUN_ID/evidence.json
```

The pilot runs a fresh GLM control and research arm, three opportunities each.
Both arms use the same clarified prompt, catalog and settings. Both prediction
sets must be jointly frozen before the existing retrospective S13 scorer runs.
Read [the protocol](you-research-pilot-protocol.md) for budgets, endpoints and
limitations. An unavailable or degraded arm is reported; do not rerun to choose
a better result.

You.com's published base rates are $0.005 per Search call and $0.001 per Contents
page. Four searches plus four pages therefore estimate $0.024, excluding GLM and
any failed-call billing differences. Actual request counts and costs are recorded
alongside each phase. No new billing or account consent is automated.
[Official billing](https://you.com/docs/administration/billing),
[Search API](https://you.com/docs/api-reference/search/v1-search), and
[Contents API](https://you.com/docs/api-reference/contents).

Use synthetic tests to verify the integration without API calls:

```sh
uv run pytest -q tests/test_you_research.py tests/test_research_evidence.py tests/test_research_pilot.py
```
