# Tech Stack & Setup

## Stack

| Layer | Choice | Why |
|---|---|---|
| Backend | Python, FastAPI | Fast to build, good async support for multiple external API calls |
| Frontend | Streamlit (fastest) or React | Streamlit if UI needs to be built quickly with limited frontend time; React if you want the tier-distinction styling in `ui-flow.md` to look polished |
| Graph | NetworkX (in-memory) or Neo4j free tier | NetworkX is zero-setup for demo scale; use Neo4j only if the team already knows it |
| Embeddings/RAG | sentence-transformers | Free, local, no API cost |
| LLM | Whatever API access you have | Must support structured/constrained output for the citation-enforcement rule in `api-contracts.md` |
| Mutation parsing | `hgvs` (PyPI), Mutalyzer API | |
| Structural prediction | AlphaFold DB lookups + RCSB PDB (primary); ColabFold only as last resort, never on demo path | |
| Binding affinity | mCSM-lig (web server, interactive use) | |

## Environment variables (placeholder names — adapt to actual implementation)

```
LLM_API_KEY=
PUBMED_API_KEY=            # optional, raises rate limit
CTGOV_API_BASE=https://clinicaltrials.gov/api/v2
CTRI_ACCESS_MODE=          # confirm actual access method first — see data-sources.md
UNIPROT_API_BASE=https://rest.uniprot.org
PUBCHEM_API_BASE=https://pubchem.ncbi.nlm.nih.gov/rest/pug
```

## Setup order

1. Clone/set up backend skeleton (FastAPI project, dependency install).
2. Test each external API in `data-sources.md` with a single manual curl/script call — confirm it actually responds before writing pipeline code around it. Do this for CIViC, ClinicalTrials.gov, CTRI, PubMed, UniProt, PubChem, AlphaFold DB, and mCSM-lig specifically, since several have documented uptime issues.
3. Locate the RareCure repository (check the paper's data-availability/methods section for the exact link) and confirm which of its six modules actually run in your environment vs. which are flagged as architectural/pending validation.
4. Build Tier 1 first, end-to-end, on a small pre-loaded CIViC subset (a handful of cancer types).
5. Wire up the intake form → Tier 1 results screen — get one full demo working before starting Tier 2.
6. Add Tier 2 pipeline per `tier2-structural-prediction.md`, testing each step (parsing → validation → AlphaMissense → structure → candidates → mCSM-lig) independently before chaining them.
7. Run `validation-plan.md` and lock in the two demo cases.
8. Pre-compute and cache both demo cases' full results (all API calls resolved, stored locally) so the live demo never depends on external API latency or uptime.

## Colab free-tier notes (if ColabFold is used as a last resort)

- Disconnects after ~90 min idle or ~12h total session.
- GPU quota is limited and shared — do not rely on it being available at demo time.
- Never run a first-time fold live on stage; if ColabFold is used at all, it should only ever populate the cached demo-case data ahead of time.

## Rate limiting to implement client-side

- PubChem PUG REST: ≤5 requests/sec, ≤400/min, ≤300s running time/min.
- mCSM-lig: no documented programmatic rate limit — throttle conservatively and avoid bulk scripted submissions; confirm acceptable use before automating.
