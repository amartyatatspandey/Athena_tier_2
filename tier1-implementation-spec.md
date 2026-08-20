# Tier 1 Implementation Spec — Evidence Retrieval & Knowledge Graph

Project: SecondLook (AI second-opinion copilot for rare and treatment-exhausted cancers, India-localized extension of RareCure). This document is the handoff spec for whoever builds Tier 1. Per `docs/checkpoint.md`, this repo's code so far (`src/secondlook/`) covers **Tier 2 only** — Tier 1 is unbuilt and is a separate scope, though several Tier 2 modules are directly reusable here (noted below). Per `docs/README.md`'s stated build order, Tier 1 is meant to be the **first thing built** — it's "the safer, more defensible V1 core" — even though Tier 2 happened to be built first in this repo; whoever picks this up should not treat Tier 2's completion as implying Tier 1 is behind schedule by design.

Source docs this spec is derived from: `docs/tier 1 docs/tier1-retrieval.md` (the Tier 1 spec itself), `docs/architecture.md` (module boundaries, data flow, orchestration), `docs/api-contracts.md` (exact schemas), `docs/ui-flow.md` (frontend), `docs/data-sources.md` and `docs/tech-stack-setup.md` (stack, endpoints, setup order), `docs/rarecure-build-reference.md` (specific, named failure modes from the prior-art codebase to avoid repeating), `docs/founder-mode-use-case-and-kg-solutions.md` (the use case this graph ultimately serves, and the graph-native capabilities it enables), `docs/checkpoint.md` (current repo state and what Tier 1 depends on from Tier 2).

---

## 1. Scope and task

Tier 1 surfaces **documented, citable** treatment options for a patient's cancer mutation: curated variant-drug evidence (CIViC), matching clinical trials (ClinicalTrials.gov + CTRI for India), and supplementary literature (PubMed). It is the primary path — Tier 2's structural prediction only runs as a fallback when Tier 1 finds nothing or a weak match.

Tier 1 also **owns the knowledge graph**, which is not merely a storage detail for CIViC data: it is the substrate that makes the system's longer-term capability possible (§2). Everything Tier 1 retrieves, and everything Tier 2 computes, lands in the same graph.

Concretely, this scope covers four retrieval components, in build-priority order:

**1. Knowledge graph — CIViC backbone.** Given a validated gene + mutation, traverse a pre-loaded CIViC graph (variant → disease → drug → evidence item, evidence level A–E per CIViC's own scale) to connected drugs and their evidence. Implementation: **FalkorDB** — see §3 for the full schema and §4 for the retrieval design. Pre-load for a defined set of cancer types — include sarcoma specifically, since it ties back to RareCure's own validation cohort and the project's origin story, and is also the intended source for the Tier 2 demo case per `validation-plan.md`.

**2. Clinical trial search.** Query ClinicalTrials.gov API v2 (global) and CTRI (ctri.nic.in, India-specific) by condition + biomarker/mutation, falling back to condition + free-text drug/gene terms where the registry doesn't support biomarker search. Every result must carry trial ID, recruiting status, location(s), and an eligibility-summary link. **CTRI's access mode is unconfirmed** — `data-sources.md` flags it as "check current access mode (API vs. scrape)" and `tier2-verification.md`'s closing note explicitly states it was "not verified in this pass... still flagged unconfirmed." Confirming this is a hard prerequisite before building the India-trial-search path, not an assumption to build around.

**3. Literature RAG.** PubMed via E-utilities, for evidence not yet in CIViC — recent case reports, off-label use. This is the one component RareCure's own reference implementation never actually built (its Module 5 is an 8-line stub returning `insufficient_evidence=True` unconditionally, despite a fully-specified config; 0 of 47,719 output options in its own committed run carried a PMID, against a README that advertised "PMID-grounded" retrieval). See Deliverable 4.

**4. Drug-gene matching.** Reused from Tier 2 where possible: the Tier 2 modules `dgidb.py`, `opentargets.py`, and `chembl.py` in this repo's `src/secondlook/` already implement a verified, working DGIdb → Open Targets → ChEMBL fallback chain with per-source exception handling. Reuse those clients directly rather than rebuilding — they were specifically hardened against the exact bug (`DgidbError`/`OpenTargetsError`/`ChemblError`, each caught at its call site) that broke RareCure's own equivalent module (`rarecure-build-reference.md` §2: RareCure's DGIdb v2 endpoint silently returned zero results because a `200 text/html` response passed `raise_for_status()` and then failed `.json()` inside a swallowed exception). Separately, the RareCure repo itself should be checked for any directly-reusable module per its data-availability/methods section — but do not reuse its `drug_match.py` DGIdb client specifically, given the confirmed defect above.

**Shared/cross-cutting scope also owned here** (per `architecture.md`'s module table and `checkpoint.md`'s "out of scope for the Tier 2 repo, tracked elsewhere" list):

- **Activation logic** — the rule deciding whether Tier 2 also runs (§6).
- **Orchestration backend (FastAPI)** — `POST /api/query`, `GET /api/results/{query_id}` (Deliverable 8).
- **LLM synthesis layer** — with citation enforcement as executable code, not a prompt instruction (Deliverable 9).
- **Frontend** — per `ui-flow.md` (Deliverable 10).

---

## 2. Why a graph, and what it is ultimately for

The immediate job (CIViC variant → drug → evidence lookup) could technically be done with a relational table. The graph earns its place because of where this system is going, and that destination should be understood before the schema is designed — retrofitting a graph schema after the fact is far more expensive than designing for it now.

**The reference use case** (`docs/founder-mode-use-case-and-kg-solutions.md`, derived from `reference_repos/founder-mode-cancer/`): GitLab co-founder Sid Sijbrandij's osteosarcoma recurrence, where standard-of-care and standard genomic workup were both exhausted. What actually worked, mechanically:

1. Standard genomic panels found nothing new after recurrence — the mutation profile was exhausted as an actionable lead.
2. The breakthrough came from a **different data modality**: single-cell RNA-seq revealed FAP overexpression. FAP is not a standard osteosarcoma biomarker; no mutation-only pipeline would ever surface it.
3. That target was matched to an **existing but unconventional modality** — FAP-targeted radioligand therapy, retargeted from prostate cancer (177Lu-PSMA) by swapping the ligand while keeping the isotope chemistry.
4. Access came through a **path most patients and physicians don't know exists** — FDA Individual Patient Expanded Access (Form 3926), approved within 48 hours. The real bottleneck was hospital IRBs, not the FDA.
5. Multiple **mechanistically-complementary therapies were layered in parallel**, not tried sequentially — checkpoint inhibitor + neoantigen vaccine + oncolytic virus + NK cells + radioligand, each hitting a different part of the immune-evasion problem. Tumor-infiltrating T cells went from 19% to 89%.
6. **AI's role was narrow and specific**: literature synthesis, cross-referencing a specific molecular profile against case reports, and structuring a transparent weighted decision — never the decision-maker. A human expert team executed and validated everything.

Every one of steps 2–5 is a **multi-hop traversal across heterogeneous entity types**: signature → target → modality-class → platform → access-pathway → precedent, where the useful path isn't known in advance (that's exactly what makes FAP a *non-obvious* target). A flat evidence list cannot express "find me options that are complementary to what this patient is already on, and not redundant with it." A pure vector-RAG retrieves *documents that sound relevant*; it does not traverse *structured relationships between entities*. Sid used plain ChatGPT and hit exactly this limitation — good at summarizing what it was pointed at, unable to do the structured reasoning that drove the outcome.

**FalkorDB specifically** because it is a property graph with a **native vector index in the same engine as Cypher**. One query can combine structured traversal with semantic search over evidence/literature text in a single retrieval pass, rather than running a graph DB and a vector DB side by side and merging results in application code. It also means the demo-scale build and any later production build share one storage engine, rather than requiring a migration off NetworkX's in-memory-only model.

**What this means for the V1 build:** build the CIViC backbone and literature RAG on a schema that already has room for the extension node types in §3.3. Do not build those extensions in V1 — but do not design a schema that forecloses them either.

---

## 3. Knowledge graph design

### 3.1 Core node types (V1 — build these)

Every node carries provenance properties: `source` (which database/API), `source_version` or `source_release` where available, and `retrieved_at` (ISO timestamp). This is non-negotiable and is a direct response to a named RareCure failure (`rarecure-build-reference.md` §4): its outputs recorded no retrieval timestamp or data version anywhere, so its published metrics could not be reconstructed from its own committed artifacts.

| Node | Key properties |
|---|---|
| `Gene` | `symbol`, `ensembl_id`, `uniprot_accession`, `hgnc_id` |
| `Variant` | `hgvs_p`, `hgvs_c`, `protein_position`, `ref_aa`, `alt_aa`, `variant_type` (missense/indel/fusion/splice), `civic_variant_id` |
| `Disease` | `name`, `doid` (Disease Ontology ID), `is_in_scope_cancer_type` |
| `Drug` | `name`, `chembl_id`, `approval_status`, `smiles`, `india_availability` (CDSCO status if wired — `ui-flow.md` names this as a results tiebreaker) |
| `EvidenceItem` | `civic_id`, `evidence_level` (A–E), `evidence_type`, `clinical_significance`, `direction`, `summary`, `citation_url` |
| `Trial` | `registry_id` (NCT or CTRI), `registry`, `status`, `phase`, `locations[]`, `country_codes[]`, `eligibility_url` |
| `Publication` | `pmid`, `title`, `journal`, `year`, `pub_type[]`, `mesh_terms[]`, `abstract`, `abstract_embedding` (vector) |

### 3.2 Core edge types (V1)

```
(Gene)-[:HAS_VARIANT]->(Variant)
(Variant)-[:OBSERVED_IN]->(Disease)
(Variant)-[:PREDICTS_RESPONSE_TO {direction: sensitive|resistant}]->(Drug)
(EvidenceItem)-[:SUPPORTS]->(Variant|Drug|Disease)
(EvidenceItem)-[:CITES]->(Publication)
(Trial)-[:RECRUITS_FOR]->(Disease)
(Trial)-[:TARGETS_BIOMARKER]->(Variant|Gene)
(Trial)-[:INVESTIGATES]->(Drug)
(Publication)-[:MENTIONS]->(Gene|Variant|Drug)
```

The `direction` property on `PREDICTS_RESPONSE_TO` is load-bearing — "this drug is associated with this variant" is useless without knowing whether the association is sensitivity or resistance. RareCure's dedup collapsed drugs by name only and lost gene attribution entirely, producing output like `"Pazopanib" | gene: "TP53"` (`rarecure-build-reference.md` §5). Key relationships on the full tuple, never on drug name alone.

### 3.3 Extension node types (design for, do NOT build in V1)

These are the graph-native capabilities from `founder-mode-use-case-and-kg-solutions.md` §3. Leave room in the schema; build only if V1 lands early or a later phase funds it.

| Node | Purpose | Why it's not V1 |
|---|---|---|
| `ModalityClass` | Radioligand-conjugate, CAR-T, oncolytic-virus, checkpoint-inhibitor, neoantigen-vaccine, small-molecule | Requires curating modality taxonomy; no existing source indexes it |
| `PlatformTechnology` | "177Lu/225Ac radioligand platform", "mRNA neoantigen platform" — with `RETARGETABLE_TO` edges encoding the swap-the-payload-keep-the-platform pattern | No data source in `data-sources.md` captures cross-target platform reusability; needs literature curation |
| `RegulatoryPathway` | FDA 3926 expanded-access precedent by modality-class and indication | Genuinely novel asset — nobody publishes this as structured data. Highest differentiation, highest curation cost |
| `ExpressionSignature` | e.g. "FAP overexpression in stromal-mimicking tumor cells" → `IMPLICATES` → `Target` | Requires an expression-data ingestion path that does not exist anywhere in either tier today |
| `Target` | Generalizes beyond `Gene` — a protein target reachable by expression, not only by mutation | Only useful once `ExpressionSignature` or `ModalityClass` exist |

Extension edges, for reference:
```
(Drug)-[:IN_MODALITY_CLASS]->(ModalityClass)
(ModalityClass)-[:COMPLEMENTARY_TO|MECHANISTICALLY_REDUNDANT_WITH|INTERACTION_RISK]->(ModalityClass)
(PlatformTechnology)-[:RETARGETABLE_TO]->(Target)
(ModalityClass)-[:ACCESSIBLE_VIA]->(RegulatoryPathway)
(ExpressionSignature)-[:IMPLICATES]->(Target)
```

### 3.4 Where Tier 2 results enter the graph

Tier 2's computed signals are graph citizens, not a separate output channel. When Tier 2 runs, its result becomes:

```
(Variant)-[:HAS_COMPUTATIONAL_SIGNAL]->(StructuralSignal)-[:PREDICTS_BINDING_CHANGE {delta, label, method}]->(Drug)
```

with `StructuralSignal` carrying `alphamissense_score`, `alphamissense_class`, `structure_source`, `structure_id`, `plddt_at_residue`, `reliability_flag`, `method` (mCSM-lig | docking), and the fixed §10 disclaimer reference. This keeps the "computational, not documented" distinction explicit in the graph itself — a `StructuralSignal` node can never be mistaken for an `EvidenceItem`, which matters because `ui-flow.md` requires these render as visually distinct sections and `architecture.md` forbids conflating them. See `tier2-implementation-spec.md` §4 for the Tier 2 side of this contract.

---

## 4. Retrieval design — hybrid Cypher + semantic

Three retrieval modes, in order of preference. **The mode used must be recorded on every returned item** so the UI and the synthesis layer can distinguish an exact structured match from a fuzzy semantic one.

**Mode 1 — Exact structured (Cypher).** Gene + exact variant → drugs + evidence. This is the strong-hit path.

```cypher
MATCH (g:Gene {symbol: $gene})-[:HAS_VARIANT]->(v:Variant {hgvs_p: $hgvs_p})
MATCH (v)-[r:PREDICTS_RESPONSE_TO]->(d:Drug)
MATCH (e:EvidenceItem)-[:SUPPORTS]->(v)
WHERE e.evidence_level IN ['A','B','C','D','E']
RETURN d.name, r.direction, e.evidence_level, e.citation_url, e.summary
ORDER BY e.evidence_level
```

**Mode 2 — Structured relaxation (Cypher).** No exact variant match → same gene, different variant; or same disease + drug class. This is the weak-hit path that triggers Tier 2. Relaxation must be **explicit and labeled** — a gene-level hit returned as if it were a variant-level hit is precisely the substitution that makes clinical output misleading.

**Mode 3 — Semantic (FalkorDB vector index).** No structured match at any relaxation level → vector search over `Publication.abstract_embedding` and `EvidenceItem.summary` embeddings. Hybrid with BM25 (§Deliverable 4) because gene symbols, `p.V600E`, NCT IDs, and drug names are lexical tokens that dense-only retrieval underperforms on.

The value of one engine: modes 2 and 3 compose in a single query — traverse structurally as far as possible, then vector-search from the nodes you landed on, without a round trip through application code.

**Hard rule carried from `api-contracts.md`:** an item with no `citation.url` must never be constructed or returned, regardless of retrieval mode.

---

## 5. Technology stack

| Layer | Technology | Notes |
|---|---|---|
| Backend | Python, FastAPI | Async support for fanning out to CIViC/CTgov/CTRI/PubMed concurrently with partial-failure tolerance |
| Frontend | Streamlit (fastest) or React | Streamlit if frontend time is limited; React if `ui-flow.md`'s tier-distinction styling (green vs. amber) needs to look polished |
| Knowledge graph | **FalkorDB** (Redis-based, Cypher-compatible, native vector index) | Not NetworkX/Neo4j — see §2. Single engine for both graph traversal and the Literature RAG's semantic search |
| Graph client | `falkordb` Python client | **Unverified in this repo** — confirm the client API and vector-index syntax against current docs before building around it (§Deliverable 1) |
| Embeddings / RAG | `sentence-transformers` for embedding generation; vectors stored/queried via FalkorDB's native index | Free, local, no API cost; combine with BM25 for hybrid retrieval |
| LLM | Whatever API access is available | Must support structured/constrained output — required for the citation-enforcement post-check |
| Mutation parsing | Reused from Tier 2: `hgvs` via `src/secondlook/mutation_validation.py` | Do not reimplement — already built, tested, hardened |
| Drug-gene matching | Reused from Tier 2: `dgidb.py`, `opentargets.py`, `chembl.py` | Already verified live; DGIdb and Open Targets are both GraphQL |

**External services:**

| Source | Purpose | Access | Auth | Status |
|---|---|---|---|---|
| CIViC | Variant-disease-drug evidence graph | Open, no license | None | **Unverified** — confirm current API/download format before building the loader; this repo's verification pass covered Tier 2 sources only |
| ClinicalTrials.gov API v2 | Global trial search | REST, free | None | Not yet verified in this repo |
| CTRI (ctri.nic.in) | India trial registry | **Unconfirmed** — API vs. scrape | Verify | Confirm first; may need a materially different integration than CTgov |
| PubMed E-utilities | Literature abstracts | REST, free | Optional API key (raises rate limit) | Not yet verified in this repo |
| RareCure modules | Reuse candidate | MIT license | None | Treat as a source of *patterns and documented failure modes* more than working code — two of six modules are stubs, several others have confirmed defects |

**Environment variables:**
```
LLM_API_KEY=
PUBMED_API_KEY=            # optional, raises rate limit
CTGOV_API_BASE=https://clinicaltrials.gov/api/v2
CTRI_ACCESS_MODE=          # confirm actual access method first
UNIPROT_API_BASE=https://rest.uniprot.org
FALKORDB_HOST=             # default localhost for local dev
FALKORDB_PORT=             # default 6379 (Redis protocol)
FALKORDB_GRAPH_NAME=       # e.g. "secondlook_tier1"
```

---

## 6. Activation logic

The rule deciding whether Tier 2 also runs. Thresholds below are **initial values** — log the actual hit-rate distribution during testing and tune before demo. Shipping untested defaults risks Tier 2 firing on every query or never firing at all.

| Condition | Behavior |
|---|---|
| **Strong hit** — ≥1 CIViC evidence item at level A/B for the *exact* variant, OR ≥1 recruiting trial matching the exact biomarker | Tier 1 result shown; Tier 2 does not run automatically |
| **Weak/partial hit** — gene has CIViC entries but not this exact amino-acid change, OR only level C–E evidence, OR literature-only (no curated DB entry), OR no matching trial | Tier 1 shows what it found, labeled with its evidence tier; Tier 2 runs automatically as a supplementary signal |
| **No hit** — nothing anywhere | Explicit "no documented evidence found" (never a blank screen); Tier 2 runs |
| **Manual override** | Doctor can request Tier 2 output even on a strong hit (UI toggle) |

The retrieval mode from §4 maps directly onto this: Mode 1 → strong, Mode 2/3 → weak, nothing → no hit.

---

## 7. Deliverables — what "done" looks like

1. **Source verification pass**, mirroring `tier2-verification.md`. A single manual call against CIViC, ClinicalTrials.gov v2, CTRI, and PubMed E-utilities, confirming actual response *shape* — not just HTTP 200. RareCure's diagnostic checked `status_code == 200` only and reported `[ OK ]` for a source that was returning HTML instead of JSON and contributing zero results. **Also verify FalkorDB itself**: stand up an instance, confirm the Python client API, create a vector index, and run one hybrid query end-to-end before building on it. Write findings into `docs/` the way `tier2-verification.md` did, including any architecture change the verification forces (e.g. if CTRI requires scraping).

2. **Knowledge graph module.** CIViC data loaded into FalkorDB per the §3 schema, for a defined, versioned (not inline-literal) set of cancer types including sarcoma. Cancer-type and gene scoping must be **config data, not Python literals** — RareCure's ontology was inline dict literals despite its README promising YAML, which made its cohort-inclusion logic unauditable. Query functions for retrieval modes 1 and 2 (§4). Output items must match the `api-contracts.md` Tier 1 result item schema exactly:
   ```json
   {
     "type": "documented",
     "source": "CIViC | ClinicalTrials.gov | CTRI | PubMed",
     "evidence_level": "string (e.g. CIViC A-E, or 'literature')",
     "citation": {"id": "string", "url": "string"},
     "summary": "string",
     "drug": "string or null",
     "trial_status": "string or null"
   }
   ```
   Hard rule: no item without `citation.url` is ever constructed or returned.

3. **Clinical trial search module.** CTgov v2 + CTRI by condition + biomarker (falling back to free-text). Adapt RareCure's ontology-expansion pattern (condition → parent ontology terms → drug/gene terms) but fix its two named defects before they recur: **(a) cache on the query string from day one** — RareCure had no caching here and made ~5,700 near-duplicate calls across a 261-patient cohort, dominating its runtime; **(b) resolve geography by country code, not substring match** — RareCure's `"United States" in ctry` gave a sample Indian patient *maximum* proximity to US-only trials, the inverse of correct behavior. For an India-localized product, CTRI results and India-located CTgov sites are first-class, not a fallback.

4. **Literature RAG module**, built to the priority order in `rarecure-build-reference.md` §3:
   - **Hybrid retrieval** (BM25 + dense/`sentence-transformers`), not pure vector — gene symbols, `p.V600E`, drug names, NCT IDs are lexical. Store dense embeddings in FalkorDB's native vector index alongside the graph, so one query can traverse structurally *and* search semantically in the same pass.
   - **Chunk by abstract, not fixed token count** — one PubMed abstract = one `Publication` node, with `(pmid, title, journal, year, pub_type, mesh_terms)` as properties, so every chunk is citable by construction and no chunk spans two unrelated abstracts.
   - **Enforce citation at the schema level** — reject any evidence summary without ≥1 PMID verified to exist in the retrieved set. Strongest available hallucination guardrail for this domain, and a direct instance of the system-wide "no evidence, no claim" rule.
   - **Filter by publication type/recency** — RCTs and systematic reviews weighted above case reports.
   - **Rerank** the fused top-k before final top-5.
   - **A ~50-pair ground-truth eval set**, curatable from CIViC citations already being fetched, wired into an actual eval loop. Without this, none of the above is tunable — and RareCure never built one, which is why its RAG config parameters were all guesses that never met real text.

5. **Drug-gene matching integration.** Import and reuse Tier 2's `dgidb.py` → `opentargets.py` → `chembl.py` chain directly; don't fork. Extend those modules if a Tier-1-specific query shape is needed.

6. **Activation logic implementation** per §6, with actual hit-rate distribution logged during testing before thresholds are locked. Manual override exposed to orchestration and surfaced in the frontend.

7. **Provenance and retrieval-mode recording.** Every graph node carries `source`/`retrieved_at`; every returned item records which retrieval mode produced it. A gene-level relaxation must never be presented as a variant-level match.

8. **Orchestration backend.** `POST /api/query` and `GET /api/results/{query_id}` per `api-contracts.md`, wired to: the reused Tier 2 mutation validator (fail-fast, exact §8 error messages, plus an **explicit documented decision** on whether gene/cancer-type-only trial and literature search still proceeds when the specific mutation fails validation — `api-contracts.md` calls for this decision to be made explicitly rather than left implicit), Tier 1's modules, the activation rule, Tier 2 routing (independently re-enforcing the missense-only gate rather than trusting Tier 1's routing), and final response assembly (`tier1_results[]`, `tier2_results[]`, `failures[]`, `synthesis`, `synthesis_citations[]`). Any module failure returns a **structured failure object** — never `null`, never a silently-empty array. This is the guard against RareCure's most serious documented failure: 260 of 261 patients silently fell back to a generic hard-coded gene panel instead of their own data, recorded only as a warning string that downstream code had to string-match to detect. Nothing here substitutes a default for missing patient-specific input without that substitution being a first-class, visible output field.

9. **LLM synthesis layer.** Structured-output call producing `synthesis` + `synthesis_citations`, with post-generation enforcement **implemented as code**: every sentence-level claim maps to a citation ID, every citation ID maps to a real item in `tier1_results`/`tier2_results`. Apply prompt sanitization (strip/cap free-text — doctor's clinical question, prior-treatment notes — before it enters any prompt, per the `sanitize_for_prompt` pattern) and JSON-retry-on-malformed-output. **This layer never computes or influences a score or ranking** — its job stops at organizing already-computed evidence into text (`rarecure-build-reference.md` §4: RareCure's LLM produced scalar scoring weights that collapsed to 2 distinct vectors across 261 patients while being reported as personalized).

10. **Frontend**, per `ui-flow.md`:
    - Screen 1 (intake): cancer type, gene, mutation (shorthand or HGVS, with format hint).
    - Screen 2 (clinical form): repeatable prior-treatments field, free-text clinical question, manual Tier 2 override toggle.
    - Screen 3 (processing): stage indicator, with an explicit **"cached result" indicator** per section whenever a demo-mode cached result is shown — never presented as live when it isn't.
    - Screen 4 (results): Section A (Tier 1, green, solid border, explicit "no documented evidence found" if empty — never a blank section); Section B (Tier 2, amber/outlined, only if it ran, persistent **"COMPUTATIONAL PLAUSIBILITY SIGNAL — NOT CLINICAL EVIDENCE"** badge, fixed disclaimer once at section top); Section C (failures, rendered whenever a failure object exists); synthesis summary at top referencing only items shown below it.
    - Copy rules everywhere: use "predicted" / "computational signal" / "plausibility" / "for physician review". Never "diagnosis" / "prescribed" / "recommended treatment" / "proven" / "will work".

11. **Demo case selection.** Per `validation-plan.md`, lock the Tier 1 demo case (strong documented CIViC/trial evidence, demonstrating the fast sourced-retrieval path) and coordinate with Tier 2 on caching both demo cases' full outputs so neither runs live on stage for the first time.

12. **Verification pass before calling this done:** every module follows the specific-exception-per-source, caught-at-call-site pattern (no bare `except`/`Exception` — repo-wide convention now, not Tier 2's alone); `pytest -m integration` against live CIViC/CTgov/CTRI/PubMed confirming the source-verification findings still hold; the citation-enforcement check **demonstrated to actually reject** a deliberately-uncited synthesis output in a test, not assumed to work from the prompt; a full walkthrough of one strong-hit and one no-hit case confirming activation logic and failure-object rendering end-to-end.

---

## 8. Open decisions to resolve explicitly

Make these calls and document them; don't let them resolve by default.

- **Invalid-mutation partial search** (Deliverable 8): does gene/cancer-type-only trial and literature search proceed when the specific mutation fails reference-residue validation? `api-contracts.md` flags this as needing an explicit decision.
- **Extension scope** (§3.3): does any post-V1 phase target the *expression-signature → non-obvious-target* path (the actual FAP-discovery mechanism, but requires expression data most patients won't have), or the *modality-agnostic options + expanded-access precedent + combination reasoning* path (buildable on what Tier 1/2 already produce, no new data-ingestion problem)? The second is the more natural next increment; the first is the more dramatic capability. `founder-mode-use-case-and-kg-solutions.md` §4 has the full framing.
- **Balanced scorecard** (if the extensions are built): `ai-role.md` documents Sid's own decision process — pillars, weights, scores, transparent. If SecondLook ever surfaces a multi-criteria ranking, the weighting must be a **disclosed, user-adjustable parameter**, never a hidden LLM output. State this as a requirement before it gets built, not as an implementation detail discovered later.
