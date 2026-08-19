# Tier 1 — Evidence Retrieval

## Purpose
Surface documented, citable treatment options for a patient's mutation/case: curated variant-drug evidence, matching clinical trials, and relevant literature. This is the safer, more defensible core of the product — build and stabilize this before Tier 2.

## Inputs
Structured case object from the doctor form (see `ui-flow.md`): cancer type, gene, mutation (validated by the shared mutation validator — see `tier2-structural-prediction.md` §2, reused here), prior treatments tried + why stopped, doctor's question (free text).

## Components

### 1. Knowledge graph — CIViC backbone
- Source: CIViC (Clinical Interpretation of Variants in Cancer), open, no license required.
- Structure: variant → disease → drug → evidence item, with an evidence level (A–E) per CIViC's own scale.
- Implementation for hackathon scale: in-memory graph (NetworkX) or Neo4j free tier, pre-loaded for a defined set of well-documented cancer types (include sarcoma, to tie back to RareCure's own validation cohort).
- Query: given gene + mutation, traverse to connected drugs and their evidence items.

### 2. Clinical trial search
- Sources: ClinicalTrials.gov API v2 (global) + CTRI — Clinical Trials Registry India (ctri.nic.in) for India-relevant trials.
- Query by condition + biomarker/mutation where the registry supports it; otherwise by condition + free-text drug/gene terms.
- Every trial result must carry: trial ID, recruiting status, location(s), eligibility summary link.

### 3. Literature RAG
- Source: PubMed via E-utilities API (free).
- Use case: supplementary evidence not yet in CIViC — recent case reports, off-label use reports.
- Implementation: vector search (sentence-transformers embeddings) over retrieved abstracts for the query terms; return top-k with PMID + citation.

### 4. Drug-gene matching (reused from RareCure where runnable)
- If the RareCure repo modules are usable as-is: reuse its drug-gene matching module rather than rebuilding.
- If not runnable in your environment: fall back to querying DGIdb/Open Targets/ChEMBL directly (see `data-sources.md`).

## Activation logic (drives whether Tier 2 also runs — see `architecture.md`)

- **Strong hit:** ≥1 CIViC evidence item at level A or B for this exact variant, OR ≥1 recruiting trial matching this exact biomarker. → Tier 1 result shown, Tier 2 does not run automatically.
- **Weak/partial hit:** gene has CIViC entries but not for this exact amino-acid change, OR only level C–E evidence, OR only literature hits (no curated DB entry), OR no matching trial. → Tier 1 shows what it found, labeled with its evidence tier, AND Tier 2 runs automatically as a supplementary signal.
- **No hit:** nothing found anywhere. → Tier 1 returns explicit "no documented evidence found" (never a blank screen), Tier 2 runs.
- **Manual override:** doctor can request Tier 2 output even on a strong hit (button in UI).

These thresholds are initial values — tune based on real query logs once the CIViC subset and trial APIs are wired up; log actual hit-rate distribution during testing before the hackathon to sanity-check the thresholds aren't firing Tier 2 on every query or never firing it at all.

## Output shape (feeds `api-contracts.md`)
Each Tier 1 result item: `{type: "documented", source: "CIViC" | "ClinicalTrials.gov" | "CTRI" | "PubMed", evidence_level: string, citation: {id, url}, summary: string}`. No item without a `citation` may be returned.

## Failure handling
- Any single source (CIViC graph, CTgov API, CTRI, PubMed) failing does not block the others — return partial results with a note on which source failed, never fail the whole query because one source timed out.
- Zero results across all sources is a valid, explicitly-labeled state, not an error.

## What NOT to build for V1
- Full NLP extraction pipeline to build the graph from raw literature — use CIViC's already-structured data instead.
- Production-scale graph DB — in-memory is fine for a demo-scale, pre-loaded cancer-type subset.
