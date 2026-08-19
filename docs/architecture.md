# Architecture

## System overview

```
                          ┌─────────────────────────┐
                          │   Doctor-facing frontend │
                          │  (intake + form + results)│
                          └────────────┬─────────────┘
                                       │ POST /api/query
                                       ▼
                          ┌─────────────────────────┐
                          │      Backend (FastAPI)   │
                          │   Orchestration layer    │
                          └────────────┬─────────────┘
                                       │
                       ┌───────────────┼────────────────┐
                       ▼                                 ▼
          ┌─────────────────────┐            ┌───────────────────────┐
          │  Entity extraction   │            │   Mutation validator   │
          │  (LLM: free text →   │            │  (reference-residue    │
          │   structured fields) │            │   check, HGVS parse)   │
          └──────────┬───────────┘            └───────────┬────────────┘
                      │                                     │
                      ▼                                     ▼
          ┌───────────────────────────────────────────────────────────┐
          │                    TIER 1 — Evidence retrieval              │
          │  Knowledge graph (CIViC) + trial search (CTgov/CTRI)        │
          │  + literature RAG (PubMed)                                  │
          └──────────────────────┬────────────────────────────────────┘
                                  │
                    hits found?  ─┼─  no hits / weak hits
                     yes │        │
                         ▼        ▼
              ┌──────────────┐  ┌───────────────────────────────────────┐
              │ Tier 1 result │  │       TIER 2 — Structural prediction    │
              │   (evidence-  │  │  AlphaMissense → structure lookup →    │
              │    grounded)  │  │  candidate drugs → mCSM-lig ΔΔ affinity │
              └──────┬────────┘  └───────────────────┬─────────────────────┘
                     │                                │
                     └───────────────┬────────────────┘
                                     ▼
                       ┌───────────────────────────┐
                       │   LLM synthesis layer       │
                       │ (organizes retrieved/computed│
                       │  evidence only — no freeform │
                       │  medical claims permitted)   │
                       └─────────────┬─────────────┘
                                     ▼
                       ┌───────────────────────────┐
                       │  Doctor-facing results UI   │
                       │  Tier 1 (green) / Tier 2     │
                       │  (amber) sections, sourced   │
                       └───────────────────────────┘
```

## Module boundaries

| Module | Responsibility | Does NOT do |
|---|---|---|
| **Frontend** | Intake form, results rendering, tier-distinction styling | Any medical reasoning or scoring |
| **Orchestration backend** | Routes query to Tier 1, checks hit quality, routes to Tier 2 on fallback condition, assembles final response | Direct DB calls (delegates to tier modules) |
| **Entity extraction** | Turns doctor's free-text question + structured fields into normalized query terms (gene, mutation, drug class, condition) | Does not decide treatment relevance |
| **Mutation validator** | Parses HGVS/shorthand, resolves canonical sequence, checks reference-residue match | Does not proceed silently on mismatch — hard stop (see `tier2-structural-prediction.md` §2) |
| **Tier 1 module** | Graph traversal (CIViC) + trial search (CTgov + CTRI) + PubMed RAG | Does not compute predictions — retrieval only |
| **Tier 2 module** | AlphaMissense lookup, structure sourcing, candidate drug generation, mCSM-lig scoring | Does not run if mutation type is out of scope (non-missense) — see input scope rules |
| **LLM synthesis layer** | Formats/explains retrieved+computed evidence into doctor-readable text | Never introduces a claim absent from the evidence payload it received — enforced by prompt AND by the API contract requiring cited `evidence[]` per claim |

## Data flow contract (see `api-contracts.md` for exact schemas)

1. Frontend sends structured case object to backend.
2. Backend runs mutation validator first — **fails fast** on reference-residue mismatch or unsupported mutation type, returns explicit error, does not proceed to any tier.
3. Backend queries Tier 1. Applies "hit quality" rule (`tier1-retrieval.md` §Activation) to decide whether Tier 2 also runs.
4. If Tier 2 runs, backend enforces the mutation-type gate (missense only in V1) independently — Tier 2 refuses non-missense input even if somehow reached.
5. Backend assembles combined result: Tier 1 items (if any) + Tier 2 items (if any) + explicit "no results" state if both are empty — never returns an empty/blank response silently.
6. LLM synthesis step runs last, over the assembled evidence only, and every output claim must resolve to an item in the evidence list or the response is rejected before returning to frontend.

## Failure/timeout policy (system-wide)

Any module failure (API down, timeout, low-confidence result) returns a structured error/low-confidence object — never `null` or an empty array with no explanation. See `tier2-structural-prediction.md` §8 and `api-contracts.md` for exact error shapes.

## Deployment notes

- Patient/case data: no persistent storage of identifiable fields in V1 (see `data-sources.md` and DPDP note in `tier2-structural-prediction.md` §10). Store only de-identified gene+mutation+context needed for the session.
- Demo-day reliability: Tier 2's external calls (mCSM-lig web server, structure DB lookups) should have cached fallbacks for pre-validated demo cases (`validation-plan.md`), never a first-time live call on stage.
