# SecondLook — AI Second-Opinion Copilot for Rare & Treatment-Exhausted Cancers

India-localized extension of RareCure (open-source, MIT-licensed, peer-reviewed pipeline). Gives oncologists a fast, evidence-ranked second opinion for cases where standard-of-care has failed or the mutation is undocumented.

## Doc index

| File | Covers |
|---|---|
| `architecture.md` | Full system architecture, both tiers, data flow, module boundaries |
| `tier1-retrieval.md` | Evidence-retrieval tier: knowledge graph + literature RAG spec |
| `tier2-structural-prediction.md` | Structural drug-response prediction tier spec (build-ready) |
| `data-sources.md` | Every external API/DB used: endpoints, auth, rate limits, licensing |
| `api-contracts.md` | Internal request/response schemas between frontend, backend, and each module |
| `ui-flow.md` | Doctor-facing UX: intake form → results screen, field-by-field |
| `validation-plan.md` | Gold-standard test cases and pass/fail criteria for both tiers |
| `tech-stack-setup.md` | Stack choices, install/setup steps, environment variables |

## Non-negotiable system-wide rule

**The LLM never states a medical fact it did not retrieve or compute in a traceable step.** Every output — Tier 1 or Tier 2 — carries a source (database ID, paper citation, or named computational method + result). This rule is enforced in the API contract (see `api-contracts.md`) and the UI (see `ui-flow.md`), not just in a prompt — every response object requires a non-empty `evidence` field or it must be rejected before reaching the UI.

## Build order (recommended)

1. `tech-stack-setup.md` — get the environment running
2. `data-sources.md` — confirm every API/DB actually works before building around it
3. `tier1-retrieval.md` — build the documented-evidence path first (it's the safer, more defensible V1 core)
4. `ui-flow.md` — wire up intake form → Tier 1 results, get one end-to-end demo working
5. `tier2-structural-prediction.md` — add the fallback path
6. `validation-plan.md` — run before touching the pitch deck; do not claim reliability you haven't checked

## Scope discipline

V1 = Tier 1 fully working + Tier 2's core path (AlphaMissense + mCSM-lig) on missense mutations only, demoed on pre-validated cases. Everything else in the specs marked "stretch" stays out of V1 unless V1 is done early.
