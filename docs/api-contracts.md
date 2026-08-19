# API Contracts

Schemas are illustrative (pseudo-JSON-Schema) — adapt field names to your actual implementation, but preserve the required fields and the "no evidence, no claim" rule below throughout.

## POST /api/query (frontend → backend)

Request body:
```json
{
  "cancer_type": "string",
  "gene": "string",
  "mutation": "string (e.g. 'R175H' or 'p.Arg175His')",
  "prior_treatments": [
    {"drug": "string", "outcome": "string (why it stopped/failed)"}
  ],
  "clinical_question": "string (free text)"
}
```

## Mutation validator output (internal, backend step 1)

```json
{
  "status": "valid | reference_mismatch | unsupported_type",
  "gene": "string",
  "hgvs_normalized": "string",
  "uniprot_accession": "string",
  "isoform_note": "string",
  "position": "integer",
  "reference_residue_expected": "string",
  "reference_residue_actual": "string",
  "mutation_type": "missense | insertion | deletion | frameshift | fusion | splice",
  "error_message": "string or null"
}
```
If `status != "valid"`, backend returns the corresponding message from `tier2-structural-prediction.md` §8 immediately and does not call Tier 1 or Tier 2 for prediction work that depends on a valid mutation object. (Tier 1's trial/literature search by gene/cancer-type alone may still proceed even if the specific mutation fails validation — decide this explicitly and document it in your implementation notes.)

## Tier 1 result item

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
**Rule: an item with no `citation.url` must not be constructed or returned.**

## Tier 2 result item

```json
{
  "type": "computational_signal",
  "mutation_validated": { "...mutation validator output above..." },
  "alphamissense": {"score": "float 0-1", "class": "likely_benign | ambiguous | likely_pathogenic"},
  "structure": {"source": "PDB | AlphaFoldDB", "id": "string", "plddt_at_residue": "float or null", "reliability_flag": "high | low"},
  "drug": "string",
  "smiles_source": "string",
  "method": "mCSM-lig | docking",
  "delta_score": "float",
  "label": "likely_reduced_binding | likely_retained_or_increased_binding | uncertain",
  "binding_site_distance_angstrom": "float or null",
  "disclaimer": "fixed text from tier2-structural-prediction.md §10"
}
```
**Rule: `label` may only be set from a delta_score using the documented, versioned threshold constant in code — never a hardcoded per-case decision. `disclaimer` is mandatory and identical on every item, sourced from a single shared constant, not re-typed per call site.**

## Failure/empty-state object (used by both tiers)

```json
{
  "type": "failure",
  "tier": "1 | 2",
  "reason": "string (one of the exact messages in tier2-structural-prediction.md §8, or an equivalent Tier 1 message)",
  "retryable": "boolean"
}
```
**Rule: this object is always rendered in the UI when produced — never filtered out silently before reaching the frontend.**

## GET /api/results/{query_id} (backend → frontend, final assembled response)

```json
{
  "query_id": "string",
  "tier1_results": ["... Tier 1 result items, or []"],
  "tier2_results": ["... Tier 2 result items, or []"],
  "failures": ["... failure objects for any component that failed, or []"],
  "synthesis": "string (LLM-generated summary text)",
  "synthesis_citations": ["array of citation/result IDs the synthesis text actually references"]
}
```
**Rule: every sentence-level claim in `synthesis` must map to at least one ID in `synthesis_citations`, and every ID in `synthesis_citations` must correspond to an actual item in `tier1_results` or `tier2_results`. Enforce this with a post-generation check (e.g., regex/structured-output parsing) before returning — don't rely on the prompt alone.**
