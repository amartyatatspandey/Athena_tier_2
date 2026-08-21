# Tier 2 — Lead Prompts for Cursor

These are the lead-owned deliverables from `tier2-implementation-spec.md` §5, written as Cursor prompts in the same style as the ones in `docs/tier2-verification.md` (Prompts A–F) that built Steps 1–5b. One prompt per deliverable, run in order. After each one lands, it gets reviewed in this session — code read in full, tests rerun independently, edge cases probed by hand — before you get a green light and a commit message, same as every prior round.

Junior-owned work (`tier2-junior-tasks.md`) is explicitly **not** in here and should not start until Prompt G below is reviewed and merged — Junior Task 1 needs a working pipeline function to run cases through.

---

## Prompt G — Step 7 pipeline orchestration

```
Build the Tier 2 pipeline orchestration function — this wires together six already-built, already-tested modules into one entry point. Read them first; do not reimplement any of their logic, only call them.

Existing modules and their exact signatures (all in src/secondlook/):

  validate_mutation(identifier: str, mutation: str, *, sequence_provider=None) -> MutationValidationResult
    — mutation_validation.py. .status is "valid" | "reference_mismatch" | "unsupported_type".
    On non-"valid", .error_message already holds the exact §8 message (OUT_OF_SCOPE_MESSAGE for
    unsupported_type, reference_mismatch_message(...) output for reference_mismatch) — reuse it,
    do not re-derive or re-type these strings.

  lookup_alphamissense(validation: MutationValidationResult, *, transcript_resolver=None, vep_client=None) -> AlphaMissenseResult
    — alphamissense.py. .status is "scored" | "unavailable". Already handles its own failures
    (EnsemblError caught internally) — never raises.

  source_structure(validation: MutationValidationResult, *, pdb_client=None, alphafold_client=None,
                    esm_client=None, allow_de_novo_folding=False, preferred_ligands=()) -> StructureResult
    — structure.py. .status is "found" | "unavailable". Leave allow_de_novo_folding at its default
    (False) — ESM Atlas must never be on this pipeline's live path, per data-sources.md's documented-
    outage caution. Already handles its own failures (RcsbError/AlphaFoldError caught internally).

  generate_candidates(validation: MutationValidationResult, *, dgidb_client=None, opentargets_client=None,
                       chembl_client=None, pubchem_client=None, max_candidates=20) -> CandidateResult
    — candidates.py. .status is "found" | "none". On "none", .error_message already holds
    ZERO_CANDIDATES_MESSAGE (candidates.py) — reuse it.

  score_binding(validation, structure, candidate: DrugCandidate, *, mcsm_client=None, vina_client=None,
                het_resolver=None, min_delay_seconds=3.0, sleeper=None) -> BindingScore
    — binding.py. .status is "scored" | "unavailable". Already handles structure.status != "found"
    internally (returns unavailable with structure_unavailable_message) and every mCSM/Vina failure
    path (returns unavailable with BINDING_UNAVAILABLE_MESSAGE) — you do not need to gate calling this
    on structure or candidate quality, it degrades safely on its own. Note min_delay_seconds=3.0 is a
    real sleep between mCSM-lig submissions (etiquette on a shared academic server, no documented rate
    limit) — when looping this over multiple candidates for one mutation, that delay is real wall-clock
    time per candidate; don't parallelize calls to this function in a way that defeats the delay.

Build one function, e.g. run_tier2_pipeline(identifier: str, mutation: str, ...) -> Tier2PipelineResult,
that does, in order:

Step 0/1 — call validate_mutation. If .status != "valid", return a single top-level failure object
  (not a list of per-candidate items) carrying the exact .error_message already on the result. Do not
  call any other module.

Step 2 — call lookup_alphamissense(validation). Always call this regardless of what Step 3/4 find —
  it's informational context per tier2-structural-prediction.md §Step 2, never a gate.

Step 3 — call source_structure(validation). If .status == "unavailable", do NOT stop the whole
  pipeline — AlphaMissense is still a valid signal on its own (this is explicitly allowed by the
  §8 "Structure fetch failed" message's own text: "AlphaMissense functional score is shown as the
  only available computational signal"). Continue to Step 4.

Step 4 — call generate_candidates(validation). If .status == "none", return a single top-level
  failure object with ZERO_CANDIDATES_MESSAGE. Do not proceed to Step 5 — there is nothing to score.

Step 5/5b — for each candidate in the Step 4 result, call score_binding(validation, structure,
  candidate). This produces one BindingScore per candidate.

Step 6 — labeling is NOT part of this prompt. The numeric threshold hasn't been decided yet (it's
  chosen empirically from a validation run that needs this pipeline to exist first — see
  tier2-implementation-spec.md §5 item 3). Do not invent a threshold. Instead:
  - Create src/secondlook/labeling.py with a function label_binding_score(delta_score: float | None) -> str
    that is a CLEARLY MARKED PLACEHOLDER: it always returns "uncertain", with a module-level constant
    THRESHOLD_DECIDED = False and a docstring/comment stating this is a placeholder pending the lead's
    validation-run decision, and pointing at tier2-implementation-spec.md §5 item 3 and
    tier2-junior-tasks.md Task 5 as where the real implementation will land.
  - The pipeline function should accept an optional labeling_fn parameter (defaulting to
    label_binding_score from labeling.py) so the real implementation can be swapped in later without
    touching this function's code at all — same dependency-injection pattern every other module in
    this repo already uses (sequence_provider, transcript_resolver, pdb_client, etc.).

Step 7 — assemble output. For each candidate, build one item matching this exact shape
  (api-contracts.md's "Tier 2 result item" schema):

  {
    "type": "computational_signal",
    "mutation_validated": {... the full MutationValidationResult, as a dict ...},
    "alphamissense": {"score": float | None, "class": str | None},
    "structure": {"source": str | None, "id": str | None, "plddt_at_residue": float | None, "reliability_flag": str | None},
    "drug": candidate.name,
    "smiles_source": candidate.smiles_source,
    "method": binding_score.method,   # "mCSM-lig" | "docking" | None
    "delta_score": binding_score.delta_score,
    "label": labeling_fn(binding_score.delta_score),
    "binding_site_distance_angstrom": binding_score.distance_angstrom,
    "disclaimer": TIER2_DISCLAIMER
  }

  Create TIER2_DISCLAIMER as a single module-level constant (this doesn't exist anywhere in the repo
  yet — check before assuming it does) holding the exact §10 text from tier2-structural-prediction.md:
  "This is a computational plausibility signal generated by structural prediction tools for research
  and hypothesis-generation only. It is not clinical evidence, a diagnosis, or a treatment
  recommendation, and it has not been clinically validated for this mutation. The underlying models
  predict physical binding/stability changes with limited accuracy and were not designed to predict
  patient response. These results must be interpreted by a qualified clinician alongside established
  evidence before any clinical decision. Do not use this output to start, stop, or change any
  treatment." — attach it from this one constant to every item; never re-type it per item.

  If binding_score.status == "unavailable" for a given candidate, still include that candidate's item
  in the output (method/delta_score/label as None, binding_score.error_message surfaced somewhere on
  the item — add a field for this, e.g. "binding_note") rather than dropping the candidate silently.
  Every candidate Step 4 found should appear in the output in some form — never silently fewer items
  than candidates found, per the project's "no evidence, no claim, and no silent gaps either" rule.

Graph integration shape (tier2-implementation-spec.md §4): also provide a function
  to_structural_signal(item: dict) -> dict returning the StructuralSignal shape:
  {alphamissense_score, alphamissense_class, structure_source, structure_id, plddt_at_residue,
   reliability_flag, method, binding_site_distance_angstrom, computed_at: <ISO timestamp of when this
   ran>, pipeline_version: <a version string constant>}. This does NOT write to any database — Tier 2
   stays dependency-light per the open-decision recommendation in tier2-implementation-spec.md §7
   (return objects, don't perform the graph write here). Just shape the data correctly so Tier 1's
   orchestration layer can consume it later.

Failure handling: same rule as every module before this one — every external call in this pipeline is
  already wrapped by the module it lives in (you're not adding new external calls here, only
  orchestrating), but the orchestration itself must not introduce a new unhandled-exception path. If
  calling into any of the five modules raises anything unexpected (it shouldn't, per their own
  contracts, but assert this rather than assume it), that's a bug in this new code, not something to
  swallow with a bare except — let it surface in tests rather than catching broadly.

Tests:
  - Unit tests with fake/injected clients for: full happy path (valid mutation, structure found,
    candidates found, binding scored) producing the correct number of correctly-shaped items with the
    disclaimer attached; reference-residue mismatch short-circuits to the single failure object with
    the exact message and calls nothing else; unsupported mutation type short-circuits the same way;
    zero candidates short-circuits to the exact ZERO_CANDIDATES_MESSAGE failure object without calling
    score_binding; structure unavailable does NOT short-circuit (AlphaMissense + candidates still run);
    a candidate whose binding_score.status is "unavailable" still appears in the output with its note.
  - Test that label_binding_score always returns "uncertain" right now (locks in the placeholder
    behavior so nobody is surprised later) and that the pipeline function's labeling_fn parameter can
    be overridden with a fake in a test (proving the injection point actually works).
  - Test to_structural_signal produces the documented shape from a sample item.
  - @pytest.mark.integration test: run the full pipeline live on TP53 R175H (identifier "P04637",
    mutation "R175H") — every individual value here has already been verified live at each step in
    prior prompts (UniProt sequence, AlphaMissense 0.9857/likely_pathogenic, PDB 9C5S structure,
    DGIdb candidates including real drugs, at least one binding score). Assert the assembled output
    has the right shape, non-empty candidates, and the disclaimer text matches exactly.

Test case to sanity-check the failure paths against: gene "NOTAREALGENE123" / mutation "R175H" should
  produce a single failure-shaped result (not a crash) — this exercises the fail-fast Step 1 gate.
```

---

## After Prompt G lands

Bring it back here for review before doing anything else — including before telling the junior to start Task 1. Once it's reviewed and green-lit, the handoff sequence is:

1. Junior Task 1 (run the 9 gold-standard cases through the now-real pipeline).
2. You (lead) pick the Step 6 cutoff from Task 1's data — `tier2-implementation-spec.md` §5 item 3.
3. Junior Task 5 replaces `labeling.py`'s placeholder with the real threshold, using the exact value you hand them.
4. Junior Tasks 2–4 and 6 can run any time after Task 1 (2, and after you lock the demo case, 3) or in parallel with the above (4, 6).
