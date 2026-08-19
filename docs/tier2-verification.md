# Tier 2 Data Source Verification — Results

Live-tested every source in `data-sources.md`'s Tier 2 table on 2026-08-20. This is the "confirm every API/DB actually works" step from `tech-stack-setup.md` (step 2), done specifically for Tier 2. Findings below change some of the original architecture assumptions — read before starting implementation.

## Summary

| Source | Status | Notes |
|---|---|---|
| UniProt REST | ✅ Pass | Exactly as documented |
| hgvs (PyPI) / Mutalyzer | ✅ Pass | Both live |
| **AlphaMissense** | ✅ Pass — **architecture change** | Use Ensembl VEP REST, not dbNSFP file/hosted lookup — see below |
| RCSB PDB | ✅ Pass | |
| AlphaFold DB | ✅ Pass | pLDDT returned directly as `globalMetricValue` |
| ESM Atlas fold API | ✅ Pass (today) | Docs' "unreliable, don't demo-path" caution stands — historical outages, not scriptable-safety |
| DGIdb | ✅ Pass — **note** | It's GraphQL, not plain REST — endpoint is `POST https://dgidb.org/api/graphql` |
| Open Targets | ✅ Pass | GraphQL, `https://api.platform.opentargets.org/api/v4/graphql` |
| ChEMBL | ✅ Pass | |
| PubChem PUG REST | ✅ Pass | |
| **mCSM-lig** | ⚠️ Confirmed risk | Plain HTML form, single mutation/ligand per submission, no API, no batch mode — see below |
| AutoDock Vina / smina | ✅ Pass (install-time) | `pip install vina` (1.2.7) + `meeko` (0.7.1) for ligand/receptor prep both resolve |
| BindingDB REST | ⚠️ Flaky | `getLigandsByPDBs` fast (1.2s); `getLigandsByUniprot` hung >25s with no response. Matches its "reference/validation only" categorization — don't build a live dependency on it |

## Detail: AlphaMissense — architecture change

Original spec assumed a dbNSFP v4.7c download or a separate hosted lookup, with license to be confirmed separately. Verified instead:

- **Ensembl VEP REST** returns AlphaMissense scores inline on any missense variant lookup, no separate call needed:
  ```
  GET https://rest.ensembl.org/vep/human/hgvs/{transcript}:{hgvs}?content-type=application/json;AlphaMissense=1
  ```
  Response includes `transcript_consequences[].alphamissense.{am_pathogenicity, am_class}` per transcript. Confirmed live: TP53 R175H → `am_pathogenicity: 0.9857, am_class: "likely_pathogenic"` (matches known pathogenic classification).
- License: AlphaMissense predictions are CC BY 4.0 (confirmed current as of this check — was CC BY-NC-SA, relicensed March 2024).
- **Implication:** Step 2 of the Tier 2 pipeline can be built as one VEP call instead of a multi-gigabyte file download or a separate lookup service. Simpler, faster, no storage cost. Update `tier2-structural-prediction.md` §Step 2 to reflect this if you want the spec doc to match.

## Detail: mCSM-lig — confirmed blocker, needs a decision

Inspected the live form at `https://biosig.lab.uq.edu.au/mcsm_lig/prediction`:

- It's a synchronous multi-part HTML form: PDB file upload *or* 4-letter PDB code, one mutation (e.g. `D30N`) + chain, one ligand per submission.
- No JSON API, no documented batch endpoint, no rate-limit docs (matches `data-sources.md`'s existing caution).
- Form posts to `https://biosig.lab.uq.edu.au/mcsm_lig/...` — technically inspectable, but scripting an undocumented endpoint on someone else's interactive-only academic server isn't something to do without asking them; not something I did or recommend automating.

**Decision: automate it.** Confirmed live (2026-08-20) that this is safe to do responsibly:

- Submitted the "Run example" form (D30N, chain A, PDB 2Z4O) and got a result back **synchronously** — a few seconds, one request/response, no job queue, no email/polling. That makes scripted automation straightforward rather than the async-job problem it could have been.
- The result page is plain server-rendered HTML with clearly labeled fields, reliably parseable:
  ```
  Predicted Affinity Change: -2.056 log(affinity fold change) - Destabilizing
  Wild-type: D
  Position: 30
  Mutant-type: N
  Chain: A
  Ligand ID: 065
  Distance to ligand: 2.814 Å
  DUET stability change: -0.087 Kcal/mol
  ```
- The submission form's real fields (read from the live DOM, not guessed): `wild` (file upload), `pdb_code` (text, alternative to file), `mutation` (text, e.g. `D30N`), `chain` (text), `lig_id` (text, ligand 3-letter code), `affin_wt` (text, optional — purpose unconfirmed, leave blank unless Cursor's testing shows it's required), submit button `name="run" value="single"`. The static `<form action=...>` attribute reads `run_example` with `method="get"`, but that's the default state before JS rewrites it on submit — the actual POST target must be discovered by driving the real page (see Prompt E), not by replaying a guessed URL.
- No documented rate limit still stands, so automation must stay conservative and serial: one submission at a time, a real delay between calls (not concurrent/bulk), and a generous timeout with a clean fallback rather than a hang. This is a shared academic server, not a rate-limited public API — scripted politeness here isn't optional.

Given the synchronous, cleanly-parseable response, the most robust automation approach is a **headless browser (Playwright)** filling the actual HTML form and reading the rendered result page — not replaying a hand-built HTTP POST — since JS owns the real submit behavior and a raw POST could silently break on any front-end change. This also means the automated path degrades the same way a human user's would if the site changes, rather than failing in some confusing lower-level way.

AutoDock Vina/smina (confirmed installable) stays as the fallback for anything mCSM-lig itself can't handle (timeout, form rejects the input, no ligand-bound structure available) — see Prompt E.

## Cursor prompts — ready to build now

These four are unblocked by the verification above and can be built independently, each testable in isolation before chaining (per `tech-stack-setup.md` step 6).

### Prompt A — Mutation parsing & validation module
```
Build the mutation validation module for Tier 2 of SecondLook (see docs/tier2-structural-prediction.md §Step 1 and §8).

Input: gene symbol/UniProt accession + mutation notation (HGVS protein format like "p.Arg175His" or shorthand like "R175H").
Steps:
1. Regex pre-parser to normalize shorthand ("R175H") into HGVS protein notation before passing to the `hgvs` package.
2. Fetch canonical sequence via UniProt REST: GET https://rest.uniprot.org/uniprotkb/{accession}.fasta — parse the FASTA, record which accession/isoform was used.
3. Hard validation gate: confirm the reference amino acid at the stated position matches the canonical sequence. On mismatch, return exactly this error (do not proceed further): "Reference residue mismatch: notation claims [X] at position [N] but canonical sequence has [Y]. This usually indicates a transcript/isoform numbering difference. Cannot proceed safely."
4. On match, apply the substitution to get the mutant sequence string.
5. Scope gate: reject (with the exact message in docs/tier2-structural-prediction.md §8) anything that isn't a single missense substitution — insertions, deletions, frameshifts, fusions, splice variants.

Write this as a standalone, independently testable module/function with unit tests covering: valid missense match, reference mismatch, and each out-of-scope mutation type. Use TP53 R175H (UniProt P04637) as the manual test case — verified working against live UniProt REST.
```

### Prompt B — AlphaMissense lookup via Ensembl VEP
```
Build the AlphaMissense pathogenicity lookup step for Tier 2 (docs/tier2-structural-prediction.md §Step 2), using Ensembl VEP REST instead of a dbNSFP download — this was verified live and returns the score inline.

Endpoint: GET https://rest.ensembl.org/vep/human/hgvs/{transcript_id}:{hgvs_c_notation}?content-type=application/json;AlphaMissense=1

Requires converting the validated protein-level mutation (from the mutation-validation module) into a transcript + coding-DNA HGVS notation VEP accepts — confirm the right transcript_id to use per gene (canonical/MANE Select transcript, not just any isoform).

Extract am_pathogenicity (float) and am_class (string) from transcript_consequences[].alphamissense in the response — there may be multiple transcript_consequences entries (one per transcript); use the canonical transcript's value and record which transcript was used.

Apply the exact thresholds from docs/tier2-structural-prediction.md §Step 2: likely_benign 0.000–0.333, ambiguous 0.334–0.564, likely_pathogenic 0.565–1.000 (use these instead of trusting am_class blindly, in case Ensembl's own thresholds ever drift from the paper's).

This score is informational only — it does not gate the rest of the pipeline. Store it for output assembly (§Step 7). Add a fallback error path (per §8) for cases where VEP returns no alphamissense field for the variant (e.g. non-missense or scoring unavailable for the position) — return "AlphaMissense functional score is shown as the only available computational signal" framing is NOT needed here since this is a pre-filter step, not the final failure — instead just record "no AlphaMissense score available" as a field and continue the pipeline; missing AlphaMissense should not block progression to Step 3.

Test case: TP53 R175H, transcript ENST00000269305, should return am_pathogenicity ≈ 0.9857, am_class "likely_pathogenic" (verified live).
```

### Prompt C — Structure sourcing module
```
Build the structure sourcing step for Tier 2 (docs/tier2-structural-prediction.md §Step 3), with this priority order, each tier verified live and working:

1. RCSB PDB — search for an experimental structure of this protein, prefer one that's ligand-bound. Use RCSB's search API (https://search.rcsb.org/rcsbsearch/v2/query) to query by UniProt accession, and https://data.rcsb.org/rest/v1/core/entry/{pdb_id} to pull entry metadata. If multiple hits, prefer ones with a bound ligand relevant to the candidate drugs from Step 4 (may require sequencing this step after Step 4, or doing a generic best-available fetch here and refining in Step 5 — decide based on how the pipeline is orchestrated).
2. AlphaFold DB — if no usable PDB structure, fetch the precomputed model: GET https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}. Extract pLDDT via the response's globalMetricValue (overall) — for the residue-level pLDDT at the specific mutated position, you'll need to parse the per-residue B-factor field from the downloadable PDB/CIF file (the summary API only gives global average — confirm this before assuming residue-level pLDDT is available from the JSON alone).
3. ESM Atlas fold API — POST https://api.esmatlas.com/foldSequence/v1/pdb/ with the raw sequence as the request body — verified working today, but per docs/data-sources.md this has documented outages. Implement it but gate it as fallback-only, never on the live demo path, and add a timeout + explicit failure message (docs/tier2-structural-prediction.md §8: "Structural analysis unavailable or low-confidence...") rather than hanging.

Apply the pLDDT < 70 threshold at the mutated residue per §Step 3 — flag structural reliability as low but do not suppress the result.

Do NOT implement folding of the mutant sequence — always use wild-type/reference structure with the mutation annotated at the relevant position, per the explicit "what NOT to build" list.

Test case: TP53 (P04637) — AlphaFold DB returns model AF-P04637-F1, globalMetricValue 75.06 (verified live).
```

### Prompt D — Candidate drug generation + SMILES fetch
```
Before writing any client code: the last two modules you built (AlphaMissense/Ensembl VEP, then structure sourcing) both shipped with the same bug. Each external API client raised a custom error on failure (EnsemblError, RcsbError, AlphaFoldError), but the orchestrating function never caught it — so a live API outage crashed the whole pipeline instead of falling through to the next source or returning the documented "unavailable" result. Both had to be fixed after review, with tests added afterward to catch the regression.

Do not repeat this. For every external call you add in this module (DGIdb, Open Targets, ChEMBL, PubChem):
1. Give each client's failure a dedicated exception type (e.g. DgidbError, OpenTargetsError, ChemblError, PubChemError), same pattern as the existing modules.
2. At the call site in the orchestrating function, catch that specific exception immediately and treat it as "this source returned nothing" — never let it propagate up uncaught. Do not use a bare except or catch Exception; catch the specific error type, the same way source_structure() now does for RcsbError/AlphaFoldError.
3. Write a test for each source's failure path before considering the step done — simulate the client raising its error and assert the function falls through to the next source (DGIdb down → tries Open Targets; both down → tries ChEMBL) or returns the exact zero-candidates message from §8 if all three fail, rather than raising.
4. The PubChem SMILES fetch is a per-candidate call in a loop — one candidate's PubChem failure must not abort the whole candidate list; catch it per-candidate and drop that candidate (or mark its SMILES unavailable) while the rest of the list still comes back.

Now the actual task. Build the candidate drug generation step for Tier 2 (docs/tier2-structural-prediction.md §Step 4). All three sources verified live:

1. DGIdb — query first, it's GraphQL not REST: POST https://dgidb.org/api/graphql with a query like:
   { genes(names: ["GENE_SYMBOL"]) { nodes { name interactions { drug { name } } } } }
2. Open Targets — if DGIdb is empty, query target-disease-drug associations via GraphQL: POST https://api.platform.opentargets.org/api/v4/graphql. Need the Ensembl gene ID (not symbol) — resolve gene symbol → Ensembl ID first (Ensembl REST or a static mapping).
3. ChEMBL — if both empty, fall back to pathway/family-level candidates: GET https://www.ebi.ac.uk/chembl/api/data/target/search.json?q={gene_symbol}, then pull related activities/mechanisms.

For each candidate drug name, fetch SMILES via PubChem PUG REST:
GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/property/CanonicalSMILES/TXT
Enforce client-side throttling: ≤5 req/sec, ≤400/min (verified endpoint working, but PubChem will rate-limit/block if you don't throttle — implement a token-bucket or simple sleep-based limiter, don't rely on catching 429s reactively).

Cap the shortlist at 10–20 candidates. Rank: exact-protein target > same family > same pathway; approved drugs before investigational within a tier (you'll need a data field indicating approval status — DGIdb interactions have an "interaction_types"/approval info in the full response; check the actual DGIdb GraphQL schema via introspection for the exact field name before assuming).

Zero-candidates case: return the exact message from §8 ("No candidate drugs targeting this protein, its family, or its pathway were found in open databases. Tier 2 cannot generate a signal for this target.") and do not proceed to Step 5.

Test case: TP53 — DGIdb returns real hits (PIRARUBICIN, OXALIPLATIN, CARBOPLATIN, etc., verified live).
```

### Prompt E — Binding affinity via automated mCSM-lig submission
```
Same rule as last time: every external call gets its own exception type, caught at the specific call site, with a test for the failure path — do not let anything here propagate uncaught. This module has more ways to fail than the previous ones (browser automation, file uploads, a slower external computation), so be thorough.

Build the binding-affinity step for Tier 2 (docs/tier2-structural-prediction.md §Step 5), automating mCSM-lig submission via headless browser (Playwright) rather than a hand-built HTTP POST — the form's actual submit behavior is owned by client-side JS (confirmed live: the static form action reads "run_example"/GET even though real submission needs the "single" button and a POST-shaped multipart body), so driving the real page is more robust than replaying a guessed request, and it degrades the same way a human user's session would if the site changes anything.

Target: https://biosig.lab.uq.edu.au/mcsm_lig/prediction

Form fields (confirmed live from the DOM):
- `pdb_code` (text) — 4-letter PDB code; use this rather than the `wild` file-upload field, since Step 3's structure sourcing already gives you a PDB ID when the source is "PDB". If the structure came from AlphaFold DB instead (no PDB ID), you'll need the `wild` file upload with the AlphaFold model's PDB text instead — handle both paths.
- `mutation` (text) — shorthand like "D30N" (wild-type single-letter + position + mutant single-letter). Build this from the validated mutation result (Prompt A), not user input directly.
- `chain` (text) — the PDB chain letter the mutation is on. You'll need to resolve which chain in the fetched structure corresponds to the UniProt-numbered position — do not assume chain "A" always; check it against the structure (SIFTS PDB-UniProt residue mapping, or parse chain identifiers from the PDB file directly).
- `lig_id` (text) — 3-letter ligand code (HET code) of the bound ligand to score against. Comes from Step 4's candidate list — you'll need the candidate drug's PDB HET code, not its PubChem SMILES; these are different identifiers. If a candidate has no known HET code (most won't, since HET codes only exist for ligands that have appeared in a deposited PDB structure), you cannot run mCSM-lig for that candidate — fall back to AutoDock Vina/smina for it instead (see below), not to skipping it.
- `affin_wt` (text) — optional, purpose unconfirmed from the live check. Leave blank first; if the result looks wrong or the form rejects blank, investigate what it expects (test manually, don't guess).
- Submit via the button with name="run" value="single".

Parse the result page (server-rendered HTML, plain labeled text, confirmed live format):
  Predicted Affinity Change: {float} log(affinity fold change) - {Destabilizing|Stabilizing}
  Wild-type: {residue}
  Position: {int}
  Mutant-type: {residue}
  Chain: {letter}
  Ligand ID: {3-letter code}
  Distance to ligand: {float} Å
  DUET stability change: {float} Kcal/mol
Extract at minimum the affinity fold-change delta and the stabilizing/destabilizing label — that's the "delta" required by §Step 5 ("mutant score minus wild-type score... report the delta, not an absolute score alone" — confirm whether mCSM-lig's reported number already is that delta, which it appears to be, or whether it needs combining with anything else).

Rate limiting / etiquette (no documented limit exists, so be conservative on your own initiative): one submission at a time, never concurrent; a real delay (a few seconds minimum) between submissions if scoring multiple candidates for the same mutation; a generous timeout (60s+, since this is a real structural computation, not a fast API) with a clean failure rather than a hang.

Failure handling — every one of these needs its own caught exception type and a fallback, not a crash:
- Playwright can't reach the page / page structure changed (selectors don't match) → treat as mCSM-lig unavailable for this candidate, fall back to AutoDock Vina/smina.
- Form submission times out → same fallback.
- No HET code available for a candidate's ligand → same fallback (this will be the common case, not an edge case — most PubChem/DGIdb candidates won't have appeared in a deposited PDB structure).
- Structure has no bound ligand at all (Step 3 returned an apo structure) → mCSM-lig cannot run for any candidate on this case; go straight to Vina/smina for all candidates, or if Vina/smina also isn't feasible for this case, return the §8 framing that only AlphaMissense + structural context is available.

AutoDock Vina/smina fallback: implement receptor/ligand prep (protonation, grid box centered on the binding site — use Step 3's structure and Step 4's SMILES) and a docking run comparing wild-type vs. mutant structure scores. This is real, nontrivial setup work per the existing docs — build it as a clearly separate function/module from the mCSM-lig client, not fused together, so each can be tested independently. If Vina/smina setup for this prompt turns out too large to do alongside mCSM-lig automation in one pass, stop and tell me — we'll split it into its own prompt rather than rushing it.

Test case: run the actual "Run example" flow (D30N, chain A, PDB 2Z4O, ligand 065) end to end through your Playwright client and confirm you parse -2.056 log(affinity fold change), Destabilizing, distance 2.814 Å — these are verified live values, use them as the integration test's expected result.
```

### Prompt F — AutoDock Vina/smina docking fallback

`src/secondlook/vina_dock.py` currently exists only as a stub — `VinaDockClient.score()` always raises `VinaError("not implemented yet")`. That was deliberate (see Prompt E): Vina/smina setup turned out to be its own real task, not something to bolt onto the mCSM-lig automation. Right now, any candidate whose ligand hasn't appeared in a deposited PDB structure — which is most PubChem/DGIdb candidates, per Prompt E's own finding — falls straight to "unavailable." This prompt replaces the stub with a working implementation.

Same rule as every module since Prompt B: every failure point gets its own caught exception, no bare `except`/`Exception`, and a test for each failure path before calling it done.

```
Implement AutoDock Vina/smina docking in src/secondlook/vina_dock.py as the WT-vs-mutant binding-affinity fallback for Tier 2 Step 5, replacing the current VinaError stub. VinaDockClient.score(pdb_text, smiles, mutation, position) must return a BindingScore with status="scored", method="docking", per docs/tier2-structural-prediction.md §Step 5.

Confirmed installable and already in the dev environment's reach: the `vina` PyPI package (1.2.7, AutoDock Vina's Python bindings) and `meeko` (0.7.1, ligand/receptor preparation for Vina). Use those rather than shelling out to a separate `smina` binary unless you hit a concrete limitation `vina`/`meeko` can't handle — if so, stop and tell me before adding a new external binary dependency.

This is real, multi-step structural work — do not approximate or skip a step silently:

1. **Receptor prep**: take the wild-type structure's PDB text (from Step 3), strip waters/heteroatoms except the target ligand, add hydrogens (protonate) at physiological pH, convert to the PDBQT format Vina needs. Use meeko's receptor preparation utilities, or openbabel/pdbfixer if meeko's coverage is insufficient — check meeko's actual API before assuming it does the whole job (research/confirm this, don't guess a method name).

2. **Mutant structure**: per docs/tier2-structural-prediction.md's explicit "what NOT to build" list, you must NOT re-fold the mutant. The correct approach is in-place side-chain mutation of the existing wild-type structure — replace the residue at the validated position with the mutant amino acid's side chain, keeping the backbone fixed, then let a rotamer-placement step (if your prep library supports one) or Vina's own local search handle the new side chain's conformation. Research what's actually available for this in your prep pipeline (PyMOL's mutagenesis wizard via its Python API, or a simpler geometric side-chain swap) before implementing — this is a real methodological choice, not a detail to wing. If you can't find a defensible way to do this without a paid/unavailable tool, stop and tell me rather than shipping something that silently produces a wrong mutant structure.

3. **Grid box**: center the docking search box on the co-crystallized ligand's coordinates in the wild-type structure (the ligand you're stripping out in step 1 tells you where it was) — this only works when Step 3's structure is ligand-bound, which is already a precondition score_binding() checks before reaching Vina (see binding.py's ligand_bound check). Pad the box with enough margin to let the candidate's SMILES-derived 3D conformer dock (a reasonable default, e.g. ligand bounding box + 5 Å per side, but confirm this is a sane default for Vina rather than assuming).

4. **Ligand prep**: generate a 3D conformer from the candidate's SMILES (RDKit, or whatever meeko expects as input) and convert to PDBQT.

5. **Dock twice**: run Vina against the prepared receptor once with the wild-type side chain, once with the mutant side chain, same grid box, same ligand. Report the delta (mutant docking score minus wild-type docking score) per §Step 5 — not either score alone.

6. **Failure handling** — each of these needs a caught, specific exception and a fallback to BindingScore(status="unavailable", ...) with an accurate message (reuse the BINDING_UNAVAILABLE_MESSAGE pattern from binding.py — this is a binding-scoring failure, not a structure-quality failure, so do not reintroduce the pLDDT-message bug that was just fixed in Step 5's review):
   - Receptor prep fails (bad protonation, missing atoms) → VinaError, caught in score_binding(), returns unavailable.
   - No defensible mutant side-chain placement method available → VinaError, do not fabricate a structure.
   - SMILES → 3D conformer generation fails (some candidate SMILES will be un-embeddable) → VinaError.
   - Vina itself errors or times out on either the WT or mutant run → VinaError. Set a real timeout; a docking run has a bounded reasonable runtime, don't let it hang indefinitely.

7. Write unit tests with fake/small inputs for each failure path, plus one real end-to-end docking test using an actual small ligand-bound PDB structure (reuse the TP53/9C5S case already verified live in Step 3, with a small real ligand from Step 4's candidate list) as an integration test — mark it `@pytest.mark.integration` like the rest, and note in the test that Vina docking scores are stochastic-ish (seed it if the library allows, otherwise assert direction/plausibility rather than an exact float).

This is very likely the largest single unit of work in the Tier 2 pipeline so far — structure prep is explicitly called out as "the hard part" in tech-stack-setup.md. If any single step here (especially step 2, mutant side-chain placement) turns out too large to finish well in one pass, stop and tell me which part, rather than shipping a version that silently produces plausible-looking but methodologically unsound docking scores — that's a worse outcome than "docking unavailable, AlphaMissense + structural context only."
```

## Not yet ready for a Cursor prompt

- **CTRI (Tier 1)** — not verified in this pass (out of scope for "Tier 2" per this session), still flagged unconfirmed in `data-sources.md`.
