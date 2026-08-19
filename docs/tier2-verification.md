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

**Decision needed before Step 5 gets built** — pick one:
1. **Demo-only path (matches validation-plan.md anyway):** manually submit the 2 locked demo cases through the form once, cache the results locally. No live dependency at all. Lowest risk, fits the "never call mCSM-lig live on stage" rule already in the docs.
2. **Semi-live path:** for any case beyond the 2 cached demos, skip binding-affinity scoring entirely and fall back to the AlphaMissense-only + binding-pocket-distance output (this fallback is already specified in `validation-plan.md`'s "if overall pass rate is below threshold" clause — you'd just be invoking it by default for non-demo cases rather than conditionally).
3. **Local docking path:** use AutoDock Vina/smina (confirmed installable) for anything beyond the cached demo cases, treating mCSM-lig as demo-only and Vina as the "live" fallback method. More setup work (structure prep is real effort per the existing docs) but gives you a scriptable live path.

Recommend **option 1 + 3**: cache mCSM-lig for the 2 demo cases (best accuracy, ρ up to 0.67), use Vina/smina for anything else so the live query path still returns *something* scored rather than immediately punting to AlphaMissense-only.

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
Build the candidate drug generation step for Tier 2 (docs/tier2-structural-prediction.md §Step 4). All three sources verified live:

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

## Not yet ready for a Cursor prompt

- **Step 5 (binding affinity)** — blocked on the mCSM-lig decision above. Once you pick an option, I'll write the corresponding prompt (manual-cache tooling for options 1, or a Vina/smina prep+dock pipeline for option 3).
- **CTRI (Tier 1)** — not verified in this pass (out of scope for "Tier 2" per this session), still flagged unconfirmed in `data-sources.md`.
