# Tier 2 Implementation Spec — Structural Drug-Response Prediction

Project: SecondLook (AI second-opinion copilot for rare and treatment-exhausted cancers, India-localized extension of RareCure). This document is the handoff spec for whoever continues Tier 2 work in this repo. It reflects the actual current state of the code (`src/secondlook/`), not just the original design — five of seven pipeline steps are built and independently verified; this spec covers what's left and how to build it consistently with what's already there.

Source docs this spec is derived from, in order of authority: `docs/checkpoint.md` (current build state — read this first), `docs/tier2-structural-prediction.md` (the build-ready pipeline spec), `docs/tier2-verification.md` (live source verification + the six Cursor prompts used to build Steps 1–5b), `docs/api-contracts.md` (the exact output schema), `docs/validation-plan.md` (the gold-standard test set), `docs/data-sources.md`, `docs/tech-stack-setup.md`, `docs/rarecure-build-reference.md` (failure modes to avoid, several directly relevant to the remaining work), `docs/founder-mode-use-case-and-kg-solutions.md` (the capability gap Tier 2's current scope leaves open, and the post-V1 extension paths), `tier1-implementation-spec.md` (the knowledge graph Tier 2's output feeds into).

## 0. Delegation model — who owns what

This document is now **lead-owned**. It covers everything requiring full-system judgment: the Step 7 architecture, the Step 6 threshold decision, and the two open decisions in §7. All of it either produces a decision that's expensive to unwind, or requires holding the whole pipeline's contracts in your head at once.

Everything mechanical-but-real — has a clear "done" state, doesn't require system-wide judgment, and shouldn't need you in the loop mid-task — is delegated to **`tier2-junior-tasks.md`**, a self-contained document the junior can work from without reading this spec end-to-end. The handoff points are marked inline below with → **[JUNIOR]**.

The dependency order matters: Deliverable 2 (Step 7, lead) has to exist before **[JUNIOR] Task 1** can run; **[JUNIOR] Task 1**'s raw data has to exist before the lead can make the Step 6 cutoff decision (Deliverable 1); the cutoff has to be decided before **[JUNIOR] Task 5** (implementing `labeling.py`) can start. `tier2-junior-tasks.md` states these dependencies explicitly so the junior knows what they're blocked on and what to flag back to you rather than guess past.

---

## 1. Scope and task

Tier 2 generates a **computational plausibility signal** for a cancer mutation when Tier 1 (documented evidence retrieval — CIViC, clinical trials, PubMed) finds nothing or only a weak match. It is not a diagnosis or a treatment recommendation. It runs a missense mutation through: pathogenicity scoring (AlphaMissense via Ensembl), protein structure sourcing (PDB/AlphaFold), candidate drug generation (DGIdb/Open Targets/ChEMBL), and predicted binding-affinity change (mCSM-lig or AutoDock Vina), then labels each candidate and attaches a fixed non-clinical disclaimer.

### 1.1 What's already built and verified

Per `checkpoint.md` — do not re-derive scope for these. They are done, unit- and integration-tested, and independently re-verified in code review (every live claim re-checked from scratch rather than taken from the build report).

| Step | Files | What it does |
|---|---|---|
| 1. Mutation parsing & validation | `mutation_validation.py`, `uniprot.py` | Normalizes shorthand (`R175H`) or HGVS (`p.Arg175His`), fetches canonical sequence from UniProt REST, hard-gates on reference-residue mismatch, rejects non-missense variant types |
| 2. AlphaMissense pathogenicity | `alphamissense.py`, `ensembl.py` | Ensembl VEP REST call (not the originally-spec'd dbNSFP download — a verified simplification found during source verification), classified with the AlphaMissense paper's fixed numeric thresholds, not Ensembl's own class label |
| 3. Structure sourcing | `structure.py`, `rcsb.py`, `alphafold.py`, `esm_atlas.py` | RCSB PDB (ligand-bound preferred) → AlphaFold DB (residue-level pLDDT via CA B-factor) → ESM Atlas (off by default, never demo path). Never folds the mutant — uses WT/reference structure with the mutation annotated |
| 4. Candidate drug generation | `candidates.py`, `dgidb.py`, `opentargets.py`, `chembl.py`, `pubchem.py` | DGIdb → Open Targets → ChEMBL fallback chain, PubChem SMILES with a token-bucket rate limiter (≤5 req/s), per-candidate failure isolation |
| 5. Binding affinity — mCSM-lig | `mcsm_lig.py`, `binding.py` | Headless-Playwright automation of the mCSM-lig web form (no API exists), resolves the real PDB chain rather than assuming `"A"` |
| 5b. Binding affinity — Vina fallback | `vina_dock.py` | For the (common) case where a candidate has no PDB HET code so mCSM-lig can't run. PDBFixer + meeko receptor prep, in-place side-chain mutation (never re-folds), grid box centered on the co-crystallized ligand, each run in its own subprocess with a timeout |

### 1.2 The failure-handling contract — mandatory for all new work

All six modules above follow one contract, and it is **non-negotiable for every new module below**: every external call gets its own dedicated exception type, caught at the specific call site (never a bare `except`/`Exception`), falling through to the next source or an explicit "unavailable" result — never crashing the pipeline, never silently returning `None`/`[]`.

This rule is explicit because it was violated three times and caught three times in review:
- **Steps 2 and 3** originally let an uncaught `EnsemblError`/`RcsbError`/`AlphaFoldError` crash the whole lookup instead of falling through to the next source.
- **Step 5** reused the pLDDT-based "structure is unreliable" message for binding-scoring failures, which produced a factually false claim (`"pLDDT = n/a, below reliability threshold"`) about a 1.01 Å experimental crystal structure that was never in question. Fixed by splitting `BINDING_UNAVAILABLE_MESSAGE` from `structure_unavailable_message()`.

That last one matters beyond its own bug: **an inaccurate explanation of *why* something is unavailable is itself a false clinical-adjacent claim**, and the system-wide rule in `README.md` ("the LLM never states a medical fact it did not retrieve or compute in a traceable step") applies to error messages too, not just to results.

### 1.3 What's in scope for this handoff (not yet built)

1. **Step 6 — Threshold labeling.** Convert each candidate's `delta_score` into one of three labels: `likely_reduced_binding`, `likely_retained_or_increased_binding`, `uncertain`. There is no established clinical cutoff — the numeric threshold must be chosen empirically by running the gold-standard cases (item 3), set as a versioned code constant, and documented explicitly as an internal heuristic (in code comments *and* the UI disclaimer), never a hardcoded per-case decision.

2. **Step 7 — Output assembly / pipeline orchestration.** The six modules above are independently callable functions, not a pipeline. Build one function running Steps 1→6 in order, looping Steps 5–6 over every Step-4 candidate, assembling results into the exact `api-contracts.md` schema (§5, Deliverable 2), attaching the fixed §10 disclaimer from a single shared constant — not re-typed at each call site.

3. **Run `validation-plan.md`'s nine gold-standard cases** through the assembled pipeline (EGFR T790M ×2 drugs, EGFR C797S, ABL1 T315I, KIT D816V, KIT V560G, BRAF V600E, ALK G1202R, ALK I1171T), checking resistance-vs-sensitive directionality against **pre-committed** pass criteria — do not adjust criteria after seeing results:
   - Pass threshold: correct directionality on ≥70% of cases.
   - Hard requirement regardless of overall score: BRAF V600E/vemurafenib and EGFR T790M/osimertinib (the two "designed-for" positive controls) must both show retained/increased binding. Either failing means the pipeline is not demo-ready as-is.
   - If below threshold: fall back to the narrower claim ("mutation is in/near the known binding pocket" + AlphaMissense flag only), dropping the binding-affinity delta/label from the UI. A legitimate, documented fallback — not a failure to hide.
   - This run also determines the Step 6 cutoff; the cutoff cannot be picked before it.

4. **Demo-case caching.** Pre-compute and locally cache the full pipeline output (all external calls resolved) for the locked Tier 2 demo case (one gold-standard case, sarcoma-related if possible per `validation-plan.md`'s tie-back to RareCure's cohort). The live demo must never depend on external API latency/uptime — especially mCSM-lig (synchronous but undocumented speed) and Vina (real compute time).

5. **Graph integration contract** — emit Tier 2 results as `StructuralSignal` nodes for the Tier 1-owned FalkorDB graph (§4). New in this revision; previously Tier 2 output had no defined destination beyond the API response.

6. **Document the Vina coverage gap honestly.** Per Prompt E/F's own finding, most PubChem/DGIdb candidates have no PDB HET code, so **Vina — not mCSM-lig — carries most real traffic** despite being the nominal "fallback." mCSM-lig has a published correlation figure (ρ up to 0.67, Pires/Blundell/Ascher 2016); Vina/smina has no equivalent figure for this task. Step 6's labeling/confidence framing must represent this honestly per `validation-plan.md`'s "cite real accuracy figures, not assumed ones" — do not imply uniform confidence across mCSM-lig-scored and Vina-scored candidates.

### 1.4 Explicitly out of scope for this repo

Owned elsewhere per `checkpoint.md`: CTRI access-mode verification, the FastAPI orchestration layer deciding Tier 1 vs. Tier 2 routing, the frontend, and the LLM synthesis layer's citation-enforcement check. All are covered in `tier1-implementation-spec.md`. Coordinate with whoever owns Tier 1 on the Step 7 exposure question (§7).

---

## 2. Scope boundaries, and the capability gap they leave

Tier 2's current scope is deliberately narrow, and the narrowness is correct for V1 — but it should be understood explicitly rather than discovered later.

**Two hard gates:**
- **Missense only.** Step 0 rejects insertions, deletions, frameshifts, fusions, and splice variants with an exact §8 message. Not a limitation to fix quietly — a degraded analysis on an out-of-scope variant type would be worse than an explicit refusal.
- **Drug-binding mechanism only.** Tier 2 answers exactly one question: *does this amino-acid substitution change how well this small molecule binds this protein?* It says nothing about expression, pathway rewiring, immune context, or any non-small-molecule modality.

**The gap this leaves** (from `docs/founder-mode-use-case-and-kg-solutions.md`): the reference case that motivates this whole project — the GitLab founder's osteosarcoma recurrence — had its breakthrough via a route Tier 2 structurally cannot reach. Standard genomic panels were exhausted; the actionable target (FAP) came from **single-cell RNA expression data, with no mutation involved at all**, and was hit by a **radioligand conjugate**, not a small molecule. Tier 2's Step 4 candidate generation (DGIdb/Open Targets/ChEMBL) indexes small molecules almost exclusively; none of those sources would surface a radioligand, oncolytic virus, or cell therapy as a candidate.

So: Tier 2 as built is a correct, verified answer to a real and narrower question. The founder-mode capability requires **additional entry points**, not a repair of the existing pipeline. Those are §6 — post-V1, explicitly not this handoff.

---

## 3. Technology stack

| Layer | Technology | Notes |
|---|---|---|
| Language / runtime | Python ≥3.11 | `pyproject.toml`, package `secondlook`, `src/` layout |
| Package/build | setuptools ≥68 (`[tool.setuptools.packages.find] where = ["src"]`) | |
| HTTP client | `httpx` ≥0.27 | All REST/GraphQL calls |
| Mutation notation | `hgvs` ≥1.5 | Structured HGVS parsing after a shorthand→HGVS regex pre-parse |
| Browser automation | `playwright` ≥1.49 | Drives the live mCSM-lig HTML form — required because the real submit target is JS-owned, not a documented API |
| Docking | `vina` ==1.2.7, `meeko` ==0.7.1 | Confirmed installable (`vina` needed a source build requiring SWIG). Don't shell out to a separate `smina` binary unless `vina`/`meeko` hit a concrete limitation |
| Cheminformatics | `rdkit` | SMILES → 3D conformer for docking |
| Structure prep | `pdbfixer` (+ `openmm`), `gemmi` | Protonation at pH 7.4, in-place side-chain mutation via `applyMutations`, PDB/CIF parsing |
| Numerical | `scipy` | meeko dependency |
| Graph client | `falkordb` (only if Step 7 writes directly — see §4) | **Unverified in this repo.** Tier 1 owns the graph; coordinate before adding this dependency |
| Testing | `pytest` ≥8.0, `pythonpath = ["src"]`, marker `integration` (deselected by default) | Existing convention: unit tests + one `@pytest.mark.integration` test against the real service per module. Keep this for Steps 6–7 |

**External services (all verified live per `tier2-verification.md`):**

| Service | Purpose | Access | Caveat |
|---|---|---|---|
| UniProt REST | Canonical protein sequence | `GET /uniprotkb/{accession}.fasta` | None documented |
| Ensembl VEP REST | AlphaMissense score, inline | `GET /vep/human/hgvs/{transcript}:{hgvs}?AlphaMissense=1` | Requires resolving canonical/MANE Select transcript per gene |
| RCSB PDB | Experimental structures | `search.rcsb.org/rcsbsearch/v2/query`, `data.rcsb.org/rest/v1/core/entry/{id}` | |
| AlphaFold DB | Precomputed models | `GET alphafold.ebi.ac.uk/api/prediction/{accession}` | Global pLDDT only in JSON; residue-level requires parsing CA B-factor from the structure file |
| ESM Atlas | De novo fold, fallback-only | `POST api.esmatlas.com/foldSequence/v1/pdb/` | Historical outages — never on demo path, gate behind a timeout |
| DGIdb | Drug-gene interactions | GraphQL, `POST dgidb.org/api/graphql` | Not REST despite older docs |
| Open Targets | Target-disease-drug | GraphQL, `POST api.platform.opentargets.org/api/v4/graphql` | Needs Ensembl gene ID, not symbol |
| ChEMBL | Pathway/family candidates | REST, `ebi.ac.uk/chembl/api/data/target/search.json` | |
| PubChem PUG REST | Candidate SMILES | `GET .../compound/name/{drug}/property/CanonicalSMILES/TXT` | ≤5 req/s, ≤400/min — token-bucket already in `pubchem.py`, reuse it |
| mCSM-lig | Primary binding method | `biosig.lab.uq.edu.au/mcsm_lig/prediction`, HTML form only | No documented rate limit — stay serial, one submission at a time, real delay between calls, generous timeout. Shared academic server |
| AutoDock Vina | Fallback binding method | Local via `vina`/`meeko` | Structure prep is the hard part; real compute time |

---

## 4. Graph integration contract

Tier 1 owns a FalkorDB knowledge graph (`tier1-implementation-spec.md` §3). Tier 2's results are graph citizens, not a parallel output channel. Step 7 must emit results in a shape the graph loader can consume:

```
(Variant)-[:HAS_COMPUTATIONAL_SIGNAL]->(StructuralSignal)
(StructuralSignal)-[:PREDICTS_BINDING_CHANGE {delta_score, label, method}]->(Drug)
```

`StructuralSignal` properties: `alphamissense_score`, `alphamissense_class`, `structure_source` (PDB | AlphaFoldDB), `structure_id`, `plddt_at_residue`, `reliability_flag`, `method` (mCSM-lig | docking), `binding_site_distance_angstrom`, `computed_at`, `pipeline_version`.

**Two rules on this boundary:**

1. **A `StructuralSignal` is never a `EvidenceItem`.** The graph must make "computed" versus "documented" structurally impossible to confuse, because `ui-flow.md` requires them to render as visually distinct sections and `architecture.md` forbids conflating them. Separate node label, separate edge type, no shared parent.

2. **`computed_at` and `pipeline_version` are mandatory.** A structural signal is a *point-in-time computation against external databases that change*. Without provenance, a cached six-month-old signal is indistinguishable from a fresh one — exactly the failure that made RareCure's published metrics unreproducible from its own artifacts (`rarecure-build-reference.md` §4: no retrieval timestamp or data version recorded anywhere in its outputs).

Whether Step 7 writes to FalkorDB directly or returns objects for the orchestration layer to persist is one of the open decisions in §7 — resolve it with Tier 1's owner before building.

---

## 5. Deliverables — what "done" looks like

Lead-owned work only. Delegated tasks are marked → **[JUNIOR]** with the corresponding task number in `tier2-junior-tasks.md`; they're listed here too so the full Step 6/7 picture stays in one place, but the actual instructions live in the junior doc.

1. **Pipeline orchestration function (Step 7) — build this first, everything else depends on it.**
   - Single entry point: given a validated case (gene, mutation notation), runs Steps 1→6, loops Steps 5–6 per candidate, returns `api-contracts.md`-shaped items:
     ```json
     {
       "type": "computational_signal",
       "mutation_validated": { "...mutation validator output..." },
       "alphamissense": {"score": "float 0-1", "class": "likely_benign | ambiguous | likely_pathogenic"},
       "structure": {"source": "PDB | AlphaFoldDB", "id": "string", "plddt_at_residue": "float or null", "reliability_flag": "high | low"},
       "drug": "string",
       "smiles_source": "string",
       "method": "mCSM-lig | docking",
       "delta_score": "float",
       "label": "likely_reduced_binding | likely_retained_or_increased_binding | uncertain",
       "binding_site_distance_angstrom": "float or null",
       "disclaimer": "fixed text, from a single shared constant"
     }
     ```
   - Every exact §8 failure message (reference-residue mismatch, out-of-scope type, structure unavailable/low-confidence, zero candidates, timeout/API error) returned as a **structured failure object** when triggered — never a silent empty list or `null`. Reuses the `BINDING_UNAVAILABLE_MESSAGE` split already fixed in `binding.py` rather than reintroducing the pLDDT-message conflation bug.
   - Emits graph-ready `StructuralSignal` objects per §4.
   - Integration test running end-to-end on TP53 R175H (values already verified at each individual step) as a smoke test.
   - This is the piece that requires holding every module's contract in your head simultaneously — not a good delegation target even though individual pieces of it (writing one more integration test, say) look small in isolation. Keep it.

2. **→ [JUNIOR] Task 1 — run the nine gold-standard cases.** Once Deliverable 1 exists, hand the pipeline function to the junior. They run all nine `validation-plan.md` cases through it and produce the raw results table (predicted delta/method per case) — no interpretation, no pass/fail call. That call is yours (item 3 below).

3. **Step 6 threshold decision — yours, informed by Task 1's data.**
   - Look at Task 1's raw results against `validation-plan.md`'s pre-committed pass criteria (≥70% correct directionality; BRAF V600E/vemurafenib and EGFR T790M/osimertinib — the two positive controls — must both show retained/increased binding, no exceptions). **Do not adjust the criteria after seeing the results.**
   - Pick the numeric delta-score cutoff that separates the three labels, as a versioned constant. → **[JUNIOR] Task 6** can hand you a literature-grounded starting range (mCSM-lig/Vina delta magnitudes treated as meaningful in comparable published validations) before you commit to a number — use it as input, not as the decision itself.
   - If the ≥70% threshold isn't met: invoke the documented fallback (drop the binding-affinity delta/label from the UI, keep AlphaMissense + binding-pocket-proximity only). This is a legitimate, pre-authorized outcome per `validation-plan.md`, not something to relitigate under demo-timeline pressure.
   - Once decided, hand the exact cutoff value + label boundaries to **[JUNIOR] Task 5**, who implements `labeling.py` to that exact spec.

4. **→ [JUNIOR] Task 2 — Vina coverage tally**, computed from Task 1's raw data.

5. **→ [JUNIOR] Task 3 — demo-case caching.** You lock *which* gold-standard case is the demo case (a judgment call — sarcoma-related if possible, per `validation-plan.md`'s tie-back to RareCure's cohort); the junior executes the caching itself.

6. **→ [JUNIOR] Task 4 — verification-pass checklist** (test suite, bare-except grep, schema field-for-field check). Review their report before merging; don't rerun it yourself unless something in it looks off.

7. **Written recommendation on the open decisions** (§7) — yours. These are architecture calls (graph write responsibility, shared-doc split) that shape how Tier 1 and Tier 2 integrate; a PR description or `docs/` addendum is enough, but the call itself isn't delegable.

8. **Final review and merge** of every junior-produced artifact before it lands — the delegation model in §0 assumes you're the last check on each handoff, not a bystander to it.

---

## 6. Post-V1 extension paths (do not build in this handoff)

From `docs/founder-mode-use-case-and-kg-solutions.md` §3. Recorded here so the §2 gap has a documented answer, and so Step 7's design doesn't accidentally foreclose them.

- **Expression-signature entry point.** A second, parallel entry into Tier 2 that takes an `ExpressionSignature` (e.g. FAP overexpression) instead of a validated missense mutation, and traverses signature → target rather than mutation → target. This is the direct computational analog of the founder-mode breakthrough. Requires an expression-data ingestion path that exists nowhere in either tier today, and expression data most patients in this project's target population won't have — which is exactly why it's not V1.
- **Modality-agnostic candidate generation.** Extend Step 4 beyond small molecules to modality *classes* (radioligand-conjugate, CAR-T, oncolytic-virus, checkpoint-inhibitor, neoantigen-vaccine). Mostly a schema and data-sourcing problem, not a structural-prediction one — note that mCSM-lig and Vina both fundamentally require a small-molecule ligand, so **non-small-molecule modalities cannot be scored by Tier 2's existing methods at all**. They'd surface as candidates with evidence and mechanism but no binding delta, which is honest and still useful.
- **Platform-retargeting and expanded-access precedent.** Graph-native, Tier 1-owned (`tier1-implementation-spec.md` §3.3) — Tier 2's only involvement would be contributing structural feasibility signals where a small-molecule analog exists.

**Design constraint for Step 7 today:** don't hard-code the assumption that every candidate has a `delta_score`. A candidate with a real evidence trail and no computable binding delta is a legitimate future output shape, and the schema already permits `null` for `binding_site_distance_angstrom` — treat `delta_score`/`label` as similarly optional in the internal model even though every V1 candidate will have both.

---

## 7. Open decisions to resolve explicitly

- **Step 7 exposure surface** — importable function only (current assumption), or a small HTTP surface of its own? Affects how the Tier 1-owned orchestration backend integrates. Decide *with* Tier 1's owner, not unilaterally.
- **Graph write responsibility** (§4) — does Step 7 write `StructuralSignal` nodes to FalkorDB directly (adds a `falkordb` dependency and a live DB requirement to this repo's tests), or return objects the orchestration layer persists (keeps Tier 2 dependency-light and testable offline)? **Recommendation: the latter** — it keeps Tier 2's test suite free of a database dependency and matches the existing pattern where every module returns structured objects rather than performing side effects.
- **Shared-doc split** — do `data-sources.md`, `api-contracts.md`, `architecture.md`, `ui-flow.md`, `validation-plan.md` get split so this repo holds only their Tier 2 sections, or stay whole (current state) because Tier 2's code depends on content in them (e.g. the result-item schema)? Flag a recommendation rather than leaving it open indefinitely.
