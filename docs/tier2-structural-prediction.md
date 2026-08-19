# Tier 2 — Structural Drug-Response Prediction

Build-ready spec. Full research/rationale behind these decisions lives in the Tier 2 research report already produced (`compass_artifact` doc) — this file states the decisions as implementable steps.

## Purpose
Generate a labeled, non-clinical "computational plausibility signal" for mutations with no documented evidence — used only when Tier 1 finds nothing or a weak hit (`tier1-retrieval.md` §Activation logic).

## Pipeline (in order)

### Step 0 — Scope gate
- **Accept:** single amino-acid substitutions (missense) only.
- **Reject with explicit message** (see §8) for: insertions, deletions, frameshifts, gene fusions, splice variants. Do not attempt a degraded analysis on these in V1.

### Step 1 — Mutation parsing & validation
- Parse input notation: HGVS protein format (`p.Arg175His`) or shorthand (`R175H`) — add a regex pre-parser to normalize shorthand into HGVS before using the `hgvs` Python package / Mutalyzer for structured parsing.
- Resolve canonical protein sequence via **UniProt REST** (`https://rest.uniprot.org/uniprotkb/{accession}.fasta`). Always record and display which UniProt accession/isoform was used.
- **Hard validation gate:** confirm the reference amino acid at the stated position matches the canonical sequence. **On mismatch: STOP, return the exact error in §8. Do not proceed.**
- Apply the mutation to the sequence (string substitution at validated position) to get the mutant sequence.

### Step 2 — Pathogenicity pre-filter (AlphaMissense)
- Look up the precomputed AlphaMissense score for this variant (do not run the model — use precomputed score files/API, e.g. via dbNSFP v4.7 "c" branch or hosted lookup).
- Classification thresholds (from the AlphaMissense paper, use exactly): likely_benign 0.000–0.333, ambiguous 0.334–0.564, likely_pathogenic 0.565–1.000.
- This score is stored and shown, but does **not** gate the rest of the pipeline — it's informational context, not a hard stop (a benign-classified mutation can still be worth a drug-binding check).
- Verify current license terms on whatever specific AlphaMissense score file/route you use before shipping (was CC BY-NC-SA, relicensed CC BY 4.0 as of March 2024 for precomputed scores — confirm on your actual source).

### Step 3 — Structure sourcing
- Priority order: (1) experimental structure from **RCSB PDB** if one exists for this protein, ideally drug-bound; (2) precomputed **AlphaFold DB** model (alphafold.ebi.ac.uk); (3) only as last resort, ColabFold/AF2 run on Colab (never on the live demo path — see `tech-stack-setup.md`).
- Record the pLDDT confidence at the mutated residue's position. **Threshold: pLDDT < 70 at that residue → flag structural reliability as low in the output (see §7), do not suppress the result, but label it clearly.**
- Do NOT attempt to fold the mutant separately expecting a different structure — this is a known invalid method (structure is near-invariant for single mutations). Use the wild-type/reference structure with the mutation noted at the relevant position for downstream binding-site distance calculations.

### Step 4 — Candidate drug generation
- Query, in order: (1) DGIdb for drugs directly targeting this gene/protein; (2) Open Targets for target–disease–drug associations if step 1 is empty; (3) ChEMBL for pathway/family-level candidates if both are empty.
- Fetch each candidate's SMILES string via **PubChem PUG REST** (`.../compound/name/{drug}/property/CanonicalSMILES/TXT`) — respect rate limits: ≤5 requests/sec, ≤400/min.
- Cap the shortlist at 10–20 candidates. Rank: exact-protein target > same family > same pathway; approved drugs before investigational within a tier.
- **Zero candidates found (even at pathway level):** return explicit message in §8, do not proceed to Step 5 for this query.

### Step 5 — Binding affinity prediction
- **Primary method: mCSM-lig** (web server: biosig.lab.uq.edu.au/mcsm_lig or structure.bioc.cam.ac.uk/mcsm_lig). Requires a protein–ligand complex structure — pair the Step 3 structure with each Step 4 candidate.
- If no complex structure is obtainable for a given drug: fall back to **AutoDock Vina/smina** docking (wild-type vs. mutant structure, compare docking scores) as the stretch/complementary method.
- **Do NOT use DeepPurpose or other pure sequence-based DTI models for this step** — they are not sensitive to single-residue changes (see research report §5 for the full justification). If already implemented, treat its output as unreliable for mutation-effect scoring and do not surface it as a mutation-specific signal.
- Report the **delta**: mutant score minus wild-type score for the same drug, not an absolute score alone.

### Step 6 — Threshold labeling
- Convert the delta into one of three labels: "likely reduced binding," "likely retained/increased binding," "uncertain / insufficient signal."
- **There is no established clinical threshold for this** — set your own numeric cutoff empirically during validation (`validation-plan.md`) and **label it explicitly as an internal heuristic, not a validated clinical cutoff**, both in code comments and in the UI disclaimer.

### Step 7 — Output assembly
Each Tier 2 result item includes (see `api-contracts.md` for exact schema): mutation (as validated), UniProt accession/isoform used, AlphaMissense score + class, structure source (PDB ID or AlphaFold DB ID) + pLDDT at mutated residue, per-drug: name, SMILES source, method used (mCSM-lig or docking), delta score, label, confidence flag, distance of mutated residue to binding site if computable, and the fixed disclaimer text (§10 below).

## §8 — Failure handling (exact messages)

| Condition | Message returned |
|---|---|
| Reference residue mismatch | "Reference residue mismatch: notation claims [X] at position [N] but canonical sequence has [Y]. This usually indicates a transcript/isoform numbering difference. Cannot proceed safely." |
| Out-of-scope mutation type | "Tier 2 supports single amino-acid substitutions (missense) only. This variant type requires evidence-based review; no computational structural signal is available." |
| Structure fetch failed / very low confidence | "Structural analysis unavailable or low-confidence for this region (pLDDT = X, below reliability threshold). No structure-based binding signal is reported. AlphaMissense functional score is shown as the only available computational signal." |
| Zero candidate drugs | "No candidate drugs targeting this protein, its family, or its pathway were found in open databases. Tier 2 cannot generate a signal for this target." |
| Timeout / external API error | "Tier 2 analysis could not complete (timeout/error). No partial or inferred result is shown to avoid misleading output. Please retry or consult evidence-based resources." |

**Rule: every failure state above must be rendered to the user explicitly. Never return an empty array or null silently — a doctor cannot distinguish "checked, found nothing concerning" from "silently failed," and that ambiguity is dangerous in a clinical-adjacent tool.**

## §10 — Fixed disclaimer text (attach to every Tier 2 result, no exceptions)

> "This is a computational plausibility signal generated by structural prediction tools for research and hypothesis-generation only. It is not clinical evidence, a diagnosis, or a treatment recommendation, and it has not been clinically validated for this mutation. The underlying models predict physical binding/stability changes with limited accuracy and were not designed to predict patient response. These results must be interpreted by a qualified clinician alongside established evidence before any clinical decision. Do not use this output to start, stop, or change any treatment."

## Data retention
No identifiable patient data stored — session operates on de-identified gene+mutation+context only (DPDP Act 2023 alignment).

## What NOT to build for V1
- De novo folding of the mutant structure as a comparison method (invalid — structures are near-invariant for point mutations).
- DeepPurpose or any pure sequence-based DTI model as the mutation-effect scorer.
- Rosetta ddG (multi-day setup time, not feasible in hackathon window).
- Any live, uncached mCSM-lig or folding call on the actual demo stage — pre-compute and cache demo cases (`validation-plan.md`, `tech-stack-setup.md`).
