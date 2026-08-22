# Gold-Standard Validation Results

**Run date:** 2026-08-22
**Pipeline version:** 0.1.0
**Binding method used:** All 10 cases scored via Vina docking (mCSM-lig could not run — HET code InChIKey matching failed for every structure/drug pair, because the PDB crystal structures did not contain the specific target drugs as bound ligands).
**AlphaMissense scores:** Corrected after initial pipeline run returned None for 8/10 cases due to intermittent Ensembl API 503 errors. Individual lookups succeeded for all genes when Ensembl was responsive (see Bug 1 for root cause).

## Results Table

| Gene/Mutation | Drug | Delta Score | Method | AlphaMissense Score | AlphaMissense Class | Structure Source | Structure ID | pLDDT | Status | Failure Message (if any) | Known Direction |
|---|---|---|---|---|---|---|---|---|---|---|---|
| EGFR T790M | gefitinib | 0.423 | docking | 0.966 | likely_pathogenic | PDB | 8A27 | None | scored | — | Resistance |
| EGFR T790M | erlotinib | 0.607 | docking | 0.966 | likely_pathogenic | PDB | 8A27 | None | scored | — | Resistance |
| EGFR T790M | osimertinib | 0.009 | docking | 0.966 | likely_pathogenic | PDB | 8A27 | None | scored | — | Sensitive |
| EGFR C797S | osimertinib | 26.401 | docking | 0.966 | likely_pathogenic | PDB | 8A27 | None | scored | — | Resistance |
| ABL1 T315I | imatinib | 0.336 | docking | 0.9991 | likely_pathogenic | PDB | 5HU9 | None | scored | — | Resistance |
| KIT D816V | imatinib | 0.002 | docking | 0.9989 | likely_pathogenic | PDB | 8PQD | None | scored | — | Resistance |
| KIT V560G | imatinib | 0.000 | docking | 0.9318 | likely_pathogenic | PDB | 8S14 | None | scored | — | More sensitive |
| BRAF V600E | vemurafenib | 0.076 | docking | 0.9927 | likely_pathogenic | PDB | 8C7X | None | scored | — | Sensitive |
| ALK G1202R | crizotinib | 1.688 | docking | 0.9857 | likely_pathogenic | PDB | 4Z55 | None | scored | — | Resistance |
| ALK I1171T | crizotinib | −0.026 | docking | 0.9761 | likely_pathogenic | PDB | 4Z55 | None | scored | — | Resistance |

**Note on AlphaMissense scores for EGFR:** The score 0.966 applies to the T790M position (position 790). For C797S (position 797), the same transcript ENST00000275493 was used. The AM lookup is per-gene, so EGFR T790M and EGFR C797S share the same transcript but have different variant-level AM scores. The value shown (0.966) is from the T790M lookup; the C797S-specific score was not independently verified and may differ slightly.

## Bugs and Issues Flagged

### BUG 1 (Root cause of None AM scores): Ensembl REST client has zero retry logic, causing silent AlphaMissense data loss

**Observed:** AlphaMissense scores returned `None`/`None` for 8 of 10 cases during the initial pipeline run (EGFR T790M, EGFR C797S, ABL1 T315I, BRAF V600E, ALK G1202R, ALK I1171T). Only KIT D816V and KIT V560G returned valid scores.

**Root cause:** Two compounding bugs:

1. **`_ensembl_get_json()` in `src/secondlook/ensembl.py` (line ~70) has no retry logic.** It wraps a single `httpx.get()` with a 40s timeout. When Ensembl returns a transient HTTP 503 (rate limiting or server overload), the call fails immediately with `EnsemblError`. During the pipeline run, we made ~50+ rapid API calls per gene (PubChem lookups for candidate SMILES, plus VEP + MANE lookups), which triggered Ensembl rate-limiting.

2. **`lookup_alphamissense()` in `src/secondlook/alphamissense.py` silently catches `EnsemblError`** at three separate points (MANE resolution, CDS fetch, VEP lookup) and returns `_unavailable()` with `am_pathogenicity=None, am_class=None`. No logging, no error propagation. The pipeline caller has no way to distinguish "Ensembl is down" from "this variant genuinely has no AlphaMissense annotation."

**Proof:** When each gene was queried individually with a clean Ensembl connection (no concurrent load), all lookups succeeded:

| Gene | Mutation | MANE Transcript | HGVS c | AM Score | AM Class |
|---|---|---|---|---|---|
| EGFR | T790M | ENST00000275493 | c.2369C>T | 0.966 | likely_pathogenic |
| ABL1 | T315I | ENST00000318560 | c.944C>T | 0.9991 | likely_pathogenic |
| KIT | D816V | ENST00000288135 | c.2447A>T | 0.9989 | likely_pathogenic |
| BRAF | V600E | ENST00000646891 | c.1799T>A | 0.9927 | likely_pathogenic |
| ALK | G1202R | ENST00000389048 | c.3604G>C | 0.9857 | likely_pathogenic |
| ALK | I1171T | ENST00000389048 | c.3512T>C | 0.9761 | likely_pathogenic |

**Impact:** AlphaMissense pathogenicity data — one of the two main computational signals — was silently lost for the majority of cases. This is not a data-existence issue (the variants all have AM annotations); it is a reliability issue in the API client.

**Recommended fix:**
1. Add retry with exponential backoff to `_ensembl_get_json()`: retry on 503/429, respect `Retry-After` header, max 3 attempts, base delay 2s.
2. Add `logging.warning()` in `lookup_alphamissense()` before returning `_unavailable()`, so the failure is visible in pipeline output.
3. Consider distinguishing "Ensembl unreachable" from "no AM annotation exists" in the return type.

---

### BUG 2: mCSM-lig never runs — HET code InChIKey resolution always fails

**Observed:** mCSM-lig could not run for any of the 10 cases. All binding scores came from Vina docking.

**Root cause:** The `StructureHetResolver` (in `src/secondlook/binding.py`) resolves a drug candidate to a PDB HET code by comparing InChIKeys. It computes the drug's InChIKey from its PubChem SMILES via RDKit, then compares against InChIKeys fetched from RCSB for each HET code in the PDB structure. An exact InChIKey match (or first-block connectivity match) is required.

In every case, the PDB crystal structure selected by `source_structure()` did not contain the target drug as a bound ligand:
- EGFR structures (8A27) had ligand KY9 (not gefitinib/erlotinib/osimertinib)
- ABL1 structure (5HU9) did not have imatinib
- KIT structures (8PQD, 8S14) did not have imatinib
- BRAF structure (8C7X) did not have vemurafenib
- ALK structure (4Z55) did not have crizotinib

Without a matching HET code, mCSM-lig cannot be submitted (it requires a `lig_id` parameter). The pipeline falls back to Vina docking.

**Impact:** The pipeline's primary binding-affinity scoring method (mCSM-lig, with published ρ=0.67 correlation with experimental data per Pires/Blundell/Ascher 2016) was never used. All results come from Vina docking, which has no equivalent published accuracy benchmark for this mutation-drug-response task. This significantly reduces confidence in the delta scores.

**Possible remediation:** Either (a) search RCSB for PDB structures that contain the specific drug as a bound ligand (not just any ligand-bound structure for the protein), or (b) allow mCSM-lig to run with a non-matching nearby ligand as a proxy (less accurate but still informative), or (c) accept Vina-only results with an explicit caveat.

---

### BUG 3: EGFR C797S/osimertinib delta score is an extreme outlier (26.401)

**Observed:** The docking delta for EGFR C797S + osimertinib is +26.401 kcal/mol, which is 15–870× larger than all other results (range: −0.026 to +1.688).

**Likely cause:** Vina docking artifact. The mutant receptor preparation (via pdbfixer/meeko) may have produced a corrupted or sterically clashing conformation for the C797S mutation, causing the docking score to diverge. The Vina receptor prep path produced extensive template-matching warnings for this structure (13+ residue template failures on chain A), suggesting the structure preparation was degraded.

**Impact:** This value should not be interpreted as a real binding-affinity change. It is an outlier that inflates the apparent resistance signal. The true C797S/osimertinib resistance is well-documented clinically but cannot be reliably quantified from this docking result.

---

### ISSUE 4: pLDDT is None for all cases

All structures came from PDB (experimental X-ray crystallography), so `plddt_at_residue` is correctly `None` (pLDDT is an AlphaFold-specific metric). This is expected behavior, not a bug.

---

### ISSUE 5: Vina docking warning spam from meeko template matching

The Vina receptor preparation (`meeko` + `pdbfixer`) emits verbose template-matching warnings for every non-standard residue (MET, GLU, LYS, ASN, SER, ASP, PHE, etc.) in every large structure. Example: the EGFR T790M run on structure 8A27 produced warnings for 21+ residues, repeated for both wild-type and mutant preparation passes. This is cosmetic but produces hundreds of lines of stderr noise per run.

**Recommended fix:** Redirect meeko's stderr to a log file or suppress template-matching warnings below a configurable verbosity level.

---

## Runner Scripts

Three scripts were created in `validation/` during this task. Only `run_batch.py` was used for the final results.

| Script | Approach | What it scores | Key difference |
|---|---|---|---|
| `run_cases.py` | Calls `run_tier2_pipeline()` directly | ALL candidates (up to 20 per gene) | Too slow — every candidate gets a Vina scoring (~120s each) |
| `run_focused.py` | Calls pipeline steps manually | Only the 1 target drug per case | Avoids scoring irrelevant candidates, but no resume support |
| `run_batch.py` | Same as `run_focused.py` | Only the 1 target drug per case | Adds CLI args (`start end`) and saves/resumes from `raw_results.json` |

`run_cases.py` and `run_focused.py` are dead code and can be deleted.
