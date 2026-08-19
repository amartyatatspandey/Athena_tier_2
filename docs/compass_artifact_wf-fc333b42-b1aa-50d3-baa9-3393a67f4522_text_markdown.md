# Tier 2 Design Decision Document: Structural Drug-Response Prediction for Undocumented Cancer Mutations (RareCure Extension)

## TL;DR

- **The originally proposed Tier 2 pipeline (ESMFold structure change + DeepPurpose binding affinity to judge a point mutation) is NOT scientifically defensible and must be re-architected** — sequence-based DTI models like DeepPurpose largely memorize protein identity and are effectively blind to single amino-acid changes, and ESMFold/AlphaFold typically return a near-identical fold for a point mutation. Building on these would produce confident-looking noise.
- **The defensible redesign** keeps the same UX but swaps the engine: AlphaMissense precomputed-score lookup as a functional pre-filter, plus a mutation-specific protein–ligand ΔΔ-affinity predictor — chiefly **mCSM-lig**, which was built for exactly this question and trained on the Platinum database (correlation up to ρ = 0.67 with experiment) — reporting a mutant-vs-wildtype **delta**, never an absolute response prediction.
- **Frame Tier 2 as a labeled "computational plausibility signal / hypothesis-generator," never a diagnosis or prescription**, validate by rediscovering textbook resistance cases (ABL1 T315I, EGFR T790M/C797S, KIT D816V, BRAF V600E, ALK G1202R), pre-cache demo results (never fold live on stage), and pre-empt CDSCO SaMD / DPDP-2023 questions.

---

## Key Findings

1. **DeepPurpose cannot detect the effect of a single amino-acid substitution.** It is a sequence-based DTI library (Huang et al., *Bioinformatics* 2020, 36(22-23):5545–5547) trained on DAVIS/KIBA/BindingDB. Sequence-based DTI/DTA models are documented to rely on protein-identity memorization and degrade sharply in cold-target splits; a single mutation is ~1 residue out of hundreds and is essentially invisible to the encoder. A 2024 attention-DTI paper explicitly states its model "is currently unable to clearly distinguish between mutated and wild-type proteins in the Davis kinase dataset."
2. **ESMFold/AlphaFold usually do not change the predicted structure for a point mutation.** Pak et al. (*PLoS ONE* 2023) and adversarial-mutation studies (bioRxiv 2026) show "nonphysical structural invariance"; ESMFold is somewhat more mutation-sensitive than AlphaFold but neither is reliable for single-residue effects. "Fold the mutant and compare" is an invalid method.
3. **mCSM-lig is the right-fit tool** — purpose-built to quantify how a single missense mutation changes protein–small-molecule affinity (Pires, Blundell & Ascher, *Sci Rep* 2016), trained on the Platinum database, correlation **up to ρ = 0.67** with experiment. Independent 2026 benchmarking (Pan, Portelli, Nguyen & Ascher, *Brief Bioinform* 27(1):bbag035) shows it and PremPLI lead on Platinum-derived data (R≈0.75–0.76) but predictions degrade 5–30% when structures are modeled rather than experimental.
4. **No validated clinical threshold exists** to convert a ΔΔ-affinity into a user label; the |ΔΔG| > ~1 kcal/mol "destabilizing" convention (used in the CAGI5 stability challenge and Rosetta/FoldX literature) is a heuristic, not a clinical boundary.
5. **The whole V1 stack runs on free tools / Colab free tier** — the only heavy step (de novo folding) should be avoided in the live path.

---

## Details

### 1. Activation logic
**DECISION: Activate on zero Tier 1 hits AND on weak/partial hits; run as a labeled fallback (not silently blended); always allow manual override.**
- **Trigger:** fire when the Tier 1 KG (CIViC + ClinicalTrials.gov/CTRI + PubMed RAG) returns zero hits OR only weak hits (no CIViC Level A/B item, no matching interventional trial, RAG similarity below your cutoff). Tier 2 exists precisely to cover the "no documented evidence" gap.
- **Parallel vs. fallback:** strictly a *labeled fallback* in the UI so a low-evidence computational signal is never conflated with high-evidence documented data. You MAY compute it silently in the background on documented cases for telemetry/validation.
- **Partially documented (gene studied, exact AA change not):** treat as **Tier 1 partial + Tier 2 run** — show gene-level context AND the separated computational signal. This is Tier 2's highest-value scenario.
- **Manual override:** yes — a "Run structural analysis anyway" button is defensible because output is explicitly non-clinical.
- **Flag:** the exact "weak hit" thresholds are arbitrary until calibrated on your own retrieval logs — treat as tunable parameters.

### 2. Input scope
**DECISION: Missense point mutations only in V1; accept HGVS protein (`p.Arg175His`) and shorthand (`R175H`); resolve via UniProt canonical isoform; hard-validate the reference residue before proceeding.**
- **Parsing (all free/open-source):** the **`hgvs` Python package** (Hart et al., *Bioinformatics* 2015; 2018 update; hgvs.readthedocs.io) is the reference implementation, installable fully locally; **Mutalyzer** (mutalyzer.nl/api, no auth; `mutalyzer-hgvs-parser` on PyPI) for syntax checking/structured models; **VariantValidator** (Freeman et al., *Human Mutation* 2018) wraps hgvs and auto-corrects. Shorthand like `R175H` is not strictly HGVS-compliant, so add a small regex pre-parser to map shorthand → `p.Arg175His`.
- **Scope justification — missense only:** every downstream mutation-effect tool you can realistically use (mCSM-lig, FoldX, DDGun, AlphaMissense) is built and validated for single substitutions; missense is the dominant actionable class; indels/frameshifts/fusions/splice each need fundamentally different modeling (truncations, novel junctions, altered reading frames).
- **Out-of-scope types:** never silently drop — return *"Tier 2 supports single amino-acid substitutions (missense) only. This variant type requires evidence-based review; no computational structural signal is available."*
- **Sequence source: UniProt REST** (`https://rest.uniprot.org/uniprotkb/{accession}.fasta` for the canonical sequence; JSON for metadata). RefSeq/Ensembl as cross-references only.
- **Pitfalls (all real):** (a) canonical isoform selection — numbering is isoform-dependent; TP53 R175H is defined on canonical isoform UniProt P04637 — always surface the accession/isoform used; (b) multiple transcripts give different residue numbers (biggest silent-error risk); (c) pin the UniProt entry version for reproducibility; (d) HGVS positions are 1-based while Python strings are 0-based — off-by-one errors are catastrophic and silent.
- **Apply + validate (non-negotiable):** fetch canonical sequence → confirm the wild-type residue at position N matches the notation's asserted residue (the "R" in R175H must be Arg at 175). **On mismatch, STOP:** *"Reference residue mismatch: notation claims Arg at position 175 but canonical sequence has [X]. This usually indicates a transcript/isoform numbering difference. Cannot proceed safely."*

### 3. Candidate drug generation
**DECISION: Exact-protein targets first, then family/pathway fallback, via DGIdb + Open Targets + ChEMBL (all genuinely free); SMILES from PubChem PUG REST; cap ~10–20 ranked by directness; filter to India-available drugs AFTER scoring.**

| Source | Gives | Free for students? |
|---|---|---|
| **DGIdb** | Drug–gene interactions (aggregated) | Yes, open |
| **Open Targets** | Target–disease–drug (GraphQL); drug data from ChEMBL | Yes, open, free API |
| **ChEMBL** | Bioactivities, mechanisms, SMILES ("Open Data") | Yes, fully open |
| **PubChem PUG REST** | SMILES, CIDs, synonyms | Yes, fully open |
| **BindingDB** | Experimental binding affinities | Yes, free (CC BY) |
| **PharmGKB** | Pharmacogenomics | Free for research; terms apply |
| **DrugBank** | Rich drug annotation | **NOT freely usable** — academic license required; avoid depending on it |
| **OncoKB** | Curated oncology variant-drug | Website free; **API needs registration + academic license**, treatment data licensed — use for building reference lists, not as a live dependency |

- **SMILES:** `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drugname}/property/CanonicalSMILES/TXT`. **Rate limits (verbatim from PubChem policy / Kim et al., *J Cheminform* 2019): no more than 5 requests/second, no more than 400/minute, no longer than 300 s running time per minute.** The `pubchempy` wrapper simplifies this; handle both `CanonicalSMILES`/`ConnectivitySMILES` field names.
- **Cap:** 10–20/query (downstream ΔΔG/docking cost scales linearly). **Rank** by directness: exact-protein interaction > same family > same pathway; approved over investigational within a tier.
- **India availability (CDSCO):** filter *after* scoring so the strongest computational candidate is visible even if import is needed; show CDSCO status as metadata.
- **Zero candidates even at pathway level:** return the explicit "no candidate drugs found" message (Section 8), never an empty card.

### 4. Structure prediction
**DECISION: Prefer experimental PDB and AlphaFold DB precomputed lookups over de novo folding; do NOT rely on ESMFold/AlphaFold to detect a point mutation's structural effect; use AlphaMissense precomputed-score lookup as a pathogenicity pre-filter.**

| Option | Access | Reality |
|---|---|---|
| **Experimental PDB (RCSB)** | Free | Best when it exists; EGFR/ABL1/BRAF/KIT have drug-bound structures |
| **AlphaFold DB lookup** | alphafold.ebi.ac.uk, free | Instant precomputed human models; preferred over self-folding |
| **ColabFold (AF2)** | Free Colab | ~5–10 min/seq; free tier disconnects ~90 min idle / ~12 h total; ~2,500-residue monomer cap on free GPU |
| **ESMFold (HF/local)** | Free | Fast single-sequence; no longer SOTA |
| **ESM Atlas fold API** (`api.esmatlas.com/foldSequence/v1/pdb/`) | Was free | **STATUS UNRELIABLE — VERIFY BEFORE DEPENDING.** Documented outages/SSL errors (facebookresearch/esm Discussion #627, "ESM fold prediction service is broken"), ~400-residue limit, uncertain maintenance. Do not put on the demo-critical path. |

- **Does folding change for a point mutation? Honest answer: usually NO.** Pak et al. (*PLoS ONE* 2023) found AlphaFold ill-suited to single-mutation stability effects; adversarial studies (bioRxiv 2026) show "nonphysical structural invariance." This directly undermines the original Tier 2 concept.
- **Confidence metrics:** **pLDDT** per-residue — convention >90 very high, 70–90 confident, 50–70 low, <50 very low/disordered (ESM Atlas treats pLDDT>0.7 and pTM>0.7 as "reasonably well structured"); **PAE** for relative positioning. **Surface these** — low pLDDT at the mutated residue means any structural claim is unreliable.
- **AlphaMissense pre-filter (recommended):** outputs a calibrated **pathogenicity score 0–1** and class using established thresholds (Cheng et al., *Science* 2023): **likely_benign 0.000–0.333, ambiguous 0.334–0.564, likely_pathogenic 0.565–1.000**, chosen for **90% precision on ClinVar** ("We classify 32% of all missense variants as likely pathogenic and 57% as likely benign using a cutoff yielding 90% precision on the ClinVar dataset"). It predicts functional disruption, **not** drug response.
  - **Licensing nuance:** predictions were originally **CC BY-NC-SA 4.0 (non-commercial)**; **as of March 2024 the precomputed scores were relicensed to CC BY 4.0** (dbNSFP v4.7 "c" branch; DeepMind GitHub also states CC BY 4.0). **Verify the license on the exact file/route you download** — for a non-commercial hackathon either is fine, but NC terms on some distributions matter if ever commercialized.
  - **Prefer precomputed-score LOOKUP** (DeepMind/EBI-VEP/Hugging Face `katielink/dm_alphamissense`/hegelab.org host all ~71M human scores) over running the model.
- **Low-confidence/failed fold:** fall back to the AlphaMissense score and mark structural analysis unavailable; never fabricate a number (Section 8).

### 5. Binding affinity prediction (the technical heart)
**DECISION: Do NOT use DeepPurpose for mutation effects. Use mutation-specific protein–ligand ΔΔG predictors — primarily mCSM-lig — as the V1 signal, with FoldX/DDGun stability ΔΔG and (stretch) AutoDock Vina/smina docking as complements. Report mutant-minus-wildtype DELTA.**

- **DeepPurpose evaluation:** takes drug SMILES + protein amino-acid sequence, outputs affinity on a **log scale (pKd/pIC50/pKi)**, trained on **DAVIS** (68 kinase inhibitors × 442 kinases, >80% of the human catalytic kinome; Davis et al. 2011), KIBA, BindingDB. **It cannot resolve a single-residue change** — sequence DTI models memorize protein identity and collapse in cold-target splits; the 2024 attention-DTI paper confirms models "cannot clearly distinguish mutated and wild-type proteins in the Davis kinase dataset." Mutant-vs-WT deltas from DeepPurpose are meaningless noise.

| Tool | Predicts | Access | Notes |
|---|---|---|---|
| **mCSM-lig** ⭐ | ΔΔ protein–small-molecule affinity upon single missense mutation | Web server (biosig.lab.uq.edu.au/mcsm_lig; structure.bioc.cam.ac.uk/mcsm_lig) | Trained on **Platinum**; graph-based signatures; **needs a protein–ligand complex structure**; ρ up to 0.67 |
| **FoldX** | Protein **stability** ΔΔG (not ligand-specific) | Academic license, downloadable | |ΔΔG|>1 kcal/mol ≈ destabilizing; weak on stabilizing mutations |
| **Rosetta ddG (cartesian/flex ddG)** | Stability / interface ΔΔG | Free academic, heavy | ~5–7 days to first successful run for non-experts — **avoid for hackathon** |
| **DDGun / DDGun3D** | Stability ΔΔG | Free, web + local | Fast simple baseline |
| **ThermoNet / PremPS / PremPLI** | Stability / ligand ΔΔG | Free | PremPLI is ligand-specific like mCSM-lig |
| **AutoDock Vina / smina** | Physics docking score per pose | Free, open (Apache-2.0) | Dock the same drug into WT vs. mutant structure and compare; CPU-fine; **smina** is often the best easy default; structure prep is the hard part |
| **DiffDock / GNINA** | ML docking pose | Free, GPU-preferred | More setup; DiffDock slow on CPU |

- **mCSM-lig access caveat (verify at build):** it is a **web server** (no advertised local install); for a hackathon use it interactively or confirm whether batch submission is permitted (respect terms/rate limits). It **requires a protein–ligand complex** — tying Section 4 (get a structure) to Section 5.
- **Recommended V1 (mostly no-GPU):** parse+validate mutation → **AlphaMissense lookup** (functional gate) → fetch experimental **PDB**/AlphaFold model (prefer drug-bound) → candidate drugs + PubChem SMILES → **mCSM-lig** ΔΔ affinity per drug (or **Vina/smina WT-vs-mutant docking** where no complex exists) → labeled plausibility signal with all caveats.
- **Stretch V2:** add flex ddG/Rosetta for select cases; ensemble mCSM-lig with docking; add DDGun stability as an orthogonal check; add resistance-specific databases.
- **Delta, not absolute:** use ΔΔ (the change attributable to the mutation cancels some systematic bias; mCSM-lig outputs this natively). **Pitfalls:** if the model is mutation-insensitive the delta is noise; subtracting two uncertain numbers inflates relative uncertainty; delta says nothing about whether the drug binds at all; sign conventions differ between tools — track them.
- **Thresholds (be explicit about what's grounded vs. arbitrary):** stability **|ΔΔG| > ~1 kcal/mol ≈ destabilizing** is a real *convention* (CAGI5 stability challenge; Rosetta/FoldX literature; some use 1.2–2 kcal/mol for "highly destabilizing") — cite as convention, not biological law (Rosetta reports Rosetta Energy Units needing ~2.94× to approximate kcal/mol). For **protein–ligand ΔΔ affinity there is NO established clinical threshold** — you must set and label an internal heuristic cutoff and validate it (Section 6). AlphaMissense: use the published 0.34/0.564 cutoffs but remember it's pathogenicity, not response.

### 6. Validation
**DECISION: Rediscovery benchmark on ~10–20 textbook cases held out as unknown; pre-commit pass/fail (≥70% correct directionality on clear-cut cases); if agreement is poor, narrow the claim to "mutation is in/near the binding pocket."**

**Gold-standard cases with known directionality (score against these):**

| Gene / mutation | Drug | Known direction |
|---|---|---|
| EGFR T790M | gefitinib/erlotinib (1st-gen) | **Resistance** (gatekeeper) |
| EGFR T790M | osimertinib (3rd-gen) | **Sensitive** (designed for T790M) |
| EGFR C797S | osimertinib | **Resistance** (abolishes covalent C797 binding) |
| ABL1 T315I | imatinib | **Resistance** (classic gatekeeper) |
| KIT D816V | imatinib | **Resistance** (kinase-domain; contrast juxtamembrane V560G, *more* sensitive) |
| BRAF V600E | vemurafenib | **Sensitive** (designed for V600E) |
| ALK G1202R | crizotinib | **Resistance** (sensitive to lorlatinib) |
| ALK I1171T | crizotinib | **Resistance** (sensitive to ceritinib/lorlatinib) |

- **Resistance databases for labels/reference:** COSMIC resistance mutations; CIViC (open); OncoKB (browsable free; API licensed — build reference lists, don't depend live); MdrDB (integrates GDSC/DepMap/Platinum/TKI).
- **Model benchmark datasets:** **Platinum** (biosig.lab.uq.edu.au/platinum; >1,000 mutations mapped to 3D complexes) — the appropriate benchmark for mCSM-lig-type tools; **SKEMPI** (protein–protein, for PPI tools); **ThermoMutDB/ProTherm** (stability, for FoldX/Rosetta/DDGun); **PSnpBind** (binding-site mutation → ligand affinity, REST API).
- **Feasible scale:** 10–20 hand-picked cases (each web-server run is minutes). **Pre-commit criteria before running:** e.g., "correct resistance-vs-sensitive directionality on ≥70% of clear-cut cases, with no confident wrong call on the designed-drug positives (BRAF V600E/vemurafenib, EGFR T790M/osimertinib)."
- **Published accuracy figures to quote (real numbers):**
  - **mCSM-lig:** verbatim from Pires, Blundell & Ascher, *Sci Rep* 2016 — "very good correlation with experimental data (**up to ρ = 0.67**) and is effective in predicting a range of chemotherapeutic, antiviral and antibiotic resistance mutations." (A finer within-5 Å sub-figure of ρ≈0.674 appears in circulated summaries but **could not be independently re-verified in this pass — check the paper's Figure 2 directly before quoting it**; the ρ = 0.67 headline is confirmed.)
  - **Pan, Portelli, Nguyen & Ascher, *Brief Bioinform* 2026, 27(1):bbag035:** predictions degrade **5–30%** on modeled vs. experimental structures; mCSM-lig & PremPLI strongest on Platinum-derived data (R≈0.75–0.76, RMSE≈1.06–1.07 kcal/mol), flex ddG weakest there (R≈0.19).
  - **Rosetta ddG (stability):** r ≈ 0.68–0.73 on ~1,210 mutations (Kellogg et al.; Rosetta Commons docs).
  - **FoldX (stability):** Pearson R ≈ 0.69–0.71 on large sets (VFSD ~2,275 mutations), but drops to R ≈ 0.30 on harder curated sets (S461) — quote both to show dataset dependence.
  - **DeepPurpose:** original paper reports "comparable" to DeepDTA/GraphDTA on DAVIS/KIBA and a generalization Pearson ≈ 0.78 on 1,000 unseen DAVIS pairs — **none of these are mutation-discrimination figures**, which is the task Tier 2 needs.
- **Fallback if agreement is poor:** narrow the claim to the structurally verifiable statement — *"This mutation lies within/near the known binding pocket of drug X (residue N, within Y Å)"* plus the AlphaMissense flag. Honest, still useful for triage, hard to get wrong.

### 7. Output and presentation
**DECISION: A visually distinct amber-outlined "Computational Plausibility Signal (NOT clinical evidence)" card, separated from Tier 1, ranked by predicted-effect × confidence × structural reliability.**
- **Fields:** input mutation (gene, HGVS + shorthand, UniProt accession/isoform used); reference-residue validation status (PASS); AlphaMissense score + class (with "pathogenicity, not drug response" note); target + structure source (PDB ID / AlphaFold DB, with pLDDT at the mutated residue); per-drug predicted ΔΔ direction (↓/↑/≈), magnitude, tool used, confidence; distance of the mutated residue to the binding site; method list + limitations; disclaimer.
- **Distinction from Tier 1:** different color (Tier 1 solid/green "Documented evidence"; Tier 2 amber/outlined "Computational signal"); never "recommended/indicated/effective" — use "predicted/computational signal/hypothesis/plausibility"; persistent badge: **"TIER 2 — COMPUTATIONAL PLAUSIBILITY SIGNAL — NOT A DIAGNOSIS OR PRESCRIPTION."**
- **Exact disclaimer:**
  > *"This is a computational plausibility signal generated by structural prediction tools for research and hypothesis-generation only. It is **not** clinical evidence, a diagnosis, or a treatment recommendation, and it has **not** been clinically validated for this mutation. The underlying models predict physical binding/stability changes with limited accuracy and were not designed to predict patient response. These results must be interpreted by a qualified clinician alongside established evidence before any clinical decision. Do not use this output to start, stop, or change any treatment."*
- **Ranking:** predicted effect magnitude × model confidence × structural reliability (pLDDT/site-distance), preferring approved + India-available in ties; always show uncertainty and allow sorting.

### 8. Failure handling
**DECISION: Every failure returns an explicit labeled message; inconclusive results are ALWAYS shown, never silently dropped.**
- **Failed/low-confidence structure:** *"Structural analysis unavailable or low-confidence for this region (pLDDT = X, below reliability threshold). No structure-based binding signal is reported. AlphaMissense functional score is shown as the only available computational signal."*
- **No candidate drugs (even pathway level):** *"No candidate drugs targeting this protein, its family, or its pathway were found in open databases (DGIdb/Open Targets/ChEMBL). Tier 2 cannot generate a signal for this target."*
- **Timeout/error:** *"Tier 2 analysis could not complete (timeout/error). No partial or inferred result is shown to avoid misleading output. Please retry or consult evidence-based resources."*
- **Why silent omission is dangerous:** in clinical decision support a physician cannot distinguish "found nothing concerning" from "failed/skipped." Silent drops create false reassurance and hide system limitations — always render the null/failure state explicitly.

### 9. Compute and timeline reality
**DECISION: V1 (lookups + web-server ΔΔG) runs on Colab free tier / a laptop; do NOT fold live on stage — pre-compute and cache demo cases while showing the live pipeline path for a pre-validated case.**

| Tool | Free hosted? | Local install | GPU/RAM | Latency | Beginner setup |
|---|---|---|---|---|---|
| hgvs parser | n/a | pip, easy | trivial CPU | ms | Low |
| UniProt/PubChem REST | Yes | n/a | none | <1 s (rate-limited) | Low |
| AlphaMissense lookup | Yes (files/VEP/HF) | download parquet | modest RAM | ms | Low–medium |
| AlphaFold DB / PDB fetch | Yes | n/a | none | <1 s | Low |
| ESMFold (HF/local) | Hosted unreliable | medium | GPU + high VRAM for big proteins | sec–min on GPU | Medium–high |
| ColabFold (AF2) | Free Colab | notebook | free Colab GPU; ~2,500-res cap | ~5–10 min/seq | Medium |
| mCSM-lig | Web server | not advertised | none (server-side) | sec–min/query | Low (interactive) |
| FoldX | License | download | CPU | sec–min | Medium |
| Rosetta ddG | Free academic | heavy | CPU-heavy | hours; ~5–7 days to first run | High (avoid) |
| AutoDock Vina/smina | Free | download, moderate | CPU fine | sec–min/pose | Medium (structure prep) |

- **Colab free-tier:** V1 path is trivially within capacity. Free Colab disconnects after ~90 min idle / ~12 h total with GPU quotas — never put a long GPU session on the demo-critical path.
- **Live vs. pre-computed — HYBRID:** pre-compute/cache your 10–20 validation/demo cases for instant, reliable on-stage results, while displaying the *actual* pipeline steps executing for a known-good case (genuinely live, not faked). Never attempt a first-time fold or cold mCSM-lig submission live.
- **Fallback demo plan:** if any API is down, switch to the cached set with a visible "using cached result (offline mode)" label; keep a fully offline copy of demo inputs (sequences, structures, SMILES, precomputed scores) on disk.

### 10. Ethical and accuracy framing
**DECISION: One consistent vocabulary everywhere; a candid validation answer; pre-empt CDSCO/DPDP.**
- **Terminology — always:** "computational plausibility signal," "hypothesis-generating," "structural prediction," "for physician review." **Never:** "diagnosis," "prescription," "recommended treatment," "proven," "accurate prediction of response."
- **If a judge asks "has this been clinically validated?":**
  > *"No — and we're deliberate about that. Tier 2 has not been clinically validated and is not a medical device. It's a computational hypothesis-generator that flags whether an undocumented mutation sits in a drug-binding region and whether physics-based tools predict it weakens or strengthens binding. We validated it technically by checking it rediscovers textbook resistance cases like ABL1 T315I/imatinib and EGFR T790M, but that's a sanity check on the method, not clinical validation. Every output is labeled as a non-clinical signal for physician review."*
- **Plain-language tool explanations (beginner-deliverable):**
  - **UniProt:** "The reference library of protein sequences — we look up the normal protein so we can apply the patient's mutation to the right spot."
  - **hgvs parser:** "A tool that reads mutation notation like 'R175H' correctly so we don't misinterpret which change is meant."
  - **AlphaMissense:** "A Google DeepMind AI that estimates how likely a single amino-acid change is to break the protein's function — a first filter, not a drug-response predictor."
  - **ESMFold/AlphaFold:** "AIs that predict a protein's 3D shape from its sequence. Like predicting a folded paper crane from the instructions — but they're weak at seeing how one tiny fold change alters things, which is why we don't rely on them for mutation effects."
  - **mCSM-lig:** "A tool built specifically to estimate how a mutation changes how tightly a drug grips its target protein — trained on real experimental measurements of exactly that."
  - **AutoDock Vina:** "A physics-based 'docking' simulator that fits a drug into a protein pocket and scores the fit — we compare the fit before and after the mutation."
  - **DeepPurpose:** "A deep-learning library that guesses drug–protein binding from sequences; we tested it and found it can't 'feel' a single-letter mutation, so we don't use it for that job."
- **Regulatory:** India's **CDSCO Draft Guidance on Medical Device Software (Oct 21, 2025)** sets a risk-based Class A–D framework where standalone clinical-decision-support / AI diagnostic software can be a regulated medical device (A/B licensed by State authorities, C/D centrally by CDSCO), with an Algorithm Change Protocol for AI updates. A deployed diagnostic version could fall in scope; a **hypothesis-generating, physician-in-the-loop research prototype making no diagnostic/treatment claim** is positioned to stay outside device classification (analogous to the US 21st Century Cures Act CDS carve-out for tools that inform rather than direct clinicians) — say this proactively. **DPDP Act 2023:** health data is sensitive personal data needing specific consent/safeguards — the prototype should store **no identifiable patient data**, operating on de-identified gene+mutation inputs only, and note alignment with DPDP and ICMR/NDHM guidance.

---

## Summary table — Recommended V1 tool stack

| Tool | Purpose | Access mode | Cost | Latency | Key limitation |
|---|---|---|---|---|---|
| **hgvs / Mutalyzer** | Parse & validate mutation notation | pip / free API | Free | ms | Shorthand needs pre-parsing |
| **UniProt REST** | Canonical protein sequence | REST | Free | <1 s | Isoform/numbering mismatches |
| **AlphaMissense (precomputed lookup)** | Pathogenicity pre-filter | File/VEP/HF download | Free (CC BY / verify route) | ms | Pathogenicity ≠ drug response |
| **RCSB PDB / AlphaFold DB** | Target structure | REST/download | Free | <1 s | May lack drug-bound complex |
| **DGIdb / Open Targets / ChEMBL** | Candidate drug generation | API | Free | seconds | Coverage gaps |
| **PubChem PUG REST** | Drug SMILES | REST | Free | <1 s (≤5/s, ≤400/min) | Naming/SMILES field changes |
| **mCSM-lig** ⭐ | Mutation → protein–ligand ΔΔ affinity | Web server | Free | sec–min | Needs complex structure; ρ≈0.67 ceiling |
| **AutoDock Vina / smina** | WT-vs-mutant docking (complement/stretch) | Local, open | Free | sec–min | Structure prep is the hard part |
| **FoldX / DDGun** | Stability ΔΔG cross-check | License / free | Free/academic | sec–min | Stability ≠ binding |

---

## Recommendations (staged, with decision thresholds)

**Stage 0 — before writing code (½ day):** Commit to the redesign (drop DeepPurpose + fold-difference as the engine). Hand-run the 8 gold-standard cases through the **mCSM-lig web server** using existing PDB complexes. *Threshold to proceed:* if mCSM-lig gets directionality right on ≥6/8, build the pipeline around it. If not, fall back immediately to the narrower "mutation-in-pocket + AlphaMissense" scope.

**Stage 1 — V1 build:** parse+validate (with the mandatory reference-residue gate) → AlphaMissense lookup → PDB/AlphaFold structure fetch → DGIdb/Open Targets candidates + PubChem SMILES → mCSM-lig ΔΔ affinity → labeled plausibility card. Pre-cache all demo cases.

**Stage 2 — validation & labels:** run the ≥10-case rediscovery benchmark; set the user-facing ΔΔ-affinity cutoff empirically and label it a heuristic. *Threshold:* publish the tool's own accuracy numbers (mCSM-lig ρ up to 0.67; 5–30% degradation on modeled structures) in the UI/deck so claims match evidence.

**Stage 3 — stretch (only if time):** add AutoDock Vina/smina WT-vs-mutant docking for targets lacking a complex; consider flex ddG/DDGun as orthogonal checks.

**Benchmarks that change the plan:** (a) if mCSM-lig batch access is disallowed or the server is down → pivot to Vina/smina docking deltas; (b) if directionality accuracy <70% on clear-cut cases → narrow to the "in-pocket + pathogenicity" claim; (c) if the ESM Atlas API is needed and unavailable → use AlphaFold DB/PDB only.

---

## Biggest risks and honest caveats

1. **The original design is unsound and must change.** DeepPurpose cannot detect point-mutation effects (protein-identity memorization; cold-split failure) and ESMFold/AlphaFold return near-invariant structures for single mutations — building on them yields confident-looking noise. This is the #1 risk and the reason for the redesign.
2. **Even the recommended tools are modest.** mCSM-lig tops out around ρ = 0.67 on experimental data and degrades 5–30% when structures are modeled. Treat every output as a weak prior.
3. **No validated clinical threshold exists** for protein–ligand ΔΔ affinity → user labels; the stability |ΔΔG| > 1 kcal/mol convention is a heuristic, not a clinical boundary. Say so.
4. **Structure dependency:** mCSM-lig needs a protein–ligand complex; where none exists you must dock first (adding error) or fall back to the narrower claim.
5. **Hosted-API fragility:** the ESM Atlas fold API is unreliable and web servers rate-limit — pre-cache everything for the demo.
6. **Numbering/isoform silent errors** are the most likely invisible bug — the reference-residue validation gate is mandatory.
7. **Regulatory/ethical:** under India's 2025 CDSCO software guidance a deployed diagnostic version could be a regulated medical device; keep the prototype firmly in "physician-reviewed hypothesis-generation, no patient data stored" territory to stay defensible under CDSCO and DPDP 2023.
8. **Honest scope wins credibility:** the strongest pitch is not "we predict the right drug" but "we honestly flag a computational plausibility signal where no evidence exists, with every limitation labeled." That framing is both more accurate and more compelling to informed judges.
9. **One figure to double-check before quoting:** the mCSM-lig within-5 Å sub-correlation (~ρ 0.674) could not be independently re-verified here — confirm against Figure 2 of the 2016 Scientific Reports paper before putting it in the deck; the ρ = 0.67 headline is confirmed.