# Tier 2 Maximum-Potential Roadmap: Computational Structural Drug-Response Plausibility Signal

## TL;DR
- **Tier 2 can be substantially improved, but its ceiling is fixed by physics and data, not engineering.** The single highest-value upgrade is replacing the arbitrary ΔΔ-affinity threshold with a *data-driven calibration model* trained on the MdrDB / Platinum / TKI resistance datasets that maps predicted ΔΔG to an empirical probability of biochemical resistance, wrapped in conformal prediction for honest confidence intervals. Every other upgrade is secondary to fixing this.
- **The best realistic accuracy for the core task (ΔΔ binding-affinity on a missense mutation) is roughly Pearson R≈0.75 / RMSE≈1.0 kcal/mol on data resembling the training set, degrading to R≈0.24–0.49 on independent targets** — meaning Tier 2 will always be a hypothesis-generator, never a predictor. Kinase conformational-state selection (DFG-in/out) is the single largest silent error source and must be handled explicitly.
- **Whole categories remain genuine open problems: fusion oncoproteins, compound mutations in cis/trans, bypass-pathway (non-binding) resistance, and any link to actual patient outcomes rather than biochemical IC50.** These should be flagged in-product as "out of computational scope," not papered over.

## Key Findings
- The V1 primary engine, mCSM-lig, achieves only ρ≈0.67 at best (Pires, Blundell & Ascher, *Scientific Reports* 2016, DOI 10.1038/srep29575: "our method provides a very good correlation with experimental data (up to ρ = 0.67)") and is web-server-only with no local install; it degrades further on non-experimental structures and should be supplemented or replaced, not trusted alone.
- Boltz-2 (MIT license, single-GPU) is the most important new capability: it co-folds protein–ligand complexes AND predicts affinity. Per Wohlwend/Passaro et al., *bioRxiv* 2025 (DOI 10.1101/2025.06.14.659707), it is "the first AI model to approach the accuracy of FEP methods for predicting binding affinities on the FEP+ benchmark … at least 1000× more computationally efficient than FEP" (open-sourced 6 June 2025 by MIT CSAIL/Jameel Clinic with Recursion). Independent benchmarks show it is not a lossless FEP replacement and handles allosteric, metal-coordinated, and mutation-delta cases poorly.
- Physics-based free-energy methods (open-source OpenFE/PMX) are the accuracy ceiling (~1–2 kcal/mol RMSE) but cost ~1 GPU-day per mutation edge — usable only for a small curated set of high-value actionable genes, not real-time.
- Kinase conformational state (DFG-in/out) determines type-I vs type-II inhibitor binding; standard AlphaFold returns one state and picking the wrong one is a major silent failure for exactly the ABL1/KIT/ALK/BRAF targets in V1's benchmark. AF2-RAVE / multi-state modeling can generate both states.
- India's CDSCO issued its final Guidance on Medical Device Software in mid-2026 with a four-tier risk classification and an Algorithm Change Protocol; a hypothesis-generating tool that outputs non-clinical signals sits at the boundary of SaMD and must be scoped carefully. Health data is "sensitive personal data" under the DPDP Act 2023.

## A. Beyond single missense — expanding mutation scope

### A1. In-frame insertions/deletions (indels)
- **Current limitation:** V1 hard-rejects all indels. Yet EGFR exon 19 deletions, EGFR exon 20 insertions, and KIT exon 11 indels are among the most clinically important actionable events in exactly the treatment-exhausted population Tier 2 targets.
- **Upgrade:** Build a physically relaxed mutant model for the indel, then score. For structure building, indels can be modeled by loop reconstruction (Rosetta remodel / MODELLER) or by AlphaFold2/3-class co-folding of the mutant sequence. For ΔΔG, the recent IFUM method explicitly handles indels (folding-stability ΔΔG): per the IFUM paper (*Nature Communications* 2026, DOI 10.1038/s41467-026-68637-4), "Within the full Mega-scale test set, IFUM achieved a PCC of 0.81, 0.80, and 0.63" for point mutants, indels, and double mutants respectively — vastly better than FoldX (PCC 0.22) and Rosetta (PCC 0.01), which are essentially non-predictive for indels.
- **Tools:** IFUM (folding-stability ΔΔG, handles indels); Rosetta cartesian_ddg (point only); MODELLER/Rosetta remodel for loop building; AlphaFold2/3 or Boltz for whole-domain remodeling.
- **Accuracy:** IFUM indel PCC 0.80 (stability, not ligand affinity). No ligand-ΔΔ-affinity predictor is validated for indels — mCSM-lig and PremPLI are single-substitution only. The CSM-lig authors state the framework could in principle extend to insertions/deletions/splicing, but this is not benchmarked.
- **Compute/licensing:** IFUM and Rosetta are academic-licensed; MODELLER academic-free. CPU-feasible for stability; whole-domain refolding needs a GPU.
- **Residual limitation:** You can get a *stability/structural plausibility* signal for indels but NOT a validated binding-affinity-change signal. Honest labeling required: "structural disruption estimate only."

### A2. Gene fusions (EML4-ALK, BCR-ABL1, COL1A1-PDGFB, SS18-SSX)
- **Current limitation:** Entirely out of scope. Fusions are structurally chimeric; the drug usually targets a retained kinase domain (e.g., the ALK or ABL1 kinase domain), so the *fusion breakpoint* is often irrelevant to binding while the *kinase-domain resistance mutation* on the fusion (e.g., ALK G1202R, BCR-ABL1 T315I) is what matters.
- **Defensible partial workaround:** Reframe — if the actionable event is a resistance mutation *within the retained druggable domain* of a fusion, Tier 2 can model that domain in isolation (as a normal kinase domain) and ignore the fusion partner. This covers EML4-ALK + G1202R and BCR-ABL1 + T315I. This is the single most useful and honest fusion capability.
- **Full fusion-structure modeling is an open problem:** AlphaFold2 models of PAX3::FOXO1, EWSR1::FLI1, EML4::ALK, SS18::SSX1 show largely low (50–70) and very-low (<50) pLDDT, indicating extensive intrinsic disorder (FusOn-pLM paper, *Nat. Commun.* 2025, DOI 10.1038/s41467-025-56745-6); predictions degrade specifically where a fusion partner is absent from the training database. FusionPDB predicted 2300+ fusion structures but only 68 of 266 curated fusions passed a pLDDT≥70 active-site quality gate. FusOn-pLM is a fusion-specific language model but predicts sequence/disorder features, not druggable 3D complexes.
- **Tools:** AlphaFold3/Boltz/Chai-1 for the retained domain; FusionPDB (knowledgebase of predicted fusion structures); FusOn-pLM (fusion language model).
- **Residual limitation:** Whole-fusion structural reasoning is not reliable. Restrict Tier 2 to resistance mutations inside retained, individually-foldable druggable domains and flag everything else as out of scope.

### A3. Frameshift/truncating variants and neoantigens
- **Current limitation:** Out of scope; correctly so for a *binding* signal.
- **Honest position:** A frameshift/truncation that removes the drug-binding domain trivially abolishes binding (no computation needed); one that truncates downstream may not affect the target pocket at all. Structural ΔΔ-affinity modeling is inapplicable. The defensible adjacent signal is *neoantigen/immunogenicity* prediction (e.g., NetMHCpan for MHC-I binding of the frameshift-derived peptide) — but this is a different therapeutic axis (immunotherapy), outside Tier 2's binding-plausibility purpose. Recommend an explicit "out of scope — consider immunotherapy pathway" flag rather than a fabricated binding signal.

### A4. Splice-site variants (MET exon 14 skipping)
- **Current limitation:** Out of scope.
- **Upgrade path:** Splice consequences are a *sequence/transcript* problem, not a structural one. The correct approach: predict the protein-level consequence (SpliceAI for splice disruption; then map the resulting in-frame deletion, e.g., MET exon 14 skipping removes the juxtamembrane region), then, if the result is an in-frame deletion in/near a druggable domain, hand off to the A1 indel path. MET exon 14 skipping does not alter the kinase pocket, so MET inhibitors (capmatinib/tepotinib) retain binding — a structural signal would (correctly) say "binding retained."
- **Tools:** SpliceAI (Illumina, open) for splice prediction; then the A1 pipeline.
- **Residual limitation:** The two-step chain compounds error; only defensible when the skip produces a clean in-frame deletion mappable to structure.

### A5. Compound/co-occurring mutations (EGFR T790M + C797S in cis vs trans)
- **Current limitation:** V1 models one mutation at a time; cis vs trans configuration is clinically decisive (T790M+C797S in *trans* remains sensitive to combination EGFR TKIs; in *cis* is resistant to all current EGFR TKIs).
- **Upgrade:** Build the double-mutant structure explicitly (FoldX BuildModel and Rosetta cartesian_ddg both accept multiple simultaneous mutations) and score the combined complex. Physics/alchemical methods (OpenFE/PMX) naturally handle multiple mutations. For stability, IFUM handles double mutants (PCC 0.63) and Rosetta (~0.78) and ESM2 (~0.70) do reasonably.
- **Critical caveat — cis/trans is not a structural property of a single monomer:** "in cis" (both mutations on one allele/protein molecule) vs "in trans" (different molecules) is a *genomic phasing* question. Tier 2 can model the cis case (both substitutions on one modeled chain) but the trans case requires modeling two separate protein populations and reasoning about combination therapy — partially out of structural scope.
- **Tools:** FoldX BuildModel, Rosetta cartesian_ddg, OpenFE/PMX (all accept multi-mutation input).
- **Residual limitation:** Can model cis compound effects on one chain; cannot infer phasing (must be supplied by Tier 1 genomics); trans + combination reasoning remains partial.

### A6. Copy-number / expression effects
- **Honest position:** Fully out of structural scope. The defensible adjacent signal is target-amplification awareness: if Tier 1 reports focal amplification of the drug target (e.g., MET or ERBB2 amplification), that is a documented resistance/sensitivity mechanism independent of binding, and Tier 2 should surface a clearly-separated non-structural note ("amplification detected — binding model does not capture dose effects") rather than a ΔΔG. No structural computation applies.

## B. Better structural fidelity

### B1. Physically relaxed mutant structures
- **Current limitation:** V1 uses the wild-type structure with the mutation merely "noted," not built. This misses side-chain repacking, steric clash relief, and local backbone relaxation — the actual physical basis of most resistance (e.g., T790M's bulky methionine, T315I's isoleucine gatekeeper).
- **Upgrade:** Build and minimize the mutant. FoldX BuildModel (fast, seconds, empirical) for triage; Rosetta cartesian_ddg (Cartesian-space local backbone relaxation, scaling factor 2.94 to convert Rosetta energy units → kcal/mol) for refinement; PyRosetta for scripting; OpenMM/GROMACS for energy minimization. Combine FoldX + Rosetta because their predictions overlap only 12–25%, so consensus is genuinely additive.
- **Accuracy:** For stability ΔΔG, Rosetta cartesian_ddg is robust to structural perturbation (works on homology models with ≥40% template identity); FoldX PCC 0.30–0.5 on independent sets with occasional large outliers (RMSE up to ~1.9 kcal/mol on the S461 set); cartesian_ddg can itself produce large errors (RMSE 3.59 on one benchmark) — neither is reliable in isolation.
- **Compute/licensing:** FoldX free for academics (licensed binary); Rosetta/PyRosetta free academic license, CPU-bound, minutes-to-hours per variant; OpenMM (MIT) and GROMACS (LGPL) free.
- **Residual limitation:** Relaxation improves realism but ΔΔG error stays ~1 kcal/mol at best; it fixes local geometry, not conformational-state or dynamics errors.

### B2. Co-folding models (Boltz-2, Chai-1, AlphaFold3, OpenFold3, Protenix, ESM3)
- **Current limitation:** V1 has no way to generate a protein–ligand *complex* structure when none exists in the PDB; it falls back to docking.
- **Upgrade:** Use an open co-folding model to jointly predict the mutant complex and (Boltz-2) its affinity.
- **Licensing (critical, verified):**
  - **Boltz-1 and Boltz-2: MIT license** — fully open, commercial use permitted.
  - **Chai-1: released Apache 2.0** (moved from an initially restrictive license); Chai-1r variant for unrestricted commercial use.
  - **Protenix (ByteDance): Apache 2.0.**
  - **OpenFold3: Apache 2.0** (academic + commercial), preview stage as of this writing.
  - **AlphaFold3: NOT open for commercial use** — code CC-BY-NC-SA 4.0, weights only to academics on request with a prohibited-use policy barring commercial/clinical activity. **Do not use AF3 in a product intended for clinical deployment.**
  - **ESM3: restricted/community license** with commercial tiers — check terms.
- **Boltz-2 affinity claim, evaluated:** Boltz-2 approaches FEP+ accuracy at ~1000× lower cost and outperformed all submitted methods on the CASP16 affinity challenge (140 complexes). BUT independent scrutiny (Rowan; PREreview by the Fraser lab) finds: predicted ligand geometry is highly variable (some poses >10 Å RMSD); allosteric binders poorly described; metal coordination, ordered waters, and multimeric cofactors not modeled well; membrane targets (GPCRs) show higher variance; and on PDBbind/DUD-E it generalizes only comparably to GNINA and end-state methods. It reports affinity as log10(IC50).
- **Compute:** Boltz-2/Chai-1 run on a single modern GPU (feasible beyond Colab-free; a ~24–40 GB GPU is realistic). AF3-class inference is heavier.
- **Residual limitation:** Co-folding gives you a plausible *static* complex and a *population-level* affinity, but it is NOT mutation-delta-sensitive in the way Tier 2 needs — it was not trained to resolve the small ΔΔG between wild-type and single-mutant. Use it for structure generation and as one ensemble member, not as the resistance oracle.

### B3. Molecular dynamics and ensemble docking
- **Current limitation:** A static single structure ignores flexibility; a mutation may act by shifting conformational equilibria.
- **Upgrade:** Short MD (OpenMM/GROMACS) to generate a conformational ensemble, then ensemble docking / MM-GBSA rescoring across snapshots. Meaningful when the pocket is flexible or the mutation is allosteric.
- **Accuracy/compute:** MD meaningfully helps for buried-water and flexible-loop cases but adds large cost — enhanced sampling (REST2, GCMC) reaches FEP-class RMSE (~0.9 kcal/mol) but at GPU-days per system. For Tier 2's compute class this is only justifiable for a precomputed library of top actionable genes.
- **Residual limitation:** Cost makes MD impractical for real-time arbitrary queries; reserve it for offline precomputation.

### B4. Kinase conformational states (DFG-in/out, αC-in/out) — THE major silent error
- **Current limitation:** V1 takes whatever PDB/AlphaFold structure it finds. AlphaFold returns predominantly the active DFG-in state. Type-II inhibitors (imatinib, nilotinib, sorafenib) bind the *inactive DFG-out* state. Docking a type-II inhibitor into a DFG-in structure — or scoring a resistance mutation in the wrong state — produces confidently wrong answers. This is a systematic, silent failure mode for exactly the kinase targets (ABL1, KIT, ALK, BRAF) in V1's benchmark.
- **Upgrade:** Explicitly detect/select the conformational state and, where needed, generate the missing state. Kincore/KLIFS classify kinase structures by DFG/αC state; AF2-RAVE and AlphaFold2 multi-state modeling (reduced-MSA-depth or state-specific templates) generate DFG-out models; DFGmodel builds DFG-out via multi-template MODELLER.
- **Accuracy:** ML classifiers label DFG state with near-perfect accuracy on crystal structures (one study mislabeled only 17/8826). AF2 multi-state modeling improves docking/virtual screening on kinases; AF2-RAVE recovers DFG-in/out relative stabilities and correctly predicts mutation-induced flipping ~2–3 orders of magnitude faster than brute-force MD. Biology caveat: DFG-out is not a perfect proxy for type-II binding (some type-I inhibitors bind DFG-out).
- **Tools:** KLIFS/OpenCADD-KLIFS (Python, open), Kincore, AF2-RAVE, DFGmodel.
- **Residual limitation:** Generating the right state adds complexity and its own error; but *ignoring* state is worse. This upgrade converts a silent error into a stated, managed uncertainty — arguably the highest-fidelity structural improvement available.

### B5. Missing loops, disorder, PTMs, cofactors, multi-chain/allosteric context
- **Current limitation:** V1 does not systematically handle missing residues, ATP/Mg²⁺ cofactors, or allosteric/multi-chain context.
- **Upgrade:** Loop building (Rosetta/MODELLER); include catalytic cofactors (ATP·Mg²⁺) explicitly since many kinase inhibitors are ATP-competitive and pocket geometry depends on Mg²⁺; flag disordered (low-pLDDT) regions as unmodelable; note allosteric/dimer context (e.g., BRAF acts as a dimer).
- **Residual limitation:** Co-folding models handle metals/cofactors poorly (per B2); PTM effects are largely unmodeled. State these as caveats.

## C. Better binding-affinity-change prediction (the analytical core)

### C1. Comparative accuracy of ΔΔ-affinity predictors
- **Current limitation:** mCSM-lig alone (ρ≈0.67 best-case), web-only, and its accuracy collapses on non-experimental structures.
- **Verified figures (kcal/mol, ΔΔG_bind):**
  - **mCSM-lig:** original paper up to ρ=0.67 (0.73–0.77 within 5 Å of the ligand after outlier removal); classification accuracy 0.76, AUC~0.8. Independent re-test (PremPLI paper): PCC 0.63, RMSE 2.19.
  - **PremPLI:** PCC 0.73, RMSE 1.07 on the mCSM-lig training set; the authors report it is "significantly better than mCSM-lig and not very sensitive to the initial input structures."
  - **Rosetta flex ddG:** competitive on independent kinase (TKI) data.
  - **Systematic 2026 benchmark** (Pan, Portelli, Nguyen & Ascher, *Briefings in Bioinformatics* 27(1):bbag035, Jan 2026, DOI 10.1093/bib/bbag035; 791 mutations / 221 complexes): On the Platinum subset (training-overlapping), mCSM-lig R=0.76/RMSE=1.07, PremPLI R=0.75/RMSE=1.06, flexddg R=0.19/RMSE=1.80. On independent TKI data all drop to R≈0.43–0.49; on AIMMS to R≈0.24–0.25. Using predicted structures instead of experimental costs 5–30%: ~5% for AlphaFold3 or experimental-receptor+docking, 10–20% for AlphaFold2+docking, >30% for low-identity homology models. Degradation is worst for interface mutations and low-molecular-weight ligands. The paper states scoring-function methods "in general did not present a good capability to identify putative drug-resistant mutations," whereas the mutation-aware methods (mCSM-lig, PremPLI, flexddg) did.
- **Physics/alchemical (accuracy ceiling):** FEP+ (commercial, Schrödinger) ~1 kcal/mol; the FEP+ benchmark median RMSE is 1.08 kcal/mol against an experimental reproducibility floor of 0.85–0.91 kcal/mol. Open-source OpenFE: public-dataset weighted RMSE 1.73 [1.53–1.96] kcal/mol across 58 systems (10 sub-kcal/mol); private-dataset 2.44 kcal/mol. PMX (GROMACS, open) performs "on par with commercial." MM-GBSA/PBSA rescoring is cheaper but weakly correlated on PDBbind.
- **Tools/licensing:** mCSM-lig (web only, no local install), PremPLI (web), Rosetta flex ddG (local, academic), OpenFE (MIT), PMX (open), Alchemlyb/BioSimSpace (open).
- **Residual limitation:** Even the ceiling (~1 kcal/mol) exceeds the ΔΔG differences that separate clinical sensitivity from resistance in some cases; the ML tools were trained largely on Platinum and do not generalize well to novel targets.

### C2. Consensus/ensemble scoring
- **Upgrade:** Combine MD-based free energy + Rosetta flex ddG; published work found the consensus beat both individual methods in ~78% of cases and dropped RMSE below 1 kcal/mol (0.82) on high-reproducibility data. For Tier 2, combine mutation-aware ML (PremPLI) + physics (flex ddG) + co-folding (Boltz-2) as an ensemble.
- **Residual limitation:** Ensembles reduce variance, not systematic bias shared across methods (all trained on similar data).

### C3. Uncertainty quantification
- **Current limitation:** V1's three-way label rests on a heuristic flag, not a calibrated interval.
- **Upgrade:** Conformal prediction (distribution-free coverage guarantee) on top of the ΔΔG predictor to emit genuine prediction intervals; ensemble variance across the C2 methods as an epistemic-uncertainty proxy; Gaussian-process regression for calibrated intervals. Conformal prediction is well-validated in drug-property prediction but its coverage guarantee degrades under covariate shift (novel chemotypes/targets) — use Mondrian/covariate-shift-corrected variants (e.g., CoDrug).
- **Residual limitation:** Conformal coverage is only valid under exchangeability; treatment-exhausted rare-cancer mutations are precisely the out-of-distribution regime where guarantees weaken. Report this explicitly.

### C4. Threshold calibration to clinical resistance (highest-value upgrade)
- **Current limitation:** The label boundary is "admittedly arbitrary." This is the weakest link in the whole system.
- **Upgrade:** Train a calibration model mapping predicted ΔΔG → probability of resistance, using resistance datasets with quantitative fold-change IC50: **MdrDB** (Yang, Ye, Qiu et al., *Communications Chemistry* 6:123, 2023, DOI 10.1038/s42004-023-00920-7, Tencent Quantum Lab) — "MdrDB contains 100,537 samples, generated from 240 proteins (5,119 total PDB structures), 2,503 mutations, and 440 drugs" — integrating Platinum, TKI, GDSC, DepMap, AIMMS, RET, KinaseMD; **Platinum** (~1000 experimental protein–ligand ΔΔ-affinity records with structures); **TKI** (144 Abl kinase mutants with ΔpIC50); **GDSC** (4.4 million IC50 values across ~1000 lines / 518 drugs); **DepMap/CCLE**. The Pan et al. benchmark already defines resistance cutoffs at log(2/5/10/20-fold) IC50 change and picks thresholds by maximizing MCC — the exact recipe to adopt.
- **Feasibility:** Yes, a calibration model CAN be trained; MdrDB explicitly demonstrates enhanced ML ΔΔG-prediction performance. This converts the arbitrary threshold into an empirically grounded, per-drug-class probability.
- **Residual limitation:** These are *biochemical/cell-line* fold-changes, not patient outcomes (see F5). The calibrated probability is "probability of biochemical resistance," which must be labeled as such — a proxy for, not a measure of, clinical failure. Note also that MdrDB is academic-use-only and its data is current only to 2022.

## D. Beyond binding — mechanism completeness

### D1. Non-binding resistance (bypass, downstream, efflux, amplification)
- **Current limitation:** V1 models only target binding. Many/most acquired resistances are non-binding: bypass-pathway activation (e.g., MET amplification bypassing EGFR), downstream mutations (KRAS/PIK3CA), drug efflux, target amplification. A "binding retained" signal can be clinically wrong if resistance is by bypass.
- **Honest scope statement (mandatory in-product):** "Tier 2 models only drug–target binding. It cannot detect bypass, downstream, efflux, or amplification-mediated resistance. A 'binding retained' result does NOT mean the drug will work."
- **Partial adjacent signal:** Network/pathway propagation over protein–protein-interaction graphs — STRING, Reactome, OmniPath — to flag whether co-occurring alterations sit in known bypass pathways; synthetic-lethality prediction for combination hypotheses.
- **Tools:** OmniPath (Python/R, open, integrates STRING/Reactome/SIGNOR), STRING (open), Reactome (open), NetworkX/igraph for propagation.
- **Residual limitation:** This is a coarse, literature-derived flag, not a quantitative prediction; keep it clearly separate from the structural signal.

### D2. Protein stability/expression effects
- **Current limitation:** A destabilizing mutation may reduce protein abundance (loss of function / reduced target) rather than block the drug — a different mechanism V1 conflates with "reduced binding."
- **Upgrade:** Compute stability ΔΔG in parallel: FoldX; DDGun/DDGun3D (untrained baseline, PCC ~0.45–0.49, well-behaved antisymmetry); ThermoNet (3D-CNN, good antisymmetry); PremPS (balanced training, strong on stabilizing mutations); Rosetta cartesian_ddg. Report stability and binding effects as distinct axes.
- **Accuracy:** Best current stability predictors reach PCC ~0.5–0.7 with RMSE ~1.0–1.5 kcal/mol; all are weaker on stabilizing than destabilizing mutations; consensus of ≥2 tools recommended (FoldX/Rosetta overlap only 12–25%).
- **Tools/licensing:** DDGun (open; web + standalone + docker), ThermoNet (open), PremPS (web), FoldX (academic), Rosetta (academic).
- **Residual limitation:** Stability ΔΔG does not directly give abundance; it is a proxy.

### D3. Allosteric-site detection / distant-mutation effects
- **Current limitation:** V1 only considers mutations at/near the pocket; allosteric mutations (distant residues that shift pocket conformation) are invisible.
- **Upgrade:** Allosteric analysis — AlloSigMA (allosteric signaling), PASSer (ML allosteric-site prediction), Allosteric Database (ASD), OHM, dynamic network analysis on MD trajectories.
- **Residual limitation:** Allosteric prediction is qualitative and low-precision; treat as a flag ("mutation is distal to pocket — possible allosteric effect, low confidence").

### D4. Cryptic/new pocket detection from a mutation
- **Upgrade:** Flag mutation-created druggable pockets — PocketMiner (GNN; per Meller et al., *Nat. Commun.* 2023, DOI 10.1038/s41467-023-36699-3, "it accurately identifies cryptic pockets (ROC-AUC: 0.87) >1,000-fold faster than existing methods" — vs simulation-based CryptoSite at ROC-AUC 0.85); fpocket (open, fast geometric); CryptoBank/sequence-PLM predictors (PR-AUC ~0.8 at >20% identity).
- **Tools/licensing:** PocketMiner (open), fpocket (open, GPL).
- **Residual limitation:** Cryptic-pocket prediction is qualitative; a "new pocket" flag is a research lead, not an actionable target.

## E. Candidate-drug-generation improvements

### E1. Structure-based virtual screening / repurposing over the mutant
- **Current limitation:** The DGIdb→Open Targets→ChEMBL cascade is target-centric and blind to the mutant structure.
- **Upgrade:** Dock approved-drug libraries against the *mutant* structure — GNINA (CNN-rescored, strongest single method: AutoDock-GNINA reported the highest median EF1% on LIT-PCBA), smina, AutoDock Vina/Vina-GPU, Uni-Dock (GPU batch), DiffDock-L (ML; best RMSD but physically-invalid poses on novel targets). Libraries: DrugBank approved subset, ChEMBL approved drugs, ZINC in-stock, Broad Drug Repurposing Hub.
- **Accuracy (verified):** On the harder PoseBusters set of unseen complexes, classical AutoDock Vina/GOLD remain best when physical plausibility (PB-valid) is required; DiffDock's headline RMSD success (38% top-1 <2 Å on PDBbind) drops sharply once PB-validity is enforced and it generalizes poorly to novel sequences. GNINA+AutoDock is the pragmatic, CPU/GPU-feasible choice (~2 ligands/sec on an RTX 3090).
- **Residual limitation:** Docking scores rank poorly as absolute affinities; use for pose generation + triage, then rescore with mutation-aware/physics methods.

### E2. Chemical-similarity analog expansion
- **Upgrade:** For a drug that a resistance mutation defeats, find analogs that may evade it — RDKit Morgan fingerprints + Tanimoto similarity; ChEMBL similarity search. Next-generation TKIs that overcome gatekeeper mutations often share scaffolds.
- **Tools:** RDKit (BSD, open), ChEMBL web services.
- **Residual limitation:** Similarity ≠ retained activity against the mutant; must re-score each analog structurally.

### E3. Combination approaches
- **Honest position:** Tier 2 should NOT propose specific drug combinations as actionable — synergy prediction is unreliable and combinations carry toxicity/interaction risk far beyond a binding model's scope. The defensible maximum is surfacing *mechanistically-motivated hypotheses* (e.g., "bypass pathway X active → combination targeting X is a documented research direction"), flagged as literature context only.

### E4. India-availability filtering (CDSCO)
- **Current limitation:** No India-availability filter; a US-approved repurposing candidate may be unavailable in India.
- **Data gap (flag explicitly):** There is no clean machine-readable CDSCO approved-drug API. Workarounds: CDSCO published approved-drug lists (PDFs requiring extraction), the Jan Aushadhi generic-medicine product list, and NPPA (National Pharmaceutical Pricing Authority) pricing data (contains marketed formulations). Build a curated, versioned local table by extracting these PDFs; treat availability as advisory metadata, not a hard gate.
- **Residual limitation:** Manual PDF extraction is brittle and will drift; needs periodic refresh and is a maintenance burden.

## F. Validation, benchmarking, calibration at scale

### F1. From 9 cases to systematic benchmarking
- **Current limitation:** 9 hand-picked rediscovery cases cannot estimate real-world performance or false-positive rate.
- **Upgrade:** Benchmark on public datasets — **Platinum, PSnpBind, MdrDB, SKEMPI 2.0** (protein–protein ΔΔG), **ThermoMutDB** (stability), **GDSC/DepMap/CCLE** (drug response + mutations), **OncoKB resistance annotations** (FDA-recognized; Levels R1/R2 for resistance), **COSMIC** resistance mutations, **CIViC** predictive evidence items (already in Tier 1).
- **Residual limitation:** Most of these are biochemical/cell-line, not patient (F5).

### F2. Temporal (retrospective) validation
- **Upgrade:** Hold out mutations first reported after a cutoff date and test whether the pipeline would have predicted them — a genuine prospective-analog. OncoKB/CIViC/COSMIC are timestamped, enabling this. Enough data exists for well-studied kinases (ABL1, EGFR, ALK, KIT) but is thin for rare targets.
- **Residual limitation:** Temporal validation is only meaningful where enough post-cutoff mutations exist — i.e., not for the rare cancers Tier 2 most wants to serve.

### F3. Metrics and sample sizes
- **Upgrade:** Report AUROC and AUPRC (AUPRC preferred given class imbalance — resistant mutations are the minority), Spearman/Pearson for ΔΔG regression, MCC for the binary call, and precision@k for the candidate-ranking use case. Follow the Pan et al. cutoff scheme (log 2/5/10/20-fold) and MCC-maximizing thresholds. Aim for hundreds of mutations per target class for stable estimates.

### F4. Negative controls (does it call everything resistant?)
- **Current limitation:** No test that the pipeline discriminates rather than defaulting to "resistant." The datasets are biased: Pan et al. report a mean ΔΔG_bind of +0.58 kcal/mol (skewed toward reduced binding), which can inflate apparent accuracy.
- **Upgrade:** Include documented *sensitive* mutations (e.g., EGFR L858R + osimertinib, which retains sensitivity) and neutral passenger mutations as negative controls; report specificity and base-rate-corrected predictive values.
- **Residual limitation:** Curated negatives are scarcer than positives; specificity estimates will be noisier.

### F5. Link to actual patient outcomes — a FUNDAMENTAL GAP
- **Honest position (flag as open problem):** There is no public dataset robustly linking a *computational structural prediction* to *individual patient treatment outcomes*. GDSC/DepMap/CCLE are cell lines; Platinum/TKI/SKEMPI are biochemical; OncoKB/CIViC/COSMIC record population-level clinical associations, not per-prediction validation. AACR Project GENIE (with BPC outcomes) links genomics to some outcomes but not to structural predictions. **Tier 2 can be validated against biochemical/cell-line resistance, never (with current public data) against whether a specific patient responded.** This is a hard ceiling on evidentiary maturity and must be stated in every output.

## G. Output, uncertainty communication, and safety

### G1. Automation bias — the core safety risk
- **Evidence:** Automation bias (over-reliance on automated advice) is documented across clinical decision support. A controlled computational-pathology study (Rosbach et al., "Automation Bias in AI-Assisted Medical Decision-Making under Time Pressure in Computational Pathology," arXiv 2411.00998; MICCAI/Springer 2026; n=28 pathology experts) found that "AI integration led to a statistically significant increase in overall performance, [but] it also resulted in a 7% automation bias rate, where initially correct evaluations were overturned by erroneous AI advice." Expertise modulates reliance: novices over-rely, experts sometimes under-rely. Fluent natural-language explanations *increase* perceived credibility and can worsen over-reliance (confirmation bias in XAI studies). JAMA (Jabbour et al. 2023) documents AI-driven CDS harm risk.
- **Design implications:** (a) Never present a bare label; always show method, confidence interval, and scope limits. (b) Avoid over-fluent "explanations" that manufacture false confidence. (c) Force an interaction that surfaces the "no clinical evidence — computational hypothesis only" framing. (d) Consider deliberately withholding a single point label in low-confidence/out-of-distribution cases.

### G2. Confidence expression
- **Recommendation:** Use qualitative bands anchored to the calibrated probability (C4) PLUS an explicit numeric interval PLUS an explicit "do NOT use for treatment decisions" statement. Avoid single deterministic labels. Show the calibration provenance (which dataset, which drug class).

### G3. Show the reasoning chain / structural evidence
- **Upgrade:** Render the mutated residue relative to the binding pocket in an embedded web viewer so the physician can judge proximity/plausibility themselves — Mol* (PDBe/RCSB, open, the current standard, pdbe-molstar plugin), NGL Viewer (MIT), 3Dmol.js (BSD). Display pLDDT/confidence coloring and the mutated-residue–ligand distance.
- **Benefit:** Supports physician judgment and counteracts automation bias by making the evidence inspectable rather than opaque.

### G4. Audit logging and reproducibility
- **Upgrade:** Version and log every input, model version, database snapshot (CIViC/OncoKB/AlphaMissense/PDB release), tool version, threshold, and random seed per prediction, so any result can be exactly regenerated. This is also a CDSCO/DPDP compliance requirement (H, I).

## H. Engineering / infrastructure maturity

### H1. Replace web-server dependencies
- **Current limitation:** mCSM-lig has no local install — a single point of failure, a privacy problem (sends patient-derived mutations to an external server, a DPDP concern), and non-reproducible.
- **Upgrade:** Self-host locally-runnable equivalents: Rosetta flex ddG / cartesian_ddg, FoldX (local binary), DDGun (docker image), OpenFE/PMX, Boltz-2, GNINA/smina/Vina, RDKit, fpocket, PocketMiner. PremPLI/mCSM-lig remain web-only — either drop them or use only as optional cross-checks with explicit privacy warnings.
- **Residual limitation:** Some best-in-class tools (mCSM-lig, PremPLI, FEP+) cannot be self-hosted; accept reduced access or replace.

### H2. Caching / precomputation
- **Upgrade:** Precompute lookup tables for known actionable genes (OncoKB curates over 165 actionable genes and 130+ associated drugs) and their documented resistance mutations, so common queries are instant and only novel mutations trigger live computation. Cache structures, docking poses, and ΔΔG results keyed by (protein, mutation, drug, tool-version).

### H3. Containerization / workflow management
- **Upgrade:** Docker + conda-lock for reproducible environments; Nextflow or Snakemake to orchestrate the structural pipeline (parse→validate→structure→state-select→relax→score→calibrate→report) with per-step provenance capture.

### H4. Compute scaling / cost
- **CPU-feasible:** parsing/validation, FoldX, Rosetta ddG, Vina/smina docking, RDKit, fpocket, network analysis.
- **GPU-needed:** Boltz-2/Chai-1 co-folding, GNINA, DiffDock, AF2-RAVE, MD, OpenFE/PMX free energy.
- **Realistic cost beyond free tier:** cloud single-GPU (e.g., A100/40 GB) roughly on the order of $1–4/GPU-hour on-demand; an OpenFE RBFE edge ≈1 GPU-day (community reports cite ~$12–50 per edge). Reserve GPU-heavy physics for offline precomputation of high-value genes; keep the real-time path CPU + light GPU.

### H5. Failure-mode monitoring / drift detection
- **Upgrade:** Monitor database-version drift (CIViC/OncoKB/AlphaMissense/PDB updates can change outputs); track input-distribution drift (queries drifting from the calibration distribution → widen conformal intervals or abstain); alert on reference-residue mismatch rates and low-pLDDT fractions.

## I. Regulatory & evidentiary maturity (India)

### I1. CDSCO Software-as-Medical-Device
- **Current framework:** CDSCO issued a Draft Guidance on Medical Device Software on 21 October 2025 and a final Guidance Document on Medical Device Software (Doc No. CDSCO/MD/GD/MDSW/01/2026) in mid-2026, under the Medical Devices Rules 2017. It distinguishes SiMD vs SaMD, applies a four-tier risk classification (Class A–D; Class A/B licensed by State authorities, C/D by CDSCO centrally, per the draft's licensing pathways), explicitly covers AI/ML-based MDSW across the lifecycle, and introduces an Algorithm Change Protocol (ACP) permitting bounded pre-specified model updates without a fresh licence.
- **Where Tier 2 sits:** Software that provides clinical decision support / prediction meeting the MDR device definition is regulated. Tier 2's deliberate framing as a *non-clinical hypothesis-generator that never recommends treatment* is the key determinant of whether it falls inside SaMD or stays a research tool — but this boundary is legally untested and the intended-use statement must be drafted carefully with regulatory counsel. If deployed to influence clinical decisions, expect Class C/D classification given "serious situation" use, requiring the full technical file, QMS (ISO 13485), IEC 62304 software-lifecycle compliance, clinical evaluation, and post-market surveillance.
- **Evidence required to become clinical-grade:** documented analytical validation (the F-section benchmarks), clinical validation against outcomes (blocked by F5's data gap), risk management (ISO 14971), and the ACP if the model updates.

### I2. Data protection (DPDP Act 2023)
- Health data is sensitive personal data; processing requires specific consent and DPDP-aligned handling, coordinated with ABDM and ICMR guidance. This directly forbids sending patient mutations to external web servers (reinforces H1) and mandates the audit logging in G4.

### I3. Model documentation / reporting standards
- **Adopt:** TRIPOD+AI (Collins et al., *BMJ* 2024; 27-item checklist for prediction-model reporting, supersedes TRIPOD 2015); DECIDE-AI (early-stage clinical evaluation of AI decision-support, *Nat. Med.* 2022); CONSORT-AI (if any interventional trial); MI-CLAIM; model cards documenting training data, intended use, out-of-scope uses, metrics, and known failure modes; PROBAST-AI for risk-of-bias self-assessment.

## Summary table of recommended tools

| Tool | Function | License | Local-runnable | Compute | Published accuracy |
|---|---|---|---|---|---|
| mCSM-lig | ΔΔ ligand-affinity (missense) | Web only | No | Web | ρ≤0.67; indep. PCC 0.63 / RMSE 2.19 |
| PremPLI | ΔΔ ligand-affinity (missense) | Web only | No | Web | PCC 0.73 / RMSE 1.07 (train set) |
| Rosetta flex ddG / cartesian_ddg | Mutant build + ΔΔG | Free academic | Yes | CPU, min–hr | Robust; RMSE ~1 (variable) |
| FoldX BuildModel | Fast empirical ΔΔG/stability | Free academic | Yes | CPU, sec | PCC 0.3–0.5; outlier-prone |
| DDGun/DDGun3D | Stability ΔΔG baseline | Open | Yes (docker) | CPU | PCC ~0.45–0.49 |
| ThermoNet | Stability ΔΔG (3D-CNN) | Open | Yes | GPU | Good antisymmetry |
| PremPS | Stability ΔΔG | Web | No | Web | Strong on stabilizing muts |
| IFUM | Stability ΔΔG incl. indels/doubles | Academic | Yes | CPU/GPU | Indel PCC 0.80; double 0.63 |
| Boltz-2 | Co-fold complex + affinity | MIT | Yes | 1 GPU | Approaches FEP; CASP16 best; geometry variable |
| Chai-1 / Chai-1r | Co-fold complex | Apache 2.0 | Yes | 1 GPU | AF3-level |
| Protenix | Co-fold complex | Apache 2.0 | Yes | GPU | Claims > AF3 on some benchmarks |
| OpenFold3 | Co-fold (AF3 reproduction) | Apache 2.0 | Yes | GPU | Preview |
| AlphaFold3 | Co-fold complex | CC-BY-NC-SA; weights gated | Restricted | GPU | SOTA but NOT commercial |
| OpenFE | Open FEP (RBFE) | MIT | Yes | GPU-days | RMSE 1.73 kcal/mol (public set) |
| PMX | Open FEP (GROMACS) | Open | Yes | GPU-days | On par with commercial |
| FEP+ | Commercial FEP | Proprietary | Yes | GPU | ~1 kcal/mol (median RMSE 1.08) |
| GNINA | CNN docking | Open (Apache/GPL) | Yes | GPU/CPU | Best single-method EF1% |
| smina / AutoDock Vina | Docking | Open | Yes | CPU | Best PB-valid on novel targets |
| DiffDock-L | ML docking | Open | Yes | GPU | 38% <2 Å RMSD; poor PB-valid/generalization |
| Uni-Dock | GPU batch docking | Open | Yes | GPU | High-throughput |
| PocketMiner | Cryptic-pocket GNN | Open | Yes | GPU | ROC-AUC 0.87 |
| fpocket | Geometric pocket detection | Open (GPL) | Yes | CPU | Fast baseline |
| AF2-RAVE | Kinase DFG state generation | Open | Yes | GPU | Recovers DFG-in/out stability |
| KLIFS / OpenCADD-KLIFS | Kinase state classification | Open | Yes | CPU | Near-perfect on crystals |
| SpliceAI | Splice consequence | Open | Yes | GPU/CPU | SOTA splice |
| FusionPDB / FusOn-pLM | Fusion structures / LM | Open | Yes | GPU | Low-confidence structures |
| OmniPath / STRING / Reactome | Pathway/network context | Open | Yes | CPU | N/A (annotation) |
| RDKit | Cheminformatics/fingerprints | BSD | Yes | CPU | N/A |
| Mol* / NGL / 3Dmol.js | Web structure viewer | MIT/open | Yes | Browser | N/A |
| MdrDB / Platinum / TKI / GDSC / DepMap | Calibration/benchmark data | Open (varies) | N/A | N/A | 100,537 / ~1k / 144 / 4.4M records |
| OncoKB / COSMIC / CIViC | Resistance annotations | Open/tiered | API/DL | N/A | FDA-recognized (OncoKB) |
| Docker / conda-lock / Nextflow / Snakemake | Reproducible pipeline | Open | Yes | N/A | N/A |

## OPEN PROBLEMS — no adequate solution exists

1. **Fusion oncoprotein whole-structure reasoning.** AlphaFold2/3 produce low-pLDDT, disordered models for fusion junctions; no tool gives reliable druggable-complex structures for chimeric proteins. *Closest workaround:* model only the retained, individually-foldable druggable domain (ALK/ABL1 kinase domain) and treat resistance mutations there as ordinary kinase-domain mutations; flag all else out of scope. *What would need to be built:* fusion-specific co-folding trained on experimentally-resolved fusion structures (which barely exist).

2. **Validated ligand-ΔΔ-affinity prediction for indels.** IFUM handles indel *stability* (PCC 0.80) but no tool is validated for indel *binding-affinity change*. *Workaround:* report a structural-disruption/stability signal only. *Needed:* an indel-capable, experimentally-benchmarked ΔΔ-affinity method and the training data to build it.

3. **Compound-mutation cis/trans phasing.** Structural tools can model a cis double-mutant on one chain but cannot determine phasing, which is genomic; the trans + combination-therapy case is only partially structural. *Workaround:* require phasing from Tier 1 genomics; model cis only. *Needed:* integration of phased sequencing with structural modeling plus combination-response models.

4. **Non-binding (bypass/downstream/efflux/amplification) resistance.** Tier 2's binding model is structurally blind to the majority mechanism of acquired resistance. *Workaround:* coarse pathway/network flags (OmniPath/STRING/Reactome) kept separate from the structural signal. *Needed:* validated multi-omic network models of resistance — an unsolved research frontier.

5. **Linking computational predictions to real patient outcomes.** No public dataset connects a structural prediction to individual patient response; all validation is biochemical/cell-line. *Workaround:* validate against biochemical/cell-line resistance and label outputs as a "biochemical-plausibility proxy." *Needed:* prospective registries linking genotype + computational prediction + treatment + outcome — none exist openly.

6. **Sub-kcal/mol ΔΔG accuracy on novel targets.** The physics ceiling (~1 kcal/mol) and ML generalization collapse (R 0.75→0.24 off-training-distribution) mean the exact rare/novel targets Tier 2 serves are where accuracy is worst. *Workaround:* conformal intervals + abstention out-of-distribution. *Needed:* fundamentally more accurate, generalizable free-energy methods or vastly larger experimental training data for rare targets.

7. **Conformational-state ground truth for arbitrary kinases.** AF2-RAVE/multi-state modeling generate DFG-out models but cannot guarantee the correct state population for a novel mutant; wrong-state scoring is a silent error. *Workaround:* generate both states, score both, report divergence as uncertainty. *Needed:* reliable per-mutant conformational-ensemble prediction.

8. **Machine-readable CDSCO drug availability.** No API for India-approved drugs; only PDFs. *Workaround:* curated, versioned local extraction from CDSCO/Jan Aushadhi/NPPA. *Needed:* an official structured CDSCO drug database.

## Recommendations (staged by capability, not time)

**Stage 1 — Fix the analytical core (do first; highest fidelity gain per unit effort):**
1. Replace the arbitrary label threshold with a **data-driven calibration model** (C4) trained on MdrDB/Platinum/TKI, emitting a calibrated *probability of biochemical resistance* per drug class.
2. Add **conformal prediction** intervals (C3) and abstention when out-of-distribution.
3. Add **explicit kinase conformational-state handling** (B4) — this eliminates the largest silent error for V1's own benchmark targets.
4. Self-host the ΔΔG engine (Rosetta flex ddG + FoldX consensus) to remove the mCSM-lig web dependency (H1) — also a DPDP requirement.
- *Threshold to proceed:* the calibration model achieves AUPRC materially above base rate on held-out TKI/AIMMS data and passes temporal validation on ABL1/EGFR.

**Stage 2 — Broaden fidelity and scope where defensible:**
5. Build relaxed mutant structures (B1); add stability ΔΔG as a separate axis (D2).
6. Add co-folding (Boltz-2, MIT) for complex generation and as an ensemble member (B2), never as sole oracle.
7. Add the **in-scope subset** of new mutation types: cis compound mutations (A5), fusion resistance mutations in retained domains (A2), MET-exon-14-style clean in-frame skips (A4) — each with explicit scope labels.
8. Add structure visualization (Mol*) and the mandatory scope/automation-bias-aware output design (G).
- *Threshold:* each new mutation class is validated on its own held-out benchmark before enabling; disable any class that cannot beat a naive baseline.

**Stage 3 — Systematic validation and infrastructure:**
9. Full benchmarking (F1–F4) with AUPRC/MCC/Spearman + negative controls; temporal validation; document per TRIPOD+AI and DECIDE-AI.
10. Containerize (Docker/Snakemake), precompute actionable-gene tables (H2–H3), add drift monitoring (H5), full audit logging (G4).
11. Add coarse pathway-context flags (D1), cryptic-pocket flags (D4), and India-availability metadata (E4) — all clearly secondary/advisory.

**Stage 4 — Regulatory positioning:**
12. Draft the intended-use statement with counsel to keep Tier 2 defensibly a non-clinical hypothesis-generator; if clinical deployment is intended, plan for CDSCO Class C/D (IEC 62304 QMS, clinical evaluation, ACP, post-market surveillance) and DPDP compliance.

## Caveats
- Accuracy figures are drawn from method papers and independent benchmarks; method-paper numbers (e.g., Boltz-2's "approaches FEP") are best-case and self-reported — I have paired them with independent critiques (Rowan, Fraser-lab PREreview, PoseBusters) throughout.
- The Pan et al. non-experimental per-tool R/RMSE breakdown and exact AUPRC/MCC values are reported here in aggregate; the paper's Figure 3 and supplementary tables (not fully accessible during research) hold the per-input-type numbers, and the "single best tool" and "beats a naive baseline" statements could not be verified verbatim.
- CDSCO guidance is very recent (Oct 2025 draft; mid-2026 final per secondary legal sources); confirm the current text directly with CDSCO before any regulatory action.
- "Open" licenses vary in detail (MIT vs Apache vs GPL vs academic-only vs CC-BY-NC); verify each tool's current license before commercial/clinical use — AlphaFold3 and ESM3 in particular are restricted, and MdrDB is academic-use-only.
- The fundamental ceiling stands: Tier 2 is a hypothesis-generator producing a biochemical-plausibility proxy, never a clinical predictor, and no upgrade in current technology changes that.