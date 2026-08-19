# Data Sources Reference

Confirm every row below actually works in your environment before building around it — several hosted APIs in this space have unstable uptime. Test with a single manual call first.

## Tier 1 sources

| Source | Purpose | Access | Auth | Notes |
|---|---|---|---|---|
| CIViC | Variant-disease-drug evidence graph backbone | Open, no license | None | Structured, community-curated; check current API/download format |
| ClinicalTrials.gov API v2 | Global trial search | REST, free | None | |
| CTRI (ctri.nic.in) | India trial registry | Check current access mode (API vs. scrape) | Verify | May require different integration approach than CTgov — confirm early |
| PubMed E-utilities | Literature abstracts | REST, free | Optional API key (higher rate limit) | |
| RareCure modules | Drug-gene matching, trial screening (reuse where runnable) | MIT license, open source | None | **Verify exact repo link from paper's data-availability section before assuming any module runs as-is** |

## Tier 2 sources

| Source | Purpose | Access | Auth | Rate limits / notes |
|---|---|---|---|---|
| UniProt REST | Canonical protein sequences | REST, free | None | `/uniprotkb/{accession}.fasta` |
| hgvs (PyPI) / Mutalyzer | Mutation notation parsing | pip / free API | None | Mutalyzer: mutalyzer.nl/api |
| AlphaMissense precomputed scores | Pathogenicity pre-filter | Download (dbNSFP v4.7c) or hosted lookup | None for download | **Verify license on the exact file/route used** — was CC BY-NC-SA, precomputed scores relicensed CC BY 4.0 as of March 2024; confirm at build time |
| RCSB PDB | Experimental structures | REST/download, free | None | |
| AlphaFold DB | Precomputed structure models | alphafold.ebi.ac.uk, free | None | Prefer over de novo folding |
| ESM Atlas fold API | De novo folding (fallback only) | Historically free | None | **Status unreliable — documented outages. Do not put on demo-critical path. Verify live before any dependency.** |
| ColabFold | De novo folding (last resort) | Free Colab notebook | Google account | ~5–10 min/sequence; free tier disconnects after ~90 min idle / ~12h total |
| DGIdb | Drug-gene interactions | REST, free | None | |
| Open Targets | Target-disease-drug associations | GraphQL, free | None | |
| ChEMBL | Bioactivities, mechanisms, SMILES | REST, "Open Data" | None | |
| PubChem PUG REST | Drug SMILES | REST, free | None | **≤5 req/sec, ≤400/min, ≤300s running time/min** — enforce client-side throttling |
| mCSM-lig | Mutation → protein-ligand ΔΔ affinity | Web server only (biosig.lab.uq.edu.au/mcsm_lig or structure.bioc.cam.ac.uk/mcsm_lig) | None advertised | **No documented batch/local install — confirm interactive-use terms before scripting bulk calls; respect server load** |
| AutoDock Vina / smina | Docking (stretch/fallback) | Free, open source, local install | None | Structure prep (protonation, grid box setup) is the hard part — budget real time for this |
| BindingDB | Experimental binding affinity (reference/validation) | Free, CC BY | None | |

## Explicitly avoid depending on

| Source | Why |
|---|---|
| DrugBank | Requires academic license for programmatic use — not freely usable for a hackathon build |
| OncoKB API | Requires registration + academic license for treatment-level data; browse the website for reference only, don't build a live dependency |
| DeepPurpose | Not mutation-sensitive (see `tier2-structural-prediction.md`) — do not use as the affinity-change scorer |
| Rosetta ddG | Days-long setup for non-experts — infeasible in hackathon window |

## Validation/benchmark datasets (for `validation-plan.md`)

| Dataset | Use |
|---|---|
| Platinum | Benchmark for mCSM-lig-type protein-ligand affinity change tools |
| SKEMPI | Protein-protein interaction benchmark (not directly needed here) |
| ThermoMutDB / ProTherm | Stability ΔΔG benchmarks (for optional FoldX/DDGun cross-check) |
| COSMIC resistance mutations | Reference for known resistance cases |
| CIViC | Also serves as ground truth for the known gold-standard cases in `validation-plan.md` |
