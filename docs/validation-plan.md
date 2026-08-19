# Validation Plan

Run this before the pitch deck is finalized. Do not claim reliability you haven't checked here.

## Method
For each case below, run the mutation through Tier 2 **as if it were undocumented** (even though the clinical answer is known) and check whether the pipeline's predicted label matches the known direction. This is a rediscovery test, not proof of general accuracy — state it as such.

## Gold-standard cases (known clinical directionality)

| Gene / mutation | Drug | Known direction | Notes |
|---|---|---|---|
| EGFR T790M | gefitinib / erlotinib (1st-gen) | **Resistance** | Gatekeeper mutation |
| EGFR T790M | osimertinib (3rd-gen) | **Sensitive** | Designed for T790M |
| EGFR C797S | osimertinib | **Resistance** | Abolishes covalent binding at C797 |
| ABL1 T315I | imatinib | **Resistance** | Classic gatekeeper mutation |
| KIT D816V | imatinib | **Resistance** | Kinase-domain mutation |
| KIT V560G | imatinib | **More sensitive** (contrast case) | Juxtamembrane mutation — use alongside D816V to test the pipeline distinguishes direction, not just "mutation = resistance" |
| BRAF V600E | vemurafenib | **Sensitive** | Drug designed for this mutation |
| ALK G1202R | crizotinib | **Resistance** | Sensitive to lorlatinib instead |
| ALK I1171T | crizotinib | **Resistance** | Sensitive to ceritinib/lorlatinib instead |

## Pre-committed pass/fail criteria (set BEFORE running, do not adjust after seeing results)

- **Pass threshold:** correct resistance-vs-sensitive directionality on ≥70% of the cases above.
- **Hard requirement regardless of overall score:** no confidently-wrong call on the two "designed-for" positive-control cases (BRAF V600E/vemurafenib should show retained/increased binding; EGFR T790M/osimertinib should show retained/increased binding). If either of these fails, treat the pipeline as not ready to demo as-is.
- **If overall pass rate is below threshold:** fall back to the narrower claim — "mutation is in/near the known binding pocket" plus AlphaMissense pathogenicity flag only, dropping the binding-affinity delta/label from the UI. Document this fallback decision plainly if you have to invoke it; it's still a legitimate, honest V1.

## Demo case selection

Pick exactly two cases for the live demo:
1. **Tier 1 path** — a case with strong documented CIViC/trial evidence, to show the "fast, sourced retrieval" experience.
2. **Tier 2 path** — one gold-standard case from the table above (sarcoma-related if possible, ties back to RareCure's own validation cohort and the project's origin story) where Tier 2 correctly rediscovers the known direction.

Both must be **pre-computed and cached** (see `tech-stack-setup.md`) — do not run either live on stage for the first time.

## What to report to judges (accurate framing)

State plainly: "We validated Tier 2 by checking whether it rediscovers known resistance/sensitivity patterns in textbook cases like ABL1 T315I and EGFR T790M — it passed [X/9]. This is a sanity check on the method, not clinical validation, and every result is labeled accordingly in the UI."

## Cite real accuracy figures, not assumed ones

Quote the underlying models' own published figures rather than implying higher confidence:
- mCSM-lig: correlation up to ρ = 0.67 with experimental data (Pires, Blundell & Ascher, *Sci Rep* 2016).
- Independent benchmarking shows accuracy degrades 5–30% when using modeled (vs. experimental) structures.
- AlphaMissense: thresholds chosen for 90% precision on ClinVar (Cheng et al., *Science* 2023) — this is a pathogenicity classifier, not a drug-response predictor; don't conflate the two in your pitch.
