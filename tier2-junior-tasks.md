# Tier 2 — Junior Tasks

Six tasks, each self-contained. You shouldn't need to read `tier2-implementation-spec.md` end-to-end to do any of these — everything you need is below. Where a task needs something from the lead first, it says so explicitly under **Waiting on**.

**General rule for all six tasks:** if you hit something not covered here, or a result looks surprising (a case fails where it shouldn't, a test that should pass doesn't), **stop and flag it — don't guess and don't quietly work around it.** This whole codebase has a rule that failures get surfaced explicitly, never silently patched over or hidden; that applies to your process too, not just the code. A wrong guess costs more of the lead's time to untangle later than a question costs now.

Report back per task, not all six at once — each one unblocks something different downstream.

---

## Task 1 — Run the nine gold-standard cases

**Waiting on:** the lead's Step 7 pipeline function to exist (`tier2-implementation-spec.md` Deliverable 1). Ask before starting if you're not sure it's ready.

**What you're doing:** running nine known mutation/drug cases through the pipeline and recording exactly what comes back. You are **not** deciding pass/fail — just producing an accurate data table. The interpretation is the lead's call.

**The nine cases** (from `docs/validation-plan.md`):

| Gene / mutation | Drug | Known direction (for later reference — don't use this to "correct" a result) |
|---|---|---|
| EGFR T790M | gefitinib or erlotinib | Resistance |
| EGFR T790M | osimertinib | Sensitive |
| EGFR C797S | osimertinib | Resistance |
| ABL1 T315I | imatinib | Resistance |
| KIT D816V | imatinib | Resistance |
| KIT V560G | imatinib | More sensitive |
| BRAF V600E | vemurafenib | Sensitive |
| ALK G1202R | crizotinib | Resistance |
| ALK I1171T | crizotinib | Resistance |

**Steps:**
1. For each row, call the Step 7 pipeline function with the gene + mutation + drug.
2. Record, per case: the raw `delta_score`, the `method` used (`mCSM-lig` or `docking`), the `alphamissense` score/class, the `structure` source/id/pLDDT, and whether the call succeeded or returned a failure object (§8 message — if so, record which one verbatim).
3. Put this in a table — markdown is fine, e.g. `validation/results.md` — one row per case, one column per field above, plus the "known direction" column for reference.
4. **Do not** add a "pass/fail" column yourself, and do not adjust anything if a result looks "wrong" against the known direction — record exactly what the pipeline returned. If a case errors out entirely (timeout, crash, unhandled exception), record that as-is and flag it immediately rather than retrying silently — an unhandled exception anywhere in Steps 1–7 is itself a bug to report, not something to route around.

**Done when:** all nine rows are recorded with real pipeline output (not partial/mocked), and you've handed the table to the lead.

**Hand back to:** the lead, directly — this feeds their Step 6 cutoff decision, which is the next thing blocked on your output.

---

## Task 2 — Vina vs. mCSM-lig coverage tally

**Waiting on:** Task 1's completed results table.

**What you're doing:** simple counting from data you already have.

**Steps:**
1. From Task 1's table, count how many of the nine cases (and, if any case scored multiple candidates, how many total candidate-scores) used `method: mCSM-lig` vs. `method: docking` (Vina).
2. Compute the percentage split.
3. Write a short note (a few sentences, in `validation/results.md` right under the table, or a separate short file — your call) stating the split plainly, e.g. "X of Y candidates scored via mCSM-lig (Z%), remainder via Vina docking."
4. One sentence of context to include: mCSM-lig has a published accuracy figure for this exact task (correlation ρ up to 0.67, Pires/Blundell/Ascher 2016); Vina/smina does not have an equivalent published figure here. You don't need to interpret what this means for the product — just state the split and that fact; the lead uses it when writing the confidence framing.

**Done when:** the tally and note exist and the numbers are traceable back to Task 1's table (i.e., someone could recount them from the table and get the same answer).

**Hand back to:** the lead.

---

## Task 3 — Demo-case caching

**Waiting on:** the lead telling you which of the nine gold-standard cases is the locked demo case. Don't guess this — ask if it's not obvious.

**What you're doing:** running the pipeline once for the locked case and saving its full output so the live demo never has to make a real API call.

**Steps:**
1. Run the Step 7 pipeline for the locked demo case exactly as in Task 1.
2. Serialize the **complete** output — not just the final result object, but everything the pipeline needed along the way if that's feasible (structure file contents, intermediate API responses) so a re-run from cache genuinely needs zero external network calls. If full intermediate caching isn't practical, at minimum cache the final assembled result object and say so explicitly in your notes — don't silently cache a partial thing and imply it's complete.
3. Save it in whatever format the frontend/orchestration layer expects to consume (check with the lead if this isn't already decided — it may not be built yet, in which case save it as clean JSON and note that the consuming format is still open).
4. Confirm — or if the consuming UI doesn't exist yet, note for later — that wherever this cached result gets shown, it's marked as cached rather than presented as live. `docs/ui-flow.md` Screen 3 requires a visible "cached result" indicator; this is a hard requirement, not a nice-to-have, per the project's rule against presenting anything as live when it isn't.
5. Write a short script (a few lines is fine) that regenerates the cache on demand, rather than a one-off manual process — so it can be refreshed if the pipeline changes before the actual demo.

**Done when:** the cache file exists, the regeneration script runs cleanly from a fresh checkout, and you've noted anywhere the "cached, not live" indicator still needs wiring up on the consuming side.

**Hand back to:** the lead.

---

## Task 4 — Verification-pass checklist

**Waiting on:** nothing — this can run against the current state of the repo at any point, but is most useful after Deliverable 1 (Step 7) and Task 5 (labeling) both exist.

**What you're doing:** a checklist pass confirming the codebase is actually in the state it's supposed to be, not assumed to be. This is mechanical — you're not exercising judgment about *whether* something is a good idea, just verifying it's true.

**Checklist:**
1. Run the full test suite two ways and record both results in full (not just "passed"):
   ```bash
   pytest
   pytest -m integration
   ```
   The second one hits live external services (UniProt, Ensembl, RCSB, AlphaFold DB, DGIdb, Open Targets, ChEMBL, PubChem, mCSM-lig, and now the Step 7/labeling additions) — expect it to be slower and occasionally flaky on a genuinely-down external service. If something fails, record exactly what and don't rerun-until-green without noting the first failure.
2. Search the codebase for bare exception handling, which this project explicitly forbids everywhere:
   ```bash
   grep -rn "except:" src/
   grep -rn "except Exception" src/
   ```
   Any hit is a violation of the project's error-handling rule (every external call needs its own specific exception type, caught at the call site) — report each hit's file/line rather than fixing it yourself, unless the lead has told you fixing-on-sight is fine for this pass.
3. For every new external-facing failure path introduced in Steps 6–7 (the labeling function's edge cases, Step 7's failure-object paths), confirm there's an actual test exercising it — not just a happy-path test. List any you find that are missing a failure-path test.
4. Take one real assembled Step 7 output (any of the nine cases from Task 1 works) and check it field-by-field against the schema in `docs/api-contracts.md`'s "Tier 2 result item" section. Confirm every field is present, correctly typed, and that `label` is one of exactly the three defined strings (`likely_reduced_binding`, `likely_retained_or_increased_binding`, `uncertain`) — never anything else.

**Done when:** you've produced a short written report (pass/fail per checklist item, with specifics for anything that failed) — not just a verbal "looks good."

**Hand back to:** the lead, as a written report — they decide what (if anything) needs fixing before merge; you don't need to fix things yourself unless asked.

---

## Task 5 — Implement `labeling.py`

**Waiting on:** the lead's Step 6 cutoff decision (`tier2-implementation-spec.md` §5 item 3). Do not start this by picking your own threshold — the number itself is the lead's call, informed by Task 1's data. You're implementing to a spec they hand you, not deriving the spec.

**What you're building:** a small, pure function. Once you have the exact cutoff value and the three label strings from the lead, this is mechanical.

**Steps:**
1. Write a function taking a `delta_score` (float) and returning one of exactly three strings: `likely_reduced_binding`, `likely_retained_or_increased_binding`, or `uncertain`. The cutoff(s) must be **named constants at module level**, not inline numbers in the function body — e.g. `REDUCED_BINDING_THRESHOLD = <value the lead gave you>`, not a bare number in an `if` statement.
2. Add a comment above the constant(s) stating: the exact value, that it's an internal heuristic (not a validated clinical threshold — this exact phrase needs to appear, it's a project-wide requirement from `docs/tier2-structural-prediction.md` §Step 6), and which validation run produced it (point at Task 1's results table).
3. Write unit tests. At minimum:
   - One case per label (a delta clearly in each range).
   - **Boundary values** — exactly at each cutoff, and one increment on either side of it. This matters more than it looks: a previous bug in a different part of this codebase (RareCure's `clamp_weights`, documented in `docs/rarecure-build-reference.md` §4) shipped with a test that was written to match the buggy output instead of the intended contract. Don't do that — write the test to assert what the function is *supposed* to guarantee, and if the implementation doesn't meet it, that's a bug to report, not a test to loosen.
   - A **post-condition check**: for a large batch of varied/random delta values, assert the function's output is always one of exactly the three defined strings, never anything else (`None`, an unexpected string, an exception).
4. If the lead's spec mentions method-aware confidence (an mCSM-lig-derived score vs. a Vina-derived score should not carry identical confidence framing) and gives you enough detail to implement it, do so; if they haven't specified how yet, implement the basic three-label function first and flag that the confidence-framing piece is still open rather than guessing at it.

**Done when:** the function exists as a named-constant-driven pure function, tests cover the three labels + boundaries + the post-condition check, and the code comment documents the threshold's provenance.

**Hand back to:** the lead for review before merge.

---

## Task 6 — Background research: mCSM-lig / docking delta-score thresholds in published work (optional, informs the lead's Task 5-blocking decision)

**Waiting on:** nothing — can run in parallel with anything else, and is most useful *before* the lead needs to pick the Step 6 cutoff, so do this early if you have spare time.

**What you're doing:** a bounded literature search. You are gathering input for the lead's decision, not making the decision.

**Steps:**
1. Look at how mCSM-lig's own validation work (Pires, Blundell & Ascher, *Scientific Reports* 2016 — the paper `docs/validation-plan.md` already cites for the ρ = 0.67 correlation figure) characterizes what counts as a "meaningfully destabilizing" vs. "negligible" predicted affinity change in `log(affinity fold change)` units — the same units the pipeline's `delta_score` uses for mCSM-lig-sourced results.
2. Check whether the **Platinum database** (the training set mCSM-lig itself was built on, also referenced in `docs/data-sources.md`'s validation-datasets table) has any documented convention for a "significant change" cutoff.
3. For the Vina/docking side, look for any published convention on what magnitude of docking-score delta (kcal/mol) is treated as meaningful for a single point mutation's effect on binding — this is a different scale/units than mCSM-lig's output, so don't try to force one shared number across both methods; note that they may need separate cutoffs.
4. Write up findings as a short doc (half a page is plenty) — what range(s) other published work treats as meaningful, with citations, and an explicit note on the units mismatch between the two methods. **Do not recommend a specific final cutoff for this project** — that's the lead's decision informed by the actual gold-standard validation run (Task 1), not by literature precedent alone, since this project's methodology (mCSM-lig on drug-response mutations specifically) isn't identical to what these papers validated.

**Done when:** the short write-up exists with citations, framed as background input rather than a recommendation.

**Hand back to:** the lead, before they finalize the Task 5 spec.
