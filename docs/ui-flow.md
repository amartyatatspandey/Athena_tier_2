# UI Flow — Doctor-Facing Interface

## Screen 1 — Patient/case intake

Fields:
- Cancer type (dropdown/autocomplete, scoped to your pre-loaded CIViC subset for V1)
- Gene (text/autocomplete)
- Mutation (text field, accepts shorthand or HGVS — show format hint: "e.g. R175H or p.Arg175His")
- Optional: report upload (PDF) — **stretch feature, not V1 core**

## Screen 2 — Doctor's clinical form

Fields:
- Prior treatments tried (repeatable: drug name + outcome/reason stopped)
- Clinical question (free text — what the doctor actually wants answered)
- Manual override toggle: "Also run structural/computational analysis" (available even if Tier 1 finds strong evidence)

Submit → calls `POST /api/query` (see `api-contracts.md`).

## Screen 3 — Processing state

- Show which stage is active: "Searching documented evidence…" → (if triggered) "No strong match found — running structural analysis…" → "Preparing results."
- If any component is using a cached/pre-computed result (demo mode), do not hide this — show a small "cached result" indicator per section rather than pretending it's live, to keep the demo honest to yourselves and avoid confusion if a judge asks.

## Screen 4 — Results screen

**Section A — Documented evidence (Tier 1)**
- Visual style: solid border, green accent.
- Each card: drug/trial name, evidence level, one-line summary, source link (trial ID / CIViC entry / PubMed citation), clickable.
- If empty: explicit "No documented evidence found for this exact mutation" message — never an empty section with no explanation.

**Section B — Computational signal (Tier 2)**, shown only if it ran
- Visual style: amber/outlined border, clearly separated from Section A, persistent header badge: **"COMPUTATIONAL PLAUSIBILITY SIGNAL — NOT CLINICAL EVIDENCE."**
- Each card: drug name, predicted label (reduced/retained/uncertain), method used (mCSM-lig or docking), structure source + confidence flag, binding-site distance if available, AlphaMissense pathogenicity context.
- Fixed disclaimer text (from `tier2-structural-prediction.md` §10) shown once at the top of this section, not per-card (avoid clutter, but must be unmissable).
- Sort order: predicted effect magnitude × confidence × structural reliability, with India-availability/approval status as a tiebreaker if you've wired in CDSCO data.

**Section C — Failures/inconclusive states**, shown if any occurred
- Render the exact message from `tier2-structural-prediction.md` §8 or the Tier 1 equivalent — never omit this section if a failure object was returned.

**Synthesis summary** (top of results screen, above A/B/C)
- Short LLM-generated paragraph referencing only items shown below it, each claim traceable via the citation/source shown in the relevant card.

## Copy rules (apply everywhere in this UI)

- Use: "predicted," "computational signal," "plausibility," "for physician review."
- Never use: "diagnosis," "prescribed," "recommended treatment," "proven," "will work."
- Every Tier 2 element visually distinct from Tier 1 at a glance (color, border style, badge) — a doctor scanning quickly should never mistake one for the other.

## Demo script screen (internal, not shown to judges as a UI screen — for your team's rehearsal)

See `validation-plan.md` for the two required demo cases (one Tier 1 path, one Tier 2 path) and how to present the difference live.
