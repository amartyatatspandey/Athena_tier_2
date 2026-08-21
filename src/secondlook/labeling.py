"""Step 6 binding-score labeling.

PLACEHOLDER: this always returns "uncertain". The numeric threshold has not been
decided yet — it is chosen empirically from a validation run that needs the
pipeline to exist first. See:

  - tier2-implementation-spec.md §5 item 3 (lead's cutoff decision)
  - tier2-junior-tasks.md Task 5 (where the real implementation will land)

Do not invent a threshold here. THRESHOLD_DECIDED stays False until Task 5.
"""

from __future__ import annotations

# Placeholder pending the lead's validation-run decision.
# See tier2-implementation-spec.md §5 item 3 and tier2-junior-tasks.md Task 5.
THRESHOLD_DECIDED = False


def label_binding_score(delta_score: float | None) -> str:
    """Always returns "uncertain".

    Placeholder pending the lead's validation-run decision. The real cutoff
    lands in junior Task 5 (tier2-junior-tasks.md) after the lead picks a
    threshold from the gold-standard run (tier2-implementation-spec.md §5 item 3).
    """
    return "uncertain"
