#!/usr/bin/env python3
"""Run the nine (plus one) gold-standard mutation/drug cases through the Tier 2 pipeline."""

import json
import sys
import time
from pathlib import Path

# Ensure src is on the path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from secondlook.pipeline import run_tier2_pipeline

# Test cases: (gene_symbol, mutation, drug_name, known_direction)
TEST_CASES = [
    ("EGFR", "T790M", "gefitinib", "Resistance"),
    ("EGFR", "T790M", "erlotinib", "Resistance"),
    ("EGFR", "T790M", "osimertinib", "Sensitive"),
    ("EGFR", "C797S", "osimertinib", "Resistance"),
    ("ABL1", "T315I", "imatinib", "Resistance"),
    ("KIT", "D816V", "imatinib", "Resistance"),
    ("KIT", "V560G", "imatinib", "More sensitive"),
    ("BRAF", "V600E", "vemurafenib", "Sensitive"),
    ("ALK", "G1202R", "crizotinib", "Resistance"),
    ("ALK", "I1171T", "crizotinib", "Resistance"),
]


def find_drug_item(items, drug_name):
    """Find a pipeline result item matching a specific drug name (case-insensitive)."""
    target = drug_name.lower()
    for item in items:
        if (item.get("drug") or "").lower() == target:
            return item
    return None


def run_all():
    results = []

    for gene, mutation, drug, known_dir in TEST_CASES:
        label = f"{gene} {mutation} / {drug}"
        print(f"Running: {label} ...", flush=True)
        start = time.time()

        try:
            result = run_tier2_pipeline(
                gene,
                mutation,
                min_delay_seconds=3.0,
                max_candidates=20,
            )
            elapsed = time.time() - start
            print(f"  done in {elapsed:.1f}s", flush=True)

            if result.failure is not None:
                results.append({
                    "gene_mutation": f"{gene} {mutation}",
                    "drug": drug,
                    "delta_score": None,
                    "method": None,
                    "am_score": None,
                    "am_class": None,
                    "structure_source": None,
                    "structure_id": None,
                    "plddt": None,
                    "status": "pipeline_failure",
                    "failure_message": result.failure.get("reason", "unknown"),
                    "known_direction": known_dir,
                    "all_drugs_found": [],
                })
                print(f"  PIPELINE FAILURE: {result.failure.get('reason')}")
                continue

            # Find the specific drug in items
            item = find_drug_item(result.items, drug)
            all_drug_names = [it.get("drug", "?") for it in result.items]

            if item is None:
                results.append({
                    "gene_mutation": f"{gene} {mutation}",
                    "drug": drug,
                    "delta_score": None,
                    "method": None,
                    "am_score": None,
                    "am_class": None,
                    "structure_source": None,
                    "structure_id": None,
                    "plddt": None,
                    "status": "drug_not_in_candidates",
                    "failure_message": f"Drug '{drug}' not found among candidates: {all_drug_names}",
                    "known_direction": known_dir,
                    "all_drugs_found": all_drug_names,
                })
                print(f"  Drug '{drug}' NOT in candidates: {all_drug_names}")
                continue

            am = item.get("alphamissense") or {}
            struct = item.get("structure") or {}
            results.append({
                "gene_mutation": f"{gene} {mutation}",
                "drug": drug,
                "delta_score": item.get("delta_score"),
                "method": item.get("method"),
                "am_score": am.get("score"),
                "am_class": am.get("class"),
                "structure_source": struct.get("source"),
                "structure_id": struct.get("id"),
                "plddt": struct.get("plddt_at_residue"),
                "status": "scored" if item.get("delta_score") is not None else "binding_unavailable",
                "failure_message": item.get("binding_note"),
                "known_direction": known_dir,
                "all_drugs_found": all_drug_names,
            })
            ds = item.get("delta_score")
            meth = item.get("method")
            print(f"  delta_score={ds}, method={meth}, binding_note={item.get('binding_note')}")

        except Exception as exc:
            elapsed = time.time() - start
            results.append({
                "gene_mutation": f"{gene} {mutation}",
                "drug": drug,
                "delta_score": None,
                "method": None,
                "am_score": None,
                "am_class": None,
                "structure_source": None,
                "structure_id": None,
                "plddt": None,
                "status": "unhandled_exception",
                "failure_message": f"{type(exc).__name__}: {exc}",
                "known_direction": known_dir,
                "all_drugs_found": [],
            })
            print(f"  UNHANDLED EXCEPTION: {type(exc).__name__}: {exc}")

    # Write JSON results
    out_path = Path(__file__).resolve().parent / "raw_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults written to {out_path}")
    return results


if __name__ == "__main__":
    run_all()
