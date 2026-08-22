#!/usr/bin/env python3
"""Run the ten gold-standard cases through the Tier 2 pipeline, scoring only the target drug."""

import json
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from secondlook.alphamissense import lookup_alphamissense
from secondlook.binding import score_binding, StructureHetResolver
from secondlook.candidates import generate_candidates
from secondlook.mutation_validation import validate_mutation
from secondlook.pipeline import _result_item, label_binding_score
from secondlook.structure import source_structure
from secondlook.uniprot import UniProtLookupError

# (gene, mutation, drug, known_direction)
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


def find_candidate_by_name(candidates, drug_name):
    """Find a DrugCandidate with matching name (case-insensitive)."""
    target = drug_name.upper().strip()
    for c in candidates:
        if c.name.upper().strip() == target:
            return c
    return None


def run_case(gene, mutation, drug, known_dir):
    """Run one test case and return a result dict."""
    label = f"{gene} {mutation} / {drug}"
    print(f"\n{'='*60}", flush=True)
    print(f"Running: {label}", flush=True)
    print(f"{'='*60}", flush=True)
    t0 = time.time()

    result = {
        "gene_mutation": f"{gene} {mutation}",
        "drug": drug,
        "known_direction": known_dir,
        "delta_score": None,
        "method": None,
        "am_score": None,
        "am_class": None,
        "structure_source": None,
        "structure_id": None,
        "plddt": None,
        "status": None,
        "failure_message": None,
    }

    # Step 1: Validate mutation
    try:
        validation = validate_mutation(gene, mutation)
        print(f"  [Step 1] Mutation: status={validation.status}, gene={validation.gene}, pos={validation.position}")
    except UniProtLookupError as e:
        result["status"] = "pipeline_failure"
        result["failure_message"] = f"UniProtLookupError: {e}"
        print(f"  FAIL (UniProt): {e}")
        return result

    if validation.status != "valid":
        result["status"] = "pipeline_failure"
        result["failure_message"] = validation.error_message or "validation failed"
        print(f"  FAIL (validation): {validation.error_message}")
        return result

    # Step 2: AlphaMissense
    try:
        alphamissense = lookup_alphamissense(validation)
        result["am_score"] = alphamissense.am_pathogenicity
        result["am_class"] = alphamissense.am_class
        print(f"  [Step 2] AlphaMissense: score={alphamissense.am_pathogenicity}, class={alphamissense.am_class}")
    except Exception as e:
        print(f"  WARN (AlphaMissense): {type(e).__name__}: {e}")

    # Step 3: Structure
    try:
        structure = source_structure(validation)
        result["structure_source"] = structure.source
        result["structure_id"] = structure.id
        result["plddt"] = structure.plddt_at_residue
        print(f"  [Step 3] Structure: source={structure.source}, id={structure.id}, ligand_bound={structure.ligand_bound}, pLDDT={structure.plddt_at_residue}")
    except Exception as e:
        print(f"  FAIL (structure): {type(e).__name__}: {e}")
        result["status"] = "structure_failure"
        result["failure_message"] = f"{type(e).__name__}: {e}"
        return result

    # Step 4: Generate candidates (to find our target drug)
    try:
        candidates_result = generate_candidates(validation, max_candidates=50)
        print(f"  [Step 4] Candidates: status={candidates_result.status}, count={len(candidates_result.candidates)}")
    except Exception as e:
        print(f"  FAIL (candidates): {type(e).__name__}: {e}")
        result["status"] = "candidates_failure"
        result["failure_message"] = f"{type(e).__name__}: {e}"
        return result

    if candidates_result.status == "none":
        result["status"] = "no_candidates"
        result["failure_message"] = candidates_result.error_message or "no candidates found"
        print(f"  FAIL (no candidates): {candidates_result.error_message}")
        return result

    # Find our target drug
    target = find_candidate_by_name(candidates_result.candidates, drug)
    if target is None:
        all_names = [c.name for c in candidates_result.candidates]
        result["status"] = "drug_not_in_candidates"
        result["failure_message"] = f"Drug '{drug}' not found among {len(all_names)} candidates: {all_names}"
        print(f"  FAIL (drug not found): '{drug}' not in candidates")
        return result

    print(f"  [Step 4] Target drug found: {target.name} (source={target.source}, smiles={bool(target.smiles)})")

    # Step 5: Score binding
    try:
        binding = score_binding(
            validation, structure, target,
            min_delay_seconds=1.0,
        )
        elapsed = time.time() - t0
        result["delta_score"] = binding.delta_score
        result["method"] = binding.method
        result["status"] = "scored" if binding.delta_score is not None else "binding_unavailable"
        result["failure_message"] = binding.error_message
        print(f"  [Step 5] Binding: method={binding.method}, delta={binding.delta_score}, error={binding.error_message}")
    except Exception as e:
        elapsed = time.time() - t0
        tb = traceback.format_exc()
        result["status"] = "unhandled_exception"
        result["failure_message"] = f"{type(e).__name__}: {e}"
        print(f"  UNHANDLED EXCEPTION: {type(e).__name__}: {e}")
        print(tb)

    elapsed = time.time() - t0
    print(f"  Total time: {elapsed:.1f}s")
    return result


def main():
    all_results = []
    total_start = time.time()

    for gene, mutation, drug, known_dir in TEST_CASES:
        r = run_case(gene, mutation, drug, known_dir)
        all_results.append(r)

    total_elapsed = time.time() - total_start
    print(f"\n\nTotal time for all cases: {total_elapsed:.1f}s")

    # Write JSON results
    out_path = Path(__file__).resolve().parent / "raw_results.json"
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"Results written to {out_path}")

    # Print summary
    print("\n--- SUMMARY ---")
    for r in all_results:
        ds = r["delta_score"]
        meth = r["method"]
        status = r["status"]
        print(f"  {r['gene_mutation']} / {r['drug']}: status={status}, method={meth}, delta={ds}")


if __name__ == "__main__":
    main()
