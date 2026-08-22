#!/usr/bin/env python3
"""Run a batch of gold-standard cases through the Tier 2 pipeline."""

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
from secondlook.structure import source_structure
from secondlook.uniprot import UniProtLookupError

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
    target = drug_name.upper().strip()
    for c in candidates:
        if c.name.upper().strip() == target:
            return c
    return None


def run_case(gene, mutation, drug, known_dir):
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

    try:
        validation = validate_mutation(gene, mutation)
        print(f"  [Step 1] status={validation.status}, gene={validation.gene}, pos={validation.position}")
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

    try:
        alphamissense = lookup_alphamissense(validation)
        result["am_score"] = alphamissense.am_pathogenicity
        result["am_class"] = alphamissense.am_class
        print(f"  [Step 2] AM: score={alphamissense.am_pathogenicity}, class={alphamissense.am_class}")
    except Exception as e:
        print(f"  WARN (AlphaMissense): {type(e).__name__}: {e}")

    try:
        structure = source_structure(validation)
        result["structure_source"] = structure.source
        result["structure_id"] = structure.id
        result["plddt"] = structure.plddt_at_residue
        print(f"  [Step 3] Struct: source={structure.source}, id={structure.id}, ligand={structure.ligand_bound}, pLDDT={structure.plddt_at_residue}")
    except Exception as e:
        result["status"] = "structure_failure"
        result["failure_message"] = f"{type(e).__name__}: {e}"
        print(f"  FAIL (structure): {e}")
        return result

    try:
        candidates_result = generate_candidates(validation, max_candidates=50)
        print(f"  [Step 4] Candidates: count={len(candidates_result.candidates)}")
    except Exception as e:
        result["status"] = "candidates_failure"
        result["failure_message"] = f"{type(e).__name__}: {e}"
        print(f"  FAIL (candidates): {e}")
        return result

    if candidates_result.status == "none":
        result["status"] = "no_candidates"
        result["failure_message"] = candidates_result.error_message
        print(f"  FAIL (no candidates)")
        return result

    target = find_candidate_by_name(candidates_result.candidates, drug)
    if target is None:
        all_names = [c.name for c in candidates_result.candidates]
        result["status"] = "drug_not_in_candidates"
        result["failure_message"] = f"'{drug}' not in {len(all_names)} candidates: {all_names}"
        print(f"  FAIL (drug not found)")
        return result

    print(f"  [Step 4] Found: {target.name}, smiles={bool(target.smiles)}")

    try:
        binding = score_binding(validation, structure, target, min_delay_seconds=1.0)
        result["delta_score"] = binding.delta_score
        result["method"] = binding.method
        result["status"] = "scored" if binding.delta_score is not None else "binding_unavailable"
        result["failure_message"] = binding.error_message
        print(f"  [Step 5] method={binding.method}, delta={binding.delta_score}")
        if binding.error_message:
            print(f"           note: {binding.error_message[:120]}")
    except Exception as e:
        result["status"] = "unhandled_exception"
        result["failure_message"] = f"{type(e).__name__}: {e}"
        print(f"  UNHANDLED: {type(e).__name__}: {e}")
        traceback.print_exc()

    elapsed = time.time() - t0
    print(f"  Time: {elapsed:.1f}s")
    return result


def main():
    start_idx = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    end_idx = int(sys.argv[2]) if len(sys.argv) > 2 else len(TEST_CASES)

    cases = TEST_CASES[start_idx:end_idx]
    print(f"Running cases {start_idx}..{end_idx-1} ({len(cases)} cases)")

    out_path = Path(__file__).resolve().parent / "raw_results.json"

    # Load existing results if resuming
    existing = []
    if out_path.exists():
        with open(out_path) as f:
            existing = json.load(f)
        print(f"Loaded {len(existing)} existing results")

    total_start = time.time()
    for i, (gene, mutation, drug, known_dir) in enumerate(cases):
        # Skip if already run
        key = f"{gene} {mutation}|{drug}"
        already = any(
            r["gene_mutation"] == f"{gene} {mutation}" and r["drug"] == drug
            for r in existing
        )
        if already:
            print(f"\nSkipping {key} (already in results)")
            continue

        r = run_case(gene, mutation, drug, known_dir)
        existing.append(r)

        # Save after each case
        with open(out_path, "w") as f:
            json.dump(existing, f, indent=2, default=str)

    total_elapsed = time.time() - total_start
    print(f"\nTotal: {total_elapsed:.1f}s")

    print("\n--- SUMMARY ---")
    for r in existing:
        ds = r["delta_score"]
        meth = r["method"]
        status = r["status"]
        print(f"  {r['gene_mutation']} / {r['drug']}: {status} | method={meth} | delta={ds}")


if __name__ == "__main__":
    main()
