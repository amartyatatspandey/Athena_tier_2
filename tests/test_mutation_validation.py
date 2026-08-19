from dataclasses import dataclass

from secondlook.mutation_validation import (
    OUT_OF_SCOPE_MESSAGE,
    normalize_protein_notation,
    validate_mutation,
)


@dataclass(frozen=True)
class FakeProtein:
    accession: str
    sequence: str
    gene: str | None = "TP53"
    isoform_note: str = "P04637 (canonical)"


class FakeSequenceProvider:
    def __init__(self, protein: FakeProtein) -> None:
        self.protein = protein

    def fetch(self, identifier: str) -> FakeProtein:
        return self.protein


def _sequence_with_residue(position: int, residue: str, length: int = 200) -> str:
    chars = ["A"] * length
    chars[position - 1] = residue
    return "".join(chars)


def test_normalizes_single_letter_shorthand_to_hgvs_protein():
    assert normalize_protein_notation("R175H") == "p.Arg175His"


def test_leaves_hgvs_protein_notation_unchanged():
    assert normalize_protein_notation("p.Arg175His") == "p.Arg175His"


def test_valid_missense_matches_canonical_residue_and_applies_substitution():
    provider = FakeSequenceProvider(
        FakeProtein(accession="P04637", sequence=_sequence_with_residue(175, "R"))
    )

    result = validate_mutation("P04637", "R175H", sequence_provider=provider)

    assert result.status == "valid"
    assert result.mutation_type == "missense"
    assert result.hgvs_normalized == "p.Arg175His"
    assert result.uniprot_accession == "P04637"
    assert result.position == 175
    assert result.reference_residue_expected == "R"
    assert result.reference_residue_actual == "R"
    assert result.error_message is None
    assert result.wildtype_sequence[174] == "R"
    assert result.mutant_sequence[174] == "H"
    assert result.mutant_sequence[:174] == result.wildtype_sequence[:174]
    assert result.mutant_sequence[175:] == result.wildtype_sequence[175:]


def test_reference_residue_mismatch_stops_with_exact_message():
    provider = FakeSequenceProvider(
        FakeProtein(accession="P04637", sequence=_sequence_with_residue(175, "Y"))
    )

    result = validate_mutation("P04637", "R175H", sequence_provider=provider)

    assert result.status == "reference_mismatch"
    assert result.position == 175
    assert result.reference_residue_expected == "R"
    assert result.reference_residue_actual == "Y"
    assert result.mutant_sequence is None
    assert result.error_message == (
        "Reference residue mismatch: notation claims R at position 175 "
        "but canonical sequence has Y. This usually indicates a transcript/isoform "
        "numbering difference. Cannot proceed safely."
    )


def _assert_out_of_scope(notation: str, mutation_type: str) -> None:
    provider = FakeSequenceProvider(
        FakeProtein(accession="P04637", sequence=_sequence_with_residue(175, "R"))
    )
    result = validate_mutation("P04637", notation, sequence_provider=provider)
    assert result.status == "unsupported_type"
    assert result.mutation_type == mutation_type
    assert result.error_message == OUT_OF_SCOPE_MESSAGE
    assert result.mutant_sequence is None


def test_rejects_insertion_as_out_of_scope():
    _assert_out_of_scope("p.Lys2_Met3insGln", "insertion")


def test_rejects_deletion_as_out_of_scope():
    _assert_out_of_scope("p.Arg175del", "deletion")


def test_rejects_frameshift_as_out_of_scope():
    _assert_out_of_scope("p.Arg175fs", "frameshift")


def test_rejects_fusion_as_out_of_scope():
    _assert_out_of_scope("TP53::NTRK1", "fusion")


def test_rejects_splice_variant_as_out_of_scope():
    _assert_out_of_scope("c.375+1G>A", "splice")


def test_accepts_already_normalized_hgvs_protein_missense():
    provider = FakeSequenceProvider(
        FakeProtein(accession="P04637", sequence=_sequence_with_residue(175, "R"))
    )
    result = validate_mutation("P04637", "p.Arg175His", sequence_provider=provider)
    assert result.status == "valid"
    assert result.hgvs_normalized == "p.Arg175His"
    assert result.mutant_sequence[174] == "H"
