import httpx
import pytest

from secondlook.alphamissense import (
    classify_alphamissense,
    extract_alphamissense,
    lookup_alphamissense,
    protein_sub_to_hgvs_c,
)
from secondlook.ensembl import EnsemblError, EnsemblVepClient
from secondlook.mutation_validation import MutationValidationResult


class ReadTimeoutClient:
    """Simulates a connection-level failure (timeout before any response
    object exists) — distinct from an HTTP status error."""

    def get(self, url, params=None, headers=None, timeout=None):
        raise httpx.ReadTimeout("timed out")


def test_read_timeout_is_wrapped_not_raised_raw():
    client = EnsemblVepClient(client=ReadTimeoutClient())
    with pytest.raises(EnsemblError):
        client.lookup_hgvs("ENST00000269305", "c.524G>A")


def test_classifies_paper_thresholds_inclusive_on_upper_bound():
    assert classify_alphamissense(0.0) == "likely_benign"
    assert classify_alphamissense(0.333) == "likely_benign"
    assert classify_alphamissense(0.334) == "ambiguous"
    assert classify_alphamissense(0.564) == "ambiguous"
    assert classify_alphamissense(0.565) == "likely_pathogenic"
    assert classify_alphamissense(1.0) == "likely_pathogenic"


VEP_TP53_R175H = [
    {
        "transcript_consequences": [
            {
                "transcript_id": "ENST00000359597",
                "alphamissense": {"am_pathogenicity": 0.1, "am_class": "likely_benign"},
                "codons": "cGc/cAc",
                "cds_start": 524,
                "amino_acids": "R/H",
            },
            {
                "transcript_id": "ENST00000269305",
                "alphamissense": {"am_pathogenicity": 0.9857, "am_class": "likely_pathogenic"},
                "codons": "cGc/cAc",
                "cds_start": 524,
                "amino_acids": "R/H",
            },
        ]
    }
]


def test_extracts_alphamissense_from_canonical_transcript_not_other_isoforms():
    extracted = extract_alphamissense(VEP_TP53_R175H, transcript_id="ENST00000269305")
    assert extracted is not None
    assert extracted["am_pathogenicity"] == 0.9857
    assert extracted["am_class"] == "likely_pathogenic"
    assert extracted["transcript_id"] == "ENST00000269305"
    assert extracted["cds_start"] == 524
    assert extracted["codons"] == "cGc/cAc"


def test_extract_returns_none_when_canonical_transcript_has_no_alphamissense():
    payload = [
        {
            "transcript_consequences": [
                {
                    "transcript_id": "ENST00000269305",
                    "codons": "cGc/cAc",
                    "cds_start": 524,
                    "amino_acids": "R/H",
                }
            ]
        }
    ]
    assert extract_alphamissense(payload, transcript_id="ENST00000269305") is None


def test_converts_protein_substitution_to_coding_hgvs():
    cds = "ATG" + "AAA" * 173 + "CGC" + "AAA" * 20
    assert protein_sub_to_hgvs_c(cds, position=175, ref_aa="R", alt_aa="H") == "c.524G>A"


def _valid_tp53() -> MutationValidationResult:
    return MutationValidationResult(
        status="valid",
        gene="TP53",
        hgvs_normalized="p.Arg175His",
        uniprot_accession="P04637",
        isoform_note="P04637 (canonical)",
        position=175,
        reference_residue_expected="R",
        reference_residue_actual="R",
        mutation_type="missense",
        error_message=None,
        wildtype_sequence="R" * 393,
        mutant_sequence="R" * 174 + "H" + "R" * 218,
    )


class FakeTranscriptResolver:
    def resolve_mane_select(self, gene_symbol: str) -> str:
        assert gene_symbol == "TP53"
        return "ENST00000269305"

    def fetch_cds(self, transcript_id: str) -> str:
        assert transcript_id == "ENST00000269305"
        return "ATG" + "AAA" * 173 + "CGC" + "AAA" * 20


class FakeVepClient:
    def __init__(self, payload: list[dict]) -> None:
        self.payload = payload
        self.calls: list[tuple[str, str]] = []

    def lookup_hgvs(self, transcript_id: str, hgvs: str) -> list[dict]:
        self.calls.append((transcript_id, hgvs))
        return self.payload


def test_lookup_uses_mane_transcript_and_paper_class_not_ensembl_class():
    vep = FakeVepClient(VEP_TP53_R175H)
    result = lookup_alphamissense(
        _valid_tp53(),
        transcript_resolver=FakeTranscriptResolver(),
        vep_client=vep,
    )
    assert vep.calls == [("ENST00000269305", "c.524G>A")]
    assert result.status == "scored"
    assert result.am_pathogenicity == 0.9857
    assert result.am_class == "likely_pathogenic"
    assert result.transcript_id == "ENST00000269305"
    assert result.hgvs_c == "c.524G>A"
    assert result.note is None


def test_missing_alphamissense_records_unavailable_and_does_not_block():
    payload = [
        {
            "transcript_consequences": [
                {
                    "transcript_id": "ENST00000269305",
                    "codons": "cGc/cAc",
                    "cds_start": 524,
                    "amino_acids": "R/H",
                }
            ]
        }
    ]
    result = lookup_alphamissense(
        _valid_tp53(),
        transcript_resolver=FakeTranscriptResolver(),
        vep_client=FakeVepClient(payload),
    )
    assert result.status == "unavailable"
    assert result.am_pathogenicity is None
    assert result.am_class is None
    assert result.note == "no AlphaMissense score available"
    assert result.transcript_id == "ENST00000269305"


class FailingTranscriptResolver:
    def resolve_mane_select(self, gene_symbol: str) -> str:
        raise EnsemblError(f"No canonical/MANE Select transcript found for gene {gene_symbol!r}")

    def fetch_cds(self, transcript_id: str) -> str:  # pragma: no cover - not reached
        raise AssertionError("fetch_cds should not be called when transcript resolution fails")


def test_unresolvable_gene_returns_unavailable_not_an_exception():
    result = lookup_alphamissense(
        _valid_tp53(),
        transcript_resolver=FailingTranscriptResolver(),
        vep_client=FakeVepClient(VEP_TP53_R175H),
    )
    assert result.status == "unavailable"
    assert result.transcript_id is None
    assert result.note == "no AlphaMissense score available"


class FailingVepClient:
    def lookup_hgvs(self, transcript_id: str, hgvs: str) -> list[dict]:
        raise EnsemblError(f"Ensembl request failed for {transcript_id}:{hgvs}")


def test_vep_request_failure_returns_unavailable_not_an_exception():
    result = lookup_alphamissense(
        _valid_tp53(),
        transcript_resolver=FakeTranscriptResolver(),
        vep_client=FailingVepClient(),
    )
    assert result.status == "unavailable"
    assert result.transcript_id == "ENST00000269305"
    assert result.hgvs_c == "c.524G>A"
    assert result.note == "no AlphaMissense score available"


@pytest.mark.integration
def test_live_vep_tp53_r175h_alphamissense():
    from secondlook.mutation_validation import validate_mutation

    validation = validate_mutation("P04637", "R175H")
    result = lookup_alphamissense(validation)
    assert result.status == "scored"
    assert result.transcript_id is not None
    assert result.transcript_id.startswith("ENST00000269305")
    assert result.hgvs_c == "c.524G>A"
    assert result.am_pathogenicity == pytest.approx(0.9857, abs=1e-4)
    assert result.am_class == "likely_pathogenic"
    assert result.note is None
