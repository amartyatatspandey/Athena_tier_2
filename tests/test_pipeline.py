import pytest

from secondlook.binding import BindingScore, McsmLigResult
from secondlook.candidates import ZERO_CANDIDATES_MESSAGE
from secondlook.mutation_validation import OUT_OF_SCOPE_MESSAGE, reference_mismatch_message
from secondlook.uniprot import UniProtLookupError

VEP_TP53_R175H = [
    {
        "transcript_consequences": [
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

MINI_PDB = """\
ATOM      1  N   ARG A 175      0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  ARG A 175      1.000   0.000   0.000  1.00 96.62           C
ATOM      3  C   ARG A 175      2.000   0.000   0.000  1.00 11.00           C
"""

TIER2_RESULT_KEYS = {
    "type",
    "mutation_validated",
    "alphamissense",
    "structure",
    "drug",
    "smiles_source",
    "method",
    "delta_score",
    "label",
    "binding_site_distance_angstrom",
    "disclaimer",
    "binding_note",
}

STRUCTURAL_SIGNAL_KEYS = {
    "alphamissense_score",
    "alphamissense_class",
    "structure_source",
    "structure_id",
    "plddt_at_residue",
    "reliability_flag",
    "method",
    "binding_site_distance_angstrom",
    "computed_at",
    "pipeline_version",
}


class FakeProtein:
    def __init__(self, *, residue_at_175: str = "R") -> None:
        chars = ["A"] * 393
        chars[174] = residue_at_175
        self.accession = "P04637"
        self.sequence = "".join(chars)
        self.gene = "TP53"
        self.isoform_note = "P04637 (canonical)"


class FakeSequenceProvider:
    def __init__(self, protein: FakeProtein | None = None, error: Exception | None = None) -> None:
        self.protein = protein or FakeProtein()
        self.error = error
        self.calls: list[str] = []

    def fetch(self, identifier: str) -> FakeProtein:
        self.calls.append(identifier)
        if self.error:
            raise self.error
        return self.protein


class FakeTranscriptResolver:
    def __init__(self) -> None:
        self.calls: list[str] = []

    def resolve_mane_select(self, gene_symbol: str) -> str:
        self.calls.append(gene_symbol)
        return "ENST00000269305"

    def fetch_cds(self, transcript_id: str) -> str:
        return "ATG" + "AAA" * 173 + "CGC" + "AAA" * 20


class FakeVepClient:
    def __init__(self, payload: list[dict] | None = None) -> None:
        self.payload = payload or VEP_TP53_R175H
        self.calls: list[tuple[str, str]] = []

    def lookup_hgvs(self, transcript_id: str, hgvs: str) -> list[dict]:
        self.calls.append((transcript_id, hgvs))
        return self.payload


class FakePdb:
    def __init__(self, hit=None) -> None:
        self.hit = hit
        self.calls: list[str] = []

    def search_by_uniprot(self, accession: str, preferred_ligands: tuple[str, ...] = ()) -> dict | None:
        self.calls.append(accession)
        return self.hit


class FakeAlphaFold:
    def __init__(self) -> None:
        self.calls: list[str] = []

    def fetch_models(self, accession: str) -> list[dict]:
        self.calls.append(accession)
        return []


class FakeEsm:
    def __init__(self) -> None:
        self.sequences: list[str] = []

    def fold_sequence(self, sequence: str) -> str:
        self.sequences.append(sequence)
        raise AssertionError("ESM Atlas must never be on this pipeline's live path")


class FakeDgidb:
    def __init__(self, rows=None) -> None:
        self.rows = list(rows or [])
        self.calls: list[str] = []

    def fetch_drugs(self, gene_symbol: str) -> list[dict]:
        self.calls.append(gene_symbol)
        return list(self.rows)


class FakeOpenTargets:
    def __init__(self) -> None:
        self.calls: list[str] = []

    def fetch_drugs(self, gene_symbol: str) -> list[dict]:
        self.calls.append(gene_symbol)
        return []


class FakeChembl:
    def __init__(self) -> None:
        self.calls: list[str] = []

    def fetch_drugs(self, gene_symbol: str) -> list[dict]:
        self.calls.append(gene_symbol)
        return []


class FakePubChem:
    def __init__(self, smiles_by_name: dict[str, str] | None = None) -> None:
        self.smiles_by_name = smiles_by_name or {"OXALIPLATIN": "C1CCC1", "CARBOPLATIN": "CC"}
        self.calls: list[str] = []

    def fetch_smiles(self, drug_name: str) -> str:
        self.calls.append(drug_name)
        return self.smiles_by_name[drug_name]


class FakeMcsm:
    def __init__(self, result: McsmLigResult | None = None) -> None:
        self.result = result or _mcsm_result()
        self.calls: list[dict] = []

    def submit(self, **kwargs) -> McsmLigResult:
        self.calls.append(kwargs)
        return self.result


class FakeVina:
    def __init__(self, result: BindingScore | None = None, error: Exception | None = None) -> None:
        self.result = result
        self.error = error
        self.calls: list[dict] = []

    def score(self, **kwargs) -> BindingScore:
        self.calls.append(kwargs)
        if self.error:
            raise self.error
        assert self.result is not None
        return self.result


class FakeHet:
    def __init__(self, code: str | None = "065") -> None:
        self.code = code
        self.calls: list[str] = []

    def resolve(self, drug_name: str, structure) -> str | None:
        self.calls.append(drug_name)
        return self.code


def _mcsm_result() -> McsmLigResult:
    return McsmLigResult(
        affinity_change=-2.056,
        affinity_class="Destabilizing",
        wild_type="R",
        position=175,
        mutant_type="H",
        chain="A",
        ligand_id="065",
        distance_angstrom=2.814,
        duet_stability_kcal=-0.087,
    )


def _pdb_hit() -> dict:
    return {
        "pdb_id": "9C5S",
        "ligand_bound": True,
        "pdb_text": MINI_PDB,
    }


def _deps(**overrides):
    deps = {
        "sequence_provider": FakeSequenceProvider(),
        "transcript_resolver": FakeTranscriptResolver(),
        "vep_client": FakeVepClient(),
        "pdb_client": FakePdb(hit=_pdb_hit()),
        "alphafold_client": FakeAlphaFold(),
        "esm_client": FakeEsm(),
        "dgidb_client": FakeDgidb(
            [
                {"name": "OXALIPLATIN", "approved": True, "score": 0.5},
                {"name": "CARBOPLATIN", "approved": True, "score": 0.4},
            ]
        ),
        "opentargets_client": FakeOpenTargets(),
        "chembl_client": FakeChembl(),
        "pubchem_client": FakePubChem(),
        "mcsm_client": FakeMcsm(),
        "vina_client": FakeVina(error=RuntimeError("vina should not be needed on happy path")),
        "het_resolver": FakeHet(),
        "sleeper": lambda _: None,
        "min_delay_seconds": 0.0,
    }
    deps.update(overrides)
    return deps


def test_happy_path_assembles_one_item_per_candidate_with_disclaimer():
    from secondlook.pipeline import TIER2_DISCLAIMER, run_tier2_pipeline

    deps = _deps()
    result = run_tier2_pipeline("P04637", "R175H", **deps)

    assert result.failure is None
    assert len(result.items) == 2
    drugs = [item["drug"] for item in result.items]
    assert drugs == ["OXALIPLATIN", "CARBOPLATIN"]
    for item in result.items:
        assert set(item) >= TIER2_RESULT_KEYS
        assert item["type"] == "computational_signal"
        assert item["disclaimer"] is TIER2_DISCLAIMER or item["disclaimer"] == TIER2_DISCLAIMER
        assert item["alphamissense"] == {"score": 0.9857, "class": "likely_pathogenic"}
        assert item["structure"] == {
            "source": "PDB",
            "id": "9C5S",
            "plddt_at_residue": None,
            "reliability_flag": "high",
        }
        assert item["smiles_source"] == "PubChem"
        assert item["method"] == "mCSM-lig"
        assert item["delta_score"] == pytest.approx(-2.056)
        assert item["label"] == "uncertain"
        assert item["binding_site_distance_angstrom"] == pytest.approx(2.814)
        assert item["binding_note"] is None
        validated = item["mutation_validated"]
        assert validated["status"] == "valid"
        assert validated["gene"] == "TP53"
        assert validated["position"] == 175
        assert validated["hgvs_normalized"] == "p.Arg175His"
        assert "wildtype_sequence" in validated
        assert "uniprot_accession" in validated
    assert deps["esm_client"].sequences == []
    assert deps["mcsm_client"].calls  # binding ran


def test_reference_mismatch_short_circuits_with_exact_message_and_calls_nothing_else():
    from secondlook.pipeline import run_tier2_pipeline

    deps = _deps(sequence_provider=FakeSequenceProvider(FakeProtein(residue_at_175="Y")))
    result = run_tier2_pipeline("P04637", "R175H", **deps)

    assert result.items == ()
    assert result.failure is not None
    assert result.failure["type"] == "failure"
    assert result.failure["tier"] == "2"
    assert result.failure["reason"] == reference_mismatch_message("R", 175, "Y")
    assert result.failure["retryable"] is False
    assert deps["transcript_resolver"].calls == []
    assert deps["vep_client"].calls == []
    assert deps["pdb_client"].calls == []
    assert deps["dgidb_client"].calls == []
    assert deps["mcsm_client"].calls == []


def test_unsupported_mutation_type_short_circuits_the_same_way():
    from secondlook.pipeline import run_tier2_pipeline

    deps = _deps()
    result = run_tier2_pipeline("P04637", "p.Arg175del", **deps)

    assert result.items == ()
    assert result.failure is not None
    assert result.failure["type"] == "failure"
    assert result.failure["reason"] == OUT_OF_SCOPE_MESSAGE
    assert result.failure["retryable"] is False
    assert deps["sequence_provider"].calls == []
    assert deps["transcript_resolver"].calls == []
    assert deps["vep_client"].calls == []
    assert deps["pdb_client"].calls == []
    assert deps["dgidb_client"].calls == []
    assert deps["mcsm_client"].calls == []


def test_zero_candidates_short_circuits_without_calling_score_binding():
    from secondlook.pipeline import run_tier2_pipeline

    deps = _deps(dgidb_client=FakeDgidb([]))
    result = run_tier2_pipeline("P04637", "R175H", **deps)

    assert result.items == ()
    assert result.failure is not None
    assert result.failure["reason"] == ZERO_CANDIDATES_MESSAGE
    assert result.failure["retryable"] is False
    assert deps["vep_client"].calls  # AlphaMissense still ran
    assert deps["pdb_client"].calls  # structure still ran
    assert deps["dgidb_client"].calls == ["TP53"]
    assert deps["mcsm_client"].calls == []
    assert deps["het_resolver"].calls == []


def test_structure_unavailable_does_not_short_circuit():
    from secondlook.pipeline import run_tier2_pipeline

    deps = _deps(pdb_client=FakePdb(hit=None))
    result = run_tier2_pipeline("P04637", "R175H", **deps)

    assert result.failure is None
    assert len(result.items) == 2
    assert deps["vep_client"].calls
    assert deps["dgidb_client"].calls == ["TP53"]
    assert deps["esm_client"].sequences == []
    for item in result.items:
        assert item["alphamissense"]["score"] == 0.9857
        assert item["structure"]["source"] is None
        assert item["structure"]["id"] is None
        assert item["drug"] in {"OXALIPLATIN", "CARBOPLATIN"}
        assert item["method"] is None
        assert item["delta_score"] is None
        assert item["label"] is None
        assert item["binding_note"]


def test_unavailable_binding_score_still_appears_with_note():
    from secondlook.pipeline import run_tier2_pipeline
    from secondlook.vina_dock import VinaError

    deps = _deps(
        het_resolver=FakeHet(code=None),
        vina_client=FakeVina(error=VinaError("docking failed")),
    )
    result = run_tier2_pipeline("P04637", "R175H", **deps)

    assert result.failure is None
    assert len(result.items) == 2
    for item in result.items:
        assert item["drug"] in {"OXALIPLATIN", "CARBOPLATIN"}
        assert item["method"] is None
        assert item["delta_score"] is None
        assert item["label"] is None
        assert item["binding_note"]


def test_labeling_fn_injection_point_can_be_overridden():
    from secondlook.pipeline import run_tier2_pipeline

    result = run_tier2_pipeline(
        "P04637",
        "R175H",
        labeling_fn=lambda delta: "likely_reduced_binding",
        **_deps(),
    )
    assert result.failure is None
    assert {item["label"] for item in result.items} == {"likely_reduced_binding"}


def test_to_structural_signal_shape_from_sample_item():
    from secondlook.pipeline import PIPELINE_VERSION, to_structural_signal

    item = {
        "type": "computational_signal",
        "alphamissense": {"score": 0.9857, "class": "likely_pathogenic"},
        "structure": {
            "source": "PDB",
            "id": "9C5S",
            "plddt_at_residue": None,
            "reliability_flag": "high",
        },
        "method": "mCSM-lig",
        "binding_site_distance_angstrom": 2.814,
    }
    signal = to_structural_signal(item)
    assert set(signal) == STRUCTURAL_SIGNAL_KEYS
    assert signal["alphamissense_score"] == 0.9857
    assert signal["alphamissense_class"] == "likely_pathogenic"
    assert signal["structure_source"] == "PDB"
    assert signal["structure_id"] == "9C5S"
    assert signal["plddt_at_residue"] is None
    assert signal["reliability_flag"] == "high"
    assert signal["method"] == "mCSM-lig"
    assert signal["binding_site_distance_angstrom"] == 2.814
    assert signal["pipeline_version"] == PIPELINE_VERSION
    assert isinstance(signal["computed_at"], str)
    assert "T" in signal["computed_at"]


def test_unknown_gene_produces_failure_shaped_result_not_a_crash():
    from secondlook.pipeline import run_tier2_pipeline

    deps = _deps(
        sequence_provider=FakeSequenceProvider(
            error=UniProtLookupError(
                "No reviewed human UniProt entry found for gene symbol 'NOTAREALGENE123'"
            )
        )
    )
    result = run_tier2_pipeline("NOTAREALGENE123", "R175H", **deps)
    assert result.items == ()
    assert result.failure is not None
    assert result.failure["type"] == "failure"
    assert result.failure["tier"] == "2"
    assert isinstance(result.failure["reason"], str)
    assert result.failure["reason"]
    assert deps["vep_client"].calls == []
    assert deps["dgidb_client"].calls == []
    assert deps["mcsm_client"].calls == []


def test_invalid_amino_acid_letter_in_shorthand_produces_failure_not_a_crash():
    from secondlook.pipeline import run_tier2_pipeline

    deps = _deps()
    for bad_notation in ("Z175H", "R175Z"):
        result = run_tier2_pipeline("P04637", bad_notation, **deps)
        assert result.items == ()
        assert result.failure is not None
        assert result.failure["reason"] == OUT_OF_SCOPE_MESSAGE
    assert deps["dgidb_client"].calls == []


@pytest.mark.integration
def test_live_pipeline_tp53_r175h_shape_candidates_and_disclaimer():
    from secondlook.pipeline import TIER2_DISCLAIMER, run_tier2_pipeline

    result = run_tier2_pipeline("P04637", "R175H", max_candidates=1)
    assert result.failure is None
    assert result.items
    item = result.items[0]
    assert item["type"] == "computational_signal"
    assert item["disclaimer"] == TIER2_DISCLAIMER
    assert item["mutation_validated"]["status"] == "valid"
    assert item["mutation_validated"]["uniprot_accession"] == "P04637"
    assert item["mutation_validated"]["position"] == 175
    assert item["alphamissense"]["score"] == pytest.approx(0.9857, abs=1e-4)
    assert item["alphamissense"]["class"] == "likely_pathogenic"
    assert item["structure"]["source"] == "PDB"
    assert item["structure"]["id"]
    assert item["drug"]
    if item["method"] is None:
        assert item["label"] is None
        assert item["binding_note"]
    else:
        assert item["label"] == "uncertain"
    assert set(item) >= TIER2_RESULT_KEYS
