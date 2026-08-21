import pytest

from secondlook.binding import BindingScore
from secondlook.vina_dock import (
    LigandPrepError,
    MutantPlacementError,
    NoBindingSiteError,
    ReceptorPrepError,
    VinaDockClient,
    VinaError,
    VinaRunError,
    VinaTimeoutError,
    grid_box_from_pdb,
    ligand_hetatm_coords,
    parse_mutation_shorthand,
    protein_atoms_only,
)

LIGAND_PDB = """\
ATOM      1  CA  ASP A  30      1.000   0.000   0.000  1.00 20.00           C
HETATM    2  C1  LIG A  99      0.000   0.000   0.000  1.00 20.00           C
HETATM    3  C2  LIG A  99     10.000   0.000   0.000  1.00 20.00           C
HETATM    4  O   HOH A 100     99.000  99.000  99.000  1.00 20.00           O
END
"""

APO_PDB = """\
ATOM      1  CA  ASP A  30      1.000   0.000   0.000  1.00 20.00           C
END
"""

PEPTIDE_PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00 20.00           O
ATOM      5  CB  ALA A   1       2.009  -0.771  -1.214  1.00 20.00           C
ATOM      6  N   ASP A   2       3.458   1.654   0.000  1.00 20.00           N
ATOM      7  CA  ASP A   2       4.200   2.900   0.000  1.00 20.00           C
ATOM      8  C   ASP A   2       5.700   2.700   0.000  1.00 20.00           C
ATOM      9  O   ASP A   2       6.200   1.600   0.000  1.00 20.00           O
ATOM     10  CB  ASP A   2       3.700   3.700   1.200  1.00 20.00           C
ATOM     11  CG  ASP A   2       2.200   4.000   1.200  1.00 20.00           C
ATOM     12  OD1 ASP A   2       1.400   3.200   0.600  1.00 20.00           O
ATOM     13  OD2 ASP A   2       1.900   5.100   1.700  1.00 20.00           O
ATOM     14  N   ALA A   3       6.400   3.800   0.000  1.00 20.00           N
ATOM     15  CA  ALA A   3       7.850   3.800   0.000  1.00 20.00           C
ATOM     16  C   ALA A   3       8.400   5.200   0.000  1.00 20.00           C
ATOM     17  O   ALA A   3       7.650   6.200   0.000  1.00 20.00           O
ATOM     18  CB  ALA A   3       8.400   3.000  -1.200  1.00 20.00           C
HETATM   19  C1  LIG A  99      4.200   2.900   5.000  1.00 20.00           C
END
"""


class FakePrep:
    def __init__(self, error: Exception | None = None) -> None:
        self.error = error
        self.calls: list[str] = []

    def to_pdbqt(self, pdb_text: str) -> str:
        self.calls.append(pdb_text)
        if self.error:
            raise self.error
        return f"PDBQT:{hash(pdb_text)}"


class FakeMutant:
    def __init__(self, error: Exception | None = None, pdb: str = "MUTANT") -> None:
        self.error = error
        self.pdb = pdb
        self.calls: list[tuple[str, str, int]] = []

    def place_sidechain(self, pdb_text: str, mutation: str, position: int) -> str:
        self.calls.append((pdb_text, mutation, position))
        if self.error:
            raise self.error
        return self.pdb


class FakeLigand:
    def __init__(self, error: Exception | None = None) -> None:
        self.error = error
        self.calls: list[str] = []

    def to_pdbqt(self, smiles: str) -> str:
        self.calls.append(smiles)
        if self.error:
            raise self.error
        return "LIGAND_PDBQT"


class FakeDock:
    def __init__(self, scores: list[float] | None = None, error: Exception | None = None) -> None:
        self.scores = list(scores or [])
        self.error = error
        self.calls: list[dict] = []

    def dock(self, receptor_pdbqt: str, ligand_pdbqt: str, box, timeout_seconds: float, seed: int) -> float:
        self.calls.append(
            {
                "receptor": receptor_pdbqt,
                "ligand": ligand_pdbqt,
                "box": box,
                "timeout": timeout_seconds,
                "seed": seed,
            }
        )
        if self.error:
            raise self.error
        return self.scores.pop(0)


def test_parse_mutation_shorthand():
    assert parse_mutation_shorthand("D30N") == ("D", 30, "N")


def test_ligand_coords_ignore_water_and_use_hetatm():
    coords = ligand_hetatm_coords(LIGAND_PDB)
    assert coords == ((0.0, 0.0, 0.0), (10.0, 0.0, 0.0))


def test_grid_box_is_ligand_bounds_plus_five_angstrom_padding():
    box = grid_box_from_pdb(LIGAND_PDB, padding=5.0)
    assert box.center == pytest.approx((5.0, 0.0, 0.0))
    assert box.size == pytest.approx((20.0, 10.0, 10.0))


def test_grid_box_without_ligand_raises_no_binding_site():
    with pytest.raises(NoBindingSiteError):
        grid_box_from_pdb(APO_PDB)


def test_protein_atoms_only_strips_heteroatoms():
    protein = protein_atoms_only(LIGAND_PDB)
    assert "HETATM" not in protein
    assert "CA" in protein


ALTLOC_PDB = """\
ATOM      1  N  AARG A  30      0.000   0.000   0.000  0.60 20.00           N
ATOM      2  N  BARG A  30      0.100   0.100   0.100  0.40 20.00           N
ATOM      3  CA  ARG A  30      1.000   0.000   0.000  1.00 20.00           C
"""


def test_protein_atoms_only_keeps_one_altloc_not_both():
    # Regression: real crystal structures (e.g. PDB 8PQD, KIT D816V) have
    # alternate-conformation atoms (altLoc A/B) for the same position. Feeding
    # both copies to meeko's automatic bond/valence perception produced
    # "Explicit valence ... greater than permitted" — meeko can't tell the two
    # overlapping N atoms apart. Keep only the blank/'A' altLoc copy.
    protein = protein_atoms_only(ALTLOC_PDB)
    n_lines = [line for line in protein.splitlines() if line.startswith("ATOM") and line[12:16].strip() == "N"]
    assert len(n_lines) == 1
    assert n_lines[0][16] == "A"


def test_score_reports_mutant_minus_wildtype_delta():
    dock = FakeDock(scores=[-8.0, -6.5])
    client = VinaDockClient(
        receptor_preparer=FakePrep(),
        mutant_builder=FakeMutant(),
        ligand_preparer=FakeLigand(),
        dock_engine=dock,
    )
    result = client.score(pdb_text=LIGAND_PDB, smiles="CCO", mutation="D30N", position=30)
    assert result.status == "scored"
    assert result.method == "docking"
    assert result.delta_score == pytest.approx(1.5)
    assert len(dock.calls) == 2
    assert dock.calls[0]["receptor"] != dock.calls[1]["receptor"]
    assert dock.calls[0]["ligand"] == "LIGAND_PDBQT"


def test_apo_structure_raises_no_binding_site_not_a_fake_box():
    client = VinaDockClient(
        receptor_preparer=FakePrep(),
        mutant_builder=FakeMutant(),
        ligand_preparer=FakeLigand(),
        dock_engine=FakeDock(scores=[-1.0, -1.0]),
    )
    with pytest.raises(NoBindingSiteError):
        client.score(pdb_text=APO_PDB, smiles="CCO", mutation="D30N", position=30)


def test_receptor_prep_failure_raises_receptor_prep_error():
    client = VinaDockClient(
        receptor_preparer=FakePrep(error=ReceptorPrepError("bad protonation")),
        mutant_builder=FakeMutant(),
        ligand_preparer=FakeLigand(),
        dock_engine=FakeDock(scores=[-1.0, -1.0]),
    )
    with pytest.raises(ReceptorPrepError):
        client.score(pdb_text=LIGAND_PDB, smiles="CCO", mutation="D30N", position=30)


def test_mutant_placement_failure_raises_mutant_placement_error():
    client = VinaDockClient(
        receptor_preparer=FakePrep(),
        mutant_builder=FakeMutant(error=MutantPlacementError("no rotamer")),
        ligand_preparer=FakeLigand(),
        dock_engine=FakeDock(scores=[-1.0, -1.0]),
    )
    with pytest.raises(MutantPlacementError):
        client.score(pdb_text=LIGAND_PDB, smiles="CCO", mutation="D30N", position=30)


def test_unembeddable_smiles_raises_ligand_prep_error():
    client = VinaDockClient(
        receptor_preparer=FakePrep(),
        mutant_builder=FakeMutant(),
        ligand_preparer=FakeLigand(error=LigandPrepError("embed failed")),
        dock_engine=FakeDock(scores=[-1.0, -1.0]),
    )
    with pytest.raises(LigandPrepError):
        client.score(pdb_text=LIGAND_PDB, smiles="not-a-smiles", mutation="D30N", position=30)


def test_vina_timeout_raises_timeout_error():
    client = VinaDockClient(
        receptor_preparer=FakePrep(),
        mutant_builder=FakeMutant(),
        ligand_preparer=FakeLigand(),
        dock_engine=FakeDock(error=VinaTimeoutError("dock hung")),
    )
    with pytest.raises(VinaTimeoutError):
        client.score(pdb_text=LIGAND_PDB, smiles="CCO", mutation="D30N", position=30)


def test_vina_run_failure_raises_vina_run_error():
    client = VinaDockClient(
        receptor_preparer=FakePrep(),
        mutant_builder=FakeMutant(),
        ligand_preparer=FakeLigand(),
        dock_engine=FakeDock(error=VinaRunError("maps failed")),
    )
    with pytest.raises(VinaRunError):
        client.score(pdb_text=LIGAND_PDB, smiles="CCO", mutation="D30N", position=30)


def test_specific_errors_are_vina_errors_for_orchestrator_catch():
    for exc in (
        ReceptorPrepError("x"),
        MutantPlacementError("x"),
        LigandPrepError("x"),
        NoBindingSiteError("x"),
        VinaTimeoutError("x"),
        VinaRunError("x"),
    ):
        assert isinstance(exc, VinaError)


def test_pdbfixer_places_mutant_sidechain_without_moving_backbone():
    from secondlook.vina_dock import PdbFixerMutantBuilder

    builder = PdbFixerMutantBuilder()
    mutant = builder.place_sidechain(PEPTIDE_PDB, "D2N", 2)
    ca_wt = [line[30:54] for line in PEPTIDE_PDB.splitlines() if line.startswith("ATOM") and line[12:16].strip() == "CA" and int(line[22:26]) == 2][0]
    ca_mut = [line[30:54] for line in mutant.splitlines() if line.startswith("ATOM") and line[12:16].strip() == "CA" and int(line[22:26]) == 2][0]
    assert ca_mut == ca_wt
    res2 = [line for line in mutant.splitlines() if line.startswith("ATOM") and int(line[22:26]) == 2]
    assert any("ASN" in line[17:20] for line in res2)
    assert any(line[12:16].strip() == "ND2" for line in res2)


def test_meeko_ligand_prep_rejects_unembeddable_smiles():
    from secondlook.vina_dock import MeekoLigandPreparer

    with pytest.raises(LigandPrepError):
        MeekoLigandPreparer().to_pdbqt("not-a-smiles")


def test_score_binding_maps_vina_errors_to_binding_unavailable_message():
    from secondlook.binding import BINDING_UNAVAILABLE_MESSAGE, score_binding
    from secondlook.candidates import DrugCandidate
    from secondlook.mutation_validation import MutationValidationResult
    from secondlook.structure import StructureResult

    validation = MutationValidationResult(
        status="valid",
        gene="TEST",
        hgvs_normalized="p.Asp30Asn",
        uniprot_accession="P00000",
        isoform_note="canonical",
        position=30,
        reference_residue_expected="D",
        reference_residue_actual="D",
        mutation_type="missense",
        error_message=None,
        wildtype_sequence="D" + "A" * 40,
        mutant_sequence="N" + "A" * 40,
    )
    structure = StructureResult(
        status="found",
        source="PDB",
        id="2Z4O",
        plddt_at_residue=90.0,
        plddt_global=90.0,
        reliability_flag="high",
        ligand_bound=True,
        annotated_position=30,
        pdb_text=LIGAND_PDB,
    )
    candidate = DrugCandidate(
        name="EXAMPLE",
        source="DGIdb",
        target_tier="exact_protein",
        approved=True,
        smiles="CCO",
        smiles_source="PubChem",
    )

    class RaisingVina:
        def score(self, **kwargs) -> BindingScore:
            raise ReceptorPrepError("bad atoms")

    result = score_binding(
        validation,
        structure,
        candidate,
        mcsm_client=type("M", (), {"submit": staticmethod(lambda **k: (_ for _ in ()).throw(Exception("unused")))})(),
        vina_client=RaisingVina(),
        het_resolver=type("H", (), {"resolve": staticmethod(lambda *a: None)})(),
        sleeper=lambda _: None,
    )
    assert result.status == "unavailable"
    assert result.error_message == BINDING_UNAVAILABLE_MESSAGE


@pytest.mark.integration
def test_9c5s_tp53_r175h_cannot_place_uniprot_residue_on_fragment():
    import httpx

    pdb_text = httpx.get("https://files.rcsb.org/download/9C5S.pdb", timeout=60.0).text
    client = VinaDockClient(timeout_seconds=30.0, seed=1, exhaustiveness=1)
    with pytest.raises((MutantPlacementError, NoBindingSiteError)):
        client.score(pdb_text=pdb_text, smiles="Fc1c[nH]c(=O)[nH]c1=O", mutation="R175H", position=175)


@pytest.mark.integration
def test_live_vina_delta_on_2z4o_with_small_ligand():
    import httpx

    # 9C5S (Step 3's TP53 hit) is a 14-residue fragment without UniProt 175 and
    # without a drug-like HET ligand, so it cannot host a defensible R175H dock.
    # 2Z4O is the verified ligand-bound complex (Prompt E example) and is the
    # smallest real structure where WT-vs-mutant Vina can actually run.
    pdb_text = httpx.get("https://files.rcsb.org/download/2Z4O.pdb", timeout=60.0).text
    client = VinaDockClient(timeout_seconds=180.0, seed=1, exhaustiveness=1)
    result = client.score(pdb_text=pdb_text, smiles="Fc1c[nH]c(=O)[nH]c1=O", mutation="D30N", position=30)
    assert result.status == "scored"
    assert result.method == "docking"
    assert result.delta_score is not None
    assert result.delta_score == result.delta_score
