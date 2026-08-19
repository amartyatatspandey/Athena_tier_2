"""Binding affinity scoring via mCSM-lig (Playwright) with Vina fallback (Tier 2 Step 5)."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Literal, Protocol

from secondlook.candidates import DrugCandidate
from secondlook.mutation_validation import MutationValidationResult
from secondlook.structure import StructureResult, structure_unavailable_message

# Distinct from structure_unavailable_message(): that one specifically claims the
# *structure* is low-pLDDT/unreliable, which would be false here — the structure
# (and AlphaMissense) may be perfectly good; only binding-affinity scoring itself
# came up empty (no HET code match, apo structure, mCSM-lig/Vina both failed).
BINDING_UNAVAILABLE_MESSAGE = (
    "No structure-based binding-affinity signal could be computed for this candidate "
    "(mCSM-lig and docking were both unavailable or inapplicable). AlphaMissense "
    "functional score and structural context remain available."
)


class McsmPageError(RuntimeError):
    """Playwright cannot reach the page or expected form selectors are missing."""


class McsmTimeoutError(RuntimeError):
    """mCSM-lig form submission exceeded the timeout."""


class McsmNoHetCodeError(RuntimeError):
    """Candidate has no PDB HET code, so mCSM-lig cannot run."""


class McsmApoStructureError(RuntimeError):
    """Structure has no bound ligand; mCSM-lig cannot run for any candidate."""


class McsmParseError(RuntimeError):
    """Result page did not contain the expected labeled fields."""


@dataclass(frozen=True)
class McsmLigResult:
    affinity_change: float
    affinity_class: str
    wild_type: str
    position: int
    mutant_type: str
    chain: str
    ligand_id: str
    distance_angstrom: float
    duet_stability_kcal: float


@dataclass(frozen=True)
class BindingScore:
    status: Literal["scored", "unavailable"]
    method: Literal["mCSM-lig", "docking"] | None
    delta_score: float | None
    affinity_class: str | None
    distance_angstrom: float | None
    ligand_id: str | None
    chain: str | None
    error_message: str | None = None


class McsmClient(Protocol):
    def submit(
        self,
        *,
        pdb_code: str | None,
        pdb_text: str | None,
        mutation: str,
        chain: str,
        lig_id: str,
    ) -> McsmLigResult: ...


class VinaClient(Protocol):
    def score(
        self,
        *,
        pdb_text: str,
        smiles: str,
        mutation: str,
        position: int,
    ) -> BindingScore: ...


class HetCodeResolver(Protocol):
    def resolve(self, drug_name: str, structure: StructureResult) -> str | None: ...


_AFFINITY = re.compile(
    r"Predicted Affinity Change:\s*([+-]?\d+(?:\.\d+)?)\s*log\(affinity fold change\)\s*-\s*(Destabilizing|Stabilizing)",
    re.IGNORECASE,
)
_WILD = re.compile(r"Wild-type:\s*([A-Z])", re.IGNORECASE)
_POSITION = re.compile(r"Position:\s*(\d+)", re.IGNORECASE)
_MUTANT = re.compile(r"Mutant-type:\s*([A-Z])", re.IGNORECASE)
_CHAIN = re.compile(r"Chain:\s*([A-Za-z0-9])", re.IGNORECASE)
_LIGAND = re.compile(r"Ligand ID:\s*([A-Za-z0-9]+)", re.IGNORECASE)
_DISTANCE = re.compile(r"Distance to ligand:\s*([+-]?\d+(?:\.\d+)?)\s*Å", re.IGNORECASE)
_DUET = re.compile(r"DUET stability change:\s*([+-]?\d+(?:\.\d+)?)\s*Kcal/mol", re.IGNORECASE)


def parse_mcsm_result(html: str) -> McsmLigResult:
    text = re.sub(r"<[^>]+>", " ", html)
    text = re.sub(r"\s+", " ", text)

    def _req(pattern: re.Pattern[str], field: str) -> re.Match[str]:
        match = pattern.search(text)
        if match is None:
            raise McsmParseError(f"mCSM-lig result page missing {field}")
        return match

    affinity = _req(_AFFINITY, "Predicted Affinity Change")
    return McsmLigResult(
        affinity_change=float(affinity.group(1)),
        affinity_class=affinity.group(2).title(),
        wild_type=_req(_WILD, "Wild-type").group(1).upper(),
        position=int(_req(_POSITION, "Position").group(1)),
        mutant_type=_req(_MUTANT, "Mutant-type").group(1).upper(),
        chain=_req(_CHAIN, "Chain").group(1).upper(),
        ligand_id=_req(_LIGAND, "Ligand ID").group(1),
        distance_angstrom=float(_req(_DISTANCE, "Distance to ligand").group(1)),
        duet_stability_kcal=float(_req(_DUET, "DUET stability change").group(1)),
    )


def mutation_shorthand(validation: MutationValidationResult) -> str:
    if validation.position is None or not validation.reference_residue_expected:
        raise ValueError("validated mutation is missing residue/position")
    alt = (
        validation.mutant_sequence[validation.position - 1]
        if validation.mutant_sequence
        else ""
    )
    if not alt:
        raise ValueError("validated mutation is missing mutant residue")
    return f"{validation.reference_residue_expected}{validation.position}{alt}"


def chain_for_residue(pdb_text: str, residue_number: int) -> str:
    covering: list[str] = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        if resseq != residue_number:
            continue
        chain = line[21].strip() or "A"
        if chain not in covering:
            covering.append(chain)
    if not covering:
        raise McsmPageError(f"No chain in structure covers residue {residue_number}")
    if "A" in covering:
        return "A"
    return covering[0]


def het_codes_from_pdb(pdb_text: str) -> tuple[str, ...]:
    codes: list[str] = []
    skip = {"HOH", "WAT", "DOD", "NA", "CL", "SO4", "PO4", "GOL", "EDO"}
    for line in pdb_text.splitlines():
        if not line.startswith(("HETATM", "HET ")):
            continue
        resname = line[17:20].strip() if line.startswith("HETATM") else line[7:10].strip()
        if not resname or resname.upper() in skip:
            continue
        if resname not in codes:
            codes.append(resname)
    return tuple(codes)


class StructureHetResolver:
    """Resolve a candidate to a PDB HET code; None means mCSM-lig cannot run."""

    def resolve(self, drug_name: str, structure: StructureResult) -> str | None:
        if not structure.pdb_text:
            return None
        codes = het_codes_from_pdb(structure.pdb_text)
        if not codes:
            return None
        needle = re.sub(r"[^A-Z0-9]", "", drug_name.upper())
        for code in codes:
            if code.upper() in needle or needle.endswith(code.upper()):
                return code
        return None


def score_binding(
    validation: MutationValidationResult,
    structure: StructureResult,
    candidate: DrugCandidate,
    *,
    mcsm_client: McsmClient | None = None,
    vina_client: VinaClient | None = None,
    het_resolver: HetCodeResolver | None = None,
    min_delay_seconds: float = 3.0,
    sleeper=None,
) -> BindingScore:
    from secondlook.vina_dock import VinaError

    if sleeper is None:
        import time

        sleeper = time.sleep
    if mcsm_client is None:
        from secondlook.mcsm_lig import McsmLigPlaywrightClient

        mcsm_client = McsmLigPlaywrightClient()
    if vina_client is None:
        from secondlook.vina_dock import VinaDockClient

        vina_client = VinaDockClient()
    if het_resolver is None:
        het_resolver = StructureHetResolver()

    if structure.status != "found":
        return _unavailable(structure_unavailable_message(structure.plddt_at_residue))
    if validation.position is None or not structure.pdb_text:
        return _unavailable(BINDING_UNAVAILABLE_MESSAGE)

    try:
        if not structure.ligand_bound:
            raise McsmApoStructureError("Structure has no bound ligand")
        het = het_resolver.resolve(candidate.name, structure)
        if not het:
            raise McsmNoHetCodeError(f"No PDB HET code for {candidate.name}")
        chain = chain_for_residue(structure.pdb_text, validation.position)
        pdb_code = (
            structure.id
            if structure.source == "PDB" and structure.id and len(structure.id) == 4
            else None
        )
        sleeper(min_delay_seconds)
        parsed = mcsm_client.submit(
            pdb_code=pdb_code,
            pdb_text=None if pdb_code else structure.pdb_text,
            mutation=mutation_shorthand(validation),
            chain=chain,
            lig_id=het,
        )
    except (
        McsmApoStructureError,
        McsmNoHetCodeError,
        McsmPageError,
        McsmTimeoutError,
        McsmParseError,
    ):
        return _vina_or_unavailable(
            validation, structure, candidate, vina_client, VinaError
        )

    return BindingScore(
        status="scored",
        method="mCSM-lig",
        delta_score=parsed.affinity_change,
        affinity_class=parsed.affinity_class,
        distance_angstrom=parsed.distance_angstrom,
        ligand_id=parsed.ligand_id,
        chain=parsed.chain,
    )


def _vina_or_unavailable(
    validation: MutationValidationResult,
    structure: StructureResult,
    candidate: DrugCandidate,
    vina_client: VinaClient,
    vina_error: type[BaseException],
) -> BindingScore:
    if not candidate.smiles or not structure.pdb_text or validation.position is None:
        return _unavailable(BINDING_UNAVAILABLE_MESSAGE)
    try:
        return vina_client.score(
            pdb_text=structure.pdb_text,
            smiles=candidate.smiles,
            mutation=mutation_shorthand(validation),
            position=validation.position,
        )
    except vina_error:
        return _unavailable(BINDING_UNAVAILABLE_MESSAGE)


def _unavailable(message: str) -> BindingScore:
    return BindingScore(
        status="unavailable",
        method=None,
        delta_score=None,
        affinity_class=None,
        distance_angstrom=None,
        ligand_id=None,
        chain=None,
        error_message=message,
    )
