"""AutoDock Vina docking fallback — WT vs mutant delta, separate from mCSM-lig."""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Protocol

from secondlook.binding import BindingScore, McsmPageError, chain_for_residue

DEFAULT_TIMEOUT_SECONDS = 120.0
DEFAULT_SEED = 1
DEFAULT_EXHAUSTIVENESS = 8
DEFAULT_PADDING_ANGSTROM = 5.0
PHYSIOLOGICAL_PH = 7.4

_MUTATION = re.compile(r"^([A-Z])(\d+)([A-Z])$")
_SKIP_HET = frozenset(
    {"HOH", "WAT", "DOD", "NA", "CL", "SO4", "PO4", "GOL", "EDO", "ACE", "NH2", "NME"}
)

ONE_TO_THREE = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


class VinaError(RuntimeError):
    """Vina/smina could not produce a WT-vs-mutant docking delta for this candidate."""


class ReceptorPrepError(VinaError):
    """Receptor protonation or PDBQT conversion failed."""


class MutantPlacementError(VinaError):
    """In-place side-chain mutation could not be applied defensibly."""


class LigandPrepError(VinaError):
    """SMILES could not be converted to a 3D PDBQT ligand."""


class NoBindingSiteError(VinaError):
    """No co-crystallized ligand coordinates to center a docking box."""


class VinaTimeoutError(VinaError):
    """A Vina run exceeded the configured timeout."""


class VinaRunError(VinaError):
    """Vina itself failed on the WT or mutant docking run."""


@dataclass(frozen=True)
class GridBox:
    center: tuple[float, float, float]
    size: tuple[float, float, float]


class ReceptorPreparer(Protocol):
    def to_pdbqt(self, pdb_text: str) -> str: ...


class MutantBuilder(Protocol):
    def place_sidechain(self, pdb_text: str, mutation: str, position: int) -> str: ...


class LigandPreparer(Protocol):
    def to_pdbqt(self, smiles: str) -> str: ...


class DockEngine(Protocol):
    def dock(
        self,
        receptor_pdbqt: str,
        ligand_pdbqt: str,
        box: GridBox,
        timeout_seconds: float,
        seed: int,
    ) -> float: ...


class Protonator(Protocol):
    def protonate(self, pdb_text: str) -> str: ...


def parse_mutation_shorthand(mutation: str) -> tuple[str, int, str]:
    match = _MUTATION.match(mutation.strip())
    if match is None:
        raise MutantPlacementError(f"Cannot parse mutation shorthand: {mutation}")
    return match.group(1), int(match.group(2)), match.group(3)


def ligand_hetatm_coords(pdb_text: str) -> tuple[tuple[float, float, float], ...]:
    groups: dict[tuple[str, str, str], list[tuple[float, float, float]]] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM") or len(line) < 54:
            continue
        resname = line[17:20].strip().upper()
        if resname in _SKIP_HET:
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        key = (line[21], line[22:26], resname)
        groups.setdefault(key, []).append((x, y, z))
    if not groups:
        return ()
    largest = max(groups.values(), key=len)
    return tuple(largest)


def grid_box_from_pdb(pdb_text: str, padding: float = DEFAULT_PADDING_ANGSTROM) -> GridBox:
    coords = ligand_hetatm_coords(pdb_text)
    if not coords:
        raise NoBindingSiteError(
            "No co-crystallized ligand coordinates in the structure; cannot center a Vina grid box"
        )
    xs, ys, zs = zip(*coords)
    center = (
        (min(xs) + max(xs)) / 2.0,
        (min(ys) + max(ys)) / 2.0,
        (min(zs) + max(zs)) / 2.0,
    )
    size = (
        (max(xs) - min(xs)) + 2.0 * padding,
        (max(ys) - min(ys)) + 2.0 * padding,
        (max(zs) - min(zs)) + 2.0 * padding,
    )
    return GridBox(center=center, size=size)


_PREFERRED_ALTLOCS = (" ", "A")


def protein_atoms_only(pdb_text: str) -> str:
    kept = [
        line
        for line in pdb_text.splitlines()
        if line.startswith(("TER", "END"))
        or (line.startswith("ATOM") and (line[16] in _PREFERRED_ALTLOCS))
    ]
    if not any(line.startswith("ATOM") for line in kept):
        raise ReceptorPrepError("PDB text has no ATOM records for receptor prep")
    return "\n".join(kept) + "\n"


def chain_atoms_only(pdb_text: str, chain: str) -> str:
    """Keep ATOM/HETATM records for one polymer chain (crystallographic copy)."""
    kept: list[str] = []
    for line in pdb_text.splitlines():
        if line.startswith(("TER", "END")):
            kept.append(line)
            continue
        if line.startswith(("ATOM", "HETATM")) and (line[21].strip() or "A") == chain:
            kept.append(line)
    if not any(line.startswith("ATOM") for line in kept):
        raise ReceptorPrepError(f"PDB text has no ATOM records for chain {chain}")
    return "\n".join(kept) + "\n"


def _is_hydrogen_atom_line(line: str) -> bool:
    if not line.startswith("ATOM"):
        return False
    element = line[76:78].strip() if len(line) >= 78 else ""
    name = line[12:16].strip()
    return element == "H" or (not element and name.startswith("H"))


def strip_hydrogens(pdb_text: str) -> str:
    """Drop hydrogen ATOM records. Meeko treats hydrogens as optional and adds its own."""
    kept = [line for line in pdb_text.splitlines() if not _is_hydrogen_atom_line(line)]
    return "\n".join(kept) + "\n"


class PdbFixerProtonator:
    def protonate(self, pdb_text: str) -> str:
        try:
            return _pdbfixer_complete(pdb_text, mutations=None, chain=None)
        except (ReceptorPrepError, MutantPlacementError):
            raise
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise ReceptorPrepError(f"Receptor protonation failed: {exc}") from exc


class PdbFixerMutantBuilder:
    """In-place side-chain swap via OpenMM PDBFixer templates; backbone kept, no re-fold."""

    def place_sidechain(self, pdb_text: str, mutation: str, position: int) -> str:
        wt, pos, mut = parse_mutation_shorthand(mutation)
        if pos != position:
            raise MutantPlacementError(f"Mutation {mutation} does not match position {position}")
        try:
            chain = chain_for_residue(pdb_text, position)
        except McsmPageError as exc:
            raise MutantPlacementError(f"Cannot resolve chain for residue {position}: {exc}") from exc
        expected = ONE_TO_THREE.get(wt)
        target = ONE_TO_THREE.get(mut)
        if expected is None or target is None:
            raise MutantPlacementError(f"Non-standard residue in mutation {mutation}")
        actual = _residue_name_at(pdb_text, chain, position)
        if actual != expected:
            raise MutantPlacementError(
                f"PDB residue {chain}:{position} is {actual}, not {expected} from {mutation}"
            )
        try:
            return _pdbfixer_complete(
                pdb_text,
                mutations=[f"{expected}-{position}-{target}"],
                chain=chain,
            )
        except MutantPlacementError:
            raise
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
            raise MutantPlacementError(f"PDBFixer could not place {mutation}: {exc}") from exc


def _update_h_positions_stable_indices(mol, indices_to_update: list[int]) -> None:
    """Place template hydrogens without relying on RDKit atom-object indices.

    meeko 0.7.1's update_H_positions keeps Atom objects across RWMol.RemoveAtom +
    Chem.AddHs. On RDKit 2026 those objects' GetIdx() no longer track the rewritten
    mol, so AddHs creates the hydrogens but only one of them is copied back
    ("Updated 1 H positions but deleted N"). Index through the heavy-atom order
    instead. Meeko docs: input hydrogens are optional; this only assigns coords
    for hydrogens the template already declared.
    """
    from rdkit import Chem

    if not indices_to_update:
        return
    conf = mol.GetConformer()
    parent_orig: dict[int, int] = {}
    for h_index in indices_to_update:
        atom = mol.GetAtomWithIdx(h_index)
        if atom.GetAtomicNum() != 1:
            raise RuntimeError("only H positions can be updated")
        parents = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]
        if len(parents) != 1:
            raise RuntimeError(f"hydrogens must have 1 non-H neighbor, got {len(parents)}")
        parent_orig[h_index] = parents[0]
    orig_heavy = [i for i in range(mol.GetNumAtoms()) if mol.GetAtomWithIdx(i).GetAtomicNum() != 1]
    orig_to_new = {orig: new for new, orig in enumerate(orig_heavy)}
    heavy = Chem.RWMol(mol)
    hs = [i for i, atom in enumerate(heavy.GetAtoms()) if atom.GetAtomicNum() == 1]
    for i in sorted(hs, reverse=True):
        heavy.RemoveAtom(i)
    heavy_mol = heavy.GetMol()
    heavy_mol.UpdatePropertyCache(strict=False)
    for orig_p in parent_orig.values():
        atom = heavy_mol.GetAtomWithIdx(orig_to_new[orig_p])
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
        atom.SetNoImplicit(True)
    only = sorted({orig_to_new[parent] for parent in parent_orig.values()})
    added = Chem.AddHs(heavy_mol, onlyOnAtoms=only, addCoords=True)
    added_conf = added.GetConformer()
    n_heavy = heavy_mol.GetNumAtoms()
    used: set[int] = set()
    for h_index, orig_p in parent_orig.items():
        new_p = orig_to_new[orig_p]
        for neigh in added.GetAtomWithIdx(new_p).GetNeighbors():
            if neigh.GetAtomicNum() == 1 and neigh.GetIdx() >= n_heavy and neigh.GetIdx() not in used:
                conf.SetAtomPosition(h_index, added_conf.GetAtomPosition(neigh.GetIdx()))
                used.add(neigh.GetIdx())
                break
    if len(used) != len(indices_to_update):
        raise RuntimeError(
            f"Updated {len(used)} H positions but deleted {len(indices_to_update)}"
        )


def _patch_meeko_hydrogen_placement() -> None:
    import meeko.polymer as polymer

    current = polymer.update_H_positions
    if getattr(current, "_secondlook_stable_indices", False):
        return
    polymer.update_H_positions = _update_h_positions_stable_indices
    _update_h_positions_stable_indices._secondlook_stable_indices = True  # type: ignore[attr-defined]


class MeekoReceptorPreparer:
    def to_pdbqt(self, pdb_text: str) -> str:
        _patch_meeko_hydrogen_placement()
        pdb_text = strip_hydrogens(pdb_text)
        try:
            from meeko import MoleculePreparation, PDBQTWriterLegacy, Polymer, ResidueChemTemplates
            from meeko.polymer import PolymerCreationError
        except ImportError as exc:
            raise ReceptorPrepError("meeko is not installed") from exc
        try:
            templates = ResidueChemTemplates.create_from_defaults()
            mk_prep = MoleculePreparation()
            polymer = Polymer.from_pdb_string(
                pdb_text,
                templates,
                mk_prep,
                allow_bad_res=True,
                default_altloc="A",
            )
            rigid, _flex = PDBQTWriterLegacy.write_string_from_polymer(polymer)
        except (PolymerCreationError, RuntimeError, ValueError, KeyError) as exc:
            raise ReceptorPrepError(f"meeko receptor prep failed: {exc}") from exc
        if not rigid.strip():
            raise ReceptorPrepError("meeko produced an empty receptor PDBQT")
        return rigid


class MeekoLigandPreparer:
    def to_pdbqt(self, smiles: str) -> str:
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from meeko import MoleculePreparation, PDBQTWriterLegacy
        except ImportError as exc:
            raise LigandPrepError("rdkit/meeko is not installed") from exc
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise LigandPrepError(f"Cannot parse SMILES: {smiles}")
        mol = Chem.AddHs(mol)
        embed = AllChem.EmbedMolecule(mol, randomSeed=DEFAULT_SEED)
        if embed != 0:
            raise LigandPrepError(f"Cannot generate a 3D conformer from SMILES: {smiles}")
        AllChem.UFFOptimizeMolecule(mol)
        setups = MoleculePreparation().prepare(mol)
        if not setups:
            raise LigandPrepError("meeko produced no ligand setup")
        pdbqt, ok, error_msg = PDBQTWriterLegacy.write_string(setups[0])
        if not ok or not pdbqt.strip():
            raise LigandPrepError(error_msg or "meeko failed to write ligand PDBQT")
        return pdbqt


class VinaEngine:
    def __init__(self, exhaustiveness: int = DEFAULT_EXHAUSTIVENESS) -> None:
        self.exhaustiveness = exhaustiveness

    def dock(
        self,
        receptor_pdbqt: str,
        ligand_pdbqt: str,
        box: GridBox,
        timeout_seconds: float,
        seed: int,
    ) -> float:
        src_root = str(Path(__file__).resolve().parents[1])
        with tempfile.TemporaryDirectory() as tmp:
            rec = Path(tmp) / "receptor.pdbqt"
            lig = Path(tmp) / "ligand.pdbqt"
            rec.write_text(receptor_pdbqt)
            lig.write_text(ligand_pdbqt)
            payload = json.dumps(
                {
                    "receptor_path": str(rec),
                    "ligand_path": str(lig),
                    "center": list(box.center),
                    "size": list(box.size),
                    "seed": seed,
                    "exhaustiveness": self.exhaustiveness,
                }
            )
            env = os.environ.copy()
            env["PYTHONPATH"] = src_root + os.pathsep + env.get("PYTHONPATH", "")
            try:
                completed = subprocess.run(
                    [sys.executable, "-c", _VINA_SUBPROCESS, payload],
                    capture_output=True,
                    text=True,
                    timeout=timeout_seconds,
                    env=env,
                    check=False,
                )
            except subprocess.TimeoutExpired as exc:
                raise VinaTimeoutError(f"Vina docking exceeded {timeout_seconds}s") from exc
        if completed.returncode != 0:
            detail = (completed.stderr or completed.stdout or "unknown Vina error").strip()
            raise VinaRunError(detail)
        try:
            return float(completed.stdout.strip().splitlines()[-1])
        except (ValueError, IndexError) as exc:
            raise VinaRunError(f"Vina returned no score: {completed.stdout!r}") from exc


class VinaDockClient:
    """Receptor/ligand prep and WT-vs-mutant docking with AutoDock Vina."""

    def __init__(
        self,
        *,
        receptor_preparer: ReceptorPreparer | None = None,
        mutant_builder: MutantBuilder | None = None,
        ligand_preparer: LigandPreparer | None = None,
        dock_engine: DockEngine | None = None,
        protonator: Protonator | None = None,
        timeout_seconds: float = DEFAULT_TIMEOUT_SECONDS,
        seed: int = DEFAULT_SEED,
        exhaustiveness: int = DEFAULT_EXHAUSTIVENESS,
    ) -> None:
        self._receptor_preparer = receptor_preparer or MeekoReceptorPreparer()
        self._mutant_builder = mutant_builder or PdbFixerMutantBuilder()
        self._ligand_preparer = ligand_preparer or MeekoLigandPreparer()
        self._dock_engine = dock_engine or VinaEngine(exhaustiveness=exhaustiveness)
        # Meeko's receptor templates add hydrogens themselves (docs: hydrogens are
        # optional on input). PDBFixer protonation before meeko is the path that
        # produced "Updated 1 H positions but deleted N" / padding errors on real
        # crystal structures, so it is not the default.
        self._protonator = protonator
        self.timeout_seconds = timeout_seconds
        self.seed = seed

    def score(
        self,
        *,
        pdb_text: str,
        smiles: str,
        mutation: str,
        position: int,
    ) -> BindingScore:
        box = grid_box_from_pdb(pdb_text)
        chain = chain_for_residue(pdb_text, position)
        protein = protein_atoms_only(chain_atoms_only(pdb_text, chain))
        wildtype_pdb = self._protonator.protonate(protein) if self._protonator else protein
        mutant_pdb = self._mutant_builder.place_sidechain(protein, mutation, position)
        wildtype_pdbqt = self._receptor_preparer.to_pdbqt(wildtype_pdb)
        mutant_pdbqt = self._receptor_preparer.to_pdbqt(mutant_pdb)
        ligand_pdbqt = self._ligand_preparer.to_pdbqt(smiles)
        wildtype_score = self._dock_engine.dock(
            wildtype_pdbqt, ligand_pdbqt, box, self.timeout_seconds, self.seed
        )
        mutant_score = self._dock_engine.dock(
            mutant_pdbqt, ligand_pdbqt, box, self.timeout_seconds, self.seed
        )
        return BindingScore(
            status="scored",
            method="docking",
            delta_score=mutant_score - wildtype_score,
            affinity_class=None,
            distance_angstrom=None,
            ligand_id=None,
            chain=None,
        )


def _residue_name_at(pdb_text: str, chain: str, position: int) -> str:
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if (line[21].strip() or "A") != chain:
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        if resseq == position:
            return line[17:20].strip().upper()
    raise MutantPlacementError(f"No residue {chain}:{position} in PDB")


def _pdbfixer_complete(pdb_text: str, mutations: list[str] | None, chain: str | None) -> str:
    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except ImportError as exc:
        raise ReceptorPrepError("pdbfixer/openmm is not installed") from exc
    fixer = PDBFixer(pdbfile=StringIO(pdb_text))
    if mutations and chain is not None:
        fixer.applyMutations(mutations, chain)
        # Filling truncated sidechains on *unrelated* residues (e.g. 8PQD A:761/A:763)
        # makes meeko treat a chain-break as a real polymer bond and fail with
        # "Expected 2 paddings ... but got 0". Only complete atoms on the residue
        # we actually mutated. Hydrogens stay off — meeko adds them from templates.
        # findMissingResidues() must still run: findMissingAtoms() reads
        # fixer.missingResidues even when we refuse to model the gaps.
        fixer.findMissingResidues()
        fixer.missingResidues = {}
        mutated_resseq = {part.split("-")[1] for part in mutations}
        fixer.findMissingAtoms()
        fixer.missingAtoms = {
            residue: atoms
            for residue, atoms in fixer.missingAtoms.items()
            if str(residue.id) in mutated_resseq
        }
        if hasattr(fixer, "missingTerminals"):
            fixer.missingTerminals = {
                residue: groups
                for residue, groups in fixer.missingTerminals.items()
                if str(residue.id) in mutated_resseq
            }
        fixer.addMissingAtoms()
    else:
        fixer.findMissingResidues()
        fixer.missingResidues = {}
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(PHYSIOLOGICAL_PH)
    out = StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, out, keepIds=True)
    return out.getvalue()


_VINA_SUBPROCESS = """
import json, sys
from secondlook.vina_dock import run_vina_from_files
args = json.loads(sys.argv[1])
print(run_vina_from_files(**args))
"""


def run_vina_from_files(
    *,
    receptor_path: str,
    ligand_path: str,
    center: list[float],
    size: list[float],
    seed: int,
    exhaustiveness: int,
) -> float:
    try:
        from vina import Vina
    except ImportError as exc:
        raise VinaRunError("vina is not installed") from exc
    engine = Vina(sf_name="vina", cpu=1, seed=seed, verbosity=0)
    engine.set_receptor(receptor_path)
    engine.set_ligand_from_file(ligand_path)
    engine.compute_vina_maps(center=center, box_size=size)
    engine.dock(exhaustiveness=exhaustiveness, n_poses=1)
    energies = engine.energies(n_poses=1)
    return float(energies[0][0])
