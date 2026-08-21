"""RCSB PDB search and entry metadata for structure sourcing."""

from __future__ import annotations

import os

from urllib.parse import quote

import httpx

DEFAULT_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DEFAULT_DATA_URL = "https://data.rcsb.org/rest/v1"
DEFAULT_FILES_URL = "https://files.rcsb.org/download"


class RcsbError(RuntimeError):
    pass


def inchikey_from_chemcomp(data: dict) -> str | None:
    """Pull InChIKey from a RCSB chemcomp JSON payload.

    Confirmed field on live responses (e.g. ligand 065): rcsb_chem_comp_descriptor.InChIKey.
    pdbx_chem_comp_descriptor entries with type InChIKey are the fallback.
    """
    desc = data.get("rcsb_chem_comp_descriptor") or {}
    key = desc.get("InChIKey") if isinstance(desc, dict) else None
    if isinstance(key, str) and key.strip():
        return key.strip()
    for item in data.get("pdbx_chem_comp_descriptor") or []:
        if not isinstance(item, dict):
            continue
        if item.get("type") == "InChIKey":
            descriptor = item.get("descriptor")
            if isinstance(descriptor, str) and descriptor.strip():
                return descriptor.strip()
    return None


def _covers_residue(pdb_text: str, position: int) -> bool:
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        if resseq == position:
            return True
    return False


class RcsbPdbClient:
    def __init__(
        self,
        search_url: str | None = None,
        data_url: str | None = None,
        client: httpx.Client | None = None,
    ) -> None:
        self.search_url = search_url or os.environ.get("RCSB_SEARCH_URL") or DEFAULT_SEARCH_URL
        self.data_url = (data_url or os.environ.get("RCSB_DATA_URL") or DEFAULT_DATA_URL).rstrip("/")
        self.files_url = (os.environ.get("RCSB_FILES_URL") or DEFAULT_FILES_URL).rstrip("/")
        self._client = client

    def search_by_uniprot(
        self,
        accession: str,
        preferred_ligands: tuple[str, ...] = (),
        *,
        position: int | None = None,
    ) -> dict | None:
        ligand_ids = self._search(accession, ligand_bound_only=True)
        if ligand_ids:
            hit = self._best_covering_hit(
                ligand_ids, ligand_bound=True, preferred_ligands=preferred_ligands, position=position
            )
            if hit is not None:
                return hit
            # None of the ligand-bound candidates cover the target residue (or
            # there was no target residue filter and this branch is unreached) —
            # fall through and also check apo structures rather than giving up;
            # a real apo structure covering the residue still beats no
            # structure at all (source_structure falls back to AlphaFold DB
            # next if this also comes back empty).

        apo_ids = self._search(accession, ligand_bound_only=False)
        if not apo_ids:
            return None
        return self._best_covering_hit(
            apo_ids, ligand_bound=False, preferred_ligands=preferred_ligands, position=position
        )

    def _best_covering_hit(
        self,
        pdb_ids: list[str],
        *,
        ligand_bound: bool,
        preferred_ligands: tuple[str, ...],
        position: int | None,
    ) -> dict | None:
        # Best-resolution-first order from _search is preserved by iterating in
        # place. If position is given, a structure that doesn't cover it is
        # useless for downstream binding scoring (chain_for_residue etc. all
        # assume the returned pdb_text actually spans the mutated residue) —
        # skip it rather than returning a structure that will silently fail
        # binding scoring for an unrelated reason later in the pipeline.
        first_covering: dict | None = None
        for pdb_id in pdb_ids:
            hit = self._entry_hit(pdb_id, ligand_bound=ligand_bound)
            if position is not None and not _covers_residue(hit["pdb_text"], position):
                continue
            if first_covering is None:
                first_covering = hit
            if preferred_ligands:
                ligands = {ligand.upper() for ligand in hit.get("ligands") or ()}
                if ligands & {ligand.upper() for ligand in preferred_ligands}:
                    return hit
                continue
            return hit
        return first_covering

    def _search(self, accession: str, *, ligand_bound_only: bool) -> list[str]:
        nodes: list[dict] = [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                    "operator": "exact_match",
                    "value": accession,
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                    "operator": "exact_match",
                    "value": "UniProt",
                },
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.structure_determination_methodology",
                    "operator": "exact_match",
                    "value": "experimental",
                },
            },
        ]
        if ligand_bound_only:
            nodes.append(
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.deposited_nonpolymer_entity_instance_count",
                        "operator": "greater",
                        "value": 0,
                    },
                }
            )
        payload = {
            "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": 0, "rows": 25},
                "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
            },
        }
        response = self._request("POST", self.search_url, json=payload)
        if response.status_code == 204:
            return []
        data = response.json()
        return [hit["identifier"] for hit in data.get("result_set") or []]

    def _entry_hit(self, pdb_id: str, *, ligand_bound: bool) -> dict:
        data = self._request("GET", f"{self.data_url}/core/entry/{pdb_id}").json()
        info = data.get("rcsb_entry_info") or {}
        identifiers = data.get("rcsb_entry_container_identifiers") or {}
        resolution = info.get("resolution_combined")
        if isinstance(resolution, list) and resolution:
            resolution = resolution[0]
        ligands = tuple(identifiers.get("non_polymer_entity_ids") or ())
        pdb_text = self._fetch_pdb_text(pdb_id)
        return {
            "pdb_id": pdb_id,
            "ligand_bound": ligand_bound or bool(identifiers.get("non_polymer_entity_ids")),
            "resolution": resolution,
            "ligands": ligands,
            "pdb_text": pdb_text,
        }

    def _fetch_pdb_text(self, pdb_id: str) -> str:
        # The data.rcsb.org entry endpoint above is metadata only — it has no
        # atom coordinates. The actual structure file lives at files.rcsb.org.
        response = self._request("GET", f"{self.files_url}/{pdb_id}.pdb")
        text = response.text
        if not text.strip():
            raise RcsbError(f"RCSB returned an empty structure file for {pdb_id}")
        return text

    def fetch_chemcomp_inchikey(self, het_code: str) -> str:
        url = f"{self.data_url}/core/chemcomp/{quote(het_code)}"
        response = self._request("GET", url)
        try:
            data = response.json()
        except ValueError as exc:
            raise RcsbError(f"RCSB chemcomp response for {het_code} was not JSON") from exc
        if not isinstance(data, dict):
            raise RcsbError(f"RCSB chemcomp response for {het_code} was not an object")
        key = inchikey_from_chemcomp(data)
        if not key:
            raise RcsbError(f"RCSB chemcomp response for {het_code} has no InChIKey")
        return key

    def _request(self, method: str, url: str, json: dict | None = None) -> httpx.Response:
        try:
            if self._client is not None:
                response = self._client.request(method, url, json=json, timeout=40.0)
            else:
                response = httpx.request(method, url, json=json, timeout=40.0, follow_redirects=True)
            if response.status_code != 204:
                response.raise_for_status()
        except httpx.HTTPError as exc:
            raise RcsbError(f"RCSB request failed for {url}") from exc
        return response
