"""RCSB PDB search and entry metadata for structure sourcing."""

from __future__ import annotations

import os

import httpx

DEFAULT_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DEFAULT_DATA_URL = "https://data.rcsb.org/rest/v1"


class RcsbError(RuntimeError):
    pass


class RcsbPdbClient:
    def __init__(
        self,
        search_url: str | None = None,
        data_url: str | None = None,
        client: httpx.Client | None = None,
    ) -> None:
        self.search_url = search_url or os.environ.get("RCSB_SEARCH_URL") or DEFAULT_SEARCH_URL
        self.data_url = (data_url or os.environ.get("RCSB_DATA_URL") or DEFAULT_DATA_URL).rstrip("/")
        self._client = client

    def search_by_uniprot(self, accession: str, preferred_ligands: tuple[str, ...] = ()) -> dict | None:
        ligand_ids = self._search(accession, ligand_bound_only=True)
        apo_ids = [] if ligand_ids else self._search(accession, ligand_bound_only=False)
        candidates = ligand_ids or apo_ids
        if not candidates:
            return None

        chosen = candidates[0]
        if preferred_ligands:
            for pdb_id in candidates:
                meta = self._entry_hit(pdb_id, ligand_bound=bool(ligand_ids))
                ligands = {ligand.upper() for ligand in meta.get("ligands") or ()}
                if ligands & {ligand.upper() for ligand in preferred_ligands}:
                    return meta
        return self._entry_hit(chosen, ligand_bound=bool(ligand_ids))

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
        return {
            "pdb_id": pdb_id,
            "ligand_bound": ligand_bound or bool(identifiers.get("non_polymer_entity_ids")),
            "resolution": resolution,
            "ligands": ligands,
        }

    def _request(self, method: str, url: str, json: dict | None = None) -> httpx.Response:
        if self._client is not None:
            response = self._client.request(method, url, json=json, timeout=40.0)
        else:
            response = httpx.request(method, url, json=json, timeout=40.0, follow_redirects=True)
        try:
            if response.status_code != 204:
                response.raise_for_status()
        except httpx.HTTPError as exc:
            raise RcsbError(f"RCSB request failed for {url}") from exc
        return response
