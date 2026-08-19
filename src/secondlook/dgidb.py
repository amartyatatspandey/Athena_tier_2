"""DGIdb GraphQL client for exact-protein drug interactions."""

from __future__ import annotations

import os

import httpx

DEFAULT_DGIDB_URL = "https://dgidb.org/api/graphql"

_DRUGS_QUERY = """
query GeneDrugs($names: [String!]!) {
  genes(names: $names) {
    nodes {
      name
      interactions {
        interactionScore
        drug { name approved }
      }
    }
  }
}
"""


class DgidbError(RuntimeError):
    pass


class DgidbClient:
    def __init__(self, url: str | None = None, client: httpx.Client | None = None) -> None:
        self.url = url or os.environ.get("DGIDB_GRAPHQL_URL") or DEFAULT_DGIDB_URL
        self._client = client

    def fetch_drugs(self, gene_symbol: str) -> list[dict]:
        payload = self._post({"query": _DRUGS_QUERY, "variables": {"names": [gene_symbol]}})
        nodes = ((payload.get("data") or {}).get("genes") or {}).get("nodes") or []
        rows: list[dict] = []
        for node in nodes:
            for interaction in node.get("interactions") or []:
                drug = interaction.get("drug") or {}
                name = drug.get("name")
                if not name:
                    continue
                rows.append(
                    {
                        "name": name,
                        "approved": bool(drug.get("approved")),
                        "score": float(interaction.get("interactionScore") or 0),
                    }
                )
        return rows

    def _post(self, body: dict) -> dict:
        try:
            if self._client is not None:
                response = self._client.post(self.url, json=body, timeout=40.0)
            else:
                response = httpx.post(self.url, json=body, timeout=40.0, follow_redirects=True)
            response.raise_for_status()
            payload = response.json()
        except (httpx.HTTPError, ValueError) as exc:
            raise DgidbError(f"DGIdb request failed for {self.url}") from exc
        if payload.get("errors"):
            raise DgidbError(f"DGIdb GraphQL error: {payload['errors']}")
        return payload
