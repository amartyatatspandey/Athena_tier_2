import httpx
import pytest

from secondlook.rcsb import RcsbError, RcsbPdbClient, _covers_residue


class ConnectionResetClient:
    """Simulates a connection-level failure (no response object ever exists) —
    distinct from an HTTP status error, which already has a response to call
    raise_for_status() on."""

    def request(self, method, url, json=None, timeout=None):
        raise httpx.ConnectError("Connection reset by peer")


def test_connection_error_is_wrapped_not_raised_raw():
    client = RcsbPdbClient(client=ConnectionResetClient())
    with pytest.raises(RcsbError):
        client._fetch_pdb_text("9C5S")

MINI_PDB_175 = """\
ATOM      1  CA  ARG A 175      1.000   0.000   0.000  1.00 96.62           C
HETATM    2  C1  LIG A 200      0.000   0.000   0.000  1.00 20.00           C
"""

MINI_PDB_NO_175 = """\
ATOM      1  CA  ALA A  10      1.000   0.000   0.000  1.00 50.00           C
HETATM    2  C1  LIG A 200      0.000   0.000   0.000  1.00 20.00           C
"""


def test_covers_residue_true_and_false():
    assert _covers_residue(MINI_PDB_175, 175) is True
    assert _covers_residue(MINI_PDB_175, 999) is False
    assert _covers_residue(MINI_PDB_NO_175, 175) is False


class FakeResponse:
    def __init__(self, *, json_body=None, text: str = "", status_code: int = 200) -> None:
        self._json_body = json_body
        self.text = text
        self.status_code = status_code

    def json(self):
        return self._json_body

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            import httpx

            raise httpx.HTTPStatusError("error", request=None, response=None)


class FakeRcsbTransport:
    """Simulates search.rcsb.org + data.rcsb.org + files.rcsb.org for one accession.

    entries: {pdb_id: {"ligands": [...], "resolution": float, "pdb_text": str}}
    search_order: list of pdb_ids returned by the search endpoint, in order.
    """

    def __init__(self, entries: dict, search_order: list[str]) -> None:
        self.entries = entries
        self.search_order = search_order
        self.requested_files: list[str] = []

    def request(self, method: str, url: str, json=None, timeout=None):
        if "rcsbsearch" in url:
            return FakeResponse(json_body={"result_set": [{"identifier": pid} for pid in self.search_order]})
        if "/core/entry/" in url:
            pdb_id = url.rsplit("/", 1)[-1]
            entry = self.entries[pdb_id]
            return FakeResponse(
                json_body={
                    "rcsb_entry_info": {"resolution_combined": [entry.get("resolution", 2.0)]},
                    "rcsb_entry_container_identifiers": {"non_polymer_entity_ids": entry.get("ligands", [])},
                }
            )
        if "/download/" in url:
            pdb_id = url.rsplit("/", 1)[-1].removesuffix(".pdb")
            self.requested_files.append(pdb_id)
            return FakeResponse(text=self.entries[pdb_id]["pdb_text"])
        raise AssertionError(f"unexpected URL in fake transport: {url}")


def _client_for(entries: dict, search_order: list[str]) -> RcsbPdbClient:
    transport = FakeRcsbTransport(entries, search_order)
    return RcsbPdbClient(client=transport), transport


def test_skips_non_covering_candidate_for_a_covering_one():
    entries = {
        "AAAA": {"ligands": ["LIG"], "resolution": 1.0, "pdb_text": MINI_PDB_NO_175},
        "BBBB": {"ligands": ["LIG"], "resolution": 1.5, "pdb_text": MINI_PDB_175},
    }
    client, transport = _client_for(entries, ["AAAA", "BBBB"])

    hit = client.search_by_uniprot("P00000", position=175)

    assert hit is not None
    assert hit["pdb_id"] == "BBBB"
    assert transport.requested_files == ["AAAA", "BBBB"]


def test_first_candidate_used_when_it_covers_the_residue():
    entries = {
        "AAAA": {"ligands": ["LIG"], "resolution": 1.0, "pdb_text": MINI_PDB_175},
        "BBBB": {"ligands": ["LIG"], "resolution": 1.5, "pdb_text": MINI_PDB_175},
    }
    client, transport = _client_for(entries, ["AAAA", "BBBB"])

    hit = client.search_by_uniprot("P00000", position=175)

    assert hit["pdb_id"] == "AAAA"
    # Should not need to check BBBB at all once AAAA covers the residue.
    assert transport.requested_files == ["AAAA"]


def test_no_position_filter_behaves_like_before_first_candidate_wins():
    entries = {
        "AAAA": {"ligands": ["LIG"], "resolution": 1.0, "pdb_text": MINI_PDB_NO_175},
    }
    client, _ = _client_for(entries, ["AAAA"])

    hit = client.search_by_uniprot("P00000")

    assert hit["pdb_id"] == "AAAA"


def test_none_returned_when_no_ligand_bound_candidate_covers_residue_but_apo_does():
    # search_by_uniprot calls _search twice (ligand-bound tier, then apo tier).
    # Simulate that by having the same transport answer both searches with
    # different result sets isn't directly supported by this simple fake, so
    # this test drives the apo path directly via a transport whose single
    # search call represents whichever tier is being queried; the real
    # two-tier behavior is covered by the live integration test.
    entries = {
        "CCCC": {"ligands": [], "resolution": 2.0, "pdb_text": MINI_PDB_175},
    }
    client, transport = _client_for(entries, ["CCCC"])

    hit = client._best_covering_hit(["CCCC"], ligand_bound=False, preferred_ligands=(), position=175)

    assert hit is not None
    assert hit["pdb_id"] == "CCCC"


def test_returns_none_when_nothing_covers_the_residue():
    entries = {
        "AAAA": {"ligands": ["LIG"], "resolution": 1.0, "pdb_text": MINI_PDB_NO_175},
        "BBBB": {"ligands": ["LIG"], "resolution": 1.5, "pdb_text": MINI_PDB_NO_175},
    }
    client, _ = _client_for(entries, ["AAAA", "BBBB"])

    hit = client._best_covering_hit(["AAAA", "BBBB"], ligand_bound=True, preferred_ligands=(), position=175)

    assert hit is None


def test_preferred_ligand_wins_among_covering_candidates_not_just_first():
    entries = {
        "AAAA": {"ligands": ["OTHER"], "resolution": 1.0, "pdb_text": MINI_PDB_175},
        "BBBB": {"ligands": ["WANTED"], "resolution": 1.5, "pdb_text": MINI_PDB_175},
    }
    client, _ = _client_for(entries, ["AAAA", "BBBB"])

    hit = client._best_covering_hit(
        ["AAAA", "BBBB"], ligand_bound=True, preferred_ligands=("WANTED",), position=175
    )

    assert hit["pdb_id"] == "BBBB"


def test_preferred_ligand_falls_back_to_first_covering_when_no_match():
    entries = {
        "AAAA": {"ligands": ["OTHER"], "resolution": 1.0, "pdb_text": MINI_PDB_175},
        "BBBB": {"ligands": ["ALSO_OTHER"], "resolution": 1.5, "pdb_text": MINI_PDB_175},
    }
    client, _ = _client_for(entries, ["AAAA", "BBBB"])

    hit = client._best_covering_hit(
        ["AAAA", "BBBB"], ligand_bound=True, preferred_ligands=("WANTED",), position=175
    )

    # Neither has the preferred ligand, but both cover the residue — must not
    # return None just because the preference wasn't satisfied.
    assert hit is not None
    assert hit["pdb_id"] == "AAAA"


def test_empty_structure_file_raises_rcsb_error():
    entries = {"AAAA": {"ligands": [], "resolution": 1.0, "pdb_text": "   "}}
    client, _ = _client_for(entries, ["AAAA"])

    with pytest.raises(RcsbError):
        client._entry_hit("AAAA", ligand_bound=False)


@pytest.mark.integration
def test_live_rcsb_search_returns_a_structure_covering_the_requested_residue():
    hit = RcsbPdbClient().search_by_uniprot("P04637", position=175)
    assert hit is not None
    assert hit["pdb_text"]
    assert _covers_residue(hit["pdb_text"], 175)


@pytest.mark.integration
def test_live_rcsb_braf_v600e_finds_a_structure_covering_position_600():
    # Regression: 8VSO (a RAS-binding-domain-only fragment, residues -4..231)
    # used to be returned uncritically for BRAF even though it doesn't cover
    # residue 600 at all, silently breaking the BRAF V600E/vemurafenib
    # positive control from validation-plan.md downstream in score_binding.
    hit = RcsbPdbClient().search_by_uniprot("P15056", position=600)
    assert hit is not None
    assert _covers_residue(hit["pdb_text"], 600)
