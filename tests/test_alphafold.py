import httpx
import pytest

from secondlook.alphafold import AlphaFoldDbClient, AlphaFoldError


class ConnectTimeoutClient:
    """Simulates a connection-level failure (no response object ever exists) —
    distinct from an HTTP status error, which already has a response to call
    raise_for_status() on."""

    def get(self, url, timeout=None):
        raise httpx.ConnectTimeout("timed out")


def test_connect_timeout_is_wrapped_not_raised_raw():
    client = AlphaFoldDbClient(client=ConnectTimeoutClient())
    with pytest.raises(AlphaFoldError):
        client.fetch_models("P04637")


@pytest.mark.integration
def test_live_alphafold_p04637_models():
    client = AlphaFoldDbClient()
    models = client.fetch_models("P04637")
    assert any(model["entry_id"] == "AF-P04637-F1" for model in models)
