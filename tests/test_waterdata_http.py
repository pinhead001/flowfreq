"""Water Data API key and 429/503 backoff (flowfreq.waterdata.request)."""

from __future__ import annotations

from typing import Any, Dict, List
from unittest.mock import Mock

import pytest
import requests

from flowfreq import waterdata


def _response(status: int, retry_after: str = "") -> Mock:
    r = Mock()
    r.status_code = status
    r.headers = {"Retry-After": retry_after} if retry_after else {}
    if status >= 400:
        r.raise_for_status.side_effect = requests.HTTPError(f"{status} error")
    else:
        r.raise_for_status.return_value = None
    return r


@pytest.fixture
def calls(monkeypatch):
    """Record every requests.get call and every sleep; clear any configured key."""
    seen: Dict[str, List[Any]] = {"get": [], "sleep": []}
    monkeypatch.delenv(waterdata.API_KEY_ENV, raising=False)
    waterdata.set_api_key(None)
    monkeypatch.setattr("time.sleep", lambda s: seen["sleep"].append(s))
    yield seen
    waterdata.set_api_key(None)


def _serve(monkeypatch, calls, responses):
    it = iter(responses)

    def fake_get(url, params=None, timeout=None, **kwargs):
        calls["get"].append(kwargs)
        return next(it)

    monkeypatch.setattr(waterdata.requests, "get", fake_get)


def test_no_key_sends_no_header(monkeypatch, calls):
    _serve(monkeypatch, calls, [_response(200)])
    waterdata.request("u", {}, 5)
    assert calls["get"] == [{}]


def test_key_from_environment_is_sent(monkeypatch, calls):
    monkeypatch.setenv(waterdata.API_KEY_ENV, " abc123 ")
    _serve(monkeypatch, calls, [_response(200)])
    waterdata.request("u", {}, 5)
    assert calls["get"] == [{"headers": {"X-Api-Key": "abc123"}}]


def test_set_api_key_overrides_the_environment(monkeypatch, calls):
    monkeypatch.setenv(waterdata.API_KEY_ENV, "from-env")
    waterdata.set_api_key("explicit")
    _serve(monkeypatch, calls, [_response(200)])
    waterdata.request("u", {}, 5)
    assert calls["get"][0]["headers"]["X-Api-Key"] == "explicit"


def test_429_is_retried_honouring_retry_after(monkeypatch, calls):
    _serve(monkeypatch, calls, [_response(429, "7"), _response(503), _response(200)])
    r = waterdata.request("u", {}, 5)
    assert r.status_code == 200
    assert calls["sleep"] == [7.0, 4.0]  # Retry-After, then 2**2 backoff


def test_retry_after_is_capped(monkeypatch, calls):
    _serve(monkeypatch, calls, [_response(429, "3600"), _response(200)])
    waterdata.request("u", {}, 5)
    assert calls["sleep"] == [waterdata.MAX_BACKOFF_S]


def test_persistent_429_names_the_api_key(monkeypatch, calls):
    _serve(monkeypatch, calls, [_response(429)] * (waterdata.MAX_RETRIES + 1))
    with pytest.raises(requests.HTTPError, match=waterdata.API_KEY_ENV):
        waterdata.request("u", {}, 5)
    assert len(calls["sleep"]) == waterdata.MAX_RETRIES


def test_rejected_key_is_not_retried(monkeypatch, calls):
    waterdata.set_api_key("bogus")
    _serve(monkeypatch, calls, [_response(403)])
    with pytest.raises(requests.HTTPError, match="rejected the API key"):
        waterdata.request("u", {}, 5)
    assert calls["sleep"] == []


def test_other_errors_are_not_retried(monkeypatch, calls):
    _serve(monkeypatch, calls, [_response(404)])
    with pytest.raises(requests.HTTPError, match="404"):
        waterdata.request("u", {}, 5)
    assert calls["sleep"] == []


def test_peak_backend_goes_through_the_shared_request(monkeypatch, calls):
    """The peaks backend gets the key and backoff too, not its own bare GET."""
    from flowfreq.peak_sources import WaterDataApiBackend

    waterdata.set_api_key("k")
    ok = _response(200)
    ok.json.return_value = {"features": [], "links": []}
    _serve(monkeypatch, calls, [_response(429), ok])
    payload = WaterDataApiBackend()._get_page("u", {}, "03606500")
    assert payload == {"features": [], "links": []}
    assert calls["get"][1]["headers"]["X-Api-Key"] == "k"
    assert len(calls["sleep"]) == 1
