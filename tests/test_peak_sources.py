"""Tests for flowfreq.peak_sources."""

from __future__ import annotations

import copy
import json
import logging
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest
import requests

from flowfreq.peak_codes import PEAK_CODE_DESCRIPTIONS, PeakTreatment, parse_codes, peak_interval
from flowfreq.peak_sources import (
    DEFAULT_BACKEND,
    PEAK_COLUMNS,
    QUALIFIER_TOKEN_CODES,
    WATERDATA_PEAKS_URL,
    LegacyNwisBackend,
    PeakDataBackend,
    WaterDataApiBackend,
    get_backend,
    qualifiers_to_codes,
    validate_peak_frame,
)


def _frame(**overrides):
    base = {
        "water_year": [2000, 2001],
        "peak_date": pd.to_datetime(["2000-03-01", "2001-04-02"]),
        "peak_flow_cfs": [100.0, 200.0],
        "qualification_code": ["", "2"],
        "extra": [1, 2],
    }
    base.update(overrides)
    return pd.DataFrame(base)


def test_validate_restricts_and_orders_columns():
    assert tuple(validate_peak_frame(_frame()).columns) == PEAK_COLUMNS


@pytest.mark.parametrize(
    "overrides, match",
    [
        ({"water_year": [2000, 2000]}, "Duplicate"),
        ({"peak_flow_cfs": [-1.0, 2.0]}, "Negative"),
    ],
)
def test_validate_rejects(overrides, match):
    with pytest.raises(ValueError, match=match):
        validate_peak_frame(_frame(**overrides))


def test_validate_missing_column():
    with pytest.raises(ValueError, match="missing"):
        validate_peak_frame(_frame().drop(columns="qualification_code"))


def test_default_is_waterdata_and_protocol_holds():
    b = get_backend()
    assert b.name == DEFAULT_BACKEND == WaterDataApiBackend.name
    assert isinstance(b, PeakDataBackend)


def test_legacy_backend_reads_the_rdb_parser(monkeypatch):
    seen = []

    def _rdb(site_no, timeout=30):
        seen.append(site_no)
        return _frame(), "NAME", 1.0

    monkeypatch.setattr("flowfreq.usgs._download_peak_flow_rdb", _rdb)
    out = LegacyNwisBackend().fetch_peaks("1234567")
    assert list(out["peak_flow_cfs"]) == [100.0, 200.0]
    assert tuple(out.columns) == PEAK_COLUMNS
    assert seen == ["01234567"]  # zero-padded, as USGSgage pads it


def test_legacy_backend_does_not_recurse_through_download_peak_flow(monkeypatch):
    """download_peak_flow routes through this module, so the backend must not call it."""
    from tests.fixtures.nwis_rdb import PEAK_PARTIAL_DATES

    def _boom(self, backend="x"):
        raise AssertionError("LegacyNwisBackend called USGSgage.download_peak_flow")

    response = MagicMock()
    response.text = PEAK_PARTIAL_DATES
    monkeypatch.setattr("flowfreq.usgs.USGSgage.download_peak_flow", _boom)
    monkeypatch.setattr("flowfreq.usgs.requests.get", lambda *a, **k: response)
    out = LegacyNwisBackend().fetch_peaks("03606500")
    assert list(out["water_year"]) == [1897, 1919, 1927, 1930, 1931, 1932]


def test_legacy_backend_rejects_an_empty_response(monkeypatch):
    response = MagicMock()
    response.text = "# only comments\n"
    monkeypatch.setattr("flowfreq.usgs.requests.get", lambda *a, **k: response)
    with pytest.raises(ValueError, match="No peak flow data"):
        LegacyNwisBackend().fetch_peaks("03606500")


# --- Water Data OGC API backend ----------------------------------------------

FIXTURE = Path(__file__).parent / "fixtures" / "waterdata_peaks_03606500.json"


@pytest.fixture
def pages():
    """Two real (trimmed) pages for Big Sandy: page1 has rel=next, page2 does not."""
    data = json.loads(FIXTURE.read_text(encoding="utf-8"))
    return data["page1"], data["page2"]


class _FakeGet:
    """Stand-in for ``requests.get`` that serves payloads in order and records calls."""

    def __init__(self, *payloads, status=200):
        self.payloads = list(payloads)
        self.status = status
        self.calls: list[tuple[str, dict | None]] = []

    def __call__(self, url, params=None, timeout=None):
        assert timeout is not None, "requests must carry a timeout"
        self.calls.append((url, params))
        resp = MagicMock()
        resp.json.return_value = self.payloads.pop(0)
        if self.status >= 400:
            resp.raise_for_status.side_effect = requests.HTTPError(f"{self.status} Error")
        return resp


def _fetch(monkeypatch, *payloads, site="03606500", status=200):
    fake = _FakeGet(*payloads, status=status)
    monkeypatch.setattr("flowfreq.peak_sources.requests.get", fake)
    return WaterDataApiBackend().fetch_peaks(site), fake


def test_waterdata_follows_next_link_across_pages(monkeypatch, pages):
    df, fake = _fetch(monkeypatch, *pages)
    assert len(fake.calls) == 2
    next_href = next(lk["href"] for lk in pages[0]["links"] if lk["rel"] == "next")
    # The second request is the next link verbatim; it already carries the query.
    assert fake.calls[1] == (next_href, None)
    assert list(df["water_year"]) == [1897, 1919, 1927, 1930, 1934, 2010, 2020]
    assert tuple(df.columns) == PEAK_COLUMNS
    validate_peak_frame(df)


def test_waterdata_sends_prefix_and_discharge_filter(monkeypatch, pages):
    _, fake = _fetch(monkeypatch, *pages)
    url, params = fake.calls[0]
    assert url == WATERDATA_PEAKS_URL
    assert params["monitoring_location_id"] == "USGS-03606500"
    assert params["parameter_code"] == "00060"
    assert params["f"] == "json"
    assert int(params["limit"]) > 10  # the server default of 10 would page needlessly


def test_waterdata_does_not_double_the_prefix(monkeypatch, pages):
    _, fake = _fetch(monkeypatch, *pages, site="USGS-03606500")
    assert fake.calls[0][1]["monitoring_location_id"] == "USGS-03606500"


def test_waterdata_values_dates_and_codes(monkeypatch, pages):
    df, _ = _fetch(monkeypatch, *pages)
    row = df.set_index("water_year")
    assert row.loc[1897, "peak_flow_cfs"] == 25000.0
    assert row.loc[1919, "peak_flow_cfs"] == 21000.0
    assert row.loc[1927, "peak_flow_cfs"] == 18500.0
    for wy in (1897, 1919, 1927):
        # DAYUNKNOWN,HISTORIC -> "7" only. Never the characters O/C of the tokens.
        assert row.loc[wy, "qualification_code"] == "7"
        interval = peak_interval(row.loc[wy, "peak_flow_cfs"], row.loc[wy, "qualification_code"])
        assert interval.is_historic and interval.treatment == PeakTreatment.SYSTEMATIC
    assert row.loc[2010, "qualification_code"] == "1"
    assert row.loc[1930, "qualification_code"] == ""
    # Placeholder date for an unknown day, as the legacy parser gives.
    assert row.loc[1927, "peak_date"] == pd.Timestamp("1926-12-01")


def test_waterdata_water_year_taken_from_its_own_field(monkeypatch, pages):
    """The API's water_year wins even where the date would say otherwise."""
    page1, page2 = copy.deepcopy(pages)
    page2["features"][0]["properties"]["water_year"] = 1901  # time stays 1926-12-01
    df, _ = _fetch(monkeypatch, page1, page2)
    assert 1901 in set(df["water_year"])
    assert 1927 not in set(df["water_year"])
    # And a Dec-1933 peak is WY 1934 because the API says so.
    assert df.set_index("water_year").loc[1934, "peak_date"] == pd.Timestamp("1933-12-18")


def test_waterdata_output_is_sorted(monkeypatch, pages):
    df, _ = _fetch(monkeypatch, *pages)
    assert df["water_year"].is_monotonic_increasing


def test_waterdata_empty_result_raises(monkeypatch):
    empty = {"type": "FeatureCollection", "numberReturned": 0, "features": [], "links": []}
    with pytest.raises(ValueError, match="No discharge peaks.*03606500"):
        _fetch(monkeypatch, empty)


def test_waterdata_http_error_names_site_and_url(monkeypatch, pages):
    with pytest.raises(requests.RequestException, match="03606500.*api.waterdata.usgs.gov"):
        _fetch(monkeypatch, pages[0], status=400)


def test_waterdata_skips_gage_height_rows(monkeypatch, pages, caplog):
    page1, page2 = copy.deepcopy(pages)
    gh = copy.deepcopy(page2["features"][0])
    gh["properties"].update(parameter_code="00065", unit_of_measure="ft", value="30.1")
    page2["features"].append(gh)
    with caplog.at_level(logging.WARNING, logger="flowfreq.peak_sources"):
        df, _ = _fetch(monkeypatch, page1, page2)
    assert len(df) == 7  # no duplicate WY 1927
    assert "00065" in caplog.text


@pytest.mark.parametrize(
    "tokens, expected",
    [
        (None, ""),
        ([], ""),
        (["DAYUNKNOWN", "HISTORIC"], "7"),
        (["MONTHUNKNOWN", "DAYUNKNOWN", "HISTORIC"], "7"),
        (["MAXDAILYMEAN"], "1"),
        (["ESTIMATED", "HISTORIC"], "2,7"),
        (["LESSTHAN"], "4"),
        (["GREATERTHAN"], "8"),
        (["DAMFAILURE"], "3"),
        (["UNKNOWNREGULATION"], "5"),
        (["REGULATED", "URBAN"], "6,C"),
        (["EVENT"], "9"),
        (["OPPORTUNISTIC"], "O"),
        (["REVISED"], "R"),
        (["GHNOTASSCPKQ", "ESTIMATED"], "2"),
        ("HISTORIC", "7"),
    ],
)
def test_qualifiers_to_codes(tokens, expected):
    out = qualifiers_to_codes(tokens)
    assert out == expected
    assert parse_codes(out) == frozenset(expected.split(",")) - {""}


def test_every_token_maps_to_a_known_peak_code():
    assert set(QUALIFIER_TOKEN_CODES.values()) <= set(PEAK_CODE_DESCRIPTIONS)


def test_censoring_codes_survive_translation():
    assert peak_interval(100.0, qualifiers_to_codes(["LESSTHAN"])).treatment == (
        PeakTreatment.LESS_THAN
    )
    assert peak_interval(100.0, qualifiers_to_codes(["GREATERTHAN"])).treatment == (
        PeakTreatment.GREATER_THAN
    )
    # A historic peak must not be removed by the O and C in "DAYUNKNOWN".
    hist = peak_interval(100.0, qualifiers_to_codes(["DAYUNKNOWN", "HISTORIC"]))
    assert hist.treatment == PeakTreatment.SYSTEMATIC and hist.is_historic


def test_unknown_token_warns_and_is_never_split(caplog):
    with caplog.at_level(logging.WARNING, logger="flowfreq.peak_sources"):
        out = qualifiers_to_codes(["NEWFLAGXYZ", "HISTORIC"], site_no="03606500", water_year=1897)
    assert out == "7"  # not "N,E,W,..." and not the token itself
    assert "NEWFLAGXYZ" in caplog.text
    assert "03606500" in caplog.text
    assert any(r.levelno == logging.WARNING for r in caplog.records)


def test_gage_height_token_is_not_a_warning(caplog):
    with caplog.at_level(logging.DEBUG, logger="flowfreq.peak_sources"):
        assert qualifiers_to_codes(["GHNOTASSCPKQ"]) == ""
    assert not any(r.levelno >= logging.WARNING for r in caplog.records)
    assert "GHNOTASSCPKQ" in caplog.text


def test_waterdata_backend_is_registered():
    assert isinstance(get_backend("waterdata-ogc"), WaterDataApiBackend)
    assert isinstance(get_backend("nwis-legacy"), LegacyNwisBackend)


@pytest.mark.requires_network
class TestLiveWaterData:
    """Live Water Data API checks. Deselect with -m 'not requires_network'."""

    def test_big_sandy_matches_fixture(self) -> None:
        from tests.fixtures.big_sandy import HISTORICAL_PEAKS, SYSTEMATIC_PEAKS

        df = WaterDataApiBackend().fetch_peaks("03606500")
        row = df.set_index("water_year")
        for wy, q in {**SYSTEMATIC_PEAKS, **HISTORICAL_PEAKS}.items():
            assert row.loc[wy, "peak_flow_cfs"] == pytest.approx(q), wy
        for wy in HISTORICAL_PEAKS:
            assert row.loc[wy, "qualification_code"] == "7", wy
        assert df["water_year"].is_monotonic_increasing
        assert not df["water_year"].duplicated().any()


def test_unknown_backend():
    with pytest.raises(KeyError):
        get_backend("nope")
