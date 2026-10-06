"""Location lookup of the regional skew: ``regional_skew_at_huc``/``_at``/``_for_site``.

The 12-digit HUs below are the real ones of each named USGS gage, taken from
the Watershed Boundary Dataset point query at the gage's NWIS coordinates on
2026-09-27 (for 12340500 also equal to the Water Data API's own
``hydrologic_unit_code``). Offline tests use them directly; network calls are
mocked except in the ``requires_network`` tests.
"""

from __future__ import annotations

import contextlib
from typing import Any, Dict, List

import pandas as pd
import pytest
import requests

import flowfreq.regional_skew as rs_mod
from flowfreq.regional_skew import (
    HUC8_STATES_TABLE,
    RegionalSkewUnavailable,
    load_huc_table,
    regional_skew_at,
    regional_skew_at_huc,
    regional_skew_for_site,
)

PNW_SKEW, PNW_MSE = -0.07, 0.18

PNW_SITES = {
    # site: (HU12, state)
    "12048000": ("171100200306", "WA"),  # Dungeness River near Sequim, WA
    "14306500": ("171002050403", "OR"),  # Alsea River near Tidewater, OR
    "13336500": ("170603020403", "ID"),  # Selway River near Lowell, ID
    "12340500": ("170102040104", "MT"),  # Clark Fork above Missoula, MT
    "12354000": ("170102040712", "MT"),  # St. Regis River near St. Regis, MT
}


# ----------------------------------------------------------------------------
# Tables
# ----------------------------------------------------------------------------


def test_packaged_tables_load():
    hucs = load_huc_table()
    assert set(hucs["rule"]) == {"member", "excluded", "boundary", "unresolved"}
    members = hucs[hucs["rule"] == "member"]
    assert sorted(members["huc"]) == [f"17{n:02d}" for n in range(1, 12)]
    assert (members["valid_states"] == "ID;MT;OR;WA").all()
    assert sorted(hucs.loc[hucs["rule"] == "unresolved", "huc"]) == ["1712", "1801"]
    srp = hucs[hucs["skew_region"] == "Snake River Plain"]
    assert srp["huc"].str.len().eq(12).all() and srp["huc"].str[:4].isin(["1704", "1705"]).all()
    assert srp["source"].str.contains("SIR 2016-5083 p. 52").all()
    states = pd.read_csv(HUC8_STATES_TABLE, dtype=str)
    assert states["huc8"].str.len().eq(8).all()


@pytest.mark.parametrize(
    "mutate, match",
    [
        (lambda df: df.assign(rule="maybe"), "unknown rule"),
        (lambda df: df.drop(columns=["source"]), "missing columns"),
        (lambda df: pd.concat([df, df.iloc[:1]]), "duplicate HUC"),
        (lambda df: df.assign(huc="170"), "malformed"),
    ],
)
def test_huc_table_validation(tmp_path, mutate, match):
    bad = mutate(load_huc_table())
    path = tmp_path / "hucs.csv"
    bad.to_csv(path, index=False)
    with pytest.raises(ValueError, match=match):
        load_huc_table(path)


# ----------------------------------------------------------------------------
# regional_skew_at_huc -- offline
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("site", sorted(PNW_SITES))
def test_pnw_sites(site, caplog):
    huc, state = PNW_SITES[site]
    rs = regional_skew_at_huc(huc, state)
    assert (rs.skew_region, rs.skew, rs.skew_mse) == ("Pacific Northwest", PNW_SKEW, PNW_MSE)
    assert rs.effective_record_years == 41.0
    assert rs.state == state
    assert "sir2016" in rs.citation


def test_state_inferred_from_huc8():
    # Alsea HUC8 17100205 is Oregon-only.
    rs = regional_skew_at_huc(PNW_SITES["14306500"][0])
    assert rs.state == "OR" and rs.advisory == ""


def test_multistate_huc8_inside_the_study():
    # Lower Selway HUC8 touches ID and MT, both in the study: no state needed,
    # but MT has not adopted it, so the advisory is set.
    rs = regional_skew_at_huc(PNW_SITES["13336500"][0])
    assert rs.state == "ID,MT" and "MT" in rs.advisory


def test_huc8_touching_canada_needs_state():
    # Dungeness-Elwha HUC8 17110020 is tagged CN,WA in the WBD.
    with pytest.raises(RegionalSkewUnavailable, match="Pass state="):
        regional_skew_at_huc(PNW_SITES["12048000"][0])


def test_western_montana_carries_the_advisory(caplog):
    """USGS MT has not adopted the PNW study (MT row notes); the lookup says so."""
    rs = regional_skew_at_huc("170102040104", "MT")
    assert rs.skew == PNW_SKEW
    assert "MT has no verified row" in rs.advisory
    assert "no verified row" in caplog.text


def test_snake_river_plain_raises():
    # Snake River near Blackfoot, ID (13069500): an HU12 wholly on the Plain.
    with pytest.raises(RegionalSkewUnavailable, match="Snake River Plain.*p. 52"):
        regional_skew_at_huc("170402061001", "ID")
    with pytest.raises(RegionalSkewUnavailable, match="Snake River Plain"):
        regional_skew_at_huc("170402061001")


@pytest.mark.parametrize(
    "huc, point",
    [
        ("170402090301", (42.7675, -112.8794444)),  # Snake River at Neeley, ID (13077000)
        ("170402180803", (43.5822222, -113.2705556)),  # Big Lost River nr Arco, ID (13132500)
    ],
)
def test_snake_river_plain_boundary_gages_raise(monkeypatch, huc, point):
    # Both HU12s are only partly Plain; both gages are on it (live point
    # query, 2026-09-28).
    monkeypatch.setattr(rs_mod, "_on_snake_river_plain", lambda lat, lon, timeout=60: True)
    with pytest.raises(RegionalSkewUnavailable, match="Snake River Plain"):
        regional_skew_at_huc(huc, "ID", point=point)


def test_snake_river_plain_boundary_needs_a_point():
    # Medicine Lodge Creek near Small, ID (13116500): a PNW study gage in an
    # HU12 that is mostly Plain. Its code alone cannot decide.
    with pytest.raises(RegionalSkewUnavailable, match="regional_skew_at"):
        regional_skew_at_huc("170402150206", "ID")


@pytest.mark.parametrize("on_plain", [True, False])
def test_snake_river_plain_boundary_decided_by_point(monkeypatch, on_plain):
    calls: List[Any] = []

    def fake(lat: float, lon: float, *, timeout: int = 60) -> bool:
        calls.append((lat, lon))
        return on_plain

    monkeypatch.setattr(rs_mod, "_on_snake_river_plain", fake)
    point = (44.25916667, -112.4102778)
    if on_plain:
        with pytest.raises(RegionalSkewUnavailable, match="Snake River Plain"):
            regional_skew_at_huc("170402150206", "ID", point=point)
    else:
        assert regional_skew_at_huc("170402150206", "ID", point=point).skew == PNW_SKEW
    assert calls == [point]


@pytest.mark.parametrize(
    "huc, what",
    [
        ("100301021602", "Missouri River at Fort Benton, MT (06090800)"),
        ("100700041006", "Yellowstone River at Billings, MT (06214500)"),
        ("090400010404", "St. Mary River near Babb, MT (05017500; Hudson Bay)"),
    ],
)
def test_eastern_montana_raises(huc, what):
    with pytest.raises(RegionalSkewUnavailable, match="No verified regional skew study"):
        regional_skew_at_huc(huc, "MT")


def test_state_with_no_study_raises():
    # Big Sandy River at Bruceton, TN (03606500).
    with pytest.raises(RegionalSkewUnavailable, match="No verified regional skew study"):
        regional_skew_at_huc("060400050505", "TN")


@pytest.mark.parametrize(
    "huc, state, match",
    [
        ("171200030108", "OR", "1712.*unresolved"),  # Donner und Blitzen, OR (10396000)
        ("180102010605", "OR", "1801.*unresolved"),  # Williamson River, OR (11502500)
        ("170401030103", "WY", "WY is outside"),  # Cache Creek near Jackson, WY (13018300)
        ("170401030103", None, "WY"),  # Greys-Hoback HUC8 touches ID and WY
        ("17020017", None, "CN, outside"),  # Headwaters Columbia River, Canada only
        ("170402100403", None, "Pass state="),  # ID/UT border HU of study gage 13078000
    ],
)
def test_edges_raise(huc, state, match):
    with pytest.raises(RegionalSkewUnavailable, match=match):
        regional_skew_at_huc(huc, state)


def test_border_hu_resolves_with_state():
    assert regional_skew_at_huc("170402100403", "ID").skew == PNW_SKEW


def test_coarse_huc_partly_excluded_raises():
    # HUC8 17040206 (American Falls) is part Plain, part not.
    with pytest.raises(RegionalSkewUnavailable, match="only partly covered"):
        regional_skew_at_huc("17040206", "ID")


def test_coarse_huc_wholly_inside_resolves():
    assert regional_skew_at_huc("1709").skew_region == "Pacific Northwest"  # Willamette


@pytest.mark.parametrize("huc", ["17010", "abc", "", "1701020401041"])
def test_malformed_huc(huc):
    with pytest.raises(ValueError, match="2-12 digits"):
        regional_skew_at_huc(huc)


# ----------------------------------------------------------------------------
# regional_skew_at / regional_skew_for_site -- mocked network
# ----------------------------------------------------------------------------


class _Resp:
    def __init__(self, payload: Dict[str, Any]) -> None:
        self._payload = payload

    def raise_for_status(self) -> None:
        return None

    def json(self) -> Dict[str, Any]:
        return self._payload


def _wbd(huc: str, states: str):
    def get(url: str, params: Dict[str, Any], timeout: int) -> _Resp:
        assert url == rs_mod.WBD_HU12_QUERY_URL
        assert params["geometry"] == "-113.9321194,46.87676389"
        return _Resp({"features": [{"attributes": {"huc12": huc, "states": states}}]})

    return get


def test_regional_skew_at_mocked(monkeypatch):
    monkeypatch.setattr(rs_mod.requests, "get", _wbd("170102040104", "MT"))
    rs = regional_skew_at(46.87676389, -113.9321194)
    assert (rs.skew, rs.state) == (PNW_SKEW, "MT")


def test_regional_skew_at_no_hu(monkeypatch):
    monkeypatch.setattr(rs_mod.requests, "get", lambda *a, **k: _Resp({"features": []}))
    with pytest.raises(RegionalSkewUnavailable, match="No hydrologic unit"):
        regional_skew_at(46.87676389, -113.9321194)


def test_regional_skew_at_service_error(monkeypatch):
    monkeypatch.setattr(rs_mod.requests, "get", lambda *a, **k: _Resp({"error": {"code": 500}}))
    with pytest.raises(rs_mod.requests.RequestException, match="WBD HU12"):
        regional_skew_at(46.87676389, -113.9321194)


def test_regional_skew_at_bad_coordinates():
    with pytest.raises(ValueError, match="latitude"):
        regional_skew_at(146.0, -113.9)


def test_regional_skew_for_site_mocked(monkeypatch):
    import flowfreq.waterdata as wd

    monkeypatch.setattr(
        wd,
        "fetch_monitoring_location",
        lambda site, timeout=60: {"hydrologic_unit_code": "170603020403", "state_code": "16"},
    )
    rs = regional_skew_for_site("13336500")
    assert (rs.skew, rs.state) == (PNW_SKEW, "ID")


def test_regional_skew_for_site_without_huc(monkeypatch):
    import flowfreq.waterdata as wd

    monkeypatch.setattr(wd, "fetch_monitoring_location", lambda site, timeout=60: {})
    with pytest.raises(RegionalSkewUnavailable, match="no hydrologic_unit_code"):
        regional_skew_for_site("00000000")


# ----------------------------------------------------------------------------
# Live
# ----------------------------------------------------------------------------


@contextlib.contextmanager
def _skip_on_service_outage():
    """Skip, not fail, when a remote service is down (5xx, connection, timeout).

    The weekly live run failed on 2026-10-05 on a 504 Gateway Time-out from
    hydro.nationalmap.gov (WBD). An outage is not a regression in this code; a
    4xx or a wrong answer still fails.
    """
    try:
        yield
    except requests.HTTPError as exc:
        status = exc.response.status_code if exc.response is not None else 0
        if status >= 500:
            pytest.skip(f"remote service unavailable: {exc}")
        raise
    except (requests.ConnectionError, requests.Timeout) as exc:
        pytest.skip(f"remote service unreachable: {exc}")


@pytest.mark.requires_network
def test_live_point_lookups():
    with _skip_on_service_outage():
        assert regional_skew_at(46.87676389, -113.9321194).skew == PNW_SKEW  # 12340500
        assert regional_skew_at(44.25916667, -112.4102778).skew == PNW_SKEW  # 13116500
        with pytest.raises(RegionalSkewUnavailable, match="Snake River Plain"):
            regional_skew_at(43.12527778, -112.5188889)  # 13069500
        with pytest.raises(RegionalSkewUnavailable):
            regional_skew_at(45.800119, -108.468031)  # 06214500


@pytest.mark.requires_network
def test_live_site_lookup():
    with _skip_on_service_outage():
        assert regional_skew_for_site("12048000").state == "WA"
        with pytest.raises(RegionalSkewUnavailable):
            regional_skew_for_site("06214500")


def test_outage_helper_skips_5xx_but_not_4xx():
    resp = requests.Response()
    resp.status_code = 504
    with pytest.raises(pytest.skip.Exception):
        with _skip_on_service_outage():
            raise requests.HTTPError("504", response=resp)
    resp4 = requests.Response()
    resp4.status_code = 404
    with pytest.raises(requests.HTTPError):
        with _skip_on_service_outage():
            raise requests.HTTPError("404", response=resp4)
