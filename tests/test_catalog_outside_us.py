"""Catalog build: a located site with no US state is skipped, not filed under the queried
state (#206).

Gage 15388960 (PORCUPINE R NR INTERNATIONAL BOUNDARY YT, in the Yukon) came back from the
DC peak query. Checked against the Water Data API on 2026-10-30, its location record carries
``state_code = "11"`` (Yukon's province code under ``country_code = "CA"``; 11 is also DC's
FIPS code), so ``location_fields`` must read ``country_code`` and give a non-US site no state.
``build_rows`` then skips a located site with no state instead of falling back to
``default_state``, which is for a site with *no* location record at all.
"""

from __future__ import annotations

import logging

import pandas as pd
import pytest

from flowfreq.catalog import build_rows, location_fields


def _peaks(n: int = 12) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "water_year": list(range(2000, 2000 + n)),
            "peak_flow_cfs": [1000.0 + 10 * i for i in range(n)],
            "qualification_code": [""] * n,
        }
    )


def _loc(
    site: str, state_code: str | None, lon: float, lat: float, country: str | None = None
) -> dict:
    props = {"monitoring_location_name": f"SITE {site}", "drainage_area": 100.0}
    if state_code is not None:
        props["state_code"] = state_code
    if country is not None:
        props["country_code"] = country
    return location_fields(
        {"id": f"USGS-{site}", "properties": props,
         "geometry": {"type": "Point", "coordinates": [lon, lat]}}
    )  # fmt: skip


def test_location_without_state_reads_as_none() -> None:
    assert _loc("15388960", None, -140.891111, 67.424167)["state"] is None
    assert _loc("15388960", "", -140.891111, 67.424167)["state"] is None


def test_foreign_country_reads_as_none_even_with_a_province_code() -> None:
    """The real 15388960 record: province code 11 under country CA, not a US state."""
    assert _loc("15388960", "11", -140.891111, 67.424167, country="CA")["state"] is None


def test_us_country_code_keeps_the_state() -> None:
    assert _loc("01648000", "11", -77.04, 38.9725, country="US")["state"] == "DC"
    assert _loc("01648000", "11", -77.04, 38.9725)["state"] == "DC"  # no country: US


def test_foreign_located_site_is_skipped_by_build_rows() -> None:
    peaks = {"01648000": _peaks(), "15388960": _peaks()}
    locs = {
        "01648000": _loc("01648000", "11", -77.04, 38.9725, country="US"),
        "15388960": _loc("15388960", "11", -140.891111, 67.424167, country="CA"),
    }
    df = build_rows(peaks, locs, default_state="DC")
    assert list(df["site_no"]) == ["01648000"]


def test_located_site_without_state_is_skipped_and_logged(
    caplog: pytest.LogCaptureFixture,
) -> None:
    peaks = {"01648000": _peaks(), "15388960": _peaks()}
    locs = {
        "01648000": _loc("01648000", "11", -77.04, 38.9725),
        "15388960": _loc("15388960", None, -140.891111, 67.424167),
    }
    with caplog.at_level(logging.WARNING, logger="flowfreq.catalog"):
        df = build_rows(peaks, locs, default_state="DC")
    assert list(df["site_no"]) == ["01648000"]
    assert list(df["state"]) == ["DC"]
    msgs = [r.getMessage() for r in caplog.records if r.name == "flowfreq.catalog"]
    assert any("15388960" in m for m in msgs), msgs


def test_site_without_location_record_still_takes_default_state() -> None:
    df = build_rows({"01648010": _peaks()}, {}, default_state="DC")
    assert list(df["site_no"]) == ["01648010"] and list(df["state"]) == ["DC"]


def test_site_without_location_or_default_is_dropped() -> None:
    assert build_rows({"01648010": _peaks()}, {}).empty


def test_located_site_keeps_its_own_state_over_default() -> None:
    locs = {"12048000": _loc("12048000", "53", -123.0, 47.0)}
    df = build_rows({"12048000": _peaks()}, locs, default_state="DC")
    assert list(df["state"]) == ["WA"]
