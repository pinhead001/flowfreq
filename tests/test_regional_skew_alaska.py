"""Alaska's regional skew areas (SIR 2016-5024, #204).

Table 5 defines Regional Skew Areas 1 and 2 by HUC8, or by HUC10 where a HUC8 is split;
everything else in HU 19 is the exclusion area, where only station skew applies. Table 6
limits each area to a drainage-area range, outside which station skew applies as well.
HUC12s below are built inside the listed units (lookup is by prefix).
"""

from __future__ import annotations

import pytest

from flowfreq.regional_skew import (
    DRAINAGE_AREA_RANGES,
    RegionalSkewUnavailable,
    regional_skew_at_huc,
)

RSA1, RSA2 = "Regional Skew Area 1", "Regional Skew Area 2"


@pytest.mark.parametrize(
    "huc, region, skew",
    [
        ("190803010101", RSA1, 0.54),  # whole HUC8 19080301
        ("190201040201", RSA1, 0.54),  # 1902010402, a listed HUC10 of split 19020104
        ("190203011001", RSA1, 0.54),  # 1902030110, the last Area 1 HUC10 of 19020301
        ("190101020101", RSA2, 0.18),  # whole HUC8 19010102
        ("190201040501", RSA2, 0.18),  # 1902010405, an Area 2 HUC10 of split 19020104
        ("190201041801", RSA2, 0.18),  # 1902010418, printed '19020104218' in Table 5
        ("190203011101", RSA2, 0.18),  # 1902030111, an Area 2 HUC10 of split 19020301
    ],
)
def test_areas_by_hydrologic_unit(huc: str, region: str, skew: float) -> None:
    rs = regional_skew_at_huc(huc, "AK", drainage_area_sqmi=50.0)
    assert rs.skew_region == region and rs.skew == skew
    assert not rs.advisory


@pytest.mark.parametrize(
    "huc",
    [
        "190501010101",  # HU 1905: outside both areas
        "190201041901",  # 1902010419 (Martin River): in neither area's HUC10 list
        "190101060401",  # 1901010604: split 19010106, not listed
    ],
)
def test_exclusion_area_uses_station_skew(huc: str) -> None:
    with pytest.raises(RegionalSkewUnavailable, match="Use station skew"):
        regional_skew_at_huc(huc, "AK")


def test_a_split_huc8_alone_is_ambiguous() -> None:
    with pytest.raises(RegionalSkewUnavailable, match="only partly covered"):
        regional_skew_at_huc("19020104", "AK")


def test_drainage_area_ranges_are_table_6() -> None:
    assert DRAINAGE_AREA_RANGES[RSA1] == (1.2, 25560.0)
    assert DRAINAGE_AREA_RANGES[RSA2] == (1.7, 123.0)


@pytest.mark.parametrize(
    "huc, area", [("190803010101", 0.5), ("190803010101", 30000.0), ("190101020101", 500.0)]
)
def test_outside_the_drainage_area_range_uses_station_skew(huc: str, area: float) -> None:
    with pytest.raises(RegionalSkewUnavailable, match="outside.*applicable.*station skew"):
        regional_skew_at_huc(huc, "AK", drainage_area_sqmi=area)


def test_without_a_drainage_area_the_range_is_an_advisory() -> None:
    rs = regional_skew_at_huc("190101020101", "AK")
    assert rs.skew_region == RSA2
    assert "1.7" in rs.advisory and "123" in rs.advisory and "station skew" in rs.advisory


def test_drainage_area_is_ignored_where_no_range_applies() -> None:
    # Henry Fork, NC: a state fallback study with no drainage-area limit.
    rs = regional_skew_at_huc("030501020000", "NC", drainage_area_sqmi=1e6)
    assert rs.skew_region == "Southeastern United States"


def test_for_site_uses_the_monitoring_location_drainage_area(monkeypatch) -> None:
    import flowfreq.waterdata as wd
    from flowfreq.regional_skew import regional_skew_for_site

    def fake(site, timeout=60):
        area = {"15000001": 50.0, "15000002": 900.0}[site]
        return {"hydrologic_unit_code": "190101020101", "state_code": "02", "drainage_area": area}

    monkeypatch.setattr(wd, "fetch_monitoring_location", fake)
    rs = regional_skew_for_site("15000001")
    assert rs.skew_region == RSA2 and not rs.advisory
    with pytest.raises(RegionalSkewUnavailable, match="outside"):
        regional_skew_for_site("15000002")
    # An explicit drainage area overrides the record's.
    assert regional_skew_for_site("15000002", drainage_area_sqmi=40.0).skew_region == RSA2


def test_at_passes_the_drainage_area(monkeypatch) -> None:
    import flowfreq.regional_skew as rs_mod

    monkeypatch.setattr(rs_mod, "huc12_at", lambda lat, lon, timeout=60: ("190101020101", ["AK"]))
    with pytest.raises(RegionalSkewUnavailable, match="outside"):
        rs_mod.regional_skew_at(55.4, -131.2, drainage_area_sqmi=900.0)
    assert rs_mod.regional_skew_at(55.4, -131.2, drainage_area_sqmi=40.0).skew == 0.18


def test_a_nan_drainage_area_is_treated_as_unknown() -> None:
    rs = regional_skew_at_huc("190101020101", "AK", drainage_area_sqmi=float("nan"))
    assert rs.skew_region == RSA2 and "123" in rs.advisory


def test_for_site_ignores_a_non_numeric_drainage_area(monkeypatch) -> None:
    import flowfreq.waterdata as wd
    from flowfreq.regional_skew import regional_skew_for_site

    monkeypatch.setattr(
        wd,
        "fetch_monitoring_location",
        lambda site, timeout=60: {
            "hydrologic_unit_code": "190101020101", "state_code": "02", "drainage_area": "n/a"
        },
    )  # fmt: skip
    assert "123" in regional_skew_for_site("15000003").advisory
