"""Tests for flowfreq.catalog and the generated coverage report."""

from __future__ import annotations

import importlib.util
import json

import pandas as pd
import pytest

from flowfreq.catalog import (
    CATALOG_PATH,
    COLUMNS,
    MIN_YEARS,
    SEED_PATH,
    build_rows,
    catalog_row,
    load_catalog,
    load_catalog_meta,
    location_fields,
    peaks_by_site,
    select,
    sites_needing_location,
    write_catalog,
)
from tests.fixtures.paths import REPO_ROOT


def test_seed_loads_with_all_columns():
    df = load_catalog(SEED_PATH)
    assert tuple(df.columns) == COLUMNS and len(df) == 3
    assert df["site_no"].str.len().min() >= 8


def _write(tmp_path, body, header="site_no,site_name,state,drainage_area_sqmi\n"):
    p = tmp_path / "c.csv"
    p.write_text(header + body)
    return p


@pytest.mark.parametrize(
    "body, match",
    [
        ("12345678,A,WA,10\n12345678,B,WA,10\n", "duplicate"),
        ("1234,A,WA,10\n", "digits"),
        ("12345678,A,ZZ,10\n", "unknown state"),
        ("12345678,A,WA,0\n", "positive"),
    ],
)
def test_validation(tmp_path, body, match):
    with pytest.raises(ValueError, match=match):
        load_catalog(_write(tmp_path, body))


def test_missing_required_column(tmp_path):
    with pytest.raises(ValueError, match="required"):
        load_catalog(_write(tmp_path, "12345678,A\n", header="site_no,site_name\n"))


def test_catalog_row_and_select():
    peaks = pd.DataFrame(
        {
            "water_year": [1990, 1991, 1995],
            "peak_flow_cfs": [1.0, 2.0, 3.0],
            "qualification_code": ["", "6", ""],
        }
    )
    row = catalog_row("12345678", "wa", peaks, site_name="Test", drainage_area_sqmi=5.0)
    assert (row["n_peaks"], row["first_water_year"], row["last_water_year"]) == (3, 1990, 1995)
    assert row["regulation_class"] == "regulated" and row["state"] == "WA"
    df = pd.DataFrame(
        [row, {**row, "site_no": "87654321", "regulation_class": "unknown", "n_peaks": 30}]
    )
    assert len(select(df, state="WA", exclude_regulated=True)) == 1
    assert len(select(df, min_years=10)) == 1


def test_catalog_row_counts_years_with_discharge_only():
    peaks = pd.DataFrame(
        {
            "water_year": [1990, 1990, 1991, 1992],
            "peak_flow_cfs": [1.0, 2.0, None, 3.0],  # 1991: gage height only
            "qualification_code": ["", "", "", ""],
        }
    )
    row = catalog_row("12048000", "WA", peaks)
    assert (row["n_peaks"], row["first_water_year"], row["last_water_year"]) == (2, 1990, 1992)
    assert row["regulation_class"] == "reference"  # GAGES-II Ref, from the #32 screen


# ------------------------------------------------- bulk build (offline) -----

CAPTURE = REPO_ROOT / "tests" / "fixtures" / "catalog_wa_capture.json"


@pytest.fixture(scope="module")
def capture():
    return json.loads(CAPTURE.read_text(encoding="utf-8"))


def test_bulk_rows_from_capture(capture):
    peaks = peaks_by_site(capture["peaks"])
    assert set(peaks) == {"12048000", "12045500", "12027000"}
    locs = {f["site_no"]: f for f in map(location_fields, capture["locations"])}
    assert locs["12048000"]["state"] == "WA" and locs["12048000"]["huc8"] == "17110020"
    assert sites_needing_location(peaks, {}) == ["12045500", "12048000"]
    assert sites_needing_location(peaks, locs) == []
    df = build_rows(peaks, locs, regions={"12048000": "GC1752"})
    assert list(df["site_no"]) == ["12045500", "12048000"]  # 12027000 has < 10 years
    by = df.set_index("site_no")
    assert by.loc["12048000", "regulation_class"] == "reference"
    assert by.loc["12045500", "regulation_class"] == "regulated"  # code-6 peaks
    assert by.loc["12048000", "regression_region"] == "GC1752"
    assert by.loc["12048000", "n_peaks"] >= 90
    assert (df["n_peaks"] >= MIN_YEARS).all()


def test_write_catalog_is_deterministic_and_valid(tmp_path, capture):
    peaks = peaks_by_site(capture["peaks"])
    locs = {f["site_no"]: f for f in map(location_fields, capture["locations"])}
    df = build_rows(peaks, locs)
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    # Same file name: gzip stores it in the header.
    a, b = tmp_path / "a" / "c.csv.gz", tmp_path / "b" / "c.csv.gz"
    write_catalog(df, a)
    write_catalog(df, b)
    assert a.read_bytes() == b.read_bytes()  # no gzip timestamp
    back = load_catalog(a)
    assert list(back["site_no"]) == list(df["site_no"]) and back["huc8"].str.len().eq(8).all()


def test_catalog_rejects_bad_class_and_huc(tmp_path):
    with pytest.raises(ValueError, match="regulation_class"):
        load_catalog(
            _write(
                tmp_path, "12345678,A,WA,altered\n", "site_no,site_name,state,regulation_class\n"
            )
        )
    with pytest.raises(ValueError, match="huc8"):
        load_catalog(_write(tmp_path, "12345678,A,WA,1701\n", "site_no,site_name,state,huc8\n"))


def _tool():
    spec = importlib.util.spec_from_file_location(
        "build_gage_catalog", REPO_ROOT / "tools" / "build_gage_catalog.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_assign_regions_point_in_polygon():
    tool = _tool()
    square = [[[-1, -1], [1, -1], [1, 1], [-1, 1], [-1, -1]]]
    polys = [
        {"gridcode": "gc1750", "rings": square},
        {"gridcode": "gc10001,gc730,gc731", "rings": square},  # helper code dropped
        {"gridcode": "gc0", "rings": square},  # National Urban: not a region
    ]
    sites = {
        "in": {"latitude": 0.0, "longitude": 0.0},
        "out": {"latitude": 5.0, "longitude": 5.0},
        "nopoint": {"latitude": None, "longitude": None},
    }
    assert tool.assign_regions(sites, polys) == {"in": "GC1750;GC730;GC731"}


def test_nss_region_rings_feed_assign_regions():
    """NSS GeoJSON (Multi)Polygons -> rounded rings, holes and parts counted even-odd."""
    tool = _tool()
    outer = [[0, 0], [4, 0], [4, 4], [0, 4], [0, 0]]
    hole = [[1, 1], [2, 1], [2, 2], [1, 2], [1, 1]]
    part = [[10.123456789, 10], [11, 10], [11, 11], [10, 11], [10.123456789, 10]]
    multi = {"type": "MultiPolygon", "coordinates": [[outer, hole], [part]]}
    rings = tool.nss_region_rings(multi)
    assert len(rings) == 3 and rings[2][0] == [10.12346, 10.0]
    assert tool.nss_region_rings({"type": "Point", "coordinates": [0, 0]}) == []
    whole = {"type": "Polygon", "coordinates": [[[-1, -1], [20, -1], [20, 20], [-1, 20], [-1, -1]]]}
    polys = [
        {"gridcode": "GC1621", "rings": rings},
        {"gridcode": "GC1618", "rings": tool.nss_region_rings(whole)},  # statewide, like AZ 1
    ]
    sites = {
        "a": {"latitude": 0.5, "longitude": 0.5},
        "hole": {"latitude": 1.5, "longitude": 1.5},
        "b": {"latitude": 10.5, "longitude": 10.5},
    }
    assert tool.assign_regions(sites, polys) == {
        "a": "GC1618;GC1621",
        "hole": "GC1618",
        "b": "GC1618;GC1621",
    }


def test_nss_region_polygons_uses_cache(tmp_path, monkeypatch):
    """A cached region is never refetched (resumable builds); one is fetched otherwise."""
    tool = _tool()
    monkeypatch.setattr(tool, "NSS_REGION_IDS", {"ZZ": (1, 2)})
    cached = {"gridcode": "GC1", "rings": [[[0, 0], [1, 0], [1, 1], [0, 0]]]}
    (tmp_path / "nss_regions").mkdir()
    (tmp_path / "nss_regions" / "1.json").write_text(json.dumps(cached))

    class Resp:
        def raise_for_status(self):
            pass

        def json(self):
            return {
                "code": "GC2",
                "location": {"geometry": {"type": "Polygon", "coordinates": cached["rings"]}},
            }

    calls = []

    def fake_get(url, params, timeout):
        calls.append(url)
        return Resp()

    monkeypatch.setattr(tool.requests, "get", fake_get)
    got = tool.nss_region_polygons("ZZ", tmp_path, refresh=False)
    assert [g["gridcode"] for g in got] == ["GC1", "GC2"]
    assert calls == [f"{tool.NSS_REGRESSION_REGIONS}/2"]


# ------------------------------------------------- packaged catalog ---------


@pytest.fixture(scope="module")
def packaged():
    return load_catalog(CATALOG_PATH)


def test_packaged_catalog(packaged):
    meta = load_catalog_meta()
    assert len(packaged) == meta["rows"]
    assert (packaged["n_peaks"] >= meta["min_years"]).all()
    assert packaged["regulation_class"].value_counts().to_dict() == meta["regulation_class_counts"]
    located = packaged["state"].isin(["WA", "OR", "ID", "MT", "CO", "UT", "AZ"])
    filled = packaged["regression_region"].notna()
    assert not (filled & ~located).any()  # Wave 1, and Wave 2's CO, UT and AZ, only
    for state, share in {
        "WA": 0.8,
        "OR": 0.8,
        "MT": 0.8,
        "CO": 0.99,
        "UT": 0.99,
        "AZ": 0.99,
    }.items():
        assert filled[packaged["state"] == state].mean() > share, state
    assert packaged["regression_region"].dropna().str.fullmatch(r"GC\d+(;GC\d+)*").all()
    assert meta["regression_region_filled"] == int(filled.sum())


def test_packaged_catalog_wave2_regions(packaged):
    """CO and UT regions are disjoint; AZ lists High Elevation (GC1618) as a candidate."""
    reg = packaged.dropna(subset=["regression_region"]).set_index("site_no")
    co = reg[reg["state"] == "CO"]["regression_region"]
    assert set(co) <= {"GC1204", "GC1205", "GC1206", "GC1207", "GC1208", "GC1698"}
    ut = reg[reg["state"] == "UT"]["regression_region"]
    assert set(ut) <= {f"GC{n}" for n in range(960, 967)}
    az = reg[reg["state"] == "AZ"]["regression_region"]
    assert set(az) == {"GC1618;GC1619", "GC1618;GC1620", "GC1618;GC1621", "GC1618;GC1622"}
    # Known gages: Crystal River above Avalanche Creek (Roaring Fork basin, CO
    # Northwest region) and East Canyon Creek near Morgan, UT. When built, 36 of 36
    # sampled CO/UT/AZ rows matched NSS's own bylocation answer for the gage point.
    assert reg.loc["09081600", "regression_region"] == "GC1205"
    assert reg.loc["10134500", "regression_region"].startswith("GC96")
    # No verified point-in-region source for WY, NV and NM: blank.
    assert (
        not packaged[packaged["state"].isin(["WY", "NV", "NM"])]["regression_region"].notna().any()
    )


def test_packaged_catalog_known_sites(packaged):
    by = packaged.set_index("site_no")
    dung = by.loc["12048000"]  # Dungeness River near Sequim, WA
    assert dung["state"] == "WA" and dung["regulation_class"] == "reference"
    assert dung["huc8"] == "17110020" and dung["regression_region"].startswith("GC175")
    assert dung["first_water_year"] <= 1924 and dung["n_peaks"] >= 90
    assert by.loc["03606500", "state"] == "TN"  # Big Sandy at Bruceton: outside Wave 1
    assert pd.isna(by.loc["03606500", "regression_region"])


# ------------------------------------------------- live ---------------------


@pytest.mark.requires_network
def test_waterdata_bulk_filters_are_live():
    """The two bulk queries the tool depends on still accept their filters."""
    from flowfreq import waterdata

    base = waterdata.WATERDATA_BASE_URL
    peaks = waterdata.request(
        f"{base}/peaks/items",
        {
            "f": "json",
            "state_code": "53",
            "parameter_code": "00060",
            "properties": "monitoring_location_id,water_year,value,qualifier",
            "skipGeometry": "true",
            "limit": 5,
        },
        60,
    ).json()["features"]
    assert len(peaks) == 5
    assert set(peaks[0]["properties"]) == {
        "monitoring_location_id",
        "water_year",
        "value",
        "qualifier",
    }
    locs = waterdata.request(
        f"{base}/monitoring-locations/items",
        {"f": "json", "state_code": "53", "site_type_code": "ST", "limit": 2},
        60,
    ).json()["features"]
    f = location_fields(locs[0])
    assert f["state"] == "WA" and f["latitude"] is not None


@pytest.mark.requires_network
def test_nss_wave2_region_ids_are_live():
    """The NSS ids the tool fetches still name the peak-flow regions it expects."""
    import requests

    tool = _tool()
    for state, ids in tool.NSS_REGION_IDS.items():
        listed = requests.get(
            f"https://streamstats.usgs.gov/nssservices/regions/{state}/regressionregions",
            params={"statisticgroups": 2},
            timeout=60,
        ).json()
        assert set(ids) <= {r["id"] for r in listed}, state
    # Wyoming's peak regions still carry no geometry (why WY stays blank).
    wy = requests.get(
        f"{tool.NSS_REGRESSION_REGIONS}/746", params={"includeGeometry": "true"}, timeout=60
    )
    assert not ((wy.json().get("location") or {}).get("geometry"))


def test_regression_coverage_doc_is_current():
    """docs/REGRESSION_COVERAGE.md must match what the generator produces now."""
    spec = importlib.util.spec_from_file_location(
        "gen_cov", REPO_ROOT / "tools" / "gen_regression_coverage.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    committed = (REPO_ROOT / "docs" / "REGRESSION_COVERAGE.md").read_text(encoding="utf-8")
    assert committed == mod.render(), "stale: run python tools/gen_regression_coverage.py"
