"""Tests for flowfreq.regulation (issue #32): the GAGES-II screen and the B17C refusal."""

from __future__ import annotations

import gzip
import json
import logging

import numpy as np
import pandas as pd
import pytest
import requests

from flowfreq import regulation as reg
from flowfreq.regulation import (
    SCREEN_COLUMNS,
    RegulatedRecordError,
    SiteClass,
    classify_site,
    load_screen,
    require_unregulated,
)
from flowfreq.streamstats import WatershedCharacteristics
from flowfreq.workflow import run_ffa
from tests.fixtures.paths import REPO_ROOT

# Rules of tools/build_regulation_screen.py, restated so the packaged table is
# checked against them rather than trusted.
REGULATED_STORAGE_DAYS = 127.8  # Dudley and others (2018), 75th percentile
URBAN_IMPERVIOUS_PCT = 5.0  # Mastin and others (2016), SIR 2016-5118 p. 23

REFERENCE_SITE = "12048000"  # Dungeness River near Sequim, WA (GAGES-II Ref)
REGULATED_SITE = "01019000"  # GAGES-II Non-ref, ~415 days of dam storage
URBAN_SITE = "01069700"  # 5.6 percent impervious
UNKNOWN_SITE = "03606500"  # Big Sandy at Bruceton: Non-ref, below both thresholds
ABSENT_SITE = "99999999"


@pytest.fixture(scope="module")
def table() -> pd.DataFrame:
    return load_screen()


def test_packaged_table_shape(table):
    assert tuple(table.columns) == SCREEN_COLUMNS
    assert table["runave7100_mm"].notna().sum() > 9000
    assert len(table) == 9322  # every GAGES-II gage
    assert table["site_no"].str.fullmatch(r"\d{8,15}").all()
    assert set(table["source"]) == {"GAGES-II"}
    counts = table["site_class"].value_counts()
    assert counts["reference"] == 2057  # Falcone (2011): 2,057 Ref gages
    assert set(counts.index) == {c.value for c in SiteClass}


def test_packaged_table_follows_its_rules(table):
    ref = table["gagesii_class"] == "Ref"
    assert (table.loc[ref, "site_class"] == "reference").all()
    nonref = table[~ref]
    storage = nonref["norm_storage_days"]
    imperv = nonref["imperv_pct_2006"]
    expect = np.where(
        storage > REGULATED_STORAGE_DAYS,
        "regulated",
        np.where(imperv > URBAN_IMPERVIOUS_PCT, "urban", "unknown"),
    )
    assert (nonref["site_class"].to_numpy() == expect).all()
    # norm_storage_days is STOR_NID_2009 / RUNAVE7100 * 365 (stored rounded).
    ok = table["runave7100_mm"] > 0
    recomputed = table.loc[ok, "stor_nid_2009_ml_km2"] / table.loc[ok, "runave7100_mm"] * 365
    assert np.allclose(recomputed, table.loc[ok, "norm_storage_days"], rtol=1e-3, atol=0.2)


@pytest.mark.parametrize(
    "site, cls",
    [
        (REFERENCE_SITE, SiteClass.REFERENCE),
        (REGULATED_SITE, SiteClass.REGULATED),
        (URBAN_SITE, SiteClass.URBAN),
        (UNKNOWN_SITE, SiteClass.UNKNOWN),
    ],
)
def test_classify_from_table(site, cls):
    c = classify_site(site)
    assert c.site_class is cls and c.in_table
    assert "GAGES-II 2011-09-30" in c.basis[0]


def test_peak_code_6_makes_regulated_and_code_5_does_not():
    c = classify_site(REFERENCE_SITE, ["", "6", "2,6", "5", None])
    assert c.site_class is SiteClass.REGULATED
    assert (c.n_peaks_code6, c.n_peaks_code5, c.n_peaks) == (2, 1, 5)
    only5 = classify_site(REFERENCE_SITE, ["5", "5", ""])
    assert only5.site_class is SiteClass.REFERENCE
    assert any("code 5" in b for b in only5.basis)


def test_absent_site_is_unknown_and_proceeds_with_warning(caplog):
    c = classify_site(ABSENT_SITE)
    assert c.site_class is SiteClass.UNKNOWN and not c.in_table and c.attributes == {}
    with caplog.at_level(logging.WARNING, logger="flowfreq.regulation"):
        prov = require_unregulated(c)
    assert prov["regulation_override"] is False
    assert "no regulation classification" in caplog.text


def test_require_unregulated_refuses_then_records_override(caplog):
    c = classify_site(REGULATED_SITE)
    with pytest.raises(RegulatedRecordError, match="allow_regulated=True"):
        require_unregulated(c)
    with caplog.at_level(logging.WARNING, logger="flowfreq.regulation"):
        prov = require_unregulated(c, allow_regulated=True)
    assert prov["regulation_override"] is True and prov["site_class"] == "regulated"
    assert "allow_regulated=True" in caplog.text
    # A non-regulated class never records an override, even when allowed.
    assert (
        require_unregulated(classify_site(URBAN_SITE), allow_regulated=True)["regulation_override"]
        is False
    )


def test_load_screen_rejects_bad_table(tmp_path):
    bad = tmp_path / "s.csv"
    pd.DataFrame([{c: "x" for c in SCREEN_COLUMNS} | {"site_no": "01234567"}]).to_csv(
        bad, index=False
    )
    with pytest.raises(ValueError, match="site_class"):
        load_screen(bad)


# ---------------------------------------------------------------- run_ffa ----

PEAKS = np.array([1200.0, 3400, 2100, 5600, 1800, 2900, 4100, 2500, 3300, 1500, 2700, 3900])
YEARS = np.arange(2000, 2000 + len(PEAKS))


def test_run_ffa_without_site_no_is_unchanged():
    r = run_ffa(PEAKS, YEARS, station_skew_only=True)
    assert r["error"] is None and r["parameters"]["site_classification"] is None


def test_run_ffa_refuses_regulated_site():
    with pytest.raises(RegulatedRecordError):
        run_ffa(PEAKS, YEARS, station_skew_only=True, site_no=REGULATED_SITE)
    with pytest.raises(RegulatedRecordError):
        run_ffa(
            PEAKS,
            YEARS,
            station_skew_only=True,
            site_no=REFERENCE_SITE,
            peak_codes=["6"] + [""] * (len(PEAKS) - 1),
        )


def test_run_ffa_override_is_recorded():
    r = run_ffa(PEAKS, YEARS, station_skew_only=True, site_no=REGULATED_SITE, allow_regulated=True)
    assert r["error"] is None
    sc = r["parameters"]["site_classification"]
    assert sc["site_class"] == "regulated" and sc["regulation_override"] is True


def test_run_ffa_reference_site_records_classification():
    r = run_ffa(PEAKS, YEARS, station_skew_only=True, site_no=REFERENCE_SITE)
    assert r["parameters"]["site_classification"]["site_class"] == "reference"


# ------------------------------------------------------------ analyze_gage ---


class _FakeGage:
    """Stands in for USGSgage: peaks only, no network."""

    def __init__(self, site_no):
        self.site_no = site_no
        self.site_name = "fake"
        self.peak_data = pd.DataFrame()

    def download_daily_flow(self):
        raise RuntimeError("offline")

    def download_peak_flow(self):
        self.peak_data = pd.DataFrame(
            {"water_year": YEARS, "peak_flow_cfs": PEAKS, "qualification_code": [""] * len(PEAKS)}
        )


def test_analyze_gage_refuses_regulated_before_fitting(monkeypatch, tmp_path):
    import flowfreq

    monkeypatch.setattr(flowfreq, "USGSgage", _FakeGage)

    def boom(*a, **k):  # the fit must never be reached
        raise AssertionError("fitted a refused record")

    monkeypatch.setattr(flowfreq, "_gage_analysis", boom)
    with pytest.raises(RegulatedRecordError):
        flowfreq.analyze_gage(REGULATED_SITE, output_dir=str(tmp_path))
    with pytest.raises(AssertionError, match="fitted a refused record"):
        flowfreq.analyze_gage(REGULATED_SITE, output_dir=str(tmp_path), allow_regulated=True)


# ---------------------------------------------------------------- live -------


@pytest.mark.requires_network
def test_gagesii_release_still_published():
    """The ScienceBase release the table is built from still lists its zip."""
    item = requests.get(
        "https://www.sciencebase.gov/catalog/item/631405bbd34e36012efa304a",
        params={"format": "json"},
        timeout=60,
    ).json()
    names = {f["name"]: f.get("size") for f in item["files"]}
    assert names.get("basinchar_and_report_sept_2011.zip") == 55147470
    assert "10.5066/P96CPHOT" in item.get("citation", "")


@pytest.mark.requires_network
def test_nid_national_csv_is_live():
    """NID's national download answers (used by load_nid)."""
    try:
        r = requests.get("https://nid.sec.usace.army.mil/api/nation/csv", stream=True, timeout=60)
    except (requests.ConnectionError, requests.Timeout) as exc:
        # A local resolver/network failure is not a NID change (seen 2026-10-04:
        # the host resolves via public DNS but not from one developer machine).
        pytest.skip(f"NID unreachable from here: {exc}")
    try:
        assert r.status_code == 200
        head = next(r.iter_content(200)).decode("utf-8", "replace")
    finally:
        r.close()
    assert head.startswith("Data Last Updated")


# ------------------------------------------- current NID / NLCD (opt-in) -----
# Trimmed live captures (2026-10-03/04): StreamStats watersheds with their
# polygon coordinates rounded to 1e-5 degrees (about 1 m), and the NID national
# CSV cut to the dams near those basins plus two associated structures and one
# row without coordinates.

FIX = REPO_ROOT / "tests" / "fixtures" / "regulation"
PARLEYS = "10171500"  # Parleys Creek nr SLC, UT: not in GAGES-II; Mountain Dell Reservoir
EAST_CANYON = "10134500"  # East Canyon Creek nr Morgan, UT: GAGES-II regulated
SWEETWATER = "02337000"  # Sweetwater Creek nr Austell, GA: GAGES-II urban


def _ws(site: str) -> WatershedCharacteristics:
    with gzip.open(FIX / f"ws_{site}.json.gz", "rt", encoding="utf-8") as fh:
        return WatershedCharacteristics.from_dict(json.load(fh))


def _nid_text() -> str:
    with gzip.open(FIX / "nid_capture.csv.gz", "rt", encoding="utf-8") as fh:
        return fh.read()


@pytest.fixture(scope="module")
def nid() -> reg.NIDInventory:
    dams, updated = reg.parse_nid_csv(_nid_text())
    return reg.NIDInventory(dams, updated, "2026-10-04")


def test_parse_nid_csv_drops_associated_and_unlocated():
    dams, updated = reg.parse_nid_csv(_nid_text())
    assert updated == "2026-09-23"
    assert tuple(dams.columns) == reg.NID_COLUMNS
    assert len(dams) == 31  # 34 rows less 2 associated structures and 1 without coordinates
    assert dams["nid_id"].is_unique and dams["latitude"].notna().all()


def test_parse_nid_csv_rejects_other_text():
    with pytest.raises(ValueError, match="last-updated"):
        reg.parse_nid_csv("Dam Name,NID ID\nx,y\n")


def test_load_nid_caches_and_falls_back(tmp_path, monkeypatch, caplog):
    calls = []

    def fake_download(timeout, retries):
        calls.append(1)
        return _nid_text()

    monkeypatch.setattr(reg, "_download_nid", fake_download)
    inv = reg.load_nid(tmp_path)
    assert len(inv.dams) == 31 and inv.data_last_updated == "2026-09-23" and calls == [1]
    assert (tmp_path / reg.NID_CACHE_FILE).exists()
    assert len(reg.load_nid(tmp_path).dams) == 31 and calls == [1]  # served from cache

    def down(timeout, retries):
        raise reg.NIDUnavailableError("host unreachable")

    monkeypatch.setattr(reg, "_download_nid", down)
    with caplog.at_level(logging.WARNING, logger="flowfreq.regulation"):
        stale = reg.load_nid(tmp_path, refresh=True)
    assert len(stale.dams) == 31 and "stale NID cache" in caplog.text
    with pytest.raises(reg.NIDUnavailableError):
        reg.load_nid(tmp_path / "empty")


def test_nid_cache_dir_honours_env(monkeypatch, tmp_path):
    monkeypatch.setenv("FLOWFREQ_CACHE", str(tmp_path))
    assert reg.nid_cache_dir() == tmp_path / "nid"


def test_points_in_geometry_with_hole_and_multipolygon():
    square = [[0, 0], [4, 0], [4, 4], [0, 4], [0, 0]]
    hole = [[1, 1], [2, 1], [2, 2], [1, 2], [1, 1]]
    far = [[10, 10], [11, 10], [11, 11], [10, 11], [10, 10]]
    geom = {"type": "MultiPolygon", "coordinates": [[square, hole], [far]]}
    lon = np.array([0.5, 1.5, 3.0, 10.5, 5.0, -1.0])
    lat = np.array([0.5, 1.5, 3.0, 10.5, 5.0, 2.0])
    got = reg.points_in_geometry(lon, lat, geom, chunk=2)
    assert got.tolist() == [True, False, True, True, False, False]
    with pytest.raises(ValueError):
        reg.points_in_geometry(lon, lat, {"type": "Point", "coordinates": [0, 0]})


def test_nid_storage_in_basin_matches_gagesii(nid, table):
    st = reg.nid_storage_in_basin(_ws(EAST_CANYON).polygon_geojson["geometry"], nid)
    assert st.n_dams == 8 and st.data_last_updated == "2026-09-23"
    assert st.storage_acft == pytest.approx(58691.7)
    assert st.area_km2 == pytest.approx(152.28 * reg.SQ_MI_TO_KM2, rel=0.01)
    # GAGES-II (NID 2009): 9 dams, 159.67 ML/km2. Today's NID: 8 dams, 15% more storage.
    assert st.storage_ml_km2 == pytest.approx(183.56, abs=0.01)
    assert st.storage_ml_km2 == pytest.approx(
        table.loc[EAST_CANYON, "stor_nid_2009_ml_km2"], rel=0.2
    )


def test_runoff_normalization_is_storage_over_a_days_flow():
    area_km2, q_cfs, storage_ml = 100.0, 50.0, 12_345.0
    runoff = reg.runoff_mm_from_flow(q_cfs, area_km2)
    days = storage_ml / area_km2 / runoff * 365
    assert days == pytest.approx(storage_ml * 1000 / (q_cfs * reg.CFS_TO_M3S * 86400))


def test_current_nid_on_gagesii_gage_uses_runave(nid):
    c = classify_site(EAST_CANYON, use_current_nid=True, watershed=_ws(EAST_CANYON), nid=nid)
    assert c.site_class is SiteClass.REGULATED
    cur = c.current
    assert cur["nid_runoff_source"] == "GAGES-II RUNAVE7100" and cur["nid_runoff_mm"] == 212.9
    assert cur["nid_norm_storage_days"] == pytest.approx(314.7)
    assert "imperv_pct" not in cur  # not requested
    assert any("current NID (2026-09-23): 8 dams" in b for b in c.basis)


def test_current_nid_classifies_gage_outside_gagesii(nid):
    plain = classify_site(PARLEYS)
    assert plain.site_class is SiteClass.UNKNOWN and not plain.in_table
    c = classify_site(
        PARLEYS,
        use_current_nid=True,
        use_current_impervious=True,
        watershed=_ws(PARLEYS),
        nid=nid,
        mean_flow_cfs=20.0,
    )
    assert c.site_class is SiteClass.REGULATED and c.has_current_evidence
    assert c.current["nid_n_dams"] == 2 and c.current["nid_runoff_source"].startswith("mean flow")
    assert c.current["imperv_year"] == 2011 and c.current["imperv_pct"] == 0.48
    prov = require_unregulated(c, allow_regulated=True)
    assert prov["current"]["nid_ids"] == ["UT00755", "UT00221"]


def test_current_nid_mean_flow_from_daily(nid, monkeypatch):
    import flowfreq.waterdata as wd

    idx = pd.date_range("2000-10-01", "2002-09-30", freq="D")
    full = pd.DataFrame({"flow_cfs": 20.0, "qualification_code": "A"}, index=idx)
    partial = pd.DataFrame(
        {"flow_cfs": [999.0], "qualification_code": ["A"]},
        index=pd.DatetimeIndex(["2003-01-01"]),
    )  # a lone day of water year 2003, which is ignored
    daily = pd.concat([full, partial])
    monkeypatch.setattr(wd, "download_daily", lambda site_no: daily)
    assert reg.mean_flow_from_daily(PARLEYS) == (20.0, 2)
    c = classify_site(PARLEYS, use_current_nid=True, watershed=_ws(PARLEYS), nid=nid)
    assert c.current["nid_runoff_source"] == "mean daily flow, 2 complete water years"

    monkeypatch.setattr(wd, "download_daily", lambda site_no: daily.iloc[:0])
    with pytest.raises(ValueError, match="no daily discharge"):
        reg.mean_flow_from_daily(PARLEYS)
    c = classify_site(PARLEYS, use_current_nid=True, watershed=_ws(PARLEYS), nid=nid)
    assert c.site_class is SiteClass.UNKNOWN
    assert any("no mean runoff" in b for b in c.basis)


def test_current_impervious_time_series_and_year():
    c = classify_site(SWEETWATER, use_current_impervious=True, watershed=_ws(SWEETWATER))
    assert c.current["imperv_series"] == {"2006": 9.33, "2011": 10.559, "2019": 12.94}
    assert (c.current["imperv_year"], c.current["imperv_pct"]) == (2019, 12.94)
    assert c.site_class is SiteClass.URBAN
    # GAGES-II's NLCD 2006 value and StreamStats' LC06IMP agree for this basin.
    assert c.attributes["imperv_pct_2006"] == pytest.approx(9.33, abs=0.1)


def test_current_impervious_from_unverified_region_is_not_used():
    ws = _ws(SWEETWATER)
    ws.region = "CO"
    c = classify_site(SWEETWATER, use_current_impervious=True, watershed=ws)
    assert "imperv_pct" not in c.current and "imperv_series_unverified" in c.current
    assert c.site_class is SiteClass.URBAN  # falls back to GAGES-II NLCD 2006
    assert any("not used" in b for b in c.basis)


def test_region_without_impervious_is_recorded():
    ws = _ws(SWEETWATER)
    ws.characteristics = {}
    c = classify_site(PARLEYS, use_current_impervious=True, watershed=ws)
    assert c.site_class is SiteClass.UNKNOWN and not c.has_current_evidence
    assert any("computes no NLCD impervious" in b for b in c.basis)


def test_reference_class_is_retained_with_note():
    c = classify_site(REFERENCE_SITE, use_current_impervious=True, watershed=_ws(SWEETWATER))
    assert c.site_class is SiteClass.REFERENCE
    assert any("GAGES-II Ref class retained" in b for b in c.basis)


def test_streamstats_regulated_warning_is_decisive(monkeypatch):
    import flowfreq.streamstats as ss

    def refuse(*a, **k):
        raise ss.DegenerateDelineationError(
            "StreamStats delineation warning at (32.8, -83.7) in region 'GA': River is "
            "regulated, streamflow characteristics  not valid."
        )

    monkeypatch.setattr(ss, "delineate_and_get_characteristics", refuse)
    c = classify_site(ABSENT_SITE, use_current_impervious=True, region="GA", lat=32.8, lon=-83.7)
    assert c.site_class is SiteClass.REGULATED and "streamstats_regulated" in c.current
    with pytest.raises(RegulatedRecordError, match="River is"):
        require_unregulated(c)


def test_delineation_failure_leaves_class_and_warns(monkeypatch, caplog):
    import flowfreq.streamstats as ss

    def fail(*a, **k):
        raise ss.UnsnappablePointError("could not snap")

    monkeypatch.setattr(ss, "delineate_and_get_characteristics", fail)
    c = classify_site(ABSENT_SITE, use_current_nid=True, region="GA", lat=32.8, lon=-83.7)
    assert c.site_class is SiteClass.UNKNOWN and not c.has_current_evidence
    assert any("delineation failed" in b for b in c.basis)
    with caplog.at_level(logging.WARNING, logger="flowfreq.regulation"):
        require_unregulated(c)
    assert "no regulation classification" in caplog.text


def test_peak_code_6_still_overrides_current_data():
    c = classify_site(SWEETWATER, ["6", ""], use_current_impervious=True, watershed=_ws(SWEETWATER))
    assert c.site_class is SiteClass.REGULATED


@pytest.mark.requires_network
def test_streamstats_still_lists_the_impervious_codes():
    """The region lists NLCD_IMPERVIOUS_CODES was verified against (2026-10-03)."""
    expect = {"GA": {"LC06IMP", "LC11IMP", "LC19IMP"}, "WY": {"LC16IMP"}, "NV": {"LC21IMP"}}
    for region, codes in expect.items():
        listed = requests.get(
            f"https://streamstats.usgs.gov/ss-hydro/v1/basin-characteristics/{region}",
            timeout=60,
        ).json()
        assert codes <= {x["code"] for x in listed}
        assert codes <= set(reg.NLCD_IMPERVIOUS_CODES)


@pytest.mark.requires_network
def test_current_classification_live(tmp_path):
    """End to end: NID download (skipped if the host is unreachable) and StreamStats."""
    try:
        inv = reg.load_nid(tmp_path, retries=2)
    except reg.NIDUnavailableError as exc:
        pytest.skip(f"NID unreachable: {exc}")
    assert len(inv.dams) > 80_000
    c = classify_site(
        PARLEYS, use_current_nid=True, use_current_impervious=True, nid=inv, mean_flow_cfs=20.0
    )
    assert c.current["nid_n_dams"] >= 2 and c.current["imperv_year"] == 2011
    assert c.site_class is SiteClass.REGULATED
