"""Tests for flowfreq.regulation (issue #32): the GAGES-II screen and the B17C refusal."""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import pytest
import requests

from flowfreq.regulation import (
    SCREEN_COLUMNS,
    RegulatedRecordError,
    SiteClass,
    classify_site,
    load_screen,
    require_unregulated,
)
from flowfreq.workflow import run_ffa

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
    """NID's national download answers (not yet used; see the build tool)."""
    r = requests.get("https://nid.sec.usace.army.mil/api/nation/csv", stream=True, timeout=60)
    try:
        assert r.status_code == 200
        head = next(r.iter_content(200)).decode("utf-8", "replace")
    finally:
        r.close()
    assert head.startswith("Data Last Updated")
