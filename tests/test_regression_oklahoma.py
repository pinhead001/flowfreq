"""Oklahoma's verified equations (SIR 2019-5143 ver. 1.1 and SIR 2015-5134)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "flowfreq" / "data" / "regression"
AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]
NSS = {
    "Region 1": "GC1487",
    "Region 1 NRCS": "GC1488",
    "Region 2": "GC1791",
    "Region 2 NRCS": "GC1792",
    "Panhandle": "GC1629",
    "Panhandle NRCS": "GC1749",
}
REMOVED = {"07325860", "07332600", "07311200"}  # SIR 2019-5143, p. 15


def _sig(x: float, figures: int = 3) -> float:
    return float(f"{x:.{figures - 1}e}")


def _pct_to_log(p: float) -> float:
    return math.sqrt(math.log(1 + (p / 100) ** 2)) / math.log(10)


@pytest.fixture(scope="module")
def lib():
    return load_state("OK")


@pytest.fixture(scope="module")
def gages():
    path = ROOT / "tests" / "fixtures" / "regression" / "oklahoma_gages.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _q(lib, region, aep, **kw):
    return evaluate(lib.equation(region, aep), kw, allow_extrapolation=True).flow_cfs


def test_oklahoma_is_verified(lib):
    assert lib.status == "verified"
    assert sorted(lib.regions) == sorted(NSS)
    assert len(lib.equations) == 42
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
        assert [e.aep for e in lib.equations if e.region_code == region] == AEPS


def test_equations_as_printed(lib):
    """A second entry of the constants, independent of the parser that built the file."""
    printed = {
        "Region 1": [
            (1.2, 0.54, 1.44),
            (1.82, 0.56, 1.55),
            (1.82, 0.55, 1.68),
            (1.90, 0.54, 1.81),
            (1.92, 0.54, 1.91),
            (1.95, 0.53, 1.98),
            (2.09, 0.51, 2.13),
        ],
        "Region 2": [
            (61.6, 0.40, 0.75),
            (97.7, 0.44, 0.77),
            (126, 0.46, 0.78),
            (174, 0.47, 0.78),
            (204, 0.50, 0.79),
            (240, 0.50, 0.79),
            (363, 0.51, 0.80),
        ],
        "Panhandle": [
            (77.6, 0.34),
            (209, 0.35),
            (339, 0.36),
            (562, 0.36),
            (759, 0.37),
            (1000, 0.38),
            (1620, 0.39),
        ],
    }
    for region, rows in printed.items():
        for suffix in ("", " NRCS"):
            for aep, row in zip(AEPS, rows):
                eq = lib.equation(region + suffix, aep)
                assert 10**eq.intercept == pytest.approx(row[0], rel=1e-12)
                assert tuple(eq.coefficients) == row[1:]


def test_regulated_sets_use_the_unregulated_area(lib):
    assert [v.code for v in lib.equation("Region 1", 0.01).variables] == ["CONTDA", "PRECIP"]
    assert [v.code for v in lib.equation("Region 1 NRCS", 0.01).variables] == [
        "DAUNREG",
        "PRECIP",
    ]
    assert [v.code for v in lib.equation("Region 2", 0.01).variables] == ["CSL10_85fm", "CONTDA"]
    assert [v.code for v in lib.equation("Region 2 NRCS", 0.01).variables] == [
        "CSL10_85fm",
        "DAUNREG",
    ]
    assert [v.code for v in lib.equation("Panhandle", 0.01).variables] == ["CONTDA"]
    assert [v.code for v in lib.equation("Panhandle NRCS", 0.01).variables] == ["DAUNREG"]


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "OK.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-12)
        assert list(e.coefficients) == s["coefficients"]
        # Pinned difference: NSS rounds the Panhandle's Sp to one decimal (147.8 for 147.82).
        tol = 1.2e-4 if e.region_code.startswith("Panhandle") else 1e-6
        assert e.sep_log == pytest.approx(s["sep_log"], abs=tol)
        for v, n in zip(e.variables, s["variables"]):
            assert v.code == n["code"]
            assert v.maximum == n["maximum"]
            if e.region_code.startswith("Panhandle"):
                # Pinned difference: the report prints no minimum. 1.21 is Table 1's
                # lowest gage; NSS has 1.0.
                assert (v.minimum, n["minimum"]) == (1.21, 1.0)
            else:
                assert v.minimum == n["minimum"]


def test_table_5_kiamichi_river(lib):
    """SIR 2019-5143 Table 5: 07335700, 54.87 ft/mi, 39.63 mi2, equations 11-17."""
    printed = [4830, 9670, 14000, 20200, 27600, 32500, 53100]
    flows = [_q(lib, "Region 2", a, CSL10_85fm=54.87, CONTDA=39.63) for a in AEPS]
    exact = [_sig(f) == p for f, p in zip(flows, printed)]
    # Two cells miss by one unit of the third figure: 9,676 and 27,653. The report
    # evaluated unrounded coefficients; the rest round to the printed value.
    assert exact == [True, False, True, True, False, True, True]
    for f, p in zip(flows, printed):
        assert abs(f - p) <= 100 if p >= 10000 else abs(f - p) <= 10


def test_table_5_weighting_row_for_10_percent_is_pinned():
    """Equation 18 from Table 5's own row gives 19,900 where it prints 20,000."""
    q_s, n, q_r, e = 20800, 52, 14000, 8
    assert _sig((q_s * n + q_r * e) / (n + e)) == 19900
    q_s, q_r, e, printed = 41400, 32500, 12, 39700
    assert _sig((q_s * n + q_r * e) / (n + e)) == printed


def test_ungaged_example_uses_the_gage_slope(lib):
    """Pinned: the example states 42.0 ft/mi, but its 19,300 comes from the gage's 54.87."""
    assert _sig(_q(lib, "Region 2", 0.01, CSL10_85fm=54.87, CONTDA=20.5)) == 19300
    assert _sig(_q(lib, "Region 2", 0.01, CSL10_85fm=42.0, CONTDA=20.5)) == 16900


def test_urban_example(lib):
    q1 = _q(lib, "Region 2", 0.01, CSL10_85fm=11.5, CONTDA=25.0)
    q50 = _q(lib, "Region 2", 0.5, CSL10_85fm=11.5, CONTDA=25.0)
    assert q1 == pytest.approx(10300, abs=100)
    assert _sig(q50) == 1830
    # Equation 26 with RL = 3, on the example's own printed flows.
    assert _sig(2.72 * 2 * 1830 + 0.167 * 4 * 10300) == 16800


def test_regulated_example_uncle_johns_creek(lib):
    """Region 1, DAUNREG 65.1 mi2, PRECIP 31.0 in., equation 33."""
    assert _sig(_q(lib, "Region 1 NRCS", 0.01, DAUNREG=65.1, PRECIP=31.0)) == 16000


def test_panhandle_gage_estimates(lib, gages):
    """SIR 2015-5134 Table 1: all 224 printed regression estimates.

    211 round to the printed three figures. The other 13 sit on a rounding boundary,
    under 0.6 of a unit of the third figure away, as rounded constants and a rounded
    CONTDA would leave them; the report does not say what it evaluated.
    """
    assert len(gages["panhandle"]) == 32
    missed = 0
    for g in gages["panhandle"]:
        for aep, printed in zip(AEPS, g["regression"]):
            flow = _q(lib, "Panhandle", aep, CONTDA=g["CONTDA"])
            unit = 10 ** (math.floor(math.log10(printed)) - 2)
            assert abs(flow - printed) < 0.6 * unit, (g["site"], aep)
            missed += _sig(flow) != printed
    assert missed == 13
    table = np.array([g["regression"] for g in gages["panhandle"]])
    assert table.mean(axis=0) == pytest.approx(gages["panhandle_average_regression"], rel=0.005)


def test_panhandle_rmse_row_is_not_an_rmse(gages):
    """Pinned: the printed RMSE exceeds any RMS difference of the table's own columns."""
    reg = np.array([g["regression"] for g in gages["panhandle"]])
    comp = np.array([g["computed"] for g in gages["panhandle"]])
    printed = np.array(gages["panhandle_rmse"])
    rms = np.sqrt(((reg - comp) ** 2).mean(axis=0))
    assert np.all(printed > 2.5 * rms)
    # It sits within 1 percent of the root of the summed squared estimates.
    assert np.sqrt((reg**2).sum(axis=0)) == pytest.approx(printed, rel=0.01)


@pytest.mark.parametrize("region", ["1", "2"])
def test_statewide_equations_are_centered_on_table_1(lib, gages, region):
    """SIR 2019-5143 prints no per-gage estimates; check against its at-site flows."""
    rows = [g for g in gages["statewide"] if g["region"] == region and g["station"] not in REMOVED]
    assert len(rows) == lib.equation(f"Region {region}", 0.01).n_sites
    limit = {"1": (0.003, 0.033), "2": (-0.012, 0.012)}[region]
    for i, aep in enumerate(AEPS):
        res = []
        for g in rows:
            kw = {"CONTDA": g["CONTDA"]}
            kw["PRECIP" if region == "1" else "CSL10_85fm"] = g[
                "PRECIP" if region == "1" else "CSL10_85fm"
            ]
            res.append(
                math.log10(g["at_site"][i]) - math.log10(_q(lib, f"Region {region}", aep, **kw))
            )
        assert limit[0] <= float(np.mean(res)) <= limit[1], aep


def test_table_1_counts_and_removed_gages(gages):
    rows = gages["statewide"]
    assert len(rows) == 212
    assert sum(g["region"] == "1" for g in rows) == 57
    assert sum(g["region"] == "2" for g in rows) == 155
    assert {g["station"]: g["region"] for g in rows if g["station"] in REMOVED} == {
        "07325860": "1",
        "07332600": "2",
        "07311200": "2",
    }


def test_printed_ranges_and_the_gages_that_exceed_them(lib, gages):
    """Report wins: the printed limits are stored, though Table 1's gages do not span them."""
    eq1, eq2 = lib.equation("Region 1", 0.01), lib.equation("Region 2", 0.01)
    assert {v.code: (v.minimum, v.maximum) for v in (*eq1.variables, *eq2.variables)} == {
        "CONTDA": (0.1, 2510.0),
        "PRECIP": (16.6, 62.1),
        "CSL10_85fm": (1.98, 342.0),
    }
    rows = gages["statewide"]
    assert min(g["CONTDA"] for g in rows) == 0.08
    assert max(g["CONTDA"] for g in rows) == 2515.63
    assert min(g["PRECIP"] for g in rows) == 22.0
    assert max(g["PRECIP"] for g in rows if g["region"] == "1") == 36.0
    assert max(g["CSL10_85fm"] for g in rows) == 182.11
    pan = lib.equation("Panhandle", 0.01).variables[0]
    assert (pan.minimum, pan.maximum) == (1.21, 2060.0)
    assert min(g["CONTDA"] for g in gages["panhandle"]) == 1.21
    assert max(g["CONTDA"] for g in gages["panhandle"]) == 2057.39


def test_table_4_error_statistics(lib):
    table = {
        "Region 1": [
            (42.99, 41.90, 2),
            (34.09, 32.97, 5),
            (32.52, 31.28, 8),
            (41.34, 39.56, 9),
            (37.48, 35.93, 11),
            (40.34, 38.65, 12),
            (51.20, 49.32, 12),
        ],
        "Region 2": [
            (46.88, 46.02, 2),
            (36.23, 35.38, 5),
            (35.03, 33.79, 8),
            (39.88, 38.82, 9),
            (37.07, 35.93, 11),
            (39.90, 38.65, 12),
            (50.70, 49.32, 12),
        ],
    }
    for region, rows in table.items():
        for suffix in ("", " NRCS"):
            for aep, (sp, sme, ey) in zip(AEPS, rows):
                eq = lib.equation(region + suffix, aep)
                assert eq.sep_log == pytest.approx(_pct_to_log(sp), abs=1e-6)
                assert eq.model_error_variance == pytest.approx(_pct_to_log(sme) ** 2, abs=1e-6)
                assert eq.equivalent_years == ey
                assert eq.covariance is None
    panhandle = [
        (147.82, 138.48),
        (103.99, 97.52),
        (92.64, 86.40),
        (85.76, 79.24),
        (84.01, 77.09),
        (83.72, 76.33),
        (91.87, 83.62),
    ]
    for suffix in ("", " NRCS"):
        for aep, (sp, sme) in zip(AEPS, panhandle):
            eq = lib.equation("Panhandle" + suffix, aep)
            assert eq.sep_log == pytest.approx(_pct_to_log(sp), abs=1e-6)
            assert eq.model_error_variance == pytest.approx(_pct_to_log(sme) ** 2, abs=1e-6)
            assert eq.n_sites == 32
            assert eq.equivalent_years is None


def test_figure_4_disagrees_with_table_4_and_table_4_is_stored(lib):
    """Pinned: Figure 4 (Region 2, 10 percent) shows Sp 33.89, SME 32.97; Table 4 35.03, 33.79."""
    eq = lib.equation("Region 2", 0.1)
    assert eq.sep_log == pytest.approx(_pct_to_log(35.03), abs=1e-6)
    assert eq.sep_log != pytest.approx(_pct_to_log(33.89), abs=1e-3)
    assert eq.coefficients == (0.46, 0.78)  # Figure 4: 0.465, 0.780, intercept 2.10
    assert eq.intercept == pytest.approx(2.10, abs=0.005)


def test_table_4_model_error_repeats_across_regions_is_pinned(lib):
    """Pinned: both regions print the same SME at 2, 1 and 0.2 percent, with different Sp."""
    for aep in (0.02, 0.01, 0.002):
        one, two = lib.equation("Region 1", aep), lib.equation("Region 2", aep)
        assert one.model_error_variance == two.model_error_variance
        assert one.sep_log != two.sep_log


def test_intervals_use_the_sep(lib):
    est = evaluate(lib.equation("Region 2", 0.01), {"CSL10_85fm": 54.87, "CONTDA": 39.63})
    assert est.interval_method == "average standard error of prediction"
    assert est.interval[0] < est.flow_cfs < est.interval[1]
