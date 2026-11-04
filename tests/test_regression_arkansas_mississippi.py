"""Arkansas's and Mississippi's verified equations.

Arkansas: SIR 2016-5081 Tables 7-11 and SIR 2021-5046. Mississippi: SIR 2018-5148 Tables 4-5 with
its data release, and SIR 2021-5046. None of the three reports prints a worked example, so the
per-gage regression flows are the check. They are compared with an interval, not a percentage: the
flow is printed to a few figures and so are the characteristics (Arkansas prints SOILINDEX to one
decimal, which alone moves a region B2 Q50% by up to 5.7 percent), so each published flow must lie
inside the range the printed precision of its inputs and of itself allows.
"""

from __future__ import annotations

import json
import math
from dataclasses import replace
from decimal import Decimal
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from flowfreq.regression import evaluate
from flowfreq.regression.equations import Variable
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
FIXTURES = ROOT / "fixtures" / "regression"
NSS_AR = {
    "A": "GC1715",
    "B1": "GC1716",
    "B2": "GC1717",
    "C": "GC1718",
    "D": "GC1719",
    "MAP": "GC1917",
}
NSS_MS = {"A": "GC1772", "B": "GC1773", "C": "GC1774", "MAP": "GC1912"}
AEP7 = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]
AEP8 = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
AR_GAGES = json.loads((FIXTURES / "arkansas_gages.json").read_text("utf-8"))
MS_GAGES = json.loads((FIXTURES / "mississippi_gages.json").read_text("utf-8"))
MAP_GAGES = json.loads((FIXTURES / "alluvial_plain_gages.json").read_text("utf-8"))


def _snapshot(state):
    snap = json.loads((DATA / "nss_snapshot" / f"{state}.json").read_text(encoding="utf-8"))
    return {(s["region_code"], s["aep"]): s for s in snap["equations"]}


@pytest.fixture(scope="module")
def ar():
    return load_state("AR")


@pytest.fixture(scope="module")
def ms():
    return load_state("MS")


def _half_unit(printed: str) -> float:
    """Half a unit in the last place of a value as printed."""
    return 0.5 * 10.0 ** Decimal(printed).as_tuple().exponent


def _term_bounds(v: Variable, coefficient: float, printed: str):
    """The range of one regression term over the printed value's rounding interval."""
    x, h = float(printed), _half_unit(printed)
    ends = [
        coefficient
        * (math.log10(v.scale * a + v.offset) if v.transform == "log10" else v.scale * a + v.offset)
        for a in (x - h, x + h)
        if v.transform != "log10" or v.scale * a + v.offset > 0
    ]
    return min(ends), max(ends)


def _flow_half_unit(flow: float, figures: int) -> float:
    """Half a unit at the printed figures, and never less than half a cubic foot per second."""
    magnitude = math.floor(math.log10(flow))
    return max(0.5 * 10.0 ** (magnitude - figures + 1), 0.5)


def _inside(eq, printed_by_code, flow, figures):
    """Whether a published flow lies in the range the printed precision allows.

    Inputs are rounded where they are printed, so each term is taken over its rounding interval.
    The equation is rounded too: coefficients to three decimals, and the constant to three decimals
    of its exponent (SIR 2016-5081) or to a whole number (SIR 2018-5148, SIR 2021-5046).
    """
    bounds = [
        _term_bounds(v, c, printed_by_code[v.code]) for v, c in zip(eq.variables, eq.coefficients)
    ]
    constant = eq.citation.year != 2016
    slack = 0.5 / 10**eq.intercept / math.log(10) if constant else 5e-4
    for v in eq.variables:
        term = v.value({v.code: float(printed_by_code[v.code])})
        slack += 5e-4 * abs(float(term))
    lo = 10 ** (eq.intercept + sum(b[0] for b in bounds) - slack)
    hi = 10 ** (eq.intercept + sum(b[1] for b in bounds) + slack)
    half = _flow_half_unit(flow, figures)
    return lo <= flow + half and flow - half <= hi


# ---- the two library files


def test_arkansas_is_verified(ar):
    assert ar.status == "verified"
    assert ar.regions == ["A", "B1", "B2", "C", "D", "MAP"]
    assert len(ar.equations) == 43
    for region, code in NSS_AR.items():
        assert region_for_nss_code(ar, code) == region
    for region in ("A", "B1", "B2", "C", "D"):
        assert ar.aeps(region) == sorted(AEP7, reverse=True) or set(ar.aeps(region)) == set(AEP7)
    assert set(ar.aeps("MAP")) == set(AEP8)


def test_mississippi_is_verified(ms):
    assert ms.status == "verified"
    assert ms.regions == ["A", "B", "C", "MAP"]
    assert len(ms.equations) == 32
    for region, code in NSS_MS.items():
        assert region_for_nss_code(ms, code) == region


def test_region_b_gate_is_the_reports_500_ft(ar):
    assert ar.region_gates["B1"].admits(500.0) and not ar.region_gates["B1"].admits(499.9)
    assert ar.region_gates["B2"].admits(499.9) and not ar.region_gates["B2"].admits(500.0)


def test_site_counts(ar, ms):
    assert {r: ar.equation(r, 0.01).n_sites for r in ar.regions} == {
        "A": 68,
        "B1": 36,
        "B2": 56,
        "C": 83,
        "D": 38,
        "MAP": 58,
    }
    assert sum(ar.equation(r, 0.01).n_sites for r in ("A", "B1", "B2", "C", "D")) == 281
    assert {r: ms.equation(r, 0.01).n_sites for r in ms.regions} == {
        "A": 27,
        "B": 134,
        "C": 120,
        "MAP": 58,
    }


# ---- Arkansas: Tables 8 and 9


def test_arkansas_table_9_reproduces_at_every_gage(ar):
    assert len(AR_GAGES) == 281
    misfits = {}
    for site, gage in AR_GAGES.items():
        for i, aep in enumerate(AEP7):
            if not _inside(ar.equation(gage["region"], aep), gage, gage["flows"][i], figures=4):
                misfits.setdefault(site, []).append(aep)
    assert not misfits


def test_the_interval_check_catches_a_wrong_coefficient(ar):
    """Moving one exponent by 0.005 has to show up, or the check above proves nothing."""
    gages = [(s, g) for s, g in AR_GAGES.items() if g["region"] == "D"]
    eq = ar.equation("D", 0.01)
    wrong = replace(eq, coefficients=(eq.coefficients[0] + 0.005, *eq.coefficients[1:]))
    misfits = [s for s, g in gages if not _inside(wrong, g, g["flows"][5], figures=4)]
    assert len(misfits) > len(gages) / 3


# ---- Arkansas: the covariance


def _mean_prediction_variance(eq, gages):
    u = np.array(eq.covariance)
    total = 0.0
    for gage in gages:
        x = np.array([1.0] + [float(v.value({v.code: float(gage[v.code])})) for v in eq.variables])
        total += float(x @ u @ x)
    return eq.model_error_variance + total / len(gages)


def test_mev_plus_mean_xux_closes_on_the_printed_avp(ar):
    for region in ("A", "B1", "B2", "C", "D"):
        gages = [g for g in AR_GAGES.values() if g["region"] == region]
        assert len(gages) == ar.equation(region, 0.01).n_sites
        for aep in AEP7:
            eq = ar.equation(region, aep)
            assert _mean_prediction_variance(eq, gages) == pytest.approx(eq.avp, abs=1.1e-3)


def test_every_matrix_is_symmetric_and_positive_definite(ar):
    for eq in ar.equations:
        if eq.region_code == "MAP":
            assert eq.covariance is None
            continue
        u = np.array(eq.covariance)
        assert np.allclose(u, u.T, atol=2e-10), (eq.region_code, eq.aep)
        assert np.linalg.eigvalsh(u)[0] > 0, (eq.region_code, eq.aep)


#: Table 11 as printed, where its cells cannot be right. Each row: region, AEP, the cell, what is
#: printed, what is stored, and the mean prediction variance from the printed matrix.
TABLE_11_MISPRINTS = [
    ("B2", 0.5, (3, 0), -0.010697946, -0.10697946, 0.618),
    ("B2", 0.01, (5, 0), -0.01477342, -0.011477342, 0.0213),
    ("B2", 0.002, (3, 0), -0.017874274, -0.1787427, 0.526),
]


@pytest.mark.parametrize("region,aep,cell,printed,stored,bad", TABLE_11_MISPRINTS)
def test_table_11_misprints_are_settled_by_the_printed_avp(
    ar, region, aep, cell, printed, stored, bad
):
    eq = ar.equation(region, aep)
    gages = [g for g in AR_GAGES.values() if g["region"] == region]
    i, j = cell
    assert eq.covariance[i][j] == stored
    assert _mean_prediction_variance(eq, gages) == pytest.approx(eq.avp, abs=1e-3)

    u = [list(r) for r in eq.covariance]
    u[i][j] = printed
    if aep == 0.5:  # this one is printed the same way in both cells
        u[j][i] = printed
    x = [
        np.array([1.0] + [float(v.value({v.code: float(g[v.code])})) for v in eq.variables])
        for g in gages
    ]
    as_printed = eq.model_error_variance + np.mean([float(r @ np.array(u) @ r) for r in x])
    assert as_printed == pytest.approx(bad, rel=0.01)
    assert abs(as_printed - eq.avp) > 3e-3


def test_table_11_t_values_follow_from_the_site_counts(ar):
    """Spot checks of Table 11's printed 90-percent t, at n - p degrees of freedom."""
    for region, aep, printed in [
        ("A", 0.5, 1.668),
        ("A", 0.01, 1.669),
        ("B1", 0.01, 1.692),
        ("B2", 0.01, 1.676),
        ("D", 0.01, 1.691),
    ]:
        eq = ar.equation(region, aep)
        dof = eq.n_sites - len(eq.variables) - 1
        assert stats.t.ppf(0.95, dof) == pytest.approx(printed, abs=1.5e-3)


def test_sp_matches_avp_but_for_one_printed_cell(ar):
    """Table 7's B1 Q0.2% prints Sp 46.86 percent. Its AVP 0.042 is 49.9 percent, and MEV + mean
    xUx there gives 0.0420, so the Sp is the odd one out. It stays as printed."""
    for eq in ar.equations:
        if eq.region_code == "MAP":
            continue
        close = eq.sep_log**2 == pytest.approx(eq.avp, abs=1.1e-3)
        assert close == ((eq.region_code, eq.aep) != ("B1", 0.002)), (eq.region_code, eq.aep)
    assert ar.equation("B1", 0.002).sep_log ** 2 == pytest.approx(0.0374, abs=1e-4)


def test_covariance_sets_the_interval_for_the_flood_regions(ar):
    est = evaluate(ar.equation("D", 0.01), {"DRNAREA": 374.0, "CSLBlue": 0.81, "BSHAPELFP": 13.39})
    assert est.interval_method.startswith("covariance") or "covariance" in est.interval_method


# ---- Arkansas against NSS


def test_arkansas_equals_nss_once_nss_folds_its_percent_terms(ar):
    snap = _snapshot("AR")
    assert {(NSS_AR[e.region_code], e.aep) for e in ar.equations} == set(snap)
    for e in ar.equations:
        s = snap[NSS_AR[e.region_code], e.aep]
        by_code = {n["code"]: c for n, c in zip(s["variables"], s["coefficients"])}
        assert set(by_code) == {v.code for v in e.variables}
        intercept = e.intercept
        for v, c in zip(e.variables, e.coefficients):
            if v.scale != 1.0:  # Table 7's c(0.01*X+1) is NSS's c/100 * X, plus c in the intercept
                intercept += c * v.offset
                assert by_code[v.code] == pytest.approx(c * v.scale, abs=1e-12)
            else:
                assert by_code[v.code] == pytest.approx(c, abs=1e-12)
        if e.region_code == "MAP":
            assert intercept == pytest.approx(s["intercept"], abs=1e-12)
        else:
            assert intercept == pytest.approx(s["intercept"], abs=1e-9)


#: NSS's Sp for regions A-D is Table 7's "Values from Hodge & Tasker (1995)" column, in percent.
NSS_SP_IS_HODGE_AND_TASKER = {
    "A": [42.0, 32.0, 31.0, 32.0, 34.0, 37.0, 45.0],
    "B1": [42.0, 39.0, 38.0, 38.0, 39.0, 40.0, 43.0],
    "B2": [42.0, 39.0, 38.0, 38.0, 39.0, 40.0, 43.0],
    "C": [42.0, 38.0, 37.0, 36.0, 35.0, 35.0, 37.0],
    "D": [33.0, 31.0, 33.0, 37.0, 41.0, 45.0, 56.0],
}
#: Table 7 Sp, as printed.
TABLE_7_SP = {
    "A": [42.82, 33.98, 32.76, 33.18, 34.98, 36.77, 44.21],
    "B1": [37.80, 38.12, 39.10, 41.39, 42.43, 44.67, 46.86],
    "B2": [37.54, 32.56, 32.81, 33.32, 34.95, 37.66, 42.02],
    "C": [52.01, 43.18, 40.16, 39.71, 40.45, 42.66, 45.76],
    "D": [33.06, 33.31, 36.91, 42.45, 46.64, 50.69, 59.53],
}


def _pct(sep_log):
    return 100 * math.sqrt(math.exp((sep_log * math.log(10)) ** 2) - 1)


def test_nss_serves_the_previous_study_standard_errors(ar):
    snap = _snapshot("AR")
    for region, values in NSS_SP_IS_HODGE_AND_TASKER.items():
        for aep, nss, printed in zip(AEP7, values, TABLE_7_SP[region]):
            assert _pct(snap[NSS_AR[region], aep]["sep_log"]) == pytest.approx(nss, abs=0.01)
            assert _pct(ar.equation(region, aep).sep_log) == pytest.approx(printed, abs=0.01)


#: Where NSS's limits round Table 10, (region, code): (NSS, Table 10).
NSS_RANGE_ROUNDING = {
    ("A", "ELEV"): ((380.0, 1840.0), (380.51, 1832.27)),
    ("B1", "BSHAPELFP"): ((2.94, 15.1), (2.94, 15.06)),
    ("B2", "PRECIP"): ((49.3, 56.6), (49.36, 56.55)),
    ("C", "BSHAPELFP"): ((1.97, 15.1), (1.97, 15.09)),
    ("C", "LC11PAST"): ((0.0, 87.6), (0.0, 87.51)),
    ("D", "CSLBlue"): ((0.36, 31.8), (0.36, 31.74)),
    ("D", "BSHAPELFP"): ((3.12, 76.0), (3.12, 75.97)),
}


def test_arkansas_ranges_are_table_10_and_nss_rounds_them(ar):
    snap = _snapshot("AR")
    seen = set()
    for e in ar.equations:
        if e.region_code == "MAP":
            continue
        s = {n["code"]: n for n in snap[NSS_AR[e.region_code], e.aep]["variables"]}
        for v in e.variables:
            nss = (s[v.code]["minimum"], s[v.code]["maximum"])
            key = (e.region_code, v.code)
            if key in NSS_RANGE_ROUNDING:
                seen.add(key)
                assert (nss, (v.minimum, v.maximum)) == NSS_RANGE_ROUNDING[key]
            else:
                assert nss == (v.minimum, v.maximum), key
    assert seen == set(NSS_RANGE_ROUNDING)


def test_table_10_differs_slightly_from_the_table_8_gages():
    """Table 10 is stored; these are the printed gages it leaves outside."""

    def values(region, code):
        return [float(g[code]) for g in AR_GAGES.values() if g["region"] == region]

    assert min(values("B1", "DRNAREA")) == 0.06  # Table 10: 0.10
    assert max(values("B2", "DRNAREA")) == 2092.0  # Table 10: 2,090
    assert max(values("D", "DRNAREA")) == 1622.0  # Table 10: 1,620


# ---- Mississippi


#: The release's RRE and weighted rows at these two are not theirs: 02484500's implies 76 mi2, the
#: area of its neighbour 02484600, whose rows repeat 02484750's (0.11 mi2).
MS_RELEASE_MISALIGNED = {"02484500", "02484600"}


def test_mississippi_release_flows_reproduce_at_every_gage(ms):
    assert len(MS_GAGES) == 251
    misfits = set()
    for site, gage in MS_GAGES.items():
        for i, aep in enumerate(AEP8):
            if not _inside(ms.equation(gage["region"], aep), gage, gage["flows"][i], figures=3):
                misfits.add(site)
    assert misfits == MS_RELEASE_MISALIGNED


def test_the_misaligned_release_rows_belong_to_their_neighbours(ms):
    eq = ms.equation("C", 0.01)
    assert evaluate(eq, {"DRNAREA": 76.4, "LC11DEV": 0.0}).flow_cfs == pytest.approx(
        MS_GAGES["02484500"]["flows"][5], rel=0.006
    )
    assert MS_GAGES["02484600"]["flows"] == MS_GAGES["02484750"]["flows"]


def test_mississippi_area_limits_are_table_5(ms):
    assert {
        r: (ms.equation(r, 0.01).variables[0].minimum, ms.equation(r, 0.01).variables[0].maximum)
        for r in "ABC"
    } == {
        "A": (1.41, 612.0),
        "B": (0.15, 1750.0),
        "C": (0.05, 1010.0),
    }
    for region in "ABC":
        areas = [float(g["DRNAREA"]) for s, g in MS_GAGES.items() if g["region"] == region]
        limits = ms.equation(region, 0.01).variables[0]
        assert (min(areas), max(areas)) == (limits.minimum, limits.maximum)


def test_mississippi_equals_nss_but_for_the_development_limit(ms):
    snap = _snapshot("MS")
    assert {(NSS_MS[e.region_code], e.aep) for e in ms.equations} == set(snap)
    for e in ms.equations:
        s = snap[NSS_MS[e.region_code], e.aep]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-12)
        assert e.coefficients == pytest.approx(s["coefficients"], abs=1e-12)
        assert _pct(e.sep_log) == pytest.approx(_pct(s["sep_log"]), abs=0.01)
        for v, n in zip(e.variables, s["variables"]):
            assert (v.minimum, v.maximum) == (n["minimum"], n["maximum"])
        if e.region_code == "MAP":
            assert not e.applicability and "applicability" not in s
        else:
            assert [(a.code, a.maximum) for a in e.applicability] == [("LC11DEV", 20.0)]
            nss_limit = s.get("applicability")
            if e.aep == 0.5:
                assert nss_limit[0]["maximum"] == 19.9999
            else:
                assert nss_limit is None


def test_mississippi_avp_matches_its_sep(ms):
    """SEP is printed as a whole percent and AVP to three decimals, about one percent between them."""
    for e in ms.equations:
        assert _pct(math.sqrt(e.avp)) == pytest.approx(_pct(e.sep_log), abs=1.1), (
            e.region_code,
            e.aep,
        )


# ---- The Mississippi Alluvial Plain, shared by AR, MS and LA


def test_alluvial_plain_table_3_reproduces_where_printed(ar, ms):
    """Table 3's title says 58 gages; it prints 50, all of them reproducing."""
    assert len(MAP_GAGES) == 50
    for lib, slope in ((ar, "CSL1085LFP"), (ms, "CSL10_85fm")):
        misfits = set()
        for site, gage in MAP_GAGES.items():
            printed = {**gage, slope: gage["SLOPE"]}
            for i, aep in enumerate(AEP8):
                if not _inside(lib.equation("MAP", aep), printed, gage["flows"][i], figures=3):
                    misfits.add(site)
        assert not misfits, misfits


def test_the_alluvial_plain_equation_is_the_same_in_arkansas_mississippi_and_louisiana(ar, ms):
    la = _snapshot("LA")
    for aep in AEP8:
        a, m = ar.equation("MAP", aep), ms.equation("MAP", aep)
        assert a.intercept == m.intercept and a.coefficients == m.coefficients
        assert a.avp == m.avp and a.sep_log == m.sep_log
        assert [(v.minimum, v.maximum) for v in a.variables] == [
            (v.minimum, v.maximum) for v in m.variables
        ]
        s = la["GC1888", aep]
        assert s["intercept"] == pytest.approx(a.intercept, abs=1e-12)
        assert s["coefficients"] == pytest.approx(a.coefficients, abs=1e-12)
        assert _pct(s["sep_log"]) == pytest.approx(_pct(a.sep_log), abs=0.01)
    assert [v.code for v in ar.equation("MAP", 0.01).variables][1] == "CSL1085LFP"
    assert [v.code for v in ms.equation("MAP", 0.01).variables][1] == "CSL10_85fm"
    assert [n["code"] for n in la["GC1888", 0.01]["variables"]][1] == "CSL10_85"


def test_alluvial_plain_avp_matches_its_sep(ar):
    for e in (ar.equation("MAP", a) for a in AEP8):
        assert _pct(math.sqrt(e.avp)) == pytest.approx(_pct(e.sep_log), abs=1.1), e.aep


def test_region_d_and_the_alluvial_plain_stay_side_by_side(ar):
    """Co-located in eastern Arkansas, no preference: SIR 2021-5046 does not say it supersedes
    region D, and NSS keeps D enabled for StreamStats while GC1917 is enabled for NSS alone."""
    assert ar.region_preference is None
    assert "MAP" not in ar.region_gates and "D" not in ar.region_gates
    assert "status 3" in ar.notes
