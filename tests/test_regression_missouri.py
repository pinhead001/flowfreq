"""Missouri's verified equations (SIR 2014-5165 Tables 7-10, 14; SIR 2010-5073 Tables 5, 7-9).

Neither report has a worked example, so the checks are per gage: the regression columns of each
report's Table 4 and the regression variances of the 2014 report's Table 12.
"""

from __future__ import annotations

import json
import math
from itertools import permutations
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {"1": "GC1483", "2": "GC1484", "3": "GC1485", "Urban": "GC1486"}
#: Inputs (2014 Table 5, 2010 Table 2), regression flows (Table 4), variances (2014 Table 12).
GAGES = json.loads((ROOT / "fixtures" / "regression" / "missouri_gages.json").read_text("utf-8"))
RURAL_AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
URBAN_AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]


@pytest.fixture(scope="module")
def lib():
    return load_state("MO")


def _rural_x(gage):
    bshape = [math.log10(gage["BSHAPE"])] if "BSHAPE" in gage else []
    return np.array([1.0, math.log10(gage["DRNAREA"]), *bshape])


def _urban_x(gage):
    return np.array([1.0, math.log10(gage["DRNAREA"]), gage["IMPNLCD01"]])


def test_missouri_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["1", "2", "3", "Urban"]
    assert len(lib.equations) == 31
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
    assert lib.aeps("1") == lib.aeps("2") == lib.aeps("3") == RURAL_AEPS
    assert lib.aeps("Urban") == URBAN_AEPS


def test_rural_regression_flows_reproduce_at_every_gage(lib):
    """Table 4's RRE column from Table 5's rounded inputs, the rounding the report used.

    The largest difference is 0.69 percent, a printed 71 cfs against 70.5 at 06921400.
    """
    assert len(GAGES["rural"]) == 278
    assert sum(g["region"] == "1" for g in GAGES["rural"].values()) == 131
    assert sum(g["region"] == "2" for g in GAGES["rural"].values()) == 135
    assert sum(g["region"] == "3" for g in GAGES["rural"].values()) == 12
    for gage in GAGES["rural"].values():
        for i, aep in enumerate(RURAL_AEPS):
            eq = lib.equation(gage["region"], aep)
            inputs = {v.code: gage[v.code] for v in eq.variables}
            est = evaluate(eq, inputs, allow_extrapolation=True)
            assert est.flow_cfs == pytest.approx(gage["flows"][i], rel=0.007)


def test_urban_regression_flows_reproduce_at_every_gage(lib):
    """Table 4's R column, three significant figures."""
    assert len(GAGES["urban"]) == 35
    for gage in GAGES["urban"].values():
        for i, aep in enumerate(URBAN_AEPS):
            eq = lib.equation("Urban", aep)
            est = evaluate(eq, {"DRNAREA": gage["DRNAREA"], "IMPNLCD01": gage["IMPNLCD01"]})
            assert est.flow_cfs == pytest.approx(gage["flows"][i], rel=0.006)


@pytest.mark.parametrize(
    "gage,flow",
    [
        ("06894000", 30100),
        ("06910230", 14300),
        ("06936475", 18000),
        ("07010022", 7470),
        ("07010035", 1940),
        ("07052000", 12700),
        ("07063200", 519),
    ],
)
def test_urban_table_11(lib, gage, flow):
    g = GAGES["urban"][gage]
    est = evaluate(
        lib.equation("Urban", 0.01), {"DRNAREA": g["DRNAREA"], "IMPNLCD01": g["IMPNLCD01"]}
    )
    assert est.flow_cfs == pytest.approx(flow, rel=0.005)


def test_rural_covariance_reproduces_table_12_at_every_gage(lib):
    """MEV + xUx against the printed RRE variance. Both MEV and the variance print to 0.001."""
    for gage in GAGES["rural"].values():
        x = _rural_x(gage)
        for i, aep in enumerate(RURAL_AEPS):
            eq = lib.equation(gage["region"], aep)
            variance = eq.model_error_variance + float(x @ np.array(eq.covariance) @ x)
            assert variance == pytest.approx(gage["variance"][i], abs=1e-3)


def test_only_the_stated_order_reproduces_table_12():
    """Region 1 settles it: the next best order misses by 0.0027, the others by 0.055 or more."""
    gages = [g for g in GAGES["rural"].values() if g["region"] == "1"]
    lib = load_state("MO")
    worst = {}
    for perm in permutations(range(3)):
        misses = []
        for gage in gages:
            x = _rural_x(gage)[list(perm)]
            for i, aep in enumerate(RURAL_AEPS):
                eq = lib.equation("1", aep)
                u = np.array(eq.covariance)
                misses.append(abs(eq.model_error_variance + x @ u @ x - gage["variance"][i]))
        worst[perm] = max(misses)
    assert worst[(0, 1, 2)] < 1e-3
    assert all(miss > 2.5e-3 for perm, miss in worst.items() if perm != (0, 1, 2))


@pytest.mark.parametrize("region", ["1", "2", "3"])
def test_rural_covariance_reproduces_printed_avp(lib, region):
    gages = [g for g in GAGES["rural"].values() if g["region"] == region]
    for aep in RURAL_AEPS:
        eq = lib.equation(region, aep)
        u = np.array(eq.covariance)
        mean = np.mean([eq.model_error_variance + _rural_x(g) @ u @ _rural_x(g) for g in gages])
        assert mean == pytest.approx(eq.avp, abs=1e-3)


def test_urban_covariance_reproduces_printed_average(lib):
    """Table 8 prints six decimals, coarse for the impervious term, and Table 7 three."""
    for aep in URBAN_AEPS:
        eq = lib.equation("Urban", aep)
        u = np.array(eq.covariance)
        gages = GAGES["urban"].values()
        mean = np.mean([eq.model_error_variance + _urban_x(g) @ u @ _urban_x(g) for g in gages])
        assert mean == pytest.approx(eq.avp, abs=1e-3)


def test_urban_table_7_prints_variances_labelled_logarithm_units(lib):
    """0.013 "logarithm units" is printed as -23.1 to +30.0 percent, which is its square root."""
    root = math.sqrt(lib.equation("Urban", 0.5).avp)
    assert 100 * (10**-root - 1) == pytest.approx(-23.1, abs=0.05)
    assert 100 * (10**root - 1) == pytest.approx(30.0, abs=0.05)


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "MO.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-9)
        assert list(e.coefficients) == s["coefficients"]
        # The rural SEPs print in percent. NSS's urban SEPs are inside Table 7's three-decimal
        # variances (0.09986 squared rounds to the printed 0.010), so they differ by up to 1.5e-4.
        tolerance = 1.5e-4 if e.region_code == "Urban" else 1e-6
        assert e.sep_log == pytest.approx(s["sep_log"], abs=tolerance)
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.transform, v.minimum, v.maximum) == (
                n["code"],
                n["transform"],
                n["minimum"],
                n["maximum"],
            )


def test_interval_uses_the_covariance(lib):
    eq = lib.equation("2", 0.01)
    x = np.array([1.0, math.log10(150.0), math.log10(8.0)])
    est = evaluate(eq, {"DRNAREA": 150.0, "BSHAPE": 8.0})
    half = stats.t.ppf(0.95, 132) * math.sqrt(
        eq.model_error_variance + x @ np.array(eq.covariance) @ x
    )
    assert est.interval[1] == pytest.approx(est.flow_cfs * 10**half, rel=1e-9)


def test_printed_t_against_n_sites(lib):
    """Regions 1, 2 and Urban agree. Region 3 prints 1.782, which is t at 12 degrees of freedom
    where n - p is 10, so flowfreq's region 3 intervals use 1.812 rather than the printed value."""
    for region, printed, n_p in [("1", 1.66, 128), ("2", 1.66, 132), ("Urban", 1.694, 32)]:
        assert stats.t.ppf(0.95, n_p) == pytest.approx(printed, abs=0.005)
        eq = lib.equation(region, 0.01)
        assert eq.n_sites - len(eq.variables) - 1 == n_p
    eq = lib.equation("3", 0.01)
    assert eq.n_sites == 12
    assert stats.t.ppf(0.95, 12) == pytest.approx(1.782, abs=5e-4)
    assert stats.t.ppf(0.95, eq.n_sites - len(eq.variables) - 1) == pytest.approx(1.812, abs=5e-4)


def test_urban_table_9_mean_impervious_is_a_misprint():
    """Table 9 prints means of 28.8 mi2 and 6.5 percent; Table 2 averages 28.8 and 26.5, and
    Table 10's sensitivities (+18.6 percent at 50 percent AEP for +20 percent) need 26.5."""
    gages = GAGES["urban"].values()
    assert np.mean([g["DRNAREA"] for g in gages]) == pytest.approx(28.8, abs=0.05)
    assert np.mean([g["IMPNLCD01"] for g in gages]) == pytest.approx(26.5, abs=0.05)
    assert 100 * (10 ** (0.014 * 0.2 * 26.5) - 1) == pytest.approx(18.6, abs=0.1)


def test_extrapolation_is_flagged_past_the_ranges(lib):
    est = evaluate(lib.equation("3", 0.01), {"DRNAREA": 3000.0}, allow_extrapolation=True)
    assert est.extrapolated
    with pytest.raises(ValueError):
        evaluate(lib.equation("Urban", 0.01), {"DRNAREA": 200.0, "IMPNLCD01": 10.0})
