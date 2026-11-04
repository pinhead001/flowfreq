"""Louisiana's verified equations: SIR 2024-5031 Tables 3-4 with its release, and SIR 2021-5046.

Neither report prints a worked example, so the check is per gage: the release's GLS flows and
variances at 211 Coastal Plain gages, and Table 3's regression flows at 50 Alluvial Plain gages.
"""

from __future__ import annotations

import itertools
import json
import math
from decimal import Decimal
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
GAGES = json.loads((ROOT / "fixtures" / "regression" / "louisiana_gages.json").read_text("utf-8"))
CP = "Coastal Plain"
NSS = {"MAP": "GC1888", CP: "GC1889"}
SUPERSEDED = ("GC879", "GC880", "GC881")
AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
#: The release's row for this gage prints the 0.2 percent value second in both GLS and VARgls.
MISALIGNED = "08031000"


@pytest.fixture(scope="module")
def lib():
    return load_state("LA")


def _aligned(site, values):
    """A release row in AEP order, moving the misaligned gage's 0.2 percent value back."""
    return [values[0], *values[2:], values[1]] if site == MISALIGNED else values


def _half_unit(printed):
    return 0.5 * 10.0 ** Decimal(printed).as_tuple().exponent


def _inside(eq, gage, flow, half_unit, constant_slack, floor=0.05):
    """Whether a printed flow lies in what the printed inputs and coefficients allow."""
    lo = hi = eq.intercept
    slack = constant_slack
    for v, c in zip(eq.variables, eq.coefficients):
        x = float(gage[v.code])
        h = half_unit(gage[v.code])
        ends = [c * math.log10(x - h), c * math.log10(x + h)]
        lo += min(ends)
        hi += max(ends)
        slack += 5e-4 * abs(math.log10(x))
    half = max(0.5 * 10.0 ** (math.floor(math.log10(flow)) - 2), floor)
    return 10 ** (lo - slack) <= flow + half and flow - half <= 10 ** (hi + slack)


def _misfits(lib, region, half_unit, constant_slack, floor=0.05):
    out = set()
    for site, gage in GAGES[region].items():
        flows = _aligned(site, gage["flows"]) if region == CP else gage["flows"]
        for aep, flow in zip(AEPS, flows):
            eq = lib.equation(region, aep)
            if not _inside(eq, gage, flow, half_unit, constant_slack(eq), floor):
                out.add((site, aep))
    return out


def _cp_misfits(lib):
    """The release prints at most two decimals, trimming trailing zeros."""
    return _misfits(lib, CP, lambda printed: 0.005, lambda eq: 5e-4)


def _map_misfits(lib):
    """Table 2 prints the constant to a whole number, and Table 3 flows to a whole number."""
    return _misfits(
        lib, "MAP", _half_unit, lambda eq: 0.5 / 10**eq.intercept / math.log(10), floor=0.5
    )


def test_louisiana_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == [CP, "MAP"]
    assert len(lib.equations) == 16
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
        assert set(lib.aeps(region)) == set(AEPS)


def test_the_superseded_regions_are_left_out(lib):
    for code in SUPERSEDED:
        assert region_for_nss_code(lib, code) is None
    assert "superseded" in lib.notes


def test_the_coastal_plain_flows_reproduce_at_every_release_gage(lib):
    assert len(GAGES[CP]) == 211
    assert not _cp_misfits(lib)


def test_the_alluvial_plain_flows_reproduce_at_every_gage_table_3_prints(lib):
    assert len(GAGES["MAP"]) == 50
    assert not _map_misfits(lib)


def test_the_misaligned_release_row_is_misaligned_as_pinned(lib):
    gage = GAGES[CP][MISALIGNED]
    assert gage["flows"][:3] == [2540.0, 16900.0, 4600.0]
    eq = lib.equation(CP, 0.002)
    inputs = {v.code: float(gage[v.code]) for v in eq.variables}
    assert evaluate(eq, inputs).flow_cfs == pytest.approx(16900, rel=0.006)


@pytest.mark.parametrize("region", ["MAP", CP])
def test_a_small_change_to_an_exponent_is_caught(lib, region):
    """The per-gage check has to be able to fail."""
    eq = lib.equation(region, 0.01)
    gages = GAGES[region]
    half_unit = _half_unit if region == "MAP" else (lambda printed: 0.005)
    slack = 0.5 / 10**eq.intercept / math.log(10) if region == "MAP" else 5e-4

    class Moved:
        intercept = eq.intercept
        variables = eq.variables
        coefficients = [eq.coefficients[0] + 0.005, *eq.coefficients[1:]]

    i = AEPS.index(0.01)
    misses = [
        site
        for site, g in gages.items()
        if site != MISALIGNED and not _inside(Moved, g, g["flows"][i], half_unit, slack)
    ]
    assert misses


def _rows():
    return {
        site: np.array([1.0, math.log10(float(g["DRNAREA"])), math.log10(float(g["CSL10_85"]))])
        for site, g in GAGES[CP].items()
    }


def test_mev_plus_xux_is_the_releases_variance_at_every_gage(lib):
    """Equation 14. Table 4 prints U to four decimals, MEV and VARgls also to four."""
    rows = _rows()
    for i, aep in enumerate(AEPS):
        eq = lib.equation(CP, aep)
        u = np.array(eq.covariance)
        assert (u == u.T).all()
        assert (np.linalg.eigvalsh(u) > 0).all()
        for site, x in rows.items():
            printed = _aligned(site, GAGES[CP][site]["variance"])[i]
            rounding = 5e-5 * float(np.abs(np.outer(x, x)).sum()) + 1e-4
            assert abs(eq.model_error_variance + x @ u @ x - printed) <= rounding, (site, aep)


def test_the_mean_variance_is_the_printed_avp(lib):
    rows = _rows()
    for i, aep in enumerate(AEPS):
        eq = lib.equation(CP, aep)
        u = np.array(eq.covariance)
        mean = np.mean([eq.model_error_variance + x @ u @ x for x in rows.values()])
        assert mean == pytest.approx(eq.avp, abs=2e-4), aep
        assert math.sqrt(eq.avp) == pytest.approx(eq.sep_log, abs=2e-4), aep


def test_no_other_variable_order_returns_the_avp(lib):
    eq = lib.equation(CP, 0.01)
    u = np.array(eq.covariance)
    x = np.array(list(_rows().values()))
    for order in itertools.permutations(range(3)):
        z = x[:, order]
        mean = float(np.mean(eq.model_error_variance + np.einsum("ij,jk,ik->i", z, u, z)))
        if order == (0, 1, 2):
            assert mean == pytest.approx(eq.avp, abs=2e-4)
        else:
            assert abs(mean - eq.avp) > 2e-3, order


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "LA.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    kept = {key for key in by_key if key[0] not in SUPERSEDED}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == kept
    assert len(by_key) - len(kept) == 22
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-12)
        assert list(e.coefficients) == s["coefficients"]
        assert e.sep_log == pytest.approx(s["sep_log"], abs=1e-6)
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.minimum, v.maximum) == (n["code"], n["minimum"], n["maximum"])


def test_the_coastal_plain_keeps_the_ranges_of_page_9(lib):
    """The Summary says 0.011 and 1.16; no release gage is at 1.16, and the lowest is 1.2."""
    eq = lib.equation(CP, 0.01)
    assert [(v.minimum, v.maximum) for v in eq.variables] == [(0.01, 1899.0), (1.2, 273.6)]
    assert min(float(g["CSL10_85"]) for g in GAGES[CP].values()) == 1.2
    gage = GAGES[CP]["07372200"]
    implied = (gage["flows"][5] / 10**eq.intercept / float(gage["DRNAREA"]) ** 0.832) ** (1 / 0.467)
    assert implied == pytest.approx(1.2, abs=0.01)


def test_the_alluvial_plain_has_no_covariance_so_intervals_use_its_avp(lib):
    eq = lib.equation("MAP", 0.01)
    assert eq.covariance is None
    assert eq.avp == 0.049
    est = evaluate(eq, {"DRNAREA": 309.0, "CSL10_85": 1.10, "BSHAPELFP": 6.38})
    assert est.flow_cfs == pytest.approx(378 * 309.0**0.683 * 1.10**0.455 * 6.38**-0.326)
