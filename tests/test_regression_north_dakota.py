"""North Dakota's verified equations (SIR 2015-5096, Tables 1, 2, 4 and 5)."""

from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {"A": "GC138", "B": "GC139", "C": "GC140"}
AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]
#: Characteristics and regression flows from Table 1, at the 230 gages that print the flows.
GAGES = json.loads(
    (ROOT / "fixtures" / "regression" / "north_dakota_gages.json").read_text("utf-8")
)
#: Table 4's printed standard errors of prediction, in percent, in AEP order.
SP = {
    "A": [57.55, 51.03, 51.09, 52.11, 53.85, 55.02, 58.84],
    "B": [75.32, 60.30, 57.95, 58.78, 60.69, 63.53, 70.14],
    "C": [82.90, 59.45, 54.24, 52.85, 54.05, 55.23, 59.85],
}
#: The report's text prints each zone's range of SP, a second place to check Table 4.
SP_RANGE_IN_TEXT = {"A": (51.03, 58.84), "B": (57.95, 75.32), "C": (52.85, 82.90)}
#: Table 4's limits, and NSS's, which round the same values outward.
LIMITS = {
    "A": {"DRNAREA": (4.328, 3393.031), "CSL1085LFP": (1.016, 17.522)},
    "B": {"DRNAREA": (0.110, 8342.387), "RUGGED": (68.018, 7819.737), "COMPRAT": (1.414, 3.480)},
    "C": {"DRNAREA": (0.132, 2811.637), "RUGGED": (21.309, 2264.002)},
}
NSS_LIMITS = {
    "A": {"DRNAREA": (4.32, 3394.0), "CSL1085LFP": (1.01, 17.6)},
    "B": {"DRNAREA": (0.11, 8343.0), "RUGGED": (68.0, 7820.0), "COMPRAT": (1.4, 3.48)},
    "C": {"DRNAREA": (0.13, 2812.0), "RUGGED": (21.3, 2265.0)},
}
#: Table 5, with the average it prints.
TABLE_5 = {
    "A": ([0.499, 0.525, 0.544, 0.566, 0.582, 0.597, 0.628], 0.563),
    "B": ([0.609, 0.597, 0.591, 0.587, 0.584, 0.581, 0.577], 0.589),
    "C": ([0.438, 0.475, 0.493, 0.511, 0.521, 0.529, 0.543], 0.501),
}


@pytest.fixture(scope="module")
def lib():
    return load_state("ND")


def test_north_dakota_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["A", "B", "C"]
    assert len(lib.equations) == 21
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
        assert lib.aeps(region) == sorted(AEPS, reverse=True)


def test_published_regression_flows_reproduce_at_every_gage(lib):
    """Table 1 prints them unrounded, computed with Table 4's rounded coefficients."""
    assert len(GAGES) == 230
    worst = 0.0
    for gage in GAGES.values():
        for aep, flow in zip(AEPS, gage["flows"]):
            eq = lib.equation(gage["region"], aep)
            est = evaluate(
                eq, {v.code: gage[v.code] for v in eq.variables}, allow_extrapolation=True
            )
            worst = max(worst, abs(est.flow_cfs / flow - 1))
    assert worst < 3e-5


@pytest.mark.parametrize("region", ["A", "B", "C"])
@pytest.mark.parametrize("aep", AEPS)
def test_fitting_table_1_recovers_table_4(lib, region, aep):
    """A second entry of every coefficient, read from Table 1 and not from Table 4."""
    eq = lib.equation(region, aep)
    gages = [g for g in GAGES.values() if g["region"] == region]
    x = np.array([[1.0] + [v.apply(g[v.code]) for v in eq.variables] for g in gages])
    y = np.log10([g["flows"][AEPS.index(aep)] for g in gages])
    fit = np.linalg.lstsq(x, y, rcond=None)[0]
    assert fit == pytest.approx([eq.intercept, *eq.coefficients], abs=1e-6)


def test_zone_b_ruggedness_limit_is_the_value_regressed(lib):
    """Table 1 prints 7819.70 at 06337000, and its regression flow implies Table 4's 7819.737."""
    gage = GAGES["06337000"]
    assert gage["RUGGED"] == pytest.approx(7819.7004, abs=1e-3)
    eq = lib.equation("B", 0.01)
    rest = (
        eq.intercept
        + eq.coefficients[0] * math.log10(gage["DRNAREA"])
        + eq.coefficients[2] * gage["COMPRAT"]
    )
    implied = 10 ** ((math.log10(gage["flows"][5]) - rest) / eq.coefficients[1])
    assert implied == pytest.approx(7819.737, abs=1e-3)
    assert eq.variables[1].maximum == 7819.737


def test_every_equation_equals_the_nss_snapshot_but_for_its_limits(lib):
    snap = json.loads((DATA / "nss_snapshot" / "ND.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == s["intercept"]
        assert list(e.coefficients) == s["coefficients"]
        # NSS serves SP rounded to 0.1 percent.
        assert e.sep_log == pytest.approx(s["sep_log"], abs=2e-4)
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.transform) == (n["code"], n["transform"])
            assert (v.minimum, v.maximum) == LIMITS[e.region_code][v.code]
            assert (n["minimum"], n["maximum"]) == NSS_LIMITS[e.region_code][v.code]


@pytest.mark.parametrize("region", ["A", "B", "C"])
def test_error_statistics_agree_inside_table_4(lib, region):
    equations = [lib.equation(region, aep) for aep in AEPS]
    assert (min(SP[region]), max(SP[region])) == SP_RANGE_IN_TEXT[region]
    for eq, sp in zip(equations, SP[region]):
        # The report's equation 1, from AVP back to the printed SP.
        assert 100 * math.sqrt(10 ** (math.log(10) * eq.avp) - 1) == pytest.approx(sp, abs=0.006)
        assert eq.sep_log == pytest.approx(math.sqrt(eq.avp), abs=2e-5)
        assert 0 < eq.model_error_variance < eq.avp
        assert eq.covariance is None


def test_intervals_use_the_average_variance_of_prediction(lib):
    """With no covariance, the 90-percent interval is Student's t on n - p - 1 times sqrt(AVP)."""
    eq = lib.equation("B", 0.01)
    est = evaluate(eq, {"DRNAREA": 556.6, "RUGGED": 2096.5, "COMPRAT": 2.459})
    factor = 10 ** (stats.t.ppf(0.95, 123 - 3 - 1) * math.sqrt(eq.avp))
    assert est.interval[1] / est.flow_cfs == pytest.approx(factor, rel=1e-9)


def test_site_counts_reach_the_231_the_report_prints(lib):
    """Table 1's flows put 38, 123 and 69 gages in the zones; 06469600 is zone C's 70th."""
    assert "06469600" not in GAGES
    assert Counter(g["region"] for g in GAGES.values()) == {"A": 38, "B": 123, "C": 69}
    counts = {r: lib.equation(r, 0.01).n_sites for r in NSS}
    assert counts == {"A": 38, "B": 123, "C": 70}
    assert sum(counts.values()) == 231


@pytest.mark.parametrize("region", ["A", "B", "C"])
def test_table_5_is_consistent_and_kept_in_the_notes(lib, region):
    """The schema has no field for an area-ratio exponent, so the notes carry it."""
    exponents, average = TABLE_5[region]
    assert sum(exponents) / 7 == pytest.approx(average, abs=5e-4)
    text = ", ".join(f"{b:.3f}" for b in exponents)
    assert f"zone {region} {text} (average {average:.3f})" in lib.notes
