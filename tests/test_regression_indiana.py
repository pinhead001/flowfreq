"""Indiana's verified equations.

Rao (2006), FHWA/IN/JTRP-2005/18, Chapter 3, with the worked examples and gage lists of Knipe and
Rao (2005), FHWA/IN/JTRP-2005/1. Neither report prints a regression flow for each gage, so the gage
checks use the percent errors they print (Tables 3.7.1 and 3.7.2).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest
from scipy import stats

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {str(r): f"GC{1004 + r}" for r in range(1, 9)}
AEPS = [0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
#: Characteristics (Knipe and Rao Appendix A), at-site Q100 (Appendix B) and Table 3.7.1's errors.
GAGES = json.loads((ROOT / "fixtures" / "regression" / "indiana_gages.json").read_text("utf-8"))


@pytest.fixture(scope="module")
def lib():
    return load_state("IN")


def inputs(eq, gage):
    return {v.code: gage[v.code] for v in eq.variables}


def rounding_bound(eq, values, flow_digits=None):
    """Relative flow error the printed rounding allows: constant and exponents to their last digit."""
    constant = 10.0**eq.intercept
    bound = 0.05 / constant
    for v in eq.variables:
        bound += 0.0005 * abs(math.log(values[v.code] + (v.offset or 0.0)))
    if flow_digits:
        bound += 0.5 * 10.0 ** flow_digits[0] / flow_digits[1]
    return bound


def test_indiana_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == list(NSS)
    assert len(lib.equations) == 48
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
        assert lib.aeps(region) == sorted(AEPS, reverse=True)


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "IN.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-12)
        assert list(e.coefficients) == s["coefficients"]
        assert len(e.variables) == len(s["variables"])
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.minimum, v.maximum) == (n["code"], n["minimum"], n["maximum"])
            # NSS writes the printed (%U + 1) as log10_plus1; it is stored as printed.
            assert (v.offset == 1.0) == (n["transform"] == "log10_plus1")


def test_nss_region_id_limit_selects_nothing_and_is_left_out(lib):
    """HIGHREG is 1005-1012 in every region, so it cannot choose one, and no report has it."""
    snap = json.loads((DATA / "nss_snapshot" / "IN.json").read_text(encoding="utf-8"))
    limits = [s["applicability"] for s in snap["equations"] if s.get("applicability")]
    assert len(limits) == 8
    assert all(
        (a[0]["code"], a[0]["minimum"], a[0]["maximum"]) == ("HIGHREG", 1005.0, 1012.0)
        for a in limits
    )
    assert all(not e.applicability for e in lib.equations)


#: %Error, Tables 3.5.2-3.5.9, for the 10- to 500-year equations.
PERCENT_ERROR = {
    "1": [27.1, 27.8, 28.3, 29.5, 30.6, 32.2],
    "2": [35.5, 35.6, 36.0, 36.8, 37.7, 39.4],
    "3": [20.9, 19.7, 19.9, 20.4, 22.1, 24.8],
    "4": [23.1, 22.5, 22.4, 22.4, 22.7, 23.5],
    "5": [26.7, 27.7, 28.3, 28.8, 33.5, 30.0],
    "6": [37.8, 38.7, 39.4, 40.1, 40.8, 41.7],
    "7": [41.7, 44.4, 46.2, 44.3, 49.4, 51.3],
    "8": [41.0, 40.4, 39.9, 39.5, 39.1, 38.6],
}


def test_percent_error_is_the_model_error_of_equation_3_5_1(lib):
    """Every printed %Error follows from the printed variance alone, within its rounding."""
    for region, row in PERCENT_ERROR.items():
        for aep, pct in zip(AEPS, row):
            mev = lib.equation(region, aep).model_error_variance
            lo = 100 * math.sqrt(math.exp(5.3019 * (mev - 0.0005)) - 1)
            hi = 100 * math.sqrt(math.exp(5.3019 * (mev + 0.0005)) - 1)
            assert lo - 0.05 <= pct <= hi + 0.05, (region, aep)


def test_table_3_7_2_column_4_is_the_one_percent_percent_error():
    """Table 3.7.2 calls that column the average model error, from equation 3.5.1."""
    column_4 = [29.5, 36.8, 20.4, 22.4, 28.8, 40.1, 44.3, 39.5]
    assert [PERCENT_ERROR[r][3] for r in sorted(PERCENT_ERROR)] == column_4


def test_no_interval_because_no_prediction_error_is_printed(lib):
    for e in lib.equations:
        assert e.model_error_variance is not None
        assert e.equivalent_years is not None
        assert e.sep_log is None and e.avp is None and e.covariance is None
    est = evaluate(lib.equation("3", 0.01), {"CONTDA": 122.37, "CSL10_85": 10.82})
    assert est.interval is None


EXAMPLE_1 = {"CONTDA": 122.37, "CSL10_85": 10.82}
EXAMPLE_3 = {"CONTDA": 169.78, "CSL10_85": 21.64, "URBAN": 0.79}
EXAMPLE_3_URBAN = {"CONTDA": 169.78, "CSL10_85": 21.64, "URBAN": 22.67}
#: Tables 6.1, 6.5 and 6.7: flow, lower and upper limit for 10-, 25-, 50-, 100-, 200- and 500-year.
PRINTED = {
    "1": (
        "3",
        EXAMPLE_1,
        [14400, 18300, 21300, 24600, 27900, 32800],
        [9810, 12700, 14700, 16700, 18400, 20600],
        [21200, 26400, 31000, 36100, 42300, 52200],
    ),
    "3": (
        "4",
        EXAMPLE_3,
        [17800, 22700, 26300, 29800, 33400, 38200],
        [11700, 14900, 17300, 19500, 21700, 24400],
        [27300, 34500, 40100, 45500, 51500, 59700],
    ),
    "3 urban": (
        "4",
        EXAMPLE_3_URBAN,
        [21900, 27800, 32100, 36200, 40500, 46000],
        [14300, 18200, 21000, 23500, 26100, 29100],
        [33600, 42500, 49200, 55700, 62900, 72600],
    ),
}


@pytest.mark.parametrize("name", PRINTED)
def test_worked_examples_1_and_3(lib, name):
    """The spreadsheet behind the examples held more digits than the tables print."""
    region, values, flows, _, _ = PRINTED[name]
    for aep, printed in zip(AEPS, flows):
        eq = lib.equation(region, aep)
        flow = evaluate(eq, values).flow_cfs
        tol = rounding_bound(eq, values, flow_digits=(2, printed))
        assert flow == pytest.approx(printed, rel=tol), (name, aep)


@pytest.mark.parametrize(
    "contda, slope, wetland, region_7, region_8",
    [
        (5.13, 21.06, 2.41, 548, 227),
        (12.87, 12.34, 3.62, 771, 399),
        (15.5, 12.02, 3.64, 891, 465),
        (21.7, 6.58, 2.95, 1091, 691),
    ],
)
def test_worked_example_2_bigler_ditch(lib, contda, slope, wetland, region_7, region_8):
    """Tables 6.2 and 6.3: one stream, both equations, the same Q100 to 1 cfs but for rounding."""
    seven = {"CONTDA": contda, "CSL10_85": slope, "WETLAND": wetland}
    eight = {"CONTDA": contda, "WETLAND": wetland}
    assert evaluate(lib.equation("7", 0.01), seven).flow_cfs == pytest.approx(region_7, abs=1.0)
    assert evaluate(lib.equation("8", 0.01), eight).flow_cfs == pytest.approx(region_8, abs=1.0)


@pytest.mark.parametrize("name", PRINTED)
def test_printed_limits_are_90_percent_though_called_95(lib, name):
    """At 95 percent they would imply a prediction variance below the model error variance."""
    region, _, _, lower, upper = PRINTED[name]
    n = {"3": 24, "4": 72}[region]
    p = {"3": 3, "4": 4}[region]
    for aep, lo, hi in zip(AEPS, lower, upper):
        half = (math.log10(hi) - math.log10(lo)) / 2
        mev = lib.equation(region, aep).model_error_variance
        at_90 = (half / stats.t.ppf(0.95, n - p)) ** 2
        at_95 = (half / stats.t.ppf(0.975, n - p)) ** 2
        assert at_95 < mev < at_90
        assert 0.001 < at_90 - mev < 0.005


def test_gage_list_matches_the_region_counts(lib):
    """Region 7 is 21 in Appendix A and 22 in Table 3.5.1; the text's total of 223 needs 21."""
    assert len(GAGES) == 223
    for region in NSS:
        count = sum(g["region"] == region for g in GAGES.values())
        assert count == lib.equation(region, 0.01).n_sites
    assert lib.equation("7", 0.01).n_sites == 21


def test_table_3_7_1_percent_errors_reproduce(lib):
    """Printed against the full regression at 44 gages, within the rounding of Appendix A's inputs."""
    checked = [g for g in GAGES.values() if "percent_error_table_371" in g]
    assert len(checked) == 44
    for g in checked:
        eq = lib.equation(g["region"], 0.01)
        values = inputs(eq, g)
        est = evaluate(eq, values, allow_extrapolation=True).flow_cfs
        pct = 100 * abs(est - g["q100_table_371"]) / g["q100_table_371"]
        # Appendix A prints slope to 0.1 and drainage area to about 3 significant figures.
        bound = rounding_bound(eq, values)
        for v, b in zip(eq.variables, eq.coefficients):
            x = values[v.code] + (v.offset or 0.0)
            step = 0.05 if v.code != "CONTDA" else 0.5 * 10 ** (math.floor(math.log10(x)) - 2)
            bound += abs(b) * step / x
        assert pct == pytest.approx(
            g["percent_error_table_371"],
            abs=100 * bound * (est / g["q100_table_371"]) + 0.05,
        )


#: Table 3.7.2, column 3: mean absolute percent error of the full regression, 100-year.
TABLE_372_COLUMN_3 = {
    "1": 21.8,
    "2": 33.8,
    "3": 25.0,
    "4": 23.0,
    "5": 23.5,
    "6": 43.7,
    "7": 45.0,
    "8": 34.7,
}
#: Unexplained: the mean over Appendix A and B's 24 region 3 gages is 24.4, not 25.0.
REGION_3_MEAN_ERROR = 24.4


def regional_errors(lib, region):
    eq = lib.equation(region, 0.01)
    errors = []
    for g in GAGES.values():
        if g["region"] == region:
            est = evaluate(eq, inputs(eq, g), allow_extrapolation=True).flow_cfs
            errors.append(100 * abs(est - g["q100"]) / g["q100"])
    return errors


@pytest.mark.parametrize("region", sorted(set(NSS) - {"3"}))
def test_table_3_7_2_regional_errors_reproduce(lib, region):
    errors = regional_errors(lib, region)
    assert sum(errors) / len(errors) == pytest.approx(TABLE_372_COLUMN_3[region], abs=0.15)


def test_region_3_regional_error_differs_from_table_3_7_2(lib):
    errors = regional_errors(lib, "3")
    assert sum(errors) / len(errors) == pytest.approx(REGION_3_MEAN_ERROR, abs=0.05)
    assert TABLE_372_COLUMN_3["3"] == 25.0


def test_region_3_area_limit_is_as_printed_though_a_listed_gage_exceeds_it(lib):
    """Tables 3.5.10 and 4.12 and NSS print 284 mi2; Appendix A has 03303000 at 476 in region 3."""
    assert lib.equation("3", 0.01).variables[0].maximum == 284.0
    assert GAGES["03303000"]["region"] == "3" and GAGES["03303000"]["CONTDA"] == 476.0
    others = [g["CONTDA"] for k, g in GAGES.items() if g["region"] == "3" and k != "03303000"]
    assert max(others) == 248.0


def test_ranges_match_the_gage_list_but_for_the_pinned_cases(lib):
    pinned = {("1", "CSL10_85", "maximum"), ("3", "CONTDA", "maximum")}
    for region in NSS:
        gages = [g for g in GAGES.values() if g["region"] == region]
        for v in lib.equation(region, 0.01).variables:
            low = min(g[v.code] for g in gages)
            high = max(g[v.code] for g in gages)
            assert v.minimum == low
            if (region, v.code, "maximum") not in pinned:
                assert v.maximum == high
    assert max(g["CSL10_85"] for g in GAGES.values() if g["region"] == "1") == 79.5
    assert lib.equation("1", 0.01).variables[1].maximum == 79.0


def test_model_error_variances_that_break_their_trends_are_as_printed(lib):
    """Their %Error, and Table 3.7.2's column 4, agree with these values."""
    assert lib.equation("5", 0.005).model_error_variance == 0.020
    assert lib.equation("7", 0.01).model_error_variance == 0.034
