"""Ohio's verified equations.

SIR 2019-5018 (Tables 1, 3, 4, 5, appendix tables 1.1 and 1.8, equations 12-15, and the
report's model archive) and WSP 2432 (Tables 6 and 8 and its worked example).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {
    "A simple": "GC1785",
    "B simple": "GC1786",
    "C simple": "GC1787",
    "A full": "GC1788",
    "B full": "GC1789",
    "C full": "GC1790",
    "Urban": "GC1235",
}
AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]
GAGES = json.loads((ROOT / "fixtures" / "regression" / "ohio_gages.json").read_text("utf-8"))

#: AEP: a, b, c, d, MEV, VPavg, SEP as Table 4 prints them.
TABLE_4 = {
    0.5: (1.894, 0.159, 0.397, 0.705, 0.031, 0.031, 42.40),
    0.2: (2.124, 0.182, 0.410, 0.672, 0.028, 0.029, 40.49),
    0.1: (2.243, 0.197, 0.418, 0.655, 0.030, 0.030, 41.84),
    0.04: (2.368, 0.214, 0.429, 0.639, 0.031, 0.032, 43.10),
    0.02: (2.448, 0.226, 0.438, 0.629, 0.033, 0.034, 44.61),
    0.01: (2.518, 0.237, 0.446, 0.620, 0.034, 0.035, 45.20),
    0.002: (2.657, 0.261, 0.466, 0.604, 0.035, 0.037, 46.51),
}
#: AEP: a, b, c, d, e, f, MEV, VPavg, SEP as Table 5 prints them.
TABLE_5 = {
    0.5: (1.668, 0.126, 0.326, 0.781, 0.162, -0.130, 0.027, 0.028, 40.15),
    0.2: (1.845, 0.144, 0.320, 0.767, 0.202, -0.174, 0.024, 0.024, 37.16),
    0.1: (1.939, 0.157, 0.319, 0.760, 0.222, -0.197, 0.024, 0.025, 37.62),
    0.04: (2.040, 0.172, 0.321, 0.752, 0.240, -0.222, 0.024, 0.026, 38.08),
    0.02: (2.106, 0.183, 0.324, 0.748, 0.251, -0.238, 0.025, 0.027, 37.79),
    0.01: (2.165, 0.194, 0.328, 0.743, 0.259, -0.251, 0.026, 0.028, 39.60),
    0.002: (2.285, 0.219, 0.340, 0.734, 0.272, -0.277, 0.027, 0.028, 40.27),
}
#: The bracket values printed under each table, regions A, B, C.
MANUAL_SIMPLE = {
    0.5: (112.870, 78.295, 195.427),
    0.2: (202.416, 132.977, 341.680),
    0.1: (275.548, 175.038, 458.759),
    0.04: (382.469, 233.610, 627.901),
    0.02: (472.189, 280.585, 768.952),
    0.01: (569.494, 329.858, 921.547),
    0.002: (827.547, 453.687, 1325.772),
}
MANUAL_FULL = {
    0.5: (62.267, 46.553, 98.542),
    0.2: (97.551, 69.945, 146.193),
    0.1: (124.558, 86.822, 181.115),
    0.04: (162.920, 109.576, 229.548),
    0.02: (194.552, 127.520, 269.046),
    0.01: (228.645, 146.165, 311.296),
    0.002: (318.795, 192.726, 422.053),
}
#: Table 3: DA max and min, SL max and min, W max, mean DA, mean SL, mean W, n.
TABLE_3 = {
    "A": (5989, 0.04, 516, 1.53, 25.35, 195, 48.2, 2.31, 175),
    "B": (6309, 0.04, 457, 1.21, 7.10, 242, 25.2, 1.44, 68),
    "C": (2514, 0.26, 131, 3.24, 1.23, 170, 37.1, 0.33, 32),
}
#: Where appendix table 1.1 differs from the WREG inputs behind Table 1, beyond rounding.
TABLE_1_1_DIFFERS = {
    "03022500": {"DRNAREA": 629.0},
    "03106000": {"DRNAREA": 354.0},
    "03112000": {"CSL1085LFP": 9.84},
    "03113700": {"CSL1085LFP": 127.0},
    "03114500": {"CSL1085LFP": 3.65},
    "03114550": {"CSL1085LFP": 131.0},
    "03114600": {"CSL1085LFP": 126.0},
    "03154500": {"DRNAREA": 79.7, "CSL1085LFP": 4.86, "LC92STOR": 0.06},
    "03155500": {"LC92STOR": 0.15},
    "03201480": {"CSL1085LFP": 118.0},
    "03204500": {"CSL1085LFP": 4.74},
    "03216563": {"DRNAREA": 0.91},
    "03238030": {"DRNAREA": 1.89},
    "03274880": {"DRNAREA": 0.709},
    "03275900": {"DRNAREA": 5.42},
    "04099060": {"DRNAREA": 1.39, "CSL1085LFP": 18.4},
    "04177720": {"DRNAREA": 37.5},
    "04179500": {"DRNAREA": 87.5},
}


@pytest.fixture(scope="module")
def lib():
    return load_state("OH")


def _inputs(gage, **overrides):
    values = {
        "OHREGA": float(gage["region"] == "A"),
        "OHREGC": float(gage["region"] == "C"),
        "DRNAREA": gage["DRNAREA"],
        "CSL1085LFP": gage["CSL1085LFP"],
        "LC92STOR": gage["LC92STOR"],
    }
    values.update(overrides)
    return values


def _x_row(gage):
    return np.array(
        [
            1.0,
            float(gage["region"] == "A"),
            float(gage["region"] == "C"),
            math.log10(gage["DRNAREA"]),
            math.log10(gage["CSL1085LFP"]),
            math.log10(gage["LC92STOR"] + 1.0),
        ]
    )


def _sep_percent(sep_log):
    return 100.0 * math.sqrt(math.exp((sep_log * math.log(10.0)) ** 2) - 1.0)


def test_ohio_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == sorted(NSS)
    assert len(lib.equations) == 48
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region


def test_the_superseded_2006_set_is_left_out(lib):
    assert region_for_nss_code(lib, "GC1233") is None
    assert region_for_nss_code(lib, "GC1234") is None
    assert "SIR 2006-5312" in lib.notes
    assert not any("2006" in e.citation.publication for e in lib.equations)


@pytest.mark.parametrize("aep", AEPS)
def test_stored_coefficients_round_to_the_printed_tables(lib, aep):
    """The archive's unrounded fit is stored, and it rounds to Tables 4 and 5."""
    for region in "ABC":
        full = lib.equation(f"{region} full", aep)
        coefficients = [full.coefficients[i] for i in range(5)]
        printed = TABLE_5[aep]
        assert round(full.intercept, 3) == printed[0]
        assert [round(c, 3) for c in coefficients] == list(printed[1:6])
        assert round(full.model_error_variance + 1e-9, 3) == printed[6]
        assert round(full.avp + 1e-9, 3) == printed[7]
        simple = lib.equation(f"{region} simple", aep)
        printed = TABLE_4[aep]
        assert round(simple.intercept, 3) == printed[0]
        assert [round(c, 3) for c in simple.coefficients] == list(printed[1:4])
        assert round(simple.model_error_variance + 1e-9, 3) == printed[4]
        assert round(simple.avp + 1e-9, 3) == printed[5]
        assert _sep_percent(simple.sep_log) == pytest.approx(printed[6], abs=5e-3)
        assert simple.n_sites == full.n_sites == 275


@pytest.mark.parametrize("aep", AEPS)
def test_manual_computation_coefficients_reproduce(lib, aep):
    """Only the unrounded coefficients give these six figures."""
    for model, table in (("full", MANUAL_FULL), ("simple", MANUAL_SIMPLE)):
        for i, (region, printed) in enumerate(zip("ABC", table[aep])):
            eq = lib.equation(f"{region} {model}", aep)
            bracket = 10.0 ** (
                eq.intercept
                + (eq.coefficients[0] if region == "A" else 0.0)
                + (eq.coefficients[1] if region == "C" else 0.0)
            )
            assert bracket == pytest.approx(printed, abs=5.01e-4), (region, i)


def test_worked_example_mill_creek(lib):
    """Equations 14-15: an ungaged site near gage 03220000, region A, the 1-percent AEP."""
    eq = lib.equation("A full", 0.01)
    inputs = {
        "OHREGA": 1.0,
        "OHREGC": 0.0,
        "DRNAREA": 167.0,
        "CSL1085LFP": 4.83,
        "LC92STOR": 0.77,
    }
    flow = evaluate(eq, inputs).flow_cfs
    assert round(flow) == 13357
    # Equation 14 prints four decimals, which alone give 13,354, and the sentence
    # before it says 13,562; the equation's result, 13,357, carries the example.
    printed = 10 ** (2.1648 + 0.1943) * 167**0.7431 * 4.83**0.2588 * 1.77**-0.2508
    assert round(printed) == 13354
    ratio = 1.2411
    assert 17500 / 14100 == pytest.approx(ratio, abs=5e-5)
    adjusted = 13357 * (ratio - 2 * (178 - 167) * (ratio - 1) / 178)
    assert round(adjusted) == 16179
    assert float(f"{adjusted:.3g}") == 16200
    bellepoint = GAGES["regression"]["03220000"]
    assert bellepoint["flows"][AEPS.index(0.01)] == 14100
    assert evaluate(eq, _inputs(bellepoint)).flow_cfs == pytest.approx(14100, rel=0.005)


def test_table_1_regression_flows_reproduce_at_every_regression_gage(lib):
    """1,925 flows, all equal to the three figures Table 1 prints."""
    assert len(GAGES["regression"]) == 275
    cells = 0
    for gage in GAGES["regression"].values():
        for i, aep in enumerate(AEPS):
            eq = lib.equation(f"{gage['region']} full", aep)
            flow = evaluate(eq, _inputs(gage), allow_extrapolation=True).flow_cfs
            assert float(f"{flow:.3g}") == gage["flows"][i]
            cells += 1
    assert cells == 1925


def test_table_1_limits_check_the_covariance_at_every_regression_gage(lib):
    """Table 1's 95-percent limits, from t(269) and the site's MEV + xUx.

    Three figures at 100 or more, whole numbers below 100: half a unit in the last place,
    plus up to 0.15 percent from the four-decimal model error variance and printed t.
    """
    for site, gage in GAGES["regression"].items():
        for i, aep in enumerate(AEPS):
            eq = lib.equation(f"{gage['region']} full", aep)
            est = evaluate(eq, _inputs(gage), allow_extrapolation=True, interval_level=0.95)
            assert est.interval_method.startswith("site-specific")
            for value, printed in zip(est.interval, (gage["lower"][i], gage["upper"][i])):
                unit = 1.0 if printed < 1000 else 10 ** (math.floor(math.log10(printed)) - 2)
                assert abs(value - printed) <= 0.5 * unit + 0.0015 * printed, (site, aep)


def test_gages_outside_the_regression_use_the_average_variance(lib):
    """58 gages from table 1.1's rounded inputs; limits use VPavg, as the report says."""
    assert len(GAGES["table_1_1_only"]) == 58
    for site, gage in GAGES["table_1_1_only"].items():
        for i, aep in enumerate(AEPS):
            eq = lib.equation(f"{gage['region']} full", aep)
            flow = evaluate(eq, _inputs(gage), allow_extrapolation=True).flow_cfs
            assert flow == pytest.approx(gage["flows"][i], rel=0.0065), (site, aep)
            factor = 10 ** (1.969 * math.sqrt(eq.avp))
            assert flow / factor == pytest.approx(gage["lower"][i], rel=0.0075), (site, aep)
            assert flow * factor == pytest.approx(gage["upper"][i], rel=0.0075), (site, aep)


@pytest.mark.parametrize("aep", AEPS)
def test_covariance_gives_the_printed_average_error(lib, aep):
    """``sqrt(mean(MEV + xUx))`` over the 275 gages.

    The model error variance carries four decimals, so the average is good to about 2e-4.
    """
    eq = lib.equation("A full", aep)
    u = np.array(eq.covariance)
    assert u.shape == (6, 6)
    assert np.allclose(u, u.T, rtol=1e-6, atol=0)
    variance = [
        eq.model_error_variance + float(_x_row(g) @ u @ _x_row(g))
        for g in GAGES["regression"].values()
    ]
    assert float(np.mean(variance)) == pytest.approx(eq.avp, abs=1e-4)
    assert math.sqrt(float(np.mean(variance))) == pytest.approx(eq.sep_log, abs=2e-4)


def test_covariance_order_from_the_archive_standard_errors(lib):
    """coefs100.txt: the standard errors by name, against the matrix's diagonal."""
    eq = lib.equation("A full", 0.01)
    assert [v.code for v in eq.variables] == [
        "OHREGA",
        "OHREGC",
        "DRNAREA",
        "CSL1085LFP",
        "LC92STOR",
    ]
    standard_errors = [0.08745191, 0.03336044, 0.05255852, 0.02335179, 0.05050128, 0.04662956]
    assert np.sqrt(np.diag(eq.covariance)) == pytest.approx(standard_errors, rel=1e-5)


def test_simple_model_has_no_covariance_and_uses_the_average_variance(lib):
    eq = lib.equation("A simple", 0.01)
    assert eq.covariance is None
    inputs = {"OHREGA": 1.0, "OHREGC": 0.0, "DRNAREA": 167.0}
    est = evaluate(eq, inputs)
    assert est.interval_method == "average variance of prediction"
    assert est.flow_cfs == pytest.approx(569.494 * 167**0.620287, rel=1e-5)


def test_table_5_sep_at_aep_0_02_is_a_misprint_and_nss_serves_it(lib):
    """37.79 percent is WREG's standard model error; the average SEP is 38.81.

    Table 5's own VPavg of 0.027 at that row agrees with 38.81, not with 37.79, and
    the row would otherwise fall below both of its neighbours.
    """
    snap = json.loads((DATA / "nss_snapshot" / "OH.json").read_text(encoding="utf-8"))
    served = next(
        s["sep_log"] for s in snap["equations"] if s["region_code"] == "GC1788" and s["aep"] == 0.02
    )
    assert _sep_percent(served) == pytest.approx(37.79, abs=5e-3)
    for region in "ABC":
        eq = lib.equation(f"{region} full", 0.02)
        assert _sep_percent(eq.sep_log) == pytest.approx(38.81, abs=5e-3)
        assert eq.avp == 0.0265
        assert _sep_percent(math.sqrt(eq.avp)) == pytest.approx(38.81, abs=0.05)
    for aep, printed in TABLE_5.items():
        if aep != 0.02:
            eq = lib.equation("A full", aep)
            assert _sep_percent(eq.sep_log) == pytest.approx(printed[8], abs=5e-3)


def test_table_3_limits_and_summary_statistics(lib):
    """Table 3 from the 275 gages, and in every equation of the region."""
    for region, printed in TABLE_3.items():
        da_max, da_min, sl_max, sl_min, w_max, da_mean, sl_mean, w_mean, n = printed
        gages = [g for g in GAGES["regression"].values() if g["region"] == region]
        assert len(gages) == n
        area = [g["DRNAREA"] for g in gages]
        slope = [g["CSL1085LFP"] for g in gages]
        storage = [g["LC92STOR"] for g in gages]
        assert (round(max(area)), round(min(area), 2)) == (da_max, da_min)
        assert (round(max(slope)), round(min(slope), 2)) == (sl_max, sl_min)
        assert (round(max(storage), 2), round(min(storage), 2)) == (w_max, 0.0)
        assert round(float(np.mean(area))) == da_mean
        assert round(float(np.mean(slope)), 1) == sl_mean
        assert round(float(np.mean(storage)), 2) == w_mean
        for aep in AEPS:
            variables = {v.code: v for v in lib.equation(f"{region} full", aep).variables}
            assert (variables["DRNAREA"].maximum, variables["DRNAREA"].minimum) == (
                da_max,
                da_min,
            )
            assert (variables["CSL1085LFP"].maximum, variables["CSL1085LFP"].minimum) == (
                sl_max,
                sl_min,
            )
            assert variables["LC92STOR"].maximum == w_max
            assert variables["LC92STOR"].offset == 1.0
            simple = {v.code: v for v in lib.equation(f"{region} simple", aep).variables}
            assert (simple["DRNAREA"].maximum, simple["DRNAREA"].minimum) == (da_max, da_min)


def test_table_1_was_computed_from_the_archive_inputs_not_table_1_1(lib):
    """Table 1.1's value is printed to its own figures, and it is not what Table 1 used."""
    for site, printed in TABLE_1_1_DIFFERS.items():
        gage = GAGES["regression"][site]
        for code, value in printed.items():
            assert float(f"{gage[code]:.3g}") != value, (site, code)
    gage = GAGES["regression"]["03274880"]
    assert gage["DRNAREA"] == pytest.approx(1.418, abs=1e-3)
    assert gage["flows"][0] == 215
    eq = lib.equation("C full", 0.5)
    assert evaluate(eq, _inputs(gage)).flow_cfs == pytest.approx(215, rel=0.005)
    halved = _inputs(gage, DRNAREA=0.709)
    assert evaluate(eq, halved).flow_cfs == pytest.approx(125.3, abs=0.1)


def test_urban_worked_example_toledo(lib):
    """WSP 2432: 0.89 mi2, 31.6 in., BDF of 9."""
    inputs = {"DRNAREA": 0.89, "PRECIP": 31.6, "BDF": 9.0}
    assert round(evaluate(lib.equation("Urban", 0.04), inputs).flow_cfs) == 204
    assert round(evaluate(lib.equation("Urban", 0.01), inputs).flow_cfs) == 265


def test_urban_table_8_and_its_limits(lib):
    printed = {
        0.5: (155, 0.68, 0.50, -0.50, 34.3),
        0.2: (200, 0.71, 0.63, -0.44, 34.8),
        0.1: (228, 0.74, 0.68, -0.41, 36.0),
        0.04: (265, 0.76, 0.72, -0.37, 37.6),
        0.02: (293, 0.78, 0.74, -0.35, 38.8),
        0.01: (321, 0.79, 0.76, -0.33, 40.1),
    }
    assert lib.aeps("Urban") == list(printed)
    for aep, (constant, b, c, d, sep) in printed.items():
        eq = lib.equation("Urban", aep)
        assert 10**eq.intercept == pytest.approx(constant, rel=1e-12)
        assert list(eq.coefficients) == [b, c, d]
        assert _sep_percent(eq.sep_log) == pytest.approx(sep, abs=5e-3)
        assert eq.covariance is None and eq.model_error_variance is None
        assert eq.n_sites == 30
        area, precip, bdf = eq.variables
        assert (area.code, area.minimum, area.maximum) == ("DRNAREA", 0.026, 4.09)
        assert (precip.code, precip.minimum, precip.maximum) == ("PRECIP", 31.5, 41.2)
        assert (precip.offset, precip.scale) == (-30.0, 1.0)
        assert (bdf.code, bdf.minimum, bdf.maximum) == ("BDF", 0.0, 12.0)
        assert (bdf.offset, bdf.scale) == (13.0, -1.0)


def test_urban_equations_reproduce_the_printed_standard_error_of_regression(lib):
    """Table 6's 30 sites against Table 8's SER, in percent (no per-gage estimates are printed).

    Table 8 prints two-decimal exponents and three-figure constants, so 0.25 points.
    """
    sites = list(GAGES["urban"].values())
    assert len(sites) == 30
    assert float(np.mean([s["BDF"] for s in sites])) == pytest.approx(7.97, abs=5e-3)
    assert float(np.mean([s["DRNAREA"] for s in sites])) == pytest.approx(0.778, abs=2e-3)
    assert float(np.mean([s["PRECIP"] for s in sites])) == pytest.approx(37.3, abs=0.05)
    printed = [32.3, 32.8, 33.7, 35.0, 35.9, 36.9]
    for i, aep in enumerate(lib.aeps("Urban")):
        eq = lib.equation("Urban", aep)
        residuals = [
            math.log10(s["flows"][i])
            - evaluate(eq, {c: s[c] for c in ("DRNAREA", "PRECIP", "BDF")}).log_flow
            for s in sites
        ]
        standard_error = math.sqrt(sum(r * r for r in residuals) / (30 - 4))
        assert _sep_percent(standard_error) == pytest.approx(printed[i], abs=0.25)


def test_every_equation_equals_the_nss_snapshot_but_for_its_errors(lib):
    """NSS carries the archive's fit to four decimals, so flows differ by 0.05 percent at most."""
    snap = json.loads((DATA / "nss_snapshot" / "OH.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    harvested = {key for key in by_key if key[0] not in ("GC1233", "GC1234")}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == harvested
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        if e.region_code == "Urban":
            assert e.intercept == pytest.approx(s["intercept"], abs=1e-12)
            assert list(e.coefficients) == s["coefficients"]
            assert s.get("sep_log") is None
        else:
            assert s["intercept"] == round(e.intercept, 4)
            assert s["coefficients"] == [round(c, 4) for c in e.coefficients]
            if not (e.region_code.endswith("full") and e.aep == 0.02):
                assert e.sep_log == pytest.approx(s["sep_log"], abs=2.5e-5)
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.minimum, v.maximum) == (n["code"], n["minimum"], n["maximum"])
