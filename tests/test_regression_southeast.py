"""Georgia, North and South Carolina verified equations.

SIR 2023-5006 (Tables 1-5, with Kolb and others' data release and Weaver and others' model archive)
and SIR 2014-5030 ver. 1.1 (Tables 5, 7-11, and its spreadsheet tool).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.equations import OutOfRangeError
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
NSS = {
    "GA": {
        "HR1-5": "GC1934",
        "HR1": "GC1572",
        "HR3": "GC1573",
        "HR4": "GC1574",
        "HR5": "GC1575",
    },
    "NC": {"HR1-5": "GC1936", "HR1": "GC1580", "HR3": "GC1581", "HR4": "GC1582"},
    "SC": {"HR1-5": "GC1935", "HR1": "GC1587", "HR3": "GC1588", "HR4": "GC1589"},
}
STATES = list(NSS)
#: Kolb table 3 follows the WREG input (SiteInfo.txt) at these, not Kolb table 1 (DA, pct 1-5).
SITEINFO = {
    "02103500": (456.0, [0.4, 0.0, 86.2, 13.4, 0.0]),
    "02114450": (42.9, [100.0, 0.0, 0.0, 0.0, 0.0]),
    "02134500": (1229.0, [0.0, 0.0, 0.0, 100.0, 0.0]),
    "02143000": (83.6, [38.5, 61.5, 0.0, 0.0, 0.0]),
}
#: NSS's minimum LC06IMP in HR4 against Table 11's.
NSS_HR4_IMP_MIN, REPORT_HR4_IMP_MIN = 0.2, 0.02
GAGES = json.loads((ROOT / "fixtures" / "regression" / "southeast_gages.json").read_text("utf-8"))


@pytest.fixture(scope="module", params=STATES)
def state(request):
    return request.param


@pytest.fixture(scope="module")
def ga():
    return load_state("GA")


def _inputs_2023(da, pct):
    return {
        "DRNAREA": da,
        "PCTREG1": pct[0],
        "PCTREG2": pct[1],
        "PCTREG3": pct[2],
        "PCTREG5": pct[4],
        "LC19IMP": 1.0,
    }


def _x(eq, values):
    return np.array([1.0] + [v.value(values) for v in eq.variables])


def test_status_regions_and_nss_codes(state):
    lib = load_state(state)
    assert lib.status == "verified"
    assert sorted(lib.regions) == sorted(NSS[state])
    assert len(lib.equations) == 8 * len(NSS[state])
    for region, code in NSS[state].items():
        assert region_for_nss_code(lib, code) == region
        assert lib.aeps(region) == sorted(AEPS, reverse=True)


def test_every_equation_against_the_nss_snapshot(state):
    lib = load_state(state)
    snap = json.loads((DATA / "nss_snapshot" / f"{state}.json").read_text("utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[state][e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[state][e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-9)
        assert e.sep_log == pytest.approx(s["sep_log"], abs=6e-4)
        stored = {(v.code, v.times): (c, v) for v, c in zip(e.variables, e.coefficients)}
        served = {
            (n["code"], n.get("times")): (c, n) for n, c in zip(s["variables"], s["coefficients"])
        }
        assert set(stored) == set(served)
        for key, (coef, v) in stored.items():
            nss_coef, n = served[key]
            assert coef == nss_coef
            assert v.maximum == n["maximum"]
            if key[0] == "LC06IMP" and e.region_code == "HR4":
                assert (v.minimum, n["minimum"]) == (REPORT_HR4_IMP_MIN, NSS_HR4_IMP_MIN)
            else:
                assert v.minimum == n["minimum"]
        if e.region_code == "HR1-5":
            # The report's limit 3 holds at every AEP; NSS's equations carry LC19IMP only
            # in South Carolina's 50-percent form, and its no-op PCTREG4 only at 50 percent.
            assert {v.code for v in e.applicability} == {"PCTREG2", "LC19IMP"}
            served_limits = {n["code"] for n in s["applicability"]}
            if e.aep != 0.5:
                assert served_limits == {"PCTREG2"}
            elif state == "SC":
                assert served_limits == {"PCTREG2", "LC19IMP", "PCTREG4"}
            else:
                assert served_limits == {"PCTREG2", "PCTREG4"}


def test_the_three_states_share_the_equations():
    libs = {st: load_state(st) for st in STATES}
    for e in libs["GA"].equations:
        if e.region_code == "HR5":
            continue
        for st in ("NC", "SC"):
            o = libs[st].equation(e.region_code, e.aep)
            assert (o.intercept, o.coefficients, o.covariance) == (
                e.intercept,
                e.coefficients,
                e.covariance,
            )
            assert (o.sep_log, o.avp, o.model_error_variance) == (
                e.sep_log,
                e.avp,
                e.model_error_variance,
            )


def test_2023_worked_examples_yadkin_river(ga):
    """The report's 1-percent examples: gage 02116500, and an ungaged site above it."""
    eq = ga.equation("HR1-5", 0.01)
    gage = _inputs_2023(2278.0, [74.7, 25.3, 0.0, 0.0, 0.0])
    est = evaluate(eq, gage)
    assert est.flow_cfs == pytest.approx(93525, abs=0.5)
    x = _x(eq, gage)
    assert eq.model_error_variance + float(x @ np.array(eq.covariance) @ x) == pytest.approx(
        0.0293, abs=5e-5
    )
    ungaged = _inputs_2023(1690.0, [66.0, 34.0, 0.0, 0.0, 0.0])
    assert evaluate(eq, ungaged).flow_cfs == pytest.approx(81931, abs=0.5)


def test_2023_interval_uses_the_interaction_row(ga):
    eq = ga.equation("HR1-5", 0.01)
    gage = _inputs_2023(2278.0, [74.7, 25.3, 0.0, 0.0, 0.0])
    x = _x(eq, gage)
    assert x[-1] == pytest.approx(25.3 * math.log10(2278.0))
    est = evaluate(eq, gage, interval_level=0.95)
    assert est.interval_method.startswith("site-specific")
    assert est.prediction_variance == pytest.approx(
        eq.model_error_variance + float(x @ np.array(eq.covariance) @ x)
    )


@pytest.mark.parametrize(
    "aep, hr1, hr2, b2, hr3, hr4, hr5, b",
    [
        # Table 3, read against Table 2.
        (0.5, 149, 66.1, 0.870, 41.5, 66.1, 102, 0.646),
        (0.2, 267, 132, 0.830, 75.2, 132, 223, 0.631),
        (0.1, 361, 191, 0.810, 104, 191, 340, 0.623),
        (0.04, 491, 275, 0.790, 143, 275, 520, 0.615),
        (0.02, 607, 355, 0.778, 178, 355, 697, 0.610),
        (0.01, 721, 437, 0.766, 213, 437, 889, 0.605),
        (0.005, 839, 525, 0.757, 251, 525, 1107, 0.601),
        (0.002, 995, 646, 0.747, 300, 646, 1419, 0.597),
    ],
)
def test_2023_table_3_is_table_2_at_100_percent(ga, aep, hr1, hr2, b2, hr3, hr4, hr5, b):
    eq = ga.equation("HR1-5", aep)
    pct = [
        [100, 0, 0, 0, 0],
        [0, 100, 0, 0, 0],
        [0, 0, 100, 0, 0],
        [0, 0, 0, 100, 0],
        [0, 0, 0, 0, 100],
    ]
    expected = [(hr1, b), (hr2, b2), (hr3, b), (hr4, b), (hr5, b)]
    for p, (c, slope) in zip(pct, expected):
        at_1 = evaluate(eq, _inputs_2023(1.0, p)).flow_cfs
        at_10 = evaluate(eq, _inputs_2023(10.0, p)).flow_cfs
        assert at_1 == pytest.approx(c, rel=0.0052)
        assert math.log10(at_10 / at_1) == pytest.approx(slope, abs=1e-9)


def test_2023_every_gage_flow(ga):
    """Kolb table 3's regression flows, 8 AEPs at the 807 gages in the three states."""
    worst = 0.0
    rural = GAGES["rural_2023"]
    assert len(rural) == 807
    for site, g in rural.items():
        da, pct = SITEINFO.get(site, (g["da"], g["pct"]))
        for aep, q in zip(AEPS, g["q"]):
            est = evaluate(
                ga.equation("HR1-5", aep), _inputs_2023(da, pct), allow_extrapolation=True
            )
            worst = max(worst, abs(est.flow_cfs / q - 1))
    # Rounding of the 3-figure flow and of table 1's area and percentages.
    assert worst < 0.0065


@pytest.mark.parametrize("site", sorted(SITEINFO))
def test_2023_table_1_differs_from_the_wreg_input_at_four_gages(ga, site):
    """Pinned: with table 1's characteristics these miss table 3; with the WREG input they hit it."""
    g = GAGES["rural_2023"][site]
    eq = ga.equation("HR1-5", 0.01)
    da, pct = SITEINFO[site]
    assert evaluate(eq, _inputs_2023(da, pct)).flow_cfs == pytest.approx(g["q"][5], rel=0.0065)
    assert abs(evaluate(eq, _inputs_2023(g["da"], g["pct"])).flow_cfs / g["q"][5] - 1) > 0.01


def test_2023_covariance_reproduces_every_regressed_gage_variance(ga):
    """MEV + xUx against Kolb table 3 at the 670 gages regressed, and its mean against AVP."""
    regressed = {s: g for s, g in GAGES["rural_2023"].items() if "var" in g}
    assert len(regressed) == 670
    for i, aep in enumerate(AEPS):
        eq = ga.equation("HR1-5", aep)
        u = np.array(eq.covariance)
        assert u == pytest.approx(u.T)
        diffs, printed = [], []
        for site, g in regressed.items():
            da, pct = SITEINFO.get(site, (g["da"], g["pct"]))
            x = _x(eq, _inputs_2023(da, pct))
            diffs.append(eq.model_error_variance + float(x @ u @ x) - g["var"][i])
            printed.append(g["var"][i])
        assert max(abs(d) for d in diffs) < 2e-4, aep
        assert float(np.sqrt(np.mean(np.square(diffs)))) < 6e-5, aep
        assert float(np.mean(printed)) == pytest.approx(eq.avp, abs=1e-4), aep


def test_2023_four_percent_statistics_are_misprinted(ga):
    """Table 4's AVP 0.0251 and Kolb table 9's MEV 0.0261 are wrong; WREG has 0.026 and 0.0251."""
    eq = ga.equation("HR1-5", 0.04)
    assert (eq.avp, eq.model_error_variance) == (0.0260, 0.0251)
    # Table 4's own SEP, by eq. 10, needs the 0.0260.
    sep = lambda avp: 100 * math.sqrt(10 ** (2.3026 * avp) - 1)  # noqa: E731
    assert sep(0.0260) == pytest.approx(38.4, abs=0.1)
    assert abs(sep(0.0251) - 38.4) > 0.5
    # With Kolb's MEV, MEV + xUx overshoots table 3's variance at every regressed gage.
    u = np.array(eq.covariance)
    for site, g in GAGES["rural_2023"].items():
        if "var" in g:
            da, pct = SITEINFO.get(site, (g["da"], g["pct"]))
            x = _x(eq, _inputs_2023(da, pct))
            assert 0.0261 + float(x @ u @ x) - g["var"][3] > 5e-4, site


def test_2023_limit_on_impervious_area(ga):
    for aep in AEPS:
        inputs = {**_inputs_2023(10.0, [100, 0, 0, 0, 0]), "LC19IMP": 10.5}
        with pytest.raises(OutOfRangeError, match="LC19IMP"):
            evaluate(ga.equation("HR1-5", aep), inputs)


def test_2023_needs_the_percentage_in_region_2_not_region_4(ga):
    inputs = _inputs_2023(10.0, [100, 0, 0, 0, 0])
    del inputs["PCTREG2"]
    with pytest.raises(KeyError, match="PCTREG2"):
        evaluate(ga.equation("HR1-5", 0.01), inputs)
    assert "PCTREG4" not in {c for v in ga.equation("HR1-5", 0.01).variables for c in v.input_codes}


def test_2014_worked_example_unnamed_creek(ga):
    """Region 1, 0.50 mi2, 30 percent impervious, 1 percent: Q reproduces, the printed xUx does not."""
    eq = ga.equation("HR1", 0.01)
    values = {"DRNAREA": 0.5, "LC06IMP": 30.0}
    assert evaluate(eq, values).flow_cfs == pytest.approx(509, abs=0.5)
    x = _x(eq, values)
    assert x == pytest.approx([1.0, math.log10(0.5), 30.0])
    xux = float(x @ np.array(eq.covariance) @ x)
    # Pinned: the report prints 0.001447 (S 0.1358, 272 to 952); Table 10 and the tool give 0.00196.
    assert xux == pytest.approx(0.00196, abs=1e-5)
    assert abs(xux - 0.001447) > 5e-4
    assert eq.model_error_variance == 0.017


def test_2014_hr1_at_two_tenths_percent_has_no_impervious_term(ga):
    eq = ga.equation("HR1", 0.002)
    assert [v.code for v in eq.variables] == ["DRNAREA"]
    assert np.array(eq.covariance).shape == (2, 2)
    assert evaluate(eq, {"DRNAREA": 0.5}).flow_cfs == pytest.approx(1045 * 0.5**0.8160)


def test_2014_every_gage_flow(ga):
    """Table 5's regression flows at the gages regressed under the stored columns."""
    small = GAGES["small_2014"]
    assert {r: sum(g["hr"] == r for g in small.values()) for r in ("HR1", "HR3", "HR4")} == {
        "HR1": 41,
        "HR3": 10,
        "HR4": 19,
    }
    worst = 0.0
    for g in small.values():
        for aep, q in zip(AEPS, g["q"]):
            eq = ga.equation(g["hr"], aep)
            values = {
                "DRNAREA": g["da"],
                "LC06IMP": g["imp"],
                "LC06DEV": g["dev"],
                "I24H50Y": g["p"],
            }
            est = evaluate(eq, values, allow_extrapolation=True)
            worst = max(worst, abs(est.flow_cfs / q - 1))
    # Table 7's constants have 3 figures and its exponents 4 decimals, on percentages up to 98.
    assert worst < 0.012


@pytest.mark.parametrize("region, n", [("HR1", 62), ("HR3", 34), ("HR4", 157)])
def test_2014_covariance_returns_table_9_avp(ga, region, n):
    """MEV + mean(xUx) over Table 8's gages, against Table 9's 3-decimal AVP."""
    rows = GAGES["table8"][region]
    assert len(rows) == n
    for aep in AEPS:
        eq = ga.equation(region, aep)
        assert eq.n_sites == n
        u = np.array(eq.covariance)
        xs = []
        for da, pct, p in rows:
            values = {"DRNAREA": da, "LC06IMP": pct, "LC06DEV": pct, "I24H50Y": p}
            xs.append(_x(eq, values))
        mean_xux = float(np.mean([x @ u @ x for x in xs]))
        assert eq.model_error_variance + mean_xux == pytest.approx(eq.avp, abs=1.1e-3), aep


def test_2014_hr5_against_the_spreadsheet_tool(ga):
    """No gages or worked example for HR5: the tool's 0.5 mi2 column, its only check."""
    flows = [114, 177, 230, 309, 373, 441, 516, 619]
    sp = [0.17602, 0.14607, 0.14521, None, 0.16168, 0.17439, 0.18667, 0.20757]
    for aep, q, s in zip(AEPS, flows, sp):
        eq = ga.equation("HR5", aep)
        est = evaluate(eq, {"DRNAREA": 0.5})
        assert est.flow_cfs == pytest.approx(q, rel=0.005)
        x = _x(eq, {"DRNAREA": 0.5})
        got = math.sqrt(eq.model_error_variance + float(x @ np.array(eq.covariance) @ x))
        if s is not None:
            assert got == pytest.approx(s, abs=2e-4)
    # Pinned: the tool's 4-percent row reuses the 2-percent xUx (S 0.15496); Table 10 gives 0.1527.
    eq = ga.equation("HR5", 0.04)
    x = _x(eq, {"DRNAREA": 0.5})
    got = math.sqrt(eq.model_error_variance + float(x @ np.array(eq.covariance) @ x))
    assert got == pytest.approx(0.1527, abs=2e-4)
    assert eq.n_sites is None


def test_2014_table_7_heading_against_table_11(ga):
    """Table 7's heading gives HR3 from 0.22 mi2, the urban minimum; Table 11 and Table 8 have 0.14."""
    assert ga.equation("HR3", 0.01).variables[0].minimum == 0.14
    assert min(r[0] for r in GAGES["table8"]["HR3"]) == 0.14


@pytest.mark.parametrize("region", ["HR1", "HR3", "HR4", "HR5"])
def test_2014_upper_area_is_the_small_rural_definition(ga, region):
    eq = ga.equation(region, 0.01)
    values = {"DRNAREA": 1.5, "LC06IMP": 1.0, "LC06DEV": 3.0, "I24H50Y": 8.0}
    assert eq.variables[0].maximum == 1.0
    with pytest.raises(OutOfRangeError, match="DRNAREA"):
        evaluate(eq, values)


def test_the_2023_set_comes_first_where_both_apply(state):
    assert load_state(state).region_preference.tiers[0] == ["HR1-5"]
