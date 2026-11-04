"""Rhode Island's verified equations (SIR 2012-5109 ver. 1.2, equations 2-8, Tables 8-10, 13, 15, 16).

The report's p. 25 worked example does not come from the published equations. It is pinned
as a known disagreement in :func:`test_p25_worked_example_is_not_this_equation_set`, against
Table 13's row for the same basin, gage 01123000.
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

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "flowfreq" / "data" / "regression"
GAGES = json.loads(
    (ROOT / "tests" / "fixtures" / "regression" / "rhode_island_gages.json").read_text(
        encoding="utf-8"
    )
)
PCT = ["20", "10", "4", "2", "1", "0.5", "0.2"]
ARCADIA = "01117800"  # in Table 13, but not in the regression


def aep(pct: str) -> float:
    return round(float(pct) / 100, 6)


def x_row(da: float, strden: float, stor: float) -> np.ndarray:
    return np.array([1.0, math.log10(da), math.log10(strden), math.log10(stor)])


def prediction_variance(eq, x: np.ndarray) -> float:
    return eq.model_error_variance + float(x @ np.array(eq.covariance) @ x)


def chars(row: dict) -> dict:
    return {"DRNAREA": row["DRNAREA"], "STRDEN": row["STRDEN"], "STORNHD": row["STORNHD"]}


def printed(q: float) -> float:
    """Table 13's rounding: 2 significant figures below 1,000 ft3/s, 3 at or above."""
    keep = 2 if q < 1000 else 3
    unit = 10.0 ** (math.floor(math.log10(q)) - keep + 1)
    return math.floor(q / unit + 0.5) * unit


@pytest.fixture(scope="module")
def lib():
    return load_state("RI")


def test_rhode_island_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["DA-only", "Statewide"]
    assert len(lib.equations) == 14
    assert region_for_nss_code(lib, "GC191") == "Statewide"
    assert "ver. 1.2" in lib.equation("Statewide", 0.01).citation.publication


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "RI.json").read_text(encoding="utf-8"))
    by_aep = {s["aep"]: s for s in snap["equations"]}
    stored = [e for e in lib.equations if e.region_code == "Statewide"]
    assert {e.aep for e in stored} == set(by_aep)
    for e in stored:
        s = by_aep[e.aep]
        assert e.intercept == s["intercept"]
        assert list(e.coefficients) == s["coefficients"]
        assert e.sep_log == pytest.approx(s["sep_log"], abs=1e-6)
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.minimum, v.maximum) == (n["code"], n["minimum"], n["maximum"])


def test_table_13_regression_flows_reproduce_at_all_42_gages(lib):
    assert len(GAGES["table13"]) == 42
    n = 0
    for row in GAGES["table13"]:
        for pct in PCT:
            est = evaluate(
                lib.equation("Statewide", aep(pct)), chars(row), allow_extrapolation=True
            )
            assert printed(est.flow_cfs) == pytest.approx(row["q"][pct]), (row["site"], pct)
            n += 1
    assert n == 294


def test_table_13_variances_reproduce_at_the_regression_gages(lib):
    """Equation 11 at the 41 gages; two values sit on a rounding edge and are named."""
    edge = {("01117000", "0.2"): 0.0145, ("01117420", "1"): 0.0105}
    seen = set()
    for row in GAGES["table13"]:
        if row["site"] == ARCADIA:
            continue
        x = x_row(row["DRNAREA"], row["STRDEN"], row["STORNHD"])
        for pct in PCT:
            v = prediction_variance(lib.equation("Statewide", aep(pct)), x)
            if (row["site"], pct) in edge:
                assert v == pytest.approx(edge[row["site"], pct], abs=1e-5)
                seen.add((row["site"], pct))
            else:
                assert abs(v - row["var"][pct]) <= 0.0005, (row["site"], pct, v)
    assert seen == set(edge)


def test_arcadia_variance_in_table_13_is_not_equation_11(lib):
    """Pinned report slip: it is close to equation 11 only with DA and STRDEN swapped in x."""
    row = next(r for r in GAGES["table13"] if r["site"] == ARCADIA)
    right = x_row(35.14, 1.35, 7.84)
    swapped = np.array([1.0, right[2], right[1], right[3]])
    for pct in PCT:
        eq = lib.equation("Statewide", aep(pct))
        assert prediction_variance(eq, right) < 0.02
        assert row["var"][pct] >= 0.058
        assert prediction_variance(eq, swapped) == pytest.approx(row["var"][pct], abs=0.0012)


def test_average_error_reproduces_table_9_over_the_41_regression_gages(lib):
    gages = [r for r in GAGES["table13"] if r["site"] != ARCADIA]
    assert len(gages) == 41
    for pct in PCT:
        eq = lib.equation("Statewide", aep(pct))
        mean = float(
            np.mean(
                [
                    prediction_variance(eq, x_row(g["DRNAREA"], g["STRDEN"], g["STORNHD"]))
                    for g in gages
                ]
            )
        )
        assert mean == pytest.approx(GAGES["workbook_average_sp_log_variance"][pct], abs=1e-5)
        assert mean == pytest.approx(eq.avp, abs=1e-5)
        percent = 100 * math.sqrt(math.exp(mean * math.log(10) ** 2) - 1)
        # Table 9 prints 24.5 at 2 percent AEP for the workbook's 24.57, a pinned slip.
        expected = 24.57 if pct == "2" else GAGES["table9"][pct]["sep_pct"]
        assert percent == pytest.approx(expected, abs=0.06)
        # Table 9's last column, headed average SEP, is the model error.
        model = 100 * math.sqrt(math.exp(eq.model_error_variance * math.log(10) ** 2) - 1)
        assert model == pytest.approx(GAGES["table9"][pct]["model_pct"], abs=0.06)


def test_covariance_rounds_to_table_10_in_the_stated_order(lib):
    for pct in PCT:
        eq = lib.equation("Statewide", aep(pct))
        table = GAGES["table10"][pct]
        assert round(eq.model_error_variance, 4) == table["mev"]
        assert np.round(np.array(eq.covariance), 5).tolist() == table["U"]


def test_p25_worked_example_is_not_this_equation_set(lib):
    """Pinned disagreement: p. 25 prints 3,590 ft3/s, gamma2 0.00908 and xUx 0.00361.

    Equation 6 and Table 10 give 3,418, 0.0078 and 0.00319. Table 13 prints 3,410 and
    0.011 for gage 01123000, the same basin at 29.99 mi2. Report values for the set win
    over the example: four other places in the report agree with one another.
    """
    eq = lib.equation("Statewide", 0.01)
    site = {"DRNAREA": 30.0, "STRDEN": 3.53, "STORNHD": 6.77}
    est = evaluate(eq, site, interval_level=0.95)
    assert est.flow_cfs == pytest.approx(3418, abs=1)
    assert est.flow_cfs != pytest.approx(3590, rel=0.04)
    x = x_row(30.0, 3.53, 6.77)
    xux = float(x @ np.array(eq.covariance) @ x)
    assert xux == pytest.approx(0.00319, abs=1e-5)
    assert eq.model_error_variance == pytest.approx(0.0078, abs=5e-5)
    # The same basin in Table 13.
    row = next(r for r in GAGES["table13"] if r["site"] == "01123000")
    assert row["q"]["1"] == 3410
    assert row["var"]["1"] == 0.011
    # The interval is equation 9-11 on the stored set, with 41 - 4 degrees of freedom.
    c = 10 ** (2.0262 * math.sqrt(eq.model_error_variance + xux))
    assert est.interval[0] == pytest.approx(est.flow_cfs / c, rel=1e-3)
    assert est.interval[1] == pytest.approx(est.flow_cfs * c, rel=1e-3)


def test_p39_worked_example_at_arcadia_and_hope_valley(lib):
    eq = lib.equation("Statewide", 0.01)
    arcadia = evaluate(eq, {"DRNAREA": 35.14, "STRDEN": 1.35, "STORNHD": 7.84})
    assert printed(arcadia.flow_cfs) == 1590
    hope = next(r for r in GAGES["table13"] if r["site"] == "01118000")
    assert printed(evaluate(eq, chars(hope)).flow_cfs) == 3200
    # Equation 22 with Table 16's b for the 1-percent AEP.
    b = lib.equation("DA-only", 0.01).coefficients[0]
    assert (35.14 / 74.53) ** b * 3160 == pytest.approx(1770, abs=5)


def test_table_15_short_term_gages_follow_their_headings_not_table_14(lib):
    """Flows are printed from rounded inputs, so each is held to its rounding bound.

    Table 15's heading areas differ from Table 14's; the flows follow the headings.
    """
    rows = GAGES["table15"]
    assert len(rows) == 19
    assert sum(r["DA_table15"] != r["DA_table14"] for r in rows) >= 10
    n = 0
    for r in rows:
        da = r["DA_table15"]
        da_half = 0.05 if da >= 10 else 0.005
        for pct, qp in r["q"].items():
            eq = lib.equation("Statewide", aep(pct))
            site = {"DRNAREA": da, "STRDEN": r["STRDEN"], "STORNHD": r["STORNHD"]}
            q = evaluate(eq, site, allow_extrapolation=True).flow_cfs
            b = eq.coefficients
            inputs = (
                abs(b[0]) * da_half / da
                + abs(b[1]) * 0.005 / r["STRDEN"]
                + abs(b[2]) * 0.005 / r["STORNHD"]
            )
            keep = 1 if qp < 10 else (2 if qp < 1000 else 3)
            half_unit = 0.5 * 10.0 ** (math.floor(math.log10(qp)) - keep + 1)
            assert abs(q - qp) <= half_unit + q * inputs, (r["site"], pct, qp, q)
            n += 1
    assert n == 133


def test_drainage_area_only_is_table_16_with_no_interval(lib):
    for pct, (b, constant) in GAGES["table16"].items():
        eq = lib.equation("DA-only", aep(pct))
        assert list(eq.coefficients) == [b]
        assert eq.intercept == constant
        est = evaluate(eq, {"DRNAREA": 50.0})
        assert est.flow_cfs == pytest.approx(10**constant * 50.0**b)
        assert est.interval is None
    assert region_for_nss_code(lib, "GC191") != "DA-only"


def test_a_basin_outside_the_range_is_refused_or_flagged(lib):
    eq = lib.equation("Statewide", 0.01)
    site = {"DRNAREA": 3.0, "STRDEN": 2.0, "STORNHD": 5.0}
    with pytest.raises(OutOfRangeError):
        evaluate(eq, site)
    assert evaluate(eq, site, allow_extrapolation=True).out_of_range
