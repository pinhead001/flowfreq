"""New Jersey's verified equations (SIR 2009-5167 v2.0, equations 7-13, Tables 3-6, appendix 2)."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict

import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {"GP": "GC1255", "IN": "GC1256", "OC": "GC1257", "UP": "GC1258", "VR": "GC1259"}
#: Table 3's characteristics as printed, and appendix 2's regression flows, at the 254 gages.
GAGES = json.loads((ROOT / "fixtures" / "regression" / "new_jersey_gages.json").read_text("utf-8"))
#: Where appendix 2's regression flows imply another drainage area than Table 3 prints, as
#: (Table 3's value, a value inside the range the seven printed flows allow). The other three
#: characteristics stand. NWIS agrees with Table 3 at the first four (3.74, 11.8, 44.9, 15.8),
#: so appendix 2 used a different area, and the fault is not in the equations.
TABLE_3_AREA_NOT_USED = {
    "01401200": (3.82, 6.70),
    "01396660": (11.7, 9.45),
    "01405030": (44.7, 42.15),
    "01409280": (15.8, 16.25),
    "01386500": (1.01, 0.997),
}
CODES = ("DRNAREA", "STORAGE", "CSL10_85", "POPDENS")


@pytest.fixture(scope="module")
def lib():
    return load_state("NJ")


def _half_unit(printed: str) -> float:
    """Half a unit in the last place of a printed value, at three figures for a large integer."""
    text = printed.replace(",", "")
    if "." in text:
        return 0.5 * 10.0 ** -len(text.split(".")[1])
    return 0.5 * 10.0 ** max(0, len(text) - 3)


def _value(printed: str) -> float:
    return float(printed.replace(",", ""))


def _outside_rounding(lib, gage: Dict, area: float | None = None) -> list:
    """The AEPs whose printed flow is farther off than the printed digits can explain."""
    outside = []
    for aep, printed_flow in zip(lib.aeps(gage["region"]), gage["flows"]):
        eq = lib.equation(gage["region"], aep)
        inputs = {code: _value(gage[code]) for code in CODES}
        half = {code: _half_unit(gage[code]) for code in CODES}
        if area is not None:
            inputs["DRNAREA"], half["DRNAREA"] = area, 0.0
        est = evaluate(eq, inputs, allow_extrapolation=True)
        allowed = 0.0
        for var, coefficient in zip(eq.variables, eq.coefficients):
            base = inputs[var.code] + var.offset
            allowed += abs(coefficient) * math.log10((base + half[var.code]) / base)
        flow_half = _half_unit(str(printed_flow))
        allowed += math.log10(printed_flow / (printed_flow - flow_half))
        if abs(est.log_flow - math.log10(printed_flow)) > allowed:
            outside.append(aep)
    return outside


def test_new_jersey_is_verified(lib):
    assert lib.status == "verified"
    assert sorted(lib.regions) == sorted(NSS)
    assert len(lib.equations) == 35
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region


def test_regions_share_slopes_and_differ_in_the_constant(lib):
    """Equations 7-13 and Table 5, whose columns run OC, UP, GP, VR, IN."""
    for aep in (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002):
        slopes = {lib.equation(r, aep).coefficients for r in NSS}
        assert len(slopes) == 1
    constants = {r: 10 ** lib.equation(r, 0.01).intercept for r in NSS}
    assert constants == pytest.approx(
        {"OC": 65.1, "UP": 256.0, "GP": 130.0, "VR": 106.0, "IN": 146.0}, abs=1e-9
    )


def test_one_is_added_to_storage_and_population_density(lib):
    eq = lib.equation("GP", 0.01)
    assert [(v.code, v.transform, v.offset) for v in eq.variables] == [
        ("DRNAREA", "log10", 0.0),
        ("STORAGE", "log10", 1.0),
        ("CSL10_85", "log10", 0.0),
        ("POPDENS", "log10", 1.0),
    ]


@pytest.mark.parametrize(
    "region, area, storage, slope, density, flow, precision",
    [
        # Atlantic County, ungaged, region OC.
        ("OC", 15.5, 5.5, 3.2, 75.0, 639.0, 0.5),
        # 01394500, Rahway River near Springfield.
        ("GP", 24.9, 2.73, 34.1, 3980.0, 4770.0, 5.0),
        # Ungaged site downstream of 01397500, Walnut Brook near Flemington.
        ("UP", 3.0, 15.2, 92.0, 250.0, 1760.0, 5.0),
        # Ungaged site on Beden Brook, printed to two figures.
        ("UP", 10.0, 12.9, 41.2, 430.0, 3600.0, 50.0),
    ],
)
def test_worked_examples(lib, region, area, storage, slope, density, flow, precision):
    est = evaluate(
        lib.equation(region, 0.01),
        {"DRNAREA": area, "STORAGE": storage, "CSL10_85": slope, "POPDENS": density},
    )
    assert est.flow_cfs == pytest.approx(flow, abs=precision)


def test_weighted_example_uses_table_6_years(lib):
    """Equation 14 at Rahway: 67 years and table 6's 5, for the printed 6,330 ft3/s."""
    eq = lib.equation("GP", 0.01)
    assert eq.equivalent_years == 5.0
    weighted = (67 * math.log10(6470) + eq.equivalent_years * math.log10(4770)) / (67 + 5)
    assert 10**weighted == pytest.approx(6330, abs=5)


def test_equivalent_years_are_table_6(lib):
    aeps = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002)
    assert [lib.equation("VR", a).equivalent_years for a in aeps] == [1, 2, 3, 4, 5, 5, 6]


def test_regression_flows_reproduce_at_every_gage_within_printed_rounding(lib):
    """All 1,778 flows of appendix 2's R column, but for the five gages pinned above."""
    assert len(GAGES) == 254
    outside = {}
    for site, gage in GAGES.items():
        missed = _outside_rounding(lib, gage)
        if missed:
            outside[site] = missed
    assert set(outside) == set(TABLE_3_AREA_NOT_USED)
    assert sum(map(len, outside.values())) == 34


@pytest.mark.parametrize("site", sorted(TABLE_3_AREA_NOT_USED))
def test_pinned_gages_reproduce_on_the_area_the_flows_imply(lib, site):
    printed, implied = TABLE_3_AREA_NOT_USED[site]
    assert _value(GAGES[site]["DRNAREA"]) == printed
    assert _outside_rounding(lib, GAGES[site], area=implied) == []


def test_table_4_ranges_are_table_3s_extremes(lib):
    """Table 3 prints population density to three figures, and Table 4 in full."""
    for region in NSS:
        eq = lib.equation(region, 0.01)
        rows = [g for g in GAGES.values() if g["region"] == region]
        for var in eq.variables:
            values = [_value(g[var.code]) for g in rows]
            if var.code == "POPDENS":
                assert min(values) == float(f"{var.minimum:.3g}")
                assert max(values) == float(f"{var.maximum:.3g}")
            else:
                assert (min(values), max(values)) == (var.minimum, var.maximum)


def test_every_equation_equals_the_nss_snapshot(lib):
    """NSS writes the added 1 as ``log10_plus1``, which evaluates the same."""
    snap = json.loads((DATA / "nss_snapshot" / "NJ.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-9)
        assert list(e.coefficients) == s["coefficients"]
        assert e.sep_log == pytest.approx(s["sep_log"], abs=1e-6)
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.minimum, v.maximum) == (n["code"], n["minimum"], n["maximum"])
            assert v.offset == (1.0 if n["transform"] == "log10_plus1" else 0.0)


def test_no_covariance_is_printed_so_intervals_use_the_sep(lib):
    eq = lib.equation("UP", 0.01)
    assert eq.covariance is None and eq.avp is None and eq.model_error_variance is None
    est = evaluate(eq, {"DRNAREA": 10.0, "STORAGE": 12.9, "CSL10_85": 41.2, "POPDENS": 430.0})
    assert est.prediction_variance == pytest.approx(eq.sep_log**2)
    assert est.interval is not None


def test_out_of_range_population_density_is_flagged(lib):
    eq = lib.equation("VR", 0.01)
    inputs = {"DRNAREA": 10.0, "STORAGE": 5.0, "CSL10_85": 40.0, "POPDENS": 2000.0}
    assert evaluate(eq, inputs, allow_extrapolation=True).out_of_range
    with pytest.raises(ValueError):
        evaluate(eq, inputs)
