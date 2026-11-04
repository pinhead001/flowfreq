"""Tennessee's verified equations (SIR 2024-5130; WRIR 84-4110 and 84-4182)."""

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
AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
NSS = {
    "HA1": "GC2015",
    "HA2": "GC2016",
    "HA3": "GC2017",
    "HA4": "GC2018",
    "HA1 DA only": "GC2019",
    "HA2 DA only": "GC2020",
    "HA3 DA only": "GC2021",
    "HA4 DA only": "GC2022",
    "Memphis urban": "GC350",
    "Urban statewide": "GC351",
}
#: WRI 03-4176, superseded by SIR 2024-5130 for the same areas and basin type.
SUPERSEDED = {f"GC{n}" for n in range(341, 350)}
GAGES = json.loads(
    (ROOT / "fixtures" / "regression" / "tennessee_gages.json").read_text(encoding="utf-8")
)
RURAL = GAGES["rural"]

#: Table B gages whose flows follow these drainage areas, not the shapefile's, at every AEP.
#: LC11IMP of 0.87, 1.57, 0.41 and 0.55 fits as well, so which characteristic differs is unknown.
OTHER_DRAINAGE_AREA = {"07022500": 1.72, "07023000": 212.0, "07269000": 19.3, "07279600": 98.3}
#: Two of the 3632 cells miss the rounding bound by under 0.0005, not explained further.
NEAR_ROUNDING = {("03431240", 0.002), ("07029252", 0.1)}
#: Table A/B cell that disagrees with the seven others at its gage: (gage, AEP): (printed, fitted).
MISPRINTED_FLOW = {("03316000", 0.5): (4110.0, 4508.0)}
#: Where NSS follows Table 5's HA3 misprints: (AEP): (stored, NSS) constant or exponent.
NSS_CONSTANT = {("HA3 DA only", 0.2): (540.0, 528.0)}
NSS_EXPONENT = {("HA3 DA only", 0.02): (0.643, 0.64)}


@pytest.fixture(scope="module")
def lib():
    return load_state("TN")


def _third_figure_half_unit(x: float) -> float:
    """Half a unit in the third significant figure, relative to ``x``."""
    return 0.5 * 10 ** (math.floor(math.log10(x)) - 2) / x


def test_tennessee_is_verified(lib):
    assert lib.status == "verified"
    assert sorted(lib.regions) == sorted(NSS)
    assert len(lib.equations) == 76
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
        want = AEPS[:6] if region in ("Memphis urban", "Urban statewide") else AEPS
        assert [e.aep for e in lib.equations if e.region_code == region] == want
    assert all(region_for_nss_code(lib, code) is None for code in SUPERSEDED)


def test_every_equation_equals_the_nss_snapshot_but_for_its_errors(lib):
    snap = json.loads((DATA / "nss_snapshot" / "TN.json").read_text(encoding="utf-8"))
    assert len(snap["equations"]) == 139
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {k for k in by_key if k[0] not in SUPERSEDED} == {
        (NSS[e.region_code], e.aep) for e in lib.equations
    }
    for e in lib.equations:
        key = (e.region_code, e.aep)
        s = by_key[(NSS[e.region_code], e.aep)]
        if key in NSS_CONSTANT:
            assert (10**e.intercept, 10 ** s["intercept"]) == pytest.approx(NSS_CONSTANT[key])
        else:
            assert e.intercept == pytest.approx(s["intercept"], abs=1e-12)
        if key in NSS_EXPONENT:
            assert (e.coefficients[0], s["coefficients"][0]) == NSS_EXPONENT[key]
        else:
            assert list(e.coefficients) == s["coefficients"]
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.minimum, v.maximum) == (n["code"], n["minimum"], n["maximum"])
            assert v.apply(5.0) == pytest.approx(
                math.log10(6.0) if n["transform"] == "log10_plus1" else math.log10(5.0)
            )
        if e.region_code == "Urban statewide":
            assert s.get("sep_log") is None and e.sep_log is not None
        elif e.region_code == "Memphis urban":
            assert s.get("sep_log") is None and e.sep_log is None
        else:
            assert e.sep_log == pytest.approx(s["sep_log"], abs=2.5e-6)


def test_published_regression_flows_reproduce_at_all_454_gages(lib):
    """Data release Tables A and B; the tolerance is the rounding of what is printed."""
    assert len(RURAL) == 454
    checked = 0
    for site, gage in RURAL.items():
        for i, aep in enumerate(AEPS):
            eq = lib.equation(gage["region"], aep)
            inputs = {v.code: gage[v.code] for v in eq.variables}
            if site in OTHER_DRAINAGE_AREA:
                inputs["DRNAREA"] = OTHER_DRAINAGE_AREA[site]
            est = evaluate(eq, inputs, allow_extrapolation=True).flow_cfs
            printed = gage["flows"][i]
            if (site, aep) in MISPRINTED_FLOW:
                assert (printed, est) == pytest.approx(MISPRINTED_FLOW[(site, aep)], abs=1.0)
                continue
            constant = 10**eq.intercept
            tol = _third_figure_half_unit(printed)
            tol += 0.5 / constant if constant >= 1000 else _third_figure_half_unit(constant)
            tol += sum(0.0005 * abs(v.apply(inputs[v.code])) * math.log(10) for v in eq.variables)
            if (site, aep) in NEAR_ROUNDING:
                tol += 0.0005
            assert abs(est / printed - 1) <= tol, (site, aep, est, printed, tol)
            checked += 1
    assert checked == 454 * 8 - len(MISPRINTED_FLOW)


@pytest.mark.parametrize("site", sorted(OTHER_DRAINAGE_AREA))
def test_the_four_pinned_gages_do_not_reproduce_with_the_shapefile_area(lib, site):
    gage = RURAL[site]
    eq = lib.equation(gage["region"], 0.01)
    inputs = {v.code: gage[v.code] for v in eq.variables}
    assert abs(evaluate(eq, inputs).flow_cfs / gage["flows"][5] - 1) > 0.006


def test_the_sep_is_the_mean_per_gage_variance(lib):
    """Table 4's SEP equals sqrt(mean GLS variance) over Tables A and B, to its printed 0.1."""
    for region in ("HA1", "HA2", "HA3", "HA4"):
        rows = [g["variances"] for g in RURAL.values() if g["region"] == region]
        for i, aep in enumerate(AEPS):
            eq = lib.equation(region, aep)
            mean = float(np.mean([r[i] for r in rows]))
            assert eq.n_sites == len(rows)
            assert mean == pytest.approx(eq.sep_log**2, abs=1.5e-4)


def test_avp_is_the_sep_implied_variance_where_it_is_stored(lib):
    """HA2 DA only prints 0.04 in every row against SEPs implying 0.018-0.032, so stores none."""
    for e in lib.equations:
        if e.region_code == "HA2 DA only":
            assert e.avp is None and e.sep_log is not None
        elif e.region_code in ("Memphis urban", "Urban statewide"):
            assert e.avp is None
        else:
            assert round(e.sep_log**2, 2) == e.avp


def test_ha3_has_one_model_stored_twice(lib):
    """Table 5 prints 528 and 0.64 for HA3; Table 4 and the per-gage flows give 540 and 0.643."""
    for aep in AEPS:
        full = lib.equation("HA3", aep)
        da = lib.equation("HA3 DA only", aep)
        assert (full.intercept, list(full.coefficients)) == (da.intercept, list(da.coefficients))
        assert (full.sep_log, full.avp) == (da.sep_log, da.avp)


def test_table_6_ranges_cover_the_gages(lib):
    """Rounded limits of the gages' characteristics; CLIMFAC2YR's 2.05 is 2.056 in the shapefile."""
    for region in ("HA1", "HA2", "HA3", "HA4"):
        eq = lib.equation(region, 0.01)
        for v in eq.variables:
            values = [g[v.code] for g in RURAL.values() if g["region"] == region]
            assert min(values) == pytest.approx(v.minimum, abs=max(0.005, 0.01 * v.minimum))
            assert max(values) == pytest.approx(v.maximum, rel=0.01)


def test_memphis_worked_example(lib):
    """WRIR 84-4110 p. 28: A = 0.91 mi2, channel condition 1.75, the 25-year flood."""
    eq = lib.equation("Memphis urban", 0.04)
    assert evaluate(eq, {"DRNAREA": 0.91, "CHANCOND": 1.75}).flow_cfs == pytest.approx(1950, abs=5)


def test_memphis_refit_reproduces_the_printed_equations(lib):
    """WRIR 84-4110 prints no per-gage regression flows, so refit Table 3 on Table 5 instead."""
    rows = GAGES["memphis"]
    assert len(rows) == 25
    x = np.array(
        [[1.0, math.log10(g["DRNAREA"]), math.log10(g["CHANCOND"])] for g in rows.values()]
    )
    for i, aep in enumerate(AEPS[:6]):
        y = np.log10([g["flows"][i] for g in rows.values()])
        b = np.linalg.lstsq(x, y, rcond=None)[0]
        eq = lib.equation("Memphis urban", aep)
        assert 10 ** b[0] == pytest.approx(10**eq.intercept, rel=0.003)
        assert b[1:] == pytest.approx(list(eq.coefficients), abs=0.005)


def _statewide(lib, aep, a, ia, p):
    eq = lib.equation("Urban statewide", aep)
    inputs = {"CONTDA": a, "IMPERV": ia, "I24H2Y": p}
    return evaluate(eq, inputs, allow_extrapolation=True).flow_cfs


def test_urban_statewide_worked_example(lib):
    """WRIR 84-4182 p. 22: A = 2.0 mi2, 60 percent impervious, 4.0 in, to the printed 3 figures."""
    printed = [1360, 1890, 2180, 2500, 2780, 3060]
    for aep, flow in zip(AEPS[:6], printed):
        assert _statewide(lib, aep, 2.0, 60.0, 4.0) == pytest.approx(flow, abs=5)


def test_urban_statewide_graphical_example(lib):
    """WRIR 84-4182 p. 15: 3.0 mi2, 40 percent, 3.50 in, as read from figures 8-13."""
    printed = [1010, 1520, 1870, 2240, 2650, 3040]
    for aep, flow in zip(AEPS[:6], printed):
        assert _statewide(lib, aep, 3.0, 40.0, 3.5) == pytest.approx(flow, rel=0.005)


#: Table 3's line at this station fits no rainfall with Table 2's area and imperviousness.
STATEWIDE_UNEXPLAINED = "07031777"


def _rainfall_fit(lib, gage):
    """The one 2-year 24-hour rainfall (Table 3 omits it) closest to all six flows, and the misfit."""
    logs, weights = [], []
    for i, aep in enumerate(AEPS[:6]):
        eq = lib.equation("Urban statewide", aep)
        rest = eq.intercept + eq.coefficients[0] * math.log10(gage["CONTDA"])
        rest += eq.coefficients[1] * math.log10(gage["IMPERV"])
        logs.append((math.log10(gage["flows"][i]) - rest) / eq.coefficients[2])
        weights.append(eq.coefficients[2])
    p = 10 ** float(np.average(logs, weights=np.square(weights)))
    return p, [
        _statewide(lib, aep, gage["CONTDA"], gage["IMPERV"], p) / gage["flows"][i] - 1
        for i, aep in enumerate(AEPS[:6])
    ]


def test_urban_statewide_flows_reproduce_with_one_rainfall_per_station(lib):
    rows = GAGES["urban_statewide"]
    assert len(rows) == 22
    for site, gage in rows.items():
        p, misfit = _rainfall_fit(lib, gage)
        if site == STATEWIDE_UNEXPLAINED:
            assert max(abs(m) for m in misfit) > 0.04
            continue
        assert 3.0 <= p <= 4.05, site
        assert max(abs(m) for m in misfit) <= 0.035, (site, misfit)


def test_only_the_urban_statewide_set_has_equivalent_years(lib):
    years = [e.equivalent_years for e in lib.equations if e.region_code == "Urban statewide"]
    assert years == [2, 3, 4, 6, 7, 8]
    assert {e.equivalent_years for e in lib.equations if e.region_code != "Urban statewide"} == {
        None
    }
