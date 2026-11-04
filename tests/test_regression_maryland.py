"""Maryland's verified equations (Thomas and Moglen, 2010, Appendix 3 of the Third Edition).

The Third Edition is superseded by the Sixth Edition (July 2023). It is stored because it is
what NSS and StreamStats serve; see ``notes`` in ``MD.json``.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.streamstats import RegressionRegionLocation
from flowfreq.ungaged import OverlappingRegionsError, assemble, region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
EASTERN = "Eastern Coastal Plain"
WESTERN = "Western Coastal Plain"
RURAL = "Piedmont and Blue Ridge (Rural)"
URBAN = "Piedmont (Urban)"
APPALACHIAN = "Appalachian Plateau"
NSS = {
    EASTERN: "GC1478",
    WESTERN: "GC1479",
    RURAL: "GC1480",
    URBAN: "GC1481",
    APPALACHIAN: "GC1482",
}
AEPS = [0.8, 0.667, 0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
#: Appendix 1 and 2 values at the gages the report lists for the Eastern, Western and
#: Appalachian fits (Tables A3.2 and A3.3 and the named outliers).
GAGES = json.loads((ROOT / "fixtures" / "regression" / "maryland_gages.json").read_text("utf-8"))


@pytest.fixture(scope="module")
def lib():
    return load_state("MD")


def test_maryland_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == sorted(NSS)
    assert len(lib.equations) == 50
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
        assert lib.aeps(region) == sorted(AEPS, reverse=True)


def test_notes_say_the_set_is_superseded(lib):
    raw = json.loads((DATA / "MD.json").read_text(encoding="utf-8"))
    assert raw["notes"].startswith("SUPERSEDED")
    assert "Sixth Edition" in raw["source_reports"][0]["superseded_by"]


#: Appendix 5, p. A5-17: discharge, standard error of prediction in percent and in logs,
#: equivalent years and the 90 percent limits, for 10.80 mi2, no limestone, 21 percent forest.
FLAT_RUN = {
    0.8: (498, 42.8, 0.1780, 2.74, 251, 987),
    0.667: (656, 38.2, 0.1603, 3.03, 354, 1220),
    0.5: (833, 36.2, 0.1523, 3.62, 463, 1500),
    0.2: (1520, 31.9, 0.1351, 8.71, 905, 2560),
    0.1: (2160, 31.3, 0.1329, 13.67, 1300, 3610),
    0.04: (3240, 32.8, 0.1387, 19.28, 1900, 5530),
    0.02: (4280, 35.0, 0.1477, 22.05, 2430, 7560),
    0.01: (5560, 38.1, 0.1598, 23.53, 3010, 10300),
    0.005: (7130, 41.6, 0.1736, 24.19, 3650, 13900),
    0.002: (9750, 47.0, 0.1939, 24.20, 4620, 20600),
}
FLAT_INPUTS = {"DRNAREA": 10.80, "LIME": 0.0, "FOREST_MD": 21.0}


@pytest.mark.parametrize("aep", AEPS)
def test_worked_example_flat_run_flows(lib, aep):
    """Printed to three significant figures."""
    est = evaluate(lib.equation(RURAL, aep), FLAT_INPUTS)
    assert est.flow_cfs == pytest.approx(FLAT_RUN[aep][0], rel=0.006)


def test_flat_run_prediction_error_is_the_stored_error_plus_a_constant_leverage(lib):
    """The Tasker output's SEP**2 / SE**2 is 1 + h0 (eq. 2.7), and h0 belongs to the site."""
    h0 = [FLAT_RUN[a][2] ** 2 / lib.equation(RURAL, a).sep_log ** 2 - 1 for a in AEPS]
    assert all(0.0 < h < 0.05 for h in h0)
    assert max(h0) - min(h0) < 0.02


def test_flat_run_intervals_are_slightly_narrow_because_only_the_se_is_printed(lib):
    """Pinned: no (X^T X)^-1 is printed, so the leverage (about 3 percent here) is not stored."""
    for aep in AEPS:
        est = evaluate(lib.equation(RURAL, aep), FLAT_INPUTS, interval_level=0.90)
        _, _, _, _, lower, upper = FLAT_RUN[aep]
        assert lower < est.interval[0] < lower * 1.04
        assert upper * 0.96 < est.interval[1] < upper


#: Table 2.1, p. 2-6: regression flows at 01650500, 21.2 mi2 with 20.1 percent impervious.
COLESVILLE = {0.5: 1550, 0.2: 2920, 0.1: 4260, 0.04: 6550, 0.02: 8860, 0.01: 11700, 0.002: 21600}


def test_worked_example_table_2_1_urban_piedmont(lib):
    for aep, flow in COLESVILLE.items():
        est = evaluate(lib.equation(URBAN, aep), {"DRNAREA": 21.2, "IMPERV": 20.1})
        assert est.flow_cfs == pytest.approx(flow, rel=0.005)
    assert math.log10(COLESVILLE[0.01]) == pytest.approx(4.06819, abs=5e-6)
    assert lib.equation(URBAN, 0.01).equivalent_years == 45


def test_worked_example_ungaged_site_near_the_gage(lib):
    """P. 2-9: 15.1 mi2 with 25 percent impervious, the 100-year flood."""
    est = evaluate(lib.equation(URBAN, 0.01), {"DRNAREA": 15.1, "IMPERV": 25.0})
    assert est.flow_cfs == pytest.approx(9940, rel=0.005)


def test_q500_lime_exponent_misprint_is_read_as_the_series_and_nss_have_it(lib):
    """Pp. A3-4 and A3-35 both print -0165."""
    assert [lib.equation(RURAL, a).coefficients[1] for a in (0.005, 0.002)] == [-0.159, -0.165]


#: Average skew and standard deviation of logs by region, pp. A3-10 and A3-11.
SKEW_AND_S = {
    EASTERN: (0.484, 0.295),
    WESTERN: (0.513, 0.309),
    RURAL: (0.527, 0.296),
    URBAN: (0.585, 0.324),
    APPALACHIAN: (0.489, 0.241),
}


@pytest.mark.parametrize("region", list(SKEW_AND_S))
def test_equivalent_years_follow_from_the_printed_standard_errors(lib, region):
    """Eqs. 2.2 and 2.3. This also proves the SE and years columns are read against the
    right equations, which the summary pages' text layer alone would not."""
    g, s = SKEW_AND_S[region]
    for aep in AEPS:
        eq = lib.equation(region, aep)
        k = stats.pearson3.ppf(1 - aep, g)
        r2 = 1 + g * k + 0.5 * (1 + 0.75 * g**2) * k**2
        nr = (s / eq.sep_log) ** 2 * r2
        # Printed to one decimal below 10 years and to whole years above.
        assert abs(nr - eq.equivalent_years) <= (0.15 if eq.equivalent_years < 10 else 0.55)


#: GC1481 from NSS's UPFS statistic group, which the PFS snapshot does not harvest
#: (Scenarios/Estimate, 2026-11-04). Constant, DRNAREA and (IMPERV+1) exponents, years.
NSS_URBAN = {
    0.8: (17.85, 0.652, 0.635, 3.3),
    0.667: (24.66, 0.648, 0.631, 3.8),
    0.5: (37.01, 0.635, 0.588, 4.5),
    0.2: (94.76, 0.624, 0.499, 13.0),
    0.1: (169.2, 0.622, 0.435, 24.0),
    0.04: (341.0, 0.619, 0.349, 38.0),
    0.02: (562.4, 0.619, 0.284, 44.0),
    0.01: (898.3, 0.619, 0.222, 45.0),
    0.005: (1413.0, 0.621, 0.160, 44.0),
    0.002: (2529.0, 0.623, 0.079, 40.0),
}


def test_every_equation_equals_the_nss_snapshot(lib):
    """No difference from NSS in any constant, exponent, range or equivalent years."""
    snap = json.loads((DATA / "nss_snapshot" / "MD.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    stored = [e for e in lib.equations if e.region_code != URBAN]
    assert {(NSS[e.region_code], e.aep) for e in stored} == set(by_key)
    for e in stored:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-12)
        assert list(e.coefficients) == s["coefficients"]
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.minimum, v.maximum) == (n["code"], n["minimum"], n["maximum"])
            # NSS's (X+1) is log10_plus1; the report's printed form is log10 with offset 1.
            assert v.apply(5.0) == pytest.approx(
                math.log10(5.0 + 1.0 if n["transform"] == "log10_plus1" else 5.0)
            )

    raw = json.loads(
        (ROOT.parent / "data" / "nss_snapshots" / "MD_2026-10-25.json").read_text("utf-8")
    )
    for region in raw["regions"]:
        lib_code = region_for_nss_code(lib, region["region_code"])
        for row in region["equations"]:
            assert lib.equation(lib_code, row["aep"]).equivalent_years == row["equivalent_years"]


def test_urban_set_equals_nss_upfs(lib):
    for aep, (constant, b_area, b_imperv, years) in NSS_URBAN.items():
        e = lib.equation(URBAN, aep)
        assert e.intercept == pytest.approx(math.log10(constant), abs=1e-12)
        assert list(e.coefficients) == [b_area, b_imperv]
        assert e.equivalent_years == years
    limits = [(v.code, v.minimum, v.maximum) for v in lib.equation(URBAN, 0.01).variables]
    assert limits == [("DRNAREA", 0.49, 102.05), ("IMPERV", 10.0, 37.5)]


#: Printed standard error in percent at the 10 return periods, and how close the printed
#: equations' residual standard error about Appendices 1 and 2 comes, in percentage points.
#: The Western gap grows above Q25 and is unexplained: the same data do not return the
#: printed Western exponents above about Q25. Pinned here rather than widened elsewhere.
RESIDUAL_GAP = {EASTERN: 0.45, APPALACHIAN: 0.25, WESTERN: 2.45}


@pytest.mark.parametrize("region", list(RESIDUAL_GAP))
def test_printed_equations_fit_the_appendix_data(lib, region):
    gages = [g for g in GAGES.values() if g["region"] == region]
    assert len(gages) == lib.equation(region, 0.01).n_sites
    for i, aep in enumerate(AEPS):
        eq = lib.equation(region, aep)
        residuals = [
            math.log10(g["flows"][i])
            - evaluate(
                eq, {v.code: g[v.code] for v in eq.variables}, allow_extrapolation=True
            ).log_flow
            for g in gages
        ]
        dof = len(gages) - len(eq.variables) - 1
        rmse = (sum(r**2 for r in residuals) / dof) ** 0.5
        pct = 100 * (math.exp((math.log(10) * rmse) ** 2) - 1) ** 0.5
        printed = 100 * (math.exp((math.log(10) * eq.sep_log) ** 2) - 1) ** 0.5
        assert abs(pct - printed) < RESIDUAL_GAP[region], (region, aep, pct, printed)


@pytest.mark.parametrize("region", list(RESIDUAL_GAP))
def test_printed_ranges_come_from_the_listed_gages(lib, region):
    """To the rounding of the printed ranges, for every variable."""
    gages = [g for g in GAGES.values() if g["region"] == region]
    for v in lib.equation(region, 0.01).variables:
        values = [g[v.code] for g in gages]
        tol = 0.05 if v.units == "percent" or v.code == "DRNAREA" else 5e-5
        assert min(values) == pytest.approx(v.minimum, abs=tol)
        assert max(values) == pytest.approx(v.maximum, abs=tol)


def _located(*pairs):
    return [
        RegressionRegionLocation(
            region_id=0, code=code, name=code, percent_weight=weight, area_sq_mi=1.0
        )
        for code, weight in pairs
    ]


def _assemble(lib, regions, imperv, region=None):
    characteristics = {"DRNAREA": 21.2, "LIME": 0.0, "FOREST_MD": 25.0, "IMPERV": imperv}
    characteristics["SOILCorD"] = 40.0
    return assemble(
        lib,
        regions,
        characteristics,
        aeps=[0.01],
        allow_extrapolation=False,
        weighting_space="log",
        region=region,
    )


def test_imperviousness_chooses_between_the_piedmont_sets(lib):
    """NSS locates GC1480 and GC1481 at 100 percent each (live, 01650500 and 01580000)."""
    both = _located(("GC1480", 100.0), ("GC1481", 100.0))
    urban = _assemble(lib, both, 20.1)
    rural = _assemble(lib, both, 9.9)
    assert urban.quantiles.iloc[0].flow_cfs == pytest.approx(11700, rel=0.005)
    assert re.search("GC1481", str(urban.quantiles.iloc[0].regions))
    assert re.search("GC1480", str(rural.quantiles.iloc[0].regions))
    # The 10 percent line is urban, as p. A3-31 has it.
    assert re.search("GC1481", str(_assemble(lib, both, 10.0).quantiles.iloc[0].regions))


def test_a_piedmont_basin_crossing_a_boundary_must_name_its_regions(lib):
    """The gates need a region to cover the whole basin, so the rule cannot decide."""
    across = _located(("GC1480", 60.0), ("GC1481", 60.0), ("GC1479", 40.0))
    with pytest.raises(OverlappingRegionsError):
        _assemble(lib, across, 5.0)
    done = _assemble(lib, across, 5.0, region=["GC1480", "GC1479"])
    assert np.isfinite(done.quantiles.iloc[0].flow_cfs)
