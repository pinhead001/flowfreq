"""New York's verified equations: SIR 2006-5112 Tables 1-4, checked against Tables 8 and 10."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "flowfreq" / "data" / "regression"
GAGES = json.loads(
    (ROOT / "tests" / "fixtures" / "regression" / "new_york_gages.json").read_text(encoding="utf-8")
)
FULL = {str(i + 1): f"GC{1071 + i}" for i in range(6)}
DA_ONLY = {f"{i + 1}-DA": f"GC{1065 + i}" for i in range(6)}
NSS = {**FULL, **DA_ONLY}
AEPS = [0.8, 0.667, 0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
#: Table 2, the drainage-area-only 1-percent (100-year) SE and equivalent years.
DA_SE_100 = {"1": (45.4, 4.2), "2": (53.5, 7.4), "3": (33.3, 13.4)}
DA_SE_100 |= {"4": (32.5, 6.7), "5": (48.9, 6.1), "6": (57.6, 2.2)}
TERMS = {
    "DRNAREA": "A",
    "LAGFACTOR": "LAG",
    "STORAGE": "ST",
    "FOREST": "FOR",
    "MAR": "RUNF",
    "PRECIP": "P",
    "MXSNO": "MXSNO",
    "SLOPERATIO": "SR",
    "CSL10_85": "SL",
    "EL1200": "EL12",
}
#: Footnote C gages: a basin in two regions. Only 04232000's weight is printed (example 2).
TWO_REGIONS = {
    "01428500": ("4", 0.328),
    "01434000": ("4", 0.218),
    "01438500": (None, 0.0),
    "01513831": ("4", 0.941),
    "01515000": ("4", 0.834),
    "04232000": ("5", 0.535),
}
#: Footnote D gages (preregulation statistics) where Table 10 does not give Table 8. One
#: changed characteristic reproduces all ten intervals. The report does not print these
#: values, so they are fits: they show what kind of difference this is, nothing more.
PREREGULATION = {
    "01325000": ("LAG", 7.243),
    "01417000": ("LAG", 3.08),
    "01417500": ("LAG", 5.061),
    "01421000": ("LAG", 6.394),
    "01436000": ("LAG", 0.524),
    "01437000": ("LAG", 1.813),
    "01523500": ("SL", 64.6),
}


def chars(gage):
    return {code: gage[term] for code, term in TERMS.items()}


def rounds_to(est: float, printed: int) -> bool:
    """Table 8 prints whole cfs below 100 and three figures above. A three-figure constant in
    Table 1 can itself be off by up to 0.5 percent, so the limit is near that. The largest
    differences found are 0.53 percent and 0.67 cfs."""
    return abs(est - printed) <= 0.75 if printed < 100 else abs(est / printed - 1) <= 0.0055


@pytest.fixture(scope="module")
def lib():
    return load_state("NY")


def test_new_york_is_verified(lib):
    assert lib.status == "verified"
    assert sorted(lib.regions) == sorted(NSS)
    assert len(lib.equations) == 120
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region


def test_urban_regions_stay_unparsed_and_are_not_stored(lib):
    """#212: WRI 84-4350's equations are powers of another equation set's rural peak."""
    assert [u["region_code"] for u in lib.unparsed] == ["GC1077", "GC1078"]
    assert all("#212" in u["reason"] for u in lib.unparsed)
    assert region_for_nss_code(lib, "GC1077") is None
    assert region_for_nss_code(lib, "GC1078") is None


def test_every_equation_equals_the_nss_snapshot_except_the_pinned_se(lib):
    snap = json.loads((DATA / "nss_snapshot" / "NY.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-9)
        ours = {v.code: (v, c) for v, c in zip(e.variables, e.coefficients)}
        theirs = {v["code"]: (v, c) for v, c in zip(s["variables"], s["coefficients"])}
        assert set(ours) == set(theirs)
        for code, (v, c) in ours.items():
            n, nc = theirs[code]
            assert c == nc
            assert v.minimum == n["minimum"]
            assert v.maximum == n["maximum"]
            # NSS writes (LAG + 1) as log10_plus1, the report's form as an offset.
            nss_offset = n.get("offset", 0.0) + (n["transform"] == "log10_plus1")
            assert v.offset == nss_offset
        if e.region_code in FULL:
            assert e.sep_log == pytest.approx(s["sep_log"], abs=6e-4)


def test_nss_serves_the_full_equations_se_on_the_drainage_area_only_regions(lib):
    """Table 2 prints 45.4 percent for region 1's drainage-area-only 100-year equation, and
    30.8 percent for the full one. NSS gives both 30.8. The report wins."""
    snap = json.loads((DATA / "nss_snapshot" / "NY.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    for full, da in zip(FULL, DA_ONLY):
        for aep in AEPS:
            nss_da = by_key[(DA_ONLY[da], aep)]["sep_log"]
            nss_full = by_key[(FULL[full], aep)]["sep_log"]
            assert nss_da == nss_full
            assert nss_da == pytest.approx(lib.equation(full, aep).sep_log, abs=6e-4)
            assert lib.equation(da, aep).sep_log > lib.equation(full, aep).sep_log


@pytest.mark.parametrize("region", list(DA_SE_100))
def test_table_2_at_100_years(lib, region):
    pct, years = DA_SE_100[region]
    eq = lib.equation(f"{region}-DA", 0.01)
    expected = math.sqrt(math.log(1 + (pct / 100) ** 2)) / math.log(10)
    assert eq.sep_log == pytest.approx(expected, abs=1e-5)
    assert eq.equivalent_years == years


def test_example_1_fishkill_creek(lib):
    """Region 2, Q50: 4,210 ft3/s, E = 15.8, weighted with 30 years at 3,600 gives 3,810."""
    lag = 21.00 / ((27.53 + 1) * (12.67 + 1)) ** 0.5
    assert lag == pytest.approx(1.063, abs=5e-4)
    eq = lib.equation("2", 0.02)
    q = evaluate(eq, {"DRNAREA": 57.3, "STORAGE": 4.97, "LAGFACTOR": 1.063, "MAR": 21.29}).flow_cfs
    assert round(q, -1) == 4210
    assert eq.equivalent_years == 15.8
    assert (3600 * 30 + 4210 * eq.equivalent_years) / (30 + eq.equivalent_years) == pytest.approx(
        3810, abs=5
    )


def test_example_2_genesee_river(lib):
    """Regions 5 and 6, Q50: 105,700 and 50,270, weighted 53.5 and 46.5 percent to 79,900."""
    q5 = evaluate(
        lib.equation("5", 0.02), {"DRNAREA": 2467, "CSL10_85": 8.05, "PRECIP": 33.92}
    ).flow_cfs
    q6 = evaluate(
        lib.equation("6", 0.02),
        {"DRNAREA": 2467, "STORAGE": 1.08, "MAR": 14.64, "EL1200": 58.8, "SLOPERATIO": 0.019},
    ).flow_cfs
    assert q5 == pytest.approx(105700, rel=5e-4)
    assert q6 == pytest.approx(50270, rel=5e-4)
    assert 0.535 * q5 + 0.465 * q6 == pytest.approx(79900, rel=5e-4)


def test_example_3_normans_kill(lib):
    """Region 2, Q50 = 9,560 at 131 mi2; Table 3's exponent 0.666 transfers 12,800 to 10,800."""
    q = evaluate(
        lib.equation("2", 0.02),
        {"DRNAREA": 131, "STORAGE": 2.39, "LAGFACTOR": 1.073, "MAR": 17.57},
    ).flow_cfs
    assert round(q, -1) == 9560
    b = lib.equation("2-DA", 0.02).coefficients[0]
    assert b == 0.666
    assert round((131 / 168) ** b * 12800, -2) == 10800


def test_fixture_is_the_report_s_388_gages_by_region():
    assert len(GAGES) == 388
    counts = {r: sum(g["region"] == r for g in GAGES) for r in FULL}
    assert counts == {"1": 108, "2": 63, "3": 39, "4": 49, "5": 80, "6": 49}


def test_table_8_regression_line_at_every_single_region_gage(lib):
    """375 gages, ten intervals each: Table 1 at Table 10's values gives Table 8."""
    checked = 0
    for g in GAGES:
        if g["site"] in TWO_REGIONS or g["site"] in PREREGULATION:
            continue
        for aep, printed in zip(AEPS, g["regression"]):
            est = evaluate(lib.equation(g["region"], aep), chars(g), allow_extrapolation=True)
            assert rounds_to(est.flow_cfs, printed), (g["site"], aep, est.flow_cfs, printed)
            checked += 1
    assert checked == 3750


def test_table_8_weighted_line_is_equation_3_at_single_region_gages(lib):
    """Not stored, checked anyway: (Qg N + Qr E) / (N + E), E from Table 2, within 1 percent."""
    for g in GAGES:
        if "C" in g["footnotes"]:
            continue
        for aep, qg, qr, qw in zip(AEPS, g["at_site"], g["regression"], g["weighted"]):
            e = lib.equation(g["region"], aep).equivalent_years
            est = (qg * g["years"] + qr * e) / (g["years"] + e)
            assert abs(est - qw) <= max(1.0, 0.01 * qw), (g["site"], aep, est, qw)


@pytest.mark.parametrize("site", list(TWO_REGIONS))
def test_table_8_two_region_gages(lib, site):
    """Example 2 prints 04232000's weight; the others are fitted to all ten intervals."""
    other, weight = TWO_REGIONS[site]
    g = next(x for x in GAGES if x["site"] == site)
    assert "C" in g["footnotes"]
    for aep, printed in zip(AEPS, g["regression"]):
        est = evaluate(lib.equation(g["region"], aep), chars(g), allow_extrapolation=True).flow_cfs
        if other is not None:
            second = evaluate(lib.equation(other, aep), chars(g), allow_extrapolation=True)
            est = weight * second.flow_cfs + (1 - weight) * est
        assert est == pytest.approx(printed, rel=0.004)


@pytest.mark.parametrize("site", list(PREREGULATION))
def test_table_8_preregulation_gages_differ_from_table_10_by_one_characteristic(lib, site):
    term, fitted = PREREGULATION[site]
    g = next(x for x in GAGES if x["site"] == site)
    assert "D" in g["footnotes"]
    printed = g["regression"]
    as_listed = [
        evaluate(lib.equation(g["region"], a), chars(g), allow_extrapolation=True).flow_cfs
        for a in AEPS
    ]
    assert max(abs(e / p - 1) for e, p in zip(as_listed, printed)) > 0.009
    changed = {**g, term: fitted}
    for a, p in zip(AEPS, printed):
        est = evaluate(lib.equation(g["region"], a), chars(changed), allow_extrapolation=True)
        assert est.flow_cfs == pytest.approx(p, rel=0.004)
