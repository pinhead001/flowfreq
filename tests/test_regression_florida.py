"""Florida's verified equations: SIR 2011-5034 (rural), WRI 84-4004 (Leon County), WRI 82-42 (Tampa Bay)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
GAGES = json.loads((ROOT / "fixtures" / "regression" / "florida_gages.json").read_text("utf-8"))
NSS = {
    "1": "GC1535",
    "2": "GC1536",
    "3": "GC1537",
    "4": "GC1538",
    "Leon County": "GC1024",
    "Lake Lafayette": "GC1025",
    "Tampa Bay": "GC1030",
}
AEPS8 = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]
LEON_AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.002]
TAMPA_AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01]

#: Appendix 1 cells the equations do not reproduce: (gage, AEP) to printed flow.
RURAL_CELLS = {
    ("02329558", 0.002): 682,  # also printed at 0.5 percent
    ("02330400", 0.02): 5960,
    ("02368329", 0.5): 31,
    ("02368329", 0.02): 185,
    ("02245255", 0.002): 3800,
    ("02246497", 0.1): 824,  # above the 4 percent value, 770
    ("02294068", 0.02): 1290,  # above the 1 percent value, 481
    ("02294650", 0.5): 948,
}
#: Table 16, site 8 at 10 percent: below its own 20 percent value, 448.
TAMPA_CELLS = {("8", 0.1): 374}


def _close(est: float, printed: float) -> bool:
    """Rounding of three-figure inputs and outputs, or of a whole-number small flow."""
    return abs(est - printed) <= max(0.008 * printed, 0.8)


@pytest.fixture(scope="module")
def lib():
    return load_state("FL")


def test_florida_is_verified(lib):
    assert lib.status == "verified"
    assert set(lib.regions) == set(NSS)
    assert len(lib.equations) == 52
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "FL.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-9)
        assert list(e.coefficients) == s["coefficients"]
        for v, n in zip(e.variables, s["variables"]):
            assert v.code == n["code"]
            assert v.minimum == n["minimum"]
            assert v.maximum == n["maximum"]
            assert v.scale == n.get("scale", 1.0)
            # NSS's log10_plus1 is this file's offset of 1.
            plus = 1.0 if n["transform"] == "log10_plus1" else 0.0
            assert v.offset == n.get("offset", 0.0) + plus


def test_nss_differs_only_in_leon_county_sep(lib):
    """NSS has no SEP for WRI 84-4004. Page 24 prints it, and the abstract gives 32 to 47 percent."""
    snap = json.loads((DATA / "nss_snapshot" / "FL.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        if e.region_code in ("Leon County", "Lake Lafayette"):
            assert s.get("sep_log") is None and e.sep_log is not None
        elif e.region_code == "Tampa Bay":
            assert s.get("sep_log") is None and e.sep_log is None
        else:
            assert e.sep_log == pytest.approx(s["sep_log"], abs=6e-4)


@pytest.mark.parametrize("region", ["1", "2", "3", "4"])
def test_table_3_statistics_agree(lib, region):
    """SEP is equation 6 of the AVP, to the percent Table 3 prints."""
    assert lib.aeps(region) == AEPS8
    for aep in AEPS8:
        eq = lib.equation(region, aep)
        assert eq.avp > eq.model_error_variance
        assert eq.avp**0.5 == pytest.approx(eq.sep_log, abs=0.003)


def test_rural_worked_example_new_river():
    """SIR 2011-5034 p. 12: an ungaged site above 02330400, 1-percent AEP."""
    lib = load_state("FL")
    eq = lib.equation("1", 0.01)
    ungaged = evaluate(eq, {"DRNAREA": 115.0, "LC92STOR": 30.1}).flow_cfs
    assert ungaged == pytest.approx(6240, rel=0.001)
    # Step 5's 6,640 is appendix 1's cell. The equation gives 6,622 at the printed inputs.
    gaged = evaluate(eq, {"DRNAREA": 151.0, "LC92STOR": 43.03}).flow_cfs
    assert gaged == pytest.approx(6622, abs=1)
    assert GAGES["rural"]["02330400"]["flows"][5] == 6640
    # Equation 41 with the printed values. Step 1's 7,010 is 7,090 in appendix 1.
    ratio = 2 * 36 / 151
    assert (ratio + (1 - ratio) * 7010 / 6640) * 6240 == pytest.approx(6420, abs=5)


def test_appendix_1_regression_line_reproduces_at_florida_gages(lib):
    assert len(GAGES["rural"]) == 275
    cells = 0
    for station, gage in GAGES["rural"].items():
        for aep, printed in zip(AEPS8, gage["flows"]):
            est = evaluate(
                lib.equation(gage["region"], aep),
                {"DRNAREA": gage["DRNAREA"], "LC92STOR": gage["LC92STOR"]},
                allow_extrapolation=True,
            ).flow_cfs
            if (station, aep) in RURAL_CELLS:
                assert printed == RURAL_CELLS[(station, aep)]
                assert not _close(est, printed), (station, aep, est)
            else:
                assert _close(est, printed), (station, aep, est, printed)
            cells += 1
    assert cells == 2200


def test_leon_worked_example_st_augustine_branch(lib):
    """WRI 84-4004 pp. 30-31: map 13, 2.06 mi2, 54 percent impervious."""
    printed = [1330, 1840, 2180, 2620, 2950, 3280, 4170]
    assert lib.aeps("Leon County") == LEON_AEPS
    for aep, q in zip(LEON_AEPS, printed):
        est = evaluate(lib.equation("Leon County", aep), {"DRNAREA": 2.06, "IMPERV": 54.0})
        assert est.flow_cfs == pytest.approx(q, rel=0.005)
        assert est.interval is not None


def test_leon_supplementary_data_reproduce(lib):
    assert len(GAGES["leon"]) == 13
    assert sum(g["region"] == "Lake Lafayette" for g in GAGES["leon"].values()) == 5
    for name, gage in GAGES["leon"].items():
        for aep, printed in zip(LEON_AEPS, gage["flows"]):
            est = evaluate(
                lib.equation(gage["region"], aep),
                {"DRNAREA": gage["DRNAREA"], "IMPERV": gage["IMPERV"]},
            ).flow_cfs
            assert _close(est, printed), (name, aep, est, printed)


@pytest.mark.parametrize(
    "aep,inputs,printed",
    [
        (0.04, {"DRNAREA": 2.80, "BDF": 10.0, "CSL10_85": 10.9}, 2050),  # example 1
        (0.02, {"DRNAREA": 2.02, "BDF": 0.0, "CSL10_85": 13.5}, 914),  # example 2
        (0.02, {"DRNAREA": 2.02, "BDF": 9.0, "CSL10_85": 13.5}, 1750),  # example 3
        (0.02, {"DRNAREA": 2.02, "BDF": 6.0, "CSL10_85": 13.5}, 1280),  # example 3, swales
        (0.01, {"DRNAREA": 2.65, "BDF": 2.0, "CSL10_85": 2.87}, 573),  # example 4
    ],
)
def test_tampa_worked_examples(lib, aep, inputs, printed):
    """WRI 82-42 pp. 47-49; examples 2 and 4 lie outside the calibrated range."""
    est = evaluate(lib.equation("Tampa Bay", aep), inputs, allow_extrapolation=True)
    assert est.flow_cfs == pytest.approx(printed, rel=0.005)


def test_tampa_table_16_reproduces(lib):
    assert lib.aeps("Tampa Bay") == TAMPA_AEPS
    assert len(GAGES["tampa"]) == 9
    for name, gage in GAGES["tampa"].items():
        inputs = {k: gage[k] for k in ("DRNAREA", "BDF", "CSL10_85", "DETEN")}
        for aep, printed in zip(TAMPA_AEPS, gage["flows"]):
            eq = lib.equation("Tampa Bay", aep)
            est = evaluate(eq, {v.code: inputs[v.code] for v in eq.variables}).flow_cfs
            if (name, aep) in TAMPA_CELLS:
                assert printed == TAMPA_CELLS[(name, aep)]
                assert est == pytest.approx(574, abs=1)
                assert printed < gage["flows"][1]
            else:
                assert _close(est, printed), (name, aep, est, printed)


def test_tampa_prints_no_standard_error_of_prediction(lib):
    for aep in TAMPA_AEPS:
        est = evaluate(
            lib.equation("Tampa Bay", aep),
            {"DRNAREA": 1.0, "BDF": 8.0, "CSL10_85": 10.0, "DETEN": 1.0},
        )
        assert est.interval is None


def test_tampa_forms_switch_at_four_percent(lib):
    """Equation 3 uses BDF and detention; equation 4 uses 13-BDF and no detention."""
    ten = lib.equation("Tampa Bay", 0.1)
    four = lib.equation("Tampa Bay", 0.04)
    assert [v.code for v in ten.variables] == ["DRNAREA", "BDF", "CSL10_85", "DETEN"]
    assert ten.variables[3].offset == 0.01
    assert [v.code for v in four.variables] == ["DRNAREA", "BDF", "CSL10_85"]
    assert (four.variables[1].scale, four.variables[1].offset) == (-1.0, 13.0)
