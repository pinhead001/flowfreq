"""Pennsylvania's verified equations (SIR 2019-5094 version 1.1, Tables 3 and 5, appendixes 1-2)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {"1": "GC907", "2": "GC908", "3": "GC909", "4": "GC910", "5": "GC1793"}
#: Characteristics (appendix 1) and regression flows (appendix 2) at the 285 gages regressed.
GAGES = json.loads(
    (ROOT / "fixtures" / "regression" / "pennsylvania_gages.json").read_text("utf-8")
)


@pytest.fixture(scope="module")
def lib():
    return load_state("PA")


def test_pennsylvania_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["1", "2", "3", "4", "5"]
    assert len(lib.equations) == 40
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region


def test_worked_example_tioga_river(lib):
    """Example 1: 01516350, region 1, 153 mi2, maximum elevation 2,446 ft, the 1-percent AEP."""
    est = evaluate(lib.equation("1", 0.01), {"DRNAREA": 153.0, "ELEVMAX": 2446.0})
    assert est.log_flow == pytest.approx(4.3666, abs=1e-4)
    assert est.flow_cfs == pytest.approx(23300, rel=0.005)


def test_published_regression_flows_reproduce_at_every_gage(lib):
    """Appendix 2 prints three significant figures."""
    assert len(GAGES) == 285
    for gage in GAGES.values():
        for i, aep in enumerate(lib.aeps(gage["region"])):
            eq = lib.equation(gage["region"], aep)
            inputs = {v.code: gage[v.code] for v in eq.variables}
            est = evaluate(eq, inputs, allow_extrapolation=True)
            assert est.flow_cfs == pytest.approx(gage["flows"][i], rel=0.006)


#: Where NSS differs from Table 3, confirmed by appendix 2.
NSS_COEFFICIENT_ERROR = {("2", 0.04): (-0.0813, -0.0831)}
NSS_SEP_ERROR = {("2", 0.1), ("5", 0.01)}


def test_every_equation_equals_the_nss_snapshot_but_for_its_errors(lib):
    snap = json.loads((DATA / "nss_snapshot" / "PA.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        key = (e.region_code, e.aep)
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == s["intercept"]
        assert e.coefficients[0] == s["coefficients"][0]
        if key in NSS_COEFFICIENT_ERROR:
            assert (e.coefficients[1], s["coefficients"][1]) == NSS_COEFFICIENT_ERROR[key]
        else:
            assert e.coefficients[1] == s["coefficients"][1]
        if key not in NSS_SEP_ERROR:
            assert e.sep_log == pytest.approx(s["sep_log"], abs=2.5e-4)
        for v, n in zip(e.variables, s["variables"]):
            assert (v.code, v.minimum, v.offset or None) == (
                n["code"],
                n["minimum"],
                n.get("offset"),
            )
            if not (e.region_code == "2" and v.code == "DRNAREA"):
                assert v.maximum == n["maximum"]


def test_region_2_keeps_the_reports_drainage_area_limit(lib):
    """NSS serves 1,160 mi2; Table 5 prints 1,610, and a region 2 gage is at 1,607."""
    assert lib.equation("2", 0.01).variables[0].maximum == 1610.0
    assert max(g["DRNAREA"] for g in GAGES.values() if g["region"] == "2") > 1600
