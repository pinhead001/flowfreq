"""Kansas's verified equations (SIR 2017-5063 version 1.1, Tables 3, 4, 5, 7 and 8)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from flowfreq.regression import OutOfRangeError, evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {"1": "GC1743", "2": "GC1744"}
#: Regression estimates (Tables 3 and 4) and characteristics (Table 5) at each region's gages.
GAGES = json.loads((ROOT / "fixtures" / "regression" / "kansas_gages.json").read_text("utf-8"))


@pytest.fixture(scope="module")
def lib():
    return load_state("KS")


def test_kansas_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["1", "2"]
    assert len(lib.equations) == 16
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region


@pytest.mark.parametrize("region,count,tolerance", [("1", 143, 0.015), ("2", 27, 0.005)])
def test_published_regression_estimates_reproduce(lib, region, count, tolerance):
    """The report prints no worked example. The tables print three significant figures."""
    assert len(GAGES[region]) == count
    for gage in GAGES[region].values():
        inputs = {"CONTDA": gage["contda"]}
        if region == "1":
            inputs["PRECPRIS10"] = gage["precip"]
        for i, aep in enumerate(lib.aeps(region)):
            est = evaluate(lib.equation(region, aep), inputs)
            assert est.flow_cfs == pytest.approx(gage["flows"][i], rel=tolerance)


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "KS.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-9)
        assert list(e.coefficients) == s["coefficients"]
        assert e.sep_log == pytest.approx(s["sep_log"], abs=6e-5)
        assert [v.code for v in e.variables] == [v["code"] for v in s["variables"]]


def test_region_2_keeps_the_reports_range_not_nss(lib):
    """NSS serves region 2 up to 14,600 mi2; the report's gages stop at 3,555."""
    eq = lib.equation("2", 0.01)
    assert evaluate(eq, {"CONTDA": 3555.0}).flow_cfs > 0
    with pytest.raises(OutOfRangeError):
        evaluate(eq, {"CONTDA": 5000.0})
