"""Maine's verified equations (SIR 2020-5092 Tables 1, 3 and 4, appendix tables 1.4, 1.6, 1.8)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {"1": "GC1798", "2": "GC1799"}
RI = [2, 5, 10, 25, 50, 100, 200, 500]
#: Characteristics (table 1.6: area, wetlands, eight rainfalls) and regression flows (table 1.4)
#: at the 96 gages found in both tables.
GAGES = json.loads((ROOT / "fixtures" / "regression" / "maine_gages.json").read_text("utf-8"))


@pytest.fixture(scope="module")
def lib():
    return load_state("ME")


def test_maine_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["1", "2"]
    assert len(lib.equations) == 16
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region


def test_published_regression_flows_reproduce_at_every_gage(lib):
    """The report prints no worked example. Table 1.4 prints three significant figures."""
    assert len(GAGES) == 96
    for gage in GAGES.values():
        c = gage["chars"]
        for i, aep in enumerate(lib.aeps("1")):
            inputs = {"DRNAREA": c[0], "STORAGE": c[1], f"I24H{RI[i]}Y": c[2 + i]}
            est = evaluate(lib.equation("1", aep), inputs, allow_extrapolation=True)
            assert est.flow_cfs == pytest.approx(gage["flows"][i], rel=0.01)


def test_covariance_reproduces_the_published_average_error(lib):
    """No worked example checks the covariance, so check it against Table 3's error."""
    for i, aep in enumerate(lib.aeps("1")):
        eq = lib.equation("1", aep)
        u = np.array(eq.covariance)
        rows = [
            np.array([1.0, np.log10(g["chars"][0]), g["chars"][1], g["chars"][2 + i]])
            for g in GAGES.values()
        ]
        variance = np.mean([eq.model_error_variance + x @ u @ x for x in rows])
        assert variance**0.5 == pytest.approx(eq.sep_log, abs=0.005)


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "ME.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-6)
        assert list(e.coefficients) == s["coefficients"]
        assert e.sep_log == pytest.approx(s["sep_log"], abs=6e-4)
        assert [(v.code, v.minimum, v.maximum) for v in e.variables] == [
            (v["code"], v["minimum"], v["maximum"]) for v in s["variables"]
        ]


def test_drainage_area_only_uses_the_average_error(lib):
    eq = lib.equation("2", 0.01)
    assert eq.covariance is None
    assert evaluate(eq, {"DRNAREA": 100.0}).interval_method.startswith("average standard error")
