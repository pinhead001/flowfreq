"""Vermont's verified equations (SIR 2025-5088 equations 3-10 and 21-28, Tables 2, 4, 5, 6)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

DATA = Path(__file__).resolve().parent.parent / "flowfreq" / "data" / "regression"
NSS = {"1": "GC2067", "2": "GC2068"}


@pytest.fixture(scope="module")
def lib():
    return load_state("VT")


def test_vermont_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["1", "2"]
    assert len(lib.equations) == 16
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region


def test_worked_example_wells_river(lib):
    """Appendix 7: 71.8 mi2, 6.49 percent wetland, 45.8 in., the 1-percent AEP."""
    eq = lib.equation("1", 0.01)
    est = evaluate(eq, {"DRNAREA": 71.8, "WETLAND": 6.49, "PRECPRIS20": 45.8})
    assert est.flow_cfs == pytest.approx(5730, rel=0.002)
    x = np.array([1.0, np.log10(71.8), 6.49, np.log10(45.8)])
    xux = float(x @ np.array(eq.covariance) @ x)
    assert xux == pytest.approx(0.00106, abs=1e-5)
    assert (eq.model_error_variance + xux) ** 0.5 == pytest.approx(0.1876, abs=5e-4)
    assert est.interval[1] == pytest.approx(11700, rel=0.005)


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "VT.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-9)
        assert list(e.coefficients) == s["coefficients"]
        assert e.sep_log == pytest.approx(s["sep_log"], abs=6e-4)
        for v, n in zip(e.variables, s["variables"]):
            assert v.code == n["code"]
            assert v.minimum == n["minimum"]
            assert v.maximum == pytest.approx(n["maximum"], abs=0.5)


def test_drainage_area_only_uses_the_average_error(lib):
    eq = lib.equation("2", 0.01)
    assert eq.covariance is None
    assert evaluate(eq, {"DRNAREA": 71.8}).flow_cfs == pytest.approx(224 * 71.8**0.807)
