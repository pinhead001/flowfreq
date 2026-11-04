"""Massachusetts's verified equations (SIR 2016-5156, equations 2-9, Tables 5, 6 and 8)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state

DATA = Path(__file__).resolve().parent.parent / "flowfreq" / "data" / "regression"
EXAMPLE = {"DRNAREA": 18.0, "ELEV": 689.0, "LC06STOR": 8.0}


@pytest.fixture(scope="module")
def lib():
    return load_state("MA")


def test_massachusetts_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["1"]
    assert len(lib.equations) == 8


def test_worked_example(lib):
    """pp. 24-25: 18 mi2, 689 ft, 8 percent storage, the 1-percent AEP."""
    eq = lib.equation("1", 0.01)
    est = evaluate(eq, EXAMPLE)
    assert est.flow_cfs == pytest.approx(1967, abs=0.5)
    x = np.array([1.0, np.log10(18.0), 689.0 / 3.28084 * 0.001, 8.0 * 0.01])
    xux = float(x @ np.array(eq.covariance) @ x)
    assert xux == pytest.approx(0.00298, abs=5e-6)
    assert (eq.model_error_variance + xux) ** 0.5 == pytest.approx(0.21210, abs=5e-5)
    # The report rounds t to 1.65 and prints 879 and 4,400; flowfreq uses t at 195 degrees of
    # freedom (1.653).
    low, high = est.interval
    assert low == pytest.approx(879, rel=0.005)
    assert high == pytest.approx(4400, rel=0.005)


def test_every_equation_equals_the_nss_snapshot(lib):
    """NSS folds the report's scales into its coefficients; its ELEV range is from meters."""
    snap = json.loads((DATA / "nss_snapshot" / "MA.json").read_text(encoding="utf-8"))
    by_aep = {s["aep"]: s for s in snap["equations"]}
    assert {e.aep for e in lib.equations} == set(by_aep)
    for e in lib.equations:
        s = by_aep[e.aep]
        assert e.intercept == s["intercept"]
        scaled = [c * v.scale for c, v in zip(e.coefficients, e.variables)]
        assert scaled == pytest.approx(s["coefficients"], rel=1e-9)
        assert e.sep_log == pytest.approx(s["sep_log"], abs=1e-6)
        for v, n in zip(e.variables, s["variables"]):
            assert v.code == n["code"]
            assert v.minimum == pytest.approx(n["minimum"], abs=0.5)
            assert v.maximum == pytest.approx(n["maximum"], abs=1.0)


def test_nss_region_code_finds_the_set(lib):
    from flowfreq.ungaged import region_for_nss_code

    assert region_for_nss_code(lib, "GC1741") == "1"


def test_flows_rise_as_floods_get_rarer(lib):
    flows = [evaluate(lib.equation("1", a), EXAMPLE).flow_cfs for a in lib.aeps("1")]
    assert flows == sorted(flows)
