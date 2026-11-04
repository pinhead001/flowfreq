"""Connecticut's verified equations (SIR 2020-5054 Tables 7-10)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

DATA = Path(__file__).resolve().parent.parent / "flowfreq" / "data" / "regression"
NSS = {"1": "GC1801", "2": "GC1800"}


@pytest.fixture(scope="module")
def lib():
    return load_state("CT")


def test_connecticut_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["1", "2"]
    assert len(lib.equations) == 16
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region


def test_worked_example_quinnipiac_river(lib):
    """p. 27-28: 01196500, 110 mi2, I24H100Y 8.53 in., SOILCorD 0.4207, the 1-percent AEP."""
    eq = lib.equation("1", 0.01)
    est = evaluate(eq, {"DRNAREA": 110.0, "I24H100Y": 8.53, "SSURGOCCDD": 0.4207})
    assert est.flow_cfs == pytest.approx(13300, rel=0.005)
    # Table 9 states the order constant, DRNAREA, I24H, SOILCorD. Of the 24 orders only this one
    # gives the printed 0.00286; the stated order gives 1.14.
    x = np.array([1.0, np.log10(110.0), 1.4207, 8.53])
    xux = float(x @ np.array(eq.covariance) @ x)
    assert xux == pytest.approx(0.00286, abs=1e-5)
    assert (eq.model_error_variance + xux) ** 0.5 == pytest.approx(0.15446, abs=5e-5)
    low, high = est.interval
    assert low == pytest.approx(7360, rel=0.005)
    assert high == pytest.approx(24000, rel=0.005)


def test_every_equation_equals_the_nss_snapshot(lib):
    """NSS folds SOILCorD's added 1 into its intercept."""
    snap = json.loads((DATA / "nss_snapshot" / "CT.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        folded = e.intercept + sum(c * v.offset for c, v in zip(e.coefficients, e.variables))
        assert folded == pytest.approx(s["intercept"], abs=1e-9)
        nss = {
            v["code"]: (c, v["minimum"], v["maximum"])
            for v, c in zip(s["variables"], s["coefficients"])
        }
        assert {
            v.code: (c, v.minimum, v.maximum) for v, c in zip(e.variables, e.coefficients)
        } == nss
        assert e.sep_log == pytest.approx(s["sep_log"], abs=1e-3)


def test_every_covariance_gives_a_plausible_prediction_error(lib):
    """The order fixed from the 1-percent example has to work at every AEP, at every range end."""
    for e in lib.equations:
        if e.region_code != "1":
            continue
        for pick in (0, 1):
            inputs = {v.code: (v.minimum, v.maximum)[pick] for v in e.variables}
            est = evaluate(e, inputs)
            assert e.model_error_variance < est.prediction_variance < 0.1


def test_drainage_area_only_has_no_covariance(lib):
    eq = lib.equation("2", 0.01)
    assert eq.covariance is None
    est = evaluate(eq, {"DRNAREA": 110.0})
    assert est.interval_method == "average variance of prediction"
