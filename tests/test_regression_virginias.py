"""Virginia's and West Virginia's verified equations (SIR 2025-5110 Tables 6, 7 and 8)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

DATA = Path(__file__).resolve().parent.parent / "flowfreq" / "data" / "regression"
NSS = {
    "VA": {
        "SHIGH": "GC2087",
        "SVYRD": "GC2088",
        "NHIGH": "GC2089",
        "SBLU": "GC2090",
        "CVYRD": "GC2091",
        "PDMT": "GC2092",
        "NBLU": "GC2093",
        "NVYRD": "GC2094",
        "NPDMT": "GC2095",
        "PLNS": "GC2096",
        "DLMRV": "GC2097",
    },
    "WV": {
        "WPLT": "GC2081",
        "NHIGH": "GC2082",
        "SHIGH": "GC2083",
        "SVYRD": "GC2084",
        "SBLU": "GC2085",
        "NVYRD": "GC2086",
    },
}


@pytest.mark.parametrize("state", ["VA", "WV"])
def test_state_is_verified(state):
    lib = load_state(state)
    assert lib.status == "verified"
    assert set(lib.regions) == set(NSS[state])
    assert len(lib.equations) == 8 * len(NSS[state])
    for region, code in NSS[state].items():
        assert region_for_nss_code(lib, code) == region


def test_worked_example_ellis_creek():
    """Example 3: Piedmont, 7.54 mi2, the 1-percent AEP. The report rounds at each step, to 3.46
    in log units and then 2,880 ft3/s."""
    est = evaluate(load_state("VA").equation("PDMT", 0.01), {"DRNAREA": 7.54})
    assert round(est.log_flow, 2) == 3.46
    assert est.flow_cfs == pytest.approx(2880, rel=0.012)


def test_worked_example_little_coal_river():
    """Example 1: 03199400, Southern Highlands, 318 mi2, the 1-percent AEP: 31,200 ft3/s."""
    est = evaluate(load_state("WV").equation("SHIGH", 0.01), {"DRNAREA": 318.0})
    assert est.flow_cfs == pytest.approx(31200, rel=0.01)


def test_shared_regions_are_the_same_equations():
    va, wv = load_state("VA"), load_state("WV")
    for region in set(NSS["VA"]) & set(NSS["WV"]):
        for aep in va.aeps(region):
            a, w = va.equation(region, aep), wv.equation(region, aep)
            assert (a.intercept, a.coefficients) == (w.intercept, w.coefficients)
            assert (a.avp, a.model_error_variance) == (w.avp, w.model_error_variance)


@pytest.mark.parametrize("state", ["VA", "WV"])
def test_every_equation_equals_the_nss_snapshot(state):
    lib = load_state(state)
    snap = json.loads((DATA / "nss_snapshot" / f"{state}.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[state][e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[state][e.region_code], e.aep)]
        assert e.intercept == s["intercept"]
        assert list(e.coefficients) == s["coefficients"]
        assert e.sep_log == pytest.approx(s["sep_log"], abs=2.5e-4)
        v, n = e.variables[0], s["variables"][0]
        assert v.code == n["code"]
        assert v.minimum == pytest.approx(n["minimum"], abs=0.01)
        assert v.maximum == pytest.approx(n["maximum"], abs=0.5)
