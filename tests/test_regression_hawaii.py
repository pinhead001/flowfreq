"""Hawaii's verified equations (SIR 2023-5014 Tables 5 and 6)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state

DATA = Path(__file__).resolve().parent.parent / "flowfreq" / "data" / "regression"
GC = {f"GC{1917 + i}": str(i) for i in range(1, 11)}


@pytest.fixture(scope="module")
def lib():
    return load_state("HI")


def test_hawaii_is_verified(lib):
    assert lib.status == "verified"
    assert sorted(lib.regions, key=int) == [str(i) for i in range(1, 11)]
    assert len(lib.equations) == 80


def test_worked_example_kawaikoi_stream(lib):
    """p. 34: 16010000, region 1, DRNAREA 3.82 mi2, PRECIP 133 in.: Q1% = 7,240 ft3/s."""
    eq = next(e for e in lib.equations if e.region_code == "1" and e.aep == 0.01)
    est = evaluate(eq, {"DRNAREA": 3.82, "PRECIP": 133.0})
    assert est.flow_cfs == pytest.approx(7240, rel=0.001)


def test_printed_coefficients_are_nss_rounded_to_three_decimals(lib):
    snap = json.loads((DATA / "nss_snapshot" / "HI.json").read_text(encoding="utf-8"))
    by_key = {(GC[s["region_code"]], s["aep"]): s for s in snap["equations"]}
    for e in lib.equations:
        s = by_key[(e.region_code, e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=5.0001e-4)
        nss = dict(zip([v["code"] for v in s["variables"]], s["coefficients"]))
        assert set(nss) == {v.code for v in e.variables}
        for v, c in zip(e.variables, e.coefficients):
            assert c == pytest.approx(nss[v.code], abs=5.0001e-4)
        limits = {v["code"]: (v["minimum"], v["maximum"]) for v in s["variables"]}
        for v in e.variables:
            assert v.minimum == pytest.approx(limits[v.code][0], abs=0.0101)
            assert v.maximum == pytest.approx(limits[v.code][1], rel=0.0035, abs=0.0101)


def test_region_10_has_no_drainage_area(lib):
    eqs = [e for e in lib.equations if e.region_code == "10"]
    assert {tuple(v.code for v in e.variables) for e in eqs} == {("I48H500Y",)}


def test_no_covariance_so_intervals_use_sep(lib):
    eq = next(e for e in lib.equations if e.region_code == "3" and e.aep == 0.01)
    assert eq.covariance is None and eq.n_sites == 43
    est = evaluate(eq, {"DRNAREA": 5.0, "PRECIP": 100.0, "PERM24IN": 3.0})
    assert est.interval is not None and est.interval[0] < est.flow_cfs < est.interval[1]
