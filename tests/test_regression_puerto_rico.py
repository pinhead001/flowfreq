"""Puerto Rico's verified equations (SIR 2021-5062 Table 3)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

DATA = Path(__file__).resolve().parent.parent / "flowfreq" / "data" / "regression"
NSS = {"1": "GC1863", "2": "GC1862"}


@pytest.fixture(scope="module")
def lib():
    return load_state("PR")


def test_puerto_rico_is_verified(lib):
    assert lib.status == "verified"
    assert lib.regions == ["1", "2"]
    assert len(lib.equations) == 16
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region


def test_worked_example_rio_valenciano(lib):
    """Example 1: 50056400, region 2, the 1-percent AEP: 39,758 ft3/s, printed as 39,800.

    The report rounds the drainage area to 16.5 mi2 in the text; 39,758 comes from 16.49.
    """
    est = evaluate(lib.equation("2", 0.01), {"DRNAREA": 16.49})
    assert est.flow_cfs == pytest.approx(39758, rel=2e-4)
    assert evaluate(lib.equation("2", 0.01), {"DRNAREA": 16.5}).flow_cfs == pytest.approx(
        39800, rel=0.001
    )


def test_every_equation_equals_the_nss_snapshot_but_one_error(lib):
    """NSS gives region 2 at 0.2 percent region 1's SEP, 18.7 percent; Tables 3 and 4 print 39.3."""
    snap = json.loads((DATA / "nss_snapshot" / "PR.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-9)
        assert list(e.coefficients) == s["coefficients"]
        v, n = e.variables[0], s["variables"][0]
        assert (v.code, v.minimum, v.maximum) == (n["code"], n["minimum"], n["maximum"])
        if (e.region_code, e.aep) == ("2", 0.002):
            assert s["sep_log"] == pytest.approx(0.08052, abs=1e-4)
            assert e.sep_log == pytest.approx(0.1648, abs=1e-3)
        else:
            assert e.sep_log == pytest.approx(s["sep_log"], abs=6e-5)
