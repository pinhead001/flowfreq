"""Alaska's verified equations (SIR 2016-5024 Tables 7 and 8)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.streamstats import RegressionRegionLocation
from flowfreq.ungaged import assemble

DATA = Path(__file__).resolve().parent.parent / "flowfreq" / "data" / "regression"


@pytest.fixture(scope="module")
def lib():
    return load_state("AK")


def _eq(lib, region, aep):
    return next(e for e in lib.equations if e.region_code == region and e.aep == aep)


def test_alaska_is_verified(lib):
    assert lib.status == "verified"
    assert sorted(lib.regions) == ["1", "2"]
    assert len(lib.equations) == 16


def test_worked_example_15129500(lib):
    """p. 27: DRNAREA 35.5 mi2, PRECPRIS00 170 in., the 1-percent AEP."""
    eq = _eq(lib, "1", 0.01)
    est = evaluate(eq, {"DRNAREA": 35.5, "PRECPRIS00": 170.0})
    assert est.flow_cfs == pytest.approx(7870, rel=0.001)
    x = np.array([1.0, np.log10(35.5), np.log10(170.0)])
    xux = float(x @ np.array(eq.covariance) @ x)
    assert xux == pytest.approx(0.001202, abs=5e-7)
    assert (0.082 + xux) ** 0.5 == pytest.approx(0.2884, abs=5e-5)
    low, high = est.interval
    assert low == pytest.approx(2630, rel=0.002)
    assert high == pytest.approx(23500, rel=0.002)


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "AK.json").read_text(encoding="utf-8"))
    by_key = {
        ("1" if s["region_code"] == "GC1656" else "2", s["aep"]): s for s in snap["equations"]
    }
    for e in lib.equations:
        s = by_key[(e.region_code, e.aep)]
        assert 10**e.intercept == pytest.approx(10 ** s["intercept"], rel=1e-12)
        assert list(e.coefficients) == s["coefficients"]
        assert [(v.minimum, v.maximum) for v in e.variables] == [
            (v["minimum"], v["maximum"]) for v in s["variables"]
        ]


def test_covariances_are_symmetric_and_positive_definite(lib):
    for e in lib.equations:
        u = np.array(e.covariance)
        assert np.array_equal(u, u.T)
        assert np.all(np.linalg.eigvalsh(u) > 0)


@pytest.mark.parametrize("area,code", [(35.5, "GC1656"), (5000.0, "GC1657")])
def test_drainage_area_picks_the_class(lib, area, code):
    regions = [
        RegressionRegionLocation(0, "GC1656", "small", 100.0, area),
        RegressionRegionLocation(0, "GC1657", "large", 100.0, area),
    ]
    a = assemble(
        lib,
        regions,
        {"DRNAREA": area, "PRECPRIS00": 60.0},
        aeps=[0.01],
        allow_extrapolation=False,
        weighting_space="log",
    )
    assert [r.code for r in a.used] == [code]
