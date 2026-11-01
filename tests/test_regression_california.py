"""California's verified equations (SIR 2012-5113 Tables 5, 6, 8 and 9)."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state

DATA = Path(__file__).resolve().parent.parent / "flowfreq" / "data" / "regression"
GC = {"GC1529": "1", "GC1530": "2", "GC1531": "3", "GC1532": "4", "GC1533": "5", "GC1534": "6"}


@pytest.fixture(scope="module")
def lib():
    return load_state("CA")


def _eq(lib, region, aep):
    return next(e for e in lib.equations if e.region_code == region and e.aep == aep)


def test_california_is_verified(lib):
    assert lib.status == "verified"
    assert sorted(lib.regions) == ["1", "2", "3", "4", "5", "6"]
    assert len(lib.equations) == 48


@pytest.mark.parametrize("area,precip,flow", [(27.8, 32.7, 6000), (25.0, 32.5, 5460)])
def test_worked_examples_north_coast(lib, area, precip, flow):
    """p. 25: the North Coast 1-percent AEP flow, printed to 3 significant figures."""
    est = evaluate(_eq(lib, "1", 0.01), {"DRNAREA": area, "PRECIP": precip})
    assert est.flow_cfs == pytest.approx(flow, rel=0.001)


@pytest.mark.parametrize(
    "region,t", [("1", 1.653), ("2", 1.675), ("3", 1.653), ("4", 1.662), ("5", 1.661)]
)
def test_site_counts_reproduce_the_printed_t_values(lib, region, t):
    from scipy import stats

    eq = _eq(lib, region, 0.01)
    df = eq.n_sites - len(eq.variables) - 1
    assert stats.t.ppf(0.95, df) == pytest.approx(t, abs=5e-4)


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "CA.json").read_text(encoding="utf-8"))
    by_key = {
        (GC[s["region_code"]], s["aep"]): s for s in snap["equations"] if s["region_code"] in GC
    }
    for e in lib.equations:
        s = by_key[(e.region_code, e.aep)]
        assert 10**e.intercept == pytest.approx(10 ** s["intercept"], rel=1e-12)
        assert dict(zip([v.code for v in e.variables], e.coefficients)) == dict(
            zip([v["code"] for v in s["variables"]], s["coefficients"])
        )
        limits = {v["code"]: (v["minimum"], v["maximum"]) for v in s["variables"]}
        assert {v.code: (v.minimum, v.maximum) for v in e.variables} == limits


def test_covariances_are_symmetric_and_positive_definite(lib):
    for e in lib.equations:
        if e.region_code == "6":
            assert e.covariance is None and e.model_error_variance is None
            continue
        u = np.array(e.covariance)
        assert u.shape == (len(e.variables) + 1,) * 2
        assert np.array_equal(u, u.T)
        assert np.all(np.linalg.eigvalsh(u) > 0)


def test_variance_of_prediction_at_a_table_10_gage(lib):
    """Table 10, 10291500 (map 36, Buckeye Creek, Lahontan; Table 2: 43.8 mi2, 37.466 in.):
    regression VP 0.1216 at 50 percent and 0.083296 at 1 percent."""
    x = np.array([1.0, np.log10(43.8), np.log10(37.466)])
    for aep, expected in ((0.5, 0.1216), (0.01, 0.083296)):
        e = _eq(lib, "2", aep)
        vp = e.model_error_variance + x @ np.array(e.covariance) @ x
        assert vp == pytest.approx(expected, abs=6e-4)
