"""Tests for the variance of quantile estimates and Bulletin 17C Appendix 9 weighting."""

from __future__ import annotations

import dataclasses
import math

import numpy as np
import pytest
from scipy.special import ndtri

from flowfreq.bulletin17c import Bulletin17C
from flowfreq.gage_weighting import weight_with_regression
from flowfreq.regression.equations import RegressionEquation, evaluate

AEPS = np.array([0.5, 0.1, 0.01])
COVARIANCE = [[0.01, -0.002], [-0.002, 0.001]]


@pytest.fixture(scope="module")
def ema():
    rng = np.random.default_rng(7)
    flows = 10 ** rng.normal(3.2, 0.25, 45)
    b = Bulletin17C(flows, np.arange(1975, 2020), regional_skew=-0.1, regional_skew_mse=0.12)
    b.run_analysis(method="ema")
    return b


def _equation(aep: float, **kw) -> RegressionEquation:
    base = dict(
        region_code="R2",
        aep=aep,
        intercept=1.0,
        variables=[{"code": "DRNAREA", "transform": "log10", "minimum": 1.0, "maximum": 1000.0}],
        coefficients=[0.8],
        citation={"publication": "synthetic", "table": "none"},
        model_error_variance=0.04,
        covariance=COVARIANCE,
        n_sites=80,
    )
    base.update(kw)
    return RegressionEquation.from_dict(base)


class TestQuantileVariance:
    def test_ema_column_is_positive_and_grows_toward_rare_floods(self, ema):
        cl = ema.compute_confidence_limits(AEPS)
        v = cl["variance_log"].to_numpy()
        assert (v > 0).all()
        assert v[0] < v[1] < v[2]

    def test_method_of_moments_variance_reproduces_its_symmetric_bounds(self):
        rng = np.random.default_rng(8)
        flows = 10 ** rng.normal(3.0, 0.3, 40)
        b = Bulletin17C(flows, np.arange(1980, 2020), regional_skew=0.0, regional_skew_mse=0.3)
        b.run_analysis(method="mom")
        cl = b.compute_confidence_limits(AEPS)
        half = np.log10(cl["upper_5pct"]) - np.log10(cl["flow_cfs"])
        np.testing.assert_allclose(half, ndtri(0.95) * np.sqrt(cl["variance_log"]))

    def test_variance_does_not_depend_on_the_confidence_level(self, ema):
        a = ema.compute_confidence_limits(AEPS, confidence=0.90)["variance_log"]
        b = ema.compute_confidence_limits(AEPS, confidence=0.95)["variance_log"]
        np.testing.assert_allclose(a, b, rtol=1e-6)

    def test_fortran_frames_carry_var_est(self):
        from flowfreq.fortran_engine import quantile_frames
        from flowfreq.validation.reference import ReferenceResult

        ref = ReferenceResult(
            source="test",
            parameters={"mean_log": 3.0, "std_log": 0.2},
            quantiles={0.1: 2000.0, 0.01: 3000.0},
            confidence_intervals={0.1: (1800.0, 2300.0), 0.01: (2500.0, 4000.0)},
            variance={0.1: 0.002, 0.01: 0.006},
        )
        _, cl = quantile_frames(ref)
        assert cl.set_index("aep")["variance_log"].to_dict() == {0.1: 0.002, 0.01: 0.006}


class TestPredictionVariance:
    def test_site_specific(self):
        est = evaluate(_equation(0.01), {"DRNAREA": 100.0})
        x = np.array([1.0, 2.0])
        assert est.prediction_variance == pytest.approx(0.04 + x @ np.array(COVARIANCE) @ x)

    def test_none_without_error_statistics(self):
        est = evaluate(
            _equation(0.01, model_error_variance=None, covariance=None), {"DRNAREA": 100.0}
        )
        assert est.prediction_variance is None


class TestWeightWithRegression:
    def test_inverse_variance_formula(self, ema):
        regs = [evaluate(_equation(a), {"DRNAREA": 100.0}) for a in AEPS]
        out = weight_with_regression(ema.results, regs)
        cl = ema.results.confidence_limits.set_index("aep")
        assert set(out["aep"]) == set(cl.index)
        out = out.set_index("aep")
        for reg in regs:
            aep = reg.equation.aep
            row = out.loc[aep]
            vs = cl.loc[aep, "variance_log"]
            qs = math.log10(cl.loc[aep, "flow_cfs"])
            vr = reg.prediction_variance
            expected = (vr * qs + vs * reg.log_flow) / (vs + vr)
            assert math.log10(row["weighted_flow_cfs"]) == pytest.approx(expected)
            assert row["weighted_variance_log"] == pytest.approx(vs * vr / (vs + vr))
            assert row["weight_at_site"] == pytest.approx(vr / (vs + vr))
            assert row["weighted_lower_cfs"] < row["weighted_flow_cfs"] < row["weighted_upper_cfs"]
            half = math.log10(row["weighted_upper_cfs"] / row["weighted_flow_cfs"])
            assert half == pytest.approx(ndtri(0.95) * math.sqrt(row["weighted_variance_log"]))

    def test_aep_without_a_regression_is_blank_with_a_reason(self, ema):
        regs = [evaluate(_equation(0.01), {"DRNAREA": 100.0})]
        out = weight_with_regression(ema.results, regs).set_index("aep")
        assert math.isnan(out.loc[0.5, "weighted_flow_cfs"])
        assert "no regression" in out.loc[0.5, "note"]
        assert not math.isnan(out.loc[0.01, "weighted_flow_cfs"])

    def test_extrapolated_regression_is_refused_unless_allowed(self, ema):
        reg = evaluate(_equation(0.01), {"DRNAREA": 5000.0}, allow_extrapolation=True)
        out = weight_with_regression(ema.results, [reg]).set_index("aep")
        assert math.isnan(out.loc[0.01, "weighted_flow_cfs"])
        assert "calibrated range" in out.loc[0.01, "note"]
        out = weight_with_regression(ema.results, [reg], allow_extrapolated=True)
        assert not math.isnan(out.set_index("aep").loc[0.01, "weighted_flow_cfs"])

    def test_regression_without_variance_is_blank_with_a_reason(self, ema):
        reg = evaluate(
            _equation(0.01, model_error_variance=None, covariance=None), {"DRNAREA": 100.0}
        )
        out = weight_with_regression(ema.results, [reg]).set_index("aep")
        assert math.isnan(out.loc[0.01, "weighted_flow_cfs"])
        assert "variance" in out.loc[0.01, "note"]

    def test_needs_the_variance_column(self, ema):
        stale = ema.results.confidence_limits.drop(columns=["variance_log"])
        results = dataclasses.replace(ema.results, confidence_limits=stale)
        with pytest.raises(ValueError, match="variance_log"):
            weight_with_regression(results, [evaluate(_equation(0.01), {"DRNAREA": 100.0})])

    def test_duplicate_regression_aep_raises(self, ema):
        regs = [evaluate(_equation(0.01), {"DRNAREA": 100.0})] * 2
        with pytest.raises(ValueError, match="more than one"):
            weight_with_regression(ema.results, regs)
