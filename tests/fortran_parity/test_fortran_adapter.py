"""Tests for ``flowfreq.fortran_engine``'s ``ReferenceResult`` -> ``FrequencyResults``
adapter (design doc section 4).

Runs everywhere -- built from a synthetic ``ReferenceResult``/``EmafitArrays``
pair, not a live Fortran call, so it does not need the built extension. The
one thing this module cannot check is whether the adapter reads real
``emafitpr`` output correctly; that needs ``requires_fortran``-gated
end-to-end coverage this environment cannot run (no gfortran/meson here) and
is left for a follow-up once the extension can be built.
"""

from __future__ import annotations

import numpy as np

from flowfreq.core import AnalysisMethod
from flowfreq.fortran_engine import EmafitArrays, _frequency_results_from_reference
from flowfreq.validation.reference import ReferenceResult


def _synthetic_reference() -> ReferenceResult:
    return ReferenceResult(
        source="synthetic (unit test)",
        station_name="Test Site",
        n_peaks=5,
        n_systematic=5,
        n_historical=0,
        low_outlier_count=1,
        low_outlier_threshold=100.0,
        parameters={
            "mean_log": 3.5,
            "std_log": 0.30,
            "skew_weighted": -0.10,
            "mean_log_at_site": 3.5,
            "std_log_at_site": 0.30,
            "skew_at_site": -0.05,
            "mse_skew": 0.09,
        },
        quantiles={0.5: 3000.0, 0.01: 15000.0},
        confidence_intervals={0.5: (2800.0, 3200.0), 0.01: (13000.0, 17000.0)},
        variance={0.5: 0.001, 0.01: 0.002},
    )


def _synthetic_arrays() -> EmafitArrays:
    years = np.array([2000, 2001, 2002, 2003, 2004])
    peaks = {2000: 50.0, 2001: 500.0, 2002: 800.0, 2003: 1200.0, 2004: 900.0}
    ql = qu = np.log10(np.array([peaks[y] for y in years]))
    tl = np.full(5, np.log10(1e-20))
    tu = np.full(5, np.log10(1e20))
    dtype = np.zeros(5, dtype=np.int32)
    return EmafitArrays(
        ql=ql,
        qu=qu,
        tl=tl,
        tu=tu,
        dtype=dtype,
        years=years,
        systematic_peaks=peaks,
        n_zeros=0,
        n_censored=0,
        gbthrsh0=-99.0,
    )


class TestNeverFabricates:
    """Design doc section 4's rule, and section 7's acceptance test for it."""

    def test_ema_iterations_and_converged_are_none(self):
        results = _frequency_results_from_reference(
            _synthetic_reference(), _synthetic_arrays(), regional_skew=-0.3, regional_skew_mse=0.3
        )
        assert results.ema_iterations is None
        assert results.ema_converged is None

    def test_method_is_ema(self):
        results = _frequency_results_from_reference(
            _synthetic_reference(), _synthetic_arrays(), regional_skew=-0.3, regional_skew_mse=0.3
        )
        assert results.method == AnalysisMethod.EMA


class TestSkewWeightingPolicy:
    def test_weighted_skew_reported_when_regional_supplied(self):
        results = _frequency_results_from_reference(
            _synthetic_reference(), _synthetic_arrays(), regional_skew=-0.3, regional_skew_mse=0.3
        )
        assert results.skew_weighted == -0.10
        assert results.skew_used == -0.10
        assert results.skew_regional == -0.3

    def test_no_weighted_skew_when_no_regional_supplied(self):
        """Station-only (no regional skew given): emafitpr's own 'weighted'
        column reduces to the at-site fit, and reporting it as skew_weighted
        would claim a weighting that never happened."""
        results = _frequency_results_from_reference(
            _synthetic_reference(), _synthetic_arrays(), regional_skew=None, regional_skew_mse=None
        )
        assert results.skew_weighted is None
        assert results.skew_used == results.skew_station == -0.05
        assert results.skew_regional is None


class TestDerivedFields:
    def test_pilf_flows_derived_from_raw_peaks_below_threshold(self):
        results = _frequency_results_from_reference(
            _synthetic_reference(), _synthetic_arrays(), regional_skew=None, regional_skew_mse=None
        )
        # threshold is 100.0; only the 50.0 cfs peak is below it.
        assert results.pilf_flows == [50.0]

    def test_no_pilfs_when_threshold_is_zero(self):
        ref = _synthetic_reference()
        ref.low_outlier_threshold = 0.0
        ref.low_outlier_count = 0
        results = _frequency_results_from_reference(
            ref, _synthetic_arrays(), regional_skew=None, regional_skew_mse=None
        )
        assert results.pilf_flows == []

    def test_n_censored_comes_from_the_arrays_not_the_reference(self):
        arrays = _synthetic_arrays()
        results = _frequency_results_from_reference(
            _synthetic_reference(), arrays, regional_skew=None, regional_skew_mse=None
        )
        assert results.n_censored == arrays.n_censored == 0

    def test_mgb_critical_value_is_recomputed_natively(self):
        """Not an emafitpr output -- design doc section 4's one exception to
        "never synthesise": grubbs_beck_critical_value(n) is a deterministic
        function of sample size alone, computed the same way the native
        engine already computes this exact diagnostic field for itself."""
        from flowfreq.core import grubbs_beck_critical_value

        results = _frequency_results_from_reference(
            _synthetic_reference(), _synthetic_arrays(), regional_skew=None, regional_skew_mse=None
        )
        assert results.mgb_critical_value == grubbs_beck_critical_value(5)

    def test_quantiles_and_confidence_limits_round_trip(self):
        results = _frequency_results_from_reference(
            _synthetic_reference(), _synthetic_arrays(), regional_skew=None, regional_skew_mse=None
        )
        q = results.quantiles.set_index("aep")["flow_cfs"]
        assert q.loc[0.5] == 3000.0
        assert q.loc[0.01] == 15000.0
        ci = results.confidence_limits.set_index("aep")
        assert ci.loc[0.5, "lower_5pct"] == 2800.0
        assert ci.loc[0.5, "upper_5pct"] == 3200.0
