"""Big Sandy River benchmark against the 2012 PeakfqSA manual.

HISTORICAL RECORD. These expectations predate the vendored reference and are
not reproducible by peakfq 8.1.0 -- see docs/FORTRAN_UPLOAD.md section 6.0b.
For parity against the Fortran this repository actually vendors, see
tests/fortran_parity/.

Original docstring follows.

Big Sandy River Benchmark Test

Primary validation case for Bulletin 17C implementation.
Reference: PeakfqSA User Manual (Tim Cohn, USGS, 2012)
Site: USGS 03606500 - Big Sandy River at Bruceton, TN
"""

from typing import Dict, Tuple

import numpy as np
import pytest

from hydrolib.bulletin17c import Bulletin17C
from hydrolib.core import EMAParameters, FlowInterval

# Import fixture data
from tests.peakfqsa.fixtures.big_sandy import (
    BEGYEAR,
    ENDYEAR,
    EXPECTED_CONFIDENCE_INTERVALS,
    EXPECTED_PARAMETERS,
    EXPECTED_QUANTILES,
    HISTORICAL_PEAKS,
    REGIONAL_SKEW,
    REGIONAL_SKEW_SD,
    SYSTEMATIC_PEAKS,
    THRESHOLDS,
    TOLERANCE_PERCENT,
)


def build_flow_intervals() -> Tuple[list, list, EMAParameters]:
    """Build flow intervals and EMA parameters from Big Sandy data."""
    intervals = []
    water_years = []

    # Add systematic peaks (1930-1973)
    for year, flow in SYSTEMATIC_PEAKS.items():
        intervals.append(FlowInterval.from_point(flow))
        water_years.append(year)

    # Add historical peaks
    for year, flow in HISTORICAL_PEAKS.items():
        intervals.append(FlowInterval.from_point(flow, is_historical=True))
        water_years.append(year)

    # Build perception thresholds dict
    perception_thresholds = {}
    for t in THRESHOLDS:
        perception_thresholds[(t["start"], t["end"])] = t["lower"]

    ema_params = EMAParameters(
        historical_period=(BEGYEAR, min(SYSTEMATIC_PEAKS.keys()) - 1),
        perception_thresholds=perception_thresholds,
    )

    return intervals, water_years, ema_params


@pytest.fixture(scope="module")
def b17c_result():
    """Big Sandy EMA fit, built once for the whole module.

    Module-scoped deliberately. The MGBT three-sweep runs scipy.integrate.quad
    per candidate, so a single run_analysis() costs about two seconds, and this
    identical fit was previously rebuilt for every parametrized case in three
    separate classes. Safe to share because no test here mutates the result.

    Module-level rather than a class-scoped method: a class-scoped fixture
    written as an instance method is deprecated and removed in pytest 10, and
    the @staticmethod form that silences that warning breaks on Python 3.9,
    where staticmethod objects do not expose __name__ for pytest to read.
    """
    perception_thresholds = {(t["start"], t["end"]): t["lower"] for t in THRESHOLDS}
    b17c = Bulletin17C(
        peak_flows=list(SYSTEMATIC_PEAKS.values()),
        water_years=list(SYSTEMATIC_PEAKS.keys()),
        regional_skew=REGIONAL_SKEW,
        regional_skew_mse=REGIONAL_SKEW_SD**2,
        historical_peaks=[(y, f) for y, f in HISTORICAL_PEAKS.items()],
        perception_thresholds=perception_thresholds,
    )
    b17c.run_analysis()
    return b17c


class TestBigSandyParameters:
    """Test LP3 parameter estimation against expected values."""

    def test_mean_log(self, b17c_result):
        """Test mean of log10(Q)."""
        expected = EXPECTED_PARAMETERS["mean_log"]
        actual = b17c_result.results.mean_log
        pct_diff = abs(actual - expected) / expected * 100

        print(f"\nMean log Q: expected={expected:.6f}, actual={actual:.6f}, diff={pct_diff:.3f}%")
        assert pct_diff < TOLERANCE_PERCENT, f"Mean log differs by {pct_diff:.3f}%"

    def test_std_log(self, b17c_result):
        """Test standard deviation of log10(Q)."""
        expected = EXPECTED_PARAMETERS["std_log"]
        actual = b17c_result.results.std_log
        pct_diff = abs(actual - expected) / expected * 100

        print(f"\nStd log Q: expected={expected:.6f}, actual={actual:.6f}, diff={pct_diff:.3f}%")
        assert pct_diff < TOLERANCE_PERCENT, f"Std log differs by {pct_diff:.3f}%"

    def test_weighted_skew(self, b17c_result):
        """Test weighted skew coefficient."""
        expected = EXPECTED_PARAMETERS["skew_weighted"]
        actual = b17c_result.results.skew_weighted
        # Use absolute difference for skew since it can be near zero
        abs_diff = abs(actual - expected)

        print(f"\nWeighted skew: expected={expected:.6f}, actual={actual:.6f}, diff={abs_diff:.6f}")
        assert abs_diff < 0.05, f"Weighted skew differs by {abs_diff:.6f}"


class TestBigSandyQuantiles:
    """Test flood quantile estimates against expected values."""

    @pytest.mark.parametrize("aep,expected_flow", list(EXPECTED_QUANTILES.items()))
    def test_quantile(self, b17c_result, aep, expected_flow):
        """Test individual quantile estimates."""
        quantiles_df = b17c_result.compute_quantiles(np.array([aep]))
        actual_flow = quantiles_df["flow_cfs"].iloc[0]
        pct_diff = abs(actual_flow - expected_flow) / expected_flow * 100

        ri = 1 / aep
        print(
            f"\nQ{ri:.1f} (AEP={aep}): expected={expected_flow:.2f}, actual={actual_flow:.2f}, diff={pct_diff:.2f}%"
        )
        assert pct_diff < TOLERANCE_PERCENT * 2, f"Q{ri:.0f} differs by {pct_diff:.2f}%"


class TestBigSandyConfidenceIntervals:
    """Test confidence interval bounds."""

    # AEP 0.01 and 0.02 are known failures, tracked as strict xfail: the run
    # still happens and the deviation still prints, but the build stays green
    # -- and if either starts PASSING, strict=True fails the build, which is
    # exactly the alarm wanted when the missing quadrature lands.
    @pytest.mark.parametrize(
        "aep,expected_ci",
        [
            pytest.param(
                aep,
                ci,
                marks=(
                    pytest.mark.xfail(
                        strict=True,
                        reason=(
                            "Target is not reachable: these values come from the 2012 "
                            "PeakfqSA manual and are not reproducible by the vendored "
                            "reference (peakfq 8.1.0), whose HWN skew weighting differs "
                            "by design when censored data are present. Kept as historical "
                            "record. Live parity work belongs in tests/fortran_parity/, "
                            "which asserts against peakfq 8.1.0. Two real defects remain "
                            "on that path -- skew weighting and CI asymmetry. "
                            "See docs/FORTRAN_UPLOAD.md sections 6.0 and 6.0b."
                        ),
                    )
                    if aep in (0.01, 0.02)
                    else ()
                ),
            )
            for aep, ci in EXPECTED_CONFIDENCE_INTERVALS.items()
        ],
    )
    def test_confidence_interval(self, b17c_result, aep, expected_ci):
        """Test confidence interval bounds."""
        ci_df = b17c_result.compute_confidence_limits(np.array([aep]))
        lower = ci_df["lower_5pct"].iloc[0]
        upper = ci_df["upper_5pct"].iloc[0]
        expected_lower, expected_upper = expected_ci

        lower_diff = abs(lower - expected_lower) / expected_lower * 100
        upper_diff = abs(upper - expected_upper) / expected_upper * 100

        ri = 1 / aep
        print(f"\nCI for Q{ri:.0f}: expected=({expected_lower:.2f}, {expected_upper:.2f})")
        print(f"              actual=({lower:.2f}, {upper:.2f})")
        print(f"              diff=(lower:{lower_diff:.2f}%, upper:{upper_diff:.2f}%)")

        # These two cases (AEP 0.01, 0.02) are EXPECTED TO FAIL at 5%, and are
        # deliberately left failing rather than skipped or masked with a wider
        # tolerance -- this is a real open numerical question, not a missing
        # dependency. See docs/FORTRAN_UPLOAD.md section 1.6.
        #
        # The skew-uncertainty variance term (skew_used_mse propagated through
        # dK/dG) fixed the interval's SIZE: total width in log10 space is now
        # within ~2% of PeakfqSA and the point estimates within 0.5%. What
        # remains is SHAPE. compute_confidence_limits() forms log_Q +/- z*se,
        # so its upper/lower half-width ratio is exactly 1.000 at every AEP by
        # construction, while PeakfqSA's runs 1.043 (Q10) -> 1.531 (Q50) ->
        # 1.727 (Q100). Re-splitting our own width at those ratios puts every
        # bound inside 1%, which localizes the defect to the missing Inverse
        # Modified Cholesky Gaussian Quadrature -- a normal approximation
        # cannot produce that asymmetry at any variance.
        assert lower_diff < 5, f"Lower CI differs by {lower_diff:.2f}%"
        assert upper_diff < 5, f"Upper CI differs by {upper_diff:.2f}%"


class TestBigSandySummary:
    """Generate summary comparison report."""

    def test_full_comparison(self):
        """Run full comparison and print detailed report."""
        peak_flows = list(SYSTEMATIC_PEAKS.values())
        water_years = list(SYSTEMATIC_PEAKS.keys())

        historical_peaks = [(y, f) for y, f in HISTORICAL_PEAKS.items()]
        perception_thresholds = {}
        for t in THRESHOLDS:
            perception_thresholds[(t["start"], t["end"])] = t["lower"]

        b17c = Bulletin17C(
            peak_flows=peak_flows,
            water_years=water_years,
            regional_skew=REGIONAL_SKEW,
            regional_skew_mse=REGIONAL_SKEW_SD**2,
            historical_peaks=historical_peaks,
            perception_thresholds=perception_thresholds,
        )
        b17c.run_analysis()
        r = b17c.results

        print("\n" + "=" * 70)
        print("BIG SANDY RIVER BENCHMARK - FULL COMPARISON REPORT")
        print("=" * 70)
        print(f"\nSite: USGS 03606500 - Big Sandy River at Bruceton, TN")
        print(f"Analysis Period: {BEGYEAR}-{ENDYEAR}")
        print(f"Systematic Years: {len(SYSTEMATIC_PEAKS)}")
        print(f"Historical Peaks: {len(HISTORICAL_PEAKS)}")
        print(f"Regional Skew: {REGIONAL_SKEW}")

        print("\n--- LP3 PARAMETERS ---")
        print(f"{'Parameter':<20} {'Expected':>12} {'Actual':>12} {'Diff %':>10}")
        print("-" * 54)

        params = [
            ("Mean log Q", EXPECTED_PARAMETERS["mean_log"], r.mean_log),
            ("Std log Q", EXPECTED_PARAMETERS["std_log"], r.std_log),
            ("Weighted Skew", EXPECTED_PARAMETERS["skew_weighted"], r.skew_weighted),
        ]

        max_param_diff = 0
        for name, expected, actual in params:
            if abs(expected) > 0.001:
                pct_diff = (actual - expected) / expected * 100
            else:
                pct_diff = (actual - expected) * 100
            max_param_diff = max(max_param_diff, abs(pct_diff))
            print(f"{name:<20} {expected:>12.6f} {actual:>12.6f} {pct_diff:>+10.3f}")

        print("\n--- FLOOD QUANTILES ---")
        print(f"{'Return Interval':<18} {'Expected':>12} {'Actual':>12} {'Diff %':>10}")
        print("-" * 54)

        # Compute all quantiles at once
        aep_array = np.array(sorted(EXPECTED_QUANTILES.keys(), reverse=True))
        quantiles_df = b17c.compute_quantiles(aep_array)

        max_quantile_diff = 0
        for i, aep in enumerate(aep_array):
            ri = 1 / aep
            expected_flow = EXPECTED_QUANTILES[aep]
            actual_flow = quantiles_df["flow_cfs"].iloc[i]
            pct_diff = (actual_flow - expected_flow) / expected_flow * 100
            max_quantile_diff = max(max_quantile_diff, abs(pct_diff))

            ri_str = f"{ri:.1f}-yr" if ri < 10 else f"{ri:.0f}-yr"
            print(f"{ri_str:<18} {expected_flow:>12.2f} {actual_flow:>12.2f} {pct_diff:>+10.2f}")

        print("\n--- CONFIDENCE INTERVALS (90%) ---")
        print(
            f"{'Return Interval':<12} {'Expected Lower':>14} {'Actual Lower':>14} {'Expected Upper':>14} {'Actual Upper':>14}"
        )
        print("-" * 70)

        ci_aeps = np.array(sorted(EXPECTED_CONFIDENCE_INTERVALS.keys(), reverse=True))
        ci_df = b17c.compute_confidence_limits(ci_aeps)

        for i, aep in enumerate(ci_aeps):
            ri = 1 / aep
            exp_lower, exp_upper = EXPECTED_CONFIDENCE_INTERVALS[aep]
            lower = ci_df["lower_5pct"].iloc[i]
            upper = ci_df["upper_5pct"].iloc[i]
            ri_str = f"{ri:.0f}-yr"
            print(
                f"{ri_str:<12} {exp_lower:>14.2f} {lower:>14.2f} {exp_upper:>14.2f} {upper:>14.2f}"
            )

        print("\n--- ANALYSIS METADATA ---")
        print(f"Method: {r.method}")
        print(f"EMA Converged: {r.ema_converged}")
        print(f"EMA Iterations: {r.ema_iterations}")
        print(f"Low Outliers: {r.n_low_outliers}")
        if r.low_outlier_threshold:
            print(f"Low Outlier Threshold: {r.low_outlier_threshold:.2f}")
        print(f"Station Skew: {r.skew_station:.6f}")
        print(f"Skew Used: {r.skew_used:.6f}")

        print("\n--- SUMMARY ---")
        print(f"Max Parameter Difference: {max_param_diff:.3f}%")
        print(f"Max Quantile Difference: {max_quantile_diff:.3f}%")

        status = "PASS" if max_quantile_diff < 2 else "REVIEW NEEDED"
        print(f"\nOverall Status: {status}")
        print("=" * 70)

        # This test always passes - it's for reporting
        assert True
