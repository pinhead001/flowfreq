"""Tests for flowfreq.core LP3 utility functions not already covered elsewhere.

kfactor and grubbs_beck_critical_value have existing coverage under
tests/test_bulletin17c.py::TestUtilities (added alongside bulletin17c.py).
This file covers the PeakFQ-style exact-gamma LP3 functions, which had none.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from flowfreq.core import (
    YEAR_TYPES,
    assign_year_label,
    compute_ci_lp3,
    lp3_frequency_factor_peakfq,
    lp3_quantile_peakfq,
)


class TestAssignYearLabel:
    """Tests for the year-definition labeling shared by lowflow and regime.

    Used under two public names -- flowfreq.lowflow.LOW_FLOW_YEAR_TYPES is
    an alias of YEAR_TYPES defined here -- so this is the single place these
    boundary rules are pinned down.
    """

    def test_climatic_year_labeled_by_starting_calendar_year(self) -> None:
        """Apr 1 Y - Mar 31 Y+1 is climatic year Y (majority-of-months convention)."""
        dates = pd.DatetimeIndex(
            ["1999-04-01", "1999-12-31", "2000-01-01", "2000-03-31", "2000-04-01"]
        )
        labels = assign_year_label(dates, "climatic")
        assert list(labels) == [1999, 1999, 1999, 1999, 2000]

    def test_water_year_matches_usgs_download_peak_flow_convention(self) -> None:
        """Must agree with the inline water_year logic in usgs.download_peak_flow:
        Oct 1 (Y-1) - Sep 30 Y is water year Y."""
        dates = pd.DatetimeIndex(
            ["1998-10-01", "1998-12-31", "1999-01-01", "1999-09-30", "1999-10-01"]
        )
        labels = assign_year_label(dates, "water")
        assert list(labels) == [1999, 1999, 1999, 1999, 2000]

    def test_calendar_year_is_trivial(self) -> None:
        dates = pd.DatetimeIndex(["2020-01-01", "2020-12-31"])
        assert list(assign_year_label(dates, "calendar")) == [2020, 2020]

    def test_unknown_year_type_raises(self) -> None:
        with pytest.raises(ValueError, match="year_type must be one of"):
            assign_year_label(pd.DatetimeIndex(["2020-01-01"]), "fiscal")

    def test_year_types_tuple_contents(self) -> None:
        assert set(YEAR_TYPES) == {"climatic", "water", "calendar"}


def _exact_k(p: float, skew: float) -> float:
    """Reference implementation: reflect explicitly rather than sharing code
    with lp3_frequency_factor_peakfq, so this doesn't just check the function
    agrees with itself."""
    if abs(skew) < 1e-8:
        return stats.norm.ppf(p)
    alpha = 4.0 / skew**2
    if skew > 0:
        return (skew / 2.0) * (stats.gamma.ppf(p, alpha) - alpha)
    return (skew / 2.0) * (stats.gamma.ppf(1.0 - p, alpha) - alpha)


class TestLp3FrequencyFactorPeakfq:
    """Regression tests for the negative-skew sign inversion.

    Before the fix, lp3_frequency_factor_peakfq(p, skew) evaluated
    gamma.ppf(p, alpha) unconditionally. For skew < 0 the correct quantile
    needs gamma.ppf(1 - p, alpha) instead (the Pearson III lower tail maps to
    the gamma distribution's upper tail once skew flips sign) -- so for any
    negative-skew fit the function returned a K-factor that DECREASED as the
    non-exceedance probability increased, i.e. quantiles ran backwards.
    """

    @pytest.mark.parametrize("skew", [-1.2, -0.6, -0.302, -0.05, 0.0, 0.05, 0.302, 0.6, 1.2])
    def test_matches_exact_reflection(self, skew: float) -> None:
        """Matches an independently-written reflected gamma quantile, for both signs."""
        for p in (0.01, 0.02, 0.10, 0.30, 0.50, 0.70, 0.90, 0.98, 0.99):
            assert lp3_frequency_factor_peakfq(p, skew) == pytest.approx(
                _exact_k(p, skew), abs=1e-9
            )

    @pytest.mark.parametrize("skew", [-1.2, -0.6, -0.302, -0.05, 0.0, 0.05, 0.302, 0.6, 1.2])
    def test_monotonic_increasing_in_p(self, skew: float) -> None:
        """K must increase with non-exceedance probability regardless of skew sign.

        This is the direct symptom of the bug: a negative-skew K-factor that
        decreases with p silently inverts every quantile computed from it.
        """
        p = np.linspace(0.001, 0.999, 200)
        K = np.array([lp3_frequency_factor_peakfq(pi, skew) for pi in p])
        assert np.all(np.diff(K) > 0)

    def test_negative_skew_100yr_factor_is_positive(self) -> None:
        """Concrete symptom check at the national B17C regional skew (-0.302).

        A 100-year (non-exceedance p=0.99) K-factor must be positive -- it
        multiplies std_log to produce a flow ABOVE the mean. The pre-fix code
        returned a large negative K here.
        """
        K = lp3_frequency_factor_peakfq(0.99, -0.302)
        assert K > 0

    def test_zero_skew_matches_normal_quantile(self) -> None:
        """Zero skew is the shared boundary case between the two branches."""
        K = lp3_frequency_factor_peakfq(0.98, 0.0)
        assert K == pytest.approx(stats.norm.ppf(0.98))

    def test_positive_and_negative_skew_are_mirror_images(self) -> None:
        """K(p, -G) should equal -K(1-p, G): flipping skew sign and the tail
        together is a symmetry of the Pearson III family."""
        for p in (0.02, 0.10, 0.50, 0.90, 0.98):
            assert lp3_frequency_factor_peakfq(p, -0.6) == pytest.approx(
                -lp3_frequency_factor_peakfq(1 - p, 0.6), abs=1e-9
            )


class TestLp3QuantilePeakfq:
    """Tests for the quantile (flow-value) wrapper, both skew signs."""

    def test_increasing_with_return_period_positive_skew(self) -> None:
        vals = [lp3_quantile_peakfq(3.0, 0.25, 0.4, T) for T in (2, 10, 50, 100)]
        assert vals == sorted(vals)

    def test_increasing_with_return_period_negative_skew(self) -> None:
        """The case that was broken: negative skew must still increase with T."""
        vals = [lp3_quantile_peakfq(3.0, 0.25, -0.6, T) for T in (2, 10, 50, 100)]
        assert vals == sorted(vals)

    def test_median_flow_near_mean_for_symmetric_case(self) -> None:
        """T=2 (median) with zero skew should recover 10**mu."""
        q = lp3_quantile_peakfq(3.0, 0.25, 0.0, 2)
        assert q == pytest.approx(10**3.0, rel=1e-3)


class TestComputeCiLp3:
    """Tests for confidence-interval computation, both skew signs."""

    def test_lower_below_estimate_below_upper_negative_skew(self) -> None:
        """The CI must bracket the point estimate -- this fails silently
        under the pre-fix sign bug for negative skew."""
        estimate, lower, upper = compute_ci_lp3(3.0, 0.25, -0.6, n=40, T=100)
        assert lower < estimate < upper

    def test_lower_below_estimate_below_upper_positive_skew(self) -> None:
        estimate, lower, upper = compute_ci_lp3(3.0, 0.25, 0.6, n=40, T=100)
        assert lower < estimate < upper

    def test_wider_ci_for_longer_return_period(self) -> None:
        _, lower10, upper10 = compute_ci_lp3(3.0, 0.25, -0.302, n=40, T=10)
        _, lower100, upper100 = compute_ci_lp3(3.0, 0.25, -0.302, n=40, T=100)
        assert (upper100 - lower100) > (upper10 - lower10)
