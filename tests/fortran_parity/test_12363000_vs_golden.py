"""Native EMA against peakfq 8.1.0 on USGS 12363000.

What this site contributes that the other three cannot.

**A gap.** 98 peaks across 102 water years; 1924-1927 are unmeasured. Every
other parity case is contiguous, and ``wymt_case`` refuses a site that is not.
The gap years are deliberately *not* censored -- no perception threshold
asserts what would have been recorded there -- and the case sets
``fill_missing_years=False`` to say so. That is not a formality: filling them
moves the at-site skew from +0.435 to +0.250 and Q100 from 120,064 to 119,473,
so a case built the other way would be comparing two different analyses and
calling the difference a parity result.

**detrat reached with nothing censored.** ``emafit.f:763`` only reaches the
Halloween determinant ratio once the at-site skew clears 0.04. Big Sandy's is
0.0066, so it never gets there. Cains Coulee does, but with 11 PILFs, so its
censoring and its skew move together and neither is isolated. Here the at-site
skew is +0.435 -- driven by a 176,000 cfs peak 1.7x the next largest -- while
MGBT finds no low outliers at all. The skew path is exercised on a wholly
uncensored record.

Tolerances below are set from the measured agreement, not chosen in advance:
the moments match to 1e-8 or better and every quantile to under 0.11%, so a
tolerance an order of magnitude above that still fails on any real regression
while tolerating float64 noise.

Runs everywhere, including CI, because the reference comes from the committed
golden rather than the f2py extension.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.fortran_parity.conftest import aep_index


class TestRung1SameData:
    """Are the two codes fitting the same rows at all?"""

    def test_interval_count(self, golden_12363000, native_12363000):
        """98 observations, not the 102 years they span."""
        assert golden_12363000["outputs"]["n"] == 98
        assert native_12363000._results.n_peaks == 98

    def test_no_low_outliers_detected(self, golden_12363000, native_12363000):
        """MGBT censors nothing here, so the skew path is tested on its own."""
        assert golden_12363000["outputs"]["mgbt"]["gbnlow"] == 0
        assert native_12363000._results.n_low_outliers == 0


class TestRung3Moments:
    """Same mean, variance and skew on identical data?"""

    def test_mean_log(self, golden_12363000, native_12363000):
        expected = golden_12363000["outputs"]["cmoms"][0][0]
        assert abs(native_12363000._results.mean_log - expected) < 1e-8

    def test_std_log(self, golden_12363000, native_12363000):
        expected = golden_12363000["outputs"]["cmoms"][1][0] ** 0.5
        assert abs(native_12363000._results.std_log - expected) < 1e-8

    def test_at_site_skew(self, golden_12363000, native_12363000):
        """+0.435, well clear of the 0.04 floor that gates ``detrat``."""
        expected = golden_12363000["outputs"]["cmoms"][2][1]
        assert expected > 0.04
        assert abs(native_12363000._results.skew_station - expected) < 1e-8

    def test_weighted_skew(self, golden_12363000, native_12363000):
        expected = golden_12363000["outputs"]["cmoms"][2][0]
        assert abs(native_12363000._results.skew_weighted - expected) < 1e-6


class TestRung5Quantiles:
    """Does the frequency curve agree across the full range of AEPs?"""

    @pytest.mark.parametrize("aep", [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002])
    def test_quantile(self, golden_12363000, native_12363000, aep):
        i = aep_index(golden_12363000, aep)
        expected = 10 ** golden_12363000["outputs"]["quantiles"]["yp"][i]
        actual = native_12363000.compute_quantiles(aep=np.array([aep]))["flow_cfs"].iloc[0]
        assert abs(actual - expected) / expected * 100 < 1.0

    def test_hundred_year_peak(self, golden_12363000, native_12363000):
        """Pinned explicitly: this is the number the app reports and a user reads."""
        i = aep_index(golden_12363000, 0.01)
        expected = 10 ** golden_12363000["outputs"]["quantiles"]["yp"][i]
        actual = native_12363000.compute_quantiles(aep=np.array([0.01]))["flow_cfs"].iloc[0]
        assert expected == pytest.approx(120_064, rel=1e-4)
        assert abs(actual - expected) / expected * 100 < 0.1


class TestRung6ConfidenceBounds:
    """Cohn's asymmetric interval, the piece var_emab exists to get right."""

    @pytest.mark.parametrize("aep", [0.1, 0.02, 0.01, 0.002])
    def test_bounds(self, golden_12363000, native_12363000, aep):
        i = aep_index(golden_12363000, aep)
        lo_exp = 10 ** golden_12363000["outputs"]["quantiles"]["ci_low"][i]
        hi_exp = 10 ** golden_12363000["outputs"]["quantiles"]["ci_high"][i]
        ci = native_12363000.compute_confidence_limits(aep=np.array([aep]))
        assert abs(ci["lower_5pct"].iloc[0] - lo_exp) / lo_exp * 100 < 0.1
        assert abs(ci["upper_5pct"].iloc[0] - hi_exp) / hi_exp * 100 < 0.1
