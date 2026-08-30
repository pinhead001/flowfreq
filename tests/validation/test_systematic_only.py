"""
Systematic-Only Regression Test

Validates that EMA correctly reduces to standard Method-of-Moments LP3
fitting when there is no historical record, no perception thresholds,
and no censored/low-outlier intervals. This is a pure internal-consistency
check (not a PeakfqSA comparison) that guards against regressions in the
EMA iteration and skew-weighting logic introduced while fixing the
Big Sandy historical-data benchmark.

Reference values are computed independently with numpy/scipy, not sourced
from PeakfqSA, since this scenario has no censoring for EMA to act on.

Data: Big Sandy systematic record only (1930-1973), no historical peaks.
"""

import numpy as np
import pytest

from hydrolib.bulletin17c import Bulletin17C
from hydrolib.core import kfactor_array
from tests.peakfqsa.fixtures.big_sandy import SYSTEMATIC_PEAKS

PEAK_FLOWS = list(SYSTEMATIC_PEAKS.values())
WATER_YEARS = list(SYSTEMATIC_PEAKS.keys())


def _reference_moments():
    """Independent MOM reference computation (numpy, not hydrolib)."""
    log_flows = np.log10(PEAK_FLOWS)
    n = len(log_flows)
    mean_log = np.mean(log_flows)
    std_log = np.std(log_flows, ddof=1)
    skew = n * np.sum((log_flows - mean_log) ** 3) / ((n - 1) * (n - 2) * std_log**3)
    return n, mean_log, std_log, skew


class TestSystematicOnlyEquivalence:
    """EMA with no historical/censored data must exactly match MOM."""

    @pytest.fixture
    def ema_result(self):
        b17c = Bulletin17C(peak_flows=PEAK_FLOWS, water_years=WATER_YEARS)
        b17c.run_analysis(method="ema")
        return b17c

    @pytest.fixture
    def mom_result(self):
        b17c = Bulletin17C(peak_flows=PEAK_FLOWS, water_years=WATER_YEARS)
        b17c.run_analysis(method="mom")
        return b17c

    def test_ema_converges_immediately(self, ema_result):
        """With no censored intervals, EMA should converge in a single iteration."""
        r = ema_result.results
        assert r.ema_converged
        assert r.ema_iterations == 1
        assert r.n_censored == 0
        assert r.n_historical == 0

    def test_ema_matches_reference_moments(self, ema_result):
        n, mean_ref, std_ref, skew_ref = _reference_moments()
        r = ema_result.results

        assert r.n_systematic == n
        assert r.mean_log == pytest.approx(mean_ref, abs=1e-9)
        assert r.std_log == pytest.approx(std_ref, abs=1e-9)
        assert r.skew_station == pytest.approx(skew_ref, abs=1e-9)

    def test_ema_matches_mom(self, ema_result, mom_result):
        """EMA and MOM must produce identical moments for this scenario."""
        re = ema_result.results
        rm = mom_result.results

        assert re.mean_log == pytest.approx(rm.mean_log, abs=1e-9)
        assert re.std_log == pytest.approx(rm.std_log, abs=1e-9)
        assert re.skew_station == pytest.approx(rm.skew_station, abs=1e-9)

    @pytest.mark.parametrize("aep", [0.5, 0.2, 0.1, 0.02, 0.01, 0.002])
    def test_quantiles_match_reference(self, ema_result, aep):
        n, mean_ref, std_ref, skew_ref = _reference_moments()
        K = kfactor_array(skew_ref, np.array([aep]))
        expected_flow = 10 ** (mean_ref + K[0] * std_ref)

        qdf = ema_result.compute_quantiles(np.array([aep]))
        actual_flow = qdf["flow_cfs"].iloc[0]

        assert actual_flow == pytest.approx(expected_flow, rel=1e-6)
