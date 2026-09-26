"""``kfactor`` is the exact Pearson III frequency factor, as peakfq computes it.

peakfq 8.1.0 computes its quantiles ``yp`` with ``emafit.f``'s ``qP3``
(line 3266), which inverts the incomplete gamma function and falls back to the
Wilson-Hilferty approximation only for ``|skew| < 0.001``. ``kfactor`` used
Wilson-Hilferty at every skew, so native quantiles drifted from peakfq's as
skew and return period grew, even when the fitted moments were identical:
1.35% at AEP 0.002 on USGS 06185500 (skew 0.87), whose moments match
peakfq's to 1e-13.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import stats
from scipy.special import ndtri

from flowfreq.core import MAX_ABS_SKEW, kfactor, kfactor_array, lp3_frequency_factor_peakfq
from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE

AEPS = (0.995, 0.9, 0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)


def _wilson_hilferty(skew: float, aep: float) -> float:
    """The approximation ``kfactor`` used before."""
    z = ndtri(1 - aep)
    k = skew / 6
    return (2 / skew) * ((1 + k * z - k * k) ** 3 - 1)


class TestKfactorIsExact:
    @pytest.mark.parametrize("skew", [-2.0, -1.0, -0.5, -0.1, 0.1, 0.5, 1.0, 2.0])
    @pytest.mark.parametrize("aep", AEPS)
    def test_matches_scipy_pearson3(self, skew, aep):
        """An independent exact reference: SciPy's standardized Pearson III."""
        expected = float(stats.pearson3.ppf(1 - aep, skew))
        assert kfactor(skew, aep) == pytest.approx(expected, rel=1e-10, abs=1e-12)

    @pytest.mark.parametrize("skew", [-1.5, -0.3, 0.3, 1.5])
    @pytest.mark.parametrize("aep", AEPS)
    def test_matches_the_exact_peakfq_frequency_factor(self, skew, aep):
        expected = lp3_frequency_factor_peakfq(1 - aep, skew)
        assert kfactor(skew, aep) == pytest.approx(expected, rel=1e-12, abs=1e-12)

    def test_high_skew_extreme_aep_is_not_wilson_hilferty(self):
        """The case that exposed the bug: the two differ by 0.06 in K here."""
        exact = kfactor(0.87, 0.002)
        assert exact == pytest.approx(float(stats.pearson3.ppf(0.998, 0.87)), rel=1e-10)
        assert abs(exact - _wilson_hilferty(0.87, 0.002)) > 0.02

    @pytest.mark.parametrize("aep", AEPS)
    def test_zero_skew_is_the_normal_quantile(self, aep):
        assert kfactor(0.0, aep) == pytest.approx(float(ndtri(1 - aep)), rel=1e-14, abs=1e-14)

    def test_continuous_through_the_wilson_hilferty_blend_band(self):
        """qP3 blends the two solutions over |skew| in [0.0007, 0.001]; no jump."""
        skews = np.linspace(-0.0015, 0.0015, 61)
        k = np.array([kfactor(float(g), 0.01) for g in skews])
        assert np.all(np.diff(k) > 0)
        assert np.max(np.abs(np.diff(k))) < 1e-4

    def test_skew_is_clipped(self):
        assert kfactor(-7.0, 0.01) == kfactor(-MAX_ABS_SKEW, 0.01)
        assert kfactor(7.0, 0.01) == kfactor(MAX_ABS_SKEW, 0.01)

    def test_array_form_agrees(self):
        aep = np.array(AEPS)
        assert np.array_equal(kfactor_array(0.87, aep), [kfactor(0.87, float(p)) for p in aep])

    def test_returns_a_python_float(self):
        assert type(kfactor(0.5, 0.01)) is float


# peakfq 8.1.0's yp for USGS 06185500.11 (Missouri River near Culbertson, MT,
# 1937-2022, station skew): a live emafitpr call on the vendored Fortran, via
# flowfreq.fortran_engine.run_fortran_reference with the station's peaks from
# wymt_ffa_2022A_EMPdata_7_4.csv. The record is uncensored and MGBT finds no
# low outliers, so the native moments equal these to 1e-13 and any quantile
# difference is the quantile function alone.
CULBERTSON_SITE = "06185500.11"
CULBERTSON_PEAKFQ_MOMENTS = (4.3182283986552585, 0.22427953912999424, 0.8705624473366247)
CULBERTSON_PEAKFQ_QUANTILES = {
    0.5: 19323.142463581287,
    0.2: 31005.28181893237,
    0.1: 41530.97121691836,
    0.04: 58787.45668942439,
    0.02: 75084.5348034185,
    0.01: 94873.85129649537,
    0.005: 118903.25617296723,
    0.002: 158705.93254407585,
}


@pytest.mark.requires_peakfqr_testdata
@pytest.mark.skipif(not TESTDATA_AVAILABLE, reason=SKIP_REASON)
class TestCulbertsonAgainstPeakfq:
    """The station the drift was reported on, against peakfq's own quantiles."""

    @pytest.fixture(scope="class")
    def b17c(self):
        from flowfreq.bulletin17c import Bulletin17C
        from tests.fixtures.wymt_peaks import load_site

        site = load_site(CULBERTSON_SITE)
        years = sorted(site.peaks)
        b17c = Bulletin17C(peak_flows=[site.peaks[y] for y in years], water_years=years)
        b17c.run_analysis(method="ema")
        return b17c

    def test_moments_match(self, b17c):
        r = b17c.results
        mean, std, skew = CULBERTSON_PEAKFQ_MOMENTS
        assert abs(r.mean_log - mean) < 1e-12
        assert abs(r.std_log - std) < 1e-12
        assert abs(r.skew_used - skew) < 1e-11

    def test_q500_matches(self, b17c):
        """AEP 0.002 was 1.35% high (160,843 cfs against peakfq's 158,706)."""
        q = b17c.compute_quantiles(aep=np.array([0.002]))
        assert float(q["flow_cfs"].iloc[0]) == pytest.approx(
            CULBERTSON_PEAKFQ_QUANTILES[0.002], rel=1e-7
        )

    def test_every_aep_down_to_0_002_matches(self, b17c):
        aeps = np.array(sorted(CULBERTSON_PEAKFQ_QUANTILES, reverse=True))
        q = b17c.compute_quantiles(aep=aeps)
        for aep, flow in zip(aeps, q["flow_cfs"]):
            assert float(flow) == pytest.approx(
                CULBERTSON_PEAKFQ_QUANTILES[float(aep)], rel=1e-7
            ), aep

    def test_confidence_table_centre_is_the_same_quantile(self, b17c):
        """The CI table's flow_cfs is the quantile the var_emab bounds surround."""
        aeps = np.array([0.01, 0.002])
        q = b17c.compute_quantiles(aep=aeps)
        ci = b17c.compute_confidence_limits(aep=aeps)
        assert np.allclose(ci["flow_cfs"], q["flow_cfs"], rtol=1e-15)
        assert np.all(ci["lower_5pct"] < ci["flow_cfs"])
        assert np.all(ci["flow_cfs"] < ci["upper_5pct"])
