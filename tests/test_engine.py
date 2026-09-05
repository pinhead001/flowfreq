"""Tests for flowfreq.engine.B17CEngine.

A simplified LP3 front end: fit moments, get quantiles and confidence
intervals, print a summary. It is not the Bulletin17C path -- it does no MGBT
screening, no EMA, no regional skew weighting -- and the tests below pin where
the two deliberately differ from where they differ by accident.

One difference was not deliberate: the station skew estimator, biased before
0.4.0 and now unbiased. See TestSkewEstimatorMatchesBulletin17C.
"""

from __future__ import annotations

import numpy as np
import pytest

from flowfreq.bulletin17c import Bulletin17C
from flowfreq.core import PeakRecord
from flowfreq.engine import STANDARD_RETURN_PERIODS, B17CEngine
from tests.fixtures.big_sandy import SYSTEMATIC_PEAKS


@pytest.fixture(scope="module")
def peaks():
    years = sorted(SYSTEMATIC_PEAKS)
    return np.array(years), np.array([SYSTEMATIC_PEAKS[y] for y in years], dtype=float)


@pytest.fixture
def fitted(peaks):
    years, flows = peaks
    engine = B17CEngine()
    engine.fit_from_flows(flows, years)
    return engine


class TestUnfittedState:
    """Every accessor has to behave before fit(), not just after."""

    def test_parameters_are_none(self):
        engine = B17CEngine()
        assert engine.mu is None
        assert engine.sigma is None
        assert engine.skew is None

    def test_n_is_zero(self):
        assert B17CEngine().n == 0

    def test_quantiles_refuses(self):
        with pytest.raises(RuntimeError, match="fit"):
            B17CEngine().quantiles()

    def test_quantiles_with_ci_refuses(self):
        with pytest.raises(RuntimeError, match="fit"):
            B17CEngine().quantiles_with_ci()

    def test_summary_says_so_rather_than_raising(self):
        """summary() is for humans; it reports the state instead of raising."""
        assert "not fitted" in B17CEngine().summary().lower()


class TestFit:
    def test_returns_and_stores_the_same_triple(self, peaks):
        years, flows = peaks
        engine = B17CEngine()
        returned = engine.fit_from_flows(flows, years)
        assert returned == engine.params
        assert returned == (engine.mu, engine.sigma, engine.skew)

    def test_moments_match_numpy_on_log10(self, fitted, peaks):
        _, flows = peaks
        lg = np.log10(flows)
        assert fitted.mu == pytest.approx(lg.mean())
        assert fitted.sigma == pytest.approx(lg.std(ddof=1))

    def test_mean_and_sigma_agree_with_bulletin17c(self, fitted, peaks):
        """The first two moments are the same estimator in both paths."""
        years, flows = peaks
        b = Bulletin17C(peak_flows=flows, water_years=years)
        b.run_analysis(method="mom")
        assert fitted.mu == pytest.approx(b.results.mean_log, rel=1e-9)
        assert fitted.sigma == pytest.approx(b.results.std_log, rel=1e-9)

    def test_n_counts_observations(self, fitted, peaks):
        _, flows = peaks
        assert fitted.n == len(flows)

    def test_fewer_than_three_observations_refuses(self):
        with pytest.raises(ValueError, match="[Aa]t least 3"):
            B17CEngine().fit([PeakRecord(year=2000, flow=100.0)])

    def test_none_flows_do_not_count_toward_the_minimum(self):
        """Three records but two null flows is one observation, not three."""
        records = [
            PeakRecord(year=2000, flow=100.0),
            PeakRecord(year=2001, flow=None),
            PeakRecord(year=2002, flow=None),
        ]
        with pytest.raises(ValueError, match="[Aa]t least 3"):
            B17CEngine().fit(records)

    def test_fit_from_flows_drops_nan(self):
        flows = [100.0, 200.0, float("nan"), 300.0, 400.0]
        engine = B17CEngine()
        engine.fit_from_flows(flows)
        assert engine.n == 4

    def test_fit_from_flows_synthesises_years_when_omitted(self):
        engine = B17CEngine()
        engine.fit_from_flows([100.0, 200.0, 300.0])
        assert engine.n == 3

    def test_refitting_replaces_rather_than_accumulates(self, peaks):
        """State from a previous fit must not leak into the next one."""
        years, flows = peaks
        engine = B17CEngine()
        engine.fit_from_flows(flows, years)
        engine.fit_from_flows(flows[:20], years[:20])
        assert engine.n == 20


class TestQuantiles:
    def test_default_return_periods(self, fitted):
        assert set(fitted.quantiles()) == set(STANDARD_RETURN_PERIODS)

    def test_monotonic_in_return_period(self, fitted):
        q = fitted.quantiles([2, 10, 50, 100, 500])
        values = [q[T] for T in (2, 10, 50, 100, 500)]
        assert values == sorted(values)

    def test_all_positive(self, fitted):
        assert all(v > 0 for v in fitted.quantiles().values())

    def test_custom_return_periods_are_honoured(self, fitted):
        assert set(fitted.quantiles([7, 13])) == {7, 13}


class TestQuantilesWithCI:
    def test_shape(self, fitted):
        ci = fitted.quantiles_with_ci([100])
        assert set(ci[100]) == {"estimate", "lower", "upper"}

    def test_bounds_bracket_the_estimate(self, fitted):
        for T, row in fitted.quantiles_with_ci([2, 10, 100]).items():
            assert row["lower"] < row["estimate"] < row["upper"], f"at T={T}"

    def test_estimate_matches_the_plain_quantile(self, fitted):
        """The two entry points must not disagree about the central value."""
        plain = fitted.quantiles([100])[100]
        with_ci = fitted.quantiles_with_ci([100])[100]["estimate"]
        assert with_ci == pytest.approx(plain, rel=1e-9)

    def test_a_wider_alpha_gives_a_narrower_interval(self, fitted):
        narrow = fitted.quantiles_with_ci([100], alpha=0.05)[100]
        wide = fitted.quantiles_with_ci([100], alpha=0.32)[100]
        assert (wide["upper"] - wide["lower"]) < (narrow["upper"] - narrow["lower"])


class TestFrequencyTable:
    def test_columns_and_length(self, fitted):
        df = fitted.frequency_table([2, 10, 100])
        assert len(df) == 3
        assert list(df.columns) == [
            "Return Period (yr)",
            "AEP (%)",
            "Flow (cfs)",
            "Lower 95%",
            "Upper 95%",
        ]

    def test_aep_is_the_reciprocal_as_a_percentage(self, fitted):
        df = fitted.frequency_table([100])
        assert df["AEP (%)"].iloc[0] == pytest.approx(1.0)

    def test_alpha_renames_the_bound_columns(self, fitted):
        df = fitted.frequency_table([100], alpha=0.10)
        assert "Lower 90%" in df.columns and "Upper 90%" in df.columns

    def test_values_agree_with_quantiles_with_ci(self, fitted):
        df = fitted.frequency_table([100])
        ci = fitted.quantiles_with_ci([100])[100]
        assert df["Flow (cfs)"].iloc[0] == pytest.approx(ci["estimate"])
        assert df["Lower 95%"].iloc[0] == pytest.approx(ci["lower"])


class TestSummary:
    def test_reports_the_fitted_moments(self, fitted):
        text = fitted.summary()
        assert f"{fitted.mu:.4f}" in text
        assert f"{fitted.sigma:.4f}" in text
        assert f"{fitted.skew:.4f}" in text

    def test_reports_sample_size(self, fitted):
        assert str(fitted.n) in fitted.summary()

    def test_quotes_quantiles_consistent_with_quantiles(self, fitted):
        """A summary that disagrees with the API it summarises is worse than none."""
        text = fitted.summary()
        assert f"{fitted.quantiles([100])[100]:,.0f}" in text


class TestSkewEstimatorMatchesBulletin17C:
    """engine.fit uses Bulletin 17C Eq. 7-2, the same estimator as Bulletin17C.

    Before 0.4.0 ``B17CEngine.fit`` computed ``((x - mean)**3).mean() /
    std**3`` -- the population (biased) coefficient of skewness -- while
    Bulletin 17C Eq. 7-2, and ``Bulletin17C.run_analysis``, use the unbiased
    sample estimator ``n * sum((x - mean)**3) / ((n-1)(n-2) * std**3)``.

    The two differ by a factor of ``n**2 / ((n-1)(n-2))``: 7.2% at n=44, and
    39% at n=10. Short records are ordinary in flood frequency work, and on
    Big Sandy the biased form put Q100 0.57% high and Q500 0.93% high.

    These tests are the guard against a regression to the biased form. The
    first pins agreement with the Bulletin 17C path; the second names the old
    estimator explicitly, so a silent revert fails on a specific, recognisable
    number rather than on a tolerance.
    """

    def test_station_skew_matches_bulletin17c(self, fitted, peaks):
        years, flows = peaks
        b = Bulletin17C(peak_flows=flows, water_years=years)
        b.run_analysis(method="mom")
        assert fitted.skew == pytest.approx(b.results.skew_station, rel=1e-6)

    def test_is_not_the_biased_population_coefficient(self, fitted, peaks):
        """The pre-0.4.0 form, pinned by name so a revert is recognisable."""
        _, flows = peaks
        log_flows = np.log10(flows)
        biased = float(((log_flows - log_flows.mean()) ** 3).mean() / log_flows.std(ddof=1) ** 3)
        n = fitted.n
        assert fitted.skew / biased == pytest.approx(n**2 / ((n - 1) * (n - 2)), rel=1e-9)

    def test_agrees_with_bulletin17c_on_a_short_record(self, peaks):
        """n=12, where the old bias factor was 1.31 rather than 1.07."""
        years, flows = peaks
        engine = B17CEngine()
        engine.fit_from_flows(flows[:12], years[:12])
        b = Bulletin17C(peak_flows=flows[:12], water_years=years[:12])
        b.run_analysis(method="mom")
        assert engine.skew == pytest.approx(b.results.skew_station, rel=1e-6)
