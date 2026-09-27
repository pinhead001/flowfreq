"""The native EMA applies every perception-threshold period, per year, as ``siteQT`` does.

``vendor/peakfqr/R/readInputs.R::siteQT`` gives every water year the ``(tl, tu)``
of each threshold period covering it, in order, the later period winning; a year
with no peak inside a nonzero-threshold period is censored below that threshold.
``ExpectedMomentsAlgorithm`` used to honour only periods lying entirely before
the systematic record, and collapse several of those into one period at the
smallest threshold. These tests pin the per-year behaviour, and check the
native rows against ``flowfreq.fortran_engine.build_emafit_arrays`` -- the
library's own ``siteQT`` translation, which feeds ``emafitpr`` -- so they run
without the Fortran extension. ``tests/fortran_parity/
test_live_perception_thresholds.py`` fits the same cases through ``emafitpr``.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import pytest

from flowfreq import ExpectedMomentsAlgorithm
from flowfreq.bulletin17c import _PERCEPTION_QMIN
from flowfreq.core import EMAParameters
from flowfreq.fortran_engine import build_emafit_arrays

Row = Tuple[int, float, float, float, int]


def _peaks(years: np.ndarray, seed: int = 7) -> np.ndarray:
    """Deterministic LP3-ish peaks, all well above the thresholds used below."""
    rng = np.random.default_rng(seed)
    return 10 ** (3.5 + 0.25 * rng.standard_normal(len(years)))


def _native_rows(ema: ExpectedMomentsAlgorithm) -> List[Row]:
    """``(year, log ql, log qu, log tl, dtype)`` per native interval, siteQT's clamping."""
    rows = []
    for iv in ema._build_flow_intervals(0.0):
        rows.append(
            (
                int(iv.year),
                round(float(np.log10(max(iv.lower, _PERCEPTION_QMIN))), 10),
                round(float(np.log10(max(iv.upper, _PERCEPTION_QMIN))), 10),
                round(float(np.log10(max(iv.perception_threshold, _PERCEPTION_QMIN))), 10),
                int(iv.is_historical),
            )
        )
    return rows


def _sitqt_rows(
    flows: np.ndarray,
    years: np.ndarray,
    thresholds: Dict[Tuple[int, int], float],
    historical: Optional[List[Tuple[int, float]]] = None,
) -> List[Row]:
    arrays = build_emafit_arrays(
        flows, water_years=years, historical_peaks=historical, perception_thresholds=thresholds
    )
    assert np.all(arrays.tu == 20.0)
    return [
        (int(y), round(float(a), 10), round(float(b), 10), round(float(c), 10), int(d))
        for y, a, b, c, d in zip(arrays.years, arrays.ql, arrays.qu, arrays.tl, arrays.dtype)
    ]


def _by_year(ema: ExpectedMomentsAlgorithm):
    return {iv.year: iv for iv in ema._build_flow_intervals(0.0)}


class TestGapInsideTheRecord:
    """A threshold period in the middle of the systematic record, over a gap."""

    YEARS = np.array([y for y in range(1950, 1980) if not 1960 <= y <= 1964])
    THRESHOLDS = {(1950, 1979): 0.0, (1960, 1964): 5000.0}

    def _ema(self):
        return ExpectedMomentsAlgorithm(
            _peaks(self.YEARS), water_years=self.YEARS, perception_thresholds=self.THRESHOLDS
        )

    def test_gap_years_are_censored_below_the_threshold(self):
        by_year = _by_year(self._ema())
        assert len(by_year) == 30
        for year in range(1960, 1965):
            iv = by_year[year]
            assert iv.is_censored and not iv.is_historical
            assert (iv.lower, iv.upper, iv.perception_threshold) == (
                1e-20,
                5000.0,
                5000.0,
            )  # (Qmin, tl)

    def test_peaks_outside_the_period_are_unrestricted(self):
        by_year = _by_year(self._ema())
        assert by_year[1959].perception_threshold == 0.0
        assert not by_year[1959].is_censored

    def test_rows_match_sitqt(self):
        ema = self._ema()
        flows = _peaks(self.YEARS)
        assert sorted(_native_rows(ema)) == sorted(_sitqt_rows(flows, self.YEARS, self.THRESHOLDS))

    def test_gap_changes_the_fit(self):
        """Ignoring the gap (the old behaviour) is a different, smaller data set."""
        flows = _peaks(self.YEARS)
        with_gap = ExpectedMomentsAlgorithm(
            flows, water_years=self.YEARS, perception_thresholds=self.THRESHOLDS
        ).run_analysis()
        without = ExpectedMomentsAlgorithm(flows, water_years=self.YEARS).run_analysis()
        assert with_gap.n_peaks == 30
        assert without.n_peaks == 25
        # Five years known to be below 5000 cfs pull the mean down.
        assert with_gap.mean_log < without.mean_log

    def test_a_zero_threshold_gap_carries_no_information(self):
        ema = ExpectedMomentsAlgorithm(
            _peaks(self.YEARS), water_years=self.YEARS, perception_thresholds={(1950, 1979): 0.0}
        )
        assert set(_by_year(ema)) == set(int(y) for y in self.YEARS)


class TestTwoHistoricalPeriods:
    """Two pre-systematic periods with different thresholds stay separate."""

    YEARS = np.arange(1950, 1980)
    HISTORICAL = [(1890, 20000.0), (1925, 15000.0)]
    THRESHOLDS = {(1880, 1899): 18000.0, (1920, 1949): 12000.0}

    def _ema(self):
        return ExpectedMomentsAlgorithm(
            _peaks(self.YEARS),
            water_years=self.YEARS,
            historical_peaks=self.HISTORICAL,
            perception_thresholds=self.THRESHOLDS,
        )

    def test_each_period_censors_at_its_own_threshold(self):
        by_year = _by_year(self._ema())
        for year in range(1880, 1900):
            if year != 1890:
                assert by_year[year].upper == 18000.0
                assert by_year[year].perception_threshold == 18000.0
        for year in range(1920, 1950):
            if year != 1925:
                assert by_year[year].upper == 12000.0
                assert by_year[year].perception_threshold == 12000.0

    def test_years_between_the_periods_have_no_row(self):
        """No period covers 1900-1919, so nothing is known about them."""
        by_year = _by_year(self._ema())
        assert not any(1900 <= y <= 1919 for y in by_year)
        assert len(by_year) == 20 + 30 + 30

    def test_historical_peaks_carry_their_own_period_threshold(self):
        by_year = _by_year(self._ema())
        assert by_year[1890].is_historical and by_year[1890].perception_threshold == 18000.0
        assert by_year[1925].is_historical and by_year[1925].perception_threshold == 12000.0

    def test_rows_match_sitqt(self):
        flows = _peaks(self.YEARS)
        expected = _sitqt_rows(flows, self.YEARS, self.THRESHOLDS, self.HISTORICAL)
        assert sorted(_native_rows(self._ema())) == sorted(expected)

    def test_summary_parameters_are_unchanged(self):
        """EMAParameters still reports one historical period (the min-threshold envelope)."""
        params = self._ema()._ema_params
        assert (params.historical_start, params.historical_end) == (1880, 1949)
        assert params.historical_threshold == 12000.0


class TestOverlappingPeriodsLaterWins:
    """Where periods overlap, the one listed later sets the threshold."""

    YEARS = np.array([y for y in range(1950, 1980) if y not in (1956, 1961)])

    def _ema(self, thresholds):
        return ExpectedMomentsAlgorithm(
            _peaks(self.YEARS), water_years=self.YEARS, perception_thresholds=thresholds
        )

    def test_inner_period_listed_last_wins(self):
        thresholds = {(1950, 1979): 0.0, (1955, 1965): 8000.0, (1960, 1962): 3000.0}
        by_year = _by_year(self._ema(thresholds))
        assert by_year[1956].upper == 8000.0
        assert by_year[1961].upper == 3000.0
        assert by_year[1957].perception_threshold == 8000.0
        assert by_year[1960].perception_threshold == 3000.0
        assert by_year[1950].perception_threshold == 0.0
        flows = _peaks(self.YEARS)
        assert sorted(_native_rows(self._ema(thresholds))) == sorted(
            _sitqt_rows(flows, self.YEARS, thresholds)
        )

    def test_outer_period_listed_last_overrides_the_inner_one(self):
        thresholds = {(1950, 1979): 0.0, (1960, 1962): 3000.0, (1955, 1965): 8000.0}
        by_year = _by_year(self._ema(thresholds))
        assert by_year[1961].upper == 8000.0
        assert by_year[1960].perception_threshold == 8000.0

    def test_a_later_zero_period_clears_an_earlier_threshold(self):
        thresholds = {(1955, 1965): 8000.0, (1950, 1979): 0.0}
        by_year = _by_year(self._ema(thresholds))
        assert 1956 not in by_year and 1961 not in by_year
        assert by_year[1957].perception_threshold == 0.0


class TestLowOutlierThresholdOnPerceptionRows:
    """``gbtest`` applies the low-outlier cutoff to every row, not just systematic peaks."""

    YEARS = np.array([y for y in range(1950, 1980) if y != 1965])

    def test_gap_year_below_the_cutoff_is_recoded_to_it(self):
        ema = ExpectedMomentsAlgorithm(
            _peaks(self.YEARS),
            water_years=self.YEARS,
            perception_thresholds={(1965, 1965): 500.0},
            user_low_outlier_threshold=800.0,
        )
        ema.run_analysis()
        gap = next(iv for iv in ema.intervals if iv.year == 1965)
        assert (gap.lower, gap.upper) == (1e-6, 800.0)  # gbtest: (gbtmin, cutoff)
        # ... and its perception threshold is raised to the cutoff too.
        nobs, tl, _ = ema._perception_threshold_groups()
        assert list(tl) == [pytest.approx(np.log10(800.0))]
        assert list(nobs) == [30.0]

    def test_gap_year_threshold_below_every_peak_enters_mgbt(self):
        """``gbtest`` counts such a year as a systematic "less-than" observation."""
        flows = _peaks(self.YEARS)
        ema = ExpectedMomentsAlgorithm(
            flows, water_years=self.YEARS, perception_thresholds={(1965, 1965): 1.0}
        )
        threshold, n_low, pilf = ema._multiple_grubbs_beck()
        assert n_low >= 1
        assert pilf[0] == pytest.approx(1.0)
        assert threshold > 1.0

    def test_gap_year_threshold_above_a_peak_stays_out_of_mgbt(self):
        flows = _peaks(self.YEARS)
        ema = ExpectedMomentsAlgorithm(
            flows, water_years=self.YEARS, perception_thresholds={(1965, 1965): 5000.0}
        )
        plain = ExpectedMomentsAlgorithm(flows, water_years=self.YEARS)
        assert ema._multiple_grubbs_beck() == plain._multiple_grubbs_beck()


class TestBackwardCompatibility:
    """Callers that never passed perception thresholds see the same intervals."""

    YEARS = np.arange(1950, 1980)

    def test_historical_peaks_without_thresholds_use_the_derived_period(self):
        ema = ExpectedMomentsAlgorithm(
            _peaks(self.YEARS), water_years=self.YEARS, historical_peaks=[(1930, 25000.0)]
        )
        by_year = _by_year(ema)
        assert by_year[1930].perception_threshold == 25000.0
        for year in range(1931, 1950):
            assert (by_year[year].upper, by_year[year].perception_threshold) == (25000.0,) * 2
        assert by_year[1960].perception_threshold == 0.0

    def test_explicit_ema_params_period_still_censors(self):
        params = EMAParameters(
            systematic_start=1950,
            systematic_end=1979,
            historical_start=1940,
            historical_end=1949,
            historical_threshold=9000.0,
        )
        ema = ExpectedMomentsAlgorithm(
            _peaks(self.YEARS), water_years=self.YEARS, ema_params=params
        )
        by_year = _by_year(ema)
        assert all(by_year[y].upper == 9000.0 for y in range(1940, 1950))
        assert by_year[1950].perception_threshold == 0.0

    def test_record_gap_without_thresholds_is_left_out(self):
        years = np.array([y for y in self.YEARS if y != 1960])
        ema = ExpectedMomentsAlgorithm(_peaks(years), water_years=years)
        assert 1960 not in _by_year(ema)
        assert len(_by_year(ema)) == 29
