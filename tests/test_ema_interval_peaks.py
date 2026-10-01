"""Interval peaks and upper perception thresholds in the native EMA and the Fortran builder.

The rules are peakfq 8.1.0's: rows as ``siteQT`` builds them
(``vendor/peakfqr/R/readInputs.R``), low-outlier censoring and MGBT inputs as
``gbtest`` decides them (``vendor/peakfqr/src/emafit.f`` lines 919-1075), and
the ``(nobs, tl, tu)`` groups ``compress2`` hands ``mse_ema``/``var_emab``.
Live-Fortran parity is in ``tests/fortran_parity/test_live_psf_convert.py``
(WY/MT 06328100) and ``tests/fortran_parity/test_interval_peaks_live.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from flowfreq.bulletin17c import _GBTMIN, Bulletin17C, ExpectedMomentsAlgorithm
from flowfreq.core import (
    FLOW_QMAX,
    FLOW_QMIN,
    FlowInterval,
    normalize_interval_peaks,
    perception_bounds,
)
from flowfreq.fortran_engine import build_emafit_arrays

YEARS = np.arange(1981, 2011)
FLOWS = np.array([float(200 + (37 * i) % 311) for i in range(len(YEARS))])


def _record(skip):
    keep = ~np.isin(YEARS, list(skip))
    return FLOWS[keep], YEARS[keep]


def _ema(skip=(), **kw):
    flows, years = _record(skip)
    return ExpectedMomentsAlgorithm(flows, water_years=years, **kw)


def _by_year(ema, low_threshold=0.0):
    return {iv.year: iv for iv in ema._build_flow_intervals(low_threshold)}


def _native_rows(ema):
    """``(year, ql, qu, tl, tu)`` rows in log10, the way psf_convert compares them."""
    return sorted(
        (
            iv.year,
            round(float(np.log10(max(iv.lower, FLOW_QMIN))), 10),
            round(float(np.log10(max(iv.upper, FLOW_QMIN))), 10),
            round(float(np.log10(max(iv.perception_threshold, FLOW_QMIN))), 10),
            round(float(np.log10(max(iv.perception_upper, FLOW_QMIN))), 10),
        )
        for iv in ema._build_flow_intervals(0.0)
    )


def _fortran_rows(arrays):
    return sorted(
        (int(y), round(float(a), 10), round(float(b), 10), round(float(c), 10), round(float(d), 10))
        for y, a, b, c, d in zip(arrays.years, arrays.ql, arrays.qu, arrays.tl, arrays.tu)
    )


class TestIntervalRows:
    def test_less_than_peak(self):
        ema = _ema(skip={1990}, interval_peaks=[(1990, 0.0, 150.0)])
        iv = _by_year(ema)[1990]
        assert (iv.lower, iv.upper) == (FLOW_QMIN, 150.0)
        assert iv.is_censored and iv.is_systematic
        assert (iv.perception_threshold, iv.perception_upper) == (0.0, FLOW_QMAX)
        results = ema.run_analysis()
        assert results.n_peaks == len(YEARS)
        assert results.n_censored == 1

    def test_greater_than_peak(self):
        ema = _ema(skip={1990}, interval_peaks=[(1990, 600.0, np.inf)])
        iv = _by_year(ema)[1990]
        assert (iv.lower, iv.upper) == (600.0, FLOW_QMAX)
        # A flood known to exceed 600 cfs pulls the fit up relative to one
        # at exactly 600.
        flows, years = _record({1990})
        exact = ExpectedMomentsAlgorithm(
            np.append(flows, 600.0), water_years=np.append(years, 1990)
        ).run_analysis()
        censored = ema.run_analysis()
        assert censored.mean_log > exact.mean_log

    def test_no_information_interval_has_no_row(self):
        ema = _ema(skip={1990}, interval_peaks=[(1990, 0.0, np.inf)])
        assert 1990 not in _by_year(ema)

    def test_interval_year_is_not_a_gap_year(self):
        """A perception threshold does not also censor an interval peak's year."""
        ema = _ema(
            skip={1990},
            interval_peaks=[(1990, 0.0, 150.0)],
            perception_thresholds={(1981, 2010): 0.0, (1990, 1990): 100.0},
        )
        intervals = ema._build_flow_intervals(0.0)
        assert [iv.year for iv in intervals].count(1990) == 1
        iv = _by_year(ema)[1990]
        assert (iv.lower, iv.upper, iv.perception_threshold) == (FLOW_QMIN, 150.0, 100.0)

    def test_interval_extends_the_record(self):
        """A trailing interval peak is still inside the systematic record."""
        ema = _ema(skip={2010}, interval_peaks=[(2010, 0.0, 150.0)])
        assert ema._ema_params.systematic_end == 2010
        assert ema._ema_params.historical_threshold is None

    @pytest.mark.parametrize(
        "interval_peaks, match",
        [
            ([(1990, 5.0, 5.0)], "upper > lower"),
            ([(1990, 10.0, 5.0)], "upper > lower"),
            ([(1990, -1.0, 5.0)], "invalid bound"),
            ([(1990, float("nan"), 5.0)], "invalid bound"),
            ([(1990, 0.0, 5.0), (1990, 1.0, 6.0)], "more than one"),
            ([(1990, 0.0)], "water_year, lower_cfs, upper_cfs"),
        ],
    )
    def test_invalid_interval_raises(self, interval_peaks, match):
        with pytest.raises(ValueError, match=match):
            _ema(skip={1990}, interval_peaks=interval_peaks)


class TestExactIntervalConflict:
    def test_systematic_conflict_raises(self):
        with pytest.raises(ValueError, match=r"\[1990\].*both an exact and an interval"):
            _ema(interval_peaks=[(1990, 0.0, 150.0)])

    def test_historical_conflict_raises(self):
        with pytest.raises(ValueError, match="both an exact and an interval"):
            _ema(historical_peaks=[(1950, 5000.0)], interval_peaks=[(1950, 0.0, 150.0)])

    def test_zero_flow_conflict_raises(self):
        flows = FLOWS.copy()
        flows[9] = 0.0
        with pytest.raises(ValueError, match="both an exact and an interval"):
            ExpectedMomentsAlgorithm(flows, water_years=YEARS, interval_peaks=[(1990, 0.0, 1.0)])

    def test_fortran_builder_conflict_raises(self):
        with pytest.raises(ValueError, match="both an exact and an interval"):
            build_emafit_arrays(FLOWS, water_years=YEARS, interval_peaks=[(1990, 0.0, 150.0)])

    def test_mom_refuses_interval_peaks(self):
        flows, years = _record({1990})
        b17c = Bulletin17C(flows, water_years=years, interval_peaks=[(1990, 0.0, 150.0)])
        with pytest.raises(ValueError, match="method='ema'"):
            b17c.run_analysis(method="mom")


class TestUpperThresholds:
    def test_float_threshold_is_unchanged(self):
        """A bare float still means (lower, Qmax), and fits exactly as a pair would."""
        assert perception_bounds(250.0) == (250.0, FLOW_QMAX)
        assert perception_bounds((250.0, np.inf)) == (250.0, FLOW_QMAX)
        flows, years = _record({1995, 1996})
        as_float = Bulletin17C(
            flows, water_years=years, perception_thresholds={(1981, 2010): 0.0, (1995, 1996): 250.0}
        ).run_analysis()
        as_pair = Bulletin17C(
            flows,
            water_years=years,
            perception_thresholds={(1981, 2010): 0.0, (1995, 1996): (250.0, np.inf)},
        ).run_analysis()
        assert (as_float.mean_log, as_float.std_log, as_float.skew_used) == (
            as_pair.mean_log,
            as_pair.std_log,
            as_pair.skew_used,
        )
        pd_eq = as_float.confidence_limits.equals(as_pair.confidence_limits)
        assert pd_eq

    def test_float_threshold_groups_keep_qmax(self):
        ema = _ema(skip={1995}, perception_thresholds={(1981, 2010): 0.0, (1995, 1995): 250.0})
        ema._build_flow_intervals(0.0)
        _, _, tu = ema._perception_threshold_groups()
        assert set(tu) == {20.0}

    def test_upper_threshold_changes_the_groups(self):
        thresholds = {(1981, 2010): 0.0, (1995, 1996): (100.0, 700.0)}
        ema = _ema(perception_thresholds=thresholds)
        ema._build_flow_intervals(0.0)
        nobs, tl, tu = ema._perception_threshold_groups()
        groups = {(round(a, 10), round(b, 10)): n for n, a, b in zip(nobs, tl, tu)}
        assert groups == {
            (-20.0, 20.0): len(YEARS) - 2,
            (2.0, round(float(np.log10(700.0)), 10)): 2,
        }
        for year in (1995, 1996):
            assert _by_year(ema)[year].perception_upper == 700.0

    def test_upper_threshold_carried_to_gap_years(self):
        ema = _ema(
            skip={1995}, perception_thresholds={(1981, 2010): 0.0, (1995, 1995): (100.0, 700.0)}
        )
        iv = _by_year(ema)[1995]
        assert (iv.lower, iv.upper) == (FLOW_QMIN, 100.0)
        assert (iv.perception_threshold, iv.perception_upper) == (100.0, 700.0)

    def test_upper_threshold_moves_the_confidence_bounds_only(self):
        """``tu`` enters the MSE/CI groups, not the moment fit (``p3est_ema`` reads
        only ``ql``/``qu``)."""
        flows, years = _record(set())
        plain = Bulletin17C(
            flows, water_years=years, perception_thresholds={(1995, 1996): 100.0}
        ).run_analysis()
        capped = Bulletin17C(
            flows, water_years=years, perception_thresholds={(1995, 1996): (100.0, 700.0)}
        ).run_analysis()
        assert capped.mean_log == plain.mean_log and capped.skew_used == plain.skew_used
        assert not np.allclose(
            capped.confidence_limits["upper_5pct"], plain.confidence_limits["upper_5pct"]
        )

    @pytest.mark.parametrize(
        "value, match",
        [((10.0, 5.0), "below the lower"), ((1.0, 2.0, 3.0), "pair"), ((np.nan, 5.0), "NaN")],
    )
    def test_invalid_threshold_raises(self, value, match):
        with pytest.raises(ValueError, match=match):
            _ema(perception_thresholds={(1981, 2010): value})

    def test_flow_interval_backward_compatible(self):
        iv = FlowInterval(lower=5.0, upper=5.0, year=2000)
        assert iv.perception_upper == FLOW_QMAX
        assert FlowInterval.from_censored(1.0, 2.0, 2000).perception_upper == FLOW_QMAX
        assert FlowInterval.from_historical(9.0, 1900, 3.0).perception_upper == FLOW_QMAX


class TestGbtestOnIntervalRows:
    def test_interval_below_cutoff_is_censored_to_gbtmin(self):
        ema = _ema(skip={1990, 1991}, interval_peaks=[(1990, 0.0, 150.0), (1991, 150.0, 400.0)])
        rows = _by_year(ema, low_threshold=180.0)
        assert (rows[1990].lower, rows[1990].upper) == (_GBTMIN, 180.0)
        # Its upper bound is above the cutoff: gbtest keeps it (`qu < gbcrit` only).
        assert (rows[1991].lower, rows[1991].upper) == (150.0, 400.0)

    def test_greater_than_is_never_censored(self):
        ema = _ema(skip={1990}, interval_peaks=[(1990, 10.0, np.inf)])
        assert _by_year(ema, low_threshold=1e6)[1990].upper == FLOW_QMAX

    def test_fixed_threshold_censors_interval(self):
        ema = _ema(
            skip={1990}, interval_peaks=[(1990, 0.0, 150.0)], user_low_outlier_threshold=180.0
        )
        ema.run_analysis()
        iv = next(i for i in ema.intervals if i.year == 1990)
        assert (iv.lower, iv.upper) == (_GBTMIN, 180.0)

    def test_less_than_below_every_peak_enters_mgbt(self):
        """``gbtest``: a censored dtype-0 row whose upper bound is at or below the
        smallest exact peak joins MGBT's sample at that bound."""
        ema = _ema(skip={1990, 1991}, interval_peaks=[(1990, 0.0, 150.0), (1991, 0.0, 250.0)])
        assert min(FLOWS[~np.isin(YEARS, [1990, 1991])]) == 200.0
        assert ema._less_than_uppers() == [150.0]

    def test_greater_than_never_enters_mgbt(self):
        ema = _ema(skip={1990}, interval_peaks=[(1990, 10.0, np.inf)])
        assert ema._less_than_uppers() == []

    def test_no_less_than_with_a_zero(self):
        """A zero is the smallest exact peak (Qmin); no censored row lies below it."""
        flows, years = _record({1990})
        flows = flows.copy()
        flows[0] = 0.0
        ema = ExpectedMomentsAlgorithm(flows, water_years=years, interval_peaks=[(1990, 0.0, 1.0)])
        assert ema._less_than_uppers() == []

    def test_less_than_can_be_the_mgbt_cutoff(self):
        """Several tiny less-than peaks are flagged together; the threshold is an
        exact value from the sample, the smallest non-outlier."""
        skip = {1983, 1985, 1987, 1989}
        ema = _ema(skip=skip, interval_peaks=[(y, 0.0, 1.0) for y in sorted(skip)])
        threshold, n_low, pilf = ema._multiple_grubbs_beck()
        assert n_low >= 4 and threshold >= 1.0
        assert pilf[:4] == pytest.approx([1.0] * 4)


class TestFortranBuilderAgrees:
    def test_rows_match_the_native_builder(self):
        skip = {1985, 1990, 2000}
        flows, years = _record(skip)
        kw = dict(
            historical_peaks=[(1950, 5000.0)],
            perception_thresholds={
                (1940, 1980): 3000.0,
                (1981, 2010): 0.0,
                (1985, 1985): 150.0,
                (1990, 1990): (50.0, 900.0),
            },
            interval_peaks=[(1990, 900.0, np.inf), (2000, 100.0, 300.0)],
        )
        ema = ExpectedMomentsAlgorithm(flows, water_years=years, **kw)
        arrays = build_emafit_arrays(flows, water_years=years, **kw)
        assert _native_rows(ema) == _fortran_rows(arrays)
        by_year = dict(zip(arrays.years, zip(arrays.ql, arrays.qu, arrays.tu, arrays.dtype)))
        assert by_year[1990] == (np.log10(900.0), 20.0, np.log10(900.0), 0)
        assert by_year[1985][:2] == (-20.0, np.log10(150.0))
        # 40 historical-period gap years (1940-1980 less 1950), 1985's gap row,
        # and the two interval peaks.
        assert arrays.n_censored == 43

    def test_normalize_interval_peaks(self):
        assert normalize_interval_peaks([(1990, 0, np.inf), (1991, 3, 4)]) == [
            (1990, FLOW_QMIN, FLOW_QMAX),
            (1991, 3.0, 4.0),
        ]
        assert normalize_interval_peaks(None) == []


class TestHistoricalIntervalPeaks:
    """Code 7 with 4 or 8: siteQT's censored row with dtype = 1."""

    KW = dict(
        historical_peaks=[(1950, 5000.0)],
        perception_thresholds={(1940, 1980): 3000.0, (1981, 2010): 0.0},
        historical_interval_peaks=[(1945, 3000.0, np.inf), (1960, 0.0, 120.0)],
    )

    def test_rows_are_historic_intervals(self):
        by_year = _by_year(_ema(**self.KW))
        greater = by_year[1945]
        assert greater.is_historical and greater.is_censored
        assert (greater.lower, greater.upper) == (3000.0, FLOW_QMAX)
        assert greater.perception_threshold == 3000.0
        less = by_year[1960]
        assert less.is_historical and (less.lower, less.upper) == (FLOW_QMIN, 120.0)

    def test_years_are_observed_not_gap_years(self):
        ema = _ema(**self.KW)
        assert 1945 not in ema._gap_year_thresholds()
        assert 1960 not in ema._gap_year_thresholds()

    def test_outside_mgbts_sample(self):
        """gbtest takes only dtype = 0 rows into MGBT, even below the smallest peak."""
        assert _ema(**self.KW)._less_than_uppers() == []

    def test_recoded_and_counted_below_the_cutoff(self):
        ema = _ema(**self.KW)
        by_year = {iv.year: iv for iv in ema._build_flow_intervals(150.0)}
        assert (by_year[1960].lower, by_year[1960].upper) == (_GBTMIN, 150.0)
        assert by_year[1960].is_historical
        assert (by_year[1945].lower, by_year[1945].upper) == (3000.0, FLOW_QMAX)
        assert ema._n_rows_below_cutoff == 1  # no systematic peak is below 150

    def test_no_information_interval_has_no_row(self):
        kw = dict(self.KW, historical_interval_peaks=[(1945, 0.0, np.inf)])
        assert 1945 not in _by_year(_ema(**kw))

    def test_fortran_builder_agrees_with_dtype_1(self):
        flows, years = _record(())
        ema = ExpectedMomentsAlgorithm(flows, water_years=years, **self.KW)
        arrays = build_emafit_arrays(flows, water_years=years, **self.KW)
        assert _native_rows(ema) == _fortran_rows(arrays)
        dtype = dict(zip(arrays.years, arrays.dtype))
        assert dtype[1945] == dtype[1960] == dtype[1950] == 1
        assert sum(arrays.dtype) == 3
        historic = {iv.year for iv in ema._build_flow_intervals(0.0) if iv.is_historical}
        assert historic == {1945, 1950, 1960}

    @pytest.mark.parametrize(
        "extra, match",
        [
            (dict(historical_peaks=[(1945, 4000.0)]), "historical interval"),
            (dict(interval_peaks=[(1945, 0.0, 10.0)]), "historical interval"),
        ],
    )
    def test_conflicts_raise(self, extra, match):
        kw = dict(self.KW, **extra)
        with pytest.raises(ValueError, match=match):
            _ema(**kw)
        flows, years = _record(())
        with pytest.raises(ValueError, match=match):
            build_emafit_arrays(flows, water_years=years, **kw)

    def test_systematic_year_conflict_raises(self):
        with pytest.raises(ValueError, match="historical interval"):
            _ema(historical_interval_peaks=[(1990, 0.0, 10.0)])

    def test_mom_refuses(self):
        flows, years = _record(())
        b = Bulletin17C(flows, years, historical_interval_peaks=[(1960, 0.0, 120.0)])
        with pytest.raises(ValueError, match="historical_interval_peaks"):
            b.run_analysis(method="mom")

    def test_ema_fits_and_counts(self):
        flows, years = _record(())
        r = Bulletin17C(flows, years, **self.KW).run_analysis(method="ema")
        assert r.n_historical == 3
        assert r.n_peaks == len(years) + 3 + 38  # and 1940-1980's 38 gap years
