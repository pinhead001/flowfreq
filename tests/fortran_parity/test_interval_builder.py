"""Tests for ``flowfreq.fortran_engine.build_emafit_arrays``.

Runs everywhere -- none of this needs the built f2py extension. It checks the
*translation* from a ``Bulletin17C`` input set to ``emafitpr``'s arrays, which
is pure array bookkeeping; whether the vendored Fortran agrees with the
resulting numbers is a separate, ``requires_fortran``-gated question answered
by ``test_live_vs_golden.py`` and (for the cases this module also builds)
would be answered the same way once the array construction itself is right.

Per docs/FORTRAN_ENGINE_DESIGN.md section 7's correctness plan:

* byte-equality against the existing ``build_emafit_inputs`` test-code builder
  on all four registered parity cases (``TestMatchesExistingParityCases``);
* the 12363000 gap-year switch, at the row-construction level -- the actual
  +0.435/+0.250 skew split needs a live Fortran call this environment cannot
  make, so this checks what *can* be checked here: that a declared threshold
  changes which rows exist and how they are censored, not just that the
  builder runs (``TestGapYearSwitch``);
* the declared-but-vacuous-threshold omission rule (``TestVacuousThreshold``);
* zero flows, overlapping perception-threshold periods, an omitted
  ``water_years``, and the ``gbthrsh0`` encoding.
"""

from __future__ import annotations

import numpy as np
import pytest

from flowfreq.core import EMAParameters
from flowfreq.fortran_engine import QMAX, QMIN, _gbthrsh0, build_emafit_arrays


def _sorted_rows(ql, qu, tl, tu, dtype):
    """(ql, qu, tl, tu, dtype) tuples sorted for an order-independent compare.

    Row order carries no meaning to emafitpr -- it is a sum over independent
    intervals -- and the two builders under comparison append rows in
    different block orders (systematic/historical/fill vs. one global
    year-sort), so "byte-equality" means the same multiset of rows, not the
    same iteration order.
    """
    return sorted(zip(ql, qu, tl, tu, dtype))


class TestMatchesExistingParityCases:
    """Byte-equality against ``tests/fortran_parity/cases.py::build_emafit_inputs``."""

    def test_big_sandy(self):
        from tests.fixtures.big_sandy import HISTORICAL_PEAKS, SYSTEMATIC_PEAKS, THRESHOLDS
        from tests.fortran_parity.cases import big_sandy_case, build_emafit_inputs

        expected = build_emafit_inputs(big_sandy_case())

        perception_thresholds = {(t["start"], t["end"]): t["lower"] for t in THRESHOLDS}
        arrays = build_emafit_arrays(
            peak_flows=list(SYSTEMATIC_PEAKS.values()),
            water_years=list(SYSTEMATIC_PEAKS.keys()),
            historical_peaks=list(HISTORICAL_PEAKS.items()),
            perception_thresholds=perception_thresholds,
        )

        assert arrays.n == len(expected["ql"])
        assert _sorted_rows(
            arrays.ql, arrays.qu, arrays.tl, arrays.tu, arrays.dtype
        ) == _sorted_rows(
            expected["ql"], expected["qu"], expected["tl"], expected["tu"], expected["dtype"]
        )
        # 37 censored historical gap years (TODO.md P3), nothing systematic.
        assert arrays.n_censored == 37
        assert arrays.n_zeros == 0

    @pytest.mark.parametrize(
        "case_name, site_no",
        [("powder_river_06326500", "06326500.00"), ("cains_coulee_06327450", "06327450.00")],
    )
    def test_wymt_sites(self, case_name, site_no):
        from tests.fixtures.wymt_peaks import load_site
        from tests.fortran_parity.cases import build_emafit_inputs, wymt_case

        try:
            site = load_site(site_no)
            case = wymt_case(site_no)
        except (FileNotFoundError, ValueError) as exc:
            pytest.skip(f"{case_name}: {exc}")

        expected = build_emafit_inputs(case)
        years = sorted(site.peaks)
        arrays = build_emafit_arrays(
            peak_flows=[site.peaks[y] for y in years],
            water_years=years,
        )

        assert arrays.n == len(expected["ql"])
        assert _sorted_rows(
            arrays.ql, arrays.qu, arrays.tl, arrays.tu, arrays.dtype
        ) == _sorted_rows(
            expected["ql"], expected["qu"], expected["tl"], expected["tu"], expected["dtype"]
        )
        assert arrays.n_censored == 0  # contiguous, no historic peaks

    def test_site_12363000_gap_years_omitted_by_default(self):
        """No perception threshold covers the gap -- matches the committed golden's 98 rows."""
        from tests.fixtures.site_12363000 import SYSTEMATIC_PEAKS
        from tests.fortran_parity.cases import build_emafit_inputs, site_12363000_case

        expected = build_emafit_inputs(site_12363000_case())
        years = sorted(SYSTEMATIC_PEAKS)
        arrays = build_emafit_arrays(
            peak_flows=[SYSTEMATIC_PEAKS[y] for y in years],
            water_years=years,
        )

        assert arrays.n == 98 == len(expected["ql"])
        assert _sorted_rows(
            arrays.ql, arrays.qu, arrays.tl, arrays.tu, arrays.dtype
        ) == _sorted_rows(
            expected["ql"], expected["qu"], expected["tl"], expected["tu"], expected["dtype"]
        )


class TestGapYearSwitch:
    """Design doc section 3's worked example, at the row-construction level.

    The published +0.435 (98 rows, gaps omitted) vs +0.250 (102 rows, gaps
    censored) skew split comes from feeding both row sets through the
    vendored Fortran, which this environment cannot build. What is checked
    here instead -- and is exactly what the builder is responsible for -- is
    that declaring a perception threshold over the gap changes the row set
    the *same way* the design doc's second construction describes: 102 rows,
    the four gap years present and censored at the declared lower bound, not
    98 with them silently dropped or filled with an invented exact value.

    The 5,000 cfs threshold below is a synthetic value chosen only to be
    comfortably below every peak in the record (the smallest systematic peak
    here is 19,700 cfs) -- it is not a claim about what USGS actually
    perceived at this gage in 1924-1927, unlike the 98-peak fixture itself.
    """

    YEARS = [1922, 1923, 1928, 1929, 1930]  # a slice, gap at 1924-1927
    FLOWS = [82_200, 88_000, 101_000, 69_700, 38_800]
    SYNTHETIC_THRESHOLD = 5_000.0

    def test_omitted_without_a_declared_threshold(self):
        arrays = build_emafit_arrays(peak_flows=self.FLOWS, water_years=self.YEARS)
        assert arrays.n == 5
        assert set(arrays.years) == set(self.YEARS)
        assert arrays.n_censored == 0

    def test_censored_with_a_declared_threshold(self):
        arrays = build_emafit_arrays(
            peak_flows=self.FLOWS,
            water_years=self.YEARS,
            perception_thresholds={(1922, 1930): self.SYNTHETIC_THRESHOLD},
        )
        assert arrays.n == 9  # 5 observed + 1924-1927
        assert set(arrays.years) == set(range(1922, 1931))
        assert arrays.n_censored == 4

        gap_rows = {
            year: (ql, qu, tl, tu)
            for year, ql, qu, tl, tu, dtype in zip(
                arrays.years, arrays.ql, arrays.qu, arrays.tl, arrays.tu, arrays.dtype
            )
            if year in (1924, 1925, 1926, 1927)
        }
        assert set(gap_rows) == {1924, 1925, 1926, 1927}
        log_threshold = np.log10(self.SYNTHETIC_THRESHOLD)
        for ql, qu, tl, tu in gap_rows.values():
            assert ql == pytest.approx(np.log10(QMIN))
            assert qu == pytest.approx(log_threshold)
            assert tl == pytest.approx(log_threshold)
            assert tu == pytest.approx(np.log10(QMAX))


class TestVacuousThreshold:
    """A declared-but-vacuous (<= Qmin) threshold is treated as no information.

    readInputs.R lines ~1030-1051: such a gap year gets the fully-open
    interval, then is dropped by the final ``keepNoInfo`` filter -- same as an
    undeclared gap year, not an exact peak at Qmin.
    """

    def test_vacuous_threshold_omits_the_gap_year(self):
        arrays = build_emafit_arrays(
            peak_flows=[1000.0, 1200.0, 900.0],
            water_years=[2000, 2002, 2003],  # 2001 missing
            perception_thresholds={(2000, 2003): 0.0},
        )
        assert arrays.n == 3
        assert set(arrays.years) == {2000, 2002, 2003}
        assert arrays.n_censored == 0

    def test_real_threshold_censors_the_same_gap_year(self):
        arrays = build_emafit_arrays(
            peak_flows=[1000.0, 1200.0, 900.0],
            water_years=[2000, 2002, 2003],
            perception_thresholds={(2000, 2003): 500.0},
        )
        assert arrays.n == 4
        assert 2001 in set(arrays.years)
        assert arrays.n_censored == 1

    def test_vacuous_threshold_does_not_break_a_real_peak_in_the_same_period(self):
        """A systematic peak inside a vacuous-threshold period keeps its
        exact value; only its (clamped-to-Qmin) perception threshold reflects
        the vacuous declaration, not log10(0)."""
        arrays = build_emafit_arrays(
            peak_flows=[1000.0, 1200.0, 900.0],
            water_years=[2000, 2001, 2002],
            perception_thresholds={(2000, 2002): 0.0},
        )
        assert arrays.n == 3
        idx = list(arrays.years).index(2000)
        assert arrays.ql[idx] == pytest.approx(np.log10(1000.0))
        assert arrays.tl[idx] == pytest.approx(np.log10(QMIN))


class TestZeroFlows:
    """Open question 1: a zero flow is an exact peak at Qmin, not an omission."""

    def test_zero_flow_is_an_exact_peak_at_qmin(self):
        arrays = build_emafit_arrays(peak_flows=[0.0, 500.0, 800.0], water_years=[2000, 2001, 2002])
        assert arrays.n == 3
        assert arrays.n_zeros == 1
        zero_idx = list(arrays.years).index(2000)
        assert arrays.ql[zero_idx] == arrays.qu[zero_idx] == pytest.approx(np.log10(QMIN))
        # Reported as its real value (0.0), not the internal Qmin clamp, so
        # the adapter can later tell a real zero apart from a tiny nonzero peak.
        assert arrays.systematic_peaks[2000] == 0.0

    def test_negative_flow_raises(self):
        with pytest.raises(ValueError):
            build_emafit_arrays(peak_flows=[-1.0], water_years=[2000], **{})


class TestOverlappingPerceptionPeriods:
    """readInputs.R: 'the last one specified is given priority' for overlapping years."""

    def test_later_period_wins_on_overlap(self):
        arrays = build_emafit_arrays(
            peak_flows=[1000.0] * 16,
            water_years=list(range(2000, 2016)),
            perception_thresholds={
                (2000, 2010): 1_000.0,
                (2005, 2015): 2_000.0,  # declared later -> wins 2005-2010
            },
        )
        by_year = dict(zip(arrays.years, arrays.tl))
        assert by_year[2000] == pytest.approx(np.log10(1_000.0))
        assert by_year[2005] == pytest.approx(np.log10(2_000.0))
        assert by_year[2010] == pytest.approx(np.log10(2_000.0))
        assert by_year[2015] == pytest.approx(np.log10(2_000.0))

    def test_earlier_period_wins_when_declared_last(self):
        """Same periods, reversed insertion order -> the outcome reverses too."""
        arrays = build_emafit_arrays(
            peak_flows=[1000.0] * 16,
            water_years=list(range(2000, 2016)),
            perception_thresholds={
                (2005, 2015): 2_000.0,
                (2000, 2010): 1_000.0,  # declared later -> wins 2005-2010 now
            },
        )
        by_year = dict(zip(arrays.years, arrays.tl))
        assert by_year[2005] == pytest.approx(np.log10(1_000.0))
        assert by_year[2010] == pytest.approx(np.log10(1_000.0))
        assert by_year[2015] == pytest.approx(np.log10(2_000.0))


class TestHistoricalThresholdFromEmaParams:
    """A caller relying on Bulletin17C's auto-derived historical period, not an
    explicit ``perception_thresholds`` entry, still gets the gap censored."""

    def test_ema_params_historical_period_is_honored(self):
        ema_params = EMAParameters(
            systematic_start=2000,
            systematic_end=2005,
            historical_start=1990,
            historical_end=1999,
            historical_threshold=3_000.0,
        )
        arrays = build_emafit_arrays(
            peak_flows=[1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0],
            water_years=list(range(2000, 2006)),
            historical_peaks=[(1995, 9_000.0)],
            ema_params=ema_params,
        )
        # 1990-1999 minus the one historical peak = 9 censored gap years.
        assert arrays.n_censored == 9
        assert arrays.n == 6 + 1 + 9

    def test_explicit_perception_thresholds_override_ema_params_historical(self):
        ema_params = EMAParameters(
            systematic_start=2000,
            systematic_end=2002,
            historical_start=1990,
            historical_end=1999,
            historical_threshold=3_000.0,
        )
        arrays = build_emafit_arrays(
            peak_flows=[1000.0, 1100.0, 1200.0],
            water_years=[2000, 2001, 2002],
            ema_params=ema_params,
            perception_thresholds={(1990, 1999): 0.0},  # explicitly vacuous
        )
        # The explicit (vacuous) declaration wins over ema_params' historical
        # threshold, so none of 1990-1999 appears at all.
        assert arrays.n == 3
        assert set(arrays.years) == {2000, 2001, 2002}


class TestWaterYearsOmitted:
    def test_defaults_to_recent_years(self):
        from datetime import datetime

        n = 5
        arrays = build_emafit_arrays(peak_flows=[100.0, 200.0, 300.0, 400.0, 500.0])
        assert arrays.n == n
        end_year = datetime.now().year
        assert set(arrays.years) == set(range(end_year - n + 1, end_year + 1))


class TestGbthrsh0Encoding:
    """fortranWrappers.R's LOthresh -> gbthrsh0 mapping (~lines 93-127)."""

    def test_no_override_runs_mgbt(self):
        assert _gbthrsh0(None) == -99.0

    def test_real_threshold_is_fixed_and_logged(self):
        assert _gbthrsh0(500.0) == pytest.approx(np.log10(500.0))

    def test_zero_or_negative_also_runs_mgbt(self):
        assert _gbthrsh0(0.0) == -99.0
        assert _gbthrsh0(-1.0) == -99.0


class TestMinimumRows:
    def test_fewer_than_three_rows_raises(self):
        with pytest.raises(ValueError):
            build_emafit_arrays(peak_flows=[100.0, 200.0], water_years=[2000, 2001])
