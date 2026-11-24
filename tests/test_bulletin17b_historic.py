"""Bulletin17B's outlier order and the Appendix 6 historic adjustment.

Printed values are from Appendix 12 Example 2 and the Appendix 6 Big Sandy example.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from flowfreq.bulletin17b import (
    Bulletin17B,
    Bulletin17BError,
    _outlier_order,
    confidence_deviates,
    expected_probability,
    historic_moments,
    log_moments,
    outlier_k,
    station_skew_mse,
    weighted_skew,
)
from tests.fixtures.b17b import load, record

BOTH = "both on unadjusted statistics"
# 18 ordinary years, used to build records with outliers.
BASE = [1000.0, 1300, 800, 1150, 900, 1250, 700, 1050, 1400]
BASE = BASE + [q * 1.03 for q in BASE]


def _p(text):
    return float(str(text).replace(",", ""))


def example2(**kw):
    r = record("example2")
    kw.setdefault("generalized_skew", -0.3)
    return Bulletin17B(r.peaks_cfs, r.water_years, **kw).run_analysis()


def big_sandy(flows=None, **kw):
    bs = record("appendix6_big_sandy")
    kw.setdefault("historic_peaks", bs.historic_peaks)
    kw.setdefault("historic_period", bs.historic_period)
    return Bulletin17B(flows or bs.peaks_cfs, bs.water_years, **kw).run_analysis()


def synthetic(extra, **kw):
    """The 18 ordinary years and then ``extra``, from 1950."""
    flows = BASE + list(extra)
    return Bulletin17B(flows, list(range(1950, 1950 + len(flows))), **kw)


# --- Example 2 ---------------------------------------------------------------------------


def test_example2_high_outlier_becomes_historic():
    res = example2(historic_period=(1892, 1973), hand_computation=True)
    assert res.outlier_order == BOTH
    assert res.high_outliers == [(1953, 71500.0)] and res.low_outliers == []
    h = res.historic
    assert (h.h, h.z, h.l) == (82, 1, 0) and h.w == pytest.approx(2.13158, abs=1e-5)
    assert (h.mean, h.std, h.skew) == (
        pytest.approx(3.5375, abs=1e-4),
        pytest.approx(0.4377, abs=1e-4),
        # Unrounded, 0.16535. The printed 0.1650 comes from M - M~ = 3.5212 - 3.5375 =
        # -0.0163 (exactly -0.01626): see test_example2_printed_skew_arithmetic.
        pytest.approx(0.1650, abs=4e-4),
    )
    assert h.historic_peaks == [(1953, 71500.0)]
    assert res.n_years == 39 and res.n_systematic == 38
    assert res.skew_station_mse == pytest.approx(0.073, abs=5e-4)
    # 0.07473: the printed 0.0745 carries the 0.1650.
    assert res.skew_weighted == pytest.approx(0.0745, abs=3e-4) and res.skew_used == 0.1
    got = res.quantiles.set_index("aep").flow_cfs
    rows = load("example2")["final_frequency_table_12_7"]["rows"]
    assert len(rows) == 9
    for row in rows:
        assert got[_p(row["P"])] == pytest.approx(_p(row["Q_cfs"]), rel=0.006)


def test_example2_printed_skew_arithmetic():
    """Eq. 12-22 at the printed 4-decimal figures gives the printed 0.1650."""
    res = example2(historic_period=(1892, 1973))

    def eq_6_4b(m, s, g, w, xz, m_adj, s_adj):
        n, h = 38, 82
        bracket = (
            w * (n - 1) * (n - 2) * s**3 * g / n
            + 3 * w * (n - 1) * (m - m_adj) * s**2
            + w * n * (m - m_adj) ** 3
            + (xz - m_adj) ** 3
        )
        return h * bracket / ((h - 1) * (h - 2) * s_adj**3)

    printed = eq_6_4b(3.5212, 0.4177, -0.0949, 2.13158, 4.8543, 3.5375, 0.4377)
    assert printed == pytest.approx(0.1650, abs=5e-5)
    r = record("example2")
    x = np.log10([q for q in r.peaks_cfs if q != 71500.0])
    h = res.historic
    exact = eq_6_4b(*log_moments(x), h.w, math.log10(71500.0), h.mean, h.std)
    assert exact == pytest.approx(h.skew, abs=1e-10)


def test_example2_thresholds_are_eqs_7_and_8a():
    res = example2(historic_period=(1892, 1973))
    # Eqs. 12-17 and 12-18: 3.5553 +/- 2.671 (.4642), printed as 62400 and 207 cfs.
    assert math.log10(res.high_outlier_threshold) == pytest.approx(4.7952, abs=2e-4)
    assert math.log10(res.low_outlier_threshold) == pytest.approx(2.3154, abs=2e-4)
    assert res.high_outlier_threshold == pytest.approx(62400, rel=0.002)
    assert res.low_outlier_threshold == pytest.approx(207, rel=0.003)


def test_example2_without_historic_period_keeps_the_outlier():
    res = example2()
    assert res.historic is None and res.high_outliers == [(1953, 71500.0)]
    assert res.n_systematic == 39 and any("kept" in s for s in res.steps)


def test_example2_plotting_positions_table_12_6():
    pp = example2(historic_period=(1892, 1973)).plotting_positions.set_index("year")
    assert len(pp) == 39
    for row in load("example2")["computations"]["step5_plotting_positions"]["table_12_6"]:
        got = pp.loc[int(row["year"])]
        assert got.flow_cfs == _p(row["Q_cfs"])
        # The bulletin carries W as 2.1316, so the order drifts by 2e-5 a rank.
        assert got.weighted_order == pytest.approx(_p(row["m_tilde"]), abs=2.5e-4)
        assert got.exceedance == pytest.approx(
            _p(row["PP_exceedance_probability"]), abs=6e-5
        )  # printed to 4 decimals


# --- Appendix 6, Big Sandy ---------------------------------------------------------------


def test_big_sandy_appendix6_exact():
    fig = load("appendix6_big_sandy")["figure_6_1_computations"]
    rows = fig["eq_6_5_frequency_table"]["rows"]
    assert len(rows) == 12
    res = big_sandy(generalized_skew=-0.2, aeps=[_p(r["exceedance_percent"]) / 100 for r in rows])
    assert res.outlier_order == BOTH and not res.high_outliers and not res.low_outliers
    h = res.historic
    assert (h.h, h.z, h.l) == (77, 3, 0) and h.w == pytest.approx(1.68182, abs=1e-5)
    assert res.mean_log == h.mean == pytest.approx(3.71581, abs=1e-5)
    assert res.std_log == h.std == pytest.approx(0.28898, abs=1e-5)
    assert res.skew_station == h.skew
    assert res.skew_station_mse == pytest.approx(0.07074, abs=5e-5)
    # The bulletin's own factors give an adjusted skew of 0.0419, and it prints 0.0418,
    # so its weighted skew carries that slip. 1.5e-4 keeps it in range.
    assert res.skew_weighted == pytest.approx(-0.00409, abs=1.5e-4)
    assert res.skew_used == res.skew_weighted
    # Figure 6-1's frequency table. The 2 percent row's S x K is printed 0.59289 for
    # 0.59287, but its log Q and flow use 0.59287, so the flow is compared as printed.
    for row, log_q, q in zip(rows, res.quantiles.log_flow, res.quantiles.flow_cfs):
        assert log_q == pytest.approx(_p(row["log_Q"]), abs=2e-4)
        assert q == pytest.approx(_p(row["Q_cfs"]), rel=5e-4)


def test_big_sandy_plotting_positions_figure_6_1():
    pp = big_sandy().plotting_positions
    assert len(pp) == 47 and list(pp.year) == sorted(pp.year)
    pp = pp.set_index("year")
    rows = load("appendix6_big_sandy")["figure_6_1_table"]["rows_in_rank_order"]
    assert len(rows) == 47
    for row in rows:
        got = pp.loc[int(row["year"])]
        # Figure 6-1 carries W as 1.682 and the constant as 2.387, off by up to 0.008 at
        # E = 47, and rounds by up to 0.007 more. The largest differences are 0.0145 in
        # the order (E = 45) and 0.0224 in the percent (E = 40).
        assert got.weighted_order == pytest.approx(_p(row["weighted_order_m_tilde"]), abs=0.015)
        assert 100 * got.exceedance == pytest.approx(
            _p(row["plotting_position_weibull_percent"]), abs=0.0225
        )


def test_big_sandy_limits_use_the_systematic_length():
    """p. 9-3: N for the limits is the systematic record, 44, not H or N + Z."""
    res = big_sandy(generalized_skew=-0.2)
    q = res.quantiles
    k_up, k_lo = confidence_deviates(q.K_factor.to_numpy(), 44, 1.644853626951)
    assert q.upper.to_numpy() == pytest.approx(10 ** (res.mean_log + k_up * res.std_log))
    assert q.lower.to_numpy() == pytest.approx(10 ** (res.mean_log + k_lo * res.std_log))
    assert q.expected_probability.to_numpy() == pytest.approx(
        expected_probability(q.aep.to_numpy(), 44)
    )


def test_a_peak_moved_to_z_still_counts_in_the_limits():
    """Example 2: 38 peaks in the moments, but a systematic record of 39 years."""
    res = example2(historic_period=(1892, 1973))
    q = res.quantiles
    k_up, _ = confidence_deviates(q.K_factor.to_numpy(), 39, 1.644853626951)
    assert q.upper.to_numpy() == pytest.approx(10 ** (res.mean_log + k_up * res.std_log))


# --- What goes into Z --------------------------------------------------------------------


def test_systematic_peak_above_smallest_historic_joins_z():
    bs = record("appendix6_big_sandy")
    flows = list(bs.peaks_cfs)
    flows[5] = 19000.0  # above the 18,500 historic peak
    res = big_sandy(flows)
    assert res.historic.z == 4 and res.n_systematic == 43
    assert res.n_years == 44 and not res.high_outliers
    assert (bs.water_years[5], 19000.0) in res.historic.historic_peaks
    expected = historic_moments(
        np.log10([q for i, q in enumerate(flows) if i != 5]),
        np.log10([25000.0, 21000.0, 18500.0, 19000.0]),
        h=77,
        l=0,
    )
    assert res.historic.skew == pytest.approx(expected.skew, abs=1e-12)
    pp = res.plotting_positions
    assert len(pp) == 47  # no peak is plotted twice
    assert list(pp.sort_values("weighted_order").flow_cfs[:4]) == [25000, 21000, 19000, 18500]


def test_systematic_peak_equal_to_smallest_historic_stays_systematic():
    flows = list(record("appendix6_big_sandy").peaks_cfs)
    flows[5] = 18500.0
    res = big_sandy(flows)
    assert res.historic.z == 3 and res.n_systematic == 44
    pp = res.plotting_positions.set_index("year")
    assert pp.loc[1927].weighted_order == 3.0  # the historic peak ranks first in the tie
    w = res.historic.w
    assert pp.loc[1935].weighted_order == pytest.approx(w * 4 - (w - 1) * 3.5)


def test_high_outlier_below_the_historic_peaks_also_joins_z():
    res = synthetic(
        [9000.0, 1100.0], historic_peaks=[(1920, 50000.0)], historic_period=(1901, 1969)
    ).run_analysis()
    assert res.high_outliers == [(1968, 9000.0)]
    assert (res.historic.z, res.n_systematic) == (2, 19)
    assert res.historic.historic_peaks == [(1920, 50000.0), (1968, 9000.0)]


def test_period_with_nothing_historic_makes_no_adjustment():
    r = record("example1")
    plain = Bulletin17B(r.peaks_cfs, r.water_years, generalized_skew=0.6).run_analysis()
    res = Bulletin17B(
        r.peaks_cfs, r.water_years, generalized_skew=0.6, historic_period=(1900, 1968)
    ).run_analysis()
    assert res.historic is None
    assert any("no historic adjustment" in s for s in res.steps)
    assert res.quantiles.equals(plain.quantiles)


def test_weight_of_one_leaves_the_moments_alone():
    """The whole period is the record: W = 1 and eqs. 6-2 to 6-4 reduce to 2 to 4a."""
    r = record("example2")
    res = example2(historic_period=(1935, 1973))
    assert res.historic.w == 1.0
    assert (res.mean_log, res.std_log, res.skew_station) == pytest.approx(
        log_moments(np.log10(r.peaks_cfs)), abs=1e-12
    )


def test_eq6_takes_the_adjusted_skew_and_h():
    res = big_sandy(generalized_skew=-0.2)
    assert res.skew_station_mse == station_skew_mse(res.historic.skew, 77)
    assert res.skew_weighted == weighted_skew(res.historic.skew, res.skew_station_mse, -0.2, 0.302)


# --- The order of the tests ----------------------------------------------------------------


@pytest.mark.parametrize(
    ("skew", "order"),
    [
        (0.4, BOTH),
        (-0.4, BOTH),
        (0.0, BOTH),
        (math.nextafter(0.4, 1.0), "high first"),
        (math.nextafter(-0.4, -1.0), "low first"),
        (0.401, "high first"),
        (-0.401, "low first"),
    ],
)
def test_outlier_order_at_and_beyond_0_4(skew, order):
    assert _outlier_order(skew) == order


def test_high_first_runs_the_low_test_by_eq_8b_at_systematic_n():
    """K is read at N = 20, not at H = 220 (which is off the table) nor N - 1."""
    res = synthetic([9000.0, 1100.0], historic_period=(1750, 1969)).run_analysis()
    assert res.outlier_order == "high first" and res.historic.h == 220
    unadjusted = log_moments(np.log10(BASE + [9000.0, 1100.0]))
    assert res.high_outlier_threshold == pytest.approx(
        10 ** (unadjusted[0] + outlier_k(20) * unadjusted[1])
    )
    assert res.low_outlier_threshold == pytest.approx(
        10 ** (res.historic.mean - outlier_k(20) * res.historic.std)
    )
    text = " ".join(res.steps)
    assert "eq. 8b" in text and "eq. 8a" not in text


def test_high_first_without_adjustment_runs_the_low_test_by_eq_8a():
    res = synthetic([9000.0, 1100.0]).run_analysis()
    assert res.outlier_order == "high first" and res.historic is None
    mean, std, _ = log_moments(np.log10(BASE + [9000.0, 1100.0]))
    assert res.low_outlier_threshold == pytest.approx(10 ** (mean - outlier_k(20) * std))
    text = " ".join(res.steps)
    assert "eq. 8a" in text and "eq. 8b" not in text


def test_both_branch_never_uses_eq_8b():
    text = " ".join(big_sandy().steps)
    assert "eq. 8a" in text and "eq. 8b" not in text


def test_high_first_low_outlier_found_after_the_adjustment():
    """The high outlier hides the low one until the adjustment has been made."""
    b = synthetic([20000.0, 220.0], historic_period=(1870, 1969))
    x = np.log10(BASE + [20000.0, 220.0])
    mean, std, g = log_moments(x)
    assert g > 0.4 and x[-1] > mean - outlier_k(20) * std  # not low by eq. 8a
    first = historic_moments(np.delete(x, 18), x[18:19], h=100, l=0)
    s = b._screen_and_adjust([])
    assert s.high_outliers == [(1968, 20000.0)] and s.low_outliers == [(1969, 220.0)]
    assert s.low_threshold_log == pytest.approx(first.mean - outlier_k(20) * first.std)
    again = historic_moments(x[:18], x[18:19], h=100, l=1)
    assert (s.historic.l, s.historic.z, s.n_systematic) == (1, 1, 18)
    assert (s.mean, s.std, s.skew) == pytest.approx(again[1:], abs=1e-12)
    with pytest.raises(Bulletin17BError, match="not implemented"):
        b.run_analysis()


def test_low_first_high_outlier_found_after_the_low_is_dropped():
    b = synthetic([0.5, 3500.0], historic_period=(1870, 1969))
    x = np.log10(BASE + [0.5, 3500.0])
    mean, std, g = log_moments(x)
    assert g < -0.4 and x[-1] < mean + outlier_k(20) * std  # not high by the first statistics
    m2, s2, _ = log_moments(np.delete(x, 18))
    s = b._screen_and_adjust([])
    assert s.order == "low first"
    assert s.low_outliers == [(1968, 0.5)] and s.high_outliers == [(1969, 3500.0)]
    assert s.low_threshold_log == pytest.approx(mean - outlier_k(20) * std)
    assert s.high_threshold_log == pytest.approx(m2 + outlier_k(19) * s2)
    expected = historic_moments(x[:18], x[19:], h=100, l=1)
    assert s.historic.l == 1
    assert (s.mean, s.std, s.skew) == pytest.approx(expected[1:], abs=1e-12)
    with pytest.raises(Bulletin17BError, match="not implemented"):
        b.run_analysis()


def test_both_branch_drops_low_outliers_then_adjusts():
    """Both kinds flagged on unadjusted statistics, and L counts the low one."""
    b = synthetic([40000.0, 26.0], historic_period=(1870, 1969))
    x = np.log10(BASE + [40000.0, 26.0])
    mean, std, g = log_moments(x)
    assert abs(g) <= 0.4
    s = b._screen_and_adjust([])
    assert s.order == BOTH
    assert s.high_outliers == [(1968, 40000.0)] and s.low_outliers == [(1969, 26.0)]
    assert s.high_threshold_log == pytest.approx(mean + outlier_k(20) * std)
    assert s.low_threshold_log == pytest.approx(mean - outlier_k(20) * std)
    expected = historic_moments(x[:18], x[18:19], h=100, l=1)
    assert (s.historic.l, s.historic.z) == (1, 1)
    assert s.historic.w == pytest.approx(99 / 19)
    assert (s.mean, s.std, s.skew) == pytest.approx(expected[1:], abs=1e-12)
    with pytest.raises(Bulletin17BError, match="not implemented"):
        b.run_analysis()


def test_low_outlier_without_history_is_dropped_from_the_moments():
    s = synthetic([0.5, 1100.0])._screen_and_adjust([])
    x = np.log10(BASE + [1100.0])
    assert s.historic is None and s.low_outliers == [(1968, 0.5)]
    assert (s.mean, s.std, s.skew) == pytest.approx(log_moments(x), abs=1e-12)
    assert s.n_systematic == 19


# --- Refusals ------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("peaks", "period", "match"),
    [
        ([(1890, 30000.0)], (1897, 1973), "outside the historic period"),
        ([(1980, 30000.0)], (1897, 1973), "outside the historic period"),
        ([(1897, 30000.0)], None, "historic_period"),
        ([(1940, 30000.0)], (1897, 1973), "systematic record"),
        ([(1897, 30000.0), (1897, 31000.0)], (1897, 1973), "more than once"),
        ([(1897, 0.0)], (1897, 1973), "positive"),
        ([(1897, -5.0)], (1897, 1973), "positive"),
        ([(1897, math.nan)], (1897, 1973), "finite"),
        ([(1897.5, 30000.0)], (1897, 1973), "integers"),
        ([1897], (1897, 1973), "pairs"),
        ([], (1973, 1897), "period"),
        ([], (1897,), "period"),
        ([], 1897, "period"),
        ([], (1897, 1972), "cover"),
        ([], (1931, 1973), "cover"),
    ],
)
def test_bad_historic_input_refuses_with_the_reason(peaks, period, match):
    bs = record("appendix6_big_sandy")
    with pytest.raises(Bulletin17BError, match=match):
        Bulletin17B(
            bs.peaks_cfs, bs.water_years, historic_peaks=peaks, historic_period=period
        ).run_analysis()


def test_every_systematic_peak_above_the_historic_one_refuses():
    with pytest.raises(Bulletin17BError, match="systematic record is empty"):
        big_sandy(historic_peaks=[(1900, 1.0)])
