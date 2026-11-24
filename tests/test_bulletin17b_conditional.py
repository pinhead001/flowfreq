"""Bulletin17B's conditional probability adjustment (Appendix 5) and synthetic statistics.

Printed values are from Appendix 12 Examples 3 (a low outlier) and 4 (zero flood years).
The mixed records are checked against arithmetic written out here, with the Pearson III
deviates from scipy rather than ``flowfreq.core.kfactor``.
"""

from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np
import pytest
from scipy.stats import pearson3

from flowfreq.bulletin17b import (
    Bulletin17B,
    Bulletin17BError,
    confidence_deviates,
    expected_probability,
    log_moments,
    station_skew_mse,
    synthetic_statistics,
    weighted_skew,
)
from tests.fixtures.b17b import load, record

BOTH = "both on unadjusted statistics"
# 18 ordinary years, used to build records with outliers.
BASE = [1000.0, 1300, 800, 1150, 900, 1250, 700, 1050, 1400]
BASE = BASE + [q * 1.03 for q in BASE]


def _p(text):
    return float(str(text).replace(",", ""))


def example3(**kw):
    r = record("example3")
    return Bulletin17B(r.peaks_cfs, r.water_years, generalized_skew=0.5, **kw).run_analysis()


def example4(**kw):
    r = record("example4")
    return Bulletin17B(r.peaks_cfs, r.water_years, generalized_skew=-0.3, **kw).run_analysis()


def example1(**kw):
    r = record("example1")
    return Bulletin17B(r.peaks_cfs, r.water_years, **kw).run_analysis()


def synthetic(extra, **kw):
    """The 18 ordinary years and then ``extra``, from 1950."""
    flows = BASE + list(extra)
    return Bulletin17B(flows, list(range(1950, 1950 + len(flows))), **kw).run_analysis()


def final_rows(name):
    table = "final_frequency_table_12_9" if name == "example3" else "final_frequency_table_12_11"
    rows = load(name)[table]["rows"]
    assert len(rows) == 9
    return [(_p(row["P"]), _p(row["Q_cfs"])) for row in rows]


# --- Independent arithmetic ---------------------------------------------------------------


def pearson_k(skew, p):
    return float(pearson3.isf(p, skew))


class Synthetic(NamedTuple):
    p_tilde: float
    q01: float
    q10: float
    q50: float
    mean: float
    std: float
    skew: float


def adjusted_moments(xs, xz, h, low):
    """Eqs. 6-1 to 6-4a, written out. Returns W, mean, S and G."""
    w = (h - len(xz)) / (len(xs) + low)
    m = h - w * low
    mean = (w * math.fsum(xs) + math.fsum(xz)) / m

    def central(power):
        return w * math.fsum((x - mean) ** power for x in xs) + math.fsum(
            (x - mean) ** power for x in xz
        )

    s = math.sqrt(central(2) / (m - 1))
    return w, mean, s, m / ((m - 1) * (m - 2)) * central(3) / s**3


def synthetic_by_hand(p_tilde, mean, std, skew):
    """Eqs. 5-2 to 5-5 on the conditional curve."""
    l01, l10, l50 = (mean + pearson_k(skew, p / p_tilde) * std for p in (0.01, 0.10, 0.50))
    g = -2.50 + 3.12 * (l01 - l10) / (l10 - l50)
    s = (l01 - l50) / (pearson_k(g, 0.01) - pearson_k(g, 0.50))
    return Synthetic(p_tilde, 10**l01, 10**l10, 10**l50, l50 - pearson_k(g, 0.50) * s, s, g)


def check_conditional(res, want: Synthetic, n_truncated, n_mse, n_years):
    """The adjustment, eq. 6 and the final curve against the arithmetic above."""
    c = res.conditional
    assert c.n_truncated == n_truncated and res.n_years == n_years
    assert c.p_tilde == pytest.approx(want.p_tilde, abs=1e-12)
    assert (c.q01, c.q10, c.q50) == pytest.approx((want.q01, want.q10, want.q50), rel=1e-6)
    assert (c.mean_s, c.std_s, c.skew_s) == pytest.approx(
        (want.mean, want.std, want.skew), abs=1e-6
    )
    assert (res.mean_log, res.std_log, res.skew_station) == (c.mean_s, c.std_s, c.skew_s)
    assert res.skew_station_mse == pytest.approx(station_skew_mse(want.skew, n_mse), abs=1e-6)
    q = res.quantiles
    k = np.array([pearson_k(res.skew_used, p) for p in q.aep])
    assert q.flow_cfs.to_numpy() == pytest.approx(10 ** (want.mean + k * want.std), rel=1e-5)
    up, lo = confidence_deviates(k, n_years, 1.644853626951)
    assert q.upper.to_numpy() == pytest.approx(10 ** (want.mean + up * want.std), rel=1e-5)
    assert q.lower.to_numpy() == pytest.approx(10 ** (want.mean + lo * want.std), rel=1e-5)
    assert q.expected_probability.to_numpy() == pytest.approx(
        expected_probability(q.aep.to_numpy(), n_years)
    )


# --- Example 3, a low outlier -------------------------------------------------------------


def test_example3_low_outlier_and_synthetic_statistics():
    res = example3(hand_computation=True)
    assert res.outlier_order == "low first"
    assert [y for y, _ in res.low_outliers] == [1969] and res.high_outliers == []
    assert res.low_outlier_threshold == pytest.approx(946, rel=1e-3)
    assert res.high_outlier_threshold == pytest.approx(22760, rel=2e-3)
    assert res.n_years == 38 and res.n_systematic == 37 and res.n_zeros == res.n_below_base == 0
    c = res.conditional
    assert c.n_truncated == 1 and c.p_tilde == pytest.approx(0.9737, abs=1e-4)
    assert (c.mean, c.std) == (pytest.approx(3.7488, abs=1e-4), pytest.approx(0.2296, abs=1e-4))
    assert c.skew == pytest.approx(0.6311, abs=2e-4)
    # Read off Figure 12-5 in the bulletin: 23880, 11210 and 5230.
    assert (c.q01, c.q10, c.q50) == pytest.approx((23880, 11210, 5230), rel=1e-3)
    assert c.skew_s == pytest.approx(0.5948, abs=0.02)
    assert c.std_s == pytest.approx(0.2310, abs=1e-4)
    assert c.mean_s == pytest.approx(3.7415, abs=2e-4)
    assert res.skew_station_mse == pytest.approx(0.183, abs=0.005)
    assert res.skew_weighted == pytest.approx(0.5590, abs=1e-3)
    assert res.skew_used == 0.6
    got = res.quantiles.set_index("aep").flow_cfs
    for p, q in final_rows("example3"):
        assert got[p] == pytest.approx(q, rel=0.02)


def test_example3_exact_mode_is_close():
    got = example3().quantiles.set_index("aep").flow_cfs
    for p, q in final_rows("example3"):
        assert got[p] == pytest.approx(q, rel=0.06)


def test_exact_mode_rounds_nothing():
    res = example3()
    x = np.log10([q for q in record("example3").peaks_cfs if q != 536.0])
    want = synthetic_by_hand(37 / 38, *log_moments(x))
    check_conditional(res, want, n_truncated=1, n_mse=38, n_years=38)
    assert res.skew_used == res.skew_weighted
    assert res.skew_weighted == pytest.approx(
        weighted_skew(want.skew, station_skew_mse(want.skew, 38), 0.5, 0.302), abs=1e-6
    )


def test_hand_mode_rounds_the_conditional_skew_and_the_k_skew():
    c = example3(hand_computation=True).conditional
    p = 37 / 38
    by_hand = [10 ** (c.mean + pearson_k(0.6, a / p) * c.std) for a in (0.01, 0.10, 0.50)]
    assert (c.q01, c.q10, c.q50) == pytest.approx(by_hand, rel=1e-6)
    assert (c.mean_s, c.std_s, c.skew_s) == pytest.approx(
        synthetic_statistics(c.q01, c.q10, c.q50, k_skew=0.6)
    )


# --- Example 4, zero flood years ------------------------------------------------------------


def test_example4_zero_years():
    res = example4(hand_computation=True)
    assert res.n_zeros == 6 and [y for y, _ in res.low_outliers] == [1955]
    assert res.outlier_order == "low first" and res.high_outliers == []
    assert res.low_outlier_threshold == pytest.approx(23.9, rel=2e-3)
    assert res.high_outlier_threshold == pytest.approx(41770, rel=2e-3)
    assert res.n_years == 42 and res.n_systematic == 35 and res.n_below_base == 0
    c = res.conditional
    assert c.n_truncated == 7 and c.p_tilde == pytest.approx(35 / 42, abs=1e-5)
    assert (c.mean, c.std, c.skew) == pytest.approx((3.1321, 0.5665, -0.4396), abs=1e-4)
    # Read off Figure 12-7 in the bulletin: 17940, 6000 and 1060.
    assert (c.q01, c.q10, c.q50) == pytest.approx((17940, 6000, 1060), rel=3e-3)
    assert c.skew_s == pytest.approx(-0.5287, abs=0.03)
    assert c.mean_s == pytest.approx(2.9708, abs=5e-4)
    assert c.std_s == pytest.approx(0.6564, abs=5e-4)
    assert res.skew_station_mse == pytest.approx(0.163, abs=1e-3)


def test_example4_final_table_at_the_printed_skew():
    """Table 12-11's rows from the class's synthetic statistics and the printed -0.4.

    The bulletin weights its graph-read G_s of -0.5287 to -0.4487, which rounds to -0.4.
    Evaluated exactly, the adjusted curve gives 17954, 6016 and 1059 cfs, not 17940, 6000
    and 1060, so G_s is -0.5358 and the weighted skew -0.4530, which rounds to -0.5. The
    hand-mode quantiles of the class are therefore at -0.5, and the 1 percent flood is 11
    percent below the printed 20100.
    """
    res = example4(hand_computation=True)
    assert res.skew_weighted == pytest.approx(-0.4487, abs=0.005)
    assert res.skew_weighted == pytest.approx(-0.45296, abs=1e-5)
    assert res.skew_used == -0.5
    got = res.quantiles.set_index("aep").flow_cfs
    for p, q in final_rows("example4"):
        at_printed_skew = 10 ** (res.mean_log + pearson_k(-0.4, p) * res.std_log)
        assert at_printed_skew == pytest.approx(q, rel=0.006)
        # The whole difference in the class's own rows is the skew's.
        shift = (pearson_k(-0.5, p) - pearson_k(-0.4, p)) * res.std_log
        assert got[p] == pytest.approx(q * 10**shift, rel=0.006)


def test_example4_printed_flows_give_the_printed_answer():
    """The same equations at the bulletin's graph readings give -0.4487, -0.4 and 20100."""
    mean, std, g = synthetic_statistics(17940, 6000, 1060, k_skew=-0.5)
    assert (mean, std, g) == pytest.approx((2.9708, 0.6564, -0.5287), abs=1e-4)
    assert weighted_skew(g, 0.163, -0.3, 0.302) == pytest.approx(-0.4487, abs=2e-4)
    assert 10 ** (mean + pearson_k(-0.4, 0.01) * std) == pytest.approx(20100, rel=0.003)


def test_example4_exact_mode_is_close():
    """Exact mode uses -0.4396 and -0.4775 where the bulletin rounds both to -0.4.

    Within 5 percent from AEP 0.50 to 0.005. The tails carry the unrounded skews: -17
    percent at 0.99 (an extrapolation below P~ = 0.83), and -7 percent at 0.90 and 0.002.
    """
    res = example4()
    got = res.quantiles.set_index("aep").flow_cfs
    wider = {0.99: 0.17, 0.90: 0.07, 0.002: 0.07}
    for p, q in final_rows("example4"):
        assert got[p] == pytest.approx(q, rel=wider.get(p, 0.05))
    x = np.log10([q for q in record("example4").peaks_cfs if q > 16.0])
    want = synthetic_by_hand(35 / 42, *log_moments(x))
    check_conditional(res, want, n_truncated=7, n_mse=42, n_years=42)


def test_example4_plotting_positions_give_truncated_years_ranks():
    pp = example4().plotting_positions
    assert len(pp) == 42 and list(pp.year) == sorted(pp.year)
    by_year = pp.set_index("year")
    assert by_year.loc[1955].weighted_order == 36.0 and by_year.loc[1955].flow_cfs == 16.0
    zeros = pp[pp.flow_cfs == 0.0]
    assert sorted(zeros.weighted_order) == [37, 38, 39, 40, 41, 42]
    assert by_year.loc[1958].exceedance == pytest.approx(1 / 43)


def test_steps_record_the_truncation_and_the_synthetic_statistics():
    res = example4(hand_computation=True)
    text = " ".join(res.steps)
    assert "6 zero year(s)" in text and "1947" in text
    assert "7 of 42 years truncated" in text and "above 16 cfs" in text and "115 cfs" in text
    assert "0.8333" in text and "eq. 5-1a" in text and "5-1b" not in text
    assert f"{res.conditional.skew_s:.4f}" in text and f"{res.mean_log:.4f}" in text
    assert "extrapolat" in text
    note = next(t for t in res.steps if "extrapolat" in t)
    assert note.startswith("AEP 0.99, 0.9: above P~ = 0.8333")


def test_no_extrapolation_note_when_every_aep_is_below_p_tilde():
    assert not any("extrapolat" in t for t in example4(aeps=[0.5, 0.01]).steps)


def test_a_plain_record_has_no_conditional_steps():
    res = example1()
    assert res.conditional is None
    assert not any("truncat" in t or "extrapolat" in t for t in res.steps)


# --- The 25 percent limit --------------------------------------------------------------------


def test_user_threshold_truncating_too_much_refuses():
    with pytest.raises(Bulletin17BError, match="25 percent"):
        example1(low_outlier_threshold=2190.0)


@pytest.mark.parametrize(("zeros", "refused"), [(6, False), (7, True)])
def test_exactly_25_percent_is_allowed(zeros, refused):
    """6 zeros in 24 years is 25 percent. 7 in 25 is 28."""
    flows = BASE + [0.0] * zeros
    b = Bulletin17B(flows, list(range(1950, 1950 + len(flows))))
    if refused:
        with pytest.raises(Bulletin17BError, match="7 of the 25.*25 percent"):
            b.run_analysis()
    else:
        assert b.run_analysis().conditional.p_tilde == 0.75


@pytest.mark.parametrize(("truncated", "refused"), [(6, False), (7, True)])
def test_the_limit_counts_low_outliers(truncated, refused):
    ranked = sorted(record("example1").peaks_cfs)
    assert len(ranked) == 24 and ranked[truncated - 1] < ranked[truncated]
    if refused:
        with pytest.raises(Bulletin17BError, match="7 of the 24.*25 percent"):
            example1(low_outlier_threshold=ranked[truncated])
    else:
        res = example1(low_outlier_threshold=ranked[truncated])
        assert len(res.low_outliers) == 6 and res.conditional.p_tilde == 0.75


def test_the_limit_counts_zeros_below_base_and_lows_together():
    """3 zeros, 2 below-base years and 2 lows of 25 years: 28 percent."""
    flows = BASE + [0.0] * 3 + [500.0, 510.0]
    years = list(range(1950, 1950 + len(flows)))
    with pytest.raises(Bulletin17BError, match="7 of the 25.*25 percent"):
        Bulletin17B(
            flows, years, below_base_years=[1990, 1991], low_outlier_threshold=600.0
        ).run_analysis()


def test_the_limit_is_given_before_the_ten_peak_minimum():
    with pytest.raises(Bulletin17BError, match="25 percent"):
        Bulletin17B(BASE[:9] + [0.0] * 9, list(range(1950, 1968))).run_analysis()
    with pytest.raises(Bulletin17BError, match="25 percent"):
        Bulletin17B([0.0] * 12, list(range(1950, 1962))).run_analysis()


def test_fewer_than_ten_peaks_left_after_zeros_and_below_base_years_refuses():
    """9 peaks of 11 years: 18 percent truncated, but too few left."""
    with pytest.raises(
        Bulletin17BError, match="9 of the 11.*1 zero year.*1 year.s. below.*at least 10"
    ):
        Bulletin17B(
            BASE[:9] + [0.0], list(range(1950, 1960)), below_base_years=[1960]
        ).run_analysis()


def test_ten_left_after_zeros_is_enough():
    res = Bulletin17B(BASE[:10] + [0.0], list(range(1950, 1961))).run_analysis()
    assert (res.n_systematic, res.n_years) == (10, 11)


def test_p_tilde_is_above_one_half_without_history():
    """The 25 percent limit keeps P~ at 0.75 or more, so the 0.50 point exists."""
    for zeros in range(1, 7):
        res = synthetic([0.0] * zeros)
        assert res.conditional.p_tilde >= 0.75


def test_p_tilde_at_one_half_or_less_refuses():
    """Only with history can it happen: most of the record moved to Z leaves W large."""
    rest = [100.0 + 4 * i for i in range(10)]
    above = [160.0 + 2 * i for i in range(80)]
    flows = rest + above + [0.0] * 30
    with pytest.raises(Bulletin17BError, match="0.4525"):
        Bulletin17B(
            flows,
            list(range(1851, 1971)),
            historic_peaks=[(1700, 150.0)],
            historic_period=(1671, 1970),
        ).run_analysis()


# --- Below-base years ---------------------------------------------------------------------


def test_below_base_years_are_truncated():
    r = record("example1")
    res = Bulletin17B(
        r.peaks_cfs, r.water_years, below_base_years=[r.water_years[0]]
    ).run_analysis()
    assert res.n_below_base == 1 and res.conditional.p_tilde == pytest.approx(23 / 24)
    assert (res.n_years, res.n_systematic, res.n_zeros) == (24, 23, 0)
    want = synthetic_by_hand(23 / 24, *log_moments(np.log10(r.peaks_cfs[1:])))
    check_conditional(res, want, n_truncated=1, n_mse=24, n_years=24)


def test_below_base_year_with_or_without_a_flow_is_the_same():
    """The record value is not used, whatever it is, and N for the limits is the same."""
    r = record("example1")
    y0 = r.water_years[0]
    without = Bulletin17B(r.peaks_cfs[1:], r.water_years[1:], below_base_years=[y0]).run_analysis()
    for flow in (r.peaks_cfs[0], 5.0, 0.0, 1e6):
        res = Bulletin17B(
            [flow, *r.peaks_cfs[1:]], r.water_years, below_base_years=[y0]
        ).run_analysis()
        assert res.quantiles.equals(without.quantiles)
        assert res.plotting_positions.equals(without.plotting_positions)
        assert (res.n_years, res.n_below_base, res.n_zeros) == (24, 1, 0)


def test_below_base_years_in_the_plotting_positions():
    """Ranked below every measured peak and above the zeros, with no flow."""
    flows = BASE + [0.0]
    res = Bulletin17B(flows, list(range(1950, 1969)), below_base_years=[1940, 1975]).run_analysis()
    pp = res.plotting_positions
    assert len(pp) == res.n_years == 21
    by_year = pp.set_index("year")
    assert sorted(by_year.loc[[1940, 1975]].weighted_order) == [19, 20]
    assert by_year.loc[[1940, 1975]].flow_cfs.isna().all()
    assert by_year.loc[1968].weighted_order == 21.0
    assert by_year.exceedance.min() == pytest.approx(1 / 22)


@pytest.mark.parametrize(
    ("years", "match"),
    [
        ([1990, 1990], "more than once"),
        ([1990.5], "integers"),
        (["x"], "number"),
        ([True], "integers"),
        (1990, "sequence"),
    ],
)
def test_bad_below_base_years_refuse(years, match):
    with pytest.raises(Bulletin17BError, match=match):
        example1(below_base_years=years)


def test_below_base_year_that_is_a_historic_peak_refuses():
    with pytest.raises(Bulletin17BError, match="1920.*systematic record"):
        example1(
            below_base_years=[1920],
            historic_peaks=[(1920, 30000.0)],
            historic_period=(1900, 1973),
        )


def test_the_historic_period_has_to_cover_the_below_base_years():
    with pytest.raises(Bulletin17BError, match="cover"):
        example1(below_base_years=[1890], historic_period=(1900, 1973))
    with pytest.raises(Bulletin17BError, match="cover"):
        example1(below_base_years=[1990], historic_period=(1900, 1973))


@pytest.mark.parametrize("kw", [{}, {"below_base_years": [2005]}])
def test_a_negative_flow_refuses(kw):
    flows = [1000.0 + 37 * i for i in range(12)]
    flows[5] = -1.0
    with pytest.raises(Bulletin17BError, match="negative"):
        Bulletin17B(flows, list(range(2000, 2012)), **kw)


# --- The user's threshold ---------------------------------------------------------------------


def test_user_threshold_replaces_eq_8a():
    """600 cfs drops what eq. 8a drops (1969's 536 cfs), so only the threshold differs."""
    res, plain = example3(low_outlier_threshold=600.0), example3()
    assert res.low_outlier_threshold == 600.0 and res.low_outliers == [(1969, 536.0)]
    assert res.quantiles.equals(plain.quantiles)
    # The high test is made on the record without the low, as in the plain run.
    assert res.high_outlier_threshold == plain.high_outlier_threshold
    note = next(t for t in res.steps if t.startswith("Low outlier"))
    assert "600" in note and "replaces eq. 8a" in note


def test_user_threshold_is_strict():
    """A peak equal to it is kept, although eq. 8a (946 cfs) would drop it."""
    res = example3(low_outlier_threshold=536.0)
    assert res.low_outliers == [] and res.conditional is None
    assert res.n_systematic == 38


def test_user_threshold_replaces_eq_8b():
    kw = {"historic_period": (1870, 1969)}
    b = Bulletin17B(BASE + [20000.0, 220.0], list(range(1950, 1970)), **kw)
    assert b._screen_and_adjust([]).low_outliers == [(1969, 220.0)]
    kept = synthetic([20000.0, 220.0], low_outlier_threshold=100.0, **kw)
    assert kept.low_outliers == [] and kept.conditional is None and kept.historic.l == 0
    assert any("replaces eq. 8b" in t for t in kept.steps)


@pytest.mark.parametrize(
    ("order_extra", "order"),
    [([2600.0, 1100.0], "high first"), ([400.0, 1100.0], "low first"), ([1900.0, 1100.0], BOTH)],
)
def test_user_threshold_in_every_order(order_extra, order):
    """Truncated peaks are those below it, and the moments are of the peaks kept."""
    res = synthetic(order_extra, low_outlier_threshold=720.0)
    assert res.outlier_order == order and res.low_outlier_threshold == 720.0
    kept = [q for q in BASE + order_extra if q >= 720.0]
    assert sorted(q for _, q in res.low_outliers) == sorted(
        q for q in BASE + order_extra if q < 720.0
    )
    want = synthetic_by_hand(len(kept) / 20, *log_moments(np.log10(kept)))
    check_conditional(res, want, n_truncated=20 - len(kept), n_mse=20, n_years=20)


@pytest.mark.parametrize("threshold", [0.0, -5.0, math.nan, math.inf, "x"])
def test_bad_user_threshold_refuses(threshold):
    with pytest.raises(Bulletin17BError, match="low_outlier_threshold"):
        example1(low_outlier_threshold=threshold)


# --- High and low outliers with a historic period -------------------------------------------


@pytest.mark.parametrize(
    ("extra", "order", "i_high", "i_low"),
    [
        ([6000.0, 320.0], "high first", 18, 19),
        ([0.5, 3500.0], "low first", 19, 18),
        ([2900.0, 370.0], BOTH, 18, 19),
    ],
)
def test_mixed_outliers_with_a_historic_period(extra, order, i_high, i_low):
    """Each skew band end to end: Z = 1, L = 1, H = 100 and P~ by eq. 5-1b."""
    res = synthetic(extra, historic_period=(1870, 1969))
    flows = BASE + extra
    x = np.log10(flows)
    n = 20
    mean = math.fsum(x) / n
    s = math.sqrt(math.fsum((x - mean) ** 2) / (n - 1))
    g = n * math.fsum((x - mean) ** 3) / ((n - 1) * (n - 2) * s**3)
    assert order == ("high first" if g > 0.4 else "low first" if g < -0.4 else BOTH)
    assert res.outlier_order == order
    assert res.high_outliers == [(1950 + i_high, flows[i_high])]
    assert res.low_outliers == [(1950 + i_low, flows[i_low])]
    if order == "high first":
        _, m1, s1, _ = adjusted_moments(np.delete(x, i_high), x[[i_high]], 100, 0)
        assert res.low_outlier_threshold == pytest.approx(10 ** (m1 - 2.385 * s1))
    else:
        assert res.low_outlier_threshold == pytest.approx(10 ** (mean - 2.385 * s))
    if order == "low first":
        m2, s2, _ = log_moments(np.delete(x, i_low))
        assert res.high_outlier_threshold == pytest.approx(10 ** (m2 + 2.361 * s2))
    else:
        assert res.high_outlier_threshold == pytest.approx(10 ** (mean + 2.385 * s))
    w, m_adj, s_adj, g_adj = adjusted_moments(x[:18], x[[i_high]], 100, 1)
    assert w == pytest.approx(99 / 19)
    h = res.historic
    assert (h.h, h.z, h.l, res.n_systematic) == (100, 1, 1, 18)
    assert (h.w, h.mean, h.std, h.skew) == pytest.approx((w, m_adj, s_adj, g_adj))
    want = synthetic_by_hand((100 - w * 1) / 100, m_adj, s_adj, g_adj)
    check_conditional(res, want, n_truncated=1, n_mse=100, n_years=20)
    c = res.conditional
    assert (c.mean, c.std, c.skew) == (h.mean, h.std, h.skew)
    assert any("eq. 5-1b" in t for t in res.steps)
    pp = res.plotting_positions
    assert len(pp) == 20 and pp.exceedance.max() == pytest.approx((w * 20 - (w - 1) * 1.5) / 101)


def test_zeros_with_a_historic_period():
    """Example 4 with a historic peak: L = 6 zeros + 1 low, P~ = (H - W L) / H."""
    res = example4(historic_peaks=[(1911, 12000.0)], historic_period=(1900, 1973))
    r = record("example4")
    x = np.log10([q for q in r.peaks_cfs if q > 16.0])
    w, mean, s, g = adjusted_moments(x, [math.log10(12000.0)], 74, 7)
    assert w == pytest.approx(73 / 42)
    h = res.historic
    assert (h.h, h.z, h.l, res.n_systematic, res.n_zeros) == (74, 1, 7, 35, 6)
    assert (h.w, h.mean, h.std, h.skew) == pytest.approx((w, mean, s, g))
    p_tilde = (74 - w * 7) / 74
    assert p_tilde == pytest.approx(0.835586, abs=1e-6)
    want = synthetic_by_hand(p_tilde, mean, s, g)
    check_conditional(res, want, n_truncated=7, n_mse=74, n_years=42)
    pp = res.plotting_positions
    assert len(pp) == 43  # the 42 years and the historic peak
    last = w * 43 - (w - 1) * 1.5
    assert pp.weighted_order.max() == pytest.approx(last)
    assert pp.exceedance.max() == pytest.approx(last / 75)


def test_every_pass_counts_zeros_and_below_base_years_in_l():
    """High first: the first pass has L = 3, which sets eq. 8b, and the second L = 4."""
    flows = BASE + [6000.0, 320.0, 0.0]
    years = list(range(1950, 1971))
    res = Bulletin17B(
        flows, years, below_base_years=[1945, 1971], historic_period=(1872, 1971)
    ).run_analysis()
    x = np.log10(flows[:20])
    _, m1, s1, _ = adjusted_moments(np.delete(x, 18), x[[18]], 100, 3)
    assert res.outlier_order == "high first"
    assert res.low_outlier_threshold == pytest.approx(10 ** (m1 - 2.385 * s1))
    w, mean, s, g = adjusted_moments(x[:18], x[[18]], 100, 4)
    assert w == pytest.approx(99 / 22)
    h = res.historic
    assert (h.l, h.z, h.mean, h.std, h.skew) == pytest.approx((4, 1, mean, s, g))
    want = synthetic_by_hand((100 - 4 * w) / 100, mean, s, g)
    check_conditional(res, want, n_truncated=4, n_mse=100, n_years=23)
    made = [t for t in res.steps if "historic adjustment of Appendix 6" in t]
    assert "L = 3" in made[0] and "L = 4" in made[1]
    assert len(res.plotting_positions) == 23


def test_the_historic_weight_counts_below_base_years_outside_the_flows():
    bs = record("appendix6_big_sandy")
    res = Bulletin17B(
        bs.peaks_cfs[1:],
        bs.water_years[1:],
        below_base_years=[bs.water_years[0]],
        historic_peaks=bs.historic_peaks,
        historic_period=bs.historic_period,
    ).run_analysis()
    assert (res.historic.l, res.n_years) == (1, 44)
    assert res.historic.w == pytest.approx(74 / 44)
    assert res.conditional.p_tilde == pytest.approx((77 - 74 / 44) / 77)


def test_low_outlier_without_history_uses_eq_5_1a():
    res = synthetic([0.5, 1100.0])
    assert res.historic is None and res.low_outliers == [(1968, 0.5)]
    want = synthetic_by_hand(19 / 20, *log_moments(np.log10(BASE + [1100.0])))
    check_conditional(res, want, n_truncated=1, n_mse=20, n_years=20)


def test_a_synthetic_skew_outside_eq_5_3_refuses():
    flows = [1.0 + 0.001 * i for i in range(30)] + [5.0, 60.0, 0.0]
    with pytest.raises(Bulletin17BError, match="eq. 5-3"):
        Bulletin17B(flows, list(range(1950, 1983))).run_analysis()
