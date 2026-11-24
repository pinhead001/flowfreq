"""Bulletin17B on a complete record, against the bulletin's Example 1 (Appendix 12)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from flowfreq.bulletin17b import (
    STANDARD_AEP,
    Bulletin17B,
    Bulletin17BError,
    _nearest_tenth,
    log_moments,
)
from flowfreq.bulletin17c import FloodFrequencyAnalysis
from tests.fixtures.b17b import load, record

EX1_Q = {
    0.99: 841,
    0.90: 1200,
    0.50: 2190,
    0.10: 4960,
    0.05: 6530,
    0.02: 9110,
    0.01: 11500,
    0.005: 14500,
    0.002: 19200,
}  # Table 12-3


def ex1(**kw):
    r = record("example1")
    return Bulletin17B(r.peaks_cfs, r.water_years, generalized_skew=0.6, **kw).run_analysis()


def _p(text):
    return float(str(text).replace(",", ""))


def test_standard_aep_is_the_bulletins_rows():
    assert STANDARD_AEP == (0.99, 0.90, 0.50, 0.10, 0.05, 0.02, 0.01, 0.005, 0.002)


def test_example1_hand_computation_reproduces_table_12_3():
    res = ex1(hand_computation=True)
    assert res.outlier_order == "high first" and not res.high_outliers and not res.low_outliers
    assert res.skew_station_mse == pytest.approx(0.277, abs=5e-4)
    assert res.skew_weighted == pytest.approx(0.6678, abs=1e-4) and res.skew_used == 0.7
    got = dict(zip(res.quantiles.aep, res.quantiles.flow_cfs))
    for aep, q in EX1_Q.items():
        assert got[aep] == pytest.approx(q, rel=0.006)  # printed to 3 significant figures


def test_example1_confidence_limits_table_12_4():
    res = ex1(hand_computation=True)
    rows = res.quantiles.set_index("aep")
    for row in load("example1")["confidence_limits_table_12_4"]["rows"]:
        got = rows.loc[_p(row["P"])]
        assert got.upper == pytest.approx(_p(row["upper_Q_cfs"]), rel=0.006)
        assert got.lower == pytest.approx(_p(row["lower_Q_cfs"]), rel=0.006)
    assert rows.loc[0.99].upper == pytest.approx(1100, rel=0.006)
    assert rows.loc[0.99].lower == pytest.approx(568, rel=0.006)


def test_example1_exact_mode_is_close():
    got = dict(zip(ex1().quantiles.aep, ex1().quantiles.flow_cfs))
    for aep, q in EX1_Q.items():
        assert got[aep] == pytest.approx(q, rel=0.05)


def test_exact_mode_does_not_round_the_skew():
    res = ex1()
    assert res.skew_used == res.skew_weighted != 0.7


def test_example1_expected_probability_table_12_5():
    res = ex1(hand_computation=True)
    rows = res.quantiles.set_index("aep")
    for row in load("example1")["expected_probability_table_12_5"]["rows"]:
        assert rows.loc[_p(row["P"])].expected_probability == pytest.approx(
            _p(row["expected_probability"]), rel=0.06
        )


def test_plain_record_has_no_conditional_adjustment():
    res = ex1()
    assert res.conditional is None and res.historic is None
    assert res.n_years == res.n_systematic == 24
    assert res.n_zeros == res.n_below_base == 0


def test_quantile_table_columns_and_plotting_positions():
    res = ex1()
    assert list(res.quantiles.columns) == [
        "aep",
        "flow_cfs",
        "log_flow",
        "K_factor",
        "lower",
        "upper",
        "expected_probability",
    ]
    pp = res.plotting_positions
    assert list(pp.columns) == ["year", "flow_cfs", "weighted_order", "exceedance"]
    top = pp.sort_values("exceedance").iloc[0]
    assert top.flow_cfs == max(record("example1").peaks_cfs)
    assert top.exceedance == pytest.approx(1 / 25)
    assert all(isinstance(s, str) for s in res.steps)


def test_unsorted_input_is_sorted_by_year():
    r = record("example1")
    a = ex1()
    b = Bulletin17B(r.peaks_cfs[::-1], r.water_years[::-1], generalized_skew=0.6).run_analysis()
    assert b.skew_used == pytest.approx(a.skew_used)
    assert list(b.plotting_positions.year) == list(a.plotting_positions.year)


def test_high_outliers_are_kept_without_a_historic_period():
    flows = [1000, 1100, 900, 1200, 950, 1050, 1000, 980, 1020, 1080, 990, 100000]
    res = Bulletin17B(flows, list(range(2000, 2012))).run_analysis()
    assert res.high_outliers == [(2011, 100000.0)]
    assert res.n_systematic == 12
    assert any("high outlier" in s.lower() and "kept" in s for s in res.steps)


def test_duplicate_year_refuses():
    with pytest.raises(Bulletin17BError, match="year"):
        Bulletin17B([100, 200, 300] * 4, [2000, 2000] + list(range(2001, 2011))).run_analysis()


def test_constant_record_refuses():
    with pytest.raises(Bulletin17BError):
        Bulletin17B([500.0] * 12, list(range(2000, 2012))).run_analysis()


def test_no_generalized_skew_uses_station_skew():
    r = record("example1")
    res = Bulletin17B(r.peaks_cfs, r.water_years).run_analysis()
    assert res.skew_weighted is None and res.skew_used == pytest.approx(0.7300, abs=2e-4)
    assert any("station skew" in s for s in res.steps)


@pytest.mark.parametrize(
    "kw",
    [
        {"below_base_years": [1950]},
        {"low_outlier_threshold": 500.0},
    ],
)
def test_later_task_inputs_say_not_implemented(kw):
    r = record("example1")
    with pytest.raises(Bulletin17BError, match="not implemented"):
        Bulletin17B(r.peaks_cfs, r.water_years, **kw).run_analysis()


def test_zero_flow_says_not_implemented():
    flows = [0.0] + [1000.0 + 37 * i for i in range(11)]
    with pytest.raises(Bulletin17BError, match="not implemented"):
        Bulletin17B(flows, list(range(2000, 2012))).run_analysis()


def test_low_outlier_found_says_not_implemented():
    flows = [1000, 1100, 900, 1200, 950, 1050, 1000, 980, 1020, 1080, 990, 5]
    with pytest.raises(Bulletin17BError, match="not implemented"):
        Bulletin17B(flows, list(range(2000, 2012))).run_analysis()


@pytest.mark.parametrize(
    "kw",
    [
        {"peak_flows": [1, 2, 3], "water_years": [2000, 2001]},
        {"peak_flows": [1, 2, float("nan")], "water_years": [1, 2, 3]},
        {"peak_flows": [1, -2, 3], "water_years": [1, 2, 3]},
        {"peak_flows": [1, 2, 3], "water_years": [1.5, 2, 3]},
        {"peak_flows": [1, 2, 3], "water_years": [1, 2, 3], "confidence": 1.0},
        {"peak_flows": [1, 2, 3], "water_years": [1, 2, 3], "confidence": 0.0},
        {"peak_flows": [1, 2, 3], "water_years": [1, 2, 3], "aeps": [0.0, 0.5]},
        {"peak_flows": [1, 2, 3], "water_years": [1, 2, 3], "aeps": [1.0]},
        {"peak_flows": [1, 2, 3], "water_years": [1, 2, 3], "generalized_skew_mse": 0.0},
        {"peak_flows": [1, 2, 3], "water_years": [1, 2, 3], "generalized_skew": math.inf},
    ],
)
def test_invalid_inputs_refuse(kw):
    with pytest.raises(Bulletin17BError):
        Bulletin17B(**kw).run_analysis()


def test_short_record_refuses_with_the_reason():
    with pytest.raises(Bulletin17BError, match="10"):
        Bulletin17B([100, 200, 300, 250, 400], list(range(2000, 2005))).run_analysis()


@pytest.mark.parametrize(
    ("x", "tenth"),
    [
        (0.25, 0.3),
        (-0.25, -0.3),
        (0.35, 0.4),
        (-0.35, -0.4),
        (0.65, 0.7),
        (-0.65, -0.7),
        (0.35 - 1e-15, 0.4),  # 0.35 but for floating-point noise
        (0.2499, 0.2),
        (-0.2499, -0.2),
        (0.6678, 0.7),
        (-0.04, 0.0),
    ],
)
def test_nearest_tenth_rounds_ties_away_from_zero(x, tenth):
    got = _nearest_tenth(x)
    assert got == tenth and math.copysign(1.0, got) == (1.0 if tenth >= 0 else -1.0)


@pytest.mark.parametrize("sign", [1, -1])
@pytest.mark.parametrize(("tie", "used"), [(0.25, 0.3), (0.35, 0.4), (0.65, 0.7)])
def test_hand_mode_rounds_the_weighted_skew_away_from_zero(monkeypatch, sign, tie, used):
    monkeypatch.setattr("flowfreq.bulletin17b.weighted_skew", lambda *a: sign * tie)
    assert ex1(hand_computation=True).skew_used == sign * used
    assert ex1().skew_used == sign * tie


def test_hand_mode_leaves_a_station_skew_used_alone_unrounded():
    r = record("example1")
    res = Bulletin17B(r.peaks_cfs, r.water_years, hand_computation=True).run_analysis()
    assert res.skew_used == res.skew_station == pytest.approx(0.7300, abs=2e-4)


@pytest.mark.parametrize("aeps", [[0.01, 0.01], [0.667, 0.67], [0.5, 0.43, 0.1, 0.429]])
def test_duplicate_aeps_refuse(aeps):
    r = record("example1")
    with pytest.raises(Bulletin17BError, match="same AEP"):
        Bulletin17B(r.peaks_cfs, r.water_years, aeps=aeps)


def test_bulletin17c_aeps_are_not_duplicates():
    """0.995 and 0.99 are 0.5 percent apart, but their non-exceedances differ by half."""
    r = record("example1")
    aeps = [float(a) for a in FloodFrequencyAnalysis.STANDARD_AEP]
    assert 0.995 in aeps and 0.99 in aeps
    res = Bulletin17B(r.peaks_cfs, r.water_years, aeps=aeps).run_analysis()
    assert list(res.quantiles.aep) == aeps


def test_outlier_thresholds_by_hand():
    """Eqs. 7 and 8a at N = 24: K_N = 2.467."""
    res = ex1()
    mean, std, _ = log_moments(np.log10(record("example1").peaks_cfs))
    assert res.high_outlier_threshold == pytest.approx(10 ** (mean + 2.467 * std))
    assert res.low_outlier_threshold == pytest.approx(10 ** (mean - 2.467 * std))
    assert mean == pytest.approx(3.3684, abs=1e-4) and std == pytest.approx(0.2456, abs=1e-4)
    # 10 ** (3.3684 +/- 2.467 (0.2456)) = 10 ** 3.9743 and 10 ** 2.7625
    assert res.high_outlier_threshold == pytest.approx(9425, rel=1e-3)
    assert res.low_outlier_threshold == pytest.approx(578.8, rel=1e-3)


def test_confidence_80_uses_the_90_percent_deviate():
    """Two-sided 80 percent is one-sided 0.90, z = 1.28155 (eqs. 9-3 to 9-6 by hand)."""
    res = ex1(confidence=0.80)
    row = res.quantiles.set_index("aep").loc[0.01]
    z, n, k = 1.2815515655, 24, row.K_factor
    a = 1 - z**2 / (2 * (n - 1))
    b = k**2 - z**2 / n
    root = math.sqrt(k**2 - a * b)
    assert row.upper == pytest.approx(10 ** (res.mean_log + res.std_log * (k + root) / a), rel=1e-6)
    assert row.lower == pytest.approx(10 ** (res.mean_log + res.std_log * (k - root) / a), rel=1e-6)
    wide = ex1().quantiles.set_index("aep").loc[0.01]
    assert wide.lower < row.lower < row.flow_cfs < row.upper < wide.upper
