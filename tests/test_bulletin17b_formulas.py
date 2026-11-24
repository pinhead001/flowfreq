"""Bulletin 17B formulas against the bulletin's own printed values."""

from __future__ import annotations

import csv

import numpy as np
import pytest
from scipy.special import ndtri

from flowfreq.bulletin17b import (
    KN_TABLE,
    Bulletin17BError,
    confidence_deviates,
    expected_probability,
    historic_moments,
    log_moments,
    outlier_k,
    station_skew_mse,
    synthetic_statistics,
    weighted_skew,
)
from flowfreq.core import kfactor
from tests.fixtures.b17b import FIXTURES, load, record


def test_outlier_k_reads_appendix_4():
    assert outlier_k(24) == 2.467 and outlier_k(39) == 2.671 and outlier_k(149) == 3.148
    with pytest.raises(Bulletin17BError):
        outlier_k(9)
    with pytest.raises(Bulletin17BError):
        outlier_k(150)


def test_kn_table_matches_the_fixture():
    with open(FIXTURES / "kn_table.csv", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))[1:]
    assert dict((int(n), float(k)) for n, k in rows) == dict(KN_TABLE)


def test_station_skew_mse_eq6():  # Ex. 1 p.12-8; Ex. 2-4 Table 1 values; App. 6
    assert station_skew_mse(0.7300, 24) == pytest.approx(0.277, abs=5e-4)
    assert station_skew_mse(0.1650, 82) == pytest.approx(0.073, abs=5e-4)
    assert station_skew_mse(0.5948, 38) == pytest.approx(0.183, abs=5e-4)
    assert station_skew_mse(-0.529, 42) == pytest.approx(0.163, abs=5e-4)
    assert station_skew_mse(0.0418, 77) == pytest.approx(0.07074, abs=5e-6)


def test_weighted_skew_eq5():
    assert weighted_skew(0.7300, 0.277, 0.6, 0.302) == pytest.approx(0.6678, abs=1e-4)
    assert weighted_skew(0.1650, 0.073, -0.3, 0.302) == pytest.approx(0.0745, abs=1e-4)


def test_log_moments_example1():  # p.12-7
    m, s, g = log_moments(np.log10(record("example1").peaks_cfs))
    assert (m, s, g) == (
        pytest.approx(3.3684, abs=1e-4),
        pytest.approx(0.2456, abs=1e-4),
        pytest.approx(0.7300, abs=2e-4),
    )


def test_log_moments_refuses_degenerate_input():
    with pytest.raises(Bulletin17BError):
        log_moments(np.array([3.0, 3.1]))
    with pytest.raises(Bulletin17BError):
        log_moments(np.array([3.0, 3.0, 3.0]))


def test_historic_moments_big_sandy():  # Figure 6-2
    bs = record("appendix6_big_sandy")
    hm = historic_moments(
        np.log10(bs.peaks_cfs), np.log10([q for _, q in bs.historic_peaks]), h=77, l=0
    )
    assert hm.w == pytest.approx(1.68182, abs=1e-5)
    assert (hm.mean, hm.std, hm.skew) == (
        pytest.approx(3.71581, abs=1e-5),
        pytest.approx(0.28898, abs=1e-5),
        # The bulletin prints 0.0418, but its own line, (77)(0.07485) / [(76)(75)(0.02413)],
        # evaluates to 0.04190, and full precision gives 0.04191, so the printed value is
        # a source slip. 2e-4 keeps it in range.
        pytest.approx(0.0418, abs=2e-4),
    )


def test_synthetic_statistics_examples():  # eqs. 12-29, 12-33
    assert synthetic_statistics(23880, 11210, 5230)[2] == pytest.approx(0.5948, abs=1e-4)
    assert synthetic_statistics(17940, 6000, 1060)[2] == pytest.approx(-0.5287, abs=1e-4)
    with pytest.raises(Bulletin17BError):
        synthetic_statistics(1e6, 1000, 999)


def test_synthetic_statistics_k_skew_overrides_the_deviates():
    default = synthetic_statistics(23880, 11210, 5230)
    rounded = synthetic_statistics(23880, 11210, 5230, k_skew=0.6)
    assert rounded[2] == default[2]
    assert rounded[1] != default[1]


def test_confidence_deviates_appendix9_check():  # pp. 9-3/9-4
    ku, kl = confidence_deviates(np.array([kfactor(0.2, 0.01)]), 50, 1.645)
    assert ku[0] == pytest.approx(3.026, abs=1e-3) and kl[0] == pytest.approx(2.059, abs=1e-3)


def test_expected_probability_appendix11_check():  # p.11-4: Table 11-1 gives 0.0174
    assert expected_probability(np.array([0.01]), 20)[0] == pytest.approx(0.0174, rel=0.02)


def _table_9_1_by_n():
    by_n = {}
    with open(FIXTURES / "table9_1_confidence_deviates.csv", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        pcols = [c for c in reader.fieldnames if c.startswith("P.")]
        for row in reader:
            n = int(row["N"])
            if n == 10:
                # The a/b approximation of eqs. 9-5 and 9-6 is off by up to 0.36 at N = 10,
                # while Table 9-1 is the exact non-central t: not a check on the formula.
                continue
            level = float(row["confidence_level"])
            by_n.setdefault(n, []).extend((level, float(c[1:]), float(row[c])) for c in pcols)
    return by_n


#: Largest error of the a/b approximation against Table 9-1 at each N, measured over all
#: eight levels and twelve probabilities, rounded up. Table 9-1 is the exact non-central
#: t. Eqs. 9-4 to 9-6 approximate it, and the error shrinks as N grows. Within the table
#: it is 0.005 or less at no N: the worst cell at N = 100 is 0.009 (level 0.25).
_APPROXIMATION_TOLERANCE = {
    15: 0.14,
    20: 0.075,
    25: 0.05,
    30: 0.04,
    40: 0.03,
    50: 0.02,
    60: 0.02,
    70: 0.015,
    80: 0.012,
    90: 0.011,
    100: 0.010,
}


@pytest.mark.parametrize("n", sorted(_table_9_1_by_n()))
def test_confidence_deviates_against_table_9_1(n):
    for level, p, printed in _table_9_1_by_n()[n]:
        k = np.array([ndtri(1 - p)])
        z = float(ndtri(max(level, 1 - level)))
        ku, kl = confidence_deviates(k, n, z)
        got = ku[0] if level < 0.5 else kl[0]
        assert got == pytest.approx(printed, abs=_APPROXIMATION_TOLERANCE[n]), (level, p)


def test_confidence_deviates_example1_table_12_4():
    data = load("example1")["confidence_limits_table_12_4"]
    z = 1.64485  # z_c printed on p. 12-12
    for row in data["rows"]:
        ku, kl = confidence_deviates(np.array([float(row["K_Gw_P"])]), data["N_used"], z)
        assert ku[0] == pytest.approx(float(row["K_U"]), abs=2e-3), row["P"]
        assert kl[0] == pytest.approx(float(row["K_L"]), abs=2e-3), row["P"]


def test_confidence_deviates_refuses_short_records():
    k = np.array([1.0])
    for n in (1, 2):
        with pytest.raises(Bulletin17BError):
            confidence_deviates(k, n, 1.645)
    # N = 3 works at 95 percent but eq. 9-5 gives a < 0 at 99 percent.
    ku, kl = confidence_deviates(k, 3, 1.645)
    assert np.all(np.isfinite(ku)) and np.all(np.isfinite(kl)) and kl[0] < ku[0]
    with pytest.raises(Bulletin17BError):
        confidence_deviates(k, 3, float(ndtri(0.99)))


@pytest.mark.parametrize("flow", [40.0, 70.0, 500.0, 12340.0])
@pytest.mark.parametrize("length", [3, 7, 10, 24, 40])
def test_log_moments_refuses_a_constant_record(flow, length):
    with pytest.raises(Bulletin17BError):
        log_moments(np.log10(np.full(length, flow)))
    with pytest.raises(Bulletin17BError):
        historic_moments(np.log10(np.full(length, flow)), np.log10([flow]), h=length + 10, l=0)


def test_moments_refuse_non_finite_logs():
    with pytest.raises(Bulletin17BError):
        log_moments(np.array([3.0, 3.2, -np.inf, 3.5]))
    with pytest.raises(Bulletin17BError):
        log_moments(np.array([3.0, 3.2, np.nan, 3.5]))
    with pytest.raises(Bulletin17BError):
        historic_moments(np.array([3.0, 3.2, 3.5]), np.array([-np.inf]), h=10, l=0)


def test_historic_moments_refuses_empty_record_and_short_period():
    with pytest.raises(Bulletin17BError):
        historic_moments(np.array([]), np.array([4.0]), h=10, l=0)
    # H = 5 is shorter than 3 systematic + 4 historic years: W < 1.
    with pytest.raises(Bulletin17BError):
        historic_moments(np.array([3.0, 3.2, 3.5]), np.array([4.0] * 4), h=5, l=0)


def test_synthetic_statistics_printed_std_and_mean():  # eqs. 12-30/12-31, 12-34/12-35
    # Example 3 takes K.01 and K.50 at the rounded skew 0.6, Example 4 at -0.5.
    mean3, std3, _ = synthetic_statistics(23880, 11210, 5230, k_skew=0.6)
    assert std3 == pytest.approx(0.2310, abs=5e-5) and mean3 == pytest.approx(3.7415, abs=5e-5)
    mean4, std4, _ = synthetic_statistics(17940, 6000, 1060, k_skew=-0.5)
    assert std4 == pytest.approx(0.6564, abs=5e-5) and mean4 == pytest.approx(2.9708, abs=5e-5)


def test_station_skew_mse_against_every_table_1_cell():
    with open(FIXTURES / "table1_skew_mse.csv", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 31
    for row in rows:
        g = float(row["station_skew_abs"])
        for col, printed in row.items():
            if col == "station_skew_abs":
                continue
            n = int(col[1:])
            if round(g, 1) == 0.9:
                # Table 1's 0.9 row was computed with the |G| > 0.90 branch of A
                # (-0.52 + 0.30|G|); eq. 6 as printed uses -0.33 + 0.08|G| there.
                a = -0.52 + 0.30 * g
                b = 0.94 - 0.26 * g
                expected = 10 ** (a - b * np.log10(n / 10))
                assert float(printed) == pytest.approx(expected, abs=1.5e-3), (g, n)
            else:
                assert station_skew_mse(g, n) == pytest.approx(float(printed), abs=1.5e-3), (g, n)


def test_historic_moments_with_low_values_against_eqs_6_2b_to_6_4b():
    # l > 0: two low values were left out of the systematic record but count in it.
    rng = np.random.default_rng(7)
    xs = rng.normal(3.5, 0.3, 30)
    xz = np.array([4.3, 4.1])
    h, l = 60, 2
    hm = historic_moments(xs, xz, h=h, l=l)
    n = xs.size
    m_sys, s_sys, g_sys = log_moments(xs)
    w = (h - xz.size) / (n + l)
    denom = h - w * l
    m_t = (w * n * m_sys + xz.sum()) / denom  # eq. 6-2b
    d = m_sys - m_t
    sum2 = (n - 1) * s_sys**2 + n * d**2
    sum3 = g_sys * (n - 1) * (n - 2) * s_sys**3 / n + 3 * d * (n - 1) * s_sys**2 + n * d**3
    s2 = (w * sum2 + np.sum((xz - m_t) ** 2)) / (denom - 1)  # eq. 6-3b
    g_t = (
        denom / ((denom - 1) * (denom - 2)) * (w * sum3 + np.sum((xz - m_t) ** 3)) / s2**1.5
    )  # eq. 6-4b
    assert hm.w == pytest.approx(w)
    assert hm.mean == pytest.approx(m_t)
    assert hm.std == pytest.approx(s2**0.5)
    assert hm.skew == pytest.approx(g_t)
