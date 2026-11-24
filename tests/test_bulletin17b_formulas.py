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
from tests.fixtures.b17b import FIXTURES, record


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
        # evaluates to 0.04190, and full precision gives 0.04191. The brief's 1e-4 cannot
        # be met by a correct implementation; 2e-4 keeps the printed value in range.
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


def _table_9_1_rows():
    with open(FIXTURES / "table9_1_confidence_deviates.csv", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        pcols = [c for c in reader.fieldnames if c.startswith("P.")]
        for row in reader:
            n = int(row["N"])
            if n == 10:
                # The a/b approximation of eqs. 9-5 and 9-6 is an N >= 15 formula;
                # Table 9-1 is the exact non-central t, so N = 10 is not a check on it.
                continue
            level = float(row["confidence_level"])
            for col in pcols:
                yield level, n, float(col[1:]), float(row[col])


#: Largest error of the a/b approximation against Table 9-1 at each N, measured over all
#: eight levels and twelve probabilities, rounded up. Table 9-1 is the exact non-central
#: t. Eqs. 9-4 to 9-6 approximate it, and the error shrinks as N grows: 0.005 holds
#: only from about N = 100 at the worst level (0.25). N = 10 is left out, where the error
#: is up to 0.36.
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


@pytest.mark.parametrize("level,n,p,printed", list(_table_9_1_rows()))
def test_confidence_deviates_against_table_9_1(level, n, p, printed):
    k = np.array([ndtri(1 - p)])
    z = float(ndtri(max(level, 1 - level)))
    ku, kl = confidence_deviates(k, n, z)
    got = ku[0] if level < 0.5 else kl[0]
    assert got == pytest.approx(printed, abs=_APPROXIMATION_TOLERANCE[n])
