"""``FrequencyResults.n_low_outliers`` is peakfq's ``gbnlow``, on every WY/MT station.

``gbtest`` (``emafit.f`` lines 1062-1076) counts every EMA row whose upper
bound lies below the low-outlier cutoff: exact peaks, zero-flow years,
interval and historic rows, and years with no peak that a perception
threshold below the cutoff censors. The native engine used to report MGBT's
own count instead, which misses the last kind: 17 against peakfq's 20 at
06328100, 12 against 15 at 06326960, 6 against 8 at 06177820. MGBT's count
is still reported, as ``n_mgbt_outliers``, and the Fortran engine derives the
same number from ``gbtest``'s sample.

The per-station check stops short of the EMA fit, which does not affect
either count (``run_analysis`` computes both before fitting), so all 24
stations cost seconds; three full fits confirm the public field.
"""

from __future__ import annotations

import pytest

from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE
from tests.fixtures.paths import testdata_path as _testdata

pytest.importorskip(
    "flowfreq.peakfqr",
    reason="Fortran extension not built; run python build_fortran/build.py "
    "(see docs/FORTRAN_UPLOAD.md)",
    # See tests/fortran_parity/test_fortran_oracles.py's exc_type comment.
    exc_type=ImportError,
)

pytestmark = [
    pytest.mark.requires_fortran,
    pytest.mark.requires_peakfqr_testdata,
    pytest.mark.skipif(not TESTDATA_AVAILABLE, reason=SKIP_REASON),
]

#: gbnlow on the stations where it differs from MGBT's own count.
GAP_YEARS_BELOW_CUTOFF = {"06328100.00": (17, 20), "06326960.00": (12, 15), "06177820.00": (6, 8)}


@pytest.fixture(scope="module")
def stations():
    from flowfreq.psf_convert import convert_psf

    return convert_psf(_testdata("wymt_ffa_2022A.psf"))


def _native_counts(station):
    """(n_mgbt_outliers, n_low_outliers) as ``run_analysis`` computes them, without the fit."""
    from flowfreq.bulletin17c import ExpectedMomentsAlgorithm

    kw = station.bulletin17c_kwargs("native")
    ema = ExpectedMomentsAlgorithm(
        kw["peak_flows"],
        water_years=kw["water_years"],
        regional_skew=kw["regional_skew"],
        regional_skew_mse=kw["regional_skew_mse"],
        historical_peaks=kw["historical_peaks"],
        perception_thresholds=kw["perception_thresholds"],
        user_low_outlier_threshold=kw["user_low_outlier_threshold"],
        interval_peaks=kw["interval_peaks"],
    )
    threshold, flagged, _pilf = ema._multiple_grubbs_beck(
        user_threshold=kw["user_low_outlier_threshold"]
    )
    ema._build_flow_intervals(threshold)
    return flagged, ema._n_rows_below_cutoff


def test_every_station_counts_gbnlow(stations):
    from flowfreq.fortran_engine import _mgbt_flagged

    assert len(stations) == 24
    for sid, station in stations.items():
        reference = station.fortran_reference()
        flagged, rows = _native_counts(station)
        assert rows == reference.low_outlier_count, sid
        expected_flagged = _mgbt_flagged(station.emafit_arrays(), reference.low_outlier_threshold)
        assert flagged == expected_flagged, sid
        if sid in GAP_YEARS_BELOW_CUTOFF:
            assert (flagged, rows) == GAP_YEARS_BELOW_CUTOFF[sid], sid
        else:
            assert flagged == rows, sid


@pytest.mark.parametrize("sid", sorted(GAP_YEARS_BELOW_CUTOFF))
def test_public_fields_on_both_engines(stations, sid):
    n_mgbt, n_rows = GAP_YEARS_BELOW_CUTOFF[sid]
    native = stations[sid].run("native")
    fortran = stations[sid].run("fortran")
    assert native.n_low_outliers == fortran.n_low_outliers == n_rows
    assert native.n_mgbt_outliers == fortran.n_mgbt_outliers == n_mgbt
    assert len(native.pilf_flows) == n_mgbt
