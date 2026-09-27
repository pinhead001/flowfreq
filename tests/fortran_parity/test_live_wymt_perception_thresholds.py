"""Native EMA against live ``emafitpr`` where perception thresholds reach into the record.

The WY/MT stations of ``wymt_thresholds.py`` have threshold periods that
overlap or lie inside the systematic record, or several separate periods.
``ExpectedMomentsAlgorithm`` used to drop the first kind and collapse the
second, so its fit on these stations was of different data than peakfq's.
With every period applied per year (``test_wymt_perception_rows.py`` checks
the rows), the fits agree -- compared here with a live ``emafitpr`` call on
the same inputs, through ``flowfreq.fortran_engine.run_fortran_reference``.

Tolerances are ``test_live_interval_builder.py``'s (the same calibration as
``test_live_vs_golden.py``): the EMA fixed point is ill-conditioned in the
third moment, and censored stations carry a few 1e-6 of cross-build noise.
Measured on Windows/MSYS2 gfortran, the largest differences on the passing
stations are 9e-7 on the mean, 2e-6 on the standard deviation and 3e-6 on
either skew.

06324500.00/.01 (17 MGBT low outliers each) were ``xfail(strict=True)`` on the
weighted fit until the native at-site skew MSE followed ``emafit.f:707``'s B17B
switch; they now match like the rest.

File name starts ``test_live_`` so it sorts after ``test_fortran_oracles.py``;
see ``test_live_interval_builder.py``'s docstring for why that matters.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE
from tests.fortran_parity.wymt_thresholds import THRESHOLDS, bulletin17c_inputs

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

# test_live_interval_builder.py's calibration.
ATOL_SKEW = 1e-3
ATOL_MEAN = 1e-4
ATOL_VARIANCE = 1e-5


def _weighted_param(site_no):
    return pytest.param(site_no)


@pytest.fixture(scope="module")
def fits():
    """site_no -> (native FrequencyResults, live ReferenceResult)."""
    from flowfreq.bulletin17c import Bulletin17C
    from flowfreq.fortran_engine import run_fortran_reference

    out = {}
    for site_no in sorted(THRESHOLDS):
        kwargs = bulletin17c_inputs(site_no)
        b17c = Bulletin17C(**kwargs)
        results = b17c.run_analysis(method="ema")
        reference, _ = run_fortran_reference(**kwargs)
        out[site_no] = (results, reference)
    return out


@pytest.mark.parametrize("site_no", sorted(THRESHOLDS))
def test_same_number_of_rows(fits, site_no):
    results, reference = fits[site_no]
    assert results.n_peaks == reference.n_peaks


@pytest.mark.parametrize("site_no", sorted(THRESHOLDS))
def test_mgbt_agrees(fits, site_no):
    results, reference = fits[site_no]
    assert results.n_low_outliers == reference.low_outlier_count
    if reference.low_outlier_count:
        assert results.low_outlier_threshold == pytest.approx(
            reference.low_outlier_threshold, rel=1e-9
        )


@pytest.mark.parametrize("site_no", sorted(THRESHOLDS))
def test_at_site_skew(fits, site_no):
    results, reference = fits[site_no]
    assert abs(results.skew_station - reference.parameters["skew_at_site"]) < ATOL_SKEW


@pytest.mark.parametrize("site_no", [_weighted_param(s) for s in sorted(THRESHOLDS)])
def test_final_moments(fits, site_no):
    results, reference = fits[site_no]
    p = reference.parameters
    assert abs(results.mean_log - p["mean_log"]) < ATOL_MEAN
    assert abs(results.std_log**2 - p["std_log"] ** 2) < ATOL_VARIANCE
    assert abs(results.skew_used - p["skew_weighted"]) < ATOL_SKEW


def test_measured_agreement_is_far_inside_the_tolerance(fits):
    """The tolerances above are loose on purpose; the real gap is ~1e-6."""
    worst = 0.0
    for site_no, (results, reference) in fits.items():
        p = reference.parameters
        worst = max(
            worst,
            abs(results.mean_log - p["mean_log"]),
            abs(results.std_log - p["std_log"]),
            abs(results.skew_used - p["skew_weighted"]),
            abs(results.skew_station - p["skew_at_site"]),
        )
    assert worst < 1e-4, worst
    assert np.isfinite(worst)
