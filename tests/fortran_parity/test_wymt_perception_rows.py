"""Native EMA rows match ``siteQT`` on the WY/MT stations with in-record thresholds.

``ExpectedMomentsAlgorithm._build_flow_intervals`` against
``flowfreq.fortran_engine.build_emafit_arrays`` -- the library's own ``siteQT``
translation, the one ``emafitpr`` is fed from -- on real stations whose
perception thresholds overlap or sit inside the systematic record (see
``wymt_thresholds.py``). Runs without the Fortran extension. The fits
themselves are compared with a live ``emafitpr`` call in
``test_live_wymt_perception_thresholds.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE
from tests.fortran_parity.wymt_thresholds import THRESHOLDS, bulletin17c_inputs

pytestmark = [
    pytest.mark.requires_peakfqr_testdata,
    pytest.mark.skipif(not TESTDATA_AVAILABLE, reason=SKIP_REASON),
]

QMIN = 1e-20


def _rows_native(kwargs):
    from flowfreq.bulletin17c import ExpectedMomentsAlgorithm

    ema = ExpectedMomentsAlgorithm(**kwargs)
    return sorted(
        (
            int(iv.year),
            round(float(np.log10(max(iv.lower, QMIN))), 10),
            round(float(np.log10(max(iv.upper, QMIN))), 10),
            round(float(np.log10(max(iv.perception_threshold, QMIN))), 10),
            int(iv.is_historical),
        )
        for iv in ema._build_flow_intervals(0.0)
    )


def _rows_sitqt(kwargs):
    from flowfreq.fortran_engine import build_emafit_arrays

    arrays = build_emafit_arrays(
        kwargs["peak_flows"],
        water_years=kwargs["water_years"],
        historical_peaks=kwargs["historical_peaks"],
        perception_thresholds=kwargs["perception_thresholds"],
    )
    return sorted(
        (int(y), round(float(a), 10), round(float(b), 10), round(float(c), 10), int(d))
        for y, a, b, c, d in zip(arrays.years, arrays.ql, arrays.qu, arrays.tl, arrays.dtype)
    )


@pytest.mark.parametrize("site_no", sorted(THRESHOLDS))
def test_native_rows_match_sitqt(site_no):
    kwargs = bulletin17c_inputs(site_no)
    assert _rows_native(kwargs) == _rows_sitqt(kwargs)


@pytest.mark.parametrize("site_no", sorted(THRESHOLDS))
def test_every_nonzero_period_year_has_a_row(site_no):
    kwargs = bulletin17c_inputs(site_no)
    years = {row[0] for row in _rows_native(kwargs)}
    for (start, end), lower in THRESHOLDS[site_no].items():
        if lower > 0:
            assert set(range(start, end + 1)) <= years, (start, end)
