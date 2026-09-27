"""Native interval peaks and upper thresholds against a live ``emafitpr`` call.

A synthetic record exercising every new row type at once -- two less-than
peaks below every exact peak (which ``gbtest`` adds to MGBT's sample), a
greater-than peak under a ``(30, 20000)`` perception threshold, and a bounded
interval -- fitted by both engines from the same ``Bulletin17C`` arguments.
Measured: weighted skew 3.7e-7, quantiles 4e-5 %, confidence bounds 0.04 %;
both engines flag the same 4 low outliers at the same 287 cfs cutoff.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip(
    "flowfreq.peakfqr",
    reason="Fortran extension not built; run python build_fortran/build.py "
    "(see docs/FORTRAN_UPLOAD.md)",
    # See tests/fortran_parity/test_fortran_oracles.py's exc_type comment.
    exc_type=ImportError,
)

pytestmark = [pytest.mark.requires_fortran]


def _kwargs():
    rng = np.random.default_rng(7)
    years = np.arange(1960, 2010)
    flows = np.round(10 ** rng.normal(3.0, 0.35, len(years)), 0)
    intervals = {1965: (0.0, 60.0), 1980: (0.0, 50.0), 1990: (2500.0, np.inf), 2000: (400.0, 900.0)}
    keep = ~np.isin(years, list(intervals))
    return dict(
        peak_flows=flows[keep],
        water_years=years[keep],
        interval_peaks=[(y, lo, hi) for y, (lo, hi) in intervals.items()],
        perception_thresholds={(1960, 2009): 0.0, (1985, 1995): (30.0, 20000.0)},
        regional_skew=-0.1,
        regional_skew_mse=0.3,
    )


def test_native_matches_fortran_engine():
    from flowfreq.bulletin17c import Bulletin17C

    kw = _kwargs()
    assert min(kw["peak_flows"]) > 60.0  # both less-than peaks enter MGBT
    native = Bulletin17C(**kw).run_analysis(engine="native")
    fortran = Bulletin17C(**kw).run_analysis(engine="fortran")

    assert native.n_low_outliers == fortran.n_low_outliers == 4
    assert native.low_outlier_threshold == pytest.approx(fortran.low_outlier_threshold)
    assert abs(native.mean_log - fortran.mean_log) < 1e-5
    assert abs(native.skew_used - fortran.skew_used) < 2e-4
    qn = native.quantiles.set_index("aep")["flow_cfs"]
    qf = fortran.quantiles.set_index("aep")["flow_cfs"]
    assert float((qn / qf - 1).abs().max()) < 5e-4
    cn = native.confidence_limits.set_index("aep")
    cf = fortran.confidence_limits.set_index("aep")
    for col in ("lower_5pct", "upper_5pct"):
        assert float((cn[col] / cf[col] - 1).abs().max()) < 2e-3, col
