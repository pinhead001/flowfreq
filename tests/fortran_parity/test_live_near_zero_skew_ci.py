"""Confidence bounds at near-zero skew: native against live ``emafitpr``.

``emafit.f`` section 4.2 (lines 822-850): when |weighted skew| <= skewmin
(0.06324555) the Fortran runs ``var_emab`` at -skewmin and +skewmin and
interpolates the bounds linearly in skew, rather than evaluating them at the
near-zero skew itself. The native engine evaluated at the skew directly and
was 0.77-0.96% off on these records; with the interpolation it is ~0.002%.

Synthetic, uncensored records with a tight regional skew near zero, so the
weighted skew lands inside the interpolation zone. File name starts
``test_live_`` so it sorts after ``test_fortran_oracles.py``; see
``test_live_interval_builder.py``'s docstring for why that matters.
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

AEPS = [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002]


@pytest.mark.parametrize("seed, regional_skew", [(1, 0.0), (2, 0.03), (5, -0.04)])
def test_near_zero_skew_bounds_match_live_fortran(seed, regional_skew):
    from flowfreq.bulletin17c import _CI_SKEWMIN, Bulletin17C
    from flowfreq.fortran_engine import run_fortran_reference

    rng = np.random.default_rng(seed)
    flows = 10 ** rng.normal(3.0, 0.25, 50)
    years = np.arange(1970, 2020)
    kwargs = dict(regional_skew=regional_skew, regional_skew_mse=0.02)

    native = Bulletin17C(flows, years, **kwargs)
    native.run_analysis(method="ema")
    assert abs(native.results.skew_used) < _CI_SKEWMIN  # inside the interpolation zone
    ci = native.compute_confidence_limits(np.array(AEPS))

    ref, _ = run_fortran_reference(flows, water_years=years, aeps=AEPS, **kwargs)
    for i, aep in enumerate(AEPS):
        lo, hi = ref.confidence_intervals[aep]
        assert ci["lower_5pct"].iloc[i] == pytest.approx(lo, rel=1e-4), aep
        assert ci["upper_5pct"].iloc[i] == pytest.approx(hi, rel=1e-4), aep
