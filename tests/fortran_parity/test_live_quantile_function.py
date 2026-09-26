"""Native quantiles against live ``emafitpr`` wherever the moments already agree.

``emafitpr`` reports ``yp = qP3(pq, cmoms(:,1))`` -- the exact Pearson III
inverse of the *weighted* moments (``emafit.f`` line 848). Where the native
fit reproduces those moments, the native quantiles must reproduce ``yp`` too,
at every AEP. They did not while ``flowfreq.core.kfactor`` used the
Wilson-Hilferty approximation: up to 1.35% at AEP 0.002 on 06185500.11
(skew 0.87), 0.96% on 06329200.00 (skew -0.44), 0.10% on Powder River
(skew -0.18).

These stations are uncensored, with no MGBT low outliers, so the moments
agree to 1e-9 or better and the EMA's ill-conditioning (see
``test_live_vs_golden.py``) does not arise. The largest measured gap is
4.2e-8 in log10 (AEP 0.995 on 06329200.00), and it is not the native side's:
``q_p3`` and the Fortran's own ``qp3sub`` differ by exactly that on
*identical* moments -- the accuracy of ``emafit.f``'s ``DCHIIN`` chi-square
inverse in the far lower tail of a negatively skewed fit. Elsewhere the gap is
1e-9 or less. The bound below is 1e-7 in log10.

File name starts ``test_live_`` so it sorts after ``test_fortran_oracles.py``;
see ``test_live_interval_builder.py``'s docstring for why that matters.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE

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

AEPS = (0.995, 0.99, 0.95, 0.9, 0.8, 0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)

#: site_no -> the skew each fit carries, to show the range is covered.
STATIONS = {
    "06185500.11": 0.871,  # station skew
    "06329200.00": -0.445,  # weighted
    "06326500.00": -0.184,  # weighted (Powder River)
}


@pytest.fixture(scope="module", params=sorted(STATIONS))
def fit(request):
    from flowfreq.bulletin17c import Bulletin17C
    from flowfreq.fortran_engine import run_fortran_reference
    from tests.fixtures.wymt_peaks import load_site

    site = load_site(request.param)
    years = sorted(site.peaks)
    weighted = site.regional_skew is not None
    kwargs = dict(
        peak_flows=[site.peaks[y] for y in years],
        water_years=years,
        regional_skew=site.regional_skew if weighted else None,
        regional_skew_mse=site.regional_skew_mse if weighted else None,
    )
    b17c = Bulletin17C(**kwargs)
    b17c.run_analysis(method="ema")
    reference, _ = run_fortran_reference(**kwargs, aeps=AEPS)
    return request.param, b17c, reference


def test_moments_agree(fit):
    site_no, b17c, reference = fit
    r, p = b17c.results, reference.parameters
    assert abs(r.mean_log - p["mean_log"]) < 1e-9
    assert abs(r.std_log - p["std_log"]) < 1e-9
    assert abs(r.skew_used - p["skew_weighted"]) < 1e-9
    assert r.skew_used == pytest.approx(STATIONS[site_no], abs=1e-3)


def test_quantiles_agree_at_every_aep(fit):
    site_no, b17c, reference = fit
    q = b17c.compute_quantiles(aep=np.array(AEPS))
    for aep, log_flow in zip(AEPS, q["log_flow"]):
        expected = np.log10(reference.quantiles[aep])
        assert abs(float(log_flow) - expected) < 1e-7, (site_no, aep)


def test_residual_is_the_fortran_gamma_inverse(fit):
    """Native quantiles track ``qp3sub`` as closely as ``q_p3`` itself does on the same moments."""
    from flowfreq._p3_moments import q_p3
    from flowfreq.peakfqr._emafort import qp3sub

    _, b17c, reference = fit
    p = reference.parameters
    mc = np.array([p["mean_log"], p["std_log"] ** 2, p["skew_weighted"]])
    q = b17c.compute_quantiles(aep=np.array(AEPS))
    for aep, log_flow in zip(AEPS, q["log_flow"]):
        port_gap = abs(q_p3(1 - aep, mc) - float(qp3sub(1 - aep, mc)))
        native_gap = abs(float(log_flow) - np.log10(reference.quantiles[aep]))
        assert native_gap < port_gap + 1e-9, aep
