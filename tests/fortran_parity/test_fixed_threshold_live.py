"""A user-supplied (FIXED) low-outlier threshold against a live ``emafitpr`` call.

The record: 50 lognormal peaks, a fixed 300 cfs threshold censoring 5 of
them, regional skew -0.1 with MSE 0.3. The at-site skew is 0.056 and the
regional weighting pulls it to about -0.006. Reported by the interval-peaks
work (PR #72) as a 2e-3 weighted-skew / 0.17 % quantile residual.

What these tests pin down:

* Everything up to the regional weighting matches: the at-site moments, the
  FIXED branch of ``gbtest`` (5 rows censored at 300 cfs, ADJE kept), and
  ``nG`` to 5e-5 (``Wd`` and ADJE's ``as_G_mse`` are within 3.4e-5 and 1.4e-5).
* The weighted skew does not, and cannot be made to: at a skew this close
  to zero ``emafitpr``'s own answer is not a smooth function of its inputs.
  ``mP3`` evaluates a censored interval's moments through ``DGAMDF``
  (Numerical Recipes ``GSER``/``GCF`` with a 6-term Lanczos log-gamma, error
  ~2e-10) at alpha = 4/skew**2 ~ 1e5, then expands about tau ~ -100; the
  CDF's rounding comes out as a 1e-4 to 1e-3 relative error in
  ``E[X^3 | censored]`` that jumps around as the skew moves. Perturbing
  ``r_G_mse`` by one part in 10**6 moves ``emafitpr``'s weighted skew by
  3e-4, by 1e-5 moves it 2e-3; its answers over such perturbations span
  -0.0073 to -0.0003 (the 1 % AEP flood 4241-4256 cfs). The native fit
  (-0.0037, exact incomplete gamma) sits inside that band. A native port of
  ``DGAMDF`` was tried: it reproduces ``mp3`` to 1e-13 wherever the
  Fortran is well conditioned, and here merely lands on a different point
  of the same band, because ``nG`` differs from the Fortran's by 5e-5.
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

THRESHOLD = 300.0


def _record():
    rng = np.random.default_rng(7)
    years = np.arange(1960, 2010)
    flows = np.round(10 ** rng.normal(3.0, 0.35, len(years)), 0)
    return flows, years


def _kwargs():
    flows, years = _record()
    return dict(
        peak_flows=flows,
        water_years=years,
        regional_skew=-0.1,
        regional_skew_mse=0.3,
        user_low_outlier_threshold=THRESHOLD,
    )


@pytest.fixture(scope="module")
def fits():
    from flowfreq.bulletin17c import Bulletin17C

    native = Bulletin17C(**_kwargs())
    native.run_analysis(engine="native")
    fortran = Bulletin17C(**_kwargs())
    fortran.run_analysis(engine="fortran")
    return native, fortran


def test_fixed_branch_and_at_site_fit_match(fits):
    native, fortran = fits
    n, f = native.results, fortran.results
    assert n.low_outlier_threshold == pytest.approx(THRESHOLD)
    assert f.low_outlier_threshold == pytest.approx(THRESHOLD)
    assert n.n_low_outliers == f.n_low_outliers == 5
    # A user threshold keeps ADJE (emafit.f:707: B17B only when MGBT computed it).
    assert native._analyzer._at_site_option == "ADJE"
    assert abs(n.skew_station - f.skew_station) < 1e-7


def test_regional_weighting_inputs_match(fits):
    """``nG = n * Wd * as_G_mse / r_G_mse`` agrees with ``emafitpr``'s to 1e-4."""
    from flowfreq.fortran_engine import run_fortran_reference

    native, _ = fits
    flows, years = _record()
    reference, _arrays = run_fortran_reference(
        flows,
        years,
        user_low_outlier_threshold=THRESHOLD,
        regional_skew=-0.1,
        regional_skew_mse=0.3,
    )
    p = reference.parameters
    ema = native._analyzer
    nobs, tl, tu = ema._perception_threshold_groups()
    r = native.results
    # The at-site moments the Fortran weights with are cmoms(:,2).
    at_site = (p["mean_log_at_site"], p["std_log_at_site"], p["skew_at_site"])
    n_rows = len(ema.intervals)
    wd = ema._detrat_wd(
        tuple(nobs),
        tuple(tl),
        tuple(tu),
        float(at_site[0]),
        float(at_site[1]) ** 2,
        float(at_site[2]),
        n_rows,
    )
    as_g_mse = ema._at_site_skew_mse(at_site[0], at_site[1], at_site[2], n_rows)
    assert wd == pytest.approx(p["weight_factor"], rel=1e-4)
    assert as_g_mse == pytest.approx(p["mse_skew"], rel=3e-5)
    assert r.skew_station == pytest.approx(at_site[2], abs=1e-7)


def test_emafitpr_weighted_skew_is_ill_conditioned_here():
    """The root cause: one part in 10**5 of ``r_G_mse`` moves peakfq's own skew by >1e-3."""
    from flowfreq.fortran_engine import run_fortran_reference

    flows, years = _record()
    skews = []
    for mse in (0.3, 0.3 * (1.0 + 1e-5)):
        reference, _ = run_fortran_reference(
            flows,
            years,
            user_low_outlier_threshold=THRESHOLD,
            regional_skew=-0.1,
            regional_skew_mse=mse,
            aeps=[0.01],
        )
        skews.append(reference.parameters["skew_weighted"])
    assert abs(skews[1] - skews[0]) > 1e-3


def test_mp3_oracle_is_rough_at_near_zero_skew():
    """``mP3``'s own rounding at skew 0.003: 1e-3 relative in E[X^3], vs 1e-10 at skew 0.06."""
    from flowfreq._p3_moments import m_p3
    from flowfreq.peakfqr import _emafort

    tl, tu = -6.0, np.log10(THRESHOLD)

    def rel_err(skew):
        m = np.array([2.8959, 0.0993, skew])
        mine = m_p3(tl, tu, m, 3)
        theirs = np.asarray(_emafort.mp3(tl, tu, m, 3))[:3]
        return abs(mine[2] / theirs[2] - 1.0)

    assert rel_err(0.0564) < 1e-9
    assert rel_err(0.003) > 1e-4


@pytest.mark.xfail(
    strict=True,
    reason="emafitpr's weighted skew is ill-conditioned on this record (see the module "
    "docstring): native -0.0037 vs Fortran -0.0056, inside the Fortran's own "
    "1e-5-perturbation band of -0.0073..-0.0003",
)
def test_weighted_skew_matches_emafitpr(fits):
    native, fortran = fits
    assert abs(native.results.skew_used - fortran.results.skew_used) < 1e-5


# --------------------------------------------------------------------------- #
# Big Sandy's systematic record with a 6000 cfs FIXED threshold: 29 of 44
# peaks censored. Until the fixed point followed p3est_ema's iteration limit
# this never reached the regional weighting natively (100 iterations were not
# enough; run_ffa fell back to MOM). It converges now, and exposes a
# different, pre-existing divergence: emafitpr's ADJE as_G_mse is 2.74 where
# the native one is 0.064, because the Fortran MN2MVARB stops after its 100
# Newton iterations without converging and returns that iterate unflagged --
# mc2mnvb of its answer misses the target s_mn by up to 0.08, the native
# root by 3e-15. The at-site fit agrees to 5e-8.
# --------------------------------------------------------------------------- #


def _big_sandy_6000():
    from tests.fixtures.big_sandy import SYSTEMATIC_PEAKS

    years = np.array(sorted(SYSTEMATIC_PEAKS))
    flows = np.array([SYSTEMATIC_PEAKS[y] for y in years])
    return dict(
        peak_flows=flows,
        water_years=years,
        regional_skew=-0.302,
        regional_skew_mse=0.55**2,
        user_low_outlier_threshold=6000.0,
    )


@pytest.fixture(scope="module")
def big_sandy_6000():
    from flowfreq.bulletin17c import Bulletin17C

    native = Bulletin17C(**_big_sandy_6000()).run_analysis(engine="native")
    fortran = Bulletin17C(**_big_sandy_6000()).run_analysis(engine="fortran")
    return native, fortran


def test_big_sandy_6000_at_site_fit_matches(big_sandy_6000):
    native, fortran = big_sandy_6000
    assert native.ema_converged
    assert abs(native.skew_station - fortran.skew_station) < 1e-6


def test_fortran_mn2mvarb_does_not_converge_on_big_sandy_6000():
    from flowfreq._mse_ema import mc2mnvb, mn2mvarb
    from flowfreq._p3_moments import m2mn
    from flowfreq.fortran_engine import run_fortran_reference
    from flowfreq.peakfqr import _emafort

    kw = _big_sandy_6000()
    reference, _ = run_fortran_reference(
        kw["peak_flows"],
        kw["water_years"],
        user_low_outlier_threshold=6000.0,
        regional_skew=kw["regional_skew"],
        regional_skew_mse=kw["regional_skew_mse"],
        aeps=[0.01],
    )
    p = reference.parameters
    mc = np.array([0.0, p["std_log_at_site"] ** 2, p["skew_at_site"]])
    shift = p["mean_log_at_site"]
    s_mn = _emafort.var_mom([44.0], [np.log10(6000.0) - shift], [20.0 - shift], mc)
    mn = m2mn(mc)
    f_mc, f_s_mc = _emafort.mn2mvarb(mn, s_mn)
    n_mc, n_s_mc = mn2mvarb(mn, s_mn)
    assert np.max(np.abs(mc2mnvb(f_mc, f_s_mc) - s_mn)) > 1e-2
    assert np.max(np.abs(mc2mnvb(n_mc, n_s_mc) - s_mn)) < 1e-12


@pytest.mark.xfail(
    strict=True,
    reason="emafitpr's MN2MVARB returns a non-converged iterate here (as_G_mse 2.74 vs "
    "the true root's 0.064), so its weighted skew is -0.281 against native -0.166",
)
def test_big_sandy_6000_weighted_skew_matches(big_sandy_6000):
    native, fortran = big_sandy_6000
    assert abs(native.skew_used - fortran.skew_used) < 1e-4
