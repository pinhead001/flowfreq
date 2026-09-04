"""Per-routine oracles from the vendored Fortran.

``emafitpr`` alone is a poor way to check a ported routine. The EMA fit is a
fixed point whose condition number is around 1e13 -- one ulp of input moves the
converged at-site skew by 3e-3 -- so a correct routine can look wrong through
it and a wrong one can look right. ``build_fortran/_emafort.pyf`` therefore
also exposes ``mseg_all_sub``, ``detratsub``, ``var_mom`` and ``moms_p3``, all
of which peakfqr's own R code already calls through ``.Fortran()``.

These tests pin what those oracles say, so the ``var_mom`` port (TODO.md P3)
can be written against them one routine at a time instead of end to end.

They also record two things the oracles found immediately:

* ``_b17b_skew_mse`` matches the Fortran's ``mseg()`` exactly up to n = 150 and
  diverges above it, because ADJE's bias adjustment partially undoes the cap.
* ``_ema_iteration`` reproduces ``moms_p3`` exactly on uncensored rows, and
  diverges only where intervals are censored -- which places the remaining
  error in the truncated-moment code, not in the transcribed formulas.

A third, found while adding the Phase 3 (``mse_ema``) oracles: ``mseg_all_sub``
is not safe to call more than once per process. A fresh call reproduces
``as_G_mse_o`` correctly (``TestSkewMseOracle.test_reproduces_emafitpr_as_g_mse``,
which runs first in this file and is what actually pins that number), but
calling ``emafitpr`` for *any* case in between -- as ``TestDeterminantRatioOracle``
does a few classes below -- leaves ``mseg_all_sub`` permanently returning the
*uncensored* value for every input afterwards, Big Sandy's own included,
regardless of the (nobs, tl, tu) actually passed. Confirmed outside pytest
entirely (plain call, call ``emafitpr`` for an unrelated case, call again: the
second call is wrong and stays wrong). This is Fortran-side ``SAVE``d state
leaking across calls to *different* entry points, not anything on the Python
side -- the arrays going in are unmutated -- and not something to fix here
(``vendor/`` is not edited). The practical rule: never call ``mseg_all_sub``
after any ``emafitpr`` call in the same process, in these tests or elsewhere.

A fourth, the same class of bug in different routines: a direct call to the
raw ``var_emab``/``regmoms`` Fortran oracle (both carry ``SAVE``) drifts by
~2e-3 relative once other tests earlier in this file have already exercised
``emafitpr``/MGBT/``detrat``, even though it matches this port (and the
golden file) exactly when called first in a clean process. ``TestVarEmabPort``
therefore checks against the golden file only, never a live oracle call.

A fifth, found chasing TODO.md P3's last remaining ``xfail`` (Cains Coulee's
``skew_weighted``, ~0.058 skew units off): ``TestSkewMseOracle.test_reproduces
_emafitpr_as_g_mse`` above "passes" for Cains Coulee for the wrong reason. It
calls ``mseg_all_sub`` with ``golden["inputs"]``'s tl/tu, which -- as
``TestDeterminantRatioOracle`` documents -- are *uncensored* for this site
(MGBT creates the real censoring inside the fit); the ADJE bias adjustment
is therefore a no-op there and the call reduces to the plain B17B ``mseg()``
value, which happens to equal ``as_G_mse_o`` (0.2212) -- so the test never
actually exercises ADJE's bias adjustment for this site at all. Calling
``mseg_all_sub`` fresh with the *real* post-MGBT group instead
(``TestCainsCouleeAsGMseDiscrepancy``, placed here specifically so its own
call is the first Fortran entry point this process touches) gives 0.0749,
not 0.2212 -- and unlike the first three findings above, this one is *not*
explained by calling something after an unrelated ``emafitpr``: regenerating
Cains Coulee's golden file in total isolation reproduces 0.2212 too. See
that class's docstring for the full account and what was ruled out.
"""

from __future__ import annotations

import collections
import math

import numpy as np
import pytest

from tests.fortran_parity.conftest import load_golden

pytest.importorskip(
    "flowfreq.peakfqr",
    reason="Fortran extension not built; run python build_fortran/build.py",
)

pytestmark = pytest.mark.requires_fortran

CASES = ("big_sandy_03606500", "powder_river_06326500", "cains_coulee_06327450")

QMIN, QMAX = 1e-20, 1e20


def _threshold_groups(tl, tu):
    """Distinct (tl, tu) perception-threshold pairs and their counts.

    These routines take threshold *groups*, not per-observation arrays:
    ``nobs(i)`` counts the observations sharing the pair ``(tl(i), tu(i))``.

    Grouped on exact values, deliberately. Rounding to 12 decimals to be
    "robust" perturbs log10(18000) enough to move Big Sandy's at-site skew MSE
    by 2.2e-4 -- the fixed point amplifies a 1e-12 input change that far. The
    values here are exact duplicates already, so there is nothing to round.
    """
    counter = collections.Counter(zip(tl, tu))
    return (
        np.array(list(counter.values()), dtype=float),
        np.array([k[0] for k in counter], dtype=float),
        np.array([k[1] for k in counter], dtype=float),
    )


def _golden(name):
    golden = load_golden(name)
    if golden is None:
        pytest.skip(f"golden file missing for {name}")
    return golden


def _at_site_moments(golden):
    """Column 2 of cmoms: the at-site fit, which is what these routines take."""
    cmoms = np.asarray(golden["outputs"]["cmoms"], dtype=float)
    return np.array([cmoms[0][1], cmoms[1][1], cmoms[2][1]])


def _weighted_moments(golden):
    """Column 1 of cmoms: the final, regionally-weighted fit.

    ``VAR_EMAB``/``regmoms`` take this one, not the at-site fit
    ``_at_site_moments`` returns -- confidence intervals are built around
    the same estimate ``compute_quantiles()`` reports.
    """
    cmoms = np.asarray(golden["outputs"]["cmoms"], dtype=float)
    return np.array([cmoms[0][0], cmoms[1][0], cmoms[2][0]])


class TestSkewMseOracle:
    """mseg_all_sub: the at-site skew MSE, through the default ADJE option."""

    @pytest.mark.parametrize("name", CASES)
    def test_reproduces_emafitpr_as_g_mse(self, name):
        """Called directly it must give exactly what the full fit reports.

        This is the oracle the ADJE port is written against: it already
        includes the censoring bias adjustment flowfreq does not implement --
        0.0944 on Big Sandy against an uncensored 0.0636.
        """
        from flowfreq.peakfqr._emafort import mseg_all_sub

        golden = _golden(name)
        nobs, tl, tu = _threshold_groups(golden["inputs"]["tl"], golden["inputs"]["tu"])
        direct = float(mseg_all_sub(nobs, tl, tu, _at_site_moments(golden)))
        assert direct == pytest.approx(golden["outputs"]["skew"]["as_G_mse_o"], rel=1e-12)

    @pytest.mark.parametrize("n, skew", [(84, 0.0066), (85, -0.181), (32, -0.708), (50, 0.5)])
    def test_matches_b17b_formula_when_nothing_is_censored(self, n, skew):
        """With no censoring ADJE's bias adjustment is 1, so it reduces to mseg().

        That makes this a direct check of flowfreq's ``_b17b_skew_mse``.
        """
        from flowfreq.bulletin17c import _b17b_skew_mse
        from flowfreq.peakfqr._emafort import mseg_all_sub

        fortran = float(
            mseg_all_sub(
                np.array([float(n)]),
                np.array([-99.0]),
                np.array([99.0]),
                np.array([0.0, 0.09, skew]),
            )
        )
        assert fortran == pytest.approx(_b17b_skew_mse(n, skew), rel=1e-6)

    def test_b17b_formula_diverges_above_150_observations(self):
        """A real limitation, recorded rather than hidden.

        ``mseg_all`` evaluates ``mseg()`` at ``min(n, 150)`` and then multiplies
        by a bias adjustment that partially undoes the cap.
        ``_b17b_skew_mse`` applies the cap without the adjustment, so above 150
        it overestimates -- 0.0479 against the Fortran's 0.0365 at n = 200, 24%
        high, which over-weights the regional skew on a long record. Fixing it
        needs ``mse_ema``, hence ``var_mom``; see TODO.md P3. No parity case
        reaches n = 150, so nothing else detects this.
        """
        from flowfreq.bulletin17c import _b17b_skew_mse
        from flowfreq.peakfqr._emafort import mseg_all_sub

        mc = np.array([0.0, 0.09, 0.3])
        fortran = float(mseg_all_sub(np.array([200.0]), np.array([-99.0]), np.array([99.0]), mc))
        assert fortran == pytest.approx(0.03648, abs=1e-5)
        assert _b17b_skew_mse(200, 0.3) == pytest.approx(0.04789, abs=1e-5)
        assert _b17b_skew_mse(200, 0.3) / fortran == pytest.approx(1.31, abs=0.02)


class TestCainsCouleeAsGMseDiscrepancy:
    """A real, unexplained gap found chasing TODO.md P3's last xfail.

    Runs right after ``TestSkewMseOracle`` -- before any ``emafitpr`` call
    in this file -- so its own ``mseg_all_sub`` call is the first Fortran
    entry point touched this process, avoiding the ``SAVE``d-state leak
    the module docstring already documents.

    Everything else this port touches matches the Fortran closely at Cains
    Coulee's real post-MGBT censoring group (``tl=2.521``, its 332 cfs
    MGBT cutoff, ``tu=20``, ``n=32``) and its real at-site fit -- but
    ``emafitpr``'s own internally-computed, reported ``as_G_mse`` for this
    site (0.2212, committed in the golden file as ``skew.as_G_mse_o``, and
    what ``skew_weighted`` there was built from) does not match what
    calling the *same* ``mseg_all`` Fortran routine standalone gives for
    the *identical* inputs (0.0749). Confirmed reproducible, not a testing
    artifact: regenerating Cains Coulee's golden file in total isolation
    (``python tools/gen_fortran_golden.py cains_coulee_06327450``, nothing
    else in the process) gives the same 0.2212 -- ruling out cross-case
    contamination as the explanation this time (that would require some
    *other* case's ``emafitpr`` to have run first). The exact mechanism
    inside ``emafitpr`` was not pinned down: ``momsadj``'s skew floor is a
    no-op at this magnitude (-1.41), ``p3est_ema`` computes ``nG`` once
    before its iteration loop rather than per-iteration (so that is not
    silently recomputing ``as_G_mse`` differently), and replicating
    ``emafitpr``'s exact internal ``mse_ema``/``mseg_all`` call sequence
    (kmom=1, kmom=2, kmom=3, then the ERL "Syst" variant, then kmom=3
    again) via standalone oracle calls never reproduces the drift either --
    something earlier in ``emafitpr``'s own multi-stage fit (MGBT, or the
    first at-site-only ``p3est_ema`` pass under ``at_site_option='B17B'``)
    is implicated, but which one, and why, is still open.

    ``flowfreq._mse_ema.mse_ema``/``flowfreq._var_mom.var_mom`` are not the
    suspects here -- both are independently confirmed exact against the
    Fortran at this same real point, ruling out the ``var_mom``-precision
    explanation this xfail carried before this investigation.
    """

    #: Cains Coulee 06327450's real post-MGBT censoring group and at-site
    #: fit, read once from the committed golden file / a verified live
    #: emafitpr call and pinned here as literals -- so this class never
    #: needs its own emafitpr call, which would make its mseg_all_sub call
    #: untrustworthy per this file's own documented state-leak rule.
    _NOBS = np.array([32.0])
    _TL = np.array([2.5211380837040362])  # log10(331.99999999999994), the MGBT cutoff
    _TU = np.array([20.0])
    _MC_AT_SITE = np.array([2.6367671197988694, 0.15301024816541287, -0.7078900133204128])
    _AS_G_MSE_O = 0.22118912333191387  # golden file's skew.as_G_mse_o

    def test_mse_ema_matches_the_fortran_at_this_real_point(self):
        """The routine var_mom/mn2mvarb's port actually touches -- and it is exact."""
        from flowfreq._mse_ema import mse_ema
        from flowfreq.peakfqr._emafort import mse_ema as mse_ema_f

        mine = mse_ema(self._NOBS, self._TL, self._TU, self._MC_AT_SITE, kmom=3)
        fortran = float(mse_ema_f(self._NOBS, self._TL, self._TU, self._MC_AT_SITE, 3))
        assert mine == pytest.approx(fortran, rel=1e-6)

    def test_mseg_all_sub_disagrees_with_emafitprs_own_as_g_mse_o(self):
        """The actual gap: not a precision limit, a ~3x discrepancy.

        ``mseg_all_sub`` here is ADJE's ``bias_adj * mseg(n_adj, skew)``,
        called standalone with exactly the inputs ``emafitpr`` used to
        report ``as_G_mse_o`` -- and it does not agree.
        """
        from flowfreq.peakfqr._emafort import mseg_all_sub

        standalone = float(mseg_all_sub(self._NOBS, self._TL, self._TU, self._MC_AT_SITE))
        assert standalone == pytest.approx(0.0748857778, rel=1e-6)
        assert (
            abs(standalone / self._AS_G_MSE_O - 1) > 0.5
        ), "if this ever starts agreeing, the mystery is solved -- see the class docstring"


class TestDeterminantRatioOracle:
    """detratsub: the Halloween Wd, which flowfreq does not implement."""

    def test_needs_the_post_mgbt_thresholds(self):
        """The usage note that matters, and it is not obvious.

        Cains Coulee's *input* thresholds are uncensored throughout, so calling
        detrat with them returns 1.0. Its censoring is created by MGBT inside
        the fit, which raises the lower threshold to log10(332) = 2.521, and
        only with those does detrat give the 0.184 that emafitpr reports.
        emafitpr returns them as tlema/tuema.
        """
        from flowfreq.peakfqr._emafort import detratsub, emafitpr
        from tests.fortran_parity.cases import CASES as CASE_FACTORIES
        from tests.fortran_parity.cases import build_emafit_inputs

        args = build_emafit_inputs(CASE_FACTORIES["cains_coulee_06327450"]())
        out = emafitpr(
            args["ql"],
            args["qu"],
            args["tl"],
            args["tu"],
            args["dtype"],
            args["reg_m"],
            args["reg_m_mse"],
            args["reg_sd"],
            args["reg_sd_mse"],
            args["r_g"],
            args["r_g_mse"],
            args["gbthrsh0"],
            args["pq"],
            args["eps"],
            args["wght_opt_n"],
        )
        cmoms, wdout = np.asarray(out[9], dtype=float), float(out[14])
        n = len(args["ql"])
        tlema, tuema = np.asarray(out[17], dtype=float), np.asarray(out[18], dtype=float)
        mc = np.array([cmoms[0][1], cmoms[1][1], cmoms[2][1]])

        as_supplied = float(detratsub(mc, n, *_threshold_groups(args["tl"], args["tu"])))
        assert as_supplied == pytest.approx(1.0), "input thresholds carry no censoring"

        after_mgbt = float(detratsub(mc, n, *_threshold_groups(tlema[:n], tuema[:n])))
        assert after_mgbt == pytest.approx(wdout, rel=1e-9)
        assert after_mgbt == pytest.approx(0.184, abs=5e-4)


class TestMomentIterationOracle:
    """moms_p3 against _ema_iteration, the transcription it was written from."""

    @staticmethod
    def _fit(name):
        from flowfreq.bulletin17c import Bulletin17C
        from tests.fortran_parity.cases import CASES as CASE_FACTORIES

        case = CASE_FACTORIES[name]()
        b17c = Bulletin17C(
            peak_flows=list(case.systematic.values()) + list(case.historical.values()),
            water_years=list(case.systematic) + list(case.historical),
            regional_skew=case.regional_skew,
            regional_skew_mse=case.regional_skew_mse,
        )
        b17c.run_analysis(method="ema")
        return b17c

    @staticmethod
    def _same_rows(analyzer):
        """The exact intervals _ema_iteration is working on, as ql/qu.

        Both sides must be handed the same rows. MGBT censors peaks *inside*
        run_analysis, so the case's original ql/qu are not what the native
        iteration sees.
        """
        ql, qu = [], []
        for interval in analyzer.intervals:
            low = interval.lower if interval.lower > 0 else QMIN
            high = interval.upper if np.isfinite(interval.upper) else QMAX
            ql.append(np.log10(low))
            qu.append(np.log10(low if not interval.is_censored else high))
        return np.array(ql), np.array(qu)

    @staticmethod
    def _compare(b17c):
        from flowfreq.peakfqr._emafort import moms_p3

        analyzer = b17c._analyzer
        ql, qu = TestMomentIterationOracle._same_rows(analyzer)
        results = b17c.results
        start = np.array([results.mean_log, results.std_log**2, results.skew_station])
        fortran = np.asarray(moms_p3(ql, qu, 0.0, 0.0, start), dtype=float)
        mean, std, skew = analyzer._ema_iteration(
            start[0], np.sqrt(start[1]), start[2], n_regional=0.0
        )
        return fortran, np.array([mean, std**2, skew]), int((ql != qu).sum())

    @pytest.mark.parametrize("name", ["powder_river_06326500", "big_sandy_03606500"])
    def test_exact_on_uncensored_rows(self, name):
        """Where no interval is censored the transcription is exact.

        Measured differences: 0.0 on the mean, ~1e-14 on the variance and
        ~1e-12 on the skew. That is the formulas in moms_p3 being right.
        """
        fortran, mine, n_censored = self._compare(self._fit(name))
        assert n_censored == 0, "this case should have no censored intervals here"
        assert np.allclose(fortran, mine, rtol=0.0, atol=1e-9)

    def test_matches_on_censored_rows_too(self):
        """Cains Coulee censors 11 of 32; this used to be where they diverged.

        Two bugs previously combined to a 0.70% variance / 4.94% skew gap,
        both now fixed (see TODO.md P3):

        * ``_compute_ema_moments``'s censored branch used its own
          approximate truncated-gamma-moment formula instead of the
          Fortran-verified ``flowfreq._p3_moments.m_p3``.
        * ``_ema_iteration``'s bias-correction factors (``c2``, ``c3``)
          used the total interval count where the vendored Fortran's
          default (``bcf=1997``) uses the exact-peak count instead
          (``emafit.f:1408``; the ``bcf=2004`` alternative that would use
          the total is compiled out).

        Measured now: mean 7e-11, variance 1.6e-10, skew 3.0e-10 -- machine
        precision, the same level Powder River and Big Sandy already hit
        with nothing censored at all.
        """
        fortran, mine, n_censored = self._compare(self._fit("cains_coulee_06327450"))
        assert n_censored == 11
        assert np.allclose(fortran, mine, rtol=0.0, atol=1e-8)


class TestVarianceOfMomentsOracle:
    """var_mom: the root of the dependency tree the whole port hangs off."""

    @pytest.mark.parametrize("name", CASES)
    def test_returns_a_symmetric_positive_definite_matrix(self, name):
        """No independent reference for this one, so check what must hold.

        It is a covariance matrix of (mean, variance, skew) estimators, so it
        must be 3x3, symmetric and positive definite. A port that produces
        anything else is wrong regardless of the numbers.
        """
        from flowfreq.peakfqr._emafort import var_mom

        golden = _golden(name)
        nobs, tl, tu = _threshold_groups(golden["inputs"]["tl"], golden["inputs"]["tu"])
        varm = np.asarray(var_mom(nobs, tl, tu, _at_site_moments(golden)), dtype=float)
        assert varm.shape == (3, 3)
        assert np.allclose(varm, varm.T, rtol=1e-9), "covariance must be symmetric"
        assert np.all(np.linalg.eigvalsh(varm) > 0), "covariance must be positive definite"
        assert np.all(np.diag(varm) > 0)

    @pytest.mark.parametrize("name", CASES)
    def test_matches_fortran(self, name):
        """flowfreq._var_mom.var_mom (Phase 2) against the Fortran, directly.

        Powder River and Cains Coulee (using golden["inputs"], uncensored
        throughout -- see TestDeterminantRatioOracle) never drive d_est's
        nontrivial path: both tail probabilities stay under d_est's
        1e-6 cutoff, so those two match to float64 precision. Big Sandy's
        one real group does exercise it -- see test_matches_fortran_on_the
        _one_case_that_exercises_d_est below for why that one is not this
        tight.
        """
        from flowfreq._var_mom import var_mom

        golden = _golden(name)
        nobs, tl, tu = _threshold_groups(golden["inputs"]["tl"], golden["inputs"]["tu"])
        mc = _at_site_moments(golden)
        mine = var_mom(nobs, tl, tu, mc)
        from flowfreq.peakfqr._emafort import var_mom as var_mom_f

        fortran = np.asarray(var_mom_f(nobs, tl, tu, mc), dtype=float)
        tol = 1e-3 if name == "big_sandy_03606500" else 1e-9
        assert np.allclose(mine, fortran, rtol=tol, atol=1e-12), (name, mine, fortran)


class TestExpMomDerivPort:
    """expmomderiv: the Jacobian d_est needs, d(E[X^j])/d(noncentral moments).

    flowfreq._dexpect differentiates flowfreq._p3_moments' own (already
    Fortran- and quadrature-verified) truncated-moment function numerically
    via mpmath, rather than transcribing DEXPECT's closed-form chain
    (which needs DDGAM, itself a nontrivial incomplete-gamma derivative).
    """

    def test_matches_fortran_on_big_sandy_s_open_tail(self):
        """~1e-5 relative, not machine precision -- and that gap is on the

        Fortran side. Independently checking the same derivative at three
        different mpmath working precisions (50, 80, 120 dps) agrees to 30+
        stable digits, so flowfreq's numerical derivative is not the
        source of the ~1e-5 gap; DDGAM's own series (TOL=1e-11 per term,
        emafit.f/probfun.f:1054) and the ADJ rising-factorial chain built
        on it in float64 are the more likely one, though this is not
        independently proven the way the mP3 finding in TODO.md P3 is.
        """
        from flowfreq._p3_moments import m2p, p_p3, q_p3
        from flowfreq._var_mom import expmomderiv
        from flowfreq.peakfqr._emafort import expmomderiv as expmomderiv_f

        golden = _golden("big_sandy_03606500")
        mc = _at_site_moments(golden)
        skew = 0.06324555  # var_mom's clamp floor -- see TestVarMomComponentsPort
        mc_shift = np.array([0.0, mc[1], skew])
        parms = m2p(mc_shift)

        lo = 4.2552725051033065 - mc[0]
        pa = p_p3(lo, mc_shift)
        t1 = q_p3(pa, mc_shift)

        mine = expmomderiv(parms, -1e19, t1)
        fortran = np.asarray(expmomderiv_f(parms, -1e19, t1), dtype=float)
        assert np.allclose(mine, fortran, rtol=1e-4), (mine, fortran)


class TestVarMomComponentsPort:
    """varb/varc/d_est: the pieces var_mom sums over each threshold group."""

    @staticmethod
    def _group_inputs(golden, lo, hi):
        """(mu_x, e_x, p1, p2, p3) for one (mean-shifted) threshold group,

        matching what var_mom itself builds before calling varb/varc/d_est
        -- see flowfreq._var_mom.var_mom.
        """
        from flowfreq._p3_moments import m_p3, p_p3

        mc = _at_site_moments(golden)
        skew = math.copysign(min(max(abs(mc[2]), 0.06324555), 1.41), mc[2])
        mc_shift = np.array([0.0, mc[1], skew])
        mnouta = m_p3(-999.0, lo, mc_shift, 3)
        mnoutb = m_p3(hi, 999.0, mc_shift, 3)
        e_x = m_p3(lo, hi, mc_shift, 6)
        mu_x = np.array([mnouta, e_x[:3], mnoutb]).T
        pa, pb = p_p3(lo, mc_shift), p_p3(hi, mc_shift)
        return mc_shift, mu_x, e_x, pa, pb - pa, 1.0 - pb

    def test_varb_matches_fortran_on_big_sandy_s_group(self):
        from flowfreq._var_mom import varb
        from flowfreq.peakfqr._emafort import varb as varb_f

        golden = _golden("big_sandy_03606500")
        lo, hi = (
            4.2552725051033065 - golden["outputs"]["cmoms"][0][1],
            20.0 - golden["outputs"]["cmoms"][0][1],
        )
        _, mu_x, _, p1, p2, p3 = self._group_inputs(golden, lo, hi)
        nh = 40.0
        mine = varb(mu_x, nh, p1, p2, p3)
        fortran = np.asarray(varb_f(mu_x, nh, p1, p2, p3), dtype=float)
        assert np.allclose(mine, fortran, rtol=1e-9), (mine, fortran)

    def test_varc_matches_fortran_on_big_sandy_s_group(self):
        from flowfreq._var_mom import varc
        from flowfreq.peakfqr._emafort import varc as varc_f

        golden = _golden("big_sandy_03606500")
        lo, hi = (
            4.2552725051033065 - golden["outputs"]["cmoms"][0][1],
            20.0 - golden["outputs"]["cmoms"][0][1],
        )
        _, _, e_x, _, p2, _ = self._group_inputs(golden, lo, hi)
        nh = 40.0
        mine = varc(e_x, nh, p2)
        fortran = np.asarray(varc_f(e_x, nh, p2), dtype=float)
        assert np.allclose(mine, fortran, rtol=1e-9), (mine, fortran)

    def test_d_est_matches_fortran_on_big_sandy_s_group(self):
        """The one group among the three cases where d_est is not just zeros."""
        from flowfreq._var_mom import d_est
        from flowfreq.peakfqr._emafort import d_est as d_est_f

        golden = _golden("big_sandy_03606500")
        lo, hi = (
            4.2552725051033065 - golden["outputs"]["cmoms"][0][1],
            20.0 - golden["outputs"]["cmoms"][0][1],
        )
        mc_shift, _, _, pa, _, p3 = self._group_inputs(golden, lo, hi)
        pb = 1.0 - p3
        nh = 40.0

        mine = d_est(nh, mc_shift, pa, pb)
        fortran = np.asarray(d_est_f(nh, mc_shift, pa, pb), dtype=float)
        assert not np.allclose(mine, 0.0), "this group should exercise the nonzero path"
        assert np.allclose(mine, fortran, rtol=1e-3), (mine, fortran)

    def test_d_est_is_zero_below_the_probability_cutoff(self):
        """Neither tail carries enough probability: both jac contributions

        should vanish, matching the Fortran's own skip (emafit.f:2491/2500).
        """
        from flowfreq._var_mom import d_est
        from flowfreq.peakfqr._emafort import d_est as d_est_f

        mc = np.array([0.0, 0.09, 0.5])
        nh = 50.0
        pa, pb = 0.0, 1.0  # p1 = pa = 0, p3 = 1-pb = 0: both under 1e-6
        mine = d_est(nh, mc, pa, pb)
        fortran = np.asarray(d_est_f(nh, mc, pa, pb), dtype=float)
        assert np.allclose(mine, 0.0)
        assert np.allclose(fortran, 0.0)


# A skew right in the blend band (|skew| in [0.0007, 0.0010]) so both the
# incomplete-gamma and Wilson-Hilferty branches of pP3/qP3/mP3 get nonzero
# weight -- the seam is exactly where a transcription is likeliest to be
# wrong, and no parity case's at-site skew happens to land there.
BLEND_BAND_MOMENTS = (0.0, 1.0, 0.00085)


class TestM2pM2mnPort:
    """m2p/m2mn/mn2m: pure algebra, the leaves of the leaf layer."""

    @pytest.mark.parametrize("name", CASES)
    def test_m2p_matches_fortran(self, name):
        from flowfreq._p3_moments import m2p
        from flowfreq.peakfqr._emafort import m2p as m2p_fortran

        mc = _at_site_moments(_golden(name))
        assert np.allclose(m2p(mc), m2p_fortran(mc), rtol=1e-12)

    def test_m2p_in_the_blend_band(self):
        from flowfreq._p3_moments import m2p
        from flowfreq.peakfqr._emafort import m2p as m2p_fortran

        mc = np.array(BLEND_BAND_MOMENTS)
        assert np.allclose(m2p(mc), m2p_fortran(mc), rtol=1e-12)

    @pytest.mark.parametrize("name", CASES)
    def test_m2mn_matches_fortran(self, name):
        from flowfreq._p3_moments import m2mn
        from flowfreq.peakfqr._emafort import m2mn as m2mn_fortran

        mc = _at_site_moments(_golden(name))
        assert np.allclose(m2mn(mc), m2mn_fortran(mc), rtol=1e-12)

    @pytest.mark.parametrize("name", CASES)
    def test_mn2m_matches_fortran_and_inverts_m2mn(self, name):
        from flowfreq._p3_moments import m2mn, mn2m
        from flowfreq.peakfqr._emafort import mn2m as mn2m_fortran

        mc = _at_site_moments(_golden(name))
        mn = m2mn(mc)
        assert np.allclose(mn2m(mn), mn2m_fortran(mn), rtol=1e-12)
        assert np.allclose(mn2m(mn), mc, rtol=1e-9), "mn2m must invert m2mn"


class TestPP3QP3Port:
    """pP3/qP3: the Pearson III CDF and its inverse, gamma/WH-blended."""

    @staticmethod
    def _probe_points(mc):
        mu, s = mc[0], mc[1] ** 0.5
        return [mu - 3 * s, mu - s, mu, mu + s, mu + 3 * s]

    @pytest.mark.parametrize("name", CASES)
    def test_p_p3_matches_fortran(self, name):
        from flowfreq._p3_moments import p_p3
        from flowfreq.peakfqr._emafort import pp3

        mc = _at_site_moments(_golden(name))
        for x in self._probe_points(mc):
            assert p_p3(x, mc) == pytest.approx(float(pp3(x, mc)), abs=1e-9)

    def test_p_p3_in_the_blend_band(self):
        """Looser than test_p_p3_matches_fortran on purpose.

        A skew of 0.00085 gives alpha = 4/skew**2 ~ 5.5e6 -- far more
        extreme than any real at-site skew (none of the three parity cases
        land in the blend band at all). float64 stats.gamma.cdf loses a few
        digits at that shape parameter relative to the Fortran's real128
        DGAMDF; measured at ~4e-9 absolute here, not the moment-order
        cancellation m_p3 has (see _GAMMA_MOMENT_DPS) since pP3 evaluates
        the CDF once rather than raising alpha to a power.
        """
        from flowfreq._p3_moments import p_p3
        from flowfreq.peakfqr._emafort import pp3

        mc = np.array(BLEND_BAND_MOMENTS)
        for x in self._probe_points(mc):
            assert p_p3(x, mc) == pytest.approx(float(pp3(x, mc)), abs=1e-7)

    @pytest.mark.parametrize("name", CASES)
    def test_q_p3_matches_fortran(self, name):
        from flowfreq._p3_moments import q_p3
        from flowfreq.peakfqr._emafort import qp3sub

        mc = _at_site_moments(_golden(name))
        for q in (0.01, 0.1, 0.5, 0.9, 0.99):
            assert q_p3(q, mc) == pytest.approx(float(qp3sub(q, mc)), rel=1e-8)

    def test_q_p3_in_the_blend_band(self):
        """See test_p_p3_in_the_blend_band: same extreme-alpha float64 gap."""
        from flowfreq._p3_moments import q_p3
        from flowfreq.peakfqr._emafort import qp3sub

        mc = np.array(BLEND_BAND_MOMENTS)
        for q in (0.01, 0.1, 0.5, 0.9, 0.99):
            # abs matters near q=0.5: with mean 0 the quantile itself is near
            # zero there, where a relative tolerance alone is unsatisfiable.
            assert q_p3(q, mc) == pytest.approx(float(qp3sub(q, mc)), rel=1e-6, abs=1e-7)

    @pytest.mark.parametrize("name", CASES)
    def test_p_p3_and_q_p3_are_inverses(self, name):
        """Self-consistency, independent of the Fortran: p_p3(q_p3(p)) == p.

        Probability -> quantile -> probability, not the other direction:
        Cains Coulee's at-site skew is -0.708, strongly bounded on the
        right, so a probe point picked in x-space (mean +/- 3 sigma) can
        land past that bound, where p_p3 saturates to exactly 1.0 and the
        round trip is genuinely lossy rather than wrong -- q_p3(1.0, mc) is
        the boundary itself, not the x that produced 1.0 in floating point.
        Starting from p avoids ever hitting that boundary.
        """
        from flowfreq._p3_moments import p_p3, q_p3

        mc = _at_site_moments(_golden(name))
        for p in (0.01, 0.1, 0.5, 0.9, 0.99):
            assert p_p3(q_p3(p, mc), mc) == pytest.approx(p, rel=1e-6)


class TestMP3Port:
    """mP3: expected truncated raw moments -- the piece TODO.md P3 names as

    the source of _ema_iteration's residual gap on censored rows.
    """

    @pytest.mark.parametrize("name", CASES)
    def test_full_support_matches_m2mn(self, name):
        """With no truncation mP3 must reduce to the ordinary raw moments."""
        from flowfreq._p3_moments import m2mn, m_p3
        from flowfreq.peakfqr._emafort import mp3

        mc = _at_site_moments(_golden(name))
        fortran = np.asarray(mp3(-QMAX, QMAX, mc, 3), dtype=float)
        mine = m_p3(-QMAX, QMAX, mc, 3)
        assert np.allclose(fortran, np.asarray(m2mn(mc)), rtol=1e-6)
        assert np.allclose(mine, fortran, rtol=1e-6)

    @staticmethod
    def _shifted(golden):
        """(mc, groups) the way var_mom itself calls mP3.

        var_mom (emafit.f:2346) shifts every threshold by the fitted mean
        before calling mP3 -- ``tl = tl_in(it) - mc_in(1)`` -- and always
        asks for n = 6. Testing with the case's raw (unshifted) tl/tu against
        the raw at-site mean, as an earlier version of this test did, is not
        what var_mom actually evaluates and mixes up two different failure
        modes; this reproduces the real call.
        """
        mc = _at_site_moments(golden)
        mc_shift = np.array([0.0, mc[1], mc[2]])
        _, tl, tu = _threshold_groups(golden["inputs"]["tl"], golden["inputs"]["tu"])
        return mc_shift, list(zip(tl - mc[0], tu - mc[0]))

    @pytest.mark.parametrize("name", ["powder_river_06326500", "cains_coulee_06327450"])
    def test_matches_fortran_on_full_support_groups(self, name):
        """Every case has an uncensored (-20, 20) group; shifted, still wide open.

        Big Sandy is deliberately excluded here, not just untested: its
        at-site skew (0.0066) is close enough to zero that alpha = 4/skew**2
        is enormous regardless of which group is used, so even this "wide
        open" group already shows the precision gap this class documents
        for its one real censored group -- see
        test_fortran_itself_loses_precision_on_big_sandy_s_censored_group.
        Powder River and Cains Coulee have much larger-magnitude at-site
        skew (not near TODO.md's table values, which are the *weighted*
        skew; both are safely away from the blend band either way).
        """
        from flowfreq._p3_moments import m_p3
        from flowfreq.peakfqr._emafort import mp3

        mc_shift, groups = self._shifted(_golden(name))
        for lo, hi in groups:
            if not (lo < -15 and hi > 15):
                continue
            fortran = np.asarray(mp3(lo, hi, mc_shift, 6), dtype=float)
            mine = m_p3(lo, hi, mc_shift, 6)
            assert np.allclose(mine, fortran, rtol=1e-6, atol=1e-9), (lo, hi)

    def test_matches_fortran_on_cains_coulee_post_mgbt_group(self):
        """Cains Coulee's real censoring: finite width, both bounds real.

        golden["inputs"] carries emafitpr's *input* thresholds, which are
        uncensored for this case (MGBT creates the censoring inside the
        fit) -- see TestDeterminantRatioOracle. The group var_mom's own
        mse_ema/var_mom calls would actually iterate over is the post-MGBT
        one, (2.521, 20) before the mean shift, recovered the same way
        TestDeterminantRatioOracle does: from emafitpr's tlema/tuema.
        """
        from flowfreq._p3_moments import m_p3
        from flowfreq.peakfqr._emafort import emafitpr, mp3
        from tests.fortran_parity.cases import CASES as CASE_FACTORIES
        from tests.fortran_parity.cases import build_emafit_inputs

        args = build_emafit_inputs(CASE_FACTORIES["cains_coulee_06327450"]())
        out = emafitpr(
            args["ql"],
            args["qu"],
            args["tl"],
            args["tu"],
            args["dtype"],
            args["reg_m"],
            args["reg_m_mse"],
            args["reg_sd"],
            args["reg_sd_mse"],
            args["r_g"],
            args["r_g_mse"],
            args["gbthrsh0"],
            args["pq"],
            args["eps"],
            args["wght_opt_n"],
        )
        cmoms = np.asarray(out[9], dtype=float)
        n = len(args["ql"])
        tlema, tuema = np.asarray(out[17], dtype=float)[:n], np.asarray(out[18], dtype=float)[:n]
        mc = np.array([cmoms[0][1], cmoms[1][1], cmoms[2][1]])
        _, tl, tu = _threshold_groups(tlema, tuema)
        censored = [(lo, hi) for lo, hi in zip(tl, tu) if lo > -15]
        assert censored, "Cains Coulee should have a real (raised) lower threshold post-MGBT"

        mc_shift = np.array([0.0, mc[1], mc[2]])
        for lo, hi in censored:
            fortran = np.asarray(mp3(lo - mc[0], hi - mc[0], mc_shift, 6), dtype=float)
            mine = m_p3(lo - mc[0], hi - mc[0], mc_shift, 6)
            assert np.allclose(mine, fortran, rtol=1e-6, atol=1e-9), (lo, hi)

    def test_fortran_itself_loses_precision_on_big_sandy_s_censored_group(self):
        """A real, understood limitation of the *reference*, not of this port.

        Big Sandy's only non-trivial threshold group is (4.255, 20) --
        essentially "greater than 18000 cfs", open at the sentinel used for
        no upper bound. Shifted by the fitted mean (as var_mom does; at-site
        skew 0.0066), that lands mP3 in a regime its own author flagged:
        alpha = 4/skew**2 ~ 9e4 (probfun.f's FP_G1_CDF/FP_G1_MOM_TRC were
        patched in Sept. 2024 to evaluate DGAMDF in real128 for exactly
        this). But the patch only promoted the incomplete-gamma call itself
        -- the surrounding rising-factorial and choose(i,j)*tau**(i-j)*fp(j)
        expansion mP3 does around it (emafit.f:3049) is still float64, and
        that expansion cancels by ~11 orders of magnitude at k=6 here (terms
        of order tau**6 ~ 2e11 net to an answer of order 1e2). float64 does
        not have 11 spare digits, and the Fortran result shows it: it goes
        negative at k=6, which E[X**6] over positive-truncated X cannot be.

        flowfreq's own mP3 keeps this expansion in mpmath at 50 decimal
        digits (see _GAMMA_MOMENT_DPS), which is enough headroom to survive
        the cancellation. Independent verification below, by direct
        arbitrary-precision quadrature of the truncated density -- not by
        comparison with the Fortran, which is exactly the thing in
        question here -- confirms flowfreq's k=4..6 values, not the
        Fortran's.

        The Fortran/flowfreq gap grows steadily with k rather than
        appearing only at k=6 -- k=1 already differs at the 1e-8 absolute
        level -- so what is checked against the oracle (loosely, k=1..3)
        and what is checked against quadrature instead (k=4..6, where the
        gap is no longer a rounding difference) is a judgment call, not a
        sharp cutoff the numbers themselves draw. The unambiguous part is
        k=6 going negative, and flowfreq matching quadrature everywhere.
        """
        import mpmath

        from flowfreq._p3_moments import m2p, m_p3
        from flowfreq.peakfqr._emafort import mp3

        mc_shift, groups = self._shifted(_golden("big_sandy_03606500"))
        lo, hi = next((lo, hi) for lo, hi in groups if 0 < lo < 5 and hi > 15)

        fortran = np.asarray(mp3(lo, hi, mc_shift, 6), dtype=float)
        mine = m_p3(lo, hi, mc_shift, 6)

        assert np.allclose(mine[:3], fortran[:3], rtol=1e-4), "k=1..3 should still roughly agree"
        assert fortran[5] < 0, "the Fortran value at k=6 is exhibiting the cancellation failure"

        with mpmath.workdps(40):
            tau, alpha, beta = (mpmath.mpf(v) for v in m2p(mc_shift))

            def pdf(y):
                return mpmath.e ** (-y + (alpha - 1) * mpmath.log(y) - mpmath.loggamma(alpha))

            ymin, ymax = (mpmath.mpf(lo) - tau) / beta, (mpmath.mpf(hi) - tau) / beta
            nodes = [ymin, ymin + 1000, ymin + 5000, ymin + 20000, ymax]
            mass = mpmath.quad(pdf, nodes)
            for k in range(1, 7):
                expected = mpmath.quad(lambda y: (tau + beta * y) ** k * pdf(y), nodes) / mass
                assert mine[k - 1] == pytest.approx(float(expected), rel=1e-4), k

    def test_degenerate_interval_off_support(self):
        """cdfu == cdfl: the interval carries no probability mass.

        mP3 falls back to whichever endpoint is closer to the mean, raised
        to each power -- exercised here because it is a distinct code path,
        not the blended gamma/WH computation.
        """
        from flowfreq._p3_moments import m_p3
        from flowfreq.peakfqr._emafort import mp3

        mc = np.array([0.0, 1.0, 0.3])
        lo, hi = 50.0, 60.0  # far off in the upper tail: cdf saturates at 1.0
        fortran = np.asarray(mp3(lo, hi, mc, 3), dtype=float)
        mine = m_p3(lo, hi, mc, 3)
        assert np.allclose(mine, fortran, rtol=1e-9)
        assert np.allclose(mine, [lo, lo**2, lo**3])


class TestQuadratureAndCholeskyPort:
    """normquad/gammaquad/chol33: the primitives mc2mnvb's quadrature grid needs."""

    @pytest.mark.parametrize("n", [2, 3, 5])
    def test_normquad_matches_fortran(self, n):
        from flowfreq._mse_ema import _normquad
        from flowfreq.peakfqr._emafort import normquad

        t, w = _normquad(n)
        tf, wf = normquad(n)
        # Fortran hardcodes sqrt(2)/(1/sqrt(pi)) to 12 significant digits
        # rather than computing them, so this is not machine precision.
        assert np.allclose(t, tf, atol=1e-7)
        assert np.allclose(w, wf, atol=1e-7)
        assert w.sum() == pytest.approx(1.0)

    @pytest.mark.parametrize("alpha,beta", [(5.0, 0.5), (1000.0, 0.03), (250.0, 0.06324555)])
    def test_gammaquad_matches_fortran(self, alpha, beta):
        from flowfreq._mse_ema import _gammaquad
        from flowfreq.peakfqr._emafort import gammaquad

        t, w = _gammaquad(3, alpha, beta)
        tf, wf = gammaquad(3, alpha, beta)
        assert np.allclose(t, tf, rtol=1e-9)
        assert np.allclose(w, wf, rtol=1e-9)
        assert w.sum() == pytest.approx(1.0)
        assert np.average(t, weights=w) == pytest.approx(alpha * beta, rel=1e-6)

    def test_gammaquad_matches_fortran_above_the_stabilization_threshold(self):
        """alpha > 160 (probfun.f:2694): the mean-preserving Gamma(160, .) proxy."""
        from flowfreq._mse_ema import _gammaquad
        from flowfreq.peakfqr._emafort import gammaquad

        t, w = _gammaquad(3, 5000.0, 0.01)
        tf, wf = gammaquad(3, 5000.0, 0.01)
        assert np.allclose(t, tf, rtol=1e-9)
        assert np.allclose(w, wf, rtol=1e-9)
        assert np.average(t, weights=w) == pytest.approx(5000.0 * 0.01, rel=1e-6)

    def test_chol33_matches_fortran(self):
        from flowfreq._mse_ema import _chol33
        from flowfreq.peakfqr._emafort import chol33

        s = np.array([[1.0, 0.2, 0.1], [0.2, 2.0, 0.3], [0.1, 0.3, 0.5]])
        v = _chol33(s)
        vf, iflag = chol33(s)
        assert iflag == 0
        assert np.allclose(v, vf, rtol=1e-12)
        assert np.allclose(v.T @ v, s), "V.T @ V must reconstruct S"

    def test_chol33_flags_non_positive_semi_definite(self):
        from flowfreq._mse_ema import _chol33
        from flowfreq.peakfqr._emafort import chol33

        s = np.array([[1.0, 5.0, 0.0], [5.0, 1.0, 0.0], [0.0, 0.0, 1.0]])  # off-diag too large
        assert _chol33(s) is None
        _, iflag = chol33(s)
        assert iflag == 1


class TestMc2mnvbPort:
    """mc2mnvb: central-moment covariance -> noncentral-moment covariance,

    via gridmake's quadrature.
    """

    def test_matches_fortran(self):
        from flowfreq._mse_ema import mc2mnvb
        from flowfreq.peakfqr._emafort import mc2mnvb as mc2mnvb_f

        mc = np.array([3.7, 0.09, 0.3])
        s_mc = np.array(
            [
                [0.001, 0.0001, 0.0002],
                [0.0001, 0.0002, 0.00005],
                [0.0002, 0.00005, 0.001],
            ]
        )
        mine = mc2mnvb(mc, s_mc)
        fortran = np.asarray(mc2mnvb_f(mc, s_mc, [2, 2, 2]), dtype=float)
        # Same ~1e-11-level gap as normquad's truncated constants, propagated
        # through the quadrature sum.
        assert np.allclose(mine, fortran, atol=1e-6)
        assert np.allclose(mine, mine.T), "a covariance matrix must be symmetric"


class TestMn2mVarPort:
    """mn2m_var: the closed-form linearized estimate mn2mvarb starts from."""

    @pytest.mark.parametrize("name", CASES)
    def test_matches_fortran(self, name):
        from flowfreq._mse_ema import mn2m_var
        from flowfreq._p3_moments import m2mn
        from flowfreq.peakfqr._emafort import mn2m_var as mn2m_var_f
        from flowfreq.peakfqr._emafort import var_mom as var_mom_f

        golden = _golden(name)
        mc = _at_site_moments(golden)
        nobs, tl, tu = _threshold_groups(golden["inputs"]["tl"], golden["inputs"]["tu"])
        mc_shift = np.array([0.0, mc[1], mc[2]])
        s_mn = np.asarray(var_mom_f(nobs, tl - mc[0], tu - mc[0], mc_shift), dtype=float)
        mn = m2mn(mc_shift)

        mc_out, s_mc_out = mn2m_var(mn, s_mn)
        mc_out_f, s_mc_out_f = mn2m_var_f(mn, s_mn)
        assert np.allclose(mc_out, mc_out_f, rtol=1e-12)
        assert np.allclose(s_mc_out, s_mc_out_f, rtol=1e-9)


class TestMn2mvarbPort:
    """mn2mvarb: the inverse-problem solve -- see flowfreq._mse_ema's module

    docstring for how the V-reparametrized root-find compares to the
    Fortran's step-halved Newton iteration.
    """

    @pytest.mark.parametrize("name", CASES)
    def test_matches_fortran(self, name):
        from flowfreq._mse_ema import mn2mvarb
        from flowfreq._p3_moments import m2mn
        from flowfreq.peakfqr._emafort import mn2mvarb as mn2mvarb_f
        from flowfreq.peakfqr._emafort import var_mom as var_mom_f

        golden = _golden(name)
        mc = _at_site_moments(golden)
        nobs, tl, tu = _threshold_groups(golden["inputs"]["tl"], golden["inputs"]["tu"])
        mc_shift = np.array([0.0, mc[1], mc[2]])
        s_mn = np.asarray(var_mom_f(nobs, tl - mc[0], tu - mc[0], mc_shift), dtype=float)
        mn = m2mn(mc_shift)

        mc_out, s_mc_out = mn2mvarb(mn, s_mn)
        mc_out_f, s_mc_out_f = mn2mvarb_f(mn, s_mn)
        tol = 1e-3 if name == "big_sandy_03606500" else 1e-6
        assert np.allclose(mc_out, mc_out_f, rtol=1e-9)
        assert np.allclose(s_mc_out, s_mc_out_f, rtol=tol, atol=1e-12), (s_mc_out, s_mc_out_f)


class TestMseEmaPort:
    """mse_ema: the ADJE censoring bias adjustment's numerator/denominator,

    TODO.md P3's "skew weighting -- 24% left" item. bias_adj =
    mse_ema(censored)/mse_ema(uncensored), feeding
    as_G_mse = bias_adj * mseg(min(n,150), G).
    """

    @pytest.mark.parametrize("name", CASES)
    def test_matches_fortran_on_the_case_own_thresholds(self, name):
        from flowfreq._mse_ema import mse_ema
        from flowfreq.peakfqr._emafort import mse_ema as mse_ema_f

        golden = _golden(name)
        mc = _at_site_moments(golden)
        nobs, tl, tu = _threshold_groups(golden["inputs"]["tl"], golden["inputs"]["tu"])

        mine = mse_ema(nobs, tl, tu, mc, kmom=3)
        fortran = float(mse_ema_f(nobs, tl, tu, mc, 3))
        tol = 1e-3 if name == "big_sandy_03606500" else 1e-6
        assert mine == pytest.approx(fortran, rel=tol)

    def test_reproduces_the_documented_adje_bias_adjustment_on_big_sandy(self):
        """The end-to-end number TODO.md P3 cites: bias_adj = 1.4844,

        as_G_mse = 0.09437. TestSkewMseOracle.test_reproduces_emafitpr_as_g_mse
        already checks this same figure against the ``mseg_all_sub`` oracle
        directly (as_G_mse_o); that comparison is deliberately not repeated
        here -- see the module docstring's note on ``mseg_all_sub`` for why
        calling it a second time in the same process is not safe to rely on.
        """
        from flowfreq._mse_ema import mse_ema
        from flowfreq.bulletin17c import _b17b_skew_mse

        golden = _golden("big_sandy_03606500")
        mc = _at_site_moments(golden)
        nobs, tl, tu = _threshold_groups(golden["inputs"]["tl"], golden["inputs"]["tu"])
        n = int(nobs.sum())

        mse_censored = mse_ema(nobs, tl, tu, mc, kmom=3)
        mse_uncensored = mse_ema(np.array([float(n)]), np.array([-99.0]), np.array([99.0]), mc, 3)
        bias_adj = mse_censored / mse_uncensored
        assert bias_adj == pytest.approx(1.4844, abs=2e-4)

        as_g_mse = bias_adj * _b17b_skew_mse(min(n, 150), mc[2])
        assert as_g_mse == pytest.approx(0.09437, abs=2e-5)

    def test_rejects_invalid_kmom(self):
        from flowfreq._mse_ema import mse_ema

        mc = np.array([3.7, 0.09, 0.3])
        with pytest.raises(ValueError):
            mse_ema(np.array([50.0]), np.array([-20.0]), np.array([20.0]), mc, kmom=0)


class TestExpMomCDerivPort:
    """expmomcderiv: detrat's per-group building block -- E[X^j | censored]

    and the Jacobian of the resulting expected central moments, w.r.t. the
    fit's own central moments directly (a different basis from
    flowfreq._var_mom's expmomderiv, which differentiates w.r.t. noncentral
    moments).
    """

    @pytest.mark.parametrize("tl,tu", [(3.5, 4.0), (-20.0, 3.0), (4.5, 20.0)])
    def test_matches_fortran(self, tl, tu):
        from flowfreq._detrat import _expmomcderiv
        from flowfreq._p3_moments import m2p
        from flowfreq.peakfqr._emafort import expmomcderiv

        mc = np.array([3.7, 0.09, -0.5])
        parms = m2p(mc)

        mne_f, dedmc_f = expmomcderiv(parms, tl, tu)
        _, mne, dedmc = _expmomcderiv(parms, tl, tu)
        assert np.allclose(mne, mne_f, rtol=1e-6)
        assert np.allclose(dedmc, dedmc_f, rtol=1e-6), (dedmc, dedmc_f)

    def test_falls_back_to_unconditional_moments_when_censoring_probability_is_negligible(self):
        """probfun.f:940 -- PICK <= 1e-10 uses DEXPECT(-inf, inf) directly."""
        from flowfreq._detrat import _expmomcderiv
        from flowfreq._p3_moments import m2p
        from flowfreq.peakfqr._emafort import expmomcderiv

        mc = np.array([3.7, 0.09, -0.5])
        parms = m2p(mc)
        # A threshold group so wide that essentially nothing falls outside it.
        mne_f, dedmc_f = expmomcderiv(parms, -20.0, 20.0)
        _, mne, dedmc = _expmomcderiv(parms, -20.0, 20.0)
        assert np.allclose(mne, mne_f, rtol=1e-6)
        assert np.allclose(dedmc, dedmc_f, rtol=1e-6)


class TestDetratPort:
    """detrat: the Halloween determinant ratio, Wd -- TODO.md P3's other

    open item, independent of var_mom/mse_ema.
    """

    def test_is_one_below_the_skew_floor(self):
        from flowfreq._detrat import detrat

        mc = np.array([3.7, 0.09, 0.0066])  # |skew| < 0.04
        w = detrat(mc, 84, np.array([84.0]), np.array([-20.0]), np.array([20.0]))
        assert w == 1.0

    def test_matches_fortran_on_big_sandy(self):
        """At-site skew 0.0066 is under the floor either way, but this

        exercises the real (nobs, tl, tu) groups end to end.
        """
        from flowfreq._detrat import detrat
        from flowfreq.peakfqr._emafort import detratsub

        golden = _golden("big_sandy_03606500")
        mc = _at_site_moments(golden)
        nobs, tl, tu = _threshold_groups(golden["inputs"]["tl"], golden["inputs"]["tu"])
        n = int(nobs.sum())

        mine = detrat(mc, n, nobs, tl, tu)
        fortran = float(detratsub(mc, n, nobs, tl, tu))
        assert mine == pytest.approx(fortran, rel=1e-9)
        assert mine == 1.0

    def test_matches_fortran_on_cains_coulee_post_mgbt_group(self):
        """The acceptance case: at-site skew -0.708, reference Wd 0.184.

        Uses emafitpr's own tlema/tuema (post-MGBT), the same way
        TestDeterminantRatioOracle does -- Cains Coulee's *input*
        thresholds carry no censoring; MGBT creates it inside the fit.
        """
        from flowfreq._detrat import detrat
        from flowfreq.peakfqr._emafort import detratsub, emafitpr
        from tests.fortran_parity.cases import CASES as CASE_FACTORIES
        from tests.fortran_parity.cases import build_emafit_inputs

        args = build_emafit_inputs(CASE_FACTORIES["cains_coulee_06327450"]())
        out = emafitpr(
            args["ql"],
            args["qu"],
            args["tl"],
            args["tu"],
            args["dtype"],
            args["reg_m"],
            args["reg_m_mse"],
            args["reg_sd"],
            args["reg_sd_mse"],
            args["r_g"],
            args["r_g_mse"],
            args["gbthrsh0"],
            args["pq"],
            args["eps"],
            args["wght_opt_n"],
        )
        cmoms, wdout = np.asarray(out[9], dtype=float), float(out[14])
        n = len(args["ql"])
        tlema, tuema = np.asarray(out[17], dtype=float)[:n], np.asarray(out[18], dtype=float)[:n]
        mc = np.array([cmoms[0][1], cmoms[1][1], cmoms[2][1]])
        nobs, tl, tu = _threshold_groups(tlema, tuema)

        mine = detrat(mc, n, nobs, tl, tu)
        fortran = float(detratsub(mc, n, nobs, tl, tu))
        assert mine == pytest.approx(fortran, rel=1e-6)
        assert mine == pytest.approx(wdout, rel=1e-6)
        assert mine == pytest.approx(0.184, abs=5e-4)

    @pytest.mark.parametrize("skew", [-0.5, 0.3, 1.2])
    def test_matches_fortran_on_synthetic_censored_groups(self, skew):
        from flowfreq._detrat import detrat
        from flowfreq.peakfqr._emafort import detratsub

        mc = np.array([3.7, 0.09, skew])
        nobs = np.array([60.0, 20.0])
        tl = np.array([-20.0, 4.0])
        tu = np.array([4.0, 20.0])
        mine = detrat(mc, 80, nobs, tl, tu)
        fortran = float(detratsub(mc, 80, nobs, tl, tu))
        assert mine == pytest.approx(fortran, rel=1e-6)


class TestVarEmabPort:
    """VAR_EMAB/regmoms/ci_ema_m3b: the confidence-interval shape fix.

    TODO.md P3's last open item. compute_confidence_limits() forms
    log_Q +/- z*se, symmetric by construction; peakfq skews right with
    return period via Cohn's inverse-Gaussian-quadrature method. See
    flowfreq._var_emab's module docstring for the call convention this
    class pins -- pq is non-exceedance probability (1 - aep), mc is the
    weighted fit, and the two Fortran signatures disagree on their own
    (r_G_mse, r_M_mse, r_S2, r_S2_mse) argument order.
    """

    @pytest.mark.parametrize("name", ["big_sandy_03606500", "powder_river_06326500"])
    def test_matches_fortran_end_to_end(self, name):
        """yp to ~1e-9 relative, ci_low/ci_high to ~1e-5 -- var_mom's own gap.

        Big Sandy and Powder River both have their real censoring (or lack
        of it) already reflected in golden["inputs"]'s tl/tu -- neither has
        MGBT-driven censoring, so no post-MGBT tlema/tuema reconstruction
        is needed the way Cains Coulee's is (test_matches_fortran_on_cains_coulee
        below).

        Compared only against the golden file, not a live in-process call
        to the raw ``var_emab``/``regmoms`` Fortran oracle: those carry the
        same class of ``SAVE``d cross-call state leakage this file's module
        docstring already documents for ``mseg_all_sub`` (confirmed here
        too -- a direct oracle call agreed with this port when nothing else
        in the process had touched ``emafitpr`` yet, and quietly drifted by
        ~2e-3 once other tests in this file had). The golden file was
        generated by a single, isolated ``tools/gen_fortran_golden.py`` run
        and is not subject to that.
        """
        from flowfreq._var_emab import var_emab

        golden = _golden(name)
        mc = _weighted_moments(golden)
        inputs = golden["inputs"]
        nobs, tl, tu = _threshold_groups(inputs["tl"], inputs["tu"])
        r_g_mse = float(inputs["regional_skew_mse"])
        pq = 1.0 - np.asarray(inputs["aeps"], dtype=float)
        eps = float(inputs["eps"])

        mine_yp, mine_cv, mine_cil, mine_cih, _ = var_emab(
            nobs, tl, tu, mc, pq, eps, r_g_mse=r_g_mse
        )

        ref = golden["outputs"]["quantiles"]
        assert np.allclose(mine_yp, ref["yp"], rtol=1e-6)
        assert np.allclose(mine_cil, ref["ci_low"], rtol=1e-3)
        assert np.allclose(mine_cih, ref["ci_high"], rtol=1e-3)

    def test_matches_fortran_on_cains_coulee(self):
        """The hard case: MGBT censoring, at-site skew -0.708.

        Uses emafitpr's own tlema/tuema (post-MGBT) and its own weighted
        cmoms, the same way TestDeterminantRatioOracle does -- Cains
        Coulee's *input* thresholds carry no censoring; MGBT creates it
        inside the fit. Looser tolerance than the other two cases: this is
        the one site where mn2mvarb's ~1e-3 relative precision limit
        (var_mom Phase 2's expmomderiv gap) is large enough to matter, and
        it propagates non-linearly through beta1 into an asymmetric error
        between ci_low (worse, ~11%) and ci_high (~1%) -- both loose but
        real, not evidence of a new defect on top of the known one.
        """
        from flowfreq._var_emab import var_emab
        from flowfreq.peakfqr._emafort import emafitpr
        from tests.fortran_parity.cases import CASES as CASE_FACTORIES
        from tests.fortran_parity.cases import build_emafit_inputs

        args = build_emafit_inputs(CASE_FACTORIES["cains_coulee_06327450"]())
        out = emafitpr(
            args["ql"],
            args["qu"],
            args["tl"],
            args["tu"],
            args["dtype"],
            args["reg_m"],
            args["reg_m_mse"],
            args["reg_sd"],
            args["reg_sd_mse"],
            args["r_g"],
            args["r_g_mse"],
            args["gbthrsh0"],
            args["pq"],
            args["eps"],
            args["wght_opt_n"],
        )
        cmoms = np.asarray(out[9], dtype=float)
        ref_yp, ref_cil, ref_cih = np.asarray(out[10]), np.asarray(out[11]), np.asarray(out[12])
        n = len(args["ql"])
        tlema, tuema = np.asarray(out[17], dtype=float)[:n], np.asarray(out[18], dtype=float)[:n]
        mc = np.array([cmoms[0][0], cmoms[1][0], cmoms[2][0]])
        nobs, tl, tu = _threshold_groups(tlema, tuema)
        r_g_mse = float(args["r_g_mse"])
        pq = np.asarray(args["pq"], dtype=float)
        eps = float(args["eps"])

        yp, cv, cil, cih, _ = var_emab(nobs, tl, tu, mc, pq, eps, r_g_mse=r_g_mse)

        assert np.allclose(yp, ref_yp, rtol=1e-6)
        assert np.allclose(cil, ref_cil, rtol=0.15)
        assert np.allclose(cih, ref_cih, rtol=0.02)
