"""Native EMA against peakfq 8.1.0 on two Wyoming/Montana sites.

Big Sandy was the only parity site for a long time, and it cannot answer two
questions on its own.

Does the fit agree when nothing is censored?
    Powder River is 85 contiguous water years with no historic peaks, no zero
    flows and no PILFs -- a plain systematic record fitted with a weighted
    regional skew. It answers yes, to machine precision, which is what says
    the in-loop regional-skew weighting is right rather than merely closer.

What does censoring cost?
    Cains Coulee is 32 contiguous years whose MGBT finds 11 PILFs, so the
    censoring is produced by the fit rather than supplied as input. Its
    reference ``Wd`` is 0.184 -- the first case in this repository where the
    Halloween determinant ratio is far from 1, since ``emafit.f:763`` only
    reaches ``detrat`` when the at-site skew clears 0.04 and Big Sandy's is
    0.0066. ``detrat`` and the ADJE bias adjustment are ported and wired
    (``hydrolib._detrat``, ``hydrolib._mse_ema``), and the at-site EMA
    moment iteration itself (``_compute_ema_moments``/``_ema_iteration``)
    now matches ``moms_p3`` to machine precision on censored rows too (it
    used to diverge there, from its own approximate truncated-moment
    formula and a bias-correction factor sized off the wrong count -- see
    ``tests/fortran_parity/test_fortran_oracles.py::TestMomentIterationOracle``).
    Native ``Wd`` here is 0.186, close to peakfq's 0.184; native at-site
    skew is -0.708 against peakfq's -0.708. What's left is downstream of
    all of that: ``skew_weighted`` still carries a real gap (0.058 skew
    units), and it is *not* a `var_mom`/`mn2mvarb` precision limit --
    `mse_ema` called standalone with this site's real post-MGBT group
    matches the Fortran oracle to 3e-8 relative. It traces to `emafitpr`'s
    own internally-computed `as_G_mse` for this site disagreeing with what
    the same `mseg_all` routine gives called standalone on identical
    inputs; see the xfail below for the full account.

Read together they localise the open P3 defect precisely: with no censoring
the native fit is exact, and every ported piece -- the at-site moment
iteration, the ADJE bias adjustment, `detrat`, and (independently, per
routine) `mse_ema`/`var_mom`/`mn2mvarb` -- matches the Fortran to a
demonstrated tolerance, including at this site's own real, sensitive input.
What remains is a discrepancy inside `emafitpr` itself that a clean
composition of correctly-verified routines has no principled way to
reproduce -- see the xfail below.

The peakfq 7.4 columns in the fixture CSVs are a sanity cross-check only.
Parity is against the committed goldens, generated from the vendored 8.1.0
Fortran by ``tools/gen_fortran_golden.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE
from tests.fortran_parity.conftest import load_golden

# The marker alone does not skip anything -- it only labels. Without the
# skipif the module-scope fixtures raise FileNotFoundError and every test in
# here reports as an ERROR rather than a skip when the reference tree is
# absent. Same pairing as tests/test_r_fixtures.py.
pytestmark = [
    pytest.mark.requires_peakfqr_testdata,
    pytest.mark.skipif(not TESTDATA_AVAILABLE, reason=SKIP_REASON),
]


def _native(site_no: str):
    from hydrolib.bulletin17c import Bulletin17C
    from tests.fixtures.wymt_peaks import load_site

    site = load_site(site_no)
    years = sorted(site.peaks)
    b17c = Bulletin17C(
        peak_flows=[site.peaks[y] for y in years],
        water_years=years,
        regional_skew=site.regional_skew,
        regional_skew_mse=site.regional_skew_mse,
    )
    b17c.run_analysis(method="ema")
    return b17c.results


def _reference(name: str):
    golden = load_golden(name)
    if golden is None:
        pytest.skip(f"golden file missing for {name} (run tools/gen_fortran_golden.py)")
    cmoms = np.asarray(golden["outputs"]["cmoms"], dtype=float)
    return {
        "mean_log": cmoms[0][0],
        "std_log": float(np.sqrt(cmoms[1][0])),
        "skew_weighted": cmoms[2][0],
        "skew_at_site": cmoms[2][1],
        "aeps": np.asarray(golden["inputs"]["aeps"], dtype=float),
        "quantiles": 10.0 ** np.asarray(golden["outputs"]["quantiles"]["yp"], dtype=float),
        "mgbt": golden["outputs"]["mgbt"],
        "wd": golden["outputs"]["skew"]["Wdout"],
    }


def _quantile_errors(results, ref) -> dict:
    """Percent error per AEP, for the AEPs the native fit reports."""
    q = results.quantiles
    errors = {}
    for aep, ref_q in zip(ref["aeps"], ref["quantiles"]):
        row = q[np.isclose(q["aep"], aep)]
        if not row.empty:
            errors[float(aep)] = abs(float(row["flow_cfs"].iloc[0]) / ref_q - 1) * 100
    return errors


@pytest.fixture(scope="module")
def powder_river():
    """(native results, reference) for USGS 06326500."""
    return _native("06326500.00"), _reference("powder_river_06326500")


@pytest.fixture(scope="module")
def cains_coulee():
    """(native results, reference) for USGS 06327450."""
    return _native("06327450.00"), _reference("cains_coulee_06327450")


class TestPowderRiverUncensored:
    """USGS 06326500, 1938-2022. No censoring anywhere, and it agrees exactly.

    This is the test that makes the weighting fix a claim rather than a
    direction: on a systematic record with a weighted regional skew, the native
    EMA reproduces peakfq 8.1.0 to machine precision. Measured absolute
    differences are 0.0 on the mean, 3.7e-14 on the standard deviation,
    4.5e-12 on the at-site skew and 7.5e-11 on the weighted skew. The bound
    below is 1e-6 -- four orders of magnitude looser than that, and four
    tighter than anything the open defect could hide in.
    """

    @pytest.mark.parametrize("field", ["mean_log", "std_log", "skew_weighted", "skew_at_site"])
    def test_moment_matches_to_machine_precision(self, powder_river, field):
        results, ref = powder_river
        native = {
            "mean_log": results.mean_log,
            "std_log": results.std_log,
            "skew_weighted": results.skew_weighted,
            "skew_at_site": results.skew_station,
        }[field]
        assert (
            abs(native - ref[field]) < 1e-6
        ), f"{field}: native {native!r} vs peakfq {ref[field]!r}"

    def test_no_pilfs_either_side(self, powder_river):
        """MGBT must agree that there is nothing to censor."""
        results, ref = powder_river
        assert results.n_low_outliers == ref["mgbt"]["gbnlow"] == 0

    def test_determinant_ratio_is_one_without_censoring(self, powder_river):
        """Wd = 1 here, so the missing detrat cannot be what is being measured."""
        _, ref = powder_river
        assert ref["wd"] == pytest.approx(1.0)

    def test_quantiles_agree_within_half_a_percent(self, powder_river):
        results, ref = powder_river
        errors = _quantile_errors(results, ref)
        worst = max(errors, key=errors.get)
        assert errors[worst] < 0.5, f"worst at AEP {worst}: {errors[worst]:.3f}%"


class TestCainsCouleeCensored:
    """USGS 06327450, 1991-2022. MGBT censors 11 peaks, and the fit diverges.

    Everything discrete matches: the same 11 PILFs at the same 332 cfs cut.
    The at-site skew now matches too (0.0002 gap), and so does ``Wd`` (0.002
    gap against peakfq's 0.184). What still doesn't match is
    ``skew_weighted`` -- a real residual, and (measured, not assumed) not
    a `var_mom`/`mse_ema` precision limit; see `test_weighted_skew_matches`'s
    xfail reason.
    """

    def test_mgbt_finds_the_same_pilfs(self, cains_coulee):
        """MGBT is the part verified line-by-line against the Fortran; it holds here."""
        results, ref = cains_coulee
        assert results.n_low_outliers == ref["mgbt"]["gbnlow"] == 11
        assert results.low_outlier_threshold == pytest.approx(332.0, rel=1e-3)

    def test_reference_exercises_the_determinant_ratio(self, cains_coulee):
        """The point of this case: Wd is nowhere near 1, unlike every other case."""
        _, ref = cains_coulee
        assert ref["wd"] == pytest.approx(0.184, abs=5e-4)

    def test_native_determinant_ratio_is_close(self, cains_coulee):
        """hydrolib's own Wd, computed from its now near-exact at-site fit.

        0.186 against peakfq's 0.184 -- close, and nowhere near the
        ``Wd = 1`` that ``_perception_threshold_groups`` reported before it
        accounted for MGBT's low-outlier threshold, or the 0.174 it reported
        before the at-site EMA moment iteration itself was fixed (feeding
        ``skew_station`` -0.830 into ``detrat`` rather than the now-correct
        -0.708).
        """
        from hydrolib.bulletin17c import ExpectedMomentsAlgorithm

        results, ref = cains_coulee
        nobs, tl, tu = (
            np.array([32.0]),
            np.array([np.log10(results.low_outlier_threshold)]),
            np.array([np.log10(1e20)]),
        )
        wd = ExpectedMomentsAlgorithm._detrat_wd(
            tuple(nobs),
            tuple(tl),
            tuple(tu),
            results.mean_log,
            results.std_log**2,
            results.skew_station,
            32,
        )
        assert abs(wd - ref["wd"]) < 0.02

    def test_mean_still_agrees_closely(self, cains_coulee):
        """Measured 3.1e-3 in log10, ~0.7% in flow -- was 5.9e-3/~1.4%.

        Moved when the at-site EMA moment iteration was fixed (TODO.md P3),
        even though the mean was never the part with a known defect: both
        bugs it fixed (the truncated-moment formula, the bias-correction
        count) touch every censored-row moment, mean included.
        """
        results, ref = cains_coulee
        assert abs(results.mean_log - ref["mean_log"]) < 0.01

    def test_at_site_skew_matches(self, cains_coulee):
        """Was 0.122 off; now 0.0002 -- the at-site EMA moment-iteration fix.

        ``_compute_ema_moments``'s censored branch used its own approximate
        truncated-gamma-moment formula instead of the already-ported,
        Fortran-verified ``hydrolib._p3_moments.m_p3``, and
        ``_ema_iteration``'s bias-correction factors used the total
        interval count where the vendored Fortran's actual default
        (``bcf=1997``, ``emafit.f:1408``) uses the exact-peak count. Both
        are now fixed; see
        ``tests/fortran_parity/test_fortran_oracles.py::TestMomentIterationOracle``
        for the routine-level version of this same measurement.
        """
        results, ref = cains_coulee
        assert abs(results.skew_station - ref["skew_at_site"]) < 0.02

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "skew_weighted still carries a real 0.058-skew-unit gap. Not "
            "the at-site fit (skew_station is correct to 0.0002), not "
            "detrat (native Wd is 0.186 against peakfq's 0.184), and -- "
            "measured directly, not assumed -- not mse_ema/var_mom/"
            "mn2mvarb either: called standalone with Cains Coulee's real "
            "post-MGBT group, mse_ema(kmom=3) matches the Fortran oracle "
            "to 3e-8 relative. The gap is that emafitpr's own internally-"
            "computed as_G_mse for this site (0.2212, what the golden "
            "skew_weighted was built from) does not match what calling "
            "the same mseg_all Fortran routine standalone gives for the "
            "identical (nobs, tl, tu, mc) -- 0.0749, reproduced even by a "
            "from-scratch, single-case golden regeneration, ruling out "
            "cross-case SAVE state contamination as the cause this time. "
            "The exact mechanism (likely something in emafitpr's own "
            "multi-stage internal fitting before the reported value is "
            "set) was not pinned down despite substantial investigation. "
            "hydrolib's own composition of independently-verified routines "
            "has no principled way to reproduce an unexplained number, so "
            "it is left as is. See TODO.md P3."
        ),
    )
    def test_weighted_skew_matches(self, cains_coulee):
        results, ref = cains_coulee
        assert abs(results.skew_weighted - ref["skew_weighted"]) < 0.02

    def test_quantile_error_is_bounded_and_worst_in_the_lower_tail(self, cains_coulee):
        """Recorded, not asserted away: 0.08% at best, 9.7% at worst, 1.5% at Q100.

        Was 0.64%-23.9% (6.3% at Q100) before the at-site EMA moment
        iteration was fixed (TODO.md P3); the remaining error traces
        entirely to the one open item, ``skew_weighted``'s own 0.058-skew-unit
        gap (``test_weighted_skew_matches``) -- everything upstream of it
        (the at-site fit, ``Wd``, ADJE) now matches peakfq 8.1.0 closely.
        The lower tail is where that residual bites hardest, which is the
        signature of a skew-weighting gap rather than a broken fit.
        """
        results, ref = cains_coulee
        errors = _quantile_errors(results, ref)
        assert max(errors.values()) < 12.0
        assert errors[0.01] < 2.0
        worst_aep = max(errors, key=errors.get)
        assert worst_aep > 0.5, f"worst error at AEP {worst_aep}, expected the lower tail"
