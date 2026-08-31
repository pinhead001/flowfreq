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
    0.0066. hydrolib has no ``detrat``, so it uses 1.

Read together they localise the open P3 defect precisely: with no censoring
the native fit is exact, and everything that remains is censored-path work --
the ADJE bias adjustment on the at-site skew MSE, and ``detrat`` itself.

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
    What does not match is what the censoring then does to the moments, and
    the reference ``Wd`` of 0.184 is a large part of why -- hydrolib has no
    ``detrat`` and uses 1, so it weights the regional skew far more heavily
    than peakfq does here.
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

    def test_mean_still_agrees_closely(self, cains_coulee):
        """Censoring moves the mean least; measured 5.9e-3 in log10, ~1.4% in flow."""
        results, ref = cains_coulee
        assert abs(results.mean_log - ref["mean_log"]) < 0.01

    @pytest.mark.xfail(
        strict=True,
        reason=(
            "Censored-path defect. peakfq's at-site skew MSE uses the ADJE "
            "censoring bias adjustment (from var_mom, not ported) and its Wd "
            "comes from detrat (also not ported, and 0.184 here against the 1 "
            "hydrolib assumes). Measured gap 0.122 on the at-site skew and "
            "0.098 on the weighted skew. See TODO.md P3."
        ),
    )
    @pytest.mark.parametrize("field", ["skew_at_site", "skew_weighted"])
    def test_skew_matches(self, cains_coulee, field):
        results, ref = cains_coulee
        native = {
            "skew_at_site": results.skew_station,
            "skew_weighted": results.skew_weighted,
        }[field]
        assert abs(native - ref[field]) < 0.02

    def test_quantile_error_is_bounded_and_worst_in_the_lower_tail(self, cains_coulee):
        """Recorded, not asserted away: 0.30% at best, 18.7% at worst, 2.7% at Q100.

        The lower tail is where the censoring bites, which is the signature of
        the open defect rather than of a broken fit.
        """
        results, ref = cains_coulee
        errors = _quantile_errors(results, ref)
        assert max(errors.values()) < 25.0
        assert errors[0.01] < 5.0
        worst_aep = max(errors, key=errors.get)
        assert worst_aep > 0.5, f"worst error at AEP {worst_aep}, expected the lower tail"
