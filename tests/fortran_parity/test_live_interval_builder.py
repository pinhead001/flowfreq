"""Live verification of ``flowfreq.fortran_engine.build_emafit_arrays``.

Everything in ``test_interval_builder.py`` runs without the built extension,
checking the builder's arrays against golden files and the existing
test-code builder. This module is the missing piece design doc section 7
calls for: running those same arrays -- built by the *library* builder, not
``tests/fortran_parity/cases.py::build_emafit_inputs`` -- through a live
``emafitpr`` call, and checking the result against the committed goldens.

Same cross-machine caveat as ``test_live_vs_golden.py``: the EMA fixed point
is ill-conditioned in the third moment, so bit-identical agreement across
different gfortran builds is not expected even from unchanged sources (that
file's own docstring measures a few ulps of input perturbation moving
``skew_at_site`` by 3.2e-3). Tolerances here are the same ones, not a looser
set chosen for convenience.

**Naming, not cosmetic**: this file (and ``test_live_compare_engines.py``,
``test_live_cli_compare.py``) calls ``emafitpr`` on several different cases
in one process. ``test_fortran_oracles.py``'s module docstring documents a
real ``mseg_all_sub`` ``SAVE``-state leak across calls to different
``emafitpr``-family entry points -- confirmed the hard way while adding this
file: named ``test_interval_builder_live.py`` it sorted *before*
``test_fortran_oracles.py`` alphabetically, and its multi-case calls left
``TestCainsCouleeAsGMseDiscrepancy``/``TestSkewMseOracle`` reading
contaminated state (0.2212 instead of 0.0749) even though nothing here calls
``mseg_all_sub`` directly. Renamed with a ``test_live_`` prefix so it sorts
after ``test_fortran_oracles.py``, the same place ``test_live_vs_golden.py``/
``test_native_vs_golden.py``/``test_wymt_vs_golden.py`` already safely sit.
Do not rename a live-Fortran, multi-case test file to sort before
``test_fortran_oracles.py`` without re-running the whole ``tests/fortran_parity/``
directory (not just the new file) to check for this.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.fortran_parity.conftest import load_golden

pytest.importorskip(
    "flowfreq.peakfqr",
    reason="Fortran extension not built; run python build_fortran/build.py "
    "(see docs/FORTRAN_UPLOAD.md)",
)

pytestmark = pytest.mark.requires_fortran

# Same calibration as test_live_vs_golden.py -- see that file's docstring for
# the cross-machine measurements this is based on.
ATOL_LOG10 = 1e-4
ATOL_SKEW = 1e-3
ATOL_MEAN = 1e-4
ATOL_VARIANCE = 1e-5
RTOL_DIAGNOSTIC = 1e-3


def _case_names():
    from tests.fortran_parity.cases import CASES

    return sorted(CASES)


@pytest.fixture(scope="module", params=_case_names())
def library_built_case(request):
    """(golden document, live emafitpr output via the *library* builder) for one case.

    Deliberately does not import ``build_emafit_inputs`` -- the point is to
    exercise ``flowfreq.fortran_engine.build_emafit_arrays`` on the same
    record and confirm it reaches the same live Fortran answer, not to
    re-check the test-code builder (``test_live_vs_golden.py`` already does
    that).
    """
    from flowfreq.fortran_engine import run_fortran_reference

    name = request.param
    golden = load_golden(name)
    if golden is None:
        pytest.skip(f"golden file missing for {name} (run tools/gen_fortran_golden.py)")

    try:
        peak_flows, water_years, historical_peaks, perception_thresholds = _record_for(name)
    except FileNotFoundError as exc:  # reference test data absent
        pytest.skip(f"{name}: {exc}")

    inputs = golden["inputs"]
    reference, arrays = run_fortran_reference(
        peak_flows=peak_flows,
        water_years=water_years,
        historical_peaks=historical_peaks,
        perception_thresholds=perception_thresholds,
        regional_skew=inputs["regional_skew"],
        regional_skew_mse=inputs["regional_skew_mse"],
        aeps=inputs["aeps"],
        eps=inputs["eps"],
        weight_opt=inputs["weight_opt"],
    )
    return name, golden, reference, arrays


def _record_for(name: str):
    """Raw ``Bulletin17C``-shaped inputs for a registered case, from its fixture directly.

    Returns
    -------
    tuple of (peak_flows, water_years, historical_peaks, perception_thresholds)
    """
    if name == "big_sandy_03606500":
        from tests.fixtures.big_sandy import HISTORICAL_PEAKS, SYSTEMATIC_PEAKS, THRESHOLDS

        years = sorted(SYSTEMATIC_PEAKS)
        perception_thresholds = {(t["start"], t["end"]): t["lower"] for t in THRESHOLDS}
        return (
            [SYSTEMATIC_PEAKS[y] for y in years],
            years,
            list(HISTORICAL_PEAKS.items()),
            perception_thresholds,
        )

    if name == "site_12363000":
        from tests.fixtures.site_12363000 import SYSTEMATIC_PEAKS

        years = sorted(SYSTEMATIC_PEAKS)
        return [SYSTEMATIC_PEAKS[y] for y in years], years, None, None

    # The two wymt cases share the same shape: plain contiguous systematic
    # records, no historical peaks, no perception thresholds.
    site_no = {"powder_river_06326500": "06326500.00", "cains_coulee_06327450": "06327450.00"}[name]
    from tests.fixtures.wymt_peaks import load_site

    site = load_site(site_no)
    years = sorted(site.peaks)
    return [site.peaks[y] for y in years], years, None, None


def test_interval_count_matches(library_built_case):
    _, golden, reference, _arrays = library_built_case
    assert reference.n_peaks == golden["outputs"]["n"]


def test_mgbt_matches(library_built_case):
    """Counts and the threshold are discrete, so these must match exactly."""
    _, golden, reference, _arrays = library_built_case
    assert reference.low_outlier_count == golden["outputs"]["mgbt"]["gbnlow"]


@pytest.mark.parametrize(
    "row, row_name, atol",
    [(0, "mean", ATOL_MEAN), (1, "variance", ATOL_VARIANCE), (2, "skew", ATOL_SKEW)],
)
def test_moments_match(library_built_case, row, row_name, atol):
    name, golden, reference, _arrays = library_built_case
    golden_row = np.asarray(golden["outputs"]["cmoms"], dtype=float)[row]
    # column 0: with regional info, column 1: at-site
    live_row = [
        reference.parameters[
            "mean_log" if row == 0 else "std_log" if row == 1 else "skew_weighted"
        ],
        reference.parameters[
            "mean_log_at_site" if row == 0 else "std_log_at_site" if row == 1 else "skew_at_site"
        ],
    ]
    if row == 1:
        live_row = [v**2 for v in live_row]  # cmoms row 1 is variance, parameters store std
    assert abs(live_row[0] - golden_row[0]) <= atol, f"{name}: {row_name} column 0 drifted"
    assert abs(live_row[1] - golden_row[1]) <= atol, f"{name}: {row_name} column 1 drifted"


def test_quantiles_match(library_built_case):
    name, golden, reference, _arrays = library_built_case
    aeps = golden["inputs"]["aeps"]
    golden_yp = np.asarray(golden["outputs"]["quantiles"]["yp"], dtype=float)
    live_log_q = np.array([np.log10(reference.quantiles[a]) for a in aeps])
    assert np.allclose(live_log_q, golden_yp, rtol=0.0, atol=ATOL_LOG10), (
        f"{name}: quantiles drifted from the golden by "
        f"{np.max(np.abs(live_log_q - golden_yp)):.3e} (log10)"
    )


class TestGapYearSwitchLive:
    """The design doc's worked example (section 3), against a live emafitpr call.

    Site 12363000: 98 peaks over 102 water years, gap at 1924-1927.

    Without a declared threshold, the library builder omits the gap and
    reproduces the committed golden *exactly* -- this is the same 98-row
    construction the golden was generated from, just built through
    ``build_emafit_arrays`` instead of ``build_emafit_inputs``.

    With a threshold declared over the gap, at-site skew and Q100 land on
    the design doc's published +0.250 / 119,473 -- found empirically here
    (any threshold from roughly 1,000-10,000 cfs, comfortably below every
    real peak in the record, reproduces the same result to 4 significant
    figures), not assumed. This confirms the design doc's own worked example
    against a live run, not just its two summary numbers.
    """

    THRESHOLD = 5_000.0  # synthetic; see class docstring and test_interval_builder.py

    @pytest.fixture(scope="class")
    def record(self):
        from tests.fixtures.site_12363000 import REGIONAL_SKEW, REGIONAL_SKEW_SD, SYSTEMATIC_PEAKS

        years = sorted(SYSTEMATIC_PEAKS)
        return years, [SYSTEMATIC_PEAKS[y] for y in years], REGIONAL_SKEW, REGIONAL_SKEW_SD**2

    def test_omitted_reproduces_the_golden_exactly(self, record):
        from flowfreq.fortran_engine import run_fortran_reference

        years, flows, r_g, r_g_mse = record
        golden = load_golden("site_12363000")
        aeps = golden["inputs"]["aeps"]

        reference, arrays = run_fortran_reference(
            peak_flows=flows,
            water_years=years,
            regional_skew=r_g,
            regional_skew_mse=r_g_mse,
            aeps=aeps,
        )
        assert arrays.n == 98
        golden_skew_at_site = golden["outputs"]["cmoms"][2][1]
        assert reference.parameters["skew_at_site"] == pytest.approx(
            golden_skew_at_site, abs=ATOL_SKEW
        )
        assert golden_skew_at_site == pytest.approx(0.435, abs=5e-4)

    def test_censored_reproduces_the_published_split(self, record):
        from flowfreq.fortran_engine import run_fortran_reference

        years, flows, r_g, r_g_mse = record
        reference, arrays = run_fortran_reference(
            peak_flows=flows,
            water_years=years,
            perception_thresholds={(min(years), max(years)): self.THRESHOLD},
            regional_skew=r_g,
            regional_skew_mse=r_g_mse,
            aeps=[0.01],
        )
        assert arrays.n == 102
        assert arrays.n_censored == 4
        assert reference.parameters["skew_at_site"] == pytest.approx(0.250, abs=ATOL_SKEW)
        assert reference.quantiles[0.01] == pytest.approx(119_473.0, rel=RTOL_DIAGNOSTIC)
