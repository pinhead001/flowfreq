"""Compare hydrolib's native EMA against the committed Fortran output.

Runs everywhere, including CI, because the reference values come from a golden
file rather than from the extension. See ``tools/gen_fortran_golden.py``.

Structured as the diagnostic ladder in ``docs/FORTRAN_UPLOAD.md`` section 6:
each rung feeds the next, so the *first* failure is the root cause and anything
below it is downstream noise. Two rungs are known-failing and marked
``xfail(strict=True)`` -- they are real, understood defects, and strict means
the build fails the moment either starts passing.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.fortran_parity.conftest import aep_index, asymmetry

_ASYM_REASON = (
    "compute_confidence_limits() forms log_Q +/- z*se, so its upper/lower half-width "
    "ratio is exactly 1.000 at every AEP. peakfq 8.1.0 runs 1.03 -> 1.31 -> 1.41 with "
    "return period via the Inverse Modified Cholesky Gaussian Quadrature, which hydrolib "
    "does not implement. See docs/FORTRAN_UPLOAD.md sections 1.6 and 6.0."
)


class TestRung1Censoring:
    """Are the two codes fitting the same data at all?"""

    def test_interval_count(self, golden_big_sandy, native_big_sandy):
        assert golden_big_sandy["outputs"]["n"] == native_big_sandy._results.n_peaks

    def test_no_low_outliers_detected(self, golden_big_sandy, native_big_sandy):
        """MGBT finds no PILFs in this record; the native path must agree."""
        assert golden_big_sandy["outputs"]["mgbt"]["gbnlow"] == 0
        assert native_big_sandy._results.n_low_outliers == 0


class TestRung3Moments:
    """Same mean, variance, skew on identical data?"""

    def test_mean_log(self, golden_big_sandy, native_big_sandy):
        expected = golden_big_sandy["outputs"]["cmoms"][0][0]
        actual = native_big_sandy._results.mean_log
        assert abs(actual - expected) / abs(expected) * 100 < 0.1

    def test_std_log(self, golden_big_sandy, native_big_sandy):
        expected = golden_big_sandy["outputs"]["cmoms"][1][0] ** 0.5
        actual = native_big_sandy._results.std_log
        assert abs(actual - expected) / abs(expected) * 100 < 1.5

    def test_at_site_skew(self, golden_big_sandy, native_big_sandy):
        """Compared absolutely: both sit near zero, so a ratio is meaningless."""
        expected = golden_big_sandy["outputs"]["cmoms"][2][1]
        actual = native_big_sandy._results.skew_station
        assert abs(actual - expected) < 0.05

    def test_weighted_skew(self, golden_big_sandy, native_big_sandy):
        """No longer xfail: TODO.md P3's var_mom port closes this.

        Was -0.1563 vs -0.1009, ~35%, when hydrolib used the standard
        Bulletin 17C weighting instead of peakfq's HWN "optimized adjustment
        factor when censored data are present" (fortranWrappers.R). Folding
        the regional skew into the EMA fixed point the way moms_p3 does it
        closed most of that gap; hydrolib._mse_ema.mse_ema (var_mom's
        censoring bias adjustment, ADJE's as_G_mse) closed the rest. See
        docs/FORTRAN_UPLOAD.md section 6.0b for the history and
        tests/integration/test_hybrid_workflow.py::TestBigSandyAgainstReference
        .test_weighted_skew_gap_is_closed for the current measured gap.
        """
        expected = golden_big_sandy["outputs"]["cmoms"][2][0]
        actual = native_big_sandy._results.skew_used
        assert abs(actual - expected) < 0.02


class TestRung5Quantiles:
    """Same quantiles, given that the moments differ slightly?"""

    @pytest.mark.parametrize("aep", [0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002])
    def test_quantile(self, golden_big_sandy, native_big_sandy, aep):
        i = aep_index(golden_big_sandy, aep)
        expected = 10 ** golden_big_sandy["outputs"]["quantiles"]["yp"][i]
        actual = float(native_big_sandy.compute_quantiles(np.array([aep]))["flow_cfs"].iloc[0])
        diff = abs(actual - expected) / expected * 100
        assert diff < 3.0, f"Q(AEP={aep}) differs by {diff:.2f}%"


class TestRung6ConfidenceIntervals:
    """The open defect: interval shape, not size."""

    # Only the rare events are marked. At AEP 0.1 the reference itself is nearly
    # symmetric (1.032), so log_Q +/- z*se is adequate there and the assertion
    # genuinely passes -- which is the useful finding: the defect bites at long
    # return periods, not everywhere. Marking 0.1 too would hide that.
    @pytest.mark.parametrize(
        "aep",
        [
            0.1,
            pytest.param(0.02, marks=pytest.mark.xfail(strict=True, reason=_ASYM_REASON)),
            pytest.param(0.01, marks=pytest.mark.xfail(strict=True, reason=_ASYM_REASON)),
        ],
    )
    def test_interval_asymmetry(self, golden_big_sandy, native_big_sandy, aep):
        i = aep_index(golden_big_sandy, aep)
        q = golden_big_sandy["outputs"]["quantiles"]
        expected = asymmetry(10 ** q["ci_low"][i], 10 ** q["yp"][i], 10 ** q["ci_high"][i])
        ci = native_big_sandy.compute_confidence_limits(np.array([aep]))
        actual = asymmetry(
            float(ci["lower_5pct"].iloc[0]),
            float(ci["flow_cfs"].iloc[0]),
            float(ci["upper_5pct"].iloc[0]),
        )
        assert abs(actual - expected) < 0.05

    def test_native_intervals_are_exactly_symmetric(self, native_big_sandy):
        """Pin the current behaviour, so the defect is documented rather than implied."""
        for aep in (0.1, 0.02, 0.01):
            ci = native_big_sandy.compute_confidence_limits(np.array([aep]))
            ratio = asymmetry(
                float(ci["lower_5pct"].iloc[0]),
                float(ci["flow_cfs"].iloc[0]),
                float(ci["upper_5pct"].iloc[0]),
            )
            assert abs(ratio - 1.0) < 1e-9


class TestGoldenProvenance:
    """A golden file that cannot say where it came from is not evidence."""

    def test_records_its_source(self, golden_big_sandy):
        meta = golden_big_sandy["meta"]
        assert meta["peakfq_version"] == "8.1.0"
        assert meta["site_no"] == "03606500"
        for key in ("generated_utc", "vendor_src_commit", "generator", "units"):
            assert meta.get(key), f"golden file is missing meta.{key}"

    def test_inputs_are_recorded_beside_outputs(self, golden_big_sandy):
        """Without its inputs a golden file is unfalsifiable and unreproducible."""
        inputs = golden_big_sandy["inputs"]
        n = golden_big_sandy["outputs"]["n"]
        for key in ("ql", "qu", "tl", "tu", "dtype"):
            assert (
                len(inputs[key]) == n
            ), f"inputs.{key} has {len(inputs[key])} entries, expected {n}"
        assert len(inputs["aeps"]) == len(golden_big_sandy["outputs"]["quantiles"]["yp"])
