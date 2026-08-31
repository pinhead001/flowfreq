"""Check the committed golden files still match the vendored Fortran.

Skipped wherever the f2py extension is absent. Where it is present -- the CI
``fortran`` job, or a developer machine that has run ``build_fortran/build.py``
-- this is the guard that stops a golden file from silently drifting away from
the sources that produced it. If ``vendor/peakfqr`` changes and nobody
regenerates, this is what says so.

Tolerances
----------
This file used to assert ``rtol=1e-9`` on the reasoning that it is the same
code against the same inputs, so anything larger means the reference moved.
That reasoning is wrong, and the first CI run of the ``fortran`` job proved it:
eight tests failed on a runner with *identical* gfortran (13.3.0), numpy
(2.4.6) and meson (1.12.0) versions, against an untouched ``vendor/peakfqr``.

The EMA fit is a fixed-point iteration, and it is ill-conditioned in the third
moment. Perturbing one input by exactly one ulp -- a relative change of
2.2e-16, the smallest a float64 can carry -- moves the converged outputs by:

    mean          3.9e-08        yp[0]         1.1e-06
    skew_w        6.7e-05        as_G_mse      9.4e-05
    as_G_PRL      1.1e-04        skew_at_site  3.2e-03

The differences the CI runner actually showed were 1.1e-07 to 1.6e-02, in the
same per-quantity order -- a few ulps' worth. So bit-agreement between two
machines is not available at any tolerance near 1e-9, whatever the toolchain,
and asserting it just makes the job red on arrival.

The tolerances below are therefore set per quantity, from what each one is
rather than from one global number: roughly thirty ulps of headroom over the
measured conditioning, which is still three to five orders of magnitude tighter
than any change to the Fortran that anyone would care about. Skews get an
absolute tolerance because the at-site skew is 0.0066 on this site -- near
enough zero that a relative test on it measures nothing.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip(
    "hydrolib.peakfqr",
    reason="Fortran extension not built; run python build_fortran/build.py "
    "(see docs/FORTRAN_UPLOAD.md)",
)

pytestmark = pytest.mark.requires_fortran

# log10(cfs). 1e-4 in log10 is 0.023% in discharge; the runner showed 1.6e-5.
ATOL_LOG10 = 1e-4
# Dimensionless, and legitimately close to zero, so absolute rather than relative.
ATOL_SKEW = 1e-3
# Central moments: means are O(3.7), variances O(0.08).
ATOL_MEAN = 1e-4
ATOL_VARIANCE = 1e-5
# MSEs, pseudo record length, Wd, and per-quantile variances -- all well away
# from zero, so a relative tolerance is the meaningful one here.
RTOL_DIAGNOSTIC = 1e-3


@pytest.fixture(scope="module")
def live_big_sandy():
    from tests.fortran_parity.cases import big_sandy_case, call_emafitpr

    return call_emafitpr(big_sandy_case())


def test_version_still_matches(golden_big_sandy):
    """The golden file names the peakfq version it came from; vendor/ must agree."""
    from tools.gen_fortran_golden import peakfq_version

    assert golden_big_sandy["meta"]["peakfq_version"] == peakfq_version(), (
        "vendor/peakfqr/DESCRIPTION reports a different version than the golden file "
        "records -- regenerate with tools/gen_fortran_golden.py"
    )


def test_interval_count(golden_big_sandy, live_big_sandy):
    assert live_big_sandy["n"] == golden_big_sandy["outputs"]["n"]


def test_mgbt_unchanged(golden_big_sandy, live_big_sandy):
    assert live_big_sandy["mgbt"] == golden_big_sandy["outputs"]["mgbt"]


@pytest.mark.parametrize("key", ["as_G_mse_o", "as_G_mse_Syst_o", "as_G_PRL_o", "Wdout"])
def test_skew_diagnostics_unchanged(golden_big_sandy, live_big_sandy, key):
    assert live_big_sandy["skew"][key] == pytest.approx(
        golden_big_sandy["outputs"]["skew"][key], rel=RTOL_DIAGNOSTIC
    )


# cmoms rows are mean, variance and skew, so one tolerance cannot serve all
# three: the row that varies most across machines is also the one nearest zero.
@pytest.mark.parametrize(
    "row, name, atol",
    [(0, "mean", ATOL_MEAN), (1, "variance", ATOL_VARIANCE), (2, "skew", ATOL_SKEW)],
)
def test_moments_unchanged(golden_big_sandy, live_big_sandy, row, name, atol):
    live = np.asarray(live_big_sandy["cmoms"], dtype=float)[row]
    golden = np.asarray(golden_big_sandy["outputs"]["cmoms"], dtype=float)[row]
    assert np.allclose(live, golden, rtol=0.0, atol=atol), (
        f"cmoms {name} row drifted from the golden file by "
        f"{np.max(np.abs(live - golden)):.3e}, past {atol:.0e}"
    )


@pytest.mark.parametrize(
    "key, rtol, atol",
    [
        ("yp", 0.0, ATOL_LOG10),
        ("ci_low", 0.0, ATOL_LOG10),
        ("ci_high", 0.0, ATOL_LOG10),
        ("var_est", RTOL_DIAGNOSTIC, 0.0),
    ],
)
def test_quantile_vectors_unchanged(golden_big_sandy, live_big_sandy, key, rtol, atol):
    live = np.asarray(live_big_sandy["quantiles"][key], dtype=float)
    golden = np.asarray(golden_big_sandy["outputs"]["quantiles"][key], dtype=float)
    assert np.allclose(live, golden, rtol=rtol, atol=atol), (
        f"{key} drifted from the golden file by " f"{np.max(np.abs(live - golden)):.3e}"
    )
