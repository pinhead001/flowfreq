"""Check the committed golden files still match the vendored Fortran.

Skipped wherever the f2py extension is absent, which includes CI. Where it is
present -- a developer machine that has run ``build_fortran/build.py`` -- this
is the guard that stops a golden file from silently drifting away from the
sources that produced it. If ``vendor/peakfqr`` changes and nobody regenerates,
this is what says so.

Tolerances are numerical noise only (1e-9 relative). This is the same code
against the same inputs; anything larger means the reference moved.
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

RTOL = 1e-9


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
        golden_big_sandy["outputs"]["skew"][key], rel=RTOL
    )


def test_moments_unchanged(golden_big_sandy, live_big_sandy):
    assert np.allclose(
        np.asarray(live_big_sandy["cmoms"]),
        np.asarray(golden_big_sandy["outputs"]["cmoms"]),
        rtol=RTOL,
    )


@pytest.mark.parametrize("key", ["yp", "ci_low", "ci_high", "var_est"])
def test_quantile_vectors_unchanged(golden_big_sandy, live_big_sandy, key):
    assert np.allclose(
        np.asarray(live_big_sandy["quantiles"][key]),
        np.asarray(golden_big_sandy["outputs"]["quantiles"][key]),
        rtol=RTOL,
    )
