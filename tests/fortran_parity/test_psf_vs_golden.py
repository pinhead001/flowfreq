"""The .psf converter against the committed peakfq 8.1.0 goldens (issue #31).

End to end: parse ``wymt_ffa_2022A.psf`` and its WATSTORE file, convert a
station (:mod:`flowfreq.psf_convert`), and compare with the golden files
``tools/gen_fortran_golden.py`` generated for Powder River (06326500) and
Cains Coulee (06327450). Runs without the Fortran extension.

Two things are checked separately, because the goldens were not generated from
the ``.psf``:

1. The converter's ``emafitpr`` rows and low-outlier/weighting/confidence
   arguments equal the goldens' inputs exactly.
2. The regional skew does not quite: the goldens took ``RegSkew``/``RegMSEG``
   from peakfq 7.4's ``EXPinfo`` CSV, which rounds them. Powder River's
   ``GenSkew -0.2002`` is -0.2 there, and every station's ``SkewSE 0.64``
   gives an MSE of 0.4096 (``main.R`` squares it) that the CSV prints as 0.41.
   The converter follows ``main.R`` and keeps the ``.psf`` values; the tests
   below assert the gap is exactly that rounding, then fit with the goldens'
   rounded values substituted, so the comparison measures the conversion and
   not the rounding. ``test_live_psf_convert.py`` compares against a live
   Fortran call with the ``.psf``'s own values instead.

Tolerances are the ones ``test_wymt_vs_golden.py`` already uses for these two
sites. Cains Coulee's weighted skew is not re-tested here: it carries the
known ``xfail(strict=True)`` in ``test_wymt_vs_golden.py``, for reasons that
have nothing to do with how the inputs are built (they are identical, as
``TestInputsMatchGolden`` shows).
"""

from __future__ import annotations

from collections import Counter

import numpy as np
import pytest

from tests.fixtures.paths import SKIP_REASON, TESTDATA_AVAILABLE
from tests.fixtures.paths import testdata_path as _testdata
from tests.fortran_parity.conftest import load_golden

pytestmark = [
    pytest.mark.requires_peakfqr_testdata,
    pytest.mark.skipif(not TESTDATA_AVAILABLE, reason=SKIP_REASON),
]

SITES = {
    "powder_river_06326500": "06326500.00",
    "cains_coulee_06327450": "06327450.00",
}


@pytest.fixture(scope="module")
def stations():
    from flowfreq.psf_convert import convert_psf

    return convert_psf(_testdata("wymt_ffa_2022A.psf"))


def _golden(name):
    golden = load_golden(name)
    if golden is None:
        pytest.skip(f"golden file missing for {name} (run tools/gen_fortran_golden.py)")
    return golden


def _rounded_rows(ql, qu, tl, tu, dtype):
    """(ql, qu, tl, tu, dtype) rows with the log10 bounds rounded to 10 decimals."""
    return [
        (round(a, 10), round(b, 10), round(c, 10), round(d, 10), int(e))
        for a, b, c, d, e in zip(ql, qu, tl, tu, dtype)
    ]


@pytest.mark.parametrize("name", sorted(SITES))
class TestInputsMatchGolden:
    def test_emafit_rows_identical(self, stations, name):
        golden = _golden(name)["inputs"]
        arrays = stations[SITES[name]].emafit_arrays()
        # Rows are log10 flows. Compared to 10 decimals, not bit-for-bit: numpy's
        # log10 can differ from the golden's by an ulp across platforms (this
        # passed on Windows and failed on Linux CI), and 1e-10 in log space is
        # far below anything that would change which interval a year falls in.
        ours = Counter(
            _rounded_rows(
                arrays.ql.tolist(),
                arrays.qu.tolist(),
                arrays.tl.tolist(),
                arrays.tu.tolist(),
                arrays.dtype.tolist(),
            )
        )
        theirs = Counter(
            _rounded_rows(golden["ql"], golden["qu"], golden["tl"], golden["tu"], golden["dtype"])
        )
        assert ours == theirs

    def test_scalar_arguments_identical(self, stations, name):
        golden = _golden(name)["inputs"]
        s = stations[SITES[name]]
        assert s.emafit_arrays().gbthrsh0 == golden["gbthrsh0"]  # MGBT
        assert s.confidence == golden["eps"]
        assert {"HWN": 1, "ERL": 2, "INV": 3}[s.weight_option] == golden["weight_opt"]

    def test_regional_skew_differs_only_by_csv_rounding(self, stations, name):
        golden = _golden(name)["inputs"]
        s = stations[SITES[name]]
        assert s.skew_option == "weighted"
        assert s.regional_skew_mse == pytest.approx(0.64**2, abs=1e-15)  # SkewSE**2
        assert round(s.regional_skew, 3) == golden["regional_skew"]
        assert round(s.regional_skew_mse, 2) == golden["regional_skew_mse"]


def _fit_with_golden_skew(stations, name):
    from flowfreq.bulletin17c import Bulletin17C

    golden = _golden(name)
    kwargs = stations[SITES[name]].bulletin17c_kwargs("native")
    kwargs["regional_skew"] = golden["inputs"]["regional_skew"]
    kwargs["regional_skew_mse"] = golden["inputs"]["regional_skew_mse"]
    b17c = Bulletin17C(**kwargs)
    results = b17c.run_analysis(method="ema")
    aeps = np.asarray(golden["inputs"]["aeps"], dtype=float)
    quantiles = b17c.compute_quantiles(aep=aeps)
    cmoms = np.asarray(golden["outputs"]["cmoms"], dtype=float)
    ref_q = 10.0 ** np.asarray(golden["outputs"]["quantiles"]["yp"], dtype=float)
    q_err = np.abs(quantiles["flow_cfs"].to_numpy() / ref_q - 1) * 100
    return results, cmoms, q_err, golden


class TestPowderRiverFromPsf:
    """85 uncensored years: must match to the same 1e-6 as the hand-built fixture."""

    @pytest.fixture(scope="class")
    def fit(self, stations):
        return _fit_with_golden_skew(stations, "powder_river_06326500")

    def test_moments(self, fit):
        results, cmoms, _, _ = fit
        assert abs(results.mean_log - cmoms[0][0]) < 1e-6
        assert abs(results.std_log - np.sqrt(cmoms[1][0])) < 1e-6
        assert abs(results.skew_weighted - cmoms[2][0]) < 1e-6
        assert abs(results.skew_station - cmoms[2][1]) < 1e-6

    def test_no_pilfs(self, fit):
        results, _, _, golden = fit
        assert results.n_low_outliers == golden["outputs"]["mgbt"]["gbnlow"] == 0

    def test_quantiles(self, fit):
        _, _, q_err, _ = fit
        assert q_err.max() < 0.5


class TestCainsCouleeFromPsf:
    """32 years, 11 PILFs: the parts ``test_wymt_vs_golden.py`` shows agreeing."""

    @pytest.fixture(scope="class")
    def fit(self, stations):
        return _fit_with_golden_skew(stations, "cains_coulee_06327450")

    def test_mgbt(self, fit):
        results, _, _, golden = fit
        assert results.n_low_outliers == golden["outputs"]["mgbt"]["gbnlow"] == 11
        assert results.low_outlier_threshold == pytest.approx(332.0, rel=1e-3)

    def test_at_site_skew(self, fit):
        results, cmoms, _, _ = fit
        assert abs(results.skew_station - cmoms[2][1]) < 0.02

    def test_mean(self, fit):
        results, cmoms, _, _ = fit
        assert abs(results.mean_log - cmoms[0][0]) < 0.01
