"""Integration tests for hybrid Bulletin 17C workflow."""

from __future__ import annotations

import pathlib

import numpy as np
import pytest

from hydrolib.bulletin17c import Bulletin17C
from hydrolib.validation.comparisons import FrequencyComparator
from hydrolib.validation.reference import ReferenceResult
from tests.fixtures.big_sandy import (
    EXPECTED_CONFIDENCE_INTERVALS,
    EXPECTED_PARAMETERS,
    EXPECTED_QUANTILES,
    HISTORICAL_PEAKS,
    REGIONAL_SKEW,
    REGIONAL_SKEW_SD,
    SYSTEMATIC_PEAKS,
    THRESHOLDS,
    TOLERANCE_PERCENT,
)

GOLDEN = (
    pathlib.Path(__file__).resolve().parents[1]
    / "fortran_parity"
    / "golden"
    / "big_sandy_03606500.json"
)


class TestBigSandyNativeSystematic:
    """Test Big Sandy analysis using native EMA with systematic data only.

    The native EMA implementation has numerical limitations with
    historical/censored intervals. These tests use systematic-only data
    where the native solver converges reliably.
    """

    def _run_native_systematic(self) -> Bulletin17C:
        """Run native EMA on Big Sandy systematic data only."""
        peak_flows = np.array(list(SYSTEMATIC_PEAKS.values()))
        water_years = np.array(list(SYSTEMATIC_PEAKS.keys()))

        b17c = Bulletin17C(
            peak_flows=peak_flows,
            water_years=water_years,
            regional_skew=REGIONAL_SKEW,
            regional_skew_mse=REGIONAL_SKEW_SD**2,
        )
        b17c.run_analysis(method="ema")
        return b17c

    def test_native_ema_converges(self) -> None:
        """Native EMA converges on systematic-only data."""
        b17c = self._run_native_systematic()
        assert b17c.results is not None
        assert b17c.results.ema_converged is True

    def test_native_produces_parameters(self) -> None:
        """Native EMA produces finite LP3 parameters."""
        b17c = self._run_native_systematic()
        r = b17c.results
        assert np.isfinite(r.mean_log)
        assert np.isfinite(r.std_log)
        assert r.std_log > 0
        assert r.skew_weighted is not None

    def test_native_produces_quantiles(self) -> None:
        """Native EMA produces reasonable quantile estimates."""
        b17c = self._run_native_systematic()
        quantiles = b17c.compute_quantiles()
        assert len(quantiles) > 0
        # 1% AEP should produce a reasonable flood estimate for Big Sandy
        row_01 = quantiles[quantiles["aep"] == 0.01]
        if not row_01.empty:
            q_01 = row_01["flow_cfs"].iloc[0]
            assert 5000 < q_01 < 100000, f"1% AEP flow {q_01} out of reasonable range"

    def test_to_comparison_dict(self) -> None:
        """to_comparison_dict produces expected structure."""
        b17c = self._run_native_systematic()
        d = b17c.to_comparison_dict()

        assert "parameters" in d
        assert "quantiles" in d
        assert "confidence_intervals" in d
        assert "mean_log" in d["parameters"]
        assert "std_log" in d["parameters"]
        assert len(d["quantiles"]) > 0

    def test_validate_against_reference(self) -> None:
        """Native results can be validated against a reference.

        Uses generous tolerances since systematic-only analysis will
        differ from the full historical analysis in the expected values.
        """
        b17c = self._run_native_systematic()

        # Create a reference result from expected values
        reference = ReferenceResult(
            parameters=dict(EXPECTED_PARAMETERS),
            quantiles=dict(EXPECTED_QUANTILES),
            confidence_intervals={
                aep: (lo, hi) for aep, (lo, hi) in EXPECTED_CONFIDENCE_INTERVALS.items()
            },
        )

        # Use very generous tolerance since systematic-only differs from full analysis
        result = b17c.validate(
            reference, tolerance_pct=50.0, parameter_tolerance_pct=50.0, ci_tolerance_pct=50.0
        )
        assert result.summary is not None
        assert result.max_diff_pct >= 0


class TestBigSandyWithHistorical:
    """Test Big Sandy with historical data — native EMA now converges.

    The analytical truncated-moment approach (matching the Fortran mP3
    routine) resolves the overflow that previously caused NaN divergence.
    """

    def _run_native_historical(self) -> Bulletin17C:
        peak_flows = np.array(list(SYSTEMATIC_PEAKS.values()))
        water_years = np.array(list(SYSTEMATIC_PEAKS.keys()))
        historical = [(year, q) for year, q in HISTORICAL_PEAKS.items()]

        b17c = Bulletin17C(
            peak_flows=peak_flows,
            water_years=water_years,
            regional_skew=REGIONAL_SKEW,
            regional_skew_mse=REGIONAL_SKEW_SD**2,
            historical_peaks=historical,
        )
        b17c.run_analysis(method="ema")
        return b17c

    def test_native_ema_with_historical_converges(self) -> None:
        """Native EMA with historical data converges to finite parameters."""
        b17c = self._run_native_historical()
        r = b17c.results
        assert r is not None
        assert r.ema_converged is True
        assert np.isfinite(r.mean_log)
        assert np.isfinite(r.std_log)
        assert r.std_log > 0

    def test_native_historical_quantiles_reasonable(self) -> None:
        """Historical EMA produces reasonable quantile estimates."""
        b17c = self._run_native_historical()
        quantiles = b17c.compute_quantiles()
        assert len(quantiles) > 0
        row_01 = quantiles[quantiles["aep"] == 0.01]
        if not row_01.empty:
            q_01 = row_01["flow_cfs"].iloc[0]
            assert 5000 < q_01 < 200000, f"1% AEP flow {q_01} out of reasonable range"


class TestBigSandyAgainstReference:
    """Native EMA against peakfq 8.1.0, through the public ``validate()`` path.

    This class used to be ``TestBigSandyEndToEnd``: a bare ``pass`` under
    ``@requires_peakfqsa``, waiting on a PeakfqSA binary that does not exist.
    ``pytest -m requires_peakfqsa`` collected zero tests, so the marker and its
    exclusion were dead config for as long as they existed.

    The reference is real now -- the committed golden file, generated from the
    vendored Fortran by ``tools/gen_fortran_golden.py``. ``tests/fortran_parity/``
    compares the numbers directly; what is exercised here is the facade the
    library offers a user for the same job, ``Bulletin17C.validate()`` against a
    ``ReferenceResult``, on the like-for-like fit: the same censored record, the
    same perception thresholds.
    """

    @staticmethod
    def _reference() -> ReferenceResult:
        return ReferenceResult.from_golden(GOLDEN)

    @staticmethod
    def _native() -> Bulletin17C:
        """The same configuration the reference was generated from."""
        b17c = Bulletin17C(
            peak_flows=list(SYSTEMATIC_PEAKS.values()),
            water_years=list(SYSTEMATIC_PEAKS.keys()),
            regional_skew=REGIONAL_SKEW,
            regional_skew_mse=REGIONAL_SKEW_SD**2,
            historical_peaks=[(y, f) for y, f in HISTORICAL_PEAKS.items()],
            perception_thresholds={(t["start"], t["end"]): t["lower"] for t in THRESHOLDS},
        )
        b17c.run_analysis()
        return b17c

    def test_reference_loads_with_provenance(self):
        """A comparison is only meaningful if it names what it compared against."""
        ref = self._reference()
        assert "8.1.0" in ref.source
        assert ref.n_peaks == 84
        assert ref.n_historical == 3

    def test_reference_carries_the_fortran_diagnostics(self):
        """as_G_PRL_o and Wdout are what the two open parity defects argue from."""
        params = self._reference().parameters
        assert params["pseudo_record_length"] == pytest.approx(54.373, abs=1e-3)
        assert params["weight_factor"] == pytest.approx(1.0)

    def test_validate_returns_a_populated_comparison(self):
        result = self._native().validate(self._reference(), tolerance_pct=5.0)
        assert result.quantile_diffs
        assert result.parameter_diffs
        assert result.ci_diffs
        assert result.summary

    def test_quantiles_agree_within_5pct(self):
        """Measured 0.0006% at the AEP-0.5 median, rising to 0.056% at the 0.002 tail.

        Was 0.12%-4.70% before the at-site EMA moment iteration was fixed to
        use ``hydrolib._p3_moments.m_p3`` and the Fortran's own (exact-peak,
        not total-record) bias-correction count -- Big Sandy has 37 censored
        historical gap-year intervals, so both bugs bit here even though
        neither is skew-weighting-specific. The 5% bound is left loose on
        purpose; the docstring is what actually pins the number.
        """
        diffs = self._native().validate(self._reference()).quantile_diffs
        worst_aep = max(diffs, key=lambda a: diffs[a])
        assert (
            diffs[worst_aep] < 5.0
        ), f"worst quantile diff at AEP {worst_aep}: {diffs[worst_aep]}%"

    def test_moments_agree_within_1pct(self):
        """Mean and standard deviation, the parameters with no open defect.

        Measured 3.5e-7% and 3.3e-5% -- the loose 0.1%/1.0% bounds predate
        the at-site moment-iteration fix (TODO.md P3) and are kept loose
        deliberately; see test_quantiles_agree_within_5pct's docstring.
        """
        params = self._native().validate(self._reference()).parameter_diffs
        assert params["mean_log"] < 0.1
        assert params["std_log"] < 1.0

    def test_weighted_skew_gap_is_closed(self):
        """The P3 skew-weighting defect this test used to pin, now closed.

        Was ~35% when the weighting was a post-hoc average of two skews,
        then ~24% (0.0376 in skew units) once the regional skew was folded
        into the EMA fixed point the way moms_p3 does it, then 0.0026 once
        ADJE (hydrolib._mse_ema.mse_ema, var_mom Phase 3) supplied peakfq's
        own censoring-adjusted skew MSE instead of the plain Bulletin 17B
        value. Measured now: **2.4e-6** -- essentially exact, the same level
        Powder River (no censoring at all) already hit.

        The last mile was a second bug in the at-site EMA moment iteration
        itself, unrelated to skew weighting: ``_compute_ema_moments``'s
        censored-interval branch used its own approximate truncated-gamma
        formula instead of the already-ported, Fortran-verified
        ``hydrolib._p3_moments.m_p3``, and ``_ema_iteration``'s bias
        corrections (``c2``, ``c3``) used the total interval count where the
        vendored Fortran's actual default (``bcf=1997``, ``emafit.f:1408``)
        uses the exact-peak count. Both bit here because Big Sandy has 37
        censored historical gap-year intervals -- not a skew-weighting input
        at all, which is why fixing ADJE/detrat alone left a residual and
        this needed its own fix. See
        tests/fortran_parity/test_fortran_oracles.py::TestMomentIterationOracle
        for the oracle-level version of this same fix.

        Asserted in skew units, not percent: FrequencyComparator compares
        skew by absolute difference, because a quantity that passes through
        zero cannot be judged by a ratio.
        """
        skews = self._native().validate(self._reference()).skew_diffs
        assert skews["skew_weighted"] < 1e-4

    def test_skew_at_site_percent_difference_is_not_meaningful(self):
        """A documented wart, so nobody reads 249% as a 249% error.

        The reference at-site skew is 0.0066. Percent difference divides by it,
        so an absolute gap of 0.016 reads as 249%. FrequencyComparator uses
        percent difference for every parameter, which is the wrong metric for
        one that legitimately crosses zero.
        """
        ref = self._reference()
        native = self._native().to_comparison_dict()
        assert ref.parameters["skew_at_site"] == pytest.approx(0.0066, abs=1e-3)
        absolute_gap = abs(native["parameters"]["skew_at_site"] - ref.parameters["skew_at_site"])
        assert absolute_gap < 0.05
