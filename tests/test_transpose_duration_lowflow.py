"""Tests for duration and low-flow transposition, and the FDC they rest on.

The guardrail worth naming here is
``test_flood_exponents_are_refused_for_low_flow``: the design doc calls
applying a flood exponent to 7Q10 a category error rather than an
approximation, and ``probability_kind`` is what turns that from a sentence in
a document into a ``ValueError``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq.core import LowFlowResults
from flowfreq.regime import DEFAULT_EXCEEDANCE_PCT, flow_duration_curve
from flowfreq.transpose import (
    LOW_FLOW_AREA_RATIO_RANGE,
    RegressionExponents,
    transpose_duration,
    transpose_frequency,
    transpose_low_flow,
)

DURATION_CITATION = "Example duration regression, region 3"
LOW_FLOW_CITATION = "Example low-flow regression, region 3"


@pytest.fixture()
def daily() -> pd.DataFrame:
    """A synthetic daily record with a realistic spread."""
    rng = np.random.default_rng(0)
    index = pd.date_range("2000-10-01", periods=3653, freq="D")
    flows = 200.0 * 10.0 ** rng.normal(0.0, 0.35, size=len(index))
    return pd.DataFrame({"flow_cfs": flows}, index=index)


@pytest.fixture()
def duration_exponents() -> RegressionExponents:
    """b near 1 at the wet end, easing off toward the dry end."""
    return RegressionExponents(
        aeps=[0.01, 0.10, 0.50, 0.90, 0.99],
        exponents=[0.95, 0.92, 0.87, 0.78, 0.70],
        citation=DURATION_CITATION,
        probability_kind="exceedance",
    )


@pytest.fixture()
def low_flow_exponents() -> RegressionExponents:
    return RegressionExponents(
        aeps=[0.50, 0.20, 0.10, 0.04, 0.02],
        exponents=[0.88, 0.85, 0.83, 0.81, 0.80],
        citation=LOW_FLOW_CITATION,
        probability_kind="non_exceedance",
    )


def _low_flow_results(flows, probs=(0.5, 0.2, 0.1, 0.04, 0.02), p_zero=0.0):
    return LowFlowResults(
        n_years=40,
        n_zero_years=int(round(p_zero * 40)),
        p_zero=p_zero,
        n_day=7,
        year_type="climatic",
        distribution="lp3",
        mean_log=1.5,
        std_log=0.3,
        skew_station=-0.1,
        skew_used=-0.1,
        quantiles=pd.DataFrame(
            {
                "non_exceedance_prob": list(probs),
                "return_period": [1.0 / p for p in probs],
                "flow_cfs": list(flows),
            }
        ),
    )


class TestFlowDurationCurve:
    def test_default_percentiles(self, daily):
        curve = flow_duration_curve(daily)
        assert list(curve["exceedance_pct"]) == list(DEFAULT_EXCEEDANCE_PCT)
        assert len(curve) == len(DEFAULT_EXCEEDANCE_PCT)

    def test_matches_the_plots_long_standing_arithmetic(self, daily):
        """The extraction must not move any number the plot already reported."""
        flows = daily["flow_cfs"].dropna().values
        curve = flow_duration_curve(daily)
        for pct, value in zip(curve["exceedance_pct"], curve["flow_cfs"]):
            assert value == pytest.approx(np.percentile(flows, 100 - pct), rel=0, abs=0)

    def test_plot_table_still_matches_the_standalone_function(self, daily):
        """Hydrograph.plot_flow_duration_curve now delegates here; its table
        must be unchanged in both values and column names."""
        import matplotlib

        matplotlib.use("Agg")
        from flowfreq.hydrograph import Hydrograph

        fig, stats = Hydrograph.plot_flow_duration_curve(daily)
        curve = flow_duration_curve(daily)
        assert list(stats.columns) == ["Percent Exceeded", "Flow (cfs)"]
        assert list(stats["Percent Exceeded"]) == [f"{p:g}%" for p in DEFAULT_EXCEEDANCE_PCT]
        np.testing.assert_allclose(stats["Flow (cfs)"], curve["flow_cfs"], rtol=0, atol=0)
        matplotlib.pyplot.close(fig)

    def test_wet_end_exceeds_dry_end(self, daily):
        curve = flow_duration_curve(daily)
        assert curve["flow_cfs"].iloc[0] > curve["flow_cfs"].iloc[-1]

    def test_arbitrary_percentiles(self, daily):
        curve = flow_duration_curve(daily, [2.5, 37.0])
        assert list(curve["exceedance_pct"]) == [2.5, 37.0]

    def test_missing_column_raises(self):
        with pytest.raises(KeyError, match="flow_cfs"):
            flow_duration_curve(pd.DataFrame({"q": [1.0, 2.0]}))

    def test_all_nan_raises(self):
        frame = pd.DataFrame({"flow_cfs": [np.nan, np.nan]})
        with pytest.raises(ValueError, match="no non-null flows"):
            flow_duration_curve(frame)

    def test_percentile_outside_range_raises(self, daily):
        with pytest.raises(ValueError, match="strictly between 0 and 100"):
            flow_duration_curve(daily, [0.0, 50.0])


class TestTransposeDuration:
    def test_transposes_every_point(self, daily, duration_exponents):
        curve = flow_duration_curve(daily)
        result = transpose_duration(curve, 100.0, 80.0, duration_exponents)
        assert len(result.quantiles) == len(curve)
        assert np.all(result.quantiles["flow_cfs"] < result.quantiles["donor_flow_cfs"])

    def test_arithmetic_matches_the_closed_form(self, daily, duration_exponents):
        curve = flow_duration_curve(daily)
        result = transpose_duration(curve, 100.0, 80.0, duration_exponents)
        expected = result.quantiles["donor_flow_cfs"] * 0.8 ** result.quantiles["exponent"]
        np.testing.assert_allclose(result.quantiles["flow_cfs"], expected, rtol=0, atol=0)

    def test_unit_ratio_is_the_identity(self, daily, duration_exponents):
        curve = flow_duration_curve(daily)
        result = transpose_duration(curve, 90.0, 90.0, duration_exponents)
        np.testing.assert_array_equal(
            result.quantiles["flow_cfs"].to_numpy(),
            result.quantiles["donor_flow_cfs"].to_numpy(),
        )

    def test_exponent_varies_across_the_curve(self, daily, duration_exponents):
        """The whole point of an exponent set rather than one number."""
        curve = flow_duration_curve(daily)
        result = transpose_duration(curve, 100.0, 80.0, duration_exponents)
        assert result.quantiles["exponent"].nunique() > 1
        wet = result.quantiles.iloc[0]["exponent"]
        dry = result.quantiles.iloc[-1]["exponent"]
        assert wet > dry

    def test_accepts_a_pct_indexed_curve(self, duration_exponents):
        curve = pd.DataFrame({"exceedance_pct": [10.0, 50.0], "flow_cfs": [100.0, 20.0]})
        result = transpose_duration(curve, 100.0, 80.0, duration_exponents)
        np.testing.assert_allclose(result.quantiles["exceedance_prob"], [0.10, 0.50])

    def test_zero_donor_statistic_is_nan_not_zero(self, duration_exponents):
        """0 * ratio is 0, which would assert the target is dry."""
        curve = pd.DataFrame(
            {"exceedance_prob": [0.10, 0.50, 0.99], "flow_cfs": [100.0, 20.0, 0.0]}
        )
        result = transpose_duration(curve, 100.0, 80.0, duration_exponents)
        assert np.isnan(result.quantiles["flow_cfs"].iloc[-1])
        assert result.provenance.degenerate_count == 1
        assert "could not be transposed" in result.to_markdown()

    def test_no_confidence_limits(self, daily, duration_exponents):
        curve = flow_duration_curve(daily)
        result = transpose_duration(curve, 100.0, 80.0, duration_exponents)
        assert result.confidence_limits.empty

    def test_flood_exponents_are_refused(self, daily):
        flood = RegressionExponents(
            aeps=[0.5, 0.01], exponents=[0.83, 0.76], citation="flood regression"
        )
        curve = flow_duration_curve(daily)
        with pytest.raises(ValueError, match="indexed by 'exceedance'"):
            transpose_duration(curve, 100.0, 80.0, flood)

    def test_empty_curve_raises(self, duration_exponents):
        with pytest.raises(ValueError, match="curve is empty"):
            transpose_duration(pd.DataFrame(), 100.0, 80.0, duration_exponents)

    def test_area_band_still_enforced(self, daily, duration_exponents):
        curve = flow_duration_curve(daily)
        with pytest.raises(ValueError, match="outside the applicable range"):
            transpose_duration(curve, 100.0, 40.0, duration_exponents)

    def test_markdown_uses_the_exceedance_label(self, daily, duration_exponents):
        curve = flow_duration_curve(daily)
        md = transpose_duration(curve, 100.0, 80.0, duration_exponents).to_markdown()
        assert "| Exceedance |" in md
        assert "flow duration" in md


class TestTransposeLowFlow:
    def test_transposes_and_scales_down(self, low_flow_exponents):
        donor = _low_flow_results([12.0, 8.0, 6.0, 4.5, 4.0])
        result = transpose_low_flow(
            donor,
            100.0,
            80.0,
            low_flow_exponents,
            hydrogeologic_setting="Valley and Ridge carbonate, same aquifer",
        )
        assert np.all(result.quantiles["flow_cfs"] < result.quantiles["donor_flow_cfs"])
        assert len(result.quantiles) == 5

    def test_unit_ratio_is_the_identity(self, low_flow_exponents):
        donor = _low_flow_results([12.0, 8.0, 6.0, 4.5, 4.0])
        result = transpose_low_flow(
            donor, 90.0, 90.0, low_flow_exponents, hydrogeologic_setting="same aquifer"
        )
        np.testing.assert_array_equal(
            result.quantiles["flow_cfs"].to_numpy(),
            result.quantiles["donor_flow_cfs"].to_numpy(),
        )

    def test_flood_exponents_are_refused_for_low_flow(self):
        """The category error the design doc names, made unreachable."""
        flood = RegressionExponents(
            aeps=[0.5, 0.01], exponents=[0.83, 0.76], citation="flood regression"
        )
        donor = _low_flow_results([12.0, 8.0, 6.0, 4.5, 4.0])
        with pytest.raises(ValueError, match="indexed by 'non_exceedance'"):
            transpose_low_flow(donor, 100.0, 80.0, flood, hydrogeologic_setting="same aquifer")

    def test_hydrogeologic_setting_is_required(self, low_flow_exponents):
        donor = _low_flow_results([12.0, 8.0, 6.0, 4.5, 4.0])
        with pytest.raises(ValueError, match="hydrogeologic_setting"):
            transpose_low_flow(donor, 100.0, 80.0, low_flow_exponents, hydrogeologic_setting="  ")

    def test_zero_donor_statistic_refuses(self, low_flow_exponents):
        donor = _low_flow_results([12.0, 8.0, 6.0, 4.5, 0.0])
        with pytest.raises(ValueError, match="zero or non-positive"):
            transpose_low_flow(
                donor, 100.0, 80.0, low_flow_exponents, hydrogeologic_setting="same aquifer"
            )

    def test_donor_that_goes_dry_refuses(self, low_flow_exponents):
        donor = _low_flow_results([12.0, 8.0, 6.0, 4.5, 4.0], p_zero=0.1)
        with pytest.raises(ValueError, match="goes dry"):
            transpose_low_flow(
                donor, 100.0, 80.0, low_flow_exponents, hydrogeologic_setting="same aquifer"
            )

    def test_p_zero_can_be_allowed_explicitly(self, low_flow_exponents):
        donor = _low_flow_results([12.0, 8.0, 6.0, 4.5, 4.0], p_zero=0.1)
        result = transpose_low_flow(
            donor,
            100.0,
            80.0,
            low_flow_exponents,
            hydrogeologic_setting="same aquifer",
            max_p_zero=0.25,
        )
        assert len(result.quantiles) == 5

    def test_tighter_area_band_is_the_default(self, low_flow_exponents):
        """A ratio of 0.6 is fine for floods and refused for low flow."""
        donor = _low_flow_results([12.0, 8.0, 6.0, 4.5, 4.0])
        assert LOW_FLOW_AREA_RATIO_RANGE == (0.7, 1.3)
        with pytest.raises(ValueError, match="outside the applicable range"):
            transpose_low_flow(
                donor, 100.0, 60.0, low_flow_exponents, hydrogeologic_setting="same aquifer"
            )

    def test_bfi_mismatch_is_recorded(self, low_flow_exponents):
        donor = _low_flow_results([12.0, 8.0, 6.0, 4.5, 4.0])
        result = transpose_low_flow(
            donor,
            100.0,
            80.0,
            low_flow_exponents,
            hydrogeologic_setting="same aquifer",
            donor_bfi=0.80,
            target_bfi=0.40,
        )
        assert result.provenance.bfi_difference == pytest.approx(0.40)
        assert "baseflow indices differ" in result.to_markdown()

    def test_setting_and_caveat_reach_the_report(self, low_flow_exponents):
        donor = _low_flow_results([12.0, 8.0, 6.0, 4.5, 4.0])
        md = transpose_low_flow(
            donor,
            100.0,
            80.0,
            low_flow_exponents,
            hydrogeologic_setting="Valley and Ridge carbonate",
        ).to_markdown()
        assert "Valley and Ridge carbonate" in md
        assert "| Non-exceedance |" in md
        assert "7-day low flow" in md
        assert "better tool than transposition" in md

    def test_duration_exponents_are_refused_for_low_flow(self, duration_exponents):
        donor = _low_flow_results([12.0, 8.0, 6.0, 4.5, 4.0])
        with pytest.raises(ValueError, match="indexed by 'non_exceedance'"):
            transpose_low_flow(
                donor, 100.0, 80.0, duration_exponents, hydrogeologic_setting="same aquifer"
            )


class TestKindGuardIsSymmetric:
    def test_transpose_frequency_refuses_duration_exponents(self, duration_exponents):
        from flowfreq.bulletin17c import Bulletin17C

        donor = Bulletin17C(
            peak_flows=np.array([100.0, 200.0, 150.0, 300.0, 250.0, 400.0]),
            water_years=np.arange(2000, 2006),
        ).run_analysis(method="mom")
        with pytest.raises(ValueError, match="indexed by 'aep'"):
            transpose_frequency(donor, 100.0, 80.0, duration_exponents)

    def test_unknown_kind_rejected_at_construction(self):
        with pytest.raises(ValueError, match="unknown probability_kind"):
            RegressionExponents(aeps=[0.5], exponents=[0.8], citation="x", probability_kind="vibes")


class TestMonotonicity:
    """A transposed curve that turns back on itself is physically impossible.

    This was a real defect in the first cut of ``transpose_duration``: with a
    groundwater-dominated donor (flat dry end) and an exponent set that falls
    toward the dry end -- the exact combination the PNW duration-regression
    literature describes -- the transposed curve inverted at area ratios of
    0.5 to 0.7, *inside* the supported band. It shipped because the flood path
    had a monotonicity test and the duration path did not.
    """

    #: Flat dry end: the signature of a basin sustained by storage.
    GROUNDWATER_CURVE = pd.DataFrame(
        {
            "exceedance_prob": [0.01, 0.10, 0.50, 0.90, 0.95, 0.99],
            "flow_cfs": [800.0, 260.0, 90.0, 52.0, 50.5, 50.0],
        }
    )

    @staticmethod
    def _steep_exponents() -> RegressionExponents:
        return RegressionExponents(
            aeps=[0.01, 0.10, 0.50, 0.90, 0.95, 0.99],
            exponents=[1.00, 0.95, 0.80, 0.62, 0.58, 0.55],
            citation="steep duration set, for the inversion case",
            probability_kind="exceedance",
        )

    @pytest.mark.parametrize("target_area", [50.0, 60.0, 70.0])
    def test_inversion_inside_the_supported_band_raises(self, target_area):
        with pytest.raises(ValueError, match="not monotone"):
            transpose_duration(self.GROUNDWATER_CURVE, 100.0, target_area, self._steep_exponents())

    def test_warn_mode_returns_and_records_it(self):
        result = transpose_duration(
            self.GROUNDWATER_CURVE,
            100.0,
            60.0,
            self._steep_exponents(),
            on_non_monotonic="warn",
        )
        assert result.provenance.monotonic is False

    def test_a_well_behaved_transposition_is_monotone(self):
        result = transpose_duration(self.GROUNDWATER_CURVE, 100.0, 100.0, self._steep_exponents())
        assert result.provenance.monotonic is True
        flows = result.quantiles["flow_cfs"].to_numpy()
        assert np.all(np.diff(flows) < 0)

    def test_unknown_mode_rejected(self, duration_exponents, daily):
        curve = flow_duration_curve(daily)
        with pytest.raises(ValueError, match="unknown on_non_monotonic"):
            transpose_duration(curve, 100.0, 80.0, duration_exponents, on_non_monotonic="shrug")

    def test_low_flow_direction_is_the_other_way(self, low_flow_exponents):
        """Low flows *rise* with non-exceedance probability -- 7Q2 exceeds
        7Q10 -- so the check must not simply reuse the falling rule."""
        donor = _low_flow_results([12.0, 8.0, 6.0, 4.5, 4.0])
        result = transpose_low_flow(
            donor, 100.0, 80.0, low_flow_exponents, hydrogeologic_setting="same aquifer"
        )
        assert result.provenance.monotonic is True

    def test_degenerate_rows_do_not_break_the_check(self, duration_exponents):
        """NaN rows from a zero donor statistic are skipped, not treated as
        an inversion."""
        curve = pd.DataFrame(
            {"exceedance_prob": [0.10, 0.50, 0.99], "flow_cfs": [100.0, 20.0, 0.0]}
        )
        result = transpose_duration(curve, 100.0, 80.0, duration_exponents)
        assert result.provenance.monotonic is True
        assert result.provenance.degenerate_count == 1
