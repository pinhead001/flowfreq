"""Tests for app.ffa_runner module."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Ensure app package is importable
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.ffa_runner import _low_outlier_source, format_parameters_df, format_quantile_df, run_ffa
from tests.fixtures.big_sandy import REGIONAL_SKEW, REGIONAL_SKEW_SD, SYSTEMATIC_PEAKS


def _big_sandy_arrays():
    """Extract peak_flows and water_years arrays from big_sandy fixture."""
    years = np.array(sorted(SYSTEMATIC_PEAKS.keys()))
    flows = np.array([SYSTEMATIC_PEAKS[y] for y in years])
    return flows, years


class TestRunFFA:
    def test_run_ffa_returns_expected_keys(self):
        flows, years = _big_sandy_arrays()
        result = run_ffa(flows, years)
        expected_keys = {"b17c", "converged", "method", "parameters", "quantile_df", "error"}
        assert set(result.keys()) == expected_keys

    def test_run_ffa_with_big_sandy_data(self):
        flows, years = _big_sandy_arrays()
        result = run_ffa(
            flows, years, regional_skew=REGIONAL_SKEW, regional_skew_se=REGIONAL_SKEW_SD
        )
        assert result["error"] is None
        assert result["converged"] is True
        assert result["method"] in ("ema", "mom")
        assert len(result["quantile_df"]) == 9
        assert result["b17c"] is not None

    def test_run_ffa_handles_bad_data(self):
        result = run_ffa(np.array([]), np.array([]))
        # Should return an error for empty input
        assert result["error"] is not None


class TestFormatParametersDf:
    def test_format_parameters_df_columns(self):
        params = {
            "mean_log": 3.7173,
            "std_log": 0.2892,
            "skew_station": -0.1187,
            "skew_weighted": -0.1187,
            "regional_skew": -0.302,
        }
        df = format_parameters_df(params)
        expected_cols = [
            "Mean (log10)",
            "Std Dev (log10)",
            "Station Skew",
            "Weighted Skew",
            "Regional Skew",
            "PILF Threshold (cfs)",
            "PILFs",
        ]
        assert list(df.columns) == expected_cols
        assert len(df) == 1

    def test_missing_pilf_keys_render_as_none(self):
        """Older result dicts have no PILF keys; the table must not blow up."""
        df = format_parameters_df({"mean_log": 3.7})
        assert df["PILF Threshold (cfs)"].iloc[0] == "none"
        assert df["PILFs"].iloc[0] == "0"

    def test_pilf_threshold_names_its_source(self):
        """A user needs to know whether the cut is MGBT's or their own."""
        mgbt = format_parameters_df(
            {"low_outlier_threshold": 1234.0, "n_low_outliers": 3, "low_outlier_source": "MGBT"}
        )
        assert mgbt["PILF Threshold (cfs)"].iloc[0] == "1,234 (MGBT)"
        assert mgbt["PILFs"].iloc[0] == "3"

        override = format_parameters_df(
            {"low_outlier_threshold": 500.0, "n_low_outliers": 1, "low_outlier_source": "override"}
        )
        assert override["PILF Threshold (cfs)"].iloc[0] == "500 (override)"


class TestPilfOverride:
    """The override has to change the fit, or the control is decorative."""

    OVERRIDE = 2000.0  # censors four Big Sandy peaks; EMA still converges

    def test_override_censors_more_peaks_than_mgbt(self):
        flows, years = _big_sandy_arrays()
        base = run_ffa(peak_flows=flows, water_years=years)
        forced = run_ffa(
            peak_flows=flows, water_years=years, low_outlier_threshold_override=self.OVERRIDE
        )
        assert forced["parameters"]["n_low_outliers"] > base["parameters"]["n_low_outliers"]
        assert forced["parameters"]["low_outlier_threshold"] == pytest.approx(self.OVERRIDE)
        assert forced["parameters"]["low_outlier_source"] == "override"
        assert base["parameters"]["low_outlier_source"] == "MGBT"

    def test_override_moves_the_fitted_moments(self):
        flows, years = _big_sandy_arrays()
        base = run_ffa(peak_flows=flows, water_years=years)
        forced = run_ffa(
            peak_flows=flows, water_years=years, low_outlier_threshold_override=self.OVERRIDE
        )
        assert base["parameters"]["mean_log"] != forced["parameters"]["mean_log"]

    def test_zero_and_none_both_mean_use_mgbt(self):
        flows, years = _big_sandy_arrays()
        zero = run_ffa(peak_flows=flows, water_years=years, low_outlier_threshold_override=0.0)
        none = run_ffa(peak_flows=flows, water_years=years, low_outlier_threshold_override=None)
        assert zero["parameters"]["low_outlier_source"] == "MGBT"
        assert zero["parameters"]["mean_log"] == none["parameters"]["mean_log"]

    def test_mom_fallback_reports_the_override_but_says_it_was_not_acted_on(self):
        """A high override stops EMA converging, and MOM does not censor at all.

        MethodOfMoments now reports the threshold the user asked for rather
        than substituting a Grubbs-Beck value they did not -- but it computes
        its moments from every peak regardless. The label has to say both, or
        it contradicts the number displayed beside it.
        """
        flows, years = _big_sandy_arrays()
        result = run_ffa(peak_flows=flows, water_years=years, low_outlier_threshold_override=6000.0)
        assert result["method"] == "mom"
        assert result["parameters"]["low_outlier_threshold"] == pytest.approx(6000.0)
        source = result["parameters"]["low_outlier_source"]
        assert "override" in source and "does not censor" in source

    def test_mom_fallback_moments_are_untouched_by_the_override(self):
        """The claim the label makes must actually hold."""
        flows, years = _big_sandy_arrays()
        forced = run_ffa(peak_flows=flows, water_years=years, low_outlier_threshold_override=6000.0)
        assert forced["method"] == "mom"
        plain = run_ffa(peak_flows=flows, water_years=years)
        assert forced["parameters"]["mean_log"] == pytest.approx(plain["parameters"]["mean_log"])

    def test_source_helper_covers_every_combination(self):
        assert _low_outlier_source(None, "ema") == "MGBT"
        assert _low_outlier_source(None, "mom") == "MGBT"
        assert _low_outlier_source(1000.0, "ema") == "override"
        assert "does not censor" in _low_outlier_source(1000.0, "mom")


class TestFormatQuantileDf:
    def test_format_quantile_df_formatting(self):
        raw = pd.DataFrame(
            {
                "Return Interval (yr)": [1.5, 2, 100],
                "AEP (%)": [0.667, 0.50, 0.01],
                "Flow (cfs)": [3957.5, 5284.4, 23158.7],
                "Lower 90% CI": [3000.0, 4000.0, 17000.0],
                "Upper 90% CI": [5000.0, 7000.0, 38000.0],
            }
        )
        df = format_quantile_df(raw)

        # Return intervals formatted correctly
        assert df["Return Interval (yr)"].iloc[0] == "1.5"
        assert df["Return Interval (yr)"].iloc[1] == "2"
        assert df["Return Interval (yr)"].iloc[2] == "100"

        # AEP as percentage
        assert df["AEP (%)"].iloc[2] == "1.0%"

        # Flow with commas
        assert df["Flow (cfs)"].iloc[2] == "23,159"
