"""Tests for app.ffa_runner module."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Ensure app package is importable
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.ffa_runner import format_parameters_df, format_quantile_df, run_ffa
from tests.peakfqsa.fixtures.big_sandy import (
    REGIONAL_SKEW,
    REGIONAL_SKEW_SD,
    SYSTEMATIC_PEAKS,
)


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
        ]
        assert list(df.columns) == expected_cols
        assert len(df) == 1


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
