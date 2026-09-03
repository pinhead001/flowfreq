"""Tests for app.ffa_runner -- display formatting only.

The analysis these tables render moved to :mod:`hydrolib.workflow`; its tests
live in ``tests/test_workflow.py``. Nothing here calls into the library, which
is the point: this module must stay renderable without an analysis run.
"""

import sys
from pathlib import Path

import pandas as pd

# Ensure app package is importable
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.ffa_runner import format_parameters_df, format_quantile_df


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
