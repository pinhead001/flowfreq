"""Tests for app.ffa_export module."""

import io
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Ensure app package is importable
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.ffa_export import export_comparison_csv, export_ffa_to_zip


@pytest.fixture
def mock_ffa():
    return {
        "converged": True,
        "method": "ema",
        "parameters": {
            "mean_log": 3.5,
            "std_log": 0.3,
            "skew_station": -0.1,
            "skew_weighted": -0.2,
            "skew_used": -0.2,
            "regional_skew": -0.302,
        },
        "quantile_df": pd.DataFrame(
            {
                "Return Interval (yr)": [1.5, 2, 5, 10, 25, 50, 100, 200, 500],
                "AEP (%)": [66.7, 50.0, 20.0, 10.0, 4.0, 2.0, 1.0, 0.5, 0.2],
                "Flow (cfs)": [
                    1000,
                    2000,
                    5000,
                    8000,
                    12000,
                    15000,
                    20000,
                    25000,
                    35000,
                ],
                "Lower 90% CI": [
                    800,
                    1600,
                    4000,
                    6500,
                    9500,
                    12000,
                    16000,
                    19000,
                    27000,
                ],
                "Upper 90% CI": [
                    1200,
                    2500,
                    6500,
                    10000,
                    16000,
                    20000,
                    27000,
                    35000,
                    50000,
                ],
            }
        ),
        "error": None,
    }


def _make_zip():
    buf = io.BytesIO()
    return zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED), buf


class TestExportFFAToZip:
    def test_export_ffa_creates_expected_files(self, mock_ffa):
        zf, buf = _make_zip()
        with zf:
            export_ffa_to_zip(zf, "12345678", mock_ffa)
        with zipfile.ZipFile(buf) as zr:
            names = zr.namelist()
        assert "12345678/frequency_table.csv" in names
        assert "12345678/lp3_parameters.csv" in names

    def test_frequency_table_csv_has_9_rows(self, mock_ffa):
        zf, buf = _make_zip()
        with zf:
            export_ffa_to_zip(zf, "12345678", mock_ffa)
        with zipfile.ZipFile(buf) as zr:
            content = zr.read("12345678/frequency_table.csv").decode()
        df = pd.read_csv(io.StringIO(content))
        assert len(df) == 9

    def test_parameters_csv_has_expected_columns(self, mock_ffa):
        zf, buf = _make_zip()
        with zf:
            export_ffa_to_zip(zf, "12345678", mock_ffa)
        with zipfile.ZipFile(buf) as zr:
            content = zr.read("12345678/lp3_parameters.csv").decode()
        df = pd.read_csv(io.StringIO(content))
        for col in ["mean_log", "std_log", "method", "converged"]:
            assert col in df.columns


class TestExportComparisonCSV:
    def test_export_comparison_csv_format(self, mock_ffa):
        site_results = {
            "12345678": {
                "site_name": "Test Creek",
                "drainage_area_sqmi": 100.0,
                "ffa_results": mock_ffa,
            }
        }
        zf, buf = _make_zip()
        with zf:
            export_comparison_csv(zf, site_results)
        with zipfile.ZipFile(buf) as zr:
            content = zr.read("comparison_summary.csv").decode()
        df = pd.read_csv(io.StringIO(content))
        expected_cols = [
            "Site No",
            "Site Name",
            "Drainage Area (sq mi)",
            "100-yr Flow (cfs)",
            "100-yr Flow/DA (cfs/sq mi)",
            "Weighted Skew",
            "Method",
            "Converged",
        ]
        for col in expected_cols:
            assert col in df.columns
