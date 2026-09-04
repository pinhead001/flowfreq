"""
Display formatting for the Streamlit FFA app.

The analysis lives in :mod:`flowfreq.workflow`; this module turns its numbers
into the labelled, rounded, comma-separated strings a table wants. Nothing here
computes anything -- that separation is what lets the library serve consumers
that have no use for these column headings.
"""

from __future__ import annotations

import logging
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)


def format_parameters_df(params: dict) -> pd.DataFrame:
    """Format analysis parameters as a single-row display DataFrame.

    Parameters
    ----------
    params : dict
        Parameters dict from run_ffa result.

    Returns
    -------
    pd.DataFrame
        Single-row DataFrame with formatted parameter values.
    """
    threshold = params.get("low_outlier_threshold", 0) or 0
    source = params.get("low_outlier_source", "MGBT")
    return pd.DataFrame(
        {
            "Mean (log10)": [f"{params.get('mean_log', 0):.4f}"],
            "Std Dev (log10)": [f"{params.get('std_log', 0):.4f}"],
            "Station Skew": [f"{params.get('skew_station', 0):.4f}"],
            "Weighted Skew": [f"{params.get('skew_weighted', 0):.4f}"],
            "Regional Skew": [f"{params.get('regional_skew', 0):.4f}"],
            "PILF Threshold (cfs)": [f"{threshold:,.0f} ({source})" if threshold > 0 else "none"],
            "PILFs": [f"{params.get('n_low_outliers', 0)}"],
        }
    )


def format_quantile_df(quantile_df: pd.DataFrame) -> pd.DataFrame:
    """Format quantile DataFrame for display.

    Parameters
    ----------
    quantile_df : pd.DataFrame
        Raw quantile DataFrame from run_ffa result.

    Returns
    -------
    pd.DataFrame
        Formatted DataFrame with comma-separated flows and percentage AEP.
    """
    df = quantile_df.copy()

    df["Return Interval (yr)"] = df["Return Interval (yr)"].apply(
        lambda x: "1.5" if x == 1.5 else f"{int(x)}"
    )

    df["AEP (%)"] = df["AEP (%)"].apply(lambda x: f"{x * 100:.1f}%")

    for col in ["Flow (cfs)", "Lower 90% CI", "Upper 90% CI"]:
        df[col] = df[col].apply(lambda x: f"{int(round(x)):,}")

    return df


def build_station_summary_df(
    site_no: str,
    peak_df: pd.DataFrame,
    ffa_result: dict,
    regional_skew: float,
    regional_skew_se: float,
    primary_skew_label: str = "Weighted Skew",
    latitude: Optional[float] = None,
    longitude: Optional[float] = None,
    map_skew_source: str = "B17C 2019 (Nationwide)",
) -> pd.DataFrame:
    """Build a PeakFQ-style station summary table for display.

    Parameters
    ----------
    site_no : str
        USGS site number.
    peak_df : pd.DataFrame
        Annual peak flow data (must have ``water_year`` column).
    ffa_result : dict
        Output from :func:`flowfreq.workflow.run_ffa`.
    regional_skew : float
        Regional skew input value.
    regional_skew_se : float
        Regional skew standard error.
    primary_skew_label : str
        The skew option currently selected (determines "Skew Option" field).
    latitude : float, optional
        Station latitude (decimal degrees).
    longitude : float, optional
        Station longitude (decimal degrees, negative = West).

    Returns
    -------
    pd.DataFrame
        Single-row DataFrame styled after PeakFQ station summary output.
    """
    b17c = ffa_result.get("b17c")
    r = b17c.results if b17c is not None else None

    start_year = int(peak_df["water_year"].min()) if not peak_df.empty else "N/A"
    end_year = int(peak_df["water_year"].max()) if not peak_df.empty else "N/A"
    n_sys = (r.n_systematic if r is not None else None) or len(peak_df)

    skew_option_map = {
        "Station Skew": "Station",
        "Weighted Skew": "Weighted",
        "Regional Skew": "Regional",
    }
    skew_option = skew_option_map.get(primary_skew_label, "Weighted")
    use_map_skew = "No" if primary_skew_label == "Station Skew" else "Yes"

    pilf_threshold = (
        f"{r.low_outlier_threshold:,.0f}" if r is not None and r.low_outlier_threshold > 0 else "0"
    )

    lat_str = f"{latitude:.5f}°N" if latitude is not None else "N/A"
    lon_str = f"{abs(longitude):.5f}°W" if longitude is not None else "N/A"
    mse = round(regional_skew_se**2, 4)

    return pd.DataFrame(
        {
            "Station ID": [site_no],
            "Start Year": [start_year],
            "End Year": [end_year],
            "Record Length": [n_sys],
            "Skew Option": [skew_option],
            "Use Map Skew": [use_map_skew],
            "Map Skew Source": [map_skew_source],
            "Regional Skew": [f"{regional_skew:.3f}"],
            "Reg Skew Std Err": [f"{regional_skew_se:.3f}"],
            "Mean Sqr Err": [f"{mse:.4f}"],
            "PILF (LO) Test": ["MGBT"],
            "PILF (LO) Threshold": [pilf_threshold],
            "Urban/Reg Peaks": ["No"],
            "Latitude": [lat_str],
            "Longitude": [lon_str],
        }
    )
