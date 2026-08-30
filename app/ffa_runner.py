"""
FFA analysis runner module for Streamlit app.

Wraps hydrolib Bulletin17C analysis with display-friendly formatting.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from hydrolib.bulletin17c import Bulletin17C
from hydrolib.core import kfactor_array

logger = logging.getLogger(__name__)

#: Nationwide B17C default generalized skew (England et al. 2019).
_B17C_DEFAULT_SKEW = -0.302

DISPLAY_RETURN_INTERVALS = [1.5, 2, 5, 10, 25, 50, 100, 200, 500]
DISPLAY_AEP = [1 / ri for ri in DISPLAY_RETURN_INTERVALS]
# = [0.667, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]


def run_ffa(
    peak_flows: np.ndarray,
    water_years: np.ndarray,
    regional_skew: float = -0.302,
    regional_skew_se: float = 0.55,
    perception_thresholds: Optional[List[dict]] = None,
    low_outlier_threshold_override: Optional[float] = None,
) -> dict:
    """Run Bulletin 17C flood frequency analysis.

    Parameters
    ----------
    peak_flows : np.ndarray
        Annual peak flows in cfs.
    water_years : np.ndarray
        Corresponding water years.
    regional_skew : float
        Regional skew coefficient.
    regional_skew_se : float
        Regional skew standard error.
    perception_thresholds : list of dict, optional
        Each dict has keys ``start_year``, ``end_year``, ``threshold_cfs``
        (legacy) or ``lower_cfs`` / ``upper_cfs``.  Converts to the
        ``Dict[Tuple[int,int], float]`` format expected by
        :class:`~hydrolib.bulletin17c.Bulletin17C` and passed to EMA so that
        years in each period without a recorded peak are treated as
        left-censored observations (peak < threshold).
    low_outlier_threshold_override : float, optional
        User-supplied PILF threshold (cfs).  When > 0, overrides the MGBT
        result and censors all peaks below this value.

    Returns
    -------
    dict
        Keys: b17c, converged, method, parameters, quantile_df, error.
    """
    result = {
        "b17c": None,
        "converged": False,
        "method": None,
        "parameters": {},
        "quantile_df": pd.DataFrame(),
        "error": None,
    }

    try:
        # Convert list of threshold dicts → Dict[Tuple[int,int], float]
        pt_dict: Optional[Dict[Tuple[int, int], float]] = None
        if perception_thresholds:
            pt_dict = {
                (int(t["start_year"]), int(t["end_year"])): float(t["threshold_cfs"])
                for t in perception_thresholds
                if float(t.get("threshold_cfs", 0)) > 0
            } or None

        lo_override = (
            float(low_outlier_threshold_override)
            if low_outlier_threshold_override and low_outlier_threshold_override > 0
            else None
        )

        b17c = Bulletin17C(
            peak_flows=peak_flows,
            water_years=water_years,
            regional_skew=regional_skew,
            regional_skew_mse=regional_skew_se**2,
            perception_thresholds=pt_dict,
            user_low_outlier_threshold=lo_override,
        )

        b17c.run_analysis(method="ema")
        method = "ema"
        converged = bool(b17c.results.ema_converged)

        # Only fall back to MOM when no perception thresholds are in play — MOM has no
        # mechanism to incorporate censored intervals, so we keep the (non-converged)
        # EMA result when thresholds extend the record.
        if not converged and not pt_dict:
            logger.warning("EMA did not converge, falling back to MOM")
            b17c.run_analysis(method="mom")
            method = "mom"
            converged = True

        aep = np.array(DISPLAY_AEP)
        quantiles_df = b17c.compute_quantiles(aep=aep)
        ci_df = b17c.compute_confidence_limits(aep=aep)

        quantile_df = pd.DataFrame(
            {
                "Return Interval (yr)": DISPLAY_RETURN_INTERVALS,
                "AEP (%)": aep,
                "Flow (cfs)": quantiles_df["flow_cfs"].values,
                "Lower 90% CI": ci_df["lower_5pct"].values,
                "Upper 90% CI": ci_df["upper_5pct"].values,
            }
        )

        r = b17c.results
        result.update(
            {
                "b17c": b17c,
                "converged": converged,
                "method": method,
                "parameters": {
                    "mean_log": r.mean_log,
                    "std_log": r.std_log,
                    "skew_station": r.skew_station,
                    "skew_weighted": r.skew_weighted,
                    "skew_used": r.skew_used,
                    "regional_skew": regional_skew,
                },
                "quantile_df": quantile_df,
            }
        )

    except Exception as e:
        logger.exception("FFA analysis failed")
        result["error"] = str(e)

    return result


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
    return pd.DataFrame(
        {
            "Mean (log10)": [f"{params.get('mean_log', 0):.4f}"],
            "Std Dev (log10)": [f"{params.get('std_log', 0):.4f}"],
            "Station Skew": [f"{params.get('skew_station', 0):.4f}"],
            "Weighted Skew": [f"{params.get('skew_weighted', 0):.4f}"],
            "Regional Skew": [f"{params.get('regional_skew', 0):.4f}"],
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


# ---------------------------------------------------------------------------
# Skew-variant helpers
# ---------------------------------------------------------------------------

#: Canonical skew option labels in display order.
SKEW_OPTIONS: List[str] = ["Station Skew", "Weighted Skew", "Regional Skew"]

#: Re-export so streamlit_app can import it without reaching into internals.
B17C_DEFAULT_SKEW: float = _B17C_DEFAULT_SKEW


def _skew_values_from_result(ffa_result: dict) -> Dict[str, Optional[float]]:
    """Return the three skew values stored in an ffa_result dict."""
    p = ffa_result.get("parameters", {})
    return {
        "Station Skew": p.get("skew_station"),
        "Weighted Skew": p.get("skew_weighted"),
        "Regional Skew": p.get("regional_skew"),
    }


def compute_skew_tables(
    ffa_result: dict,
    selected_labels: List[str],
) -> Dict[str, pd.DataFrame]:
    """Compute a raw quantile+CI table for each selected skew option.

    Uses the LP3 moments (mean_log, std_log) already fitted by EMA/MOM and
    substitutes the requested skew value to produce separate frequency tables
    without re-running the full analysis.

    Parameters
    ----------
    ffa_result : dict
        Output from :func:`run_ffa`.
    selected_labels : list[str]
        Subset of ``["Station Skew", "Weighted Skew", "Regional Skew"]``.

    Returns
    -------
    dict[str, pd.DataFrame]
        Maps label → DataFrame with columns:
        ``Return Interval (yr)``, ``AEP (%)``,
        ``Flow (cfs)``, ``Lower 90% CI``, ``Upper 90% CI``.
        Returns an empty dict if ffa_result has an error.
    """
    if ffa_result.get("error") or ffa_result.get("b17c") is None:
        return {}

    r = ffa_result["b17c"].results
    mean_log = r.mean_log
    std_log = r.std_log
    n = r.n_systematic or r.n_peaks

    skew_map = _skew_values_from_result(ffa_result)
    aep = np.array(DISPLAY_AEP)
    z_alpha = 1.6449  # norm.ppf(0.95): two-sided 90% CI

    tables: Dict[str, pd.DataFrame] = {}
    for label in selected_labels:
        skew_val = skew_map.get(label)
        if skew_val is None:
            continue

        K = kfactor_array(skew_val, aep)
        log_Q = mean_log + K * std_log
        Q = 10.0**log_Q

        var_factor = 1 / n + K**2 * (1 + 0.75 * skew_val**2) / (2 * (n - 1))
        se_log = std_log * np.sqrt(var_factor)
        lower = 10.0 ** (log_Q - z_alpha * se_log)
        upper = 10.0 ** (log_Q + z_alpha * se_log)

        tables[label] = pd.DataFrame(
            {
                "Return Interval (yr)": DISPLAY_RETURN_INTERVALS,
                "AEP (%)": aep,
                "Flow (cfs)": Q,
                "Lower 90% CI": lower,
                "Upper 90% CI": upper,
            }
        )

    return tables


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
        Output from :func:`run_ffa`.
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


def build_skew_curves_dict(
    ffa_result: dict,
    selected_labels: List[str],
) -> Dict[str, float]:
    """Return ``{label: skew_value}`` for the selected skew options.

    Intended for passing directly to
    :func:`hydrolib.freq_plot.plot_frequency_curve_streamlit` as the
    ``skew_curves`` argument.

    Parameters
    ----------
    ffa_result : dict
        Output from :func:`run_ffa`.
    selected_labels : list[str]
        Skew labels the user has checked.

    Returns
    -------
    dict[str, float]
        Empty dict (fall back to default) when no valid labels are found.
    """
    skew_map = _skew_values_from_result(ffa_result)
    return {lbl: skew_map[lbl] for lbl in selected_labels if skew_map.get(lbl) is not None}
