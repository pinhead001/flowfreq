"""
FFA analysis runner module for Streamlit app.

Wraps hydrolib Bulletin17C analysis with display-friendly formatting.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from hydrolib.bulletin17c import Bulletin17C
from hydrolib.core import kfactor_array

logger = logging.getLogger(__name__)

DISPLAY_RETURN_INTERVALS = [1.5, 2, 5, 10, 25, 50, 100, 200, 500]
DISPLAY_AEP = [1 / ri for ri in DISPLAY_RETURN_INTERVALS]
# = [0.667, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]


def run_ffa(
    peak_flows: np.ndarray,
    water_years: np.ndarray,
    regional_skew: float = -0.302,
    regional_skew_se: float = 0.55,
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
        b17c = Bulletin17C(
            peak_flows=peak_flows,
            water_years=water_years,
            regional_skew=regional_skew,
            regional_skew_mse=regional_skew_se**2,
        )

        b17c.run_analysis(method="ema")
        method = "ema"
        converged = bool(b17c.results.ema_converged)

        if not converged:
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
