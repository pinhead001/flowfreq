"""
hydrolib.workflow - High-level Bulletin 17C analysis entry points.

One call from annual peaks to a fitted frequency curve, plus the skew-variant
helpers that go with it. This is the layer a consumer wants when it does not
want to assemble :class:`~hydrolib.bulletin17c.Bulletin17C`, quantiles and
confidence limits by hand.

Everything here returns plain numbers and DataFrames of numbers. Turning those
into labelled, rounded, comma-separated strings is presentation and belongs to
the caller.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .bulletin17c import Bulletin17C
from .core import kfactor_array

logger = logging.getLogger(__name__)

#: Nationwide B17C default generalized skew (England et al. 2019).
B17C_DEFAULT_SKEW: float = -0.302

#: Return intervals reported by :func:`run_ffa` and :func:`compute_skew_tables`.
DEFAULT_RETURN_INTERVALS: List[float] = [1.5, 2, 5, 10, 25, 50, 100, 200, 500]

#: The same series as annual exceedance probabilities.
#: = [0.667, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
DEFAULT_AEP: List[float] = [1 / ri for ri in DEFAULT_RETURN_INTERVALS]

#: Canonical skew option labels, in the order a report presents them.
SKEW_OPTIONS: List[str] = ["Station Skew", "Weighted Skew", "Regional Skew"]


def _low_outlier_source(override: Optional[float]) -> str:
    """Describe where the reported PILF threshold came from.

    Parameters
    ----------
    override : float or None
        The user's requested threshold, if any.

    Returns
    -------
    str
        ``"MGBT"`` or ``"override"``. Both EMA and MOM censor on a supplied
        threshold now, so the label needs no method-specific caveat.
    """
    if override is None:
        return "MGBT"
    return "override"


def run_ffa(
    peak_flows: np.ndarray,
    water_years: np.ndarray,
    regional_skew: float = B17C_DEFAULT_SKEW,
    regional_skew_se: float = 0.55,
    perception_thresholds: Optional[List[dict]] = None,
    low_outlier_threshold_override: Optional[float] = None,
) -> dict:
    """Run Bulletin 17C flood frequency analysis.

    Fits by EMA, falling back to MOM only when EMA fails to converge *and* no
    perception thresholds are in play -- MOM cannot represent censored
    intervals, so a threshold-bearing record keeps its non-converged EMA fit
    rather than silently losing the thresholds.

    Errors are returned, not raised: any failure comes back as a string under
    ``error`` with the other keys left at their empty defaults. That suits an
    interactive caller that wants to show the message rather than crash. A
    caller that would rather have an exception should check ``error`` and
    raise its own.

    Parameters
    ----------
    peak_flows : np.ndarray
        Annual peak flows in cfs.
    water_years : np.ndarray
        Corresponding water years.
    regional_skew : float
        Regional skew coefficient. Defaults to the nationwide B17C value.
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
        result and censors all peaks below this value.  The threshold actually
        applied, its source and the resulting PILF count come back under
        ``parameters`` so a caller can show which cut produced the fit.

    Returns
    -------
    dict
        Keys: b17c, converged, method, parameters, quantile_df, error.

    Examples
    --------
    >>> result = run_ffa(peak_flows, water_years)
    >>> result["quantile_df"]["Flow (cfs)"]

    >>> result = run_ffa(peak_flows, water_years, regional_skew=-0.05,
    ...                  low_outlier_threshold_override=500.0)
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

        aep = np.array(DEFAULT_AEP)
        quantiles_df = b17c.compute_quantiles(aep=aep)
        ci_df = b17c.compute_confidence_limits(aep=aep)

        quantile_df = pd.DataFrame(
            {
                "Return Interval (yr)": DEFAULT_RETURN_INTERVALS,
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
                    # The low-outlier cut and where it came from. Without these
                    # a caller could offer the override but never show its effect.
                    # Both EMA and MOM censor on it now, so the source is just
                    # whether it came from MGBT or the user.
                    "low_outlier_threshold": r.low_outlier_threshold,
                    "n_low_outliers": r.n_low_outliers,
                    "low_outlier_source": _low_outlier_source(lo_override),
                },
                "quantile_df": quantile_df,
            }
        )

    except Exception as e:
        logger.exception("FFA analysis failed")
        result["error"] = str(e)

    return result


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
        Subset of :data:`SKEW_OPTIONS`.

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
    aep = np.array(DEFAULT_AEP)
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
                "Return Interval (yr)": DEFAULT_RETURN_INTERVALS,
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
        Skew labels the caller has selected.

    Returns
    -------
    dict[str, float]
        Empty dict (fall back to default) when no valid labels are found.
    """
    skew_map = _skew_values_from_result(ffa_result)
    return {lbl: skew_map[lbl] for lbl in selected_labels if skew_map.get(lbl) is not None}
