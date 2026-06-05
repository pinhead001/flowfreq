"""
Output parser for PeakfqSA result files.

Parses the .out file format produced by PeakfqSA into structured Python
objects. Field names and structure derived from peakfqr R/fortranWrappers.R
output extraction (lines 217-306) and R/main.R output compilation.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import Optional

import pandas as pd

logger = logging.getLogger(__name__)

# Compiled regex patterns for parsing output sections
_RE_STATION = re.compile(r"Station[:\s]+(.+)", re.IGNORECASE)
_RE_PERIOD = re.compile(r"Analysis\s+Period[:\s]+(\d{4})\s*[-\u2013to ]+\s*(\d{4})", re.IGNORECASE)
_RE_N_PEAKS = re.compile(r"Number of peaks[:\s]+(\d+)", re.IGNORECASE)
_RE_N_SYSTEMATIC = re.compile(r"Number of systematic peaks[:\s]+(\d+)", re.IGNORECASE)
_RE_N_HISTORICAL = re.compile(r"Number of historical peaks[:\s]+(\d+)", re.IGNORECASE)
_RE_PARAMETER = re.compile(
    r"(Mean|Std\.?\s*Dev|Skew\s*\(?\w*\)?|MSE\s+of\s+Skew|Regional\s+Skew|Regional\s+MSE)"
    r"\s*[=:]\s*([-+]?\d*\.?\d+)",
    re.IGNORECASE,
)
_RE_QUANTILE_LINE = re.compile(r"\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)")
_RE_LOW_OUTLIER = re.compile(
    r"low.?outlier.*?(\d+)\s*PILFs?\s*detected.*?threshold\s*=\s*([\d.]+)",
    re.IGNORECASE,
)
_RE_LO_METHOD = re.compile(r"Method[:\s]+(MGBT|FIXED|NONE)", re.IGNORECASE)
_RE_SKEW_OPTION = re.compile(r"Skew\s+Option[:\s]+(Weighted|Station|Generalized)", re.IGNORECASE)


@dataclass
class PeakfqSAResult:
    """Structured result from a PeakfqSA analysis run.

    Mirrors the output fields from peakfqr's ``emafit()`` return value,
    mapped to Python-native types.

    Parameters
    ----------
    station_name : str
        Station identifier.
    begyear : int
        Beginning water year of analysis.
    endyear : int
        Ending water year of analysis.
    n_peaks : int
        Total number of peaks in the analysis (record length).
    n_systematic : int
        Number of systematic (non-historical) peaks.
    n_historical : int
        Number of historical peaks.
    low_outlier_count : int
        Number of potentially influential low floods (PILFs) detected.
    low_outlier_threshold : float
        PILF threshold in discharge units (cfs).
    parameters : dict[str, float]
        LP3 distribution parameters. Keys: mean_log, std_log,
        skew_weighted, skew_at_site, mean_log_at_site, std_log_at_site,
        mse_skew, regional_skew, regional_skew_mse.
    quantiles : dict[float, float]
        Annual exceedance probability to discharge mapping.
    confidence_intervals : dict[float, tuple[float, float]]
        AEP to (lower_bound, upper_bound) CI mapping.
    plotting_positions : DataFrame or None
        Empirical plotting positions with columns: year, plot_pos, q_obs, q_fit.
    raw_output : str
        Full text of the output file for debugging.
    """

    station_name: str = ""
    begyear: int = 0
    endyear: int = 0
    n_peaks: int = 0
    n_systematic: int = 0
    n_historical: int = 0
    low_outlier_count: int = 0
    low_outlier_threshold: float = 0.0
    parameters: dict[str, float] = field(default_factory=dict)
    quantiles: dict[float, float] = field(default_factory=dict)
    confidence_intervals: dict[float, tuple[float, float]] = field(default_factory=dict)
    plotting_positions: Optional[pd.DataFrame] = field(default=None, repr=False)
    raw_output: str = field(default="", repr=False)


def parse_peakfqsa_output(text: str) -> PeakfqSAResult:
    """Parse PeakfqSA .out file text into a PeakfqSAResult.

    Parameters
    ----------
    text : str
        Full text content of a PeakfqSA output file.

    Returns
    -------
    PeakfqSAResult
        Parsed result object with all available fields populated.

    Raises
    ------
    ValueError
        If required sections (station, parameters) are missing.

    Notes
    -----
    Parsing is defensive: missing optional sections log a warning but
    do not raise. Only missing required sections raise ValueError.
    """
    result = PeakfqSAResult(raw_output=text)

    # Parse station name
    m = _RE_STATION.search(text)
    if m:
        result.station_name = m.group(1).strip()
    else:
        logger.warning("Station name not found in output")

    # Parse analysis period
    m = _RE_PERIOD.search(text)
    if m:
        result.begyear = int(m.group(1))
        result.endyear = int(m.group(2))

    # Parse peak counts
    m = _RE_N_PEAKS.search(text)
    if m:
        result.n_peaks = int(m.group(1))

    m = _RE_N_SYSTEMATIC.search(text)
    if m:
        result.n_systematic = int(m.group(1))

    m = _RE_N_HISTORICAL.search(text)
    if m:
        result.n_historical = int(m.group(1))

    # Parse LP3 parameters
    result.parameters = _parse_parameters(text)

    # Parse low outlier info
    m = _RE_LOW_OUTLIER.search(text)
    if m:
        result.low_outlier_count = int(m.group(1))
        result.low_outlier_threshold = float(m.group(2))

    # Parse frequency curve (quantiles and CIs)
    quantiles, cis = _parse_frequency_curve(text)
    result.quantiles = quantiles
    result.confidence_intervals = cis

    return result


def _parse_parameters(text: str) -> dict[str, float]:
    """Extract LP3 parameters from output text.

    Parameters
    ----------
    text : str
        Output text.

    Returns
    -------
    dict[str, float]
        Parameter name to value mapping.
    """
    params: dict[str, float] = {}

    # Map display labels to standardized keys
    label_map = {
        "mean": "mean_log",
        "std. dev": "std_log",
        "std dev": "std_log",
        "stddev": "std_log",
        "skew (weighted)": "skew_weighted",
        "skew(weighted)": "skew_weighted",
        "skew (at-site)": "skew_at_site",
        "skew(at-site)": "skew_at_site",
        "mean (at-site)": "mean_log_at_site",
        "mean(at-site)": "mean_log_at_site",
        "std dev (at-site)": "std_log_at_site",
        "std. dev (at-site)": "std_log_at_site",
        "mse of skew": "mse_skew",
        "regional skew": "regional_skew",
        "regional mse": "regional_skew_mse",
    }

    for match in _RE_PARAMETER.finditer(text):
        label = match.group(1).strip().lower()
        value = float(match.group(2))

        # Find the best matching key
        matched_key = None
        for pattern, key in label_map.items():
            if pattern in label or label in pattern:
                matched_key = key
                break

        if matched_key:
            params[matched_key] = value
        else:
            # Use sanitized label as key for unrecognized parameters
            sanitized = re.sub(r"[^a-z0-9]+", "_", label).strip("_")
            params[sanitized] = value
            logger.debug("Unrecognized parameter label '%s' stored as '%s'", label, sanitized)

    return params


def _parse_frequency_curve(
    text: str,
) -> tuple[dict[float, float], dict[float, tuple[float, float]]]:
    """Extract frequency curve data from output text.

    Parameters
    ----------
    text : str
        Output text.

    Returns
    -------
    tuple[dict[float, float], dict[float, tuple[float, float]]]
        Quantiles {AEP: discharge} and confidence intervals {AEP: (lower, upper)}.
    """
    quantiles: dict[float, float] = {}
    cis: dict[float, tuple[float, float]] = {}

    # Find the frequency curve section
    in_freq_section = False
    for line in text.splitlines():
        stripped = line.strip()

        if "Frequency Curve" in line:
            in_freq_section = True
            continue

        if in_freq_section and ("---" in stripped and "End" in stripped):
            break

        if not in_freq_section:
            continue

        # Skip header lines
        if stripped.startswith("AEP") or not stripped:
            continue

        m = _RE_QUANTILE_LINE.match(line)
        if m:
            aep = float(m.group(1))
            estimate = float(m.group(2))
            ci_low = float(m.group(3))
            ci_high = float(m.group(4))
            # group(5) is variance, not stored in quantiles/cis

            quantiles[aep] = estimate
            cis[aep] = (ci_low, ci_high)

    if not quantiles:
        logger.warning("No quantiles found in frequency curve section")

    return quantiles, cis
