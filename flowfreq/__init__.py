"""
flowfreq - Python library for hydrologic analysis

Includes:
- USGS gage data download (daily and peak flows)
- Summary hydrograph plotting
- Bulletin 17C flood frequency analysis:
  - Method of Moments (MOM)
  - Expected Moments Algorithm (EMA) with full PeakFQ parity
- Historical flood information handling
- Technical report generation
- One-call analysis workflow (flowfreq.workflow.run_ffa)
"""

import logging

from .bulletin17c import (
    Bulletin17C,
    ExpectedMomentsAlgorithm,
    FloodFrequencyAnalysis,
    MethodOfMoments,
)
from .core import (
    AnalysisMethod,
    EMAParameters,
    FlowInterval,
    FrequencyResults,
    LowFlowResults,
    SkewMethod,
    grubbs_beck_critical_value,
    kfactor,
    kfactor_array,
)
from .flowio import load_flow_frame, save_flow_frame
from .hydrograph import Hydrograph
from .lowflow import LOW_FLOW_YEAR_TYPES, LowFlowFrequency, annual_minimum_flow
from .regime import (
    BASEFLOW_METHODS,
    FlowRegime,
    baseflow_index,
    diel_variation,
    diel_variation_summary,
    monthly_flow_summary,
    richards_baker_flashiness,
    seasonal_flow_summary,
    separate_baseflow,
    tqmean,
)
from .report import HydroReport
from .usgs import (
    GageAttributes,
    NoInstantaneousDataError,
    USGSgage,
    fetch_nwis_batch,
    fetch_nwis_peaks,
)
from .workflow import (
    B17C_DEFAULT_SKEW,
    DEFAULT_AEP,
    DEFAULT_RETURN_INTERVALS,
    SKEW_OPTIONS,
    build_skew_curves_dict,
    compute_skew_tables,
    run_ffa,
)

# Alias for backwards compatibility
USGSGage = USGSgage

logger = logging.getLogger(__name__)


def analyze_gage(
    site_no: str,
    method: str = "ema",
    regional_skew: float = None,
    regional_skew_mse: float = None,
    historical_peaks: list = None,
    output_dir: str = "./output",
) -> dict:
    """
    Complete flood frequency analysis for a USGS gage.

    Parameters
    ----------
    site_no : str
        USGS site number
    method : str
        'mom' or 'ema' (default: 'ema')
    regional_skew : float, optional
        Regional skew coefficient
    regional_skew_mse : float, optional
        Mean squared error of regional skew
    historical_peaks : list of (year, flow) tuples, optional
        Historical peak observations
    output_dir : str
        Output directory
    """
    import os

    os.makedirs(output_dir, exist_ok=True)

    logger.info("Downloading data for USGS %s...", site_no)
    gage = USGSgage(site_no)

    try:
        gage.download_daily_flow()
        logger.info("Downloaded %d days of daily flow data", len(gage.daily_data))
    except Exception as e:
        logger.warning("Could not download daily flow data: %s", e)

    gage.download_peak_flow()
    logger.info("Downloaded %d annual peak flow records", len(gage.peak_data))
    logger.info("Site name: %s", gage.site_name)

    logger.info("Running Bulletin 17C analysis (method=%s)...", method.upper())

    water_years = gage.peak_data["water_year"].values

    analysis = Bulletin17C(
        gage.peak_data["peak_flow_cfs"].values,
        water_years=water_years,
        regional_skew=regional_skew,
        regional_skew_mse=regional_skew_mse,
        historical_peaks=historical_peaks,
    )

    results = analysis.run_analysis(method=method)

    logger.info("Station skew: %.4f", results.skew_station)
    if results.skew_weighted is not None:
        logger.info("Weighted skew: %.4f", results.skew_weighted)
    logger.info("Low outlier threshold: %s cfs", f"{results.low_outlier_threshold:,.0f}")

    if results.method == AnalysisMethod.EMA:
        logger.info("EMA iterations: %s", results.ema_iterations)
        logger.info("EMA converged: %s", results.ema_converged)

    logger.info("Generating report and figures...")
    report = HydroReport(gage, analysis)
    figures = report.generate_all_figures(output_dir)

    report_path = os.path.join(output_dir, "flood_frequency_report.md")
    report.save_report(report_path)
    logger.info("Report saved to: %s", report_path)

    return {
        "gage": gage,
        "analysis": analysis,
        "results": results,
        "figures": figures,
        "report_path": report_path,
    }


__version__ = "0.4.0"
__author__ = "FlowFreq"

__all__ = [
    # Core
    "AnalysisMethod",
    "SkewMethod",
    "FlowInterval",
    "EMAParameters",
    "FrequencyResults",
    "kfactor",
    "kfactor_array",
    "grubbs_beck_critical_value",
    # Low-flow frequency analysis
    "LowFlowResults",
    "LowFlowFrequency",
    "annual_minimum_flow",
    "LOW_FLOW_YEAR_TYPES",
    # Flow regime metrics
    "FlowRegime",
    "richards_baker_flashiness",
    "tqmean",
    "baseflow_index",
    "separate_baseflow",
    "monthly_flow_summary",
    "seasonal_flow_summary",
    "BASEFLOW_METHODS",
    "diel_variation",
    "diel_variation_summary",
    # USGS data retrieval
    "USGSgage",
    "USGSGage",  # Alias for backwards compatibility
    "GageAttributes",
    "NoInstantaneousDataError",
    "fetch_nwis_peaks",
    "fetch_nwis_batch",
    "save_flow_frame",
    "load_flow_frame",
    # Hydrograph
    "Hydrograph",
    # Bulletin 17C
    "Bulletin17C",
    "MethodOfMoments",
    "ExpectedMomentsAlgorithm",
    "FloodFrequencyAnalysis",
    # Report
    "HydroReport",
    # Convenience
    "analyze_gage",
]
