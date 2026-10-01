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
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _installed_version

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
from .qppq import (
    MONTHLY_SEASONS,
    SNOWMELT_SEASONS,
    FlowDurationCurve,
    QppqResult,
    apply_lag,
    center_of_timing,
    estimate_donor_lag,
    loocv_qppq,
    performance,
    qppq,
    rank_donors,
    seasonal_curves,
)
from .regime import (
    BASEFLOW_METHODS,
    DEFAULT_EXCEEDANCE_PCT,
    FlowRegime,
    baseflow_index,
    flow_duration_curve,
    monthly_flow_summary,
    richards_baker_flashiness,
    seasonal_flow_summary,
    separate_baseflow,
    tqmean,
)
from .report import HydroReport
from .subdaily import (
    DEFAULT_MAX_GAP_MULTIPLE,
    INTERVAL_LABELS,
    KNOWN_VALUE_UNITS,
    TRUE_ZERO_UNITS,
    circular_hour_statistics,
    daily_extreme_timing,
    diel_variation,
    diel_variation_summary,
    extreme_timing_summary,
    ramping_rate_summary,
    ramping_rates,
)
from .transpose import (
    DEFAULT_AREA_RATIO_RANGE,
    LOW_FLOW_AREA_RATIO_RANGE,
    PROBABILITY_KINDS,
    RegressionExponents,
    TransposedResults,
    TranspositionProvenance,
    transpose_duration,
    transpose_frequency,
    transpose_low_flow,
)
from .usgs import (
    GageAttributes,
    NoInstantaneousDataError,
    USGSgage,
    fetch_nwis_batch,
    fetch_nwis_peaks,
    join_flow_and_stage,
)
from .workflow import (
    B17C_DEFAULT_SKEW,
    B17C_DEFAULT_SKEW_SE,
    DEFAULT_AEP,
    DEFAULT_RETURN_INTERVALS,
    SKEW_OPTIONS,
    build_skew_curves_dict,
    compute_skew_tables,
    peak_code_kwargs,
    resolve_regional_skew,
    run_ffa,
)

# Alias for backwards compatibility
USGSGage = USGSgage

logger = logging.getLogger(__name__)


def _gage_analysis(
    peak_data,
    site_no: str,
    *,
    method: str,
    regional_skew: float = None,
    regional_skew_mse: float = None,
    historical_peaks: list = None,
    apply_peak_codes: bool = True,
) -> Bulletin17C:
    """Build ``analyze_gage``'s ``Bulletin17C`` from a downloaded peak frame.

    Codes go through :func:`flowfreq.workflow.peak_code_kwargs`, the same
    route :func:`~flowfreq.workflow.run_ffa` and
    :func:`~flowfreq.workflow.compare_engines` take.
    """
    codes = (
        peak_data["qualification_code"].tolist()
        if "qualification_code" in peak_data.columns
        else None
    )
    if apply_peak_codes and codes is not None and method.lower() != "ema":
        from .peak_codes import count_acted_on_codes

        acted_on = count_acted_on_codes(codes)
        if acted_on:
            raise ValueError(
                f"USGS {site_no}: peaks carry qualification codes peakfq acts on "
                f"({', '.join(f'code {c}: {n}' for c, n in acted_on.items())}), which need "
                "method='ema'; pass apply_peak_codes=False to fit them as exact values "
                "with the method of moments."
            )
    kwargs = peak_code_kwargs(
        peak_data["peak_flow_cfs"].values,
        peak_data["water_year"].values,
        codes,
        apply_peak_codes=apply_peak_codes,
        regional_skew=regional_skew,
        regional_skew_mse=regional_skew_mse,
        historical_peaks=historical_peaks,
        site_name=f"USGS {site_no}",
    )
    if kwargs is not None:
        return Bulletin17C(**kwargs)
    return Bulletin17C(
        peak_data["peak_flow_cfs"].values,
        water_years=peak_data["water_year"].values,
        regional_skew=regional_skew,
        regional_skew_mse=regional_skew_mse,
        historical_peaks=historical_peaks,
    )


def analyze_gage(
    site_no: str,
    method: str = "ema",
    regional_skew: float = None,
    regional_skew_mse: float = None,
    historical_peaks: list = None,
    output_dir: str = "./output",
    apply_peak_codes: bool = True,
    allow_regulated: bool = False,
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
    apply_peak_codes : bool, default True
        Apply the NWIS qualification codes the way peakfq 8.1.0's ``siteQT``
        does (:func:`flowfreq.workflow.peak_code_kwargs`): code 7 peaks are
        historic; 3 and O peaks, and 6 and C (regulated/urban) peaks, are
        removed and their years treated as missing data; 4/8 peaks are
        less-than/greater-than intervals. A gage with none of these codes is
        fitted exactly as before. With them the result changes, and needs
        ``method="ema"`` and no ``historical_peaks`` (a ``ValueError``
        otherwise). Code 4/8 peaks are fitted as censored intervals
        (``Bulletin17C``'s ``interval_peaks``), as peakfq fits them. Raises
        :class:`flowfreq.psf_convert.UnsupportedSpecError`, naming the years,
        when the native engine cannot express the coded record exactly (a
        code 7 peak that is also code 4 or 8).
        ``False`` restores the old behaviour -- every peak fitted as an
        exact systematic value, with ignored codes logged. It was opt-in
        (default ``False``) before.
    allow_regulated : bool, default False
        The gage is screened for regulation (:mod:`flowfreq.regulation`:
        the packaged GAGES-II screen plus code 6 in the downloaded peaks)
        before fitting. A gage classed ``regulated`` raises
        :class:`flowfreq.regulation.RegulatedRecordError` unless this is
        ``True``; the override is recorded in ``site_classification``. A gage
        with no classification proceeds, with a warning.

    Returns
    -------
    dict
        Keys ``gage``, ``analysis``, ``results``, ``figures``,
        ``report_path`` and ``site_classification`` (the screen's provenance).
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

    from .regulation import classify_site, require_unregulated

    codes = (
        gage.peak_data["qualification_code"].tolist()
        if "qualification_code" in gage.peak_data.columns
        else None
    )
    site_classification = require_unregulated(
        classify_site(site_no, codes), allow_regulated=allow_regulated
    )

    logger.info("Running Bulletin 17C analysis (method=%s)...", method.upper())

    analysis = _gage_analysis(
        gage.peak_data,
        site_no,
        method=method,
        regional_skew=regional_skew,
        regional_skew_mse=regional_skew_mse,
        historical_peaks=historical_peaks,
        apply_peak_codes=apply_peak_codes,
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
        "site_classification": site_classification,
    }


# Read from the installed distribution rather than repeated here as a
# literal. The literal drifted every time: it read "0.3.0" through the 0.4.0
# release (docs/PHASE1_RUNBOOK.md still records it saying so) and "0.4.0"
# through 0.5.0, 0.6.0 and 0.6.1, so `flowfreq.__version__` has spent most of
# this project's life reporting a version that was not the one installed.
# pyproject.toml is the single source of truth; this reads what pip actually
# installed from it, which cannot disagree.
try:
    __version__ = _installed_version("flowfreq")
except PackageNotFoundError:  # pragma: no cover - only in an uninstalled checkout
    # No dist-info: someone is importing from a source tree they never
    # installed. Saying so is better than naming a version that may be wrong.
    __version__ = "unknown"

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
    # QPPQ daily-series transfer
    "FlowDurationCurve",
    "QppqResult",
    "qppq",
    "seasonal_curves",
    "center_of_timing",
    "estimate_donor_lag",
    "apply_lag",
    "rank_donors",
    "performance",
    "loocv_qppq",
    "SNOWMELT_SEASONS",
    "MONTHLY_SEASONS",
    # Transposition to an ungaged site
    "RegressionExponents",
    "TransposedResults",
    "TranspositionProvenance",
    "transpose_frequency",
    "transpose_duration",
    "transpose_low_flow",
    "PROBABILITY_KINDS",
    "DEFAULT_AREA_RATIO_RANGE",
    "LOW_FLOW_AREA_RATIO_RANGE",
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
    "flow_duration_curve",
    "DEFAULT_EXCEEDANCE_PCT",
    "separate_baseflow",
    "monthly_flow_summary",
    "seasonal_flow_summary",
    "BASEFLOW_METHODS",
    "diel_variation",
    "diel_variation_summary",
    # Sub-daily (instantaneous-series) metrics
    "DEFAULT_MAX_GAP_MULTIPLE",
    "INTERVAL_LABELS",
    "KNOWN_VALUE_UNITS",
    "TRUE_ZERO_UNITS",
    "circular_hour_statistics",
    "daily_extreme_timing",
    "extreme_timing_summary",
    "ramping_rates",
    "ramping_rate_summary",
    # USGS data retrieval
    "USGSgage",
    "USGSGage",  # Alias for backwards compatibility
    "GageAttributes",
    "NoInstantaneousDataError",
    "fetch_nwis_peaks",
    "fetch_nwis_batch",
    "join_flow_and_stage",
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
