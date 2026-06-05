"""
FFA export utilities for writing flood frequency analysis results to ZIP archives.
"""

import io
import logging
from typing import Any, Dict, Optional

import pandas as pd

try:
    from matplotlib.figure import Figure
except ImportError:
    Figure = None  # type: ignore[assignment,misc]

logger = logging.getLogger(__name__)


def export_ffa_to_zip(
    zf: Any,
    site_no: str,
    ffa_results: Dict[str, Any],
    freq_fig: Optional[Any] = None,
) -> None:
    """
    Add flood frequency analysis results for one site to an open ZipFile.

    Parameters
    ----------
    zf : zipfile.ZipFile
        Open ZIP archive for writing.
    site_no : str
        USGS site number used as folder prefix.
    ffa_results : dict
        Result dict from ffa_runner.run_ffa containing keys like
        ``quantile_df``, ``parameters``, ``method``, ``converged``, and
        optionally ``peak_df``.
    freq_fig : matplotlib.figure.Figure or None
        Frequency curve figure to save as PNG.
    """
    prefix = f"{site_no}/"

    # 1. Frequency curve PNG
    if freq_fig is not None:
        buf = io.BytesIO()
        freq_fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
        buf.seek(0)
        zf.writestr(f"{prefix}frequency_curve.png", buf.read())

    # 2. Frequency table CSV
    quantile_df = ffa_results.get("quantile_df")
    if quantile_df is not None:
        buf = io.StringIO()
        quantile_df.to_csv(buf, index=False)
        zf.writestr(f"{prefix}frequency_table.csv", buf.getvalue())

    # 3. LP3 parameters CSV
    params = ffa_results.get("parameters", {})
    param_row = {
        "mean_log": params.get("mean_log"),
        "std_log": params.get("std_log"),
        "skew_station": params.get("skew_station"),
        "skew_weighted": params.get("skew_weighted"),
        "skew_used": params.get("skew_used"),
        "method": ffa_results.get("method"),
        "converged": ffa_results.get("converged"),
    }
    param_df = pd.DataFrame([param_row])
    buf = io.StringIO()
    param_df.to_csv(buf, index=False)
    zf.writestr(f"{prefix}lp3_parameters.csv", buf.getvalue())

    # 4. Peak flows CSV (optional)
    peak_df = ffa_results.get("peak_df")
    if peak_df is not None:
        buf = io.StringIO()
        peak_df.to_csv(buf, index=False)
        zf.writestr(f"{prefix}peak_flows.csv", buf.getvalue())


def export_comparison_csv(
    zf: Any,
    site_results: Dict[str, Dict[str, Any]],
) -> None:
    """
    Write a cross-site comparison summary CSV to the ZIP root.

    Parameters
    ----------
    zf : zipfile.ZipFile
        Open ZIP archive for writing.
    site_results : dict[str, dict]
        Keyed by site_no. Each value has ``site_name``,
        ``drainage_area_sqmi``, and ``ffa_results``.
    """
    rows = []
    for site_no, info in site_results.items():
        ffa = info.get("ffa_results") or {}
        params = ffa.get("parameters") or {}
        quantile_df = ffa.get("quantile_df")

        # Extract 100-yr flow
        flow_100 = None
        if quantile_df is not None:
            mask = quantile_df["Return Interval (yr)"] == 100.0
            matched = quantile_df.loc[mask, "Flow (cfs)"]
            if len(matched) > 0:
                flow_100 = matched.iloc[0]

        da = info.get("drainage_area_sqmi")
        if flow_100 is not None and da is not None and da > 0:
            flow_per_da = flow_100 / da
        else:
            flow_per_da = "N/A"

        rows.append(
            {
                "Site No": site_no,
                "Site Name": info.get("site_name", ""),
                "Drainage Area (sq mi)": da,
                "100-yr Flow (cfs)": flow_100,
                "100-yr Flow/DA (cfs/sq mi)": flow_per_da,
                "Weighted Skew": params.get("skew_weighted"),
                "Method": ffa.get("method"),
                "Converged": ffa.get("converged"),
            }
        )

    df = pd.DataFrame(rows)
    buf = io.StringIO()
    df.to_csv(buf, index=False)
    zf.writestr("comparison_summary.csv", buf.getvalue())
