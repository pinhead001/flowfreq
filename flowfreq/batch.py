"""
flowfreq.batch - Multi-site batch processing for flood frequency analysis
"""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import pandas as pd

from .core import PeakRecord
from .engine import STANDARD_RETURN_PERIODS, B17CEngine
from .usgs import fetch_nwis_batch

PeakInput = Union[PeakRecord, Mapping[str, Any]]


def _as_peak_records(records: Sequence[PeakInput]) -> List[PeakRecord]:
    """Accept ``PeakRecord`` objects or ``fetch_nwis_peaks``-style dicts.

    ``usgs.fetch_nwis_peaks`` (and so ``fetch_nwis_batch``) returns plain dicts
    ``{"year", "flow", "source"}`` -- its documented public contract -- while
    ``B17CEngine.fit`` reads attributes. Converting here is what lets
    :func:`analyze_sites` work on real fetched data at all.
    """
    out: List[PeakRecord] = []
    for r in records:
        if isinstance(r, PeakRecord):
            out.append(r)
        elif isinstance(r, Mapping):
            out.append(
                PeakRecord(
                    year=int(r["year"]),
                    flow=None if r.get("flow") is None else float(r["flow"]),
                    source=r.get("source"),
                )
            )
        else:
            raise TypeError(f"Expected PeakRecord or dict peak record, got {type(r).__name__}")
    return out


def run_multi_site(
    data: Mapping[str, Sequence[PeakInput]],
    return_periods: Sequence[float] = STANDARD_RETURN_PERIODS,
) -> Dict[str, Dict[str, Any]]:
    """
    Run flood frequency analysis for multiple sites.

    Parameters
    ----------
    data : dict
        Mapping of site_no to a list of :class:`PeakRecord`, or of the plain
        ``{"year", "flow", "source"}`` dicts :func:`flowfreq.usgs.fetch_nwis_peaks`
        returns
    return_periods : sequence of float
        Return periods to compute (default: standard set)

    Returns
    -------
    dict
        Mapping of site_no to analysis results containing:
        - params: (mu, sigma, skew) tuple
        - n: sample size
        - quantiles: dict of return period to flow
        - ci: dict of return period to confidence interval dict
    """
    results = {}

    for site, records in data.items():
        try:
            engine = B17CEngine()
            engine.fit(_as_peak_records(records))
            results[site] = {
                "params": engine.params,
                "n": engine.n,
                "quantiles": engine.quantiles(return_periods),
                "ci": engine.quantiles_with_ci(return_periods),
            }
        except Exception as e:
            results[site] = {"error": str(e)}

    return results


def analyze_sites(
    site_nos: List[str],
    return_periods: Sequence[float] = STANDARD_RETURN_PERIODS,
    workers: int = 6,
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, str]]:
    """
    Fetch data and run analysis for multiple USGS sites.

    Parameters
    ----------
    site_nos : list of str
        USGS site numbers
    return_periods : sequence of float
        Return periods to compute
    workers : int
        Number of parallel workers for data fetching

    Returns
    -------
    tuple
        (analysis_results, fetch_errors) where:
        - analysis_results: dict of site_no to analysis results
        - fetch_errors: dict of site_no to fetch error message
    """
    # Fetch data in parallel
    data, fetch_errors = fetch_nwis_batch(site_nos, workers=workers)

    # Run analysis on fetched data
    analysis_results = run_multi_site(data, return_periods)

    return analysis_results, fetch_errors


def batch_summary_table(
    results: Dict[str, Dict[str, Any]],
    return_periods: Sequence[float] = (10, 50, 100),
) -> "pd.DataFrame":
    """
    Generate a summary table of multi-site analysis results.

    Parameters
    ----------
    results : dict
        Output from run_multi_site or analyze_sites
    return_periods : sequence of float
        Return periods to include in table

    Returns
    -------
    pd.DataFrame
        Summary table with sites as rows and quantiles as columns
    """

    rows = []
    for site, result in results.items():
        if "error" in result:
            row = {"Site": site, "Error": result["error"]}
        else:
            mu, sigma, skew = result["params"]
            row = {
                "Site": site,
                "n": result["n"],
                "Mean (log)": mu,
                "Std (log)": sigma,
                "Skew": skew,
            }
            for T in return_periods:
                if T in result["quantiles"]:
                    row[f"Q{int(T)}"] = result["quantiles"][T]

        rows.append(row)

    return pd.DataFrame(rows)
