"""Report sections computed from a gage's daily mean flows: ``low_flow`` and ``flow_duration``.

Both take the daily record inline, as ``daily``: a first date and one value per calendar day
after it, ``null`` where the record has no value. Given a ``site_no`` and no ``daily``, the
record is downloaded and written into the options, so the saved specification rebuilds the
section without the network, as a ``gaged_ffa`` section's peaks do.
"""

from __future__ import annotations

import datetime as _dt
import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .blocks import Block, Figure, Heading, Paragraph, Table
from .sections import SECTION_KINDS, SectionContent, SectionKind, _figure_png, _version, cfs
from .spec import ReportError

#: The reference a low-flow frequency section cites.
RIGGS_1972 = (
    "Riggs, H.C., 1972, Low-flow investigations: U.S. Geological Survey Techniques of "
    "Water-Resources Investigations, book 4, chap. B1, 18 p."
)
#: The reference a flow-duration section cites.
SEARCY_1959 = (
    "Searcy, J.K., 1959, Flow-duration curves, Manual of hydrology: Part 2, Low-flow "
    "techniques: U.S. Geological Survey Water-Supply Paper 1542-A, 33 p."
)
#: Captions of the two sections' main tables.
LOW_FLOW_CAPTION = "Low-flow quantiles with 90-percent confidence limits."
FLOW_DURATION_CAPTION = "Flow-duration statistics for the period used."

_DAILY_KEYS = ("start", "flow_cfs", "source", "retrieved")
_YEAR_WORDS = {
    "climatic": "climatic year (April 1 to March 31)",
    "water": "water year (October 1 to September 30)",
    "calendar": "calendar year",
}
_STATION_OPTIONS = (
    "site_no",
    "site_name",
    "drainage_area_sqmi",
    "daily",
    "start_date",
    "end_date",
)


def download_daily(
    site_no: str, start_date: Optional[str] = None, end_date: Optional[str] = None
) -> Dict[str, Any]:
    """A gage's daily mean flows and site attributes, in the form a spec stores.

    Parameters
    ----------
    site_no : str
        USGS site number.
    start_date, end_date : str, optional
        ISO dates bounding the request. Both default to the period of record.

    Returns
    -------
    dict
        ``daily`` (the inline record: ``start``, the first day as ``YYYY-MM-DD``;
        ``flow_cfs``, one value per calendar day from ``start`` with ``None`` where the
        record has none; ``source``; ``retrieved``), ``site_name`` and
        ``drainage_area_sqmi`` (from the download, else flowfreq's gage catalog; either may
        be None).

    Raises
    ------
    ReportError
        The download failed or returned no values.
    """
    import requests

    from ..usgs import USGSgage

    gage = USGSgage(site_no)
    try:
        frame = gage.download_daily_flow(start_date=start_date, end_date=end_date)
    except (requests.RequestException, ValueError, KeyError) as exc:
        raise ReportError(
            f"daily flows for USGS {gage.site_no} could not be downloaded: {exc}"
        ) from exc
    series = pd.Series(
        frame["flow_cfs"].to_numpy(dtype=float), index=pd.DatetimeIndex(frame.index)
    ).dropna()
    if series.empty:
        raise ReportError(f"USGS {gage.site_no} has no daily flows")
    series = series[~series.index.duplicated()].sort_index()
    days = pd.date_range(series.index[0].normalize(), series.index[-1].normalize(), freq="D")
    values = series.reindex(days).to_numpy(dtype=float)
    return {
        "daily": {
            "start": days[0].date().isoformat(),
            "flow_cfs": [None if math.isnan(v) else float(v) for v in values],
            "source": f"USGS Water Data API daily mean flow, site {gage.site_no}",
            "retrieved": _dt.date.today().isoformat(),
        },
        **_station_attributes(gage),
    }


def _station_attributes(gage: Any) -> Dict[str, Any]:
    """The station's name and drainage area: the gage's own, else the shipped catalog's.

    The daily download does not always carry them, and the catalog needs no request. It
    lists gages with annual peaks, so either may still be None.
    """
    name, area = gage.site_name, gage.drainage_area
    if name is None or area is None:
        from ..catalog import load_catalog

        catalog = load_catalog()
        rows = catalog[catalog["site_no"].astype(str) == str(gage.site_no)]
        if not rows.empty:
            row = rows.iloc[0]
            if name is None and pd.notna(row["site_name"]):
                name = str(row["site_name"])
            if area is None and pd.notna(row["drainage_area_sqmi"]):
                area = float(row["drainage_area_sqmi"])
    return {"site_name": name, "drainage_area_sqmi": area}


def _date(owner: str, name: str, value: Any) -> pd.Timestamp:
    try:
        stamp = pd.Timestamp(_dt.date.fromisoformat(str(value)))
    except ValueError as exc:
        raise ReportError(f"{owner}: {name} must be a date as YYYY-MM-DD, got {value!r}") from exc
    return stamp


def _drainage_area(owner: str, options: Dict[str, Any]) -> Optional[float]:
    """The drainage area as a number; None when not given."""
    value = options.get("drainage_area_sqmi")
    if value is None:
        return None
    try:
        area = float(value)
    except (TypeError, ValueError):
        area = math.nan
    if isinstance(value, bool) or not (math.isfinite(area) and area > 0):
        raise ReportError(f"{owner}: drainage_area_sqmi must be a positive number")
    return area


def _daily_frame(owner: str, options: Dict[str, Any]) -> pd.DataFrame:
    """The inline daily record as a frame indexed by date, days without a value absent.

    Downloads the record first when the options carry none, and writes it into them.
    """
    _drainage_area(owner, options)
    site_no = str(options["site_no"]).strip() if options.get("site_no") else None
    if options.get("daily") is None:
        if site_no is None:
            raise ReportError(f"{owner}: give site_no (to download the daily flows) or daily")
        fetched = download_daily(site_no, options.get("start_date"), options.get("end_date"))
        options["daily"] = fetched["daily"]
        for key in ("site_name", "drainage_area_sqmi"):
            if options.get(key) is None and fetched.get(key) is not None:
                options[key] = fetched[key]
    daily = options["daily"]
    if not isinstance(daily, dict):
        raise ReportError(f"{owner}: daily must be an object with start and flow_cfs")
    extra = sorted(set(daily) - set(_DAILY_KEYS))
    if extra:
        raise ReportError(
            f"{owner}: unknown daily key(s) {', '.join(extra)}; allowed: {', '.join(_DAILY_KEYS)}"
        )
    values = daily.get("flow_cfs")
    if daily.get("start") is None or not isinstance(values, list) or not values:
        raise ReportError(f"{owner}: daily needs start (YYYY-MM-DD) and a flow_cfs list")
    start = _date(owner, "daily start", daily["start"])
    try:
        flows = np.array([np.nan if v is None else float(v) for v in values], dtype=float)
    except (TypeError, ValueError) as exc:
        raise ReportError(f"{owner}: daily flow_cfs is not numeric: {exc}") from exc
    if np.any(np.isinf(flows)):
        raise ReportError(f"{owner}: daily flow_cfs holds an infinite value")
    frame = pd.DataFrame(
        {"flow_cfs": flows}, index=pd.date_range(start, periods=len(flows), freq="D", name="date")
    )
    if options.get("start_date") is not None:
        frame = frame[frame.index >= _date(owner, "start_date", options["start_date"])]
    if options.get("end_date") is not None:
        frame = frame[frame.index <= _date(owner, "end_date", options["end_date"])]
    frame = frame.dropna()
    if frame.empty:
        raise ReportError(f"{owner}: the daily record has no value in the period asked for")
    return frame


def _label(options: Dict[str, Any]) -> str:
    site_no = str(options["site_no"]).strip() if options.get("site_no") else None
    parts = (f"USGS {site_no}" if site_no else None, options.get("site_name"))
    return " ".join(str(x) for x in parts if x)


def _record_table(options: Dict[str, Any], frame: pd.DataFrame) -> Table:
    first, last = frame.index[0], frame.index[-1]
    span = (last - first).days + 1
    rows: List[Tuple[str, str]] = []
    if options.get("site_no"):
        rows.append(("USGS station number", str(options["site_no"]).strip()))
    if options.get("site_name"):
        rows.append(("Station name", str(options["site_name"])))
    if options.get("drainage_area_sqmi") is not None:
        rows.append(
            ("Drainage area (square miles)", f"{float(options['drainage_area_sqmi']):,.4g}")
        )
    rows.extend(
        [
            ("Period used", f"{first.date().isoformat()} to {last.date().isoformat()}"),
            ("Days with a daily mean flow", f"{len(frame):,}"),
            ("Days without one", f"{span - len(frame):,}"),
            ("Largest daily mean flow (cfs)", cfs(frame["flow_cfs"].max())),
            ("Smallest daily mean flow (cfs)", cfs(frame["flow_cfs"].min())),
        ]
    )
    source = options["daily"].get("source")
    if source:
        retrieved = options["daily"].get("retrieved")
        rows.append(
            ("Data source", str(source) + (f", retrieved {retrieved}" if retrieved else ""))
        )
    return Table(("Item", "Value"), tuple(rows), ("left", "left"), "Station and daily-flow record.")


def _probabilities(owner: str, name: str, value: Any, top: float) -> Optional[List[float]]:
    """A list of numbers strictly between 0 and ``top``, without repeats; None when unset."""
    if value is None:
        return None
    if not isinstance(value, list) or not value:
        raise ReportError(f"{owner}: {name} must be a non-empty list of numbers")
    try:
        numbers = [float(v) for v in value]
    except (TypeError, ValueError) as exc:
        raise ReportError(f"{owner}: {name} must be a list of numbers") from exc
    if any(not (math.isfinite(v) and 0.0 < v < top) for v in numbers):
        raise ReportError(f"{owner}: every {name} value must lie strictly between 0 and {top:g}")
    if len(set(numbers)) != len(numbers):
        raise ReportError(f"{owner}: {name} repeats a value")
    return numbers


# -- low-flow frequency ---------------------------------------------------------------------


def _statistic(n_day: int, non_exceedance: float) -> str:
    """``7Q10`` for the 7-day flow at non-exceedance 0.1; empty off a whole recurrence interval."""
    years = 1.0 / non_exceedance
    return f"{n_day}Q{round(years)}" if abs(years - round(years)) < 1e-9 else ""


def _low_flow_figure(fit: Any, results: Any, label: str) -> bytes:
    """Annual minima at Weibull plotting positions with the fitted curve, on a normal axis."""
    import matplotlib.pyplot as plt
    from scipy.special import ndtr, ndtri

    minima = np.sort(fit.annual_minimums["flow_cfs"].to_numpy(dtype=float))
    n = len(minima)
    positions = np.arange(1, n + 1) / (n + 1.0)
    low = max(results.p_zero + 1e-4, min(0.005, positions[0] / 2))
    grid = ndtr(np.linspace(ndtri(low), ndtri(0.995), 200))
    curve = fit.compute_quantiles(grid)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ndtri(curve["non_exceedance_prob"]), curve["flow_cfs"], color="C0", label="Fitted")
    ax.plot(ndtri(positions), minima, "o", color="black", markersize=4, label="Annual minimum")
    ticks = [p for p in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99) if p >= low]
    ax.set_xticks(ndtri(ticks))
    ax.set_xticklabels([f"{100 * p:g}" for p in ticks])
    if minima[0] > 0:
        ax.set_yscale("log")
    ax.set_xlabel("Annual non-exceedance probability (percent)")
    ax.set_ylabel(f"Annual minimum {results.n_day}-day mean flow (cfs)")
    if label:
        ax.set_title(label, parse_math=False)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right")
    return _figure_png(fig)


def _build_low_flow(given: Dict[str, Any]) -> SectionContent:
    from ..core import YEAR_TYPES
    from ..lowflow import LowFlowFrequency

    owner = "low_flow"
    options = dict(given)
    n_day = options.get("n_day", 7)
    if isinstance(n_day, bool) or not isinstance(n_day, int) or not 1 <= n_day <= 365:
        raise ReportError(f"{owner}: n_day must be a whole number of days, 1 to 365")
    year_type = str(options.get("year_type", "climatic"))
    if year_type not in YEAR_TYPES:
        raise ReportError(f"{owner}: year_type must be one of {', '.join(YEAR_TYPES)}")
    distribution = str(options.get("distribution", "lp3"))
    if distribution not in ("lp3", "lognormal"):
        raise ReportError(f"{owner}: distribution must be lp3 or lognormal")
    asked = _probabilities(owner, "non_exceedance", options.get("non_exceedance"), 1.0)
    min_days = options.get("min_days", 350)
    if isinstance(min_days, bool) or not isinstance(min_days, int) or not 1 <= min_days <= 366:
        raise ReportError(f"{owner}: min_days must be a whole number of days from 1 to 366")
    frame = _daily_frame(owner, options)

    try:
        fit = LowFlowFrequency(
            frame, n_day=n_day, year_type=year_type, min_days=min_days, distribution=distribution
        )
        results = fit.run_analysis()
        probabilities = (
            np.sort(np.array(asked))[::-1] if asked else LowFlowFrequency.STANDARD_NONEXCEEDANCE
        )
        limits = fit.compute_confidence_limits(probabilities)
    except ValueError as exc:
        raise ReportError(f"{owner}: the low-flow analysis failed: {exc}") from exc

    label = _label(options)
    every_year = fit.annual_minimums_all_years
    used = fit.annual_minimums
    left_out = len(every_year) - len(used)
    blocks: List[Block] = [_record_table(options, frame)]

    method = [
        f"The annual minimum {n_day}-day mean flow was taken for each "
        f"{_YEAR_WORDS.get(year_type, year_type)} with at least {min_days} days of record."
    ]
    if n_day > 1:
        method.append("A period that runs from one such year into the next was not used.")
    if distribution == "lp3":
        method.append(
            "A log-Pearson Type III distribution was fitted to the logarithms of those minima "
            "by the method of moments with the station skew (Riggs, 1972)."
        )
    else:
        method.append(
            "A lognormal distribution (log-Pearson Type III with zero skew) was fitted to the "
            "logarithms of those minima by the method of moments (Riggs, 1972)."
        )
    method.append(f"The computation was made with flowfreq {_version()}.")
    if left_out:
        method.append(f"{left_out} year(s) with too few days of record were left out.")
    if results.n_zero_years:
        method.append(
            f"{results.n_zero_years} of the {results.n_years} years had a minimum of zero. The "
            "distribution was fitted to the other years, and the probabilities were adjusted "
            "for the fraction of zero years (a conditional probability adjustment); a "
            "probability at or below that fraction has a flow of zero. The confidence limits "
            "treat that fraction as known, so they are narrower than the record supports."
        )
    if results.n_years < LowFlowFrequency.RECOMMENDED_MIN_YEARS:
        method.append(
            f"**The fit rests on {results.n_years} years, fewer than the "
            f"{LowFlowFrequency.RECOMMENDED_MIN_YEARS} recommended; the rarer quantiles are "
            "supported by very few observations.**"
        )
    method.append(
        "No low-outlier test and no regional skew were applied: neither exists for low flows."
    )
    blocks.append(Paragraph(" ".join(method)))

    fitted: List[Tuple[str, str]] = [
        ("Years used", f"{results.n_years}"),
        ("Years left out as incomplete", f"{left_out}"),
        ("Years with a minimum of zero", f"{results.n_zero_years}"),
        ("Mean of log10 flow", f"{results.mean_log:.4f}"),
        ("Standard deviation of log10 flow", f"{results.std_log:.4f}"),
        ("Station skew", f"{results.skew_station:.3f}"),
        ("Skew used", f"{results.skew_used:.3f}"),
    ]
    blocks.append(
        Table(
            ("Parameter", "Value"),
            tuple(fitted),
            ("left", "right"),
            f"Fitted distribution of the annual minimum {n_day}-day mean flow.",
        )
    )
    blocks.append(
        Table(
            (
                "Annual non-exceedance probability (percent)",
                "Recurrence interval (years)",
                "Statistic",
                "Flow (cfs)",
                "Lower 90-percent limit (cfs)",
                "Upper 90-percent limit (cfs)",
            ),
            tuple(
                (
                    f"{row['non_exceedance_prob'] * 100:.3g}",
                    f"{row['return_period']:.4g}",
                    _statistic(n_day, float(row["non_exceedance_prob"])),
                    cfs(row["flow_cfs"]),
                    cfs(row["lower_90pct"]),
                    cfs(row["upper_90pct"]),
                )
                for _, row in limits.iterrows()
            ),
            ("right", "right", "left", "right", "right", "right"),
            LOW_FLOW_CAPTION,
        )
    )
    blocks.append(
        Figure(
            _low_flow_figure(fit, results, label),
            f"Annual minimum {n_day}-day mean flows at Weibull plotting positions, with the "
            "fitted distribution.",
        )
    )
    if bool(options.get("annual_table", True)):
        blocks.append(Heading(2, "Annual minima"))
        blocks.append(
            Table(
                ("Year", f"Minimum {n_day}-day mean flow (cfs)", "Last day of the period"),
                tuple(
                    (
                        f"{int(row['year'])}",
                        cfs(row["flow_cfs"]),
                        "" if pd.isna(row["window_end_date"]) else str(row["window_end_date"])[:10],
                    )
                    for _, row in used.iterrows()
                ),
                ("right", "right", "left"),
                f"Annual minimum {n_day}-day mean flows used in the analysis.",
            )
        )
    title = f"Low-flow frequency at {label}" if label else "Low-flow frequency"
    return SectionContent(title=title, blocks=blocks, options=options, references=[RIGGS_1972])


# -- flow duration --------------------------------------------------------------------------


def _build_flow_duration(given: Dict[str, Any]) -> SectionContent:
    from ..hydrograph import Hydrograph
    from ..regime import flow_duration_curve

    owner = "flow_duration"
    options = dict(given)
    asked = _probabilities(owner, "exceedance_pct", options.get("exceedance_pct"), 100.0)
    frame = _daily_frame(owner, options)
    try:
        curve = flow_duration_curve(frame, asked)
    except ValueError as exc:  # no day with a flow: every value negative
        raise ReportError(f"{owner}: {exc}") from exc
    label = _label(options)
    first, last = frame.index[0], frame.index[-1]
    missing = (last - first).days + 1 - len(frame)

    method = [
        "The flow-duration curve ranks every daily mean flow in the period used, as one "
        "population (Searcy, 1959). A statistic qN is the daily mean flow equaled or exceeded "
        "on N percent of those days, interpolated between ranked values. It describes the "
        "period as recorded and is not an annual probability.",
        f"The computation was made with flowfreq {_version()}.",
    ]
    if missing:
        method.append(
            f"{missing:,} day(s) in the period have no value and are not in the ranking; if "
            "they fall mostly in one season, the curve leans away from that season's flows."
        )
    blocks: List[Block] = [
        _record_table(options, frame),
        Paragraph(" ".join(method)),
        Table(
            ("Percent of days equaled or exceeded", "Statistic", "Flow (cfs)"),
            tuple(
                (f"{row['exceedance_pct']:g}", f"q{row['exceedance_pct']:g}", cfs(row["flow_cfs"]))
                for _, row in curve.iterrows()
            ),
            ("right", "left", "right"),
            FLOW_DURATION_CAPTION,
        ),
    ]
    figure, _ = Hydrograph.plot_flow_duration_curve(
        frame,
        site_name=options.get("site_name"),
        site_no=str(options["site_no"]).strip() if options.get("site_no") else None,
        yscale="log" if float(frame["flow_cfs"].min()) > 0 else "linear",
    )
    blocks.append(Figure(_figure_png(figure), "Flow-duration curve of the daily mean flows."))
    title = f"Flow duration at {label}" if label else "Flow duration"
    return SectionContent(title=title, blocks=blocks, options=options, references=[SEARCY_1959])


_KINDS: Sequence[SectionKind] = (
    SectionKind(
        name="low_flow",
        summary=(
            "Low-flow frequency at a gage from its daily mean flows: the annual minimum "
            "n-day mean flow fitted with a log-Pearson Type III distribution (7Q10, 7Q2 and "
            "the like), with 90-percent confidence limits, the fitted curve and the annual "
            "minima."
        ),
        required=(),
        optional=(
            *_STATION_OPTIONS,
            "n_day",
            "year_type",
            "distribution",
            "non_exceedance",
            "min_days",
            "annual_table",
        ),
        build=_build_low_flow,
    ),
    SectionKind(
        name="flow_duration",
        summary=(
            "Flow duration at a gage from its daily mean flows: the flow equaled or exceeded "
            "on a given percent of days (q10, q50, q90 and the like) and the curve."
        ),
        required=(),
        optional=(*_STATION_OPTIONS, "exceedance_pct"),
        build=_build_flow_duration,
    ),
)
SECTION_KINDS.update({kind.name: kind for kind in _KINDS})

__all__ = [
    "FLOW_DURATION_CAPTION",
    "LOW_FLOW_CAPTION",
    "RIGGS_1972",
    "SEARCY_1959",
    "download_daily",
]
