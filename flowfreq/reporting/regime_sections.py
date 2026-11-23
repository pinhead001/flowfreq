"""Report sections describing a gage's flow regime: ``monthly_statistics``,
``summary_hydrograph`` and ``baseflow``.

Computed from the daily mean flows a section carries inline, as :mod:`.daily`'s kinds are,
and from :mod:`flowfreq.regime`.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .blocks import Block, Figure, Heading, Paragraph, Table
from .daily import _STATION_OPTIONS, _daily_frame, _drainage_area, _label, _record_table
from .sections import SECTION_KINDS, SectionContent, SectionKind, _figure_png, _version, cfs
from .spec import ReportError

#: Caption of the ``monthly_statistics`` section's main table.
MONTHLY_CAPTION = "Statistics of the monthly mean flows."
#: Caption of the ``summary_hydrograph`` section's table.
SUMMARY_HYDROGRAPH_CAPTION = "Daily mean flow by date, for the 1st and 15th of each month."
#: Caption of the ``baseflow`` section's main table.
BASEFLOW_CAPTION = "Baseflow index and mean flows for the complete water years."

_SLOTO_CROUSE_1996 = (
    "Sloto, R.A., and Crouse, M.Y., 1996, HYSEP: A computer program for streamflow hydrograph "
    "separation and analysis: U.S. Geological Survey Water-Resources Investigations Report "
    "96-4040, 46 p."
)
#: What a report says of each baseflow separation method, and the reference it cites. The
#: keys are :data:`flowfreq.regime.BASEFLOW_METHODS`, in its order.
_METHOD_WORDS: Dict[str, Tuple[str, str]] = {
    "ih_smoothed_minima": (
        "The smoothed-minima method (Institute of Hydrology, 1980) takes the lowest flow of "
        "each 5-day block, keeps it as a turning point where 0.9 times it is no more than the "
        "lowest flow of either neighboring block, and joins the turning points with straight "
        "lines, never above the streamflow.",
        "Institute of Hydrology, 1980, Low flow studies report: Wallingford, United Kingdom, "
        "Institute of Hydrology.",
    ),
    "ih_smoothed_minima_log": (
        "The turning points of the smoothed-minima method (Institute of Hydrology, 1980: the "
        "lowest flow of each 5-day block, where 0.9 times it is no more than the lowest flow "
        "of either neighboring block) were joined with straight lines through the logarithms "
        "of their flows, never above the streamflow.",
        "Institute of Hydrology, 1980, Low flow studies report: Wallingford, United Kingdom, "
        "Institute of Hydrology.",
    ),
    "hysep_fixed": (
        "The HYSEP fixed-interval method (Sloto and Crouse, 1996) gives every day of each "
        "interval the lowest flow in it.",
        _SLOTO_CROUSE_1996,
    ),
    "hysep_sliding": (
        "The HYSEP sliding-interval method (Sloto and Crouse, 1996) gives each day the lowest "
        "flow in the interval centered on it.",
        _SLOTO_CROUSE_1996,
    ),
    "hysep_local_minimum": (
        "The HYSEP local-minimum method (Sloto and Crouse, 1996) keeps each day that is the "
        "lowest in the interval centered on it and joins those days with straight lines, "
        "never above the streamflow.",
        _SLOTO_CROUSE_1996,
    ),
    "lyne_hollick": (
        "The recursive digital filter of Lyne and Hollick (1979) was passed over the record "
        "three times (forward, backward, forward) with a filter parameter of 0.925.",
        "Lyne, V., and Hollick, M., 1979, Stochastic time-variable rainfall-runoff modelling: "
        "Institute of Engineers Australia National Conference, p. 89-93.",
    ),
}

#: Calendar months in water-year order.
_MONTHS: Tuple[Tuple[int, str], ...] = (
    (10, "October"),
    (11, "November"),
    (12, "December"),
    (1, "January"),
    (2, "February"),
    (3, "March"),
    (4, "April"),
    (5, "May"),
    (6, "June"),
    (7, "July"),
    (8, "August"),
    (9, "September"),
)
#: A month counts when at least this fraction of its days has a value.
_MONTH_COMPLETE = 0.9
#: Percentiles (of non-exceedance) a summary hydrograph gives for each date.
_PERCENTILES = [10, 25, 50, 75, 90]


def _valid_frame(owner: str, options: Dict[str, Any]) -> pd.DataFrame:
    """The daily record without its negative values, which are counted as missing."""
    frame = _daily_frame(owner, options)
    valid = frame[frame["flow_cfs"] >= 0]
    if valid.empty:
        raise ReportError(f"{owner}: the daily record has no value of zero or more")
    return valid


def _negative_note(owner: str, options: Dict[str, Any], frame: pd.DataFrame) -> List[str]:
    """A sentence counting the negative values left out; nothing when there were none."""
    dropped = len(_daily_frame(owner, options)) - len(frame)
    if not dropped:
        return []
    return [f"{dropped:,} negative value(s) in the record were counted as missing."]


def _flag(owner: str, name: str, value: Any) -> bool:
    if not isinstance(value, bool):
        raise ReportError(f"{owner}: {name} must be true or false")
    return value


def _flow_axis(ax: Any, low: float, high: float) -> bool:
    """Make the flow axis log, with plain numbers at 1, 2 and 5 times each power of ten.

    It stays linear, and False is returned, where the flows reach zero or span less than a
    factor of five, over which such an axis could have a single label.
    """
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    if low <= 0 or high < 5.0 * low:
        return False
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(subs=(1.0, 2.0, 5.0)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else f"{v:g}"))
    ax.yaxis.set_minor_formatter(NullFormatter())
    return True


# -- monthly statistics ----------------------------------------------------------------------


def _monthly_rows(complete: pd.DataFrame) -> List[Tuple[str, ...]]:
    """A table row for each month, from the complete months of ``monthly_flow_summary``."""
    groups = {m: complete[complete["month"] == m] for m, _ in _MONTHS}
    volumes = {
        m: float(g["mean_flow_cfs"].mean() * g["days_in_month"].mean())
        for m, g in groups.items()
        if len(g)
    }
    annual = sum(volumes.values()) if len(volumes) == len(_MONTHS) else 0.0
    rows: List[Tuple[str, ...]] = []
    for month, name in _MONTHS:
        group = groups[month]
        if group.empty:
            rows.append((name, "0", *["-"] * 6))
            continue
        means = group["mean_flow_cfs"]
        values, years = means.to_numpy(dtype=float), group["year"].to_numpy(dtype=int)
        rows.append(
            (
                name,
                f"{len(group)}",
                cfs(means.mean()),
                cfs(means.median()),
                cfs(means.std(ddof=1)),
                f"{cfs(values.max())} ({years[values.argmax()]})",
                f"{cfs(values.min())} ({years[values.argmin()]})",
                f"{100.0 * volumes[month] / annual:.1f}" if annual > 0 else "-",
            )
        )
    return rows


def _monthly_figure(complete: pd.DataFrame, label: str) -> bytes:
    """Box plots of the monthly mean flows, October to September."""
    import matplotlib.pyplot as plt

    means = [complete.loc[complete["month"] == m, "mean_flow_cfs"].to_numpy() for m, _ in _MONTHS]
    present = [i for i, values in enumerate(means) if len(values)]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.boxplot(
        [means[i] for i in present],
        positions=present,
        widths=0.6,
        medianprops={"color": "black"},
        flierprops={"markersize": 4},
    )
    ax.set_xlim(-0.6, len(_MONTHS) - 0.4)
    ax.set_xticks(range(len(_MONTHS)))
    ax.set_xticklabels([name[:3] for _, name in _MONTHS])
    _flow_axis(ax, min(means[i].min() for i in present), max(means[i].max() for i in present))
    ax.set_ylabel("Monthly mean flow (cfs)")
    if label:
        ax.set_title(label)
    ax.grid(True, axis="y", which="both", alpha=0.3)
    return _figure_png(fig)


def _annual_monthly_table(monthly: pd.DataFrame) -> Table:
    """Monthly mean flows, one row to a water year, a dash where the month is incomplete.

    A water year with no complete month has no row.
    """
    cells = {
        (int(row["year"]), int(row["month"])): row["mean_flow_cfs"]
        for _, row in monthly[monthly["complete"]].iterrows()
    }
    return Table(
        ("Water year", *(name[:3] for _, name in _MONTHS)),
        tuple(
            (f"{year}", *(cfs(cells.get((year, month), np.nan)) for month, _ in _MONTHS))
            for year in sorted({year for year, _ in cells})
        ),
        ("right",) * (len(_MONTHS) + 1),
        "Monthly mean flows by water year (cfs).",
    )


def _build_monthly_statistics(given: Dict[str, Any]) -> SectionContent:
    from ..regime import monthly_flow_summary

    owner = "monthly_statistics"
    options = dict(given)
    annual_table = _flag(owner, "annual_table", options.get("annual_table", False))
    frame = _valid_frame(owner, options)
    monthly = monthly_flow_summary(frame, min_completeness_frac=_MONTH_COMPLETE)
    complete = monthly[monthly["complete"]]
    if complete.empty:
        raise ReportError(
            f"{owner}: the record has no complete month (at least "
            f"{100 * _MONTH_COMPLETE:g} percent of its days with a value)"
        )
    label = _label(options)
    partial = int((~monthly["complete"]).sum())
    method = [
        "The mean of the daily mean flows was taken for each calendar month with a value on at "
        f"least {100 * _MONTH_COMPLETE:g} percent of its days. The statistics are of those "
        "monthly means, one for each year; the year beside a maximum or minimum is its water "
        "year (October 1 to September 30). The percent of annual volume is the month's mean "
        "flow times its length, as a share of the twelve.",
        f"The computation was made with flowfreq {_version()}.",
    ]
    if partial:
        method.append(f"{partial} month(s) with too few days of record, or none, were left out.")
    if complete.groupby("month").size().nunique() > 1:
        method.append("The months do not all rest on the same number of years.")
    method.extend(_negative_note(owner, options, frame))
    blocks: List[Block] = [
        _record_table(options, frame),
        Paragraph(" ".join(method)),
        Table(
            (
                "Month",
                "Years",
                "Mean (cfs)",
                "Median (cfs)",
                "Standard deviation (cfs)",
                "Maximum (cfs) and water year",
                "Minimum (cfs) and water year",
                "Percent of annual volume",
            ),
            tuple(_monthly_rows(complete)),
            ("left", *["right"] * 7),
            MONTHLY_CAPTION,
        ),
        Figure(
            _monthly_figure(complete, label),
            "Monthly mean flows. The box spans the middle half of the years and the line in it "
            "is the median; the whiskers reach the farthest value within 1.5 times the box's "
            "height, and the circles are values beyond them.",
        ),
    ]
    if annual_table:
        blocks.append(Heading(2, "Monthly mean flows by year"))
        blocks.append(_annual_monthly_table(monthly))
    title = f"Monthly flow statistics at {label}" if label else "Monthly flow statistics"
    return SectionContent(title=title, blocks=blocks, options=options, references=[])


# -- summary hydrograph ----------------------------------------------------------------------


def _build_summary_hydrograph(given: Dict[str, Any]) -> SectionContent:
    from ..hydrograph import Hydrograph

    owner = "summary_hydrograph"
    options = dict(given)
    frame = _valid_frame(owner, options)
    label = _label(options)
    stats = Hydrograph.get_summary_stats(frame, _PERCENTILES).set_index("day_of_water_year")
    days = Hydrograph._compute_dowy_series(pd.DatetimeIndex(frame.index))
    counts = pd.Series(days).value_counts()
    columns = ("min", *(f"p{p}" for p in _PERCENTILES), "max")
    rows: List[Tuple[str, ...]] = []
    for month, name in _MONTHS:
        for day in (1, 15):
            number = Hydrograph.day_of_water_year(month, day)
            if number in stats.index:
                cells = tuple(cfs(stats.at[number, column]) for column in columns)
            else:
                cells = ("-",) * len(columns)
            rows.append((f"{name[:3]} {day}", f"{int(counts.get(number, 0))}", *cells))

    if all(row[1] == "0" for row in rows):
        raise ReportError(f"{owner}: the record has no value on the 1st or 15th of any month")

    method = [
        "For each date of the year, the daily mean flows of every year with a value on that "
        "date were ranked and their percentiles taken by linear interpolation. They are "
        "percentiles of non-exceedance: the 10th percentile is the flow that 10 percent of the "
        "years were at or below on that date, which a flow-duration table would call q90. A "
        "curve joins one date's statistic to the next and is not the flow of any one year.",
        f"The computation was made with flowfreq {_version()}.",
    ]
    leap = int(counts.get(Hydrograph.day_of_water_year(2, 29), 0))
    if leap:
        method.append(f"February 29 rests on {leap} years, the leap years with a value on it.")
    method.extend(_negative_note(owner, options, frame))
    figure = Hydrograph.plot_summary_hydrograph(
        frame,
        site_name=options.get("site_name"),
        site_no=str(options["site_no"]).strip() if options.get("site_no") else None,
        percentiles=_PERCENTILES,
        yscale="log" if float(frame["flow_cfs"].min()) > 0 else "linear",
    )
    blocks: List[Block] = [
        _record_table(options, frame),
        Paragraph(" ".join(method)),
        Table(
            (
                "Date",
                "Years",
                "Minimum (cfs)",
                *("Median (cfs)" if p == 50 else f"{p}th percentile (cfs)" for p in _PERCENTILES),
                "Maximum (cfs)",
            ),
            tuple(rows),
            ("left", *["right"] * (len(columns) + 1)),
            SUMMARY_HYDROGRAPH_CAPTION,
        ),
        Figure(
            _figure_png(figure),
            "Summary hydrograph: the median, the middle half and the middle 80 percent of the "
            "daily mean flows on each date, with the smallest and largest of record.",
        ),
    ]
    title = f"Summary hydrograph at {label}" if label else "Summary hydrograph"
    return SectionContent(title=title, blocks=blocks, options=options, references=[])


# -- baseflow --------------------------------------------------------------------------------


def _typical_year(means: pd.Series) -> int:
    """The year whose mean flow is nearest the median of them, the earliest on a tie."""
    means = means.sort_index()
    return int(means.index[int((means - means.median()).abs().to_numpy().argmin())])


def _index(value: float) -> str:
    return f"{value:.3f}" if math.isfinite(value) else "-"


def _separated_figure(year: int, days: pd.DataFrame, label: str) -> bytes:
    """One water year's streamflow with its baseflow shaded, broken where a day has no value."""
    import matplotlib.pyplot as plt
    from matplotlib.dates import DateFormatter, MonthLocator, date2num

    first, last = pd.Timestamp(year - 1, 10, 1), pd.Timestamp(year, 9, 30)
    days = days.reindex(pd.date_range(first, last, freq="D"))
    fig, ax = plt.subplots(figsize=(8, 5))
    lowest = float(days[["flow", "base"]].min().min())
    log = _flow_axis(ax, lowest, float(days["flow"].max()))
    floor = 0.8 * lowest if log else 0.0
    ax.fill_between(
        days.index, floor, days["base"], color="C0", alpha=0.4, linewidth=0, label="Baseflow"
    )
    ax.plot(days.index, days["flow"], color="black", linewidth=0.9, label="Streamflow")
    ax.set_ylim(bottom=floor)  # after the curves: a limit set first stops the autoscaling
    ax.set_xlim(date2num(first), date2num(last))
    ax.xaxis.set_major_locator(MonthLocator())
    ax.xaxis.set_major_formatter(DateFormatter("%b"))
    ax.set_xlabel(f"Water year {year}")
    ax.set_ylabel("Daily mean flow (cfs)")
    if label:
        ax.set_title(label)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper right")
    return _figure_png(fig)


def _index_figure(annual: pd.DataFrame, pooled: float, label: str) -> bytes:
    """The baseflow index of each complete year, with the pooled value."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(annual["year"], annual["baseflow_index"], color="C0", alpha=0.7)
    ax.axhline(pooled, color="black", linewidth=1, label=f"All complete years ({pooled:.3f})")
    ax.set_ylim(0, 1)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("Water year")
    ax.set_ylabel("Baseflow index")
    if label:
        ax.set_title(label)
    ax.grid(True, axis="y", alpha=0.3)
    ax.legend(loc="upper right")
    return _figure_png(fig)


def _build_baseflow(given: Dict[str, Any]) -> SectionContent:
    from ..core import assign_year_label
    from ..regime import BASEFLOW_METHODS, FlowRegime, _hysep_interval_days

    owner = "baseflow"
    options = dict(given)
    method = str(options.get("method", "ih_smoothed_minima"))
    if method not in BASEFLOW_METHODS:
        raise ReportError(f"{owner}: method must be one of {', '.join(BASEFLOW_METHODS)}")
    min_days = options.get("min_days", 350)
    if isinstance(min_days, bool) or not isinstance(min_days, int) or not 1 <= min_days <= 366:
        raise ReportError(f"{owner}: min_days must be a whole number of days from 1 to 366")
    annual_table = _flag(owner, "annual_table", options.get("annual_table", True))
    frame = _valid_frame(owner, options)
    area = _drainage_area(owner, options)
    hysep = method.startswith("hysep")
    if hysep and area is None:
        raise ReportError(
            f"{owner}: {method} sets its interval from the drainage area: give drainage_area_sqmi"
        )

    regime = FlowRegime(frame, min_days=min_days, baseflow_method=method, drainage_area_sqmi=area)
    annual = regime.annual
    used = annual[annual["baseflow_complete"]]
    if used.empty:
        raise ReportError(
            f"{owner}: the record has no complete water year (at least {min_days} days with "
            "both a streamflow and a separated baseflow)"
        )
    pooled = float(regime.summary()["baseflow_index"])
    if not math.isfinite(pooled):
        raise ReportError(f"{owner}: every flow of the complete water years is zero")
    separated = pd.DataFrame({"flow": frame["flow_cfs"], "base": regime.baseflow_series})
    days = separated.dropna()
    days["year"] = assign_year_label(pd.DatetimeIndex(days.index), "water")
    days = days[days["year"].isin(used["year"])]
    label = _label(options)

    words, citation = _METHOD_WORDS[method]
    left_out = len(annual) - len(used)
    text = [
        "Baseflow was separated from the daily mean flows. " + words,
    ]
    if hysep and area is not None:
        text.append(
            f"The drainage area of {area:,.4g} square miles gives an interval of "
            f"{_hysep_interval_days(area)} days."
        )
    text.append(
        "The baseflow index is the volume of baseflow divided by the volume of streamflow over "
        f"the water years (October 1 to September 30) with at least {min_days} days that have "
        "both."
    )
    text.append(f"The computation was made with flowfreq {_version()}.")
    if left_out:
        text.append(f"{left_out} water year(s) with too few days, or none, were left out.")
    text.extend(_negative_note(owner, options, frame))
    cautions = (
        "The separation is a rule applied to the hydrograph, not a measurement of groundwater "
        "discharge. Baseflow indexes from different methods are not comparable, so compare "
        "this one only with others by the same method. Reservoir regulation, diversion, return "
        "flow and snowmelt all raise or lower the index without any change in groundwater."
    )
    rows: List[Tuple[str, str]] = [
        ("Separation method", method),
        ("Complete water years", f"{len(used)}"),
        ("Mean streamflow (cfs)", cfs(days["flow"].mean())),
        ("Mean baseflow (cfs)", cfs(days["base"].mean())),
        ("Baseflow index", _index(pooled)),
    ]
    if area is not None:
        inches = 365.25 * 86400.0 * 12.0 / (area * 5280.0**2)
        rows.extend(
            [
                ("Mean streamflow (inches per year)", f"{days['flow'].mean() * inches:.2f}"),
                ("Mean baseflow (inches per year)", f"{days['base'].mean() * inches:.2f}"),
            ]
        )
    year = _typical_year(days.groupby("year")["flow"].mean())
    blocks: List[Block] = [
        _record_table(options, frame),
        Paragraph(" ".join(text)),
        Paragraph(cautions),
        Table(("Item", "Value"), tuple(rows), ("left", "right"), BASEFLOW_CAPTION),
        Figure(
            _separated_figure(year, separated, label),
            f"Streamflow and separated baseflow in water year {year}, the complete year whose "
            "mean flow is nearest the median.",
        ),
        Figure(
            _index_figure(used, pooled, label),
            "Baseflow index of each complete water year.",
        ),
    ]
    if annual_table:
        blocks.append(Heading(2, "Annual baseflow index"))
        blocks.append(
            Table(
                ("Water year", "Days used", "Baseflow index"),
                tuple(
                    (f"{y}", f"{n}", _index(float(index)))
                    for y, n, index in zip(
                        used["year"],
                        days.groupby("year").size().reindex(used["year"]),
                        used["baseflow_index"],
                    )
                ),
                ("right", "right", "right"),
                "Baseflow index by water year.",
            )
        )
    title = f"Baseflow at {label}" if label else "Baseflow"
    return SectionContent(title=title, blocks=blocks, options=options, references=[citation])


_KINDS: Sequence[SectionKind] = (
    SectionKind(
        name="monthly_statistics",
        summary=(
            "Monthly flow statistics at a gage from its daily mean flows: the mean, median, "
            "standard deviation and range of each month's mean flow, its share of the annual "
            "volume, and box plots."
        ),
        required=(),
        optional=(*_STATION_OPTIONS, "annual_table"),
        build=_build_monthly_statistics,
    ),
    SectionKind(
        name="summary_hydrograph",
        summary=(
            "Summary hydrograph at a gage from its daily mean flows: the median, quartiles, "
            "10th and 90th percentiles and range of the flow on each date of the year."
        ),
        required=(),
        optional=_STATION_OPTIONS,
        build=_build_summary_hydrograph,
    ),
    SectionKind(
        name="baseflow",
        summary=(
            "Baseflow at a gage from its daily mean flows, by one named separation method: the "
            "baseflow index for the record and each year, mean streamflow and baseflow, and a "
            "separated hydrograph for a typical year."
        ),
        required=(),
        optional=(*_STATION_OPTIONS, "method", "min_days", "annual_table"),
        build=_build_baseflow,
    ),
)
SECTION_KINDS.update({kind.name: kind for kind in _KINDS})

__all__ = ["BASEFLOW_CAPTION", "MONTHLY_CAPTION", "SUMMARY_HYDROGRAPH_CAPTION"]
