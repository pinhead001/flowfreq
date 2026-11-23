"""Report sections describing a gage's flow regime: ``monthly_statistics`` and
``summary_hydrograph``.

Computed from the daily mean flows a section carries inline, as :mod:`.daily`'s kinds are,
and from :mod:`flowfreq.regime`.
"""

from __future__ import annotations

from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

from .blocks import Block, Figure, Heading, Paragraph, Table
from .daily import _STATION_OPTIONS, _daily_frame, _label, _record_table
from .sections import SECTION_KINDS, SectionContent, SectionKind, _figure_png, _version, cfs
from .spec import ReportError

#: Caption of the ``monthly_statistics`` section's main table.
MONTHLY_CAPTION = "Statistics of the monthly mean flows."
#: Caption of the ``summary_hydrograph`` section's table.
SUMMARY_HYDROGRAPH_CAPTION = "Daily mean flow by date, for the 1st and 15th of each month."

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
    if all(means[i].min() > 0 for i in present):
        ax.set_yscale("log")
    ax.set_ylabel("Monthly mean flow (cfs)")
    if label:
        ax.set_title(label)
    ax.grid(True, axis="y", which="both", alpha=0.3)
    return _figure_png(fig)


def _annual_monthly_table(monthly: pd.DataFrame) -> Table:
    """Monthly mean flows, one row to a water year, a dash where the month is incomplete."""
    cells = {
        (int(row["year"]), int(row["month"])): row["mean_flow_cfs"]
        for _, row in monthly[monthly["complete"]].iterrows()
    }
    return Table(
        ("Water year", *(name[:3] for _, name in _MONTHS)),
        tuple(
            (f"{year}", *(cfs(cells.get((year, month), np.nan)) for month, _ in _MONTHS))
            for year in sorted(int(y) for y in monthly["year"].unique())
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
        method.append(f"{partial} month(s) with too few days of record were left out.")
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
            rows.append((f"{name} {day}", f"{int(counts.get(number, 0))}", *cells))

    method = [
        "For each date of the year, the daily mean flows of every year with a value on that "
        "date were ranked and their percentiles taken by linear interpolation. They are "
        "percentiles of non-exceedance: the 10th percentile is the flow that 10 percent of the "
        "years were at or below on that date, which a flow-duration table would call q90. A "
        "curve joins one date's statistic to the next and is not the flow of any one year.",
        f"The computation was made with flowfreq {_version()}.",
    ]
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
)
SECTION_KINDS.update({kind.name: kind for kind in _KINDS})

__all__ = ["MONTHLY_CAPTION", "SUMMARY_HYDROGRAPH_CAPTION"]
