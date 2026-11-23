"""Report sections describing a gage's flow regime: ``monthly_statistics``.

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
        most, least = group.loc[means.idxmax()], group.loc[means.idxmin()]
        rows.append(
            (
                name,
                f"{len(group)}",
                cfs(means.mean()),
                cfs(means.median()),
                cfs(means.std(ddof=1)),
                f"{cfs(most['mean_flow_cfs'])} ({int(most['year'])})",
                f"{cfs(least['mean_flow_cfs'])} ({int(least['year'])})",
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
)
SECTION_KINDS.update({kind.name: kind for kind in _KINDS})

__all__ = ["MONTHLY_CAPTION"]
