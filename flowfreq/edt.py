"""
flowfreq.edt - Level 2 hydrology attribute statistics for EDT

The rating statistics and record-sufficiency checks for the hydrology
attributes of Ecosystem Diagnosis and Treatment, as defined by:

    Doyle, E.G., and L.C. Lestelle. 2021. Updated guidelines for rating
    Level 2 environmental attributes in Ecosystem Diagnosis and Treatment
    (EDT). Prepared for the Okanogan Basin Monitoring and Evaluation Program
    and Okanogan Subbasin Habitat Implementation Program, Omak, WA, by
    Confluence, Seattle, Washington.

Sections cited throughout as "D&L §15.6" and so on. This module lives in
flowfreq rather than in a basin repository because EDT is applied across the
Columbia Basin: nothing here is specific to one watershed.

What this module does NOT do
----------------------------
**It does not produce an EDT index rating.** The published ratings are
*continuous* -- D&L Table 3 rates real streams at 2.53, 1.92, 3.37 -- and
§§36.6/37.6 confirm that "index values can be identified as non-integers to
represent the lower or upper ends of a range". That interpolation belongs to
the Flow High rating tool in the EDT Excel plug-in (D&L §15.7), and a second
implementation of it here could silently disagree with the official one on a
number that drives a habitat model.

So :func:`index_band` returns the **integer band** a statistic falls in,
which is nothing more than reading the published Categorical Conclusions
table, and the continuous rating is left to EDT's own tool. A band is a
label; a rating is a number. Do not report a band as a rating.

**It does not build flow patterns.** Every flow attribute is *shaped*: the
analyst rates one month and EDT propagates it across the year with a set of
monthly multipliers (D&L Appendix A). The standard procedure for deriving
those multipliers is specified and implementable, but it is not implemented
here yet.

**The attribute names mislead, and the statistics are what count.**
Attributes 15 and 16 are both named for inter-annual *variability* and
neither is rated on a variability statistic -- 15 is the percent change in
Q2yr and 16 the percent change in the 45- or 60-day lowest average daily
flow, both changes in *magnitude* between two periods. §15.7 specifically
cautions that the coefficient of variation of annual peak flow "requires a
very long times series of data" and moves counter-intuitively, decreasing as
impervious surface rises. Do not substitute a spread statistic for these.

Sign convention
---------------
Every band table here is expressed in **signed percent change**, period 1 to
period 2, where negative means a reduction. The source mixes "decrease",
"increase" and "reduction" between attributes and between index levels of
the same attribute; normalizing to one signed scale is the only way to keep
the ladders comparable and is why this module restates them rather than
quoting them.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .lowflow import annual_minimum_flow
from .regime import tqmean

logger = logging.getLogger(__name__)

#: Averaging windows D&L §16.6 names for the low-flow statistic. The guideline
#: says "the 45 or 60-day consecutive lowest average daily flow" and does not
#: choose between them, so both are computed and reported.
LOW_FLOW_WINDOWS: Tuple[int, ...] = (45, 60)


@dataclass(frozen=True)
class AttributeSpec:
    """One EDT Level 2 hydrology attribute, as D&L defines it.

    Attributes
    ----------
    attribute_id : str
        D&L section number, e.g. ``"15"``.
    name : str
        The attribute name as published, misleading or not.
    statistic : str
        The quantity actually rated. Read this, not `name`.
    shaped : bool
        Whether the attribute is shaped across months from a single rated
        month (D&L Appendix A).
    rating_month : str
        Which month carries the rating.
    min_years_total, min_years_per_state : float, optional
        Record length the statistic needs. ``None`` where the guideline
        states none.
    needs_reference_period : bool
        Whether the rating is a change against a reference state. False for
        attribute 17, whose pristine value is fixed at index 0 a priori
        (§17.7), and for 36, which is the natural regime itself.
    citation : str
        The D&L sections this row is read from.
    """

    attribute_id: str
    name: str
    statistic: str
    shaped: bool
    rating_month: str
    min_years_total: Optional[float]
    min_years_per_state: Optional[float]
    needs_reference_period: bool
    citation: str


#: The hydrology attributes, keyed by D&L section number. Attribute 37 is
#: included and marked retired so that a caller reaching for it finds out why
#: rather than finding nothing.
EDT_HYDROLOGY_ATTRIBUTES: Dict[str, AttributeSpec] = {
    "15": AttributeSpec(
        attribute_id="15",
        name="Flow - inter-annual variability in high flows",
        statistic=(
            "percent change in Q2yr, computed from annual instantaneous peaks "
            "(NOT the highest mean daily flow -- §15.7); or percent change in "
            "average TQmean as a sanctioned alternative"
        ),
        shaped=True,
        rating_month="the month with the highest flow (Appendix A)",
        min_years_total=40.0,
        min_years_per_state=20.0,
        needs_reference_period=True,
        citation="D&L §§15.4, 15.6, 15.7, 15.8; Appendix A",
    ),
    "16": AttributeSpec(
        attribute_id="16",
        name="Flow - change in inter-annual variability in low flows",
        statistic=("percent change in the 45- or 60-day consecutive lowest average " "daily flow"),
        shaped=True,
        rating_month="the month when low flow variability between years is greatest",
        min_years_total=40.0,
        min_years_per_state=20.0,
        needs_reference_period=True,
        citation="D&L §§16.4, 16.6, 16.8",
    ),
    "17": AttributeSpec(
        attribute_id="17",
        name="Flow - intra-daily (diel) variation",
        statistic="monthly mean of the change in stage per hour, in inches",
        shaped=True,
        rating_month="the month when diel variation is greatest",
        min_years_total=None,
        min_years_per_state=None,
        needs_reference_period=False,
        citation="D&L §§17.4, 17.6, 17.7, 17.8",
    ),
    "18": AttributeSpec(
        attribute_id="18",
        name="Flow - intra-annual flow pattern",
        statistic="percent change in average TQmean (flashiness)",
        shaped=True,
        rating_month=(
            "CONTESTED -- §18.8 says the month when intra-annual flow variation "
            "is greatest, assumed to be a high-runoff month; Appendix A says "
            "the month with the lowest flow. Resolve with the rating owner."
        ),
        min_years_total=10.0,
        min_years_per_state=10.0,
        needs_reference_period=True,
        citation="D&L §§18.4, 18.6, 18.7, 18.8; Appendix A",
    ),
    "36": AttributeSpec(
        attribute_id="36",
        name="Hydrologic regime - natural",
        statistic=(
            "categorical class by dominant water source; integer only -- "
            '"only a categorical rating is acceptable for this attribute"'
        ),
        shaped=False,
        rating_month="not shaped",
        min_years_total=None,
        min_years_per_state=None,
        needs_reference_period=False,
        citation="D&L §§36.4, 36.6, 36.7, 36.8",
    ),
    "37": AttributeSpec(
        attribute_id="37",
        name="Hydrologic regime - regulated (RETIRED)",
        statistic=(
            "RETIRED -- no rule curves were developed and it computes no Level 3 "
            "survival factor. Kept only to document legacy projects. Diversions "
            "belong to §55 Water Withdrawals."
        ),
        shaped=False,
        rating_month="not shaped",
        min_years_total=None,
        min_years_per_state=None,
        needs_reference_period=False,
        citation="D&L §37",
    ),
}


@dataclass(frozen=True)
class IndexBand:
    """One index level of a Categorical Conclusions ladder.

    Half-open on its ladder's scale: ``low <= value < high``.

    Attributes
    ----------
    index : int
        The EDT index level, 0 to 4.
    low, high : float
        Band edges. The outermost bands of most ladders are opened to
        infinity so that no value is unclassifiable -- see `printed_edge`.
    label : str
        Condensed from the published Categorical Conclusions cell.
    printed_edge : float, optional
        Where the *source* printed this band's outer edge, when that differs
        from `low`/`high` above. D&L's outermost cells stop at finite values
        (attribute 18's index 4 at 45% reduction, for instance) while the
        real world does not, so the lookup extends them and
        :func:`index_band` warns when a value falls past what was actually
        published. Without this the extension would be invisible, and an
        extrapolated band would read exactly like a quoted one.
    provisional : bool
        True where one of this band's edges is a local convention rather than
        a published threshold. Only the diel ladder's 0/1 boundary is -- see
        :data:`DIEL_NO_VARIATION_IN_PER_HR`. :func:`index_band` warns, so a
        convention never passes for a citation.
    """

    index: int
    low: float
    high: float
    label: str
    printed_edge: Optional[float] = None
    provisional: bool = False


#: The 0/1 boundary of the diel ladder is **not in the source**. D&L §17.6
#: describes index 0 as "essentially no variation in discharge" and index 1 as
#: "<2 inches change in stage per hour", giving no number between them. This
#: is therefore a local convention, not a quoted threshold: a monthly mean
#: below it reads as no diel signal at all. It needs the rating owner's
#: agreement before it decides anything, and :func:`index_band` says so.
DIEL_NO_VARIATION_IN_PER_HR: float = 0.25


#: Attribute 15 via Q2yr (D&L §15.6), on signed percent change. Index 4's
#: upper edge reads ">40% and <110%+" in the source; the trailing "+" is taken
#: as open-ended, so 110 is recorded as the printed edge rather than a cap.
Q2YR_CHANGE_BANDS: Tuple[IndexBand, ...] = (
    IndexBand(0, -100.0, -40.0, "peak flows strongly reduced"),
    IndexBand(1, -40.0, -20.0, "peak flows moderately reduced"),
    IndexBand(2, -20.0, 20.0, "comparable to undisturbed"),
    IndexBand(3, 20.0, 40.0, "peak flows moderately increased"),
    IndexBand(4, 40.0, np.inf, "peak flows strongly increased", printed_edge=110.0),
)

#: Attribute 15 via TQmean (D&L §15.6, the third rating option). Bands are on
#: signed percent change, so a reduction is negative.
#:
#: **This ladder cannot reach index 0 or 1.** Those cells are defined by a
#: decrease in Q2yr or by known regulation, not by TQmean, so the lowest level
#: the TQmean option can assign is 2. A reach whose peaks are genuinely
#: *suppressed* by regulation therefore cannot be rated from TQmean alone --
#: it needs Q2yr or the qualitative option. Index 2 extends upward without
#: limit because its published cell reads "<5% reduction in average TQmean",
#: and an increase in TQmean satisfies that as written.
TQMEAN_HIGH_FLOW_BANDS: Tuple[IndexBand, ...] = (
    IndexBand(2, -5.0, np.inf, "comparable to undisturbed"),
    IndexBand(3, -15.0, -5.0, "moderately increased peaks"),
    IndexBand(4, -np.inf, -15.0, "strongly increased peaks", printed_edge=-45.0),
)

#: Attribute 18 via TQmean (D&L §18.6), on signed percent change.
TQMEAN_FLASHINESS_BANDS: Tuple[IndexBand, ...] = (
    IndexBand(0, 15.0, np.inf, "runoff response greatly slowed"),
    IndexBand(1, 5.0, 15.0, "runoff response moderately slower"),
    IndexBand(2, -5.0, 5.0, "comparable to undisturbed"),
    IndexBand(3, -15.0, -5.0, "runoff response moderately increased"),
    IndexBand(4, -np.inf, -15.0, "runoff response strongly increased", printed_edge=-45.0),
)

#: Attribute 16 (D&L §16.6), on signed percent change. Index 4's lower edge is
#: -100%, which is a physical floor rather than a printed cap: a flow cannot
#: fall by more than all of itself.
LOW_FLOW_CHANGE_BANDS: Tuple[IndexBand, ...] = (
    IndexBand(0, 75.0, np.inf, "low flows strongly increased"),
    IndexBand(1, 20.0, 75.0, "low flows moderately increased"),
    IndexBand(2, -20.0, 20.0, "comparable to undisturbed"),
    IndexBand(3, -50.0, -20.0, "low flows moderately reduced"),
    IndexBand(4, -100.0, -50.0, "low flows severely reduced"),
)

#: Attribute 17 (D&L §17.6) -- the one ladder that is not a percent change.
#: Units are inches of stage change per hour, as a monthly mean, and **index 0
#: is the pristine reference**, not index 2 as in every other flow attribute.
#: See :data:`DIEL_NO_VARIATION_IN_PER_HR` for the 0/1 boundary, which the
#: source does not supply.
DIEL_STAGE_BANDS: Tuple[IndexBand, ...] = (
    IndexBand(
        0,
        0.0,
        DIEL_NO_VARIATION_IN_PER_HR,
        "essentially no variation; pristine",
        provisional=True,
    ),
    IndexBand(1, DIEL_NO_VARIATION_IN_PER_HR, 2.0, "slight to low variation", provisional=True),
    IndexBand(2, 2.0, 6.0, "low to moderate variation"),
    IndexBand(3, 6.0, 12.0, "moderate to high variation"),
    IndexBand(4, 12.0, np.inf, "extreme variation", printed_edge=24.0),
)


@dataclass(frozen=True)
class BandResult:
    """Which index band a statistic falls in, and what to distrust about it.

    A band is **not** an EDT rating, and is not the whole part of one either.
    The integer index sits at the *midpoint* of its band's range (the footnote
    to D&L §§36.6/37.6), so a value low in band 3 rates below 3.0: flooring
    the band gives the wrong integer for six of the eight examples in D&L
    Table 3. Use :func:`q2yr_rating` for a rating on the Q2yr ladder, and
    EDT's own Flow High rating tool for the others.
    """

    index: int
    label: str
    value: float
    warnings: Tuple[str, ...] = ()


def index_band(
    value: float, bands: Sequence[IndexBand], *, ladder_name: str = "ladder"
) -> BandResult:
    """Find the index band a statistic value falls in.

    Parameters
    ----------
    value : float
        The statistic, on that ladder's own scale -- signed percent change for
        every ladder except :data:`DIEL_STAGE_BANDS`, which is inches of stage
        change per hour.
    bands : sequence of IndexBand
        One of the module's ladders.
    ladder_name : str
        Used in messages so a warning names which ladder produced it.

    Returns
    -------
    BandResult
        The matched band, with warnings for anything the caller should not
        take at face value: a value past the source's printed edge, or a band
        whose edge is a local convention.

        The band is a **label**, not a rating, and not the floor of one --
        see :class:`BandResult`.

    Raises
    ------
    ValueError
        `value` is not finite, or no band covers it. Raising beats rounding
        into a neighbouring band: :data:`TQMEAN_HIGH_FLOW_BANDS` has no level
        below 2 at all, so a value it cannot place is a real gap in the rating
        option rather than a number to nudge.

    Examples
    --------
    >>> index_band(-28.0, LOW_FLOW_CHANGE_BANDS).index
    3
    >>> index_band(3.0, TQMEAN_FLASHINESS_BANDS).index
    2
    """
    if value is None or not np.isfinite(value):
        raise ValueError(f"{ladder_name}: cannot band a non-finite value ({value!r})")

    for band in bands:
        if band.low <= value < band.high:
            notes: List[str] = []
            if band.printed_edge is not None and _beyond_printed_edge(value, band):
                notes.append(
                    f"{value:+.1f} is beyond the outer edge D&L actually printed for "
                    f"index {band.index} ({band.printed_edge:+.1f}); the band was "
                    f"extended to classify it, so treat the level as an extrapolation."
                )
            if band.provisional:
                notes.append(
                    f"index {band.index} on this ladder has an edge that is a local "
                    f"convention, not a published threshold -- confirm it with whoever "
                    f"owns the ratings before it decides anything."
                )
            return BandResult(band.index, band.label, float(value), tuple(notes))

    raise ValueError(
        f"{ladder_name}: {value:+.1f} falls in no band. This ladder covers "
        + ", ".join(f"[{b.low:+g}, {b.high:+g})" for b in bands)
        + "."
    )


def _beyond_printed_edge(value: float, band: IndexBand) -> bool:
    """Whether `value` lies outside the edge the source printed for `band`."""
    if band.printed_edge is None:
        return False
    # Which side the printed edge sits on follows the band's open side.
    if band.high == np.inf:
        return value > band.printed_edge
    if band.low == -np.inf:
        return value < band.printed_edge
    return False


Period = Tuple[pd.Timestamp, pd.Timestamp]


def split_periods(daily_data: pd.DataFrame, split: Optional[str] = None) -> Tuple[Period, Period]:
    """Divide a record into an earlier and a later period.

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily flow with a DatetimeIndex.
    split : str, optional
        Date to split on, ``YYYY-MM-DD``. The earlier period ends the day
        before; the later period starts on it. Default ``None`` splits at the
        midpoint of the record's date range. A split that leaves either period
        without data raises -- see Raises.

    Returns
    -------
    tuple
        ``((start1, end1), (start2, end2))``, inclusive bounds.

    Notes
    -----
    **A contiguous halving is a convenience, not the guideline's own
    construction.** D&L Table 3 compares *separated* windows -- the John Day
    example contrasts 1905-1935 against 1990-2019, leaving a 55-year gap -- on
    the reasoning that each window should represent one watershed development
    state rather than an arbitrary half of the record. §15.7 says so directly:
    "actual application to a watershed might involve different divisions of the
    data record to represent various states of the watershed over time."

    So pass explicit periods to :func:`tqmean_change` and
    :func:`low_flow_change` for real work. This default exists to make a first
    look cheap, and it will understate change wherever development happened
    gradually across the middle of the record.

    Raises
    ------
    ValueError
        The frame is empty, `split` falls outside the record, or `split` leaves
        either period with no data at all.

        Note what is *not* an error: a period that is non-empty but very short.
        A record starting days before the split yields a period spanning under a
        year, which slices fine and is useless. :func:`audit_record` reports
        that as ``split_applicable=False`` rather than raising, because it is a
        judgement about sufficiency, not about slicing.
    """
    if daily_data.empty:
        raise ValueError("cannot split an empty record into periods")
    index = pd.DatetimeIndex(daily_data.index)
    first, last = index.min(), index.max()

    if split is None:
        midpoint = first + (last - first) / 2
    else:
        midpoint = pd.Timestamp(split)
        if not (first < midpoint <= last):
            raise ValueError(
                f"split={split!r} is outside the record, which runs "
                f"{first.date()} to {last.date()}"
            )

    period1 = (first, midpoint - pd.Timedelta(days=1))
    period2 = (midpoint, last)

    # Defensive: for a contiguous record a split inside [first, last] cannot
    # empty either side, so this rarely fires. The case that actually bites in a
    # batch is a period that is non-empty but spans under a year -- a record
    # beginning days before the split -- and that is reported by
    # `audit_record`'s `split_applicable`, not raised here, because how short is
    # too short is a question about sufficiency rather than about slicing.
    for label, (start, end) in (("earlier", period1), ("later", period2)):
        if start > end or daily_data.loc[start:end].empty:
            raise ValueError(
                f"split={midpoint.date()} leaves the {label} period with no data; "
                f"this record runs {first.date()} to {last.date()}."
            )

    return period1, period2


def _years_spanned(start: pd.Timestamp, end: pd.Timestamp) -> float:
    """Calendar years between two timestamps, as a float."""
    return float((end - start).days) / 365.25


def audit_record(
    daily_data: pd.DataFrame,
    *,
    split: Optional[str] = None,
    site_no: str = "",
) -> pd.DataFrame:
    """Check a daily record against each EDT statistic's length requirement.

    The audit that decides which attributes are ratable empirically and which
    need a model, and therefore the one to run before any estimation work.
    D&L attaches a record-length condition to each quantitative option:
    Q2yr and the 45/60-day low flow need "~40 yrs or longer with at least 20
    yrs pertaining to a watershed development state" (§§15.6, 16.6), while
    TQmean needs only about ten years per state because it "does not display
    high inter-annual variability" (§18.7).

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily mean flow with a DatetimeIndex and a ``flow_cfs`` column.
    split : str, optional
        Development-state split date, passed to :func:`split_periods`. Without
        one the per-state clause cannot be checked and the result says so.
    site_no : str
        Carried into the output so audits for several gages concatenate.

    Returns
    -------
    pd.DataFrame
        One row per statistic:

        - ``site_no``, ``statistic``, ``attribute_id``
        - ``years_available`` : span of the record, valid days only
        - ``years_period1``, ``years_period2`` : NaN when no `split` is given
        - ``required_total``, ``required_per_state``
        - ``split_applicable`` : False when `split` leaves a period under a
          year, or when no `split` was given. Distinct from ``sufficient``:
          an inapplicable split needs a different date for this gage, where an
          insufficient record needs more data.
        - ``sufficient`` : bool, or False when it cannot be established
        - ``note`` : why, in words

    Notes
    -----
    ``years_available`` is measured from the span between the first and last
    days carrying a **valid** flow, not from the row count, and a record with
    interior gaps therefore reads as longer than the data it actually holds.
    That is the right way round for this test -- the guideline's condition is
    about how much watershed history the record brackets -- but it means a
    sparse forty-year record can pass an audit it should not. Check
    ``n_valid_days`` against ``years_available * 365`` when that matters.
    """
    if "flow_cfs" not in daily_data.columns:
        raise KeyError(
            f"daily_data has no 'flow_cfs' column; columns are {list(daily_data.columns)}"
        )

    valid = daily_data[daily_data["flow_cfs"] >= 0].dropna(subset=["flow_cfs"])
    if valid.empty:
        raise ValueError("daily_data carries no valid (non-negative, non-NaN) flow values")

    index = pd.DatetimeIndex(valid.index)
    years_available = _years_spanned(index.min(), index.max())
    n_valid_days = int(len(valid))

    if split is None:
        years1 = years2 = float("nan")
        split_applicable = False
    else:
        (s1, e1), (s2, e2) = split_periods(valid, split)
        years1, years2 = _years_spanned(s1, e1), _years_spanned(s2, e2)
        # A period under a year means the split sits at or beside one end of
        # this record: it slices without error and answers nothing. Reported
        # separately from "too short", because the remedy differs -- a
        # degenerate period needs a different split date for this gage, not a
        # longer record.
        split_applicable = min(years1, years2) >= 1.0

    rows: List[Dict[str, object]] = []
    for attribute_id, statistic in (
        ("15", "Q2yr change"),
        ("16", "45/60-day low-flow change"),
        ("18", "TQmean change"),
    ):
        spec = EDT_HYDROLOGY_ATTRIBUTES[attribute_id]
        need_total = spec.min_years_total
        need_state = spec.min_years_per_state
        assert need_total is not None and need_state is not None  # all three are set

        total_ok = years_available >= need_total
        if split is not None and not split_applicable:
            rows.append(
                {
                    "site_no": site_no,
                    "attribute_id": attribute_id,
                    "statistic": statistic,
                    "years_available": years_available,
                    "n_valid_days": n_valid_days,
                    "years_period1": years1,
                    "years_period2": years2,
                    "required_total": need_total,
                    "required_per_state": need_state,
                    "split_applicable": False,
                    "sufficient": False,
                    "note": (
                        f"split {split} does not suit this record: it leaves periods of "
                        f"{years1:.1f} and {years2:.1f} yr, one of them under a year. The "
                        f"record runs {index.min().date()} to {index.max().date()}; choose "
                        f"a split inside it, or exclude this gage from the comparison."
                    ),
                }
            )
            continue
        if split is None:
            state_ok = False
            note = (
                f"record spans {years_available:.1f} yr "
                f"({'meets' if total_ok else 'short of'} the {need_total:.0f} yr total); "
                f"the {need_state:.0f} yr-per-state clause is unchecked because no "
                f"development-state split was given, and choosing that date is a "
                f"judgement, not a data field"
            )
        else:
            state_ok = min(years1, years2) >= need_state
            note = (
                f"record spans {years_available:.1f} yr "
                f"({'meets' if total_ok else 'short of'} {need_total:.0f} yr); periods "
                f"{years1:.1f} / {years2:.1f} yr "
                f"({'both meet' if state_ok else 'short of'} {need_state:.0f} yr per state)"
            )

        rows.append(
            {
                "site_no": site_no,
                "attribute_id": attribute_id,
                "statistic": statistic,
                "years_available": years_available,
                "n_valid_days": n_valid_days,
                "years_period1": years1,
                "years_period2": years2,
                "required_total": need_total,
                "required_per_state": need_state,
                "split_applicable": split_applicable,
                "sufficient": bool(total_ok and state_ok),
                "note": note,
            }
        )

    return pd.DataFrame(rows)


@dataclass(frozen=True)
class PeriodChange:
    """A statistic compared between two periods, with its index band.

    Attributes
    ----------
    attribute_id : str
        The D&L section this serves.
    statistic : str
        What was compared.
    period1, period2 : tuple of pd.Timestamp
        Inclusive bounds actually used.
    value1, value2 : float
        The period means. NaN where a period held no complete year.
    pct_change : float
        ``100 * (value2 - value1) / value1``. Positive is an increase. NaN
        where either value is NaN or `value1` is not positive.
    n_years1, n_years2 : int
        Complete years contributing to each mean.
    band : BandResult, optional
        Where `pct_change` lands on the attribute's ladder. ``None`` when the
        change could not be computed, or when the ladder has no band for it.
    warnings : tuple of str
        Everything to carry into a report.
    """

    attribute_id: str
    statistic: str
    period1: Period
    period2: Period
    value1: float
    value2: float
    pct_change: float
    n_years1: int
    n_years2: int
    band: Optional[BandResult] = None
    warnings: Tuple[str, ...] = field(default_factory=tuple)

    def to_row(self) -> Dict[str, object]:
        """Flatten to one dict, for concatenating across gages."""
        return {
            "attribute_id": self.attribute_id,
            "statistic": self.statistic,
            "period1_start": self.period1[0],
            "period1_end": self.period1[1],
            "period2_start": self.period2[0],
            "period2_end": self.period2[1],
            "value1": self.value1,
            "value2": self.value2,
            "pct_change": self.pct_change,
            "n_years1": self.n_years1,
            "n_years2": self.n_years2,
            "index_band": self.band.index if self.band else None,
            "band_label": self.band.label if self.band else "",
            "warnings": " | ".join(self.warnings),
        }


def _period_mean(
    daily_data: pd.DataFrame,
    period: Period,
    per_year: str,
    compute,
) -> Tuple[float, int, int]:
    """Mean of a per-year statistic over the complete years inside `period`.

    The daily record is sliced first and the per-year statistic computed on
    the slice, rather than computing once and filtering year labels. That
    avoids having to map period dates onto water- or climatic-year labels --
    two different mappings in this library -- and it has the right edge
    behaviour for free: a part-year at either end of the period comes back
    marked incomplete and is excluded from the mean.
    """
    sliced = daily_data.loc[period[0] : period[1]]
    if sliced.empty:
        return float("nan"), 0, 0
    table = compute(sliced)
    complete = table[table["complete"]]
    if complete.empty:
        return float("nan"), 0, int(len(table))
    return float(complete[per_year].mean()), int(len(complete)), int(len(table))


def _change(value1: float, value2: float) -> float:
    """Signed percent change from `value1` to `value2`."""
    if not np.isfinite(value1) or not np.isfinite(value2) or value1 <= 0:
        return float("nan")
    return 100.0 * (value2 - value1) / value1


def tqmean_change(
    daily_data: pd.DataFrame,
    *,
    periods: Optional[Tuple[Period, Period]] = None,
    split: Optional[str] = None,
    year_type: str = "water",
    min_days: int = 350,
    ladder: str = "flashiness",
) -> PeriodChange:
    """Percent change in average TQmean between two periods.

    TQmean is the rating metric for attribute 18 and a sanctioned alternative
    for attribute 15 (D&L §§18.4, 15.6), and it is the statistic a short
    record can actually support: §18.7 establishes that it "does not display
    high inter-annual variability" and "can be estimated reliably from a
    relatively short (e.g., ~10 years) stream flow record", against the forty
    years Q2yr needs.

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily mean flow with a DatetimeIndex and a ``flow_cfs`` column.
    periods : tuple of two (start, end) pairs, optional
        The two windows to compare, inclusive. **Prefer this** -- D&L's own
        examples use separated windows chosen to represent watershed states.
    split : str, optional
        Used only when `periods` is None: split the record here instead of at
        its midpoint. See :func:`split_periods` on why a contiguous halving is
        a weak default.
    year_type : str
        Year definition for TQmean. Default ``"water"``, matching
        :func:`flowfreq.regime.tqmean`.
    min_days : int
        Per-year completeness threshold. Default 350.
    ladder : str
        ``"flashiness"`` for attribute 18's ladder (the default), or
        ``"high_flow"`` for attribute 15's TQmean option, which has no band
        for a TQmean increase.

    Returns
    -------
    PeriodChange

    Notes
    -----
    **TQmean is attenuated by basin size** (§18.7): "runoff patterns are
    attenuated as drainage area increases". That does not affect a
    within-gage comparison like this one, where both periods share a basin,
    but it does mean a TQmean *level* is not comparable between gages of
    different area, and any regression of TQmean onto basin characteristics
    has to carry drainage area as a predictor.

    Raises
    ------
    ValueError
        `ladder` is not a known ladder name, or the record cannot be split.
    KeyError
        No ``flow_cfs`` column.
    """
    ladders = {
        "flashiness": (TQMEAN_FLASHINESS_BANDS, "18"),
        "high_flow": (TQMEAN_HIGH_FLOW_BANDS, "15"),
    }
    if ladder not in ladders:
        raise ValueError(f"ladder must be one of {sorted(ladders)}, got {ladder!r}")
    bands, attribute_id = ladders[ladder]

    if "flow_cfs" not in daily_data.columns:
        raise KeyError(
            f"daily_data has no 'flow_cfs' column; columns are {list(daily_data.columns)}"
        )

    p1, p2 = periods if periods is not None else split_periods(daily_data, split)

    def compute(frame: pd.DataFrame) -> pd.DataFrame:
        return tqmean(frame, year_type=year_type, min_days=min_days)

    v1, n1, _ = _period_mean(daily_data, p1, "tqmean", compute)
    v2, n2, _ = _period_mean(daily_data, p2, "tqmean", compute)
    pct = _change(v1, v2)

    warnings = _period_warnings(n1, n2, EDT_HYDROLOGY_ATTRIBUTES[attribute_id])
    warnings.extend(_construction_warnings(periods, split))

    band: Optional[BandResult] = None
    if np.isfinite(pct):
        try:
            band = index_band(pct, bands, ladder_name=f"TQmean/{ladder}")
            warnings.extend(band.warnings)
        except ValueError as exc:
            warnings.append(str(exc))

    return PeriodChange(
        attribute_id=attribute_id,
        statistic=f"average TQmean ({year_type} year)",
        period1=p1,
        period2=p2,
        value1=v1,
        value2=v2,
        pct_change=pct,
        n_years1=n1,
        n_years2=n2,
        band=band,
        warnings=tuple(warnings),
    )


def low_flow_change(
    daily_data: pd.DataFrame,
    *,
    n_day: int = 45,
    periods: Optional[Tuple[Period, Period]] = None,
    split: Optional[str] = None,
    year_type: str = "climatic",
    min_days: int = 350,
) -> PeriodChange:
    """Percent change in the n-day consecutive lowest average daily flow.

    Attribute 16's rating statistic (D&L §16.6), which names "the 45 or 60-day
    consecutive lowest average daily flow" and does not choose between the two
    windows -- see :data:`LOW_FLOW_WINDOWS`.

    Note what this is *not*: despite the attribute being named "change in
    inter-annual variability in low flows", the quantitative criterion is a
    change in low-flow **magnitude**, not in its spread. A coefficient of
    variation or log-space standard deviation is the wrong statistic here.

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily mean flow with a DatetimeIndex and a ``flow_cfs`` column.
    n_day : int
        Averaging window. 45 or 60 per the guideline; any positive window is
        accepted so a caller can run a sensitivity check, but a value outside
        :data:`LOW_FLOW_WINDOWS` is warned about because it is not the rated
        statistic.
    periods, split : optional
        As :func:`tqmean_change`.
    year_type : str
        Default ``"climatic"`` (Apr 1 - Mar 31), matching
        :func:`flowfreq.lowflow.annual_minimum_flow`: the climatic year exists
        precisely so a single low-flow event is not split across two labels.
    min_days : int
        Per-year completeness threshold. Default 350.

    Returns
    -------
    PeriodChange

    Raises
    ------
    ValueError
        `n_day` is not positive, or the record cannot be split.
    KeyError
        No ``flow_cfs`` column.
    """
    if n_day < 1:
        raise ValueError(f"n_day must be >= 1, got {n_day}")
    if "flow_cfs" not in daily_data.columns:
        raise KeyError(
            f"daily_data has no 'flow_cfs' column; columns are {list(daily_data.columns)}"
        )

    p1, p2 = periods if periods is not None else split_periods(daily_data, split)

    def compute(frame: pd.DataFrame) -> pd.DataFrame:
        return annual_minimum_flow(frame, n_day=n_day, year_type=year_type, min_days=min_days)

    v1, n1, _ = _period_mean(daily_data, p1, "flow_cfs", compute)
    v2, n2, _ = _period_mean(daily_data, p2, "flow_cfs", compute)
    pct = _change(v1, v2)

    warnings = _period_warnings(n1, n2, EDT_HYDROLOGY_ATTRIBUTES["16"])
    if n_day not in LOW_FLOW_WINDOWS:
        warnings.append(
            f"n_day={n_day} is not one of D&L §16.6's windows {list(LOW_FLOW_WINDOWS)}, "
            f"so this is a sensitivity check rather than the rated statistic."
        )
    warnings.extend(_construction_warnings(periods, split))

    band: Optional[BandResult] = None
    if np.isfinite(pct):
        try:
            band = index_band(pct, LOW_FLOW_CHANGE_BANDS, ladder_name="low flow")
            warnings.extend(band.warnings)
        except ValueError as exc:
            warnings.append(str(exc))

    return PeriodChange(
        attribute_id="16",
        statistic=f"{n_day}-day lowest average daily flow ({year_type} year)",
        period1=p1,
        period2=p2,
        value1=v1,
        value2=v2,
        pct_change=pct,
        n_years1=n1,
        n_years2=n2,
        band=band,
        warnings=tuple(warnings),
    )


def _construction_warnings(
    periods: Optional[Tuple[Period, Period]], split: Optional[str]
) -> List[str]:
    """Warn about how the comparison windows were arrived at.

    Three cases, and only one of them is a halving. Conflating them was a bug:
    the warning fired on any ``periods=None`` call and told a caller who had
    supplied an explicit development-state date that their periods "were derived
    by halving the record", which was simply untrue.
    """
    if periods is not None:
        return []
    if split is not None:
        return [
            f"periods were derived by splitting at {split}, which is contiguous. D&L "
            f"Table 3 instead compares separated windows chosen to represent watershed "
            f"states, leaving the transition out of both. Pass periods= to match that "
            f"construction; a contiguous split understates change that accrued gradually "
            f"across the split date."
        ]
    return [
        "periods were derived by halving the record, which is not D&L's construction "
        "and carries no development-state meaning at all -- Table 3 compares separated "
        "windows chosen to represent watershed states. Pass split= at minimum, and "
        "periods= for reportable work."
    ]


def _period_warnings(n1: int, n2: int, spec: AttributeSpec) -> List[str]:
    """Record-length warnings shared by the two change statistics."""
    notes: List[str] = []
    if n1 == 0 or n2 == 0:
        notes.append(
            f"a period held no complete year ({n1} and {n2}), so no change could be "
            f"computed -- this is a record-length failure, not a hydrologic result"
        )
        return notes
    need = spec.min_years_per_state
    if need is not None and min(n1, n2) < need:
        notes.append(
            f"{min(n1, n2)} complete years in the shorter period, against the "
            f"~{need:.0f} D&L asks for attribute {spec.attribute_id} "
            f"({spec.citation}); the change is computable but not reportable as a rating"
        )
    return notes


#: Published percent-change / rating pairs from D&L Table 3, the only worked
#: examples in the hydrology sections. Used to calibrate and to test
#: :func:`q2yr_rating`; the stream names are kept so a failure says which row
#: broke.
TABLE3_Q2YR_EXAMPLES: Tuple[Tuple[str, float, float], ...] = (
    ("Klickitat", -2.3, 1.92),
    ("Nooksack", 10.0, 2.33),
    ("John Day", 16.0, 2.53),
    ("Naselle", 22.0, 2.74),
    ("Dungeness", 28.0, 2.92),
    ("Newaukum", 44.0, 3.17),
    ("NF Stillaguamish", 51.0, 3.26),
    ("Mercer", 60.0, 3.37),
)

#: Breakpoints of the attribute 15 rating transform, as ``(percent change,
#: index)`` pairs. Reproduces every row of :data:`TABLE3_Q2YR_EXAMPLES` to
#: within 0.013 index units -- see :func:`q2yr_rating` for how they were
#: obtained and why the index 4 anchor breaks the pattern of the other four.
Q2YR_RATING_ANCHORS: Tuple[Tuple[float, float], ...] = (
    (-70.0, 0.0),
    (-30.0, 1.0),
    (0.0, 2.0),
    (30.0, 3.0),
    (110.0, 4.0),
)

#: Largest index-unit disagreement between :func:`q2yr_rating` and D&L
#: Table 3, over the eight published examples. Asserted by the test suite, so
#: a change to the anchors that degrades the fit fails the build.
Q2YR_RATING_MAX_RESIDUAL: float = 0.014


def q2yr_rating(pct_change: float) -> float:
    """Continuous attribute 15 rating from a percent change in Q2yr.

    **Read this before using it.** EDT ratings are continuous, not the integer
    band a statistic falls in, and the two do not agree: flooring the band from
    :func:`index_band` gives the wrong integer for six of the eight examples in
    D&L Table 3. So a band is a label and this is the rating.

    How the transform was obtained
    ------------------------------
    From two things the source supplies, not by fitting an arbitrary curve.

    First the footnote to D&L §§36.6 and 37.6: "where an index value is
    associated with a range, **the integer value is assumed for modeling to be
    the midpoint**". Applied to §15.6's ladder that puts index 2 at 0% change
    (midpoint of -20..+20), index 3 at +30% (midpoint of +20..+40), index 1 at
    -30% and index 0 at -70%.

    Second, index 4. Its cell reads ">40% and <110%+", and anchoring 4.0 at the
    printed **+110** rather than at the band midpoint of +75 is what the data
    requires: the midpoint choice misses Table 3 by up to 0.30 index units,
    the printed edge by 0.013. So the outermost anchor is the number the
    source printed, and the interior anchors are midpoints.

    Linear interpolation between those five anchors reproduces all eight
    published examples to within :data:`Q2YR_RATING_MAX_RESIDUAL`, which is
    tighter than the two-significant-figure precision Table 3 reports the
    percent changes at.

    What this does not license
    --------------------------
    **Only the Q2yr ladder.** D&L publishes worked examples for no other
    hydrology statistic, so there is nothing to calibrate or check a transform
    for TQmean, the low-flow change or diel stage against. The midpoint rule
    alone is not enough -- index 4 above is the counter-example -- so those
    ladders get :func:`index_band` and nothing more, and their continuous
    rating should come from the Flow High rating tool in the EDT Excel plug-in
    (§15.7). Do not extend this function by analogy.

    Parameters
    ----------
    pct_change : float
        Signed percent change in Q2yr from the reference period to the current
        one. Positive is an increase.

    Returns
    -------
    float
        The rating, clamped to ``[0, 4]``. Values past the outermost anchors
        clamp rather than extrapolate: an index cannot exceed 4, and a
        straight-line extension of the outer segments has nothing behind it.

    Raises
    ------
    ValueError
        `pct_change` is not finite.

    Examples
    --------
    >>> round(q2yr_rating(16.0), 2)
    2.53
    >>> round(q2yr_rating(60.0), 2)
    3.38
    """
    if pct_change is None or not np.isfinite(pct_change):
        raise ValueError(f"q2yr_rating needs a finite percent change, got {pct_change!r}")
    xs = [a[0] for a in Q2YR_RATING_ANCHORS]
    ys = [a[1] for a in Q2YR_RATING_ANCHORS]
    return float(np.clip(np.interp(float(pct_change), xs, ys), 0.0, 4.0))
