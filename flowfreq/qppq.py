"""
QPPQ streamflow transfer for snowmelt-dominated basins.

QPPQ estimates a daily series at an ungaged target from a gaged donor by
going through both sites' flow-duration curves: donor flow -> donor
exceedance probability -> target flow at that same probability. The donor
supplies the timing, the target's own curve supplies the magnitude.

    Q_donor(t) --FDC_donor--> p(t) --FDC_target--> Q_target(t)

Why not a ratio: a drainage-area ratio applies one multiplier every day,
which assumes the two basins scale identically at every flow level. They do
not -- high flows scale roughly with contributing area, low flows with
storage and geology. Routing through two curves lets the effective
multiplier vary with flow level on its own.

Melt-timing phase error, and what actually fixes it
---------------------------------------------------
QPPQ's load-bearing assumption is that donor and target sit at the *same*
position in their own duration curves on the *same* day. In a snowmelt
basin that fails in a specific, systematic way: a higher, colder target
melts later than a lower donor, so on a rising-limb day the donor may be at
p = 0.10 while the target is still at p = 0.40. That is a phase error, not
noise, and it biases the timing of the freshet.

**Seasonal curves are not the fix, and a coarse season split makes it
worse.** Measured on a synthetic pair with an imposed 25-day melt offset:

===================  ========  ===========
grouping             log-NSE   dry-end bias
===================  ========  ===========
annual (one curve)     0.273        +382%
three seasons          0.040        +369%
monthly                0.605        +149%
25-day donor lag       0.895         +52%
===================  ========  ===========

A three-season split is worse than no split at all, because a 25-day offset
happens *inside* a four-month melt season: both sites are labelled
"freshet", the within-season curve is narrower than the annual one, and the
mismatch is amplified rather than absorbed. Monthly grouping helps because
its blocks are shorter than the offset.

What actually works is aligning the timing: :func:`estimate_donor_lag`
recovers the offset from rank correlation (it found 24 days for an imposed
25), and transferring the lagged donor series restores almost exactly the
no-offset skill. For a target with no record, estimate the lag from
:func:`center_of_timing` instead -- it regresses on basin elevation, so a
target's melt timing is predictable from hypsometry when its flows are not.

Use monthly grouping *and* a lag when both are available; the lag is doing
most of the work.

What this module does not do
----------------------------
- **It cannot reproduce a losing reach.** Where a channel exchanges water
  with an alluvial aquifer and goes subsurface, nothing here knows about
  it: the target's curve is estimated from donor or regional behaviour that
  has never seen that reach. Losing reaches need seepage runs and a reach
  water balance applied *afterwards*. A QPPQ series through one will
  overstate summer flow, confidently.
- **It does not naturalize.** A regulated donor transfers its regulation.
  Choose unregulated donors or naturalize first, and be consistent.
- **Donor ranking here is not the map-correlation method.**
  :func:`rank_donors` needs *some* record at the target. Selecting a donor
  for a truly ungaged site by kriging the correlation field (Archfield &
  Vogel, 2010, WRR 46, W10513) is not implemented.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Mapping, Optional, Sequence, Tuple, Union, cast

import numpy as np
import pandas as pd
from scipy.special import ndtri

logger = logging.getLogger(__name__)

#: Anything that can be read as an array of flows or probabilities.
ArrayLike = Union[Sequence[float], np.ndarray]

#: A regime-based season grouping for a snowmelt basin: accumulation,
#: freshet, recession. Deliberately not ``flowfreq.regime``'s meteorological
#: seasons (DJF/MAM/JJA/SON), which split a spring freshet across two
#: seasons and mix the melt tail with baseflow in a third.
#:
#: **Measured caveat:** for correcting melt-timing phase error between donor
#: and target this grouping is *worse than no grouping at all* -- its blocks
#: are longer than the offset being corrected. See the module docstring's
#: table. Use :data:`MONTHLY_SEASONS` plus :func:`estimate_donor_lag`
#: instead. This grouping remains useful for describing regime differences,
#: which is a different question from correcting timing.
SNOWMELT_SEASONS: Dict[str, Tuple[int, ...]] = {
    "accumulation": (11, 12, 1, 2, 3),
    "freshet": (4, 5, 6, 7),
    "recession": (8, 9, 10),
}

#: Below this, a flow is treated as zero rather than logged.
_ZERO_FLOW_CFS: float = 1e-9


def _season_of(dates: pd.DatetimeIndex, seasons: Mapping[str, Sequence[int]]) -> np.ndarray:
    """Label each date with its season, validating the grouping covers 1-12 once."""
    seen: Dict[int, str] = {}
    for name, months in seasons.items():
        for month in months:
            if month not in range(1, 13):
                raise ValueError(f"season {name!r} names month {month}, which is not 1-12")
            if month in seen:
                raise ValueError(f"month {month} is in both {seen[month]!r} and {name!r}")
            seen[month] = name
    missing = sorted(set(range(1, 13)) - set(seen))
    if missing:
        raise ValueError(f"seasons do not cover month(s) {missing}")

    lookup = np.array([""] * 13, dtype=object)
    for month, name in seen.items():
        lookup[month] = name
    return lookup[dates.month.to_numpy()]


@dataclass
class FlowDurationCurve:
    """An empirical flow-duration curve that can be read in both directions.

    Built from a daily record with Weibull plotting positions
    (``i / (n + 1)``, the convention Vogel & Fennessey use for FDCs).
    Interpolation is done in ``(z, log10 Q)`` space, where
    ``z = ndtri(1 - p)`` -- an FDC is close to a straight line there and
    violently curved in linear space, so interpolating linearly on the raw
    axes would cut corners off the tails.

    Attributes
    ----------
    flows : np.ndarray
        Positive flows, ascending.
    exceedance : np.ndarray
        Exceedance probability of each, descending to match.
    n : int
        Days in the record the curve was built from, zeros included.
    n_zero : int
        Days of zero flow.
    p_zero : float
        ``n_zero / n``. Zero-flow days occupy the wettest-probability end of
        the curve by definition: they are the smallest flows, so they sit at
        exceedance probabilities above ``1 - p_zero``.
    label : str
        For reports and error messages.
    """

    flows: np.ndarray
    exceedance: np.ndarray
    n: int
    n_zero: int
    p_zero: float
    label: str = ""

    @classmethod
    def from_daily(
        cls,
        daily_data: pd.DataFrame,
        label: str = "",
        flow_column: str = "flow_cfs",
    ) -> "FlowDurationCurve":
        """Build a curve from a daily record.

        Parameters
        ----------
        daily_data : pd.DataFrame
            Daily flows. NaNs are dropped; negatives are rejected rather
            than silently clamped, since a negative discharge is a data
            problem, not a small flow.
        label : str, optional
            Name for reports.
        flow_column : str, optional

        Returns
        -------
        FlowDurationCurve

        Raises
        ------
        KeyError, ValueError
        """
        if flow_column not in daily_data.columns:
            raise KeyError(f"daily_data has no {flow_column!r} column")

        values = daily_data[flow_column].dropna().to_numpy(dtype=float)
        if values.size == 0:
            raise ValueError(f"no non-null flows for {label or 'this record'}")
        if np.any(values < 0):
            raise ValueError(
                f"{label or 'this record'} contains negative discharge; that is a data "
                "problem rather than a low flow, and clamping it would hide it"
            )

        n = int(values.size)
        positive = np.sort(values[values > _ZERO_FLOW_CFS])
        n_zero = n - int(positive.size)

        if positive.size == 0:
            raise ValueError(f"{label or 'this record'} is entirely zero flow")

        # Ascending flow, so rank 1 is the smallest. Exceedance of the k-th
        # smallest of n total values is (n - k + 1) / (n + 1).
        ranks = np.arange(1, positive.size + 1, dtype=float) + n_zero
        exceedance = (n - ranks + 1.0) / (n + 1.0)

        return cls(
            flows=positive,
            exceedance=exceedance,
            n=n,
            n_zero=n_zero,
            p_zero=n_zero / n,
            label=label,
        )

    @property
    def is_intermittent(self) -> bool:
        """Whether the record contains any zero-flow days."""
        return self.n_zero > 0

    def probability_of(self, flow: ArrayLike) -> np.ndarray:
        """Exceedance probability of each flow, read off this curve.

        A flow outside the record's range is clamped to the endpoint
        probability rather than extrapolated: an empirical curve says
        nothing beyond the record it was built from, and inventing a tail
        here would propagate straight into the transferred series.

        Zero flows map to the midpoint of the zero block. Their probability
        is genuinely unidentifiable within that block -- every zero-flow day
        is tied -- so the midpoint is a stated convention, not a measurement.
        """
        values = np.asarray(flow, dtype=float).ravel()
        result = np.empty(values.shape, dtype=float)

        zero = values <= _ZERO_FLOW_CFS
        if np.any(zero):
            if self.p_zero > 0:
                result[zero] = 1.0 - self.p_zero / 2.0
            else:
                # The donor never went dry, so it has no zero block to place
                # this in. The driest thing this curve knows is its own
                # minimum.
                result[zero] = float(self.exceedance[0])

        if np.any(~zero):
            z_curve = ndtri(1.0 - self.exceedance)  # ascending with flow
            log_curve = np.log10(self.flows)
            log_query = np.log10(values[~zero])
            # np.interp needs ascending x; log_curve ascends with flow.
            result[~zero] = 1.0 - _ndtr(np.interp(log_query, log_curve, z_curve))

        return result

    def flow_at(self, prob: ArrayLike) -> np.ndarray:
        """Flow at each exceedance probability, read off this curve.

        Probabilities inside the zero block (above ``1 - p_zero``) return
        zero, which is the honest answer for an intermittent record.
        """
        probs = np.asarray(prob, dtype=float).ravel()
        if np.any((probs < 0.0) | (probs > 1.0)):
            raise ValueError("every exceedance probability must lie in [0, 1]")

        result = np.zeros(probs.shape, dtype=float)
        in_zero_block = (
            probs > (1.0 - self.p_zero) if self.p_zero > 0 else np.zeros_like(probs, dtype=bool)
        )

        live = ~in_zero_block
        if np.any(live):
            z_curve = ndtri(1.0 - self.exceedance)
            log_curve = np.log10(self.flows)
            z_query = ndtri(1.0 - np.clip(probs[live], 1e-12, 1.0 - 1e-12))
            result[live] = 10.0 ** np.interp(z_query, z_curve, log_curve)

        return result

    def to_frame(self, exceedance_pct: Optional[ArrayLike] = None) -> pd.DataFrame:
        """Tabulate the curve at given percent-exceeded points."""
        if exceedance_pct is None:
            exceedance_pct = (1, 5, 10, 20, 50, 80, 90, 95, 99)
        pct = np.asarray(exceedance_pct, dtype=float)
        return pd.DataFrame(
            {
                "exceedance_pct": pct,
                "exceedance_prob": pct / 100.0,
                "flow_cfs": self.flow_at(pct / 100.0),
            }
        )


#: A single curve, or one curve per season.
CurveOrCurves = Union[FlowDurationCurve, Dict[str, FlowDurationCurve]]


def _ndtr(z: np.ndarray) -> np.ndarray:
    """Standard normal CDF (scipy.special.ndtr, imported lazily by name)."""
    from scipy.special import ndtr

    return np.asarray(ndtr(z), dtype=float)


def seasonal_curves(
    daily_data: pd.DataFrame,
    seasons: Optional[Mapping[str, Sequence[int]]] = None,
    label: str = "",
    flow_column: str = "flow_cfs",
) -> Dict[str, FlowDurationCurve]:
    """One :class:`FlowDurationCurve` per season.

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily flows, DatetimeIndex.
    seasons : mapping, optional
        ``{season_name: (months,)}``. Defaults to :data:`SNOWMELT_SEASONS`.
        Must cover months 1-12 exactly once.
    label, flow_column : str, optional

    Returns
    -------
    dict of str to FlowDurationCurve

    Raises
    ------
    ValueError
        If the grouping does not partition the year, or a season has no data.
    """
    grouping = SNOWMELT_SEASONS if seasons is None else seasons
    if not isinstance(daily_data.index, pd.DatetimeIndex):
        raise TypeError("daily_data must have a DatetimeIndex to be split by season")

    season_labels = _season_of(daily_data.index, grouping)
    curves: Dict[str, FlowDurationCurve] = {}
    for name in grouping:
        subset = daily_data.loc[season_labels == name]
        if subset.empty:
            raise ValueError(f"season {name!r} has no data in this record")
        curves[name] = FlowDurationCurve.from_daily(
            subset, label=f"{label} {name}".strip(), flow_column=flow_column
        )
    return curves


@dataclass
class QppqResult:
    """A transferred daily series and what it rests on.

    Attributes
    ----------
    series : pd.DataFrame
        Indexed by date: ``donor_flow_cfs``, ``exceedance_prob``,
        ``flow_cfs`` (the estimate), ``season``.
    n_clamped_high, n_clamped_low : int
        Donor days more extreme than anything in the donor's own record, so
        their probability was clamped to the record's endpoint rather than
        extrapolated. A large count means the donor record is too short to
        describe the days being transferred.
    seasonal : bool
    donor_label, target_label : str
    """

    series: pd.DataFrame
    n_clamped_high: int = 0
    n_clamped_low: int = 0
    seasonal: bool = False
    donor_label: str = ""
    target_label: str = ""
    notes: Tuple[str, ...] = field(default_factory=tuple)


def qppq(
    donor_daily: pd.DataFrame,
    donor_curve: "CurveOrCurves",
    target_curve: "CurveOrCurves",
    *,
    seasons: Optional[Mapping[str, Sequence[int]]] = None,
    flow_column: str = "flow_cfs",
    donor_label: str = "",
    target_label: str = "",
) -> QppqResult:
    """Transfer a donor's daily series to a target through both FDCs.

    Parameters
    ----------
    donor_daily : pd.DataFrame
        The donor's daily record, DatetimeIndex.
    donor_curve, target_curve : FlowDurationCurve or dict of str to FlowDurationCurve
        Pass single curves for an annual transfer, or matching dicts keyed
        by season for a seasonal one. **Seasonal is strongly preferred in a
        snowmelt basin** -- see the module docstring on phase error.
    seasons : mapping, optional
        Required (or defaulted) when curves are dicts. Defaults to
        :data:`SNOWMELT_SEASONS`.
    flow_column, donor_label, target_label : str, optional

    Returns
    -------
    QppqResult

    Raises
    ------
    ValueError, TypeError
    """
    if flow_column not in donor_daily.columns:
        raise KeyError(f"donor_daily has no {flow_column!r} column")
    if not isinstance(donor_daily.index, pd.DatetimeIndex):
        raise TypeError("donor_daily must have a DatetimeIndex")

    seasonal = isinstance(donor_curve, dict)
    if seasonal != isinstance(target_curve, dict):
        raise ValueError(
            "donor_curve and target_curve must both be single curves or both be "
            "dicts of seasonal curves; mixing them would transfer a season's "
            "probability through an annual curve"
        )

    frame = donor_daily[[flow_column]].dropna().copy()
    if frame.empty:
        raise ValueError("donor_daily has no non-null flows to transfer")
    donor_flows = frame[flow_column].to_numpy(dtype=float)

    probs = np.empty(donor_flows.shape, dtype=float)
    estimates = np.empty(donor_flows.shape, dtype=float)

    if isinstance(donor_curve, dict) and isinstance(target_curve, dict):
        donor_map: Dict[str, FlowDurationCurve] = donor_curve
        target_map: Dict[str, FlowDurationCurve] = target_curve
        grouping = SNOWMELT_SEASONS if seasons is None else seasons
        labels = _season_of(donor_daily.index, grouping)
        absent = set(np.unique(labels)) - set(donor_map)
        if absent:
            raise ValueError(f"no donor curve for season(s) {sorted(absent)}")
        absent = set(np.unique(labels)) - set(target_map)
        if absent:
            raise ValueError(f"no target curve for season(s) {sorted(absent)}")

        for name in np.unique(labels):
            mask = labels == name
            probs[mask] = donor_map[name].probability_of(donor_flows[mask])
            estimates[mask] = target_map[name].flow_at(probs[mask])
        season_column = labels
        reference_donor = next(iter(donor_map.values()))
    else:
        assert isinstance(donor_curve, FlowDurationCurve)
        assert isinstance(target_curve, FlowDurationCurve)
        probs = donor_curve.probability_of(donor_flows)
        estimates = target_curve.flow_at(probs)
        season_column = np.array(["annual"] * len(frame), dtype=object)
        reference_donor = donor_curve

    # Days outside the donor record's own range had their probability
    # clamped rather than extrapolated; count them so a short donor record
    # is visible rather than implied.
    if isinstance(donor_curve, dict):
        n_high = 0
        n_low = 0
        for name in np.unique(season_column):
            mask = season_column == name
            curve = donor_curve[name]
            n_high += int(np.count_nonzero(donor_flows[mask] > curve.flows[-1]))
            n_low += int(
                np.count_nonzero(
                    (donor_flows[mask] < curve.flows[0]) & (donor_flows[mask] > _ZERO_FLOW_CFS)
                )
            )
    else:
        n_high = int(np.count_nonzero(donor_flows > reference_donor.flows[-1]))
        n_low = int(
            np.count_nonzero(
                (donor_flows < reference_donor.flows[0]) & (donor_flows > _ZERO_FLOW_CFS)
            )
        )

    if n_high or n_low:
        logger.warning(
            "%d donor day(s) fell outside the donor record's own range and were clamped "
            "to its endpoints rather than extrapolated (%d above, %d below).",
            n_high + n_low,
            n_high,
            n_low,
        )

    notes = []
    if not seasonal:
        notes.append(
            "Annual transfer. In a snowmelt basin a higher target melts later than a "
            "lower donor, so same-day probability equivalence carries a seasonal phase "
            "error; seasonal curves absorb most of it."
        )

    series = pd.DataFrame(
        {
            "donor_flow_cfs": donor_flows,
            "exceedance_prob": probs,
            "flow_cfs": estimates,
            "season": season_column,
        },
        index=frame.index,
    )

    return QppqResult(
        series=series,
        n_clamped_high=n_high,
        n_clamped_low=n_low,
        seasonal=seasonal,
        donor_label=donor_label,
        target_label=target_label,
        notes=tuple(notes),
    )


# =============================================================================
# Donor selection and validation
# =============================================================================


def rank_donors(
    target_daily: pd.DataFrame,
    candidates: Mapping[str, pd.DataFrame],
    *,
    min_concurrent_days: int = 365,
    flow_column: str = "flow_cfs",
) -> pd.DataFrame:
    """Rank candidate donors by concurrent-day agreement with the target.

    Correlation is computed on ``log10`` flows, because an untransformed
    correlation is dominated by the freshet and says almost nothing about
    whether the two records agree in the low-flow window that usually
    limits habitat. Spearman is reported alongside because QPPQ only ever
    uses the donor's *rank* -- a donor with a perfect rank relationship and
    a poor linear one is a good donor for this method.

    **This is not the map-correlation method.** It needs a record at the
    target. Selecting a donor for a fully ungaged site means kriging the
    correlation field (Archfield & Vogel, 2010); that is not implemented,
    and this function will not stand in for it.

    Parameters
    ----------
    target_daily : pd.DataFrame
        The target's record, however short.
    candidates : mapping of str to pd.DataFrame
        Candidate donor records.
    min_concurrent_days : int, optional
        Candidates with less overlap are returned with NaN scores rather
        than dropped, so a thin overlap is visible instead of invisible.
    flow_column : str, optional

    Returns
    -------
    pd.DataFrame
        One row per candidate, sorted best-first by Spearman: ``donor``,
        ``n_concurrent``, ``pearson_log``, ``spearman``, ``sufficient``.
    """
    if flow_column not in target_daily.columns:
        raise KeyError(f"target_daily has no {flow_column!r} column")

    target = target_daily[flow_column].dropna()
    rows = []
    for name, frame in candidates.items():
        if flow_column not in frame.columns:
            raise KeyError(f"candidate {name!r} has no {flow_column!r} column")
        joined = pd.concat(
            [target.rename("target"), frame[flow_column].dropna().rename("donor")],
            axis=1,
            join="inner",
        ).dropna()
        both_positive = joined[(joined["target"] > 0) & (joined["donor"] > 0)]
        n = int(len(both_positive))
        if n >= 2:
            log_target = np.log10(both_positive["target"])
            log_donor = np.log10(both_positive["donor"])
            pearson = float(log_target.corr(log_donor))
            spearman = float(
                both_positive["target"].corr(both_positive["donor"], method="spearman")
            )
        else:
            pearson = float("nan")
            spearman = float("nan")
        rows.append(
            {
                "donor": name,
                "n_concurrent": n,
                "pearson_log": pearson,
                "spearman": spearman,
                "sufficient": n >= min_concurrent_days,
            }
        )

    return (
        pd.DataFrame(rows)
        .sort_values(["sufficient", "spearman"], ascending=[False, False])
        .reset_index(drop=True)
    )


def performance(observed: ArrayLike, estimated: ArrayLike) -> Dict[str, float]:
    """Goodness-of-fit for a transferred series.

    Reports both NSE and log-NSE deliberately. NSE is dominated by the
    freshet and will look respectable even when summer flows are wrong by a
    factor of two; log-NSE and the dry-end bias are what speak to the
    low-flow window that limits habitat. Judging a Methow transfer on NSE
    alone would pass a series that is useless for the question being asked.

    Returns
    -------
    dict
        ``nse``, ``log_nse``, ``kge``, ``pbias`` (percent), and
        ``pbias_dry`` (percent bias over the driest 10% of observed days).
    """
    obs = np.asarray(observed, dtype=float).ravel()
    est = np.asarray(estimated, dtype=float).ravel()
    if obs.shape != est.shape:
        raise ValueError(f"shape mismatch: {obs.shape} vs {est.shape}")

    valid = np.isfinite(obs) & np.isfinite(est)
    obs, est = obs[valid], est[valid]
    if obs.size < 2:
        raise ValueError("need at least two paired values")

    def _nse(a: np.ndarray, b: np.ndarray) -> float:
        denominator = float(np.sum((a - a.mean()) ** 2))
        if denominator == 0:
            return float("nan")
        return float(1.0 - np.sum((a - b) ** 2) / denominator)

    positive = (obs > 0) & (est > 0)
    log_nse = (
        _nse(np.log10(obs[positive]), np.log10(est[positive]))
        if positive.sum() >= 2
        else float("nan")
    )

    obs_std, est_std = float(obs.std()), float(est.std())
    obs_mean, est_mean = float(obs.mean()), float(est.mean())
    if obs_std == 0 or obs_mean == 0:
        kge = float("nan")
    else:
        r = float(np.corrcoef(obs, est)[0, 1])
        kge = float(
            1.0
            - np.sqrt(
                (r - 1.0) ** 2 + (est_std / obs_std - 1.0) ** 2 + (est_mean / obs_mean - 1.0) ** 2
            )
        )

    pbias = float(100.0 * (est.sum() - obs.sum()) / obs.sum()) if obs.sum() else float("nan")

    dry_cut = float(np.quantile(obs, 0.10))
    dry = obs <= dry_cut
    pbias_dry = (
        float(100.0 * (est[dry].sum() - obs[dry].sum()) / obs[dry].sum())
        if dry.any() and obs[dry].sum() > 0
        else float("nan")
    )

    return {
        "nse": _nse(obs, est),
        "log_nse": log_nse,
        "kge": kge,
        "pbias": pbias,
        "pbias_dry": pbias_dry,
    }


def loocv_qppq(
    sites: Mapping[str, pd.DataFrame],
    *,
    seasons: Optional[Mapping[str, Sequence[int]]] = None,
    seasonal: bool = True,
    flow_column: str = "flow_cfs",
    donor_chooser: Optional[Callable[[str, pd.DataFrame, Dict[str, pd.DataFrame]], str]] = None,
) -> pd.DataFrame:
    """Leave-one-out check of the probability-transfer assumption.

    Each gaged site is treated in turn as if it were ungaged: a donor is
    chosen from the others, and QPPQ transfers that donor's series through
    the *withheld site's own observed* FDC. That isolates the assumption
    this method actually rests on -- that donor and target occupy the same
    position in their own curves on the same day -- from the separate
    question of how well the target's curve can be estimated when it truly
    is ungaged.

    Read it that way. A poor score here means probability transfer itself
    does not hold between these basins (phase error, differing storm
    coverage, regulation), and no improvement in FDC regionalization will
    rescue it. A good score here is necessary, not sufficient: the full
    ungaged chain also has to estimate the target curve.

    **This is the gate for whether added complexity earns its keep.** Run it
    seasonal and annual; if seasonal does not beat annual on ``log_nse`` and
    ``pbias_dry``, the seasonal machinery is not paying for itself here.

    Parameters
    ----------
    sites : mapping of str to pd.DataFrame
        At least two daily records with DatetimeIndexes.
    seasons : mapping, optional
    seasonal : bool, optional
        Whether to build and transfer through seasonal curves.
    flow_column : str, optional
    donor_chooser : callable, optional
        ``(withheld_name, withheld_daily, others) -> donor_name``. Defaults
        to the best Spearman rank from :func:`rank_donors`.

    Returns
    -------
    pd.DataFrame
        One row per site: ``site``, ``donor``, ``n_days``, ``seasonal``,
        and the metrics from :func:`performance`.
    """
    if len(sites) < 2:
        raise ValueError("leave-one-out needs at least two sites")

    def _default_chooser(
        _name: str, withheld: pd.DataFrame, others: Dict[str, pd.DataFrame]
    ) -> str:
        ranked = rank_donors(withheld, others, flow_column=flow_column)
        return str(ranked.iloc[0]["donor"])

    chooser = donor_chooser or _default_chooser
    rows = []

    for name, withheld in sites.items():
        others = {other: frame for other, frame in sites.items() if other != name}
        donor_name = chooser(name, withheld, others)
        donor_daily = sites[donor_name]

        if seasonal:
            donor_curve: CurveOrCurves = seasonal_curves(
                donor_daily, seasons, label=donor_name, flow_column=flow_column
            )
            target_curve: CurveOrCurves = seasonal_curves(
                withheld, seasons, label=name, flow_column=flow_column
            )
        else:
            donor_curve = FlowDurationCurve.from_daily(
                donor_daily, label=donor_name, flow_column=flow_column
            )
            target_curve = FlowDurationCurve.from_daily(
                withheld, label=name, flow_column=flow_column
            )

        transferred = qppq(
            donor_daily,
            donor_curve,
            target_curve,
            seasons=seasons,
            flow_column=flow_column,
            donor_label=donor_name,
            target_label=name,
        )

        paired = pd.concat(
            [
                withheld[flow_column].dropna().rename("observed"),
                transferred.series["flow_cfs"].rename("estimated"),
            ],
            axis=1,
            join="inner",
        ).dropna()

        if len(paired) < 2:
            logger.warning(
                "Site %s and donor %s share fewer than two days; no score.", name, donor_name
            )
            continue

        metrics = performance(
            paired["observed"].to_numpy(dtype=float), paired["estimated"].to_numpy(dtype=float)
        )
        rows.append(
            {
                "site": name,
                "donor": donor_name,
                "n_days": int(len(paired)),
                "seasonal": seasonal,
                **metrics,
            }
        )

    return pd.DataFrame(rows)


#: One curve per calendar month. Measured to work substantially better than
#: a coarse three-season split when donor and target melt at different
#: times -- see :data:`SNOWMELT_SEASONS`' own caveat and
#: ``docs/METHOW_METHODS.md``.
MONTHLY_SEASONS: Dict[str, Tuple[int, ...]] = {f"m{m:02d}": (m,) for m in range(1, 13)}


def center_of_timing(
    daily_data: pd.DataFrame,
    flow_column: str = "flow_cfs",
    min_days: int = 350,
) -> pd.DataFrame:
    """Day of the water year by which half the year's volume has passed.

    The standard summary of *when* a snowmelt basin delivers its water. It
    is the physically meaningful way to compare melt timing between two
    basins, and therefore the way to estimate the donor-to-target lag
    (:func:`estimate_donor_lag`) for a site with no flow record of its own:
    center of timing regresses on mean basin elevation, so a target's
    timing can be predicted from its hypsometry even when its flows cannot.

    Parameters
    ----------
    daily_data : pd.DataFrame
        Daily flows with a DatetimeIndex.
    flow_column : str, optional
    min_days : int, optional
        Water years with fewer valid days are reported with NaN rather than
        computed from a partial record.

    Returns
    -------
    pd.DataFrame
        Columns ``water_year``, ``n_days``, ``center_of_timing_dowy``
        (1 = Oct 1), ``complete``.
    """
    from .core import assign_year_label

    if flow_column not in daily_data.columns:
        raise KeyError(f"daily_data has no {flow_column!r} column")
    if not isinstance(daily_data.index, pd.DatetimeIndex):
        raise TypeError("daily_data must have a DatetimeIndex")

    frame = daily_data[[flow_column]].copy()
    frame["water_year"] = assign_year_label(daily_data.index, "water")
    frame = frame.dropna(subset=[flow_column])

    rows = []
    for year, group in frame.groupby("water_year"):
        # pandas-stubs types a groupby key over a mixed-dtype-capable union
        # (str/date/complex/...); this column is always int (assign_year_label's
        # contract), so int() on it is genuinely safe -- the cast documents that
        # to mypy without disabling checking on the rest of the line.
        water_year = int(cast(Any, year))
        flows = group[flow_column].to_numpy(dtype=float)
        n_days = int(flows.size)
        total = float(flows.sum())
        if n_days < min_days or total <= 0:
            rows.append(
                {
                    "water_year": water_year,
                    "n_days": n_days,
                    "center_of_timing_dowy": float("nan"),
                    "complete": False,
                }
            )
            continue
        cumulative = np.cumsum(flows)
        # +1 so the first day of the water year is day 1, not day 0.
        dowy = int(np.searchsorted(cumulative, 0.5 * total) + 1)
        rows.append(
            {
                "water_year": water_year,
                "n_days": n_days,
                "center_of_timing_dowy": float(dowy),
                "complete": True,
            }
        )

    return pd.DataFrame(rows)


def estimate_donor_lag(
    target_daily: pd.DataFrame,
    donor_daily: pd.DataFrame,
    *,
    max_lag_days: int = 45,
    flow_column: str = "flow_cfs",
) -> pd.DataFrame:
    """Find the donor lag that best aligns donor and target timing.

    Scored by **Spearman** rank correlation on concurrent days, because
    QPPQ uses only the donor's rank -- the probability it reads off the
    donor curve. Maximizing rank agreement is therefore maximizing the
    thing the method actually depends on, and it needs no FDC to evaluate.

    A positive lag means the donor's series is shifted *later* to match the
    target, i.e. the target melts later than the donor -- normally because
    it is higher or colder.

    This needs a record at the target. For a genuinely ungaged target,
    estimate the lag instead from the difference in
    :func:`center_of_timing`, predicting the target's from basin elevation.

    Parameters
    ----------
    target_daily, donor_daily : pd.DataFrame
        Daily records with DatetimeIndexes.
    max_lag_days : int, optional
        Largest shift to consider, in each direction.
    flow_column : str, optional

    Returns
    -------
    pd.DataFrame
        ``lag_days``, ``spearman``, ``n_concurrent``, sorted best-first.
        The caller gets the whole profile, not just the argmax, because a
        flat or double-peaked profile means the lag is not identified and
        that should be visible rather than reduced to a single number.
    """
    if max_lag_days < 0:
        raise ValueError("max_lag_days must be non-negative")
    target = target_daily[flow_column].dropna().rename("target")
    rows = []
    for lag in range(-max_lag_days, max_lag_days + 1):
        shifted = donor_daily[flow_column].dropna().rename("donor")
        shifted.index = shifted.index + pd.Timedelta(days=lag)
        joined = pd.concat([target, shifted], axis=1, join="inner").dropna()
        n = int(len(joined))
        spearman = (
            float(joined["target"].corr(joined["donor"], method="spearman"))
            if n >= 2
            else float("nan")
        )
        rows.append({"lag_days": lag, "spearman": spearman, "n_concurrent": n})

    return pd.DataFrame(rows).sort_values("spearman", ascending=False).reset_index(drop=True)


def apply_lag(daily_data: pd.DataFrame, lag_days: int) -> pd.DataFrame:
    """Shift a daily record forward by ``lag_days`` (positive = later)."""
    shifted = daily_data.copy()
    shifted.index = shifted.index + pd.Timedelta(days=int(lag_days))
    return shifted
