"""
flowfreq.fortran_engine - The vendored USGS Fortran (``emafitpr``) as a selectable engine.

Specified in ``docs/FORTRAN_ENGINE_DESIGN.md``. Three pieces, matching that
document's sections 3-4:

``build_emafit_arrays``
    Translates a :class:`~flowfreq.bulletin17c.Bulletin17C` input set (peaks,
    water years, historical peaks, perception thresholds, a low-outlier
    override) into the ``ql/qu/tl/tu/dtype`` arrays ``emafitpr`` takes.
    Independent of :meth:`ExpectedMomentsAlgorithm._build_flow_intervals`
    (``flowfreq.bulletin17c``) by design -- that method does not censor gap
    years against a perception threshold declared over the *systematic*
    range, only the historical range, which is exactly what the 12363000
    gap-year case (design doc section 3) exercises. Follows ``siteQT``
    (``vendor/peakfqr/R/readInputs.R``) instead.

``run_fortran_reference`` / ``run_fortran_ema``
    Build the arrays, call ``emafitpr`` through
    :meth:`flowfreq.validation.reference.ReferenceResult.from_emafit`, and (for
    ``run_fortran_ema``) adapt the result into a
    :class:`~flowfreq.core.FrequencyResults` -- never synthesising a field the
    Fortran did not report (design doc section 4's rule).

Everything here requires the built f2py extension
(``python build_fortran/build.py``); importing this module does not, but
calling ``run_fortran_reference``/``run_fortran_ema`` raises the same
``ImportError`` :mod:`flowfreq.peakfqr` raises when it is absent.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .core import AnalysisMethod, EMAParameters, FrequencyResults, grubbs_beck_critical_value
from .validation.reference import QMAX, QMIN, ReferenceResult

__all__ = [
    "EmafitArrays",
    "build_emafit_arrays",
    "run_fortran_reference",
    "run_fortran_ema",
    "quantile_frames",
]

# emafitpr's own "no regional information" sentinel for a *_mse argument --
# see CLAUDE.md's skew-MSE encoding and flowfreq._var_emab.NO_REGIONAL_INFO,
# which the confidence-interval path already uses for the same purpose.
_NO_REGIONAL_INFO = 1e10

# fortranWrappers.R's own low-outlier threshold encoding (~lines 93-127),
# verified directly rather than assumed: LOthresh >= 1e-6 is a fixed cfs
# threshold (log10'd); 1e-99 < LOthresh < 1e-6 disables the low-outlier test
# entirely (log10(Qmin), which lands in emafit.f:985's "-99 < gbthrsh0 < -6"
# no-test branch); LOthresh <= 1e-99 (no override) runs MGBT via the exact
# sentinel -99.0. Bulletin17C's own ``user_low_outlier_threshold`` is either
# None or a real cfs value, so the middle branch is unreachable from this
# library's public API today -- kept faithful to the source anyway since it
# costs nothing.
_FIXED_THRESHOLD_MIN = 1e-6
_DISABLED_THRESHOLD_MIN = 1e-99
_MGBT_SENTINEL = -99.0


@dataclass(frozen=True)
class EmafitArrays:
    """One record's ``emafitpr`` input arrays, plus what the adapter needs afterward.

    ``ql``, ``qu``, ``tl``, ``tu`` are log10(cfs); ``dtype`` is 1 only for a
    row carrying the USGS historic-peak flag. ``years`` is a parallel
    diagnostic array (not an ``emafitpr`` argument) recording which water year
    produced each row, in the same order.

    Attributes
    ----------
    systematic_peaks : dict of int to float
        Raw (real-space, not log10) systematic peak values actually used as
        *exact* observations, keyed by water year. A zero flow is recorded
        here as ``0.0`` (its real reported value), even though its row in
        ``ql``/``qu`` is clamped to ``log10(QMIN)`` -- this is what lets the
        adapter recover ``pilf_flows`` by comparing against the Fortran's
        returned threshold without reintroducing the clamp into the report.
    n_zeros : int
        Count of systematic peaks recorded as exactly zero.
    n_censored : int
        Count of rows built to fill a gap year against a declared perception
        threshold (``ql != qu``). Historic and exact-peak rows are not
        censored; MGBT-flagged low outliers are not censored *here* either --
        that censoring happens inside ``emafitpr`` itself, driven by
        ``gbthrsh0``, not by this builder.
    gbthrsh0 : float
        The low-outlier argument ``emafitpr`` expects: ``-99.0`` to run MGBT,
        or ``log10(threshold)`` to fix it, per the encoding above.
    """

    ql: np.ndarray
    qu: np.ndarray
    tl: np.ndarray
    tu: np.ndarray
    dtype: np.ndarray
    years: np.ndarray
    systematic_peaks: Dict[int, float] = field(default_factory=dict)
    n_zeros: int = 0
    n_censored: int = 0
    gbthrsh0: float = _MGBT_SENTINEL

    @property
    def n(self) -> int:
        return len(self.ql)


def _gbthrsh0(user_low_outlier_threshold: Optional[float]) -> float:
    """Map a user PILF override (cfs, or None) onto ``emafitpr``'s ``gbthrsh0``."""
    if user_low_outlier_threshold is None or user_low_outlier_threshold <= _DISABLED_THRESHOLD_MIN:
        return _MGBT_SENTINEL
    if user_low_outlier_threshold >= _FIXED_THRESHOLD_MIN:
        return float(np.log10(user_low_outlier_threshold))
    # 1e-99 < threshold < 1e-6 cfs: fortranWrappers.R's "NONE" branch, which
    # disables the low-outlier test outright. Not reachable through
    # Bulletin17C's public API (nobody declares a sub-microcfs threshold), but
    # transcribed for the same reason CLAUDE.md's other sentinel tables are:
    # so the mapping is complete, not curve-fit to the cases at hand.
    return float(np.log10(QMIN))


def build_emafit_arrays(
    peak_flows: Sequence[float],
    water_years: Optional[Sequence[int]] = None,
    historical_peaks: Optional[List[Tuple[int, float]]] = None,
    perception_thresholds: Optional[Dict[Tuple[int, int], float]] = None,
    user_low_outlier_threshold: Optional[float] = None,
    ema_params: Optional[EMAParameters] = None,
) -> EmafitArrays:
    """Translate a ``Bulletin17C`` input set into ``emafitpr``'s arrays.

    Follows ``siteQT`` (``vendor/peakfqr/R/readInputs.R``), not
    ``ExpectedMomentsAlgorithm._build_flow_intervals``: every water year that
    either carries an observation (systematic or historical) or falls inside
    a *non-vacuous* declared perception-threshold period gets a row; every
    other year -- unmeasured, with no threshold asserting what would have
    been recorded -- is omitted, matching site 12363000 (design doc section
    3's worked example).

    Parameters
    ----------
    peak_flows : sequence of float
        Annual peak flows in cfs. A value of exactly 0 is a real observation
        (a documented zero-flow year), not a gap -- ``siteQT`` gives it an
        exact interval at ``Qmin`` (``ql = qu = 1e-20``, log10 -20) rather
        than treating it as censored or missing.
    water_years : sequence of int, optional
        Water year for each entry of ``peak_flows``. Defaults to the same
        "count back from this year" convention
        ``ExpectedMomentsAlgorithm.__init__`` uses when omitted.
    historical_peaks : list of (int, float), optional
        Known historic floods outside (or gap-filling within) the systematic
        record. Every entry gets ``dtype = 1`` -- this argument *is* USGS's
        historic-peak flag in this API, there being no separate peak-code
        concept to check, so setting it for exactly these rows (and no
        others) already satisfies "dtype is 1 only for the historic flag,
        not every peak in the historical period."
    perception_thresholds : dict of (int, int) to float, optional
        Declared perception-threshold periods, lower bound only (upper is
        always ``Qmax``, matching every existing caller of this API). Applied
        in insertion order with a later period overwriting an earlier one for
        any year they both cover -- ``siteQT``'s own documented rule ("the
        last one specified is given priority"), reproduced here as a
        sequential overwrite of a year->threshold map rather than an actual
        priority sort, which is what the R implementation does too.
    user_low_outlier_threshold : float, optional
        PILF override in cfs. Mapped onto ``gbthrsh0``; does not affect the
        interval construction itself -- ``emafitpr`` does its own low-outlier
        censoring internally, driven by this argument, the same way MGBT
        does when no override is given.
    ema_params : EMAParameters, optional
        When given and it carries a historical period
        (``historical_start``/``historical_end``/``historical_threshold``),
        that period is folded in as an additional perception-threshold
        period -- lower priority than an explicit ``perception_thresholds``
        entry for the same years -- so a caller relying on
        ``Bulletin17C``'s auto-configuration (historical peaks given without
        an explicit threshold dict) gets the same historical censoring the
        native path would build.

    Returns
    -------
    EmafitArrays

    Raises
    ------
    ValueError
        Mismatched ``peak_flows``/``water_years`` lengths, a duplicate water
        year, or a negative discharge.
    """
    flows_arr = np.asarray(peak_flows, dtype=float)
    if water_years is None:
        end_year = datetime.now().year
        years_arr = np.arange(end_year - len(flows_arr) + 1, end_year + 1)
    else:
        years_arr = np.asarray(water_years, dtype=int)
    if len(flows_arr) != len(years_arr):
        raise ValueError(
            f"peak_flows has {len(flows_arr)} entries but water_years has "
            f"{len(years_arr)}; they must be parallel"
        )

    historical_peaks_list = list(historical_peaks or [])
    thresholds = dict(perception_thresholds or {})

    # --- year -> declared lower threshold, lowest priority first --------- #
    threshold_by_year: Dict[int, float] = {}
    if (
        ema_params is not None
        and ema_params.historical_start is not None
        and ema_params.historical_end is not None
        and ema_params.historical_threshold
    ):
        for year in range(ema_params.historical_start, ema_params.historical_end + 1):
            threshold_by_year[year] = float(ema_params.historical_threshold)
    for (start, end), lower in thresholds.items():
        for year in range(int(start), int(end) + 1):
            threshold_by_year[year] = float(lower)  # later entries win, per siteQT

    # --- observed rows ----------------------------------------------------#
    peak_by_year: Dict[int, float] = {}
    n_zeros = 0
    for flow_val, year_val in zip(flows_arr, years_arr):
        year = int(year_val)
        flow = float(flow_val)
        if year in peak_by_year:
            raise ValueError(f"duplicate water year {year} in peak_flows/water_years")
        if flow < 0:
            raise ValueError(f"negative discharge {flow} for water year {year}")
        if flow == 0.0:
            n_zeros += 1
        peak_by_year[year] = flow

    historical_by_year: Dict[int, float] = {}
    for hist_year, hist_flow in historical_peaks_list:
        year = int(hist_year)
        flow = float(hist_flow)
        if flow < 0:
            raise ValueError(f"negative historical discharge {flow} for water year {year}")
        historical_by_year[year] = flow

    all_years = set(peak_by_year) | set(historical_by_year) | set(threshold_by_year)

    rows: List[Tuple[int, float, float, float, float, int]] = []
    systematic_peaks: Dict[int, float] = {}
    n_censored = 0

    for year in sorted(all_years):
        if year in historical_by_year:
            flow = historical_by_year[year]
            ql = qu = QMIN if flow == 0.0 else flow
            # No declared threshold for this specific historic year: the
            # flood's own value is the smallest flow that would have been
            # noticed, the same fallback
            # ExpectedMomentsAlgorithm._build_flow_intervals uses
            # (`threshold = self._ema_params.historical_threshold or flow`).
            tl = threshold_by_year.get(year, flow if flow > 0.0 else QMIN)
            # siteQT clamps every tl below Qmin up to Qmin (readInputs.R line
            # 1046) unconditionally -- not only for gap-filled rows -- so a
            # historic peak sitting inside a declared-but-vacuous (0.0)
            # threshold period reads as "unrestricted", the same as no
            # threshold at all, rather than taking log10(0).
            tl = max(tl, QMIN)
            rows.append((year, ql, qu, tl, QMAX, 1))
        elif year in peak_by_year:
            flow = peak_by_year[year]
            ql = qu = QMIN if flow == 0.0 else flow
            tl = max(threshold_by_year.get(year, QMIN), QMIN)
            systematic_peaks[year] = flow
            rows.append((year, ql, qu, tl, QMAX, 0))
        else:
            # A gap year: no observation, but a perception threshold was
            # declared for it. siteQT (readInputs.R lines ~1030-1051):
            # a *non-vacuous* threshold (> Qmin) censors the year at its
            # lower bound; a declared-but-vacuous one (<= Qmin -- e.g. a
            # systematic-period threshold of literally 0.0) is treated
            # exactly like "no information" and dropped, same as an
            # undeclared gap year. Since undeclared years are never added to
            # `all_years` at all, only the vacuous case needs handling here.
            threshold = threshold_by_year[year]
            if threshold <= QMIN:
                continue
            rows.append((year, QMIN, threshold, threshold, QMAX, 0))
            n_censored += 1

    if len(rows) < 3:
        raise ValueError(
            f"only {len(rows)} row(s) after applying perception thresholds; "
            "emafitpr needs at least 3 to fit a skew"
        )

    rows.sort(key=lambda r: r[0])
    years_arr = np.array([r[0] for r in rows], dtype=int)
    ql = np.log10(np.array([r[1] for r in rows], dtype=float))
    qu = np.log10(np.array([r[2] for r in rows], dtype=float))
    tl = np.log10(np.array([r[3] for r in rows], dtype=float))
    tu = np.log10(np.array([r[4] for r in rows], dtype=float))
    dtype = np.array([r[5] for r in rows], dtype=np.int32)

    return EmafitArrays(
        ql=ql,
        qu=qu,
        tl=tl,
        tu=tu,
        dtype=dtype,
        years=years_arr,
        systematic_peaks=systematic_peaks,
        n_zeros=n_zeros,
        n_censored=n_censored,
        gbthrsh0=_gbthrsh0(user_low_outlier_threshold),
    )


def run_fortran_reference(
    peak_flows: Sequence[float],
    water_years: Optional[Sequence[int]] = None,
    historical_peaks: Optional[List[Tuple[int, float]]] = None,
    perception_thresholds: Optional[Dict[Tuple[int, int], float]] = None,
    user_low_outlier_threshold: Optional[float] = None,
    ema_params: Optional[EMAParameters] = None,
    regional_skew: Optional[float] = None,
    regional_skew_mse: Optional[float] = None,
    aeps: Optional[Sequence[float]] = None,
    eps: float = 0.90,
    weight_opt: int = 1,
    station_name: str = "",
) -> Tuple[ReferenceResult, EmafitArrays]:
    """Build the arrays and call the vendored Fortran through the f2py bridge.

    Parameters
    ----------
    See :func:`build_emafit_arrays` for the record-shape arguments.
    regional_skew, regional_skew_mse : float, optional
        ``None`` for either means "no regional skew supplied" -- mapped onto
        ``r_G_mse >= 1e10`` ("STATION SKEW", ``emafit.f`` lines 146/1266/2221),
        not onto 0 ("generalized, no error"), which is a different thing.
        When both are given they are passed straight through as ``r_G``/
        ``r_G_mse``: ``Bulletin17C``'s own ``regional_skew_mse`` is already a
        plain MSE in the same units and sign convention ``emafitpr`` expects
        (0 = generalized with no error, negative = generalized with
        ``MSE = -value``, positive = weighted), so no translation is needed
        beyond the None case.
    aeps : sequence of float, optional
        Annual exceedance probabilities to evaluate. Defaults to
        ``FloodFrequencyAnalysis.STANDARD_AEP``.
    eps : float, optional
        Confidence-interval coverage, 0.90 for a 90% interval.
    weight_opt : int, optional
        1=HWN, 2=ERL, 3=INV -- peakfq's default is HWN, matching
        ``ExpectedMomentsAlgorithm``'s own default weighting.
    station_name : str, optional
        Carried through to the result for reporting.

    Returns
    -------
    tuple of (ReferenceResult, EmafitArrays)
        The Fortran's raw output, and the arrays it was called with (needed
        afterward to derive ``pilf_flows``/``n_censored`` without ever
        reading them off ``ReferenceResult``, which does not report them).

    Raises
    ------
    ImportError
        The f2py extension is not built; see
        :func:`flowfreq.validation.reference.ReferenceResult.from_emafit`.
    """
    if aeps is None:
        from .bulletin17c import FloodFrequencyAnalysis

        aeps_arr = np.asarray(FloodFrequencyAnalysis.STANDARD_AEP, dtype=float)
    else:
        aeps_arr = np.asarray(aeps, dtype=float)

    arrays = build_emafit_arrays(
        peak_flows,
        water_years=water_years,
        historical_peaks=historical_peaks,
        perception_thresholds=perception_thresholds,
        user_low_outlier_threshold=user_low_outlier_threshold,
        ema_params=ema_params,
    )

    if regional_skew is None or regional_skew_mse is None:
        r_g, r_g_mse = 0.0, _NO_REGIONAL_INFO
    else:
        r_g, r_g_mse = float(regional_skew), float(regional_skew_mse)

    reference = ReferenceResult.from_emafit(
        ql=arrays.ql.tolist(),
        qu=arrays.qu.tolist(),
        tl=arrays.tl.tolist(),
        tu=arrays.tu.tolist(),
        dtype=arrays.dtype.tolist(),
        aeps=aeps_arr.tolist(),
        regional_skew=r_g,
        regional_skew_mse=r_g_mse,
        station_name=station_name,
        eps=eps,
        weight_opt=weight_opt,
        gbthrsh0=arrays.gbthrsh0,
    )
    return reference, arrays


def quantile_frames(reference: ReferenceResult) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """``(quantiles, confidence_limits)`` DataFrames, matching the native engine's shape.

    Shared by the adapter and by ``Bulletin17C``'s ``engine="fortran"``
    ``compute_quantiles``/``compute_confidence_limits`` re-invocation path
    (design doc section 5: Fortran quantiles come from ``qP3sub``'s exact
    gamma quantile at Fortran-computed AEPs, not from re-applying a K-factor
    to already-fitted moments the way the native path does for an AEP it
    was not originally fit at -- so a new AEP list means a fresh
    ``emafitpr`` call, which is what produces the ``ReferenceResult`` this
    function reads, not a recomputation from ``reference`` alone).

    Parameters
    ----------
    reference : ReferenceResult

    Returns
    -------
    tuple of (pd.DataFrame, pd.DataFrame)
        ``quantiles`` has ``aep``, ``return_period``, ``flow_cfs``,
        ``log_flow``, ``K_factor``; ``confidence_limits`` has ``aep``,
        ``return_period``, ``flow_cfs``, ``lower_5pct``, ``upper_5pct``.
    """
    p = reference.parameters
    aeps = sorted(reference.quantiles, reverse=True)
    quantiles = pd.DataFrame(
        {
            "aep": aeps,
            "return_period": [1.0 / a for a in aeps],
            "flow_cfs": [reference.quantiles[a] for a in aeps],
        }
    )
    quantiles["log_flow"] = np.log10(quantiles["flow_cfs"])
    # K = (log_flow - mean_log) / std_log is an exact algebraic rearrangement
    # of the Fortran's own reported quantile and moments, not a fabricated
    # value -- included so this DataFrame has the same columns
    # FloodFrequencyAnalysis.compute_quantiles produces.
    quantiles["K_factor"] = (quantiles["log_flow"] - p["mean_log"]) / p["std_log"]
    ci_aeps = sorted(reference.confidence_intervals, reverse=True)
    confidence_limits = pd.DataFrame(
        {
            "aep": ci_aeps,
            "return_period": [1.0 / a for a in ci_aeps],
            "flow_cfs": [reference.quantiles.get(a, np.nan) for a in ci_aeps],
            "lower_5pct": [reference.confidence_intervals[a][0] for a in ci_aeps],
            "upper_5pct": [reference.confidence_intervals[a][1] for a in ci_aeps],
        }
    )
    return quantiles, confidence_limits


def _frequency_results_from_reference(
    reference: ReferenceResult,
    arrays: EmafitArrays,
    regional_skew: Optional[float],
    regional_skew_mse: Optional[float],
) -> FrequencyResults:
    """Adapt a ``ReferenceResult`` into ``FrequencyResults``, per design doc section 4.

    **Never synthesises a field the Fortran did not report.** ``ema_iterations``
    and ``ema_converged`` are forced ``None`` -- ``emafitpr`` does not report an
    iteration count or a convergence flag, and a plausible-looking default
    (e.g. ``True``) would be a fabrication, not a report. ``mgb_critical_value``
    is the one field this recomputes rather than leaves ``None``:
    ``grubbs_beck_critical_value(n)`` is a deterministic function of sample
    size alone, already computed the same way for the *native* engine's own
    diagnostic report (``ExpectedMomentsAlgorithm.run_analysis``) rather than
    read off the fit -- so reporting it here is consistent with what the
    native engine already does with this exact field, not a new fabrication.

    One field-semantics mismatch, documented rather than silently
    "corrected": ``reference.n_systematic`` counts every ``dtype == 0`` row,
    including gap-year censoring rows this builder adds, while the *native*
    engine's ``n_systematic`` counts only uncensored rows. Passed through
    directly, as the design doc's mapping table says ("direct") -- this is a
    genuine difference in what the two engines mean by the field, not a bug
    in either one, and reinterpreting it here would be silently rewriting
    validation code outside this module's lane.

    Parameters
    ----------
    reference : ReferenceResult
        Live ``emafitpr`` output.
    arrays : EmafitArrays
        The arrays *reference* was computed from -- source of ``pilf_flows``
        (derived from the raw systematic peaks and the returned threshold,
        never returned by the Fortran itself) and ``n_censored``.
    regional_skew, regional_skew_mse : float, optional
        Echoed into ``skew_regional``/``skew_used`` bookkeeping -- ``emafitpr``
        does not return its own inputs.

    Returns
    -------
    FrequencyResults
    """
    p = reference.parameters
    has_regional = regional_skew is not None and regional_skew_mse is not None

    skew_station = p["skew_at_site"]
    # Only reported as a distinct "weighted" value when the caller actually
    # supplied regional information -- station-only (r_G_mse >= 1e10) makes
    # emafitpr's own "weighted" column numerically equal to its at-site one,
    # and reporting that as skew_weighted would claim a weighting that never
    # happened. This mirrors ExpectedMomentsAlgorithm.run_analysis's own
    # policy of leaving skew_weighted None when no regional skew is supplied.
    skew_weighted = p["skew_weighted"] if has_regional else None
    skew_used = skew_weighted if skew_weighted is not None else skew_station

    threshold = reference.low_outlier_threshold
    pilf_flows = (
        sorted(v for v in arrays.systematic_peaks.values() if v < threshold)
        if threshold > 0.0
        else []
    )

    quantiles, confidence_limits = quantile_frames(reference)

    return FrequencyResults(
        n_peaks=reference.n_peaks,
        n_systematic=reference.n_systematic,
        n_historical=reference.n_historical,
        n_censored=arrays.n_censored,
        n_low_outliers=reference.low_outlier_count,
        mean_log=p["mean_log"],
        std_log=p["std_log"],
        skew_station=skew_station,
        skew_regional=regional_skew,
        skew_weighted=skew_weighted,
        skew_used=skew_used,
        low_outlier_threshold=threshold,
        mgb_critical_value=grubbs_beck_critical_value(reference.n_peaks),
        method=AnalysisMethod.EMA,
        quantiles=quantiles,
        confidence_limits=confidence_limits,
        ema_iterations=None,
        ema_converged=None,
        skew_used_mse=p.get("mse_skew"),
        n_zeros=arrays.n_zeros,
        pilf_flows=pilf_flows,
    )


def run_fortran_ema(
    peak_flows: Sequence[float],
    water_years: Optional[Sequence[int]] = None,
    historical_peaks: Optional[List[Tuple[int, float]]] = None,
    perception_thresholds: Optional[Dict[Tuple[int, int], float]] = None,
    user_low_outlier_threshold: Optional[float] = None,
    ema_params: Optional[EMAParameters] = None,
    regional_skew: Optional[float] = None,
    regional_skew_mse: Optional[float] = None,
    aeps: Optional[Sequence[float]] = None,
    eps: float = 0.90,
    weight_opt: int = 1,
    station_name: str = "",
) -> Tuple[FrequencyResults, ReferenceResult, EmafitArrays]:
    """Run the Fortran EMA end to end: build, call, adapt.

    The thin combination ``Bulletin17C(engine="fortran")`` and
    ``compare_engines`` (``flowfreq.workflow``) both need -- the latter wants
    the ``ReferenceResult`` too, to hand to ``Bulletin17C.validate()``, and
    the former wants the raw arrays back so it can re-invoke ``emafitpr`` at a
    different AEP list later without rebuilding them.

    Returns
    -------
    tuple of (FrequencyResults, ReferenceResult, EmafitArrays)
    """
    reference, arrays = run_fortran_reference(
        peak_flows,
        water_years=water_years,
        historical_peaks=historical_peaks,
        perception_thresholds=perception_thresholds,
        user_low_outlier_threshold=user_low_outlier_threshold,
        ema_params=ema_params,
        regional_skew=regional_skew,
        regional_skew_mse=regional_skew_mse,
        aeps=aeps,
        eps=eps,
        weight_opt=weight_opt,
        station_name=station_name,
    )
    results = _frequency_results_from_reference(reference, arrays, regional_skew, regional_skew_mse)
    return results, reference, arrays
