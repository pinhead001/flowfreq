"""
flowfreq.workflow - High-level Bulletin 17C analysis entry points.

One call from annual peaks to a fitted frequency curve, plus the skew-variant
helpers that go with it. This is the layer a consumer wants when it does not
want to assemble :class:`~flowfreq.bulletin17c.Bulletin17C`, quantiles and
confidence limits by hand.

Everything here returns plain numbers and DataFrames of numbers. Turning those
into labelled, rounded, comma-separated strings is presentation and belongs to
the caller.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .bulletin17c import Bulletin17C
from .core import EMAParameters, FrequencyResults, kfactor_array

if TYPE_CHECKING:  # pragma: no cover - annotations only
    # Same layering reason as flowfreq.bulletin17c.Bulletin17C.validate: a
    # runtime import here would make every `import flowfreq.workflow` pull in
    # the validation subsystem, which flowfreq/__init__ leaves opt-in.
    # compare_engines imports what it needs lazily at call time.
    from .validation.comparisons import ComparisonResult
    from .validation.reference import ReferenceResult

logger = logging.getLogger(__name__)

#: Fallback regional skew, used only when a caller opts in with
#: ``use_default_skew=True``. **Unsourced**: Bulletin 17C (England and others,
#: 2019, p. 31) has no national skew value, and says the Bulletin 17B plate 1
#: estimates "are not recommended for use in flood frequency studies"; no
#: publication giving -0.302 has been found (0.302 is 0.55**2, plate 1's MSE).
#: Prefer a published regional skew -- see :mod:`flowfreq.regional_skew`.
B17C_DEFAULT_SKEW: float = -0.302

#: Standard error paired with :data:`B17C_DEFAULT_SKEW` under the same opt-in.
B17C_DEFAULT_SKEW_SE: float = 0.55

#: How the skew used by :func:`run_ffa` / :func:`compare_engines` was chosen.
SKEW_SOURCES: Tuple[str, ...] = ("user", "default", "station")


def resolve_regional_skew(
    regional_skew: Optional[float],
    regional_skew_se: Optional[float],
    use_default_skew: bool = False,
    station_skew_only: bool = False,
) -> Tuple[Optional[float], Optional[float], str]:
    """Settle which skew an analysis uses, refusing to pick one silently.

    Exactly one of three choices must be made: supply a regional skew and its
    standard error; opt in to the unsourced :data:`B17C_DEFAULT_SKEW`; or ask
    for the station skew alone. This mirrors peakfq 8.1.0, which stops with an
    error when a weighted or regional skew option has no ``GenSkew``/``SkewSE``
    (``vendor/peakfqr/R/main.R``) rather than defaulting.

    Parameters
    ----------
    regional_skew, regional_skew_se : float or None
        A published regional skew and its standard error. Both or neither.
    use_default_skew : bool
        Opt in to :data:`B17C_DEFAULT_SKEW` / :data:`B17C_DEFAULT_SKEW_SE`.
        Logs a warning, since the value has no located source.
    station_skew_only : bool
        Use the at-site skew with no regional weighting (peakfq's ``Station``).

    Returns
    -------
    tuple
        ``(regional_skew, regional_skew_mse, source)``; skew and MSE are
        ``None`` for station-only, and *source* is one of :data:`SKEW_SOURCES`.

    Raises
    ------
    ValueError
        When no choice, or more than one, is made, or only one of skew and
        standard error is given.
    """
    supplied = regional_skew is not None or regional_skew_se is not None
    chosen = int(supplied) + int(use_default_skew) + int(station_skew_only)
    if chosen == 0:
        raise ValueError(
            "No regional skew chosen. Pass regional_skew and regional_skew_se from a "
            "published study (see flowfreq.regional_skew.regional_skew_for), "
            "station_skew_only=True for at-site skew alone, or use_default_skew=True to "
            f"accept the unsourced fallback {B17C_DEFAULT_SKEW} (SE {B17C_DEFAULT_SKEW_SE}). "
            "Bulletin 17C gives no national default and does not recommend the 17B map."
        )
    if chosen > 1:
        raise ValueError(
            "Choose one skew source: a supplied regional_skew, use_default_skew, or "
            "station_skew_only -- not several."
        )
    if supplied:
        if regional_skew is None or regional_skew_se is None:
            raise ValueError(
                "regional_skew and regional_skew_se go together; a skew without its "
                "standard error cannot be weighted."
            )
        if regional_skew_se <= 0:
            raise ValueError(f"regional_skew_se must be positive, got {regional_skew_se}")
        return float(regional_skew), float(regional_skew_se) ** 2, "user"
    if use_default_skew:
        logger.warning(
            "Using the unsourced fallback regional skew %s (SE %s) at the caller's request; "
            "Bulletin 17C recommends a published regional skew study instead.",
            B17C_DEFAULT_SKEW,
            B17C_DEFAULT_SKEW_SE,
        )
        return B17C_DEFAULT_SKEW, B17C_DEFAULT_SKEW_SE**2, "default"
    return None, None, "station"


def peak_code_kwargs(
    peak_flows: Sequence[float],
    water_years: Optional[Sequence[int]],
    peak_codes: Optional[Sequence[object]],
    *,
    apply_peak_codes: bool = True,
    regional_skew: Optional[float] = None,
    regional_skew_mse: Optional[float] = None,
    user_low_outlier_threshold: Optional[float] = None,
    historical_peaks: Optional[Sequence[Tuple[int, float]]] = None,
    perception_thresholds: Optional[Mapping[Tuple[int, int], float]] = None,
    engines: Sequence[str] = ("native",),
    site_name: str = "site",
) -> Optional[Dict[str, Any]]:
    """``Bulletin17C`` arguments for a record whose qualification codes change it.

    The one place the analysis entry points (:func:`flowfreq.analyze_gage`,
    :func:`run_ffa`, :func:`compare_engines`, ``flowfreq compare``) apply NWIS
    peak codes, so they all apply them identically: through
    :func:`flowfreq.psf_convert.convert_peak_frame`, i.e. peakfq 8.1.0's
    ``siteQT`` rules with PeakFQ's default perception threshold for a new
    site. Code 7 peaks become historic; codes 3 and O, and 6 and C (peakfq's
    default ``Urb/Reg = No``), remove the peak, leaving its year as missing
    data; codes 4 and 8 make it a less-than / greater-than interval, passed to
    ``Bulletin17C`` as ``interval_peaks``.

    Parameters
    ----------
    peak_flows, water_years : sequence
        The record, as for :class:`~flowfreq.bulletin17c.Bulletin17C`.
    peak_codes : sequence or None
        One raw NWIS ``qualification_code`` per peak, aligned with
        ``peak_flows`` (blank/``None``/NaN for none). ``None`` means no codes.
    apply_peak_codes : bool, default True
        ``False`` fits every peak as an exact systematic value, as flowfreq
        did before codes were applied; codes that would have changed that are
        logged as ignored.
    regional_skew, regional_skew_mse, user_low_outlier_threshold
        As for ``Bulletin17C``; carried into the returned arguments.
    historical_peaks, perception_thresholds
        The caller's own historic information. With acted-on codes present
        they are refused rather than merged: the codes already say which
        peaks are historic, and merging a second account of the historic
        period into ``siteQT``'s rows is not something this reproduces.
    engines : sequence of {"native", "fortran"}
        Every engine that will fit the record. Each must build exactly
        ``siteQT``'s rows from the returned arguments
        (:meth:`~flowfreq.psf_convert.StationInputs.bulletin17c_kwargs`), so
        a comparison gives both engines the same record.
    site_name : str
        Used in messages only.

    Returns
    -------
    dict or None
        ``None`` when there is nothing to apply -- no codes, no code that
        ``siteQT`` acts on, or ``apply_peak_codes=False`` -- and the caller's
        own arguments stand unchanged. Otherwise ``Bulletin17C`` keyword
        arguments: ``peak_flows``, ``water_years``, ``historical_peaks``,
        ``perception_thresholds``, ``user_low_outlier_threshold``,
        ``regional_skew``, ``regional_skew_mse`` and ``interval_peaks`` (the
        code 4/8 peaks as ``(water_year, lower_cfs, upper_cfs)``, or ``None``).

    Raises
    ------
    ValueError
        ``peak_codes`` not aligned with the record, or acted-on codes
        together with ``historical_peaks``/``perception_thresholds``.
    flowfreq.psf_convert.UnsupportedSpecError
        A record an engine would not fit as ``siteQT`` does (the message
        names the years), such as a code 7 peak that is also code 4 or 8 -- a
        historic interval, which ``Bulletin17C`` has no argument for. Pass
        ``apply_peak_codes=False`` to fit the peaks as exact values instead,
        or use :meth:`~flowfreq.psf_convert.StationInputs.fortran_reference`.
    """
    if not engines:
        raise ValueError("engines must name at least one of 'native', 'fortran'")
    if peak_codes is None:
        return None
    codes = list(peak_codes)
    flows = np.asarray(peak_flows, dtype=float)
    if water_years is None:
        raise ValueError("peak_codes need water_years to say which year each code belongs to")
    years = np.asarray(water_years)
    if not len(codes) == len(flows) == len(years):
        raise ValueError(
            f"peak_codes has {len(codes)} entries but the record has {len(flows)} peak(s) "
            f"and {len(years)} water year(s); they must be aligned"
        )

    from .peak_codes import count_acted_on_codes

    acted_on = count_acted_on_codes(codes)
    if not acted_on:
        return None
    summary = ", ".join(f"code {c}: {n}" for c, n in acted_on.items())
    if not apply_peak_codes:
        logger.warning(
            "%s: peaks carry qualification codes peakfq would act on (%s); they are "
            "fitted as exact systematic peaks here because peak codes were turned off.",
            site_name,
            summary,
        )
        return None
    extra = [
        name
        for name, value in (
            ("historical_peaks", historical_peaks),
            ("perception_thresholds", perception_thresholds),
        )
        if value
    ]
    if extra:
        raise ValueError(
            f"{site_name}: the peaks carry qualification codes ({summary}), which say "
            f"themselves which peaks are historic or removed; do not also pass "
            f"{' or '.join(extra)}. Turn peak codes off to use your own historic "
            "information with the peaks fitted as exact values."
        )

    from .psf_convert import convert_peak_frame

    frame = pd.DataFrame(
        {"water_year": years.astype(int), "peak_flow_cfs": flows, "qualification_code": codes}
    )
    inputs = convert_peak_frame(
        frame,
        station_id=site_name,
        regional_skew=regional_skew,
        regional_skew_mse=regional_skew_mse,
        low_outlier_threshold=user_low_outlier_threshold,
    )
    # Each engine's check rebuilds its own rows; the arguments are the same.
    kwargs: Dict[str, Any] = {}
    for engine in engines:
        kwargs = inputs.bulletin17c_kwargs(engine)
    logger.info("%s: applying peak qualification codes as peakfq does (%s)", site_name, summary)
    return kwargs


#: Return intervals reported by :func:`run_ffa` and :func:`compute_skew_tables`.
DEFAULT_RETURN_INTERVALS: List[float] = [1.5, 2, 5, 10, 25, 50, 100, 200, 500]

#: The same series as annual exceedance probabilities.
#: = [0.667, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
DEFAULT_AEP: List[float] = [1 / ri for ri in DEFAULT_RETURN_INTERVALS]

#: Canonical skew option labels, in the order a report presents them.
SKEW_OPTIONS: List[str] = ["Station Skew", "Weighted Skew", "Regional Skew"]


def _low_outlier_source(override: Optional[float]) -> str:
    """Describe where the reported PILF threshold came from.

    Parameters
    ----------
    override : float or None
        The user's requested threshold, if any.

    Returns
    -------
    str
        ``"MGBT"`` or ``"override"``. Both EMA and MOM censor on a supplied
        threshold now, so the label needs no method-specific caveat.
    """
    if override is None:
        return "MGBT"
    return "override"


def run_ffa(
    peak_flows: np.ndarray,
    water_years: np.ndarray,
    regional_skew: Optional[float] = None,
    regional_skew_se: Optional[float] = None,
    perception_thresholds: Optional[List[dict]] = None,
    low_outlier_threshold_override: Optional[float] = None,
    *,
    use_default_skew: bool = False,
    station_skew_only: bool = False,
    peak_codes: Optional[Sequence[object]] = None,
    apply_peak_codes: bool = True,
) -> dict:
    """Run Bulletin 17C flood frequency analysis.

    Fits by EMA, falling back to MOM only when EMA fails to converge *and* no
    perception thresholds are in play -- MOM cannot represent censored
    intervals, so a threshold-bearing record keeps its non-converged EMA fit
    rather than silently losing the thresholds.

    Errors are returned, not raised: any failure comes back as a string under
    ``error`` with the other keys left at their empty defaults. That suits an
    interactive caller that wants to show the message rather than crash. A
    caller that would rather have an exception should check ``error`` and
    raise its own. The one exception is the skew choice: making none (or
    several) is a calling error, not a data error, and raises ``ValueError``
    before any fitting -- see :func:`resolve_regional_skew`.

    Parameters
    ----------
    peak_flows : np.ndarray
        Annual peak flows in cfs.
    water_years : np.ndarray
        Corresponding water years.
    regional_skew, regional_skew_se : float, optional
        A published regional skew and its standard error. There is no
        default: supply these, or set ``station_skew_only`` or
        ``use_default_skew``.
    perception_thresholds : list of dict, optional
        Each dict has keys ``start_year``, ``end_year``, ``threshold_cfs``
        (legacy) or ``lower_cfs`` / ``upper_cfs``.  Converts to the
        ``Dict[Tuple[int,int], float]`` format expected by
        :class:`~flowfreq.bulletin17c.Bulletin17C` and passed to EMA so that
        years in each period without a recorded peak are treated as
        left-censored observations (peak < threshold).
    low_outlier_threshold_override : float, optional
        User-supplied PILF threshold (cfs).  When > 0, overrides the MGBT
        result and censors all peaks below this value.  The threshold actually
        applied, its source and the resulting PILF count come back under
        ``parameters`` so a caller can show which cut produced the fit.
    use_default_skew : bool
        Explicitly accept the unsourced :data:`B17C_DEFAULT_SKEW`. Logs a warning.
    station_skew_only : bool
        Use the at-site skew with no regional weighting.
    peak_codes : sequence, optional
        One NWIS ``qualification_code`` per peak, aligned with
        ``peak_flows`` -- the ``qualification_code`` column of the frame
        :meth:`flowfreq.usgs.USGSgage.download_peak_flow` returns. When given,
        the codes are applied as peakfq 8.1.0 applies them (see
        :func:`peak_code_kwargs`): code 7 peaks are historic, 3/O/6/C peaks
        are removed. A record with no code peakfq acts on is fitted exactly
        as without ``peak_codes``. A code 4/8 (less-than/greater-than) peak
        cannot be fitted here and comes back as an ``error``, naming the
        years and codes.
    apply_peak_codes : bool, default True
        ``False`` ignores ``peak_codes`` (logging what was ignored) and fits
        every peak as an exact systematic value.

    Returns
    -------
    dict
        Keys: b17c, converged, method, parameters, quantile_df, error.
        ``parameters["peak_codes_applied"]`` counts, per code, the peaks
        whose treatment the codes changed (empty when none were applied).

    Raises
    ------
    ValueError
        When no skew source, or more than one, is chosen; when
        ``peak_codes`` is not aligned with ``peak_flows``; or when codes
        peakfq acts on are combined with ``perception_thresholds``.

    Examples
    --------
    >>> result = run_ffa(peak_flows, water_years, regional_skew=-0.07,
    ...                  regional_skew_se=0.36)
    >>> result["quantile_df"]["Flow (cfs)"]

    >>> result = run_ffa(peak_flows, water_years, station_skew_only=True,
    ...                  low_outlier_threshold_override=500.0)
    """
    skew, skew_mse, skew_source = resolve_regional_skew(
        regional_skew, regional_skew_se, use_default_skew, station_skew_only
    )
    if peak_codes is not None:
        if len(peak_codes) != len(peak_flows):
            raise ValueError(
                f"peak_codes has {len(peak_codes)} entries for {len(peak_flows)} peak(s); "
                "they must be aligned"
            )
        if apply_peak_codes and perception_thresholds:
            from .peak_codes import count_acted_on_codes

            acted_on = count_acted_on_codes(peak_codes)
            if acted_on:
                raise ValueError(
                    "peak_codes carrying codes peakfq acts on "
                    f"({', '.join(f'code {c}: {n}' for c, n in acted_on.items())}) cannot be "
                    "combined with perception_thresholds; pass apply_peak_codes=False to "
                    "use your thresholds with every peak fitted as an exact value."
                )
    result = {
        "b17c": None,
        "converged": False,
        "method": None,
        "parameters": {},
        "quantile_df": pd.DataFrame(),
        "error": None,
    }

    try:
        # Convert list of threshold dicts → Dict[Tuple[int,int], float]
        pt_dict: Optional[Dict[Tuple[int, int], float]] = None
        if perception_thresholds:
            pt_dict = {
                (int(t["start_year"]), int(t["end_year"])): float(t["threshold_cfs"])
                for t in perception_thresholds
                if float(t.get("threshold_cfs", 0)) > 0
            } or None

        lo_override = (
            float(low_outlier_threshold_override)
            if low_outlier_threshold_override and low_outlier_threshold_override > 0
            else None
        )

        b17c_kwargs: Dict[str, Any] = {
            "peak_flows": peak_flows,
            "water_years": water_years,
            "regional_skew": skew,
            "regional_skew_mse": skew_mse,
            "perception_thresholds": pt_dict,
            "user_low_outlier_threshold": lo_override,
        }
        coded = peak_code_kwargs(
            peak_flows,
            water_years,
            peak_codes,
            apply_peak_codes=apply_peak_codes,
            regional_skew=skew,
            regional_skew_mse=skew_mse,
            user_low_outlier_threshold=lo_override,
            perception_thresholds=pt_dict,
        )
        codes_applied: Dict[str, int] = {}
        if coded is not None and peak_codes is not None:
            from .peak_codes import count_acted_on_codes

            b17c_kwargs = coded
            codes_applied = count_acted_on_codes(list(peak_codes))
        b17c = Bulletin17C(**b17c_kwargs)

        b17c.run_analysis(method="ema")
        method = "ema"
        converged = bool(b17c.results.ema_converged)

        # Only fall back to MOM when the record is plain systematic peaks -- MOM has no
        # mechanism to incorporate censored intervals or historic peaks, so we keep the
        # (non-converged) EMA result when thresholds or applied codes extend the record.
        extended = (
            b17c_kwargs.get("perception_thresholds")
            or b17c_kwargs.get("historical_peaks")
            or b17c_kwargs.get("interval_peaks")
        )
        if not converged and not extended:
            logger.warning("EMA did not converge, falling back to MOM")
            b17c.run_analysis(method="mom")
            method = "mom"
            converged = True

        aep = np.array(DEFAULT_AEP)
        quantiles_df = b17c.compute_quantiles(aep=aep)
        ci_df = b17c.compute_confidence_limits(aep=aep)

        quantile_df = pd.DataFrame(
            {
                "Return Interval (yr)": DEFAULT_RETURN_INTERVALS,
                "AEP (%)": aep,
                "Flow (cfs)": quantiles_df["flow_cfs"].values,
                "Lower 90% CI": ci_df["lower_5pct"].values,
                "Upper 90% CI": ci_df["upper_5pct"].values,
            }
        )

        r = b17c.results
        result.update(
            {
                "b17c": b17c,
                "converged": converged,
                "method": method,
                "parameters": {
                    "mean_log": r.mean_log,
                    "std_log": r.std_log,
                    "skew_station": r.skew_station,
                    "skew_weighted": r.skew_weighted,
                    "skew_used": r.skew_used,
                    "regional_skew": skew,
                    "regional_skew_source": skew_source,
                    # The low-outlier cut and where it came from. Without these
                    # a caller could offer the override but never show its effect.
                    # Both EMA and MOM censor on it now, so the source is just
                    # whether it came from MGBT or the user.
                    "low_outlier_threshold": r.low_outlier_threshold,
                    "n_low_outliers": r.n_low_outliers,
                    "low_outlier_source": _low_outlier_source(lo_override),
                    "peak_codes_applied": codes_applied,
                },
                "quantile_df": quantile_df,
            }
        )

    except Exception as e:
        logger.exception("FFA analysis failed")
        result["error"] = str(e)

    return result


def _skew_values_from_result(ffa_result: dict) -> Dict[str, Optional[float]]:
    """Return the three skew values stored in an ffa_result dict."""
    p = ffa_result.get("parameters", {})
    return {
        "Station Skew": p.get("skew_station"),
        "Weighted Skew": p.get("skew_weighted"),
        "Regional Skew": p.get("regional_skew"),
    }


def compute_skew_tables(
    ffa_result: dict,
    selected_labels: List[str],
) -> Dict[str, pd.DataFrame]:
    """Compute a raw quantile+CI table for each selected skew option.

    Uses the LP3 moments (mean_log, std_log) already fitted by EMA/MOM and
    substitutes the requested skew value to produce separate frequency tables
    without re-running the full analysis.

    Parameters
    ----------
    ffa_result : dict
        Output from :func:`run_ffa`.
    selected_labels : list[str]
        Subset of :data:`SKEW_OPTIONS`.

    Returns
    -------
    dict[str, pd.DataFrame]
        Maps label → DataFrame with columns:
        ``Return Interval (yr)``, ``AEP (%)``,
        ``Flow (cfs)``, ``Lower 90% CI``, ``Upper 90% CI``.
        Returns an empty dict if ffa_result has an error.
    """
    if ffa_result.get("error") or ffa_result.get("b17c") is None:
        return {}

    r = ffa_result["b17c"].results
    mean_log = r.mean_log
    std_log = r.std_log
    n = r.n_systematic or r.n_peaks

    skew_map = _skew_values_from_result(ffa_result)
    aep = np.array(DEFAULT_AEP)
    z_alpha = 1.6449  # norm.ppf(0.95): two-sided 90% CI

    tables: Dict[str, pd.DataFrame] = {}
    for label in selected_labels:
        skew_val = skew_map.get(label)
        if skew_val is None:
            continue

        K = kfactor_array(skew_val, aep)
        log_Q = mean_log + K * std_log
        Q = 10.0**log_Q

        var_factor = 1 / n + K**2 * (1 + 0.75 * skew_val**2) / (2 * (n - 1))
        se_log = std_log * np.sqrt(var_factor)
        lower = 10.0 ** (log_Q - z_alpha * se_log)
        upper = 10.0 ** (log_Q + z_alpha * se_log)

        tables[label] = pd.DataFrame(
            {
                "Return Interval (yr)": DEFAULT_RETURN_INTERVALS,
                "AEP (%)": aep,
                "Flow (cfs)": Q,
                "Lower 90% CI": lower,
                "Upper 90% CI": upper,
            }
        )

    return tables


def build_skew_curves_dict(
    ffa_result: dict,
    selected_labels: List[str],
) -> Dict[str, float]:
    """Return ``{label: skew_value}`` for the selected skew options.

    Intended for passing directly to
    :func:`flowfreq.freq_plot.plot_frequency_curve` as the
    ``skew_curves`` argument.

    Parameters
    ----------
    ffa_result : dict
        Output from :func:`run_ffa`.
    selected_labels : list[str]
        Skew labels the caller has selected.

    Returns
    -------
    dict[str, float]
        Empty dict (fall back to default) when no valid labels are found.
    """
    skew_map = _skew_values_from_result(ffa_result)
    return {lbl: skew_map[lbl] for lbl in selected_labels if skew_map.get(lbl) is not None}


@dataclass
class EngineComparisonReport:
    """The result of running one record through both engines and comparing them.

    ``docs/FORTRAN_ENGINE_DESIGN.md`` section 5: "the comparison ... is the
    feature; ``engine=`` alone leaves the user diffing two analyses by hand."
    Built by :func:`compare_engines`, which reuses the already-existing
    :meth:`~flowfreq.bulletin17c.Bulletin17C.validate` /
    :class:`~flowfreq.validation.comparisons.FrequencyComparator` machinery
    rather than duplicating comparison logic.

    Attributes
    ----------
    native : FrequencyResults
        The native engine's fitted result, evaluated at the same AEPs the
        Fortran side was.
    reference : ReferenceResult
        The live ``emafitpr`` output.
    comparison : ComparisonResult
        Per-field differences and pass/fail, from
        :class:`~flowfreq.validation.comparisons.FrequencyComparator`.
    site_name : str
        Carried through for the markdown title only.
    """

    native: FrequencyResults
    reference: "ReferenceResult"
    comparison: "ComparisonResult"
    site_name: str = ""

    @property
    def max_quantile_deviation_pct(self) -> float:
        """Largest quantile percent difference across every AEP compared.

        Deliberately not :attr:`ComparisonResult.max_diff_pct` -- that also
        folds in parameter and confidence-interval differences, which answer
        a different question than "how far apart are the flood quantiles
        themselves."
        """
        return max(self.comparison.quantile_diffs.values(), default=0.0)

    def to_markdown(self) -> str:
        """A markdown report suitable for a submittal appendix.

        Parameters, skews (in skew units, not percent -- see
        :class:`~flowfreq.validation.comparisons.ComparisonResult`), quantiles
        and confidence limits, each as a table of native vs. Fortran with the
        percent (or, for skew, absolute) difference already computed by
        :class:`~flowfreq.validation.comparisons.FrequencyComparator`.
        """
        n, r, c = self.native, self.reference, self.comparison
        lines: List[str] = []

        title = "Engine comparison"
        if self.site_name:
            title += f": {self.site_name}"
        lines.append(f"# {title}")
        lines.append("")
        status = "PASS" if c.passed else "FAIL"
        lines.append(
            f"**{status}** -- native EMA vs. USGS peakfq (`emafitpr`, the vendored Fortran "
            "reference)"
        )
        lines.append("")
        lines.append(
            f"- Max quantile deviation: **{self.max_quantile_deviation_pct:.3f}%** "
            f"(tolerance {c.tolerance_pct:.2f}%)"
        )
        lines.append(f"- {c.summary}")
        lines.append("")

        native_params = {"mean_log": n.mean_log, "std_log": n.std_log}
        if c.parameter_diffs:
            lines.append("## Parameters")
            lines.append("")
            lines.append("| Parameter | Native | Fortran | Diff % |")
            lines.append("|---|---:|---:|---:|")
            for key in sorted(c.parameter_diffs):
                nat_val = native_params.get(key)
                ref_val = r.parameters.get(key)
                nat_str = f"{nat_val:.6g}" if nat_val is not None else "n/a"
                ref_str = f"{ref_val:.6g}" if ref_val is not None else "n/a"
                lines.append(f"| {key} | {nat_str} | {ref_str} | {c.parameter_diffs[key]:.3f} |")
            lines.append("")

        if c.skew_diffs:
            native_skew = {"skew_at_site": n.skew_station, "skew_weighted": n.skew_weighted}
            lines.append("## Skew (skew units, not percent)")
            lines.append("")
            lines.append("| Parameter | Native | Fortran | Abs diff |")
            lines.append("|---|---:|---:|---:|")
            for key in sorted(c.skew_diffs):
                nat_val = native_skew.get(key)
                ref_val = r.parameters.get(key)
                nat_str = f"{nat_val:.4f}" if nat_val is not None else "n/a"
                ref_str = f"{ref_val:.4f}" if ref_val is not None else "n/a"
                lines.append(f"| {key} | {nat_str} | {ref_str} | {c.skew_diffs[key]:.4f} |")
            lines.append("")

        if c.quantile_diffs:
            native_q: Dict[float, float] = {}
            if n.quantiles is not None and not n.quantiles.empty:
                native_q = dict(zip(n.quantiles["aep"], n.quantiles["flow_cfs"]))
            lines.append("## Quantiles (cfs)")
            lines.append("")
            lines.append("| AEP | Return period (yr) | Native | Fortran | Diff % |")
            lines.append("|---:|---:|---:|---:|---:|")
            for aep in sorted(c.quantile_diffs, reverse=True):
                nat_val = native_q.get(aep)
                ref_val = r.quantiles.get(aep)
                nat_str = f"{nat_val:,.0f}" if nat_val is not None else "n/a"
                ref_str = f"{ref_val:,.0f}" if ref_val is not None else "n/a"
                lines.append(
                    f"| {aep:.4g} | {1 / aep:,.1f} | {nat_str} | {ref_str} | "
                    f"{c.quantile_diffs[aep]:.3f} |"
                )
            lines.append("")

        if c.ci_diffs:
            native_ci: Dict[float, Tuple[float, float]] = {}
            if n.confidence_limits is not None and not n.confidence_limits.empty:
                for _, row in n.confidence_limits.iterrows():
                    native_ci[float(row["aep"])] = (row["lower_5pct"], row["upper_5pct"])
            lines.append("## Confidence intervals (cfs)")
            lines.append("")
            lines.append(
                "| AEP | Native lower | Native upper | Fortran lower | Fortran upper | Diff % |"
            )
            lines.append("|---:|---:|---:|---:|---:|---:|")
            for aep in sorted(c.ci_diffs, reverse=True):
                nat_lo, nat_hi = native_ci.get(aep, (None, None))
                ref_lo, ref_hi = r.confidence_intervals.get(aep, (None, None))
                nat_lo_s = f"{nat_lo:,.0f}" if nat_lo is not None else "n/a"
                nat_hi_s = f"{nat_hi:,.0f}" if nat_hi is not None else "n/a"
                ref_lo_s = f"{ref_lo:,.0f}" if ref_lo is not None else "n/a"
                ref_hi_s = f"{ref_hi:,.0f}" if ref_hi is not None else "n/a"
                lines.append(
                    f"| {aep:.4g} | {nat_lo_s} | {nat_hi_s} | {ref_lo_s} | {ref_hi_s} | "
                    f"{c.ci_diffs[aep]:.3f} |"
                )
            lines.append("")

        return "\n".join(lines)


def compare_engines(
    peak_flows: np.ndarray,
    water_years: Optional[np.ndarray] = None,
    regional_skew: Optional[float] = None,
    regional_skew_se: Optional[float] = None,
    historical_peaks: Optional[List[Tuple[int, float]]] = None,
    perception_thresholds: Optional[Dict[Tuple[int, int], float]] = None,
    user_low_outlier_threshold: Optional[float] = None,
    ema_params: Optional[EMAParameters] = None,
    aeps: Optional[np.ndarray] = None,
    site_name: str = "",
    tolerance_pct: float = 1.0,
    parameter_tolerance_pct: float = 0.5,
    ci_tolerance_pct: float = 2.0,
    *,
    use_default_skew: bool = False,
    station_skew_only: bool = False,
    peak_codes: Optional[Sequence[object]] = None,
    apply_peak_codes: bool = True,
) -> EngineComparisonReport:
    """Run one record through both engines and compare them.

    ``docs/FORTRAN_ENGINE_DESIGN.md`` section 5-- this is the feature the
    Fortran bridge exists for: a per-run comparison against USGS peakfq 8.1.0
    on the caller's own data, not just the four sites already committed as
    parity goldens.

    **Requires the built f2py extension, and raises the same actionable
    ``ImportError`` :mod:`flowfreq.peakfqr` raises when it is absent -- there
    is no golden-file fallback.** Decided, not left open (design doc section
    9, question 4): a fallback would only ever cover the four sites already
    proven to agree, which is exactly where this comparison proves the least;
    on a caller's own record there is no golden to fall back to, so the
    feature would silently work on some inputs and not others.

    Parameters
    ----------
    peak_flows, water_years, historical_peaks, perception_thresholds,
    user_low_outlier_threshold, ema_params : see :class:`~flowfreq.bulletin17c.Bulletin17C`.
    regional_skew, regional_skew_se, use_default_skew, station_skew_only
        As in :func:`run_ffa`: exactly one skew source must be chosen. Both
        engines get the same skew, so the comparison is like for like either way.
    aeps : array-like, optional
        Annual exceedance probabilities to compare at. Defaults to
        :attr:`~flowfreq.bulletin17c.FloodFrequencyAnalysis.STANDARD_AEP`.
        Both engines are evaluated at exactly this list -- the native side's
        own quantiles/confidence limits are recomputed at *aeps* even when it
        matches the default, so the two are always compared like for like
        rather than at whatever AEPs the native fit happened to be run at.
    site_name : str, optional
        Carried through to :meth:`EngineComparisonReport.to_markdown`'s title.
    tolerance_pct, parameter_tolerance_pct, ci_tolerance_pct : float
        As in :meth:`~flowfreq.bulletin17c.Bulletin17C.validate`.
    peak_codes : sequence, optional
        One NWIS ``qualification_code`` per peak, aligned with
        ``peak_flows``. Applied as in :func:`run_ffa`, through
        :func:`peak_code_kwargs`, and checked against *both* engines: each
        must build exactly peakfq's ``siteQT`` rows from the converted
        arguments, so the comparison stays like for like. Codes peakfq acts
        on cannot be combined with ``historical_peaks`` or
        ``perception_thresholds``.
    apply_peak_codes : bool, default True
        ``False`` ignores ``peak_codes`` and fits every peak as exact.

    Returns
    -------
    EngineComparisonReport

    Raises
    ------
    ValueError
        When no skew source, or more than one, is chosen, or ``peak_codes``
        cannot be applied (see :func:`peak_code_kwargs`).
    flowfreq.psf_convert.UnsupportedSpecError
        A coded record either engine would not fit as ``siteQT`` does; the
        message names the years. Code 4/8 peaks are not such a record: both
        engines get them as ``interval_peaks``.
    ImportError
        The f2py extension is not built; run
        ``python build_fortran/build.py`` (needs gfortran and meson).
    """
    skew, skew_mse, _ = resolve_regional_skew(
        regional_skew, regional_skew_se, use_default_skew, station_skew_only
    )

    # Codes first: a record that cannot be coded is the caller's to fix, and
    # saying so should not need the extension.
    coded = peak_code_kwargs(
        peak_flows,
        water_years,
        peak_codes,
        apply_peak_codes=apply_peak_codes,
        regional_skew=skew,
        regional_skew_mse=skew_mse,
        user_low_outlier_threshold=user_low_outlier_threshold,
        historical_peaks=historical_peaks,
        perception_thresholds=perception_thresholds,
        engines=("native", "fortran"),
        site_name=site_name or "site",
    )
    interval_peaks = None
    if coded is not None:
        peak_flows = coded["peak_flows"]
        water_years = coded["water_years"]
        historical_peaks = coded["historical_peaks"]
        perception_thresholds = coded["perception_thresholds"]
        user_low_outlier_threshold = coded["user_low_outlier_threshold"]
        interval_peaks = coded.get("interval_peaks")

    import flowfreq.peakfqr  # noqa: F401 -- raise before doing any native work if absent

    from .bulletin17c import FloodFrequencyAnalysis
    from .fortran_engine import run_fortran_reference

    if aeps is None:
        aeps = FloodFrequencyAnalysis.STANDARD_AEP
    aeps = np.asarray(aeps, dtype=float)

    native = Bulletin17C(
        peak_flows=peak_flows,
        water_years=water_years,
        regional_skew=skew,
        regional_skew_mse=skew_mse,
        historical_peaks=historical_peaks,
        perception_thresholds=perception_thresholds,
        ema_params=ema_params,
        user_low_outlier_threshold=user_low_outlier_threshold,
        interval_peaks=interval_peaks,
    )
    native.run_analysis(method="ema", engine="native")
    # Recompute at *aeps* explicitly: run_analysis() always fits at
    # STANDARD_AEP internally, so a caller-supplied aeps list would otherwise
    # leave native.results.quantiles at a different AEP set than the
    # reference below, and FrequencyComparator would silently compare nothing
    # (no aep keys in common) rather than raise.
    native.results.quantiles = native.compute_quantiles(aep=aeps)
    native.results.confidence_limits = native.compute_confidence_limits(aep=aeps)

    reference, _arrays = run_fortran_reference(
        peak_flows,
        water_years=water_years,
        historical_peaks=historical_peaks,
        perception_thresholds=perception_thresholds,
        user_low_outlier_threshold=user_low_outlier_threshold,
        ema_params=ema_params,
        regional_skew=skew,
        regional_skew_mse=skew_mse,
        aeps=aeps,
        station_name=site_name,
        interval_peaks=interval_peaks,
    )

    comparison = native.validate(
        reference,
        tolerance_pct=tolerance_pct,
        parameter_tolerance_pct=parameter_tolerance_pct,
        ci_tolerance_pct=ci_tolerance_pct,
    )
    return EngineComparisonReport(
        native=native.results, reference=reference, comparison=comparison, site_name=site_name
    )
