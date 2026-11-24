"""The ``b17b_ffa`` section: a Bulletin 17B (superseded) fit, for comparison with 17C.

Nothing else in a report computes Bulletin 17B. The section takes the peaks and the skew
options of a ``gaged_ffa`` section, so the two differ in procedure only, unless the caller
gives the generalized skew of an earlier study.
"""

from __future__ import annotations

import copy
import math
from contextlib import contextmanager
from typing import Any, Dict, Iterator, List, Mapping, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.special import ndtr, ndtri

from ..aep import same_aep
from ..bulletin17b import (
    AEP_RANGE,
    STANDARD_AEP,
    Bulletin17B,
    Bulletin17BError,
    Bulletin17BResults,
    _threshold_row,
    bulletin17b_kwargs,
    confidence_deviates,
    removed_peaks,
)
from ..regulation import RegulatedRecordError, classify_site, require_unregulated
from . import gaged
from ._base import SectionContent, _figure_png, _version, cfs
from .blocks import Block, Bullets, Figure, Heading, Paragraph, Table
from .gaged import _peak_frame, _resolve_skew, _skew_sentence, _threshold_rows
from .spec import ReportError

_KIND = "b17b_ffa"

#: The label on everything that shows a Bulletin 17B number.
SUPERSEDED = "Bulletin 17B (superseded)"

#: The caption of the section's main table.
B17B_CAPTION = f"{SUPERSEDED} quantiles, with 90-percent limits."

#: The reference every Bulletin 17B section cites.
BULLETIN_17B = (
    "Interagency Advisory Committee on Water Data, 1982, Guidelines for determining flood "
    "flow frequency, Bulletin 17B of the Hydrology Subcommittee (revised September 1981, "
    "editorial corrections March 1982): U.S. Geological Survey, Office of Water Data "
    "Coordination, https://pubs.usgs.gov/unnumbered/70275162/report.pdf."
)

#: The section's AEPs: the bulletin's table and the 5- and 25-year floods.
DEFAULT_AEP: Tuple[float, ...] = tuple(sorted({*STANDARD_AEP, 0.2, 0.04}, reverse=True))

#: Generalized-skew MSE when none is given (Plate I).
_PLATE_I_MSE = 0.302

_REGIONAL = ("regional_skew", "regional_skew_se")
_LOOKUPS = ("regional_skew_state", "regional_skew_region", "regional_skew_from_site")


# -- options ------------------------------------------------------------------------------


@contextmanager
def _own_errors() -> Iterator[None]:
    """Name this kind in the errors of the helpers shared with ``gaged_ffa``."""
    try:
        yield
    except ReportError as exc:
        raise ReportError(str(exc).replace("gaged_ffa: ", f"{_KIND}: ", 1)) from exc


def _number(options: Mapping[str, Any], key: str, *, positive: bool = False) -> Optional[float]:
    """An option as a finite number, or None when it is not given."""
    value = options.get(key)
    if value is None:
        return None
    try:
        number = math.nan if isinstance(value, bool) else float(value)
    except (TypeError, ValueError):
        number = math.nan
    if not math.isfinite(number):
        raise ReportError(f"{_KIND}: {key} must be a number; got {value!r}")
    if positive and number <= 0.0:
        raise ReportError(f"{_KIND}: {key} must be positive; got {number:g}")
    return number


def _period(value: Any) -> Optional[List[int]]:
    """The ``historic_period`` option as ``[first, last]``."""
    if value is None:
        return None
    if (
        not isinstance(value, (list, tuple))
        or len(value) != 2
        or any(isinstance(y, bool) or not isinstance(y, int) for y in value)
    ):
        raise ReportError(
            f"{_KIND}: historic_period must be [first year, last year]; got {value!r}"
        )
    if value[1] < value[0]:
        raise ReportError(
            f"{_KIND}: historic_period ends before it begins: {value[0]} to {value[1]}"
        )
    return list(value)


def _aeps(value: Any) -> List[float]:
    """The ``aeps`` option, or the section's default."""
    if value is None:
        return list(DEFAULT_AEP)
    if not isinstance(value, (list, tuple)) or not value:
        raise ReportError(f"{_KIND}: aeps must be a non-empty list of AEPs")
    try:
        aeps = [float(a) for a in value]
    except (TypeError, ValueError) as exc:
        raise ReportError(f"{_KIND}: aeps must be a list of numbers") from exc
    if any(not AEP_RANGE[0] <= a <= AEP_RANGE[1] for a in aeps):
        raise ReportError(f"{_KIND}: aeps must be between {AEP_RANGE[0]:g} and {AEP_RANGE[1]:.6f}")
    for i, a in enumerate(aeps):
        for b in aeps[:i]:
            if same_aep(a, b):
                raise ReportError(f"{_KIND}: aeps {b:g} and {a:g} are the same AEP given twice")
    return aeps


def _low_threshold(options: Mapping[str, Any]) -> Optional[float]:
    """The ``low_outlier_threshold`` option; zero means none, as in ``gaged_ffa``."""
    value = _number(options, "low_outlier_threshold")
    if value is not None and value < 0.0:
        raise ReportError(f"{_KIND}: low_outlier_threshold must not be negative; got {value:g}")
    return value or None


def _codes(frame: pd.DataFrame) -> List[Optional[str]]:
    """Each peak's qualification codes as text, None where it has none."""
    codes: List[Optional[str]] = []
    for year, code in zip(frame["water_year"], frame["qualification_code"]):
        if code is None or (isinstance(code, float) and math.isnan(code)):
            codes.append(None)
        elif isinstance(code, str):
            codes.append(code)
        elif (
            isinstance(code, (int, float))
            and not isinstance(code, bool)
            and float(code).is_integer()
        ):
            codes.append(str(int(code)))
        else:
            raise ReportError(f"{_KIND}: the qualification code of {year} cannot be read: {code!r}")
    return codes


def _checked_rows(thresholds: Any) -> List[Dict[str, Any]]:
    """The perception thresholds with whole years, which may be given as text, and a flow."""
    with _own_errors():
        rows = [dict(row) for row in _threshold_rows(thresholds)]
    for row in rows:
        for key in ("start_year", "end_year"):
            if isinstance(row[key], str):
                try:
                    row[key] = int(row[key])
                except ValueError as exc:
                    raise ReportError(
                        f"{_KIND}: a threshold's {key} must be a year; got {row[key]!r}"
                    ) from exc
        try:
            _threshold_row(row)
        except Bulletin17BError as exc:
            raise ReportError(f"{_KIND}: {exc}") from exc
    return rows


def _skew(options: Dict[str, Any]) -> Tuple[Optional[float], float, str, Optional[str]]:
    """The generalized skew, its MSE, the sentence that says where they come from, a citation.

    A ``generalized_skew`` wins. Without one, the regional skew of the ``gaged_ffa``
    options is resolved, in place, as that section resolves it. ``station_skew_only`` goes
    with neither.
    """
    if options.get("station_skew_only"):
        others = [
            key
            for key in ("generalized_skew", "generalized_skew_mse", *_REGIONAL, *_LOOKUPS)
            if options.get(key) is not None and options.get(key) is not False
        ]
        if others:
            raise ReportError(f"{_KIND}: station_skew_only cannot go with {', '.join(others)}")
        return None, _PLATE_I_MSE, "The station skew was used alone, as asked.", None
    generalized = _number(options, "generalized_skew")
    mse = _number(options, "generalized_skew_mse", positive=True)
    regional = _number(options, "regional_skew")
    se = _number(options, "regional_skew_se", positive=True)
    if (regional is None) != (se is None):
        raise ReportError(
            f"{_KIND}: regional_skew and regional_skew_se go together; a skew without its "
            "standard error cannot be weighted"
        )
    if generalized is not None:
        used = _PLATE_I_MSE if mse is None else mse
        sentence = (
            "The station skew was weighted with the generalized skew from the earlier study, "
            f"{generalized:.3f} (mean square error {used:.3f}"
            + (", the Plate I value" if mse is None else "")
            + ")"
        )
        if regional is not None or any(options.get(key) for key in _LOOKUPS):
            sentence += ", in place of the regional skew"
        return generalized, used, sentence + ".", None
    if mse is not None:
        raise ReportError(f"{_KIND}: generalized_skew_mse needs generalized_skew")
    with _own_errors():
        _resolve_skew(options)
    if options.get("regional_skew") is None:
        return (
            None,
            _PLATE_I_MSE,
            "No regional or generalized skew was given, so the station skew was used alone.",
            None,
        )
    sentence = _skew_sentence(options, "user") + (
        " It is the regional skew a Bulletin 17C analysis with the same options uses, so the "
        "difference is the change in procedure."
    )
    citation = options.get("regional_skew_citation")
    return (
        float(options["regional_skew"]),
        float(options["regional_skew_se"]) ** 2,
        sentence,
        str(citation) if citation else None,
    )


# -- what the section says ----------------------------------------------------------------


def _count(n: int, one: str, many: str) -> str:
    return f"{n} {one if n == 1 else many}"


def _peak_list(peaks: List[Tuple[int, float]]) -> str:
    return ", ".join(f"{y} ({cfs(q)} cfs)" for y, q in peaks) or "none"


def _record_rows(
    options: Mapping[str, Any], frame: pd.DataFrame, site_no: Optional[str]
) -> List[Tuple[str, str]]:
    """The station and the record as given, before any Bulletin 17B treatment."""
    years = frame["water_year"].to_numpy(dtype=int)
    flows = frame["peak_flow_cfs"].to_numpy(dtype=float)
    rows: List[Tuple[str, str]] = []
    if site_no:
        rows.append(("USGS station number", site_no))
    if options.get("site_name"):
        rows.append(("Station name", str(options["site_name"])))
    area = options.get("drainage_area_sqmi")
    if area is not None:
        rows.append(("Drainage area (square miles)", f"{float(area):,.4g}"))
    rows.append(("Annual peaks", f"{len(frame)}"))
    rows.append(("Water years", f"{years.min()} to {years.max()}"))
    for name, i in (("Largest", int(np.argmax(flows))), ("Smallest", int(np.argmin(flows)))):
        rows.append((f"{name} peak (cfs)", f"{cfs(flows[i])} (water year {years[i]})"))
    if options["peaks"].get("source"):
        retrieved = options["peaks"].get("retrieved")
        rows.append(
            (
                "Data source",
                str(options["peaks"]["source"]) + (f", retrieved {retrieved}" if retrieved else ""),
            )
        )
    return rows


def _removed_sentence(removed: List[Tuple[int, float, str]]) -> str:
    listed = "; ".join(
        f"{year} ({cfs(flow)} cfs, code{'s' if ',' in codes else ''} {codes})"
        for year, flow, codes in removed
    )
    return f"Peaks removed from the record for their qualification codes: {listed}."


def _threshold_sentence(rows: List[Dict[str, Any]]) -> str:
    listed = []
    for row in rows:
        start, end, threshold = _threshold_row(row)
        level = f"{cfs(threshold)} cfs" if math.isfinite(threshold) else "not observed"
        listed.append(f"{start} to {end}, {level}")
    return (
        "Perception thresholds: "
        + "; ".join(listed)
        + ". Bulletin 17B takes only a historic period from them."
    )


def _parameter_rows(res: Bulletin17BResults) -> List[Tuple[str, str]]:
    """The fitted parameters, the adjustments made and the outlier tests."""
    h, c = res.historic, res.conditional
    if c is not None:
        skew, suffix = "Synthetic skew", " (synthetic)"
    elif h is not None:
        skew, suffix = "Historically adjusted skew", " (historically adjusted)"
    else:
        skew, suffix = "Station skew", ""
    rows = [
        ("Mean of log10 flow" + suffix, f"{res.mean_log:.4f}"),
        ("Standard deviation of log10 flow" + suffix, f"{res.std_log:.4f}"),
        (skew, f"{res.skew_station:.3f}"),
        ("Mean square error of that skew", f"{res.skew_station_mse:.3f}"),
    ]
    if res.skew_generalized is not None and res.skew_generalized_mse is not None:
        rows.append(("Generalized skew", f"{res.skew_generalized:.3f}"))
        rows.append(
            ("Mean square error of the generalized skew", f"{res.skew_generalized_mse:.3f}")
        )
    if res.skew_weighted is not None:
        rows.append(("Weighted skew", f"{res.skew_weighted:.3f}"))
    rows.append(("Skew used", f"{res.skew_used:.3f}"))
    rows.append(("Years in the systematic record", f"{res.n_years}"))
    rows.append(("Systematic peaks in the moments", f"{res.n_systematic}"))
    if res.n_zeros:
        rows.append(("Zero years", f"{res.n_zeros}"))
    if res.n_below_base:
        rows.append(("Years below the gage base", f"{res.n_below_base}"))
    if h is not None:
        rows.append(("Years in the historic period, H", f"{h.h}"))
        rows.append(("Historic peaks, Z", f"{h.z}"))
        rows.append(("Weight of a systematic peak, W", f"{h.w:.4f}"))
    if c is not None:
        rows.append(("Probability of a peak above the truncation level, P~", f"{c.p_tilde:.4f}"))
    rows.append(("High-outlier threshold (cfs)", cfs(res.high_outlier_threshold)))
    rows.append(("High outliers", _peak_list(res.high_outliers)))
    rows.append(("Low-outlier threshold (cfs)", cfs(res.low_outlier_threshold)))
    rows.append(("Low outliers", _peak_list(res.low_outliers)))
    return rows


def _quantile_table(res: Bulletin17BResults) -> Table:
    return Table(
        (
            "Annual exceedance probability (percent)",
            "Recurrence interval (years)",
            "Flow (cfs)",
            "Lower 90-percent limit (cfs)",
            "Upper 90-percent limit (cfs)",
            "Expected probability (percent)",
        ),
        tuple(
            (
                f"{row['aep'] * 100:.3g}",
                f"{1.0 / row['aep']:.3g}",
                cfs(row["flow_cfs"]),
                cfs(row["lower"]),
                cfs(row["upper"]),
                f"{row['expected_probability'] * 100:.3g}",
            )
            for _, row in res.quantiles.iterrows()
        ),
        ("right",) * 6,
        B17B_CAPTION,
    )


# -- figure -------------------------------------------------------------------------------


def _frequency_figure(res: Bulletin17BResults, where: str) -> Any:
    """The fitted curve, its 90-percent limits and the peaks at their plotting positions.

    Zero and below-base years have no flow to put on a log axis and are left out.
    """
    import matplotlib.pyplot as plt

    from ..core import kfactor_array
    from ..freq_plot import _ANNOT_FONT_SIZE, _FONT_SIZE

    positions = res.plotting_positions
    shown = positions[positions["flow_cfs"] > 0.0]
    exceedance = shown["exceedance"].to_numpy(dtype=float)
    rare = min(0.002, float(res.quantiles["aep"].min()), float(exceedance.min()))
    common = max(0.999, float(res.quantiles["aep"].max()), float(exceedance.max()))

    fig, ax = plt.subplots(figsize=(10, 6))
    historic = {y for y, _ in res.historic.historic_peaks} if res.historic else set()
    low = {y for y, _ in res.low_outliers}
    groups = (
        ("Annual peaks", shown[~shown["year"].isin(historic | low)], "o", 20, "steelblue"),
        ("Historic peaks", shown[shown["year"].isin(historic)], "D", 36, "darkorange"),
        ("Low outliers", shown[shown["year"].isin(low)], "x", 40, "red"),
    )
    for label, rows, marker, size, color in groups:
        if not rows.empty:
            x = -ndtri(rows["exceedance"].to_numpy(dtype=float))
            ax.scatter(x, rows["flow_cfs"], s=size, c=color, marker=marker, zorder=5, label=label)

    x = np.linspace(-float(ndtri(common)), -float(ndtri(rare)), 300)
    k = kfactor_array(res.skew_used, ndtr(-x))
    k_up, k_lo = confidence_deviates(k, res.n_years, float(ndtri(0.95)))
    flow, upper, lower = (10.0 ** (res.mean_log + d * res.std_log) for d in (k, k_up, k_lo))
    ax.plot(x, flow, "k-", linewidth=1, label="Fitted frequency curve", zorder=4)
    ax.fill_between(x, lower, upper, color="k", alpha=0.12)
    ax.plot(x, lower, "k--", linewidth=0.7, alpha=0.55, label="90-percent limits")
    ax.plot(x, upper, "k--", linewidth=0.7, alpha=0.55)

    ax.set_yscale("log")
    flows = shown["flow_cfs"]
    exponents = range(
        math.floor(math.log10(float(flows.min()))), math.ceil(math.log10(float(flows.max()))) + 1
    )
    ax.set_yticks([10.0**e for e in exponents])
    ax.set_yticklabels([f"{10**e:,}" if e >= 0 else f"{10.0**e:g}" for e in exponents])
    ax.set_ylabel("Discharge (cfs)", fontsize=_FONT_SIZE)
    ticks = [0.999, 0.99, 0.95, 0.90, 0.80, 0.50, 0.20, 0.10, 0.04, 0.02, 0.01, 0.005, 0.002]
    ax.set_xticks([-float(ndtri(p)) for p in ticks])
    ax.set_xticklabels(
        [f"{100 * p:g}%" for p in ticks], fontsize=_ANNOT_FONT_SIZE, rotation=45, ha="right"
    )
    ax.set_xlim(x[0], x[-1])
    ax.set_xlabel("Annual Exceedance Probability", fontsize=_FONT_SIZE)
    title = f"Flood Frequency Curve, {SUPERSEDED}"
    ax.set_title(f"{title}\n{where}" if where else title, fontsize=_FONT_SIZE, fontweight="bold")
    ax.legend(loc="upper left", fontsize=_ANNOT_FONT_SIZE)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    return fig


def _figure_caption(res: Bulletin17BResults) -> str:
    positions = "Weibull" if res.historic is None else "historically weighted (Appendix 6)"
    caption = (
        f"{SUPERSEDED} frequency curve with 90-percent limits, and the annual peaks at "
        f"{positions} plotting positions."
    )
    missing = []
    if res.n_zeros:
        missing.append(_count(res.n_zeros, "zero year", "zero years"))
    if res.n_below_base:
        years = _count(res.n_below_base, "year", "years")
        missing.append(f"{years} below the gage base")
    if missing:
        caption += f" Not shown: {' and '.join(missing)}."
    return caption


# -- the section --------------------------------------------------------------------------


def _build_b17b_ffa(given: Dict[str, Any]) -> SectionContent:
    """Build the section.

    The options are those of ``gaged_ffa`` where they apply, with one difference. A record
    with code-7 (historic) peaks and ``perception_thresholds`` needs
    ``apply_peak_codes=False`` in ``gaged_ffa``, which does not take the two together. Here
    it needs the codes applied, because Bulletin 17B gets its historic peaks from them, so
    such a record takes two sets of options.
    """
    return _b17b_ffa(given)[0]


def _b17b_ffa(given: Dict[str, Any]) -> Tuple[SectionContent, Optional[Bulletin17BResults]]:
    """The section and the Bulletin 17B results behind it, None when Bulletin 17B refuses.

    A refusal (:class:`~flowfreq.bulletin17b.Bulletin17BError`, or a regulated record) is a
    sentence in the section, so the rest of the report builds. Options that cannot be read
    and a failed download or skew lookup raise :class:`ReportError`. ``given`` is not changed.
    """
    options = copy.deepcopy(given)
    _number(options, "drainage_area_sqmi", positive=True)
    threshold = _low_threshold(options)
    aeps = _aeps(options.get("aeps"))
    period = _period(options.get("historic_period"))
    if period is not None:
        options["historic_period"] = period
    rows = _checked_rows(options.get("perception_thresholds"))
    site_no = str(options["site_no"]).strip() if options.get("site_no") else None
    if options.get("peaks") is None:
        if site_no is None:
            raise ReportError(f"{_KIND}: give site_no (to download the peaks) or peaks")
        fetched = gaged.download_peaks(site_no)
        options["peaks"] = fetched["peaks"]
        for key in ("site_name", "drainage_area_sqmi"):
            if options.get(key) is None and fetched.get(key) is not None:
                options[key] = fetched[key]
    with _own_errors():
        frame = _peak_frame(options["peaks"])
    skew, skew_mse, skew_sentence, citation = _skew(options)

    codes = _codes(frame)
    has_codes = any(c is not None for c in codes)
    apply_codes = bool(options.get("apply_peak_codes", True))
    applied = codes if has_codes and apply_codes else None
    flows = frame["peak_flow_cfs"].tolist()
    years = frame["water_year"].tolist()

    where = " ".join(
        str(x) for x in (f"USGS {site_no}" if site_no else None, options.get("site_name")) if x
    )
    title = f"Flood frequency{' at ' + where if where else ''}, {SUPERSEDED}"
    record = _record_rows(options, frame, site_no)
    data: List[str] = []
    override = False
    basis = ""
    references = [BULLETIN_17B]
    refusal: Optional[str] = None
    results: Optional[Bulletin17BResults] = None
    try:
        if site_no:
            # The codes count here whether or not they are applied, as in gaged_ffa.
            classed = classify_site(site_no, codes if has_codes else None)
            basis = "; ".join(classed.basis)
            screen = require_unregulated(
                classed, allow_regulated=bool(options.get("allow_regulated", False))
            )
            record.append(
                ("Regulation screen", str(screen["site_class"]) + (f": {basis}" if basis else ""))
            )
            override = bool(screen.get("regulation_override"))
        if has_codes and not apply_codes:
            data.append(
                "NWIS qualification codes were not applied: every peak is an exact systematic "
                "value."
            )
        removed = removed_peaks(flows, years, applied)
        if removed:
            data.append(_removed_sentence(removed))
        if rows:
            data.append(_threshold_sentence(rows))
        kwargs = bulletin17b_kwargs(flows, years, applied, rows or None)
        read = [
            f"{n} {what}"
            for n, what in (
                (len(kwargs["historic_peaks"]), "historic (code 7)"),
                (len(kwargs["below_base_years"]), "below the gage base (code 4)"),
            )
            if n
        ]
        if read:
            data.append(f"Peaks by their qualification codes: {', '.join(read)}.")
        derived = kwargs["historic_period"]
        if period is not None:
            if derived is not None and list(derived) != period:
                data.append(
                    f"The historic period given, {period[0]} to {period[1]}, was used in place "
                    f"of the one the codes and perception thresholds give, {derived[0]} to "
                    f"{derived[1]}."
                )
            kwargs["historic_period"] = (period[0], period[1])
        results = Bulletin17B(
            **kwargs,
            generalized_skew=skew,
            generalized_skew_mse=skew_mse,
            low_outlier_threshold=threshold,
            aeps=aeps,
            hand_computation=bool(options.get("hand_computation", False)),
        ).run_analysis()
    except RegulatedRecordError:
        # The screen's own message is about Bulletin 17C.
        refusal = (
            f"USGS {site_no} is classed regulated ({basis}), and Bulletin 17B describes "
            "unregulated flood peaks. The section's allow_regulated option overrides this and "
            "fits the record as observed"
        )
    except Bulletin17BError as exc:
        refusal = str(exc)

    blocks: List[Block] = [
        Table(("Item", "Value"), tuple(record), ("left", "left"), "Station and peak-flow record.")
    ]
    if override:
        blocks.append(
            Paragraph(
                "**This record is classified as regulated.** Bulletin 17B describes "
                "unregulated flood peaks; the fit was made at the analyst's direction and "
                "describes the regulated record as observed."
            )
        )
    if data:
        blocks.append(Paragraph(" ".join(data)))
    if results is None:
        sentence = (
            f"Bulletin 17B gives no answer for this record: {str(refusal).rstrip().rstrip('.')}."
        )
        blocks.append(Paragraph(sentence))
        return SectionContent(title, blocks, options, references, [sentence]), None

    if citation:
        references.append(citation)
    blocks.append(Heading(2, f"Procedure, {SUPERSEDED}"))
    method = (
        "The annual peaks were fitted with a log-Pearson Type III distribution by the "
        "procedure of Bulletin 17B (Interagency Advisory Committee on Water Data, 1982), "
        "which was superseded by Bulletin 17C. The fit is given for comparison only. The "
        f"computation was made with flowfreq {_version()}. {skew_sentence}"
    )
    if options.get("hand_computation"):
        method += " The roundings of the bulletin's worked examples were made (hand computation)."
    blocks.append(Paragraph(method))
    blocks.append(Bullets(tuple(results.steps)))
    blocks.append(
        Table(
            ("Parameter", "Value"),
            tuple(_parameter_rows(results)),
            ("left", "right"),
            f"{SUPERSEDED} fitted parameters and outlier tests.",
        )
    )
    blocks.append(_quantile_table(results))
    blocks.append(
        Paragraph(
            "The expected probability is that of Bulletin 17B Appendix 11 (eq. 11-1) for each "
            "flow. It is not applied to the flows."
        )
    )
    blocks.append(Figure(_figure_png(_frequency_figure(results, where)), _figure_caption(results)))
    return SectionContent(title, blocks, options, references), results


__all__ = ["B17B_CAPTION", "BULLETIN_17B", "DEFAULT_AEP", "SUPERSEDED"]
