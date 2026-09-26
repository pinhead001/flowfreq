"""
Screening and ranking donor basins on basin similarity, not drainage area alone.

The drainage-area ratio screen in :mod:`flowfreq.transpose`
(:data:`~flowfreq.transpose.DEFAULT_AREA_RATIO_RANGE`) is necessary but not
sufficient. Two basins of similar size can still differ in the predictors a
regional regression carries besides area -- mean annual precipitation above
all in a basin with an orographic gradient -- and a transposition between them
inherits that difference as bias. Issue #13 measured it in the Methow: a donor
pair the area band calls acceptable was still 14% low, and the whole residual
was the precipitation and canopy terms the area ratio drops.

This module adds two opt-in tools. Neither changes any existing number:
drainage-area screening in :mod:`flowfreq.transpose` is untouched.

:func:`screen_donors`
    The area-ratio screen, plus an optional per-predictor ratio band for each
    basin characteristic the caller names (``PRECPRIS10``, ``CANOPY_PCT``, ...).
    It reports every donor with the criterion that failed it rather than
    raising, so a caller can see *why* a donor was set aside. With no
    predictors it is exactly the area screen.

:func:`rank_donors_by_similarity`
    A region-of-influence distance in standardized basin-characteristic
    space (Burn, 1990):

        d(t, j) = sqrt( sum_k w_k * ((x_tk - x_jk) / s_k) ** 2 )

    where ``x`` are (optionally log10-transformed) attributes, ``s_k`` is the
    standard deviation of attribute ``k`` across the pool, and ``w_k`` are
    weights -- equal unless the caller supplies them. Burn's formulation is the
    one the USGS region-of-influence regressions build on (Tasker, Hodge and
    Barks, 1996), which log-transform drainage area; that is the one transform
    applied by default here.

:func:`similarity_donor_chooser`
    Adapts the ranking to :func:`flowfreq.qppq.loocv_qppq`'s ``donor_chooser``
    hook, so a similarity-ranked donor can be checked against the default
    correlation-ranked one on the same gage network.

**Refuse rather than guess.** A target missing an attribute raises; a donor
missing one is excluded with the missing names stated in its ``reason``. No
value is ever imputed. Attributes supplied as
:class:`flowfreq.streamstats.Characteristic` objects carry a unit, and a donor
whose unit disagrees with the target's is excluded rather than compared.

**What this is not.** It is not the map-correlation method (Archfield and
Vogel, 2010), which selects a QPPQ donor for a fully ungaged site by kriging
the cross-correlation field; that remains unimplemented. The thresholds in
``similarity_band`` are local choices, as the area band is; this module ships
none for predictors.

References
----------
Burn, D.H., 1990, Evaluation of regional flood frequency analysis with a region
of influence approach: Water Resources Research, v. 26, no. 10, p. 2257-2265.

Tasker, G.D., Hodge, S.A., and Barks, C.S., 1996, Region of influence
regression for estimating the 50-year flood at ungaged sites: Journal of the
American Water Resources Association, v. 32, no. 1, p. 163-170.

Archfield, S.A., and Vogel, R.M., 2010, Map correlation method: Selection of a
reference streamgage to estimate daily streamflow at ungaged catchments: Water
Resources Research, v. 46, W10513.
"""

from __future__ import annotations

import logging
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from flowfreq.transpose import DEFAULT_AREA_RATIO_RANGE

logger = logging.getLogger(__name__)

#: A site's basin characteristics: StreamStats code to a number, or to an
#: object with ``value`` (and optionally ``unit``) such as
#: :class:`flowfreq.streamstats.Characteristic`.
SiteAttributes = Mapping[str, Any]

#: Attributes log10-transformed by default in :func:`rank_donors_by_similarity`.
#: Drainage area spans orders of magnitude and every USGS regional regression
#: (and the ROI regressions of Tasker and others, 1996) enters it in logs.
#: Nothing else is transformed unless the caller says so.
DEFAULT_LOG_ATTRIBUTES: Tuple[str, ...] = ("DRNAREA",)

#: ``screen_donors`` status values.
STATUS_PASS = "pass"
STATUS_FAIL = "fail"
STATUS_EXCLUDED = "excluded"

#: ``rank_donors_by_similarity`` status values.
STATUS_RANKED = "ranked"


def _read_attribute(
    attributes: SiteAttributes, code: str
) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    """Return ``(value, unit, problem)`` for one attribute of one site.

    ``problem`` is a human-readable reason the value cannot be used, or None.
    Never substitutes a value for a missing one.
    """
    if code not in attributes:
        return None, None, "missing"
    raw = attributes[code]
    unit: Optional[str] = None
    if hasattr(raw, "value"):
        unit = getattr(raw, "unit", None)
        raw = raw.value
    if raw is None:
        return None, unit, "missing"
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None, unit, f"not numeric ({raw!r})"
    if not np.isfinite(value):
        return None, unit, "not finite"
    return value, unit, None


def _validate_band(name: str, band: Tuple[float, float]) -> Tuple[float, float]:
    low, high = float(band[0]), float(band[1])
    if not (np.isfinite(low) and np.isfinite(high)) or low <= 0 or high < low:
        raise ValueError(
            f"band for {name!r} must be (low, high) with 0 < low <= high; got {band!r}"
        )
    return low, high


def screen_donors(
    target_area_sqmi: float,
    donor_areas: Mapping[str, float],
    *,
    target_predictors: Optional[SiteAttributes] = None,
    donor_predictors: Optional[Mapping[str, SiteAttributes]] = None,
    similarity_band: Optional[Mapping[str, Tuple[float, float]]] = None,
    area_ratio_range: Tuple[float, float] = DEFAULT_AREA_RATIO_RANGE,
) -> pd.DataFrame:
    """Screen donors on drainage-area ratio and, optionally, predictor ratios.

    Every ratio is target over donor, the convention of
    :func:`flowfreq.transpose.transpose_frequency`. A donor passes when its
    area ratio lies inside ``area_ratio_range`` and, for each predictor named
    in ``similarity_band``, its predictor ratio lies inside that predictor's
    band. Nothing raises for a failing donor: this reports, and the caller (or
    the transpose call it goes on to make) decides.

    Parameters
    ----------
    target_area_sqmi : float
        Target drainage area, square miles.
    donor_areas : mapping of str to float
        Donor name to drainage area, same units.
    target_predictors : mapping, optional
        Target basin characteristics, keyed by code (e.g. ``"PRECPRIS10"``).
        Values may be numbers or :class:`flowfreq.streamstats.Characteristic`.
        Required when ``similarity_band`` is given.
    donor_predictors : mapping of str to mapping, optional
        Donor name to that donor's basin characteristics.
    similarity_band : mapping of str to (float, float), optional
        Predictor code to the accepted ``(low, high)`` target/donor ratio. No
        default is shipped: like the area band, the threshold is local. When
        omitted, the screen is area-only and matches today's behaviour.
    area_ratio_range : tuple of float, optional
        Defaults to :data:`flowfreq.transpose.DEFAULT_AREA_RATIO_RANGE`.

    Returns
    -------
    pd.DataFrame
        One row per donor: ``donor``, ``donor_area_sqmi``, ``area_ratio``,
        ``area_ok``, one ``ratio_<CODE>`` column per screened predictor,
        ``similarity_ok`` (None when a predictor is missing), ``status``
        (``"pass"``, ``"fail"``, or ``"excluded"`` when a donor predictor is
        missing or unusable), and ``reason`` naming every criterion that failed.
        Sorted passing donors first, then by ``|log(area_ratio)|``.

    Raises
    ------
    ValueError
        If the target area or a donor area is not positive, ``donor_areas`` is
        empty, a band is malformed, or the target lacks (or has a non-positive
        value for) a predictor in ``similarity_band``.
    """
    if not np.isfinite(target_area_sqmi) or target_area_sqmi <= 0:
        raise ValueError(f"target_area_sqmi must be positive; got {target_area_sqmi!r}")
    if not donor_areas:
        raise ValueError("no donor areas supplied; nothing to screen")
    area_low, area_high = _validate_band("area ratio", area_ratio_range)

    bands: Dict[str, Tuple[float, float]] = {
        code: _validate_band(code, band) for code, band in (similarity_band or {}).items()
    }
    target_values: Dict[str, Tuple[float, Optional[str]]] = {}
    if bands:
        if target_predictors is None:
            raise ValueError("similarity_band was given but target_predictors was not")
        for code in bands:
            value, unit, problem = _read_attribute(target_predictors, code)
            if problem is not None or value is None:
                raise ValueError(
                    f"target predictor {code!r} is {problem}; cannot screen on it. "
                    "Supply it or drop it from similarity_band -- it will not be imputed."
                )
            if value <= 0:
                raise ValueError(
                    f"target predictor {code!r} is {value!r}; a ratio screen needs positive "
                    "values. Use rank_donors_by_similarity for attributes that can be zero."
                )
            target_values[code] = (value, unit)
    donor_predictors = donor_predictors or {}

    rows: List[Dict[str, Any]] = []
    for donor, area in donor_areas.items():
        area = float(area)
        if not np.isfinite(area) or area <= 0:
            raise ValueError(f"donor {donor!r}: area must be positive, got {area!r}")
        area_ratio = float(target_area_sqmi) / area
        area_ok = bool(area_low <= area_ratio <= area_high)
        reasons: List[str] = []
        if not area_ok:
            reasons.append(f"area ratio {area_ratio:.3f} outside {(area_low, area_high)}")

        row: Dict[str, Any] = {
            "donor": donor,
            "donor_area_sqmi": area,
            "area_ratio": area_ratio,
            "area_ok": area_ok,
        }
        unusable: List[str] = []
        dissimilar = False
        attrs = donor_predictors.get(donor, {})
        for code, (low, high) in bands.items():
            target_value, target_unit = target_values[code]
            value, unit, problem = _read_attribute(attrs, code)
            if problem is None and target_unit and unit and unit != target_unit:
                problem = f"in {unit!r}, target in {target_unit!r}"
            if problem is None and value is not None and value <= 0:
                problem = f"not positive ({value!r})"
            if problem is not None or value is None:
                unusable.append(f"{code} {problem}")
                row[f"ratio_{code}"] = float("nan")
                continue
            ratio = target_value / value
            row[f"ratio_{code}"] = ratio
            if not low <= ratio <= high:
                dissimilar = True
                reasons.append(f"{code} ratio {ratio:.3f} outside {(low, high)}")

        if unusable:
            row["similarity_ok"] = None
            row["status"] = STATUS_EXCLUDED
            reasons.insert(0, "donor predictor " + ", ".join(unusable))
        else:
            row["similarity_ok"] = not dissimilar
            row["status"] = STATUS_PASS if area_ok and not dissimilar else STATUS_FAIL
        row["reason"] = "; ".join(reasons)
        rows.append(row)

    frame = pd.DataFrame(rows)
    order = {STATUS_PASS: 0, STATUS_FAIL: 1, STATUS_EXCLUDED: 2}
    frame["_order"] = frame["status"].map(order)
    frame["_log_ratio"] = np.abs(np.log(frame["area_ratio"].to_numpy(dtype=float)))
    frame = frame.sort_values(["_order", "_log_ratio"], kind="stable")
    frame = frame.drop(columns=["_order", "_log_ratio"]).reset_index(drop=True)

    for _, failed in frame[frame["status"] != STATUS_PASS].iterrows():
        logger.info("donor %r %s: %s", failed["donor"], failed["status"], failed["reason"])
    return frame


def rank_donors_by_similarity(
    target_attributes: SiteAttributes,
    donor_attributes: Mapping[str, SiteAttributes],
    attributes: Sequence[str],
    *,
    weights: Optional[Mapping[str, float]] = None,
    log_attributes: Sequence[str] = DEFAULT_LOG_ATTRIBUTES,
    scale: Optional[Mapping[str, float]] = None,
) -> pd.DataFrame:
    """Rank donors by region-of-influence distance to the target (Burn, 1990).

    ``d = sqrt(sum_k w_k * ((x_target,k - x_donor,k) / s_k) ** 2)`` over the
    named ``attributes``, after log10-transforming those in ``log_attributes``.
    Smaller is more similar; an identical basin has distance zero.

    Parameters
    ----------
    target_attributes : mapping
        Target basin characteristics, keyed by code. Values may be numbers or
        :class:`flowfreq.streamstats.Characteristic`.
    donor_attributes : mapping of str to mapping
        Donor name to that donor's basin characteristics.
    attributes : sequence of str
        The characteristic codes to compare on. Required: there is no default
        attribute set, because which predictors matter is regional.
    weights : mapping of str to float, optional
        Non-negative weight per attribute. Default equal weights (1 each).
        When given, every attribute must have one -- a partial mapping is
        refused rather than filled in.
    log_attributes : sequence of str, optional
        Attributes compared as log10. Default :data:`DEFAULT_LOG_ATTRIBUTES`.
        A non-positive value of a logged attribute is unusable.
    scale : mapping of str to float, optional
        Standard deviation ``s_k`` per attribute, in transformed units. Supply
        the regional network's values (as Burn and the USGS ROI studies do) to
        make distances independent of which donors happen to be candidates.
        Default: the sample standard deviation over the target and the usable
        donors, which requires at least two usable donors.

    Returns
    -------
    pd.DataFrame
        One row per donor: ``donor``, ``distance``, ``rank`` (1 = most
        similar; NaN if excluded), ``status`` (``"ranked"`` or
        ``"excluded"``), ``reason``, and one ``z_<CODE>`` column per attribute
        holding the standardized difference ``(x_donor - x_target) / s``, so
        a caller can see which attribute drives a distance. Ranked donors
        first, by ascending distance.

    Raises
    ------
    ValueError
        If ``attributes`` is empty or repeated, no donors are supplied,
        ``weights`` or ``scale`` is malformed, the target lacks an attribute
        or has an unusable value, or the pool is too small to estimate ``scale``.
    """
    codes = list(attributes)
    if not codes:
        raise ValueError("attributes is empty; name the basin characteristics to compare on")
    if len(set(codes)) != len(codes):
        raise ValueError(f"attributes contains duplicates: {codes!r}")
    if not donor_attributes:
        raise ValueError("no donors supplied; nothing to rank")

    if weights is None:
        weight = {code: 1.0 for code in codes}
    else:
        unknown = set(weights) - set(codes)
        absent = [code for code in codes if code not in weights]
        if unknown or absent:
            raise ValueError(
                f"weights must name exactly the attributes; unknown {sorted(unknown)!r}, "
                f"missing {absent!r}"
            )
        weight = {code: float(weights[code]) for code in codes}
        if any(not np.isfinite(w) or w < 0 for w in weight.values()):
            raise ValueError(f"weights must be finite and non-negative; got {weight!r}")
        if not any(w > 0 for w in weight.values()):
            raise ValueError("at least one weight must be positive")

    logged = set(log_attributes)

    def _transformed(site: SiteAttributes, code: str) -> Tuple[Optional[float], Optional[str]]:
        value, unit, problem = _read_attribute(site, code)
        if problem is None and value is not None and code in logged:
            if value <= 0:
                return None, f"{code} not positive ({value!r}) and compared in log10"
            value = float(np.log10(value))
        if problem is not None:
            return None, f"{code} {problem}"
        return value, unit

    target_x: Dict[str, float] = {}
    target_unit: Dict[str, Optional[str]] = {}
    for code in codes:
        value, detail = _transformed(target_attributes, code)
        if value is None:
            raise ValueError(
                f"target attribute unusable: {detail}. Supply it or drop it from "
                "attributes -- it will not be imputed."
            )
        target_x[code] = value
        target_unit[code] = detail

    usable: Dict[str, Dict[str, float]] = {}
    excluded: Dict[str, str] = {}
    for donor, site in donor_attributes.items():
        values: Dict[str, float] = {}
        problems: List[str] = []
        for code in codes:
            value, detail = _transformed(site, code)
            if value is None:
                problems.append(str(detail))
                continue
            if target_unit[code] and detail and detail != target_unit[code]:
                problems.append(f"{code} in {detail!r}, target in {target_unit[code]!r}")
                continue
            values[code] = value
        if problems:
            excluded[donor] = "; ".join(problems)
            logger.info("donor %r excluded from similarity ranking: %s", donor, excluded[donor])
        else:
            usable[donor] = values

    if scale is not None:
        absent = [code for code in codes if code not in scale]
        if absent:
            raise ValueError(f"scale is missing attributes {absent!r}")
        s = {code: float(scale[code]) for code in codes}
        if any(not np.isfinite(v) or v <= 0 for v in s.values()):
            raise ValueError(f"scale values must be finite and positive; got {s!r}")
    else:
        if len(usable) < 2:
            raise ValueError(
                f"only {len(usable)} donor(s) have every attribute; estimating the "
                "standardizing scale needs the target plus at least two. Supply scale "
                "(e.g. the regional network's standard deviations)."
            )
        s = {}
        for code in codes:
            pool = np.array([target_x[code]] + [v[code] for v in usable.values()])
            s[code] = float(np.std(pool, ddof=1))
            if s[code] == 0:
                logger.info(
                    "attribute %s is identical across the pool; it contributes nothing", code
                )

    rows: List[Dict[str, Any]] = []
    for donor in donor_attributes:
        row: Dict[str, Any] = {"donor": donor}
        if donor in usable:
            total = 0.0
            for code in codes:
                diff = usable[donor][code] - target_x[code]
                z = diff / s[code] if s[code] > 0 else 0.0
                row[f"z_{code}"] = z
                total += weight[code] * z * z
            row.update(distance=float(np.sqrt(total)), status=STATUS_RANKED, reason="")
        else:
            for code in codes:
                row[f"z_{code}"] = float("nan")
            row.update(distance=float("nan"), status=STATUS_EXCLUDED, reason=excluded[donor])
        rows.append(row)

    frame = pd.DataFrame(rows)
    frame["_order"] = (frame["status"] != STATUS_RANKED).astype(int)
    frame = frame.sort_values(["_order", "distance"], kind="stable").drop(columns="_order")
    frame = frame.reset_index(drop=True)
    ranked = frame["status"] == STATUS_RANKED
    frame["rank"] = np.where(ranked, np.cumsum(ranked), np.nan)
    leading = ["donor", "distance", "rank", "status", "reason"]
    return frame[leading + [f"z_{code}" for code in codes]]


def similarity_donor_chooser(
    site_attributes: Mapping[str, SiteAttributes],
    attributes: Sequence[str],
    *,
    weights: Optional[Mapping[str, float]] = None,
    log_attributes: Sequence[str] = DEFAULT_LOG_ATTRIBUTES,
    scale: Optional[Mapping[str, float]] = None,
) -> Callable[[str, pd.DataFrame, Dict[str, pd.DataFrame]], str]:
    """A ``donor_chooser`` for :func:`flowfreq.qppq.loocv_qppq`.

    Picks, for each withheld site, the most similar of the remaining sites by
    :func:`rank_donors_by_similarity`. Running ``loocv_qppq`` with and without
    it compares similarity-ranked against correlation-ranked donors on the same
    network.

    Parameters
    ----------
    site_attributes : mapping of str to mapping
        Basin characteristics for every site passed to ``loocv_qppq``.
    attributes, weights, log_attributes, scale
        As for :func:`rank_donors_by_similarity`.

    Returns
    -------
    callable
        ``(withheld_name, withheld_daily, others) -> donor_name``. It raises
        ``KeyError`` for a site with no attributes and ``ValueError`` when no
        other site can be ranked.
    """

    def _choose(name: str, _withheld: pd.DataFrame, others: Dict[str, pd.DataFrame]) -> str:
        if name not in site_attributes:
            raise KeyError(f"no basin attributes for withheld site {name!r}")
        candidates = {other: site_attributes.get(other, {}) for other in others if other != name}
        ranked = rank_donors_by_similarity(
            site_attributes[name],
            candidates,
            attributes,
            weights=weights,
            log_attributes=log_attributes,
            scale=scale,
        )
        best = ranked[ranked["status"] == STATUS_RANKED]
        if best.empty:
            raise ValueError(f"no donor for {name!r} has every attribute in {list(attributes)!r}")
        return str(best.iloc[0]["donor"])

    return _choose
