"""Future-condition flood quantiles from published change factors.

This is the general framework from ``docs/MASTER_ROADMAP.md`` §6.3.1. A
:class:`ChangeFactorSet` is a table of per-AEP multipliers from one cited
source, for one scenario and one planning horizon. :func:`apply_change_factors`
scales a regulatory quantile table by it.

Rules this module holds to:

- The regulatory (Bulletin 17C or regression) value is never overwritten.
  The output carries both columns, labelled.
- A factor set cannot be built without a citation, scenario, and horizon.
- Factors are interpolated against the normal deviate, as
  :mod:`flowfreq.transpose` does for exponents. AEPs outside the published
  range are held at the nearest endpoint and flagged ``extrapolated``.
- State guidance overrides the national default only when it applies to
  the state in question (:func:`select_factor_set`). The override is recorded,
  and the national result is still available for comparison.

Where factor sets come from (roadmap §6.3.1, issue #36): national sources
first, then state guidance as each state's rollout wave reaches it (§6.3.2).
Both kinds live under ``flowfreq/data/future/`` as ``*.factors.json`` files,
each transcribed from a primary source with its table and page, and are
loaded by :func:`available_factor_sets`.

**No national factor set ships, and that is deliberate.** The national
sources the roadmap names were reviewed against the primary documents
(:func:`national_source_review`, ``data/future/national_sources.json``):
FHWA HEC-17 (2016) and NCHRP Project 15-61 (2019) give *methods* for deriving
site-specific factors from downscaled climate projections, plus site-specific
worked examples, but tabulate no national per-AEP flood multipliers. NOAA
Atlas 15 future precipitation is not yet published. Nothing is invented to
fill the gap; ``docs/FUTURE_FLOW_GUIDANCE.md`` records the status. A caller
with a site-specific factor set derived by one of those methods builds a
:class:`ChangeFactorSet` directly.

**No state factor set ships either.** The Wave 1 and Wave 2 states (WA, OR,
ID, MT, CO, UT, WY, NM, AZ, NV) were surveyed in 2026-10; none publishes a
tabulated per-AEP flood multiplier. Washington requires a per-site factor read
from WDFW's web tool, which is not a table; the guidance doc has each state's
findings.

Three of NCHRP 15-61's per-site procedures are mechanical once the
projection inputs exist, and are implemented here, each reproducing the
Guide's own worked example: :func:`ensemble_ratio_summary` (Ch. 6, GCM
ensemble ratios), :func:`regression_change_factors` (Ch. 7, projected inputs
in regional regression equations) and :func:`index_flood_projection` (Ch. 8,
index-flood scaling). :func:`apply_with_precedence` applies the governing set
and reports the national result alongside a state override.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from scipy.special import ndtri

from flowfreq.regression.equations import RegressionEquation, evaluate

logger = logging.getLogger(__name__)

LEGAL_STATUSES = frozenset({"required", "recommended", "informational"})
#: ``"site"`` marks a set derived for one site by a per-site procedure (for
#: example :func:`regression_change_factors` or :func:`index_flood_projection`).
#: It never takes part in national/state precedence.
LEVELS = frozenset({"national", "state", "site"})


@dataclass(frozen=True)
class ChangeFactorSet:
    """Per-AEP multipliers on flood quantiles from one source.

    Parameters
    ----------
    name : str
        Short identifier.
    aeps : sequence of float
        AEPs the source publishes factors for.
    factors : sequence of float
        Multipliers, one per AEP, all positive.
    scenario : str
        Emissions or warming scenario, e.g. "SSP5-8.5" or "+2 C".
    horizon : str
        Planning horizon, e.g. "2050" or "2041-2070".
    citation : str
        Source document and table.
    level : {"national", "state", "site"}
    geography : tuple of str
        Postal codes the set applies to. Empty means national (or, for a
        site-level set, the one site it was derived for).
    legal_status : {"required", "recommended", "informational"}
    """

    name: str
    aeps: Tuple[float, ...]
    factors: Tuple[float, ...]
    scenario: str
    horizon: str
    citation: str
    level: str = "national"
    geography: Tuple[str, ...] = ()
    legal_status: str = "informational"

    def __post_init__(self) -> None:
        for label in ("name", "scenario", "horizon", "citation"):
            if not str(getattr(self, label)).strip():
                raise ValueError(f"ChangeFactorSet requires a non-empty {label}")
        if len(self.aeps) != len(self.factors) or not self.aeps:
            raise ValueError("aeps and factors must be non-empty and equal in length")
        if any(not 0.0 < a < 1.0 for a in self.aeps):
            raise ValueError("every AEP must be in (0, 1)")
        if len(set(self.aeps)) != len(self.aeps):
            raise ValueError("duplicate AEPs")
        if any(f <= 0 for f in self.factors):
            raise ValueError("factors must be positive")
        if self.level not in LEVELS:
            raise ValueError(f"level must be one of {sorted(LEVELS)}")
        if self.legal_status not in LEGAL_STATUSES:
            raise ValueError(f"legal_status must be one of {sorted(LEGAL_STATUSES)}")
        if self.level == "state" and not self.geography:
            raise ValueError("a state-level set must name its geography")

    def applies_to(self, state: Optional[str]) -> bool:
        """Whether this set covers a given state (national and site sets cover all)."""
        if self.level != "state":
            return True
        return state is not None and state.upper() in {g.upper() for g in self.geography}

    def factor_at(self, aeps: ArrayLike) -> Tuple[np.ndarray, np.ndarray]:
        """Interpolate factors onto ``aeps``.

        Returns
        -------
        (factors, extrapolated)
            The interpolated multipliers, and a boolean mask marking the AEPs
            outside the published range, where the endpoint was held.
        """
        z_pub = ndtri(1.0 - np.asarray(self.aeps, dtype=float))
        order = np.argsort(z_pub)
        z_pub, f_pub = z_pub[order], np.asarray(self.factors, dtype=float)[order]
        z = ndtri(1.0 - np.asarray(aeps, dtype=float))
        extrap = (z < z_pub[0] - 1e-12) | (z > z_pub[-1] + 1e-12)
        return np.interp(z, z_pub, f_pub), extrap


@dataclass
class FutureQuantiles:
    """A future-condition quantile table with its provenance.

    Attributes
    ----------
    table : pandas.DataFrame
    factor_set : ChangeFactorSet
        The set that produced ``future_flow_cfs``.
    overrode : str, optional
        Name of the national set a state set replaced.
    notes : list of str
    national : FutureQuantiles, optional
        When a state set overrode a national one (:func:`apply_with_precedence`),
        the national result, reported alongside for comparison. Its numbers are
        also in ``table`` as ``national_change_factor``,
        ``national_future_flow_cfs`` and ``national_extrapolated``.
    """

    table: pd.DataFrame
    factor_set: ChangeFactorSet
    overrode: Optional[str] = None
    notes: List[str] = field(default_factory=list)
    national: Optional["FutureQuantiles"] = None

    @property
    def basis(self) -> str:
        """A label for which rule picked :attr:`factor_set`."""
        if self.overrode is not None:
            return f"state override of national default {self.overrode!r}"
        return {
            "national": "national default",
            "state": "state guidance (no national default given)",
            "site": "site-specific derivation",
        }[self.factor_set.level]

    @property
    def provenance(self) -> Dict[str, object]:
        """What produced this table, for reports."""
        fs = self.factor_set
        return {
            "factor_set": fs.name,
            "scenario": fs.scenario,
            "horizon": fs.horizon,
            "citation": fs.citation,
            "level": fs.level,
            "legal_status": fs.legal_status,
            "basis": self.basis,
            "overrode": self.overrode,
            "national_alongside": None if self.national is None else self.national.factor_set.name,
            "n_extrapolated": int(self.table["extrapolated"].sum()),
        }


def select_factor_set(
    national: ChangeFactorSet,
    state_sets: Sequence[ChangeFactorSet] = (),
    *,
    state: Optional[str] = None,
) -> Tuple[ChangeFactorSet, Optional[str]]:
    """Pick the set to apply at a site: an applicable state set wins.

    Parameters
    ----------
    national : ChangeFactorSet
        The default. Must be national-level.
    state_sets : sequence of ChangeFactorSet
        Candidate state-level sets.
    state : str, optional
        The site's postal code.

    Returns
    -------
    (ChangeFactorSet, str or None)
        The chosen set, and the name of the national set it overrode, if any.
        :func:`apply_with_precedence` applies both and labels them.

    Raises
    ------
    ValueError
        If ``national`` is not national-level, or more than one state set applies.
    """
    if national.level != "national":
        raise ValueError("the default set must be national-level")
    matches = [s for s in state_sets if s.level == "state" and s.applies_to(state)]
    if len(matches) > 1:
        raise ValueError(f"{len(matches)} state sets apply to {state}: {[s.name for s in matches]}")
    if matches:
        return matches[0], national.name
    return national, None


def apply_change_factors(
    quantiles: pd.DataFrame,
    factor_set: ChangeFactorSet,
    *,
    flow_column: str = "flow_cfs",
    overrode: Optional[str] = None,
) -> FutureQuantiles:
    """Scale a quantile table by a change-factor set.

    Parameters
    ----------
    quantiles : pandas.DataFrame
        Must have an ``aep`` column and ``flow_column`` (for example,
        ``FrequencyResults.quantiles``).
    factor_set : ChangeFactorSet
    flow_column : str, default "flow_cfs"
    overrode : str, optional
        The name of the national set this one replaced, from :func:`select_factor_set`.

    Returns
    -------
    FutureQuantiles
        The table columns are ``aep``, ``regulatory_flow_cfs``, ``change_factor``,
        ``future_flow_cfs``, and ``extrapolated``.

    Raises
    ------
    KeyError
        If a required column is missing.
    """
    for col in ("aep", flow_column):
        if col not in quantiles.columns:
            raise KeyError(f"quantiles missing column {col!r}")
    aeps = quantiles["aep"].to_numpy(dtype=float)
    factors, extrap = factor_set.factor_at(aeps)
    table = pd.DataFrame(
        {
            "aep": aeps,
            "regulatory_flow_cfs": quantiles[flow_column].to_numpy(dtype=float),
            "change_factor": factors,
        }
    )
    table["future_flow_cfs"] = table["regulatory_flow_cfs"] * table["change_factor"]
    table["extrapolated"] = extrap
    notes = []
    if extrap.any():
        msg = (
            f"{int(extrap.sum())} AEP(s) outside {factor_set.name}'s published range; "
            "endpoint factor held"
        )
        logger.warning(msg)
        notes.append(msg)
    return FutureQuantiles(table=table, factor_set=factor_set, overrode=overrode, notes=notes)


def apply_with_precedence(
    quantiles: pd.DataFrame,
    national: Optional[ChangeFactorSet],
    state_sets: Sequence[ChangeFactorSet] = (),
    *,
    state: Optional[str] = None,
    flow_column: str = "flow_cfs",
) -> FutureQuantiles:
    """Apply the governing factor set at a site, with the national result alongside.

    This is the precedence rule of roadmap §6.3.1 in one call: an applicable
    state set overrides the national default (:func:`select_factor_set`), the
    override is recorded, and the national result is still computed and
    reported next to it.

    Parameters
    ----------
    quantiles : pandas.DataFrame
        As for :func:`apply_change_factors`.
    national : ChangeFactorSet or None
        The national default. ``None`` is allowed because no national set ships
        (:func:`national_source_review`); a state set must then apply.
    state_sets : sequence of ChangeFactorSet
        Candidate state-level sets, for example
        ``available_factor_sets(level="state")``.
    state : str, optional
        The site's postal code.
    flow_column : str, default "flow_cfs"

    Returns
    -------
    FutureQuantiles
        Governed by the chosen set. When a state set overrode ``national``, the
        result's :attr:`FutureQuantiles.national` holds the national result and
        ``table`` gains the labelled columns ``national_change_factor``,
        ``national_future_flow_cfs`` and ``national_extrapolated``.
        :attr:`FutureQuantiles.basis` and ``provenance`` say which rule applied.

    Raises
    ------
    LookupError
        If ``national`` is None and no state set applies to ``state``.
    ValueError
        As :func:`select_factor_set`.
    """
    if national is None:
        matches = [s for s in state_sets if s.level == "state" and s.applies_to(state)]
        if len(matches) > 1:
            raise ValueError(
                f"{len(matches)} state sets apply to {state}: {[s.name for s in matches]}"
            )
        if not matches:
            raise LookupError(
                f"no national default was given and no state factor set applies to {state!r}"
            )
        return apply_change_factors(quantiles, matches[0], flow_column=flow_column)

    chosen, overrode = select_factor_set(national, state_sets, state=state)
    result = apply_change_factors(quantiles, chosen, flow_column=flow_column, overrode=overrode)
    if overrode is None:
        return result
    nat = apply_change_factors(quantiles, national, flow_column=flow_column)
    result.national = nat
    result.table["national_change_factor"] = nat.table["change_factor"].to_numpy()
    result.table["national_future_flow_cfs"] = nat.table["future_flow_cfs"].to_numpy()
    result.table["national_extrapolated"] = nat.table["extrapolated"].to_numpy()
    result.notes.append(
        f"{chosen.name} (state) overrides {national.name} (national); "
        "the national result is reported alongside for comparison"
    )
    return result


# --------------------------------------------------------------------------
# Shipped data: factor-set files and the national source review
# --------------------------------------------------------------------------

#: Suffix of a transcribed factor-set file under ``flowfreq/data/future/``.
FACTOR_SET_SUFFIX = ".factors.json"
_REQUIRED_KEYS = ("name", "aeps", "factors", "scenario", "horizon", "citation")


def _data_dir() -> Path:
    return Path(str(resources.files("flowfreq") / "data" / "future"))


def factor_set_from_dict(d: Dict[str, Any]) -> ChangeFactorSet:
    """Build a :class:`ChangeFactorSet` from its JSON form.

    Parameters
    ----------
    d : dict
        Keys ``name``, ``aeps``, ``factors``, ``scenario``, ``horizon`` and
        ``citation`` are required; ``level``, ``geography`` and
        ``legal_status`` are optional. Other keys (for example ``table`` or
        ``page`` notes) are ignored.

    Raises
    ------
    KeyError
        If a required key is missing.
    ValueError
        If the values fail :class:`ChangeFactorSet` validation.
    """
    missing = [k for k in _REQUIRED_KEYS if k not in d]
    if missing:
        raise KeyError(f"factor set missing key(s) {missing}")
    return ChangeFactorSet(
        name=str(d["name"]),
        aeps=tuple(float(a) for a in d["aeps"]),
        factors=tuple(float(f) for f in d["factors"]),
        scenario=str(d["scenario"]),
        horizon=str(d["horizon"]),
        citation=str(d["citation"]),
        level=str(d.get("level", "national")),
        geography=tuple(str(g) for g in d.get("geography", ())),
        legal_status=str(d.get("legal_status", "informational")),
    )


def available_factor_sets(
    directory: Optional[Union[str, Path]] = None,
    *,
    level: Optional[str] = None,
) -> Tuple[ChangeFactorSet, ...]:
    """Load every transcribed factor set shipped under ``data/future/``.

    Each ``*.factors.json`` file holds one set (an object) or several (a list).
    **None ship today** (see the module docstring), so this returns an empty
    tuple until a source with tabulated factors is transcribed.

    Parameters
    ----------
    directory : str or Path, optional
        Read from here instead of the packaged data directory.
    level : {"national", "state"}, optional
        Return only sets at this level.

    Returns
    -------
    tuple of ChangeFactorSet
        Sorted by name.

    Raises
    ------
    ValueError
        If ``level`` is not a known level, or two files define the same name.
    """
    if level is not None and level not in LEVELS:
        raise ValueError(f"level must be one of {sorted(LEVELS)}")
    root = Path(directory) if directory is not None else _data_dir()
    sets: Dict[str, ChangeFactorSet] = {}
    for path in sorted(root.glob(f"*{FACTOR_SET_SUFFIX}")):
        raw = json.loads(path.read_text(encoding="utf-8"))
        for entry in raw if isinstance(raw, list) else [raw]:
            fs = factor_set_from_dict(entry)
            if fs.name in sets:
                raise ValueError(f"duplicate factor set name {fs.name!r} in {path.name}")
            sets[fs.name] = fs
    out = sorted(sets.values(), key=lambda s: s.name)
    return tuple(s for s in out if level is None or s.level == level)


@dataclass(frozen=True)
class SourceReview:
    """What one national source was found to provide (``national_sources.json``).

    Attributes
    ----------
    id : str
    citation : str
    url : str
    provides : str
        ``"method"``, ``"factors"`` or ``"not yet published"``.
    tabulated_national_factors : bool
        Whether the source publishes a national per-AEP flow-multiplier table.
    method_summary : str
    tables_examined : tuple of dict
        Each table checked, with ``table``, ``page``, ``content`` and
        ``why_not_a_factor_set``.
    """

    id: str
    citation: str
    url: str
    provides: str
    tabulated_national_factors: bool
    method_summary: str
    tables_examined: Tuple[Dict[str, str], ...] = ()


def national_source_review(path: Optional[Union[str, Path]] = None) -> Tuple[SourceReview, ...]:
    """The review of national future-flow sources behind the empty national default.

    Parameters
    ----------
    path : str or Path, optional
        Read this file instead of the packaged ``data/future/national_sources.json``.

    Returns
    -------
    tuple of SourceReview
    """
    p = Path(path) if path is not None else _data_dir() / "national_sources.json"
    raw = json.loads(p.read_text(encoding="utf-8"))
    return tuple(
        SourceReview(
            id=s["id"],
            citation=s["citation"],
            url=s["url"],
            provides=s["provides"],
            tabulated_national_factors=bool(s["tabulated_national_factors"]),
            method_summary=s["method_summary"],
            tables_examined=tuple(dict(t) for t in s.get("tables_examined", [])),
        )
        for s in raw["sources"]
    )


# --------------------------------------------------------------------------
# NCHRP Project 15-61 per-site procedures
# --------------------------------------------------------------------------
#
# Kilgore et al. (2019), *Applying Climate Change Information to Hydrologic and
# Coastal Design of Transportation Infrastructure: Design Practices*, NCHRP
# Project 15-61 (the "Guide"). The national review found no national factors in
# it, only per-site procedures. Three of them are mechanical once the projection
# inputs are in hand, and are implemented here; each reproduces the Guide's own
# worked example (tests/test_future_flow.py). Choosing GCMs, scenarios, periods
# and grid cells, and discarding anomalous models (Guide Ch. 6, Step 8), are
# engineering judgement and stay with the caller. No climate data is fetched.

#: Citation carried into the provenance of results from the procedures below.
NCHRP_15_61 = (
    "Kilgore et al. (2019), NCHRP Project 15-61, Applying Climate Change Information to "
    "Hydrologic and Coastal Design of Transportation Infrastructure: Design Practices"
)


def ensemble_ratio_summary(
    ratios: Mapping[float, Sequence[float]],
    *,
    historical: Optional[Mapping[float, float]] = None,
    confidence: float = 0.90,
    index_aep: float = 0.1,
) -> pd.DataFrame:
    """Summarise an ensemble of future/baseline quantile ratios (Guide Ch. 6, Steps 8-9).

    Each ratio is one downscaled GCM's future quantile over its own baseline
    quantile at one grid cell (Guide Eq. 6.1). This is the mechanical part of
    Steps 8 and 9 for one grid cell, or for a watershed average the caller has
    formed: the ensemble mean, standard deviation and confidence limits of the
    ratios, and the projected quantile ``historical * ratio`` (Eq. 6.2). For
    AEPs rarer than ``index_aep`` the Guide does not recommend the models' own
    ratio and substitutes the ``index_aep`` ratio (Eq. 6.3); the
    ``applied_ratio`` column does that.

    The Guide uses this for precipitation feeding a rainfall-runoff model
    (Ch. 6), and for projecting a regression equation's precipitation variable
    (Ch. 7, Eq. 7.2). It yields precipitation, not flow; turning it into flow is
    the job of that model or of :func:`regression_change_factors`.

    Parameters
    ----------
    ratios : mapping of float to sequence of float
        AEP -> one future/baseline ratio per model, after any models the
        engineer has discarded.
    historical : mapping of float to float, optional
        AEP -> historical (observed, e.g. NOAA Atlas 14) quantile to project.
    confidence : float, default 0.90
        Two-sided level of the confidence limits.
    index_aep : float, default 0.1
        The AEP whose ratio is substituted for rarer AEPs (Eq. 6.3).

    Returns
    -------
    pandas.DataFrame
        One row per AEP, most frequent first: ``aep``, ``n_models``,
        ``mean_ratio``, ``sd_ratio``, ``lower_cl``, ``upper_cl``,
        ``substituted``, ``applied_ratio`` and, with ``historical``,
        ``historical`` and ``projected``.

    Raises
    ------
    ValueError
        If ``confidence`` is not in (0, 1), an AEP has fewer than two ratios, a
        ratio is not positive and finite, an AEP rarer than ``index_aep`` is
        given without ``index_aep`` itself, or ``historical`` lacks an AEP.

    Notes
    -----
    The standard deviation is the sample one (``ddof=1``). The Guide says the
    limits assume a normal distribution but prints no formula; the limits here
    are ``mean -/+ z * sd / sqrt(n)``, the convention that reproduces its
    Tables 6.8 and 6.9 (p. 56) from Tables 6.6 and 6.7 (p. 54). As in those
    tables, the limits are per AEP, without the Eq. 6.3 substitution.
    """
    if not 0.0 < confidence < 1.0:
        raise ValueError("confidence must be in (0, 1)")
    aeps = sorted((float(a) for a in ratios), reverse=True)
    if any(a < index_aep for a in aeps) and not any(np.isclose(a, index_aep) for a in aeps):
        raise ValueError(
            f"AEPs rarer than {index_aep} take the {index_aep}-AEP ratio (Guide Eq. 6.3), "
            "which was not given"
        )
    z = float(ndtri(0.5 + confidence / 2.0))
    by_aep = {float(a): v for a, v in ratios.items()}
    rows = []
    for a in aeps:
        r = np.asarray(by_aep[a], dtype=float)
        if r.size < 2:
            raise ValueError(f"AEP {a}: need at least two model ratios")
        if np.any(~np.isfinite(r)) or np.any(r <= 0):
            raise ValueError(f"AEP {a}: ratios must be positive and finite")
        mean, sd = float(r.mean()), float(r.std(ddof=1))
        half = z * sd / float(np.sqrt(r.size))
        rows.append(
            {
                "aep": a,
                "n_models": int(r.size),
                "mean_ratio": mean,
                "sd_ratio": sd,
                "lower_cl": mean - half,
                "upper_cl": mean + half,
            }
        )
    out = pd.DataFrame(rows)
    index_rows = [row["mean_ratio"] for row in rows if np.isclose(row["aep"], index_aep)]
    index_mean = index_rows[0] if index_rows else float("nan")
    out["substituted"] = out["aep"] < index_aep - 1e-12
    out["applied_ratio"] = np.where(out["substituted"], index_mean, out["mean_ratio"])
    if historical is not None:
        hist = {float(k): float(v) for k, v in historical.items()}
        missing = [a for a in aeps if a not in hist]
        if missing:
            raise ValueError(f"historical quantiles missing for AEP(s) {missing}")
        out["historical"] = [hist[a] for a in aeps]
        out["projected"] = out["historical"] * out["applied_ratio"]
    return out


@dataclass
class RegressionProjection:
    """Regional-regression flows under historical and projected inputs (Guide Ch. 7).

    Attributes
    ----------
    table : pandas.DataFrame
        ``aep``, ``historical_flow_cfs``, ``projected_flow_cfs`` and
        ``change_factor`` (projected / historical), most frequent AEP first.
    factor_set : ChangeFactorSet
        The per-AEP ratios as a site-level set, ready for :func:`apply_change_factors`.
    changed : dict
        Each input that differs between the two runs, as ``(historical, projected)``.
    """

    table: pd.DataFrame
    factor_set: ChangeFactorSet
    changed: Dict[str, Tuple[float, float]]


def regression_change_factors(
    equations: Sequence[RegressionEquation],
    historical: Mapping[str, float],
    projected: Mapping[str, float],
    *,
    name: str,
    scenario: str,
    horizon: str,
    citation: str,
) -> RegressionProjection:
    """Project flows by substituting projected inputs into regression equations (Guide Ch. 7).

    The Guide's four steps (Section 7.1, pp. 77-79): (1) use an equation with a
    precipitation (or other projectable climate) variable; (2-3) confirm every
    input, historical and projected, lies within the equation's calibrated
    range; (4) evaluate it with both, and compare. Step 3's projected value is
    usually ``historical * (future / baseline)`` from downscaled GCM output
    (Eq. 7.2), whose ensemble :func:`ensemble_ratio_summary` summarises.

    Parameters
    ----------
    equations : sequence of RegressionEquation
        One per AEP, for one region, as from
        :meth:`flowfreq.regression.library.StateLibrary.equation`.
    historical, projected : mapping of str to float
        Basin characteristics by code. ``projected`` need list only the inputs
        that change; the rest are taken from ``historical``.
    name, scenario, horizon, citation : str
        Provenance of the projection, carried into :attr:`RegressionProjection.factor_set`.

    Returns
    -------
    RegressionProjection

    Raises
    ------
    ValueError
        If ``equations`` is empty or repeats an AEP, if ``projected`` changes
        nothing, or if an equation uses none of the changed inputs (the Guide,
        Section 7.3: such an equation "is not equipped to estimate future
        conditions").
    KeyError
        If an equation input, or a projected input's historical value, is missing.
    flowfreq.regression.equations.OutOfRangeError
        If any input, historical or projected, is outside the calibrated range.
        There is no override: the Guide makes staying in range a condition of
        the method (Section 7.3, "prudent limits").
    """
    if not equations:
        raise ValueError("no equations given")
    aeps = [float(eq.aep) for eq in equations]
    if len(set(aeps)) != len(aeps):
        raise ValueError(f"duplicate AEPs in equations: {aeps}")
    missing_hist = [k for k in projected if k not in historical]
    if missing_hist:
        raise KeyError(f"projected inputs {missing_hist} have no historical value")
    changed = {
        k: (float(historical[k]), float(v))
        for k, v in projected.items()
        if float(historical[k]) != float(v)
    }
    if not changed:
        raise ValueError("projected inputs equal the historical ones; nothing to project")
    future = {
        **{k: float(v) for k, v in historical.items()},
        **{k: float(v) for k, v in projected.items()},
    }
    rows = []
    for eq in sorted(equations, key=lambda e: -e.aep):
        if not {v.code for v in eq.variables} & set(changed):
            raise ValueError(
                f"{eq.region_code} AEP {eq.aep} uses none of the projected inputs "
                f"{sorted(changed)}, so cannot estimate future conditions (NCHRP 15-61 §7.3)"
            )
        q_hist = evaluate(eq, historical).flow_cfs
        q_proj = evaluate(eq, future).flow_cfs
        rows.append(
            {
                "aep": float(eq.aep),
                "historical_flow_cfs": q_hist,
                "projected_flow_cfs": q_proj,
                "change_factor": q_proj / q_hist,
            }
        )
    table = pd.DataFrame(rows)
    fs = ChangeFactorSet(
        name=name,
        aeps=tuple(float(a) for a in table["aep"]),
        factors=tuple(float(f) for f in table["change_factor"]),
        scenario=scenario,
        horizon=horizon,
        citation=f"{citation}; method: {NCHRP_15_61}, Ch. 7",
        level="site",
    )
    return RegressionProjection(table=table, factor_set=fs, changed=changed)


def index_flood_projection(
    quantiles: pd.DataFrame,
    *,
    historical_index_estimate: float,
    projected_index_estimate: float,
    name: str,
    scenario: str,
    horizon: str,
    citation: str,
    index_aep: float = 0.1,
    flow_column: str = "flow_cfs",
) -> FutureQuantiles:
    """Project a flood-frequency curve by the index-flood method (Guide Ch. 8).

    The Guide's four steps (Section 8.1, pp. 85-86): (1) a historical curve by
    any method; (2) flood ratios, each quantile over the ``index_aep``
    (10-year) flood; (3) a projected index flood; (4) the projected curve,
    flood ratio times projected index flood. When the projected index flood
    comes from a different method than the curve (for example a regression
    equation, Ch. 7), Step 3 scales the curve's own index flood by that
    method's projected/historical ratio. That is the form taken here: pass the
    projecting method's historical and projected estimates of the index flood.
    If the curve itself came from the projecting method, pass the curve's own
    index flow as ``historical_index_estimate``.

    Parameters
    ----------
    quantiles : pandas.DataFrame
        The historical curve: ``aep`` and ``flow_column``, including ``index_aep``.
    historical_index_estimate, projected_index_estimate : float
        The projecting method's index flood under historical and projected inputs.
    name, scenario, horizon, citation : str
        Provenance of the projection.
    index_aep : float, default 0.1
        The Guide uses the 0.1 AEP (Section 8.1).
    flow_column : str, default "flow_cfs"

    Returns
    -------
    FutureQuantiles
        ``future_flow_cfs`` is the projected curve, and ``table`` also carries
        ``flood_ratio`` (Step 2). The factor set is site-level, with the same
        ratio at every AEP of the curve.

    Raises
    ------
    KeyError
        If a required column is missing.
    ValueError
        If ``index_aep`` is not on the curve exactly once, or an estimate is not positive.
    """
    for col in ("aep", flow_column):
        if col not in quantiles.columns:
            raise KeyError(f"quantiles missing column {col!r}")
    if not (historical_index_estimate > 0 and projected_index_estimate > 0):
        raise ValueError("index-flood estimates must be positive")
    aeps = quantiles["aep"].to_numpy(dtype=float)
    flows = quantiles[flow_column].to_numpy(dtype=float)
    at_index = np.isclose(aeps, index_aep)
    if int(at_index.sum()) != 1:
        raise ValueError(f"the curve must contain the index AEP {index_aep} exactly once")
    index_flow = float(flows[at_index][0])
    ratio = projected_index_estimate / historical_index_estimate
    fs = ChangeFactorSet(
        name=name,
        aeps=tuple(float(a) for a in aeps),
        factors=tuple([ratio] * len(aeps)),
        scenario=scenario,
        horizon=horizon,
        citation=f"{citation}; method: {NCHRP_15_61}, Ch. 8",
        level="site",
    )
    result = apply_change_factors(quantiles, fs, flow_column=flow_column)
    result.table.insert(2, "flood_ratio", flows / index_flow)
    result.notes.append(
        f"index flood ({index_aep} AEP) {index_flow:g} cfs projected to "
        f"{index_flow * ratio:g} cfs (ratio {ratio:.4f})"
    )
    return result
