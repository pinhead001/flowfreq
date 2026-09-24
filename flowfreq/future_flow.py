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

No factor values ship with this module. Each set is transcribed from its
source with a citation, under ``flowfreq/data/future/``, as each state's wave
reaches it (issue #36).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from scipy.special import ndtri

logger = logging.getLogger(__name__)

LEGAL_STATUSES = frozenset({"required", "recommended", "informational"})
LEVELS = frozenset({"national", "state"})


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
    level : {"national", "state"}
    geography : tuple of str
        Postal codes the set applies to. Empty means national.
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
        """Whether this set covers a given state (national sets cover all)."""
        if self.level == "national":
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
    """A future-condition quantile table with its provenance."""

    table: pd.DataFrame
    factor_set: ChangeFactorSet
    overrode: Optional[str] = None
    notes: List[str] = field(default_factory=list)

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
            "overrode": self.overrode,
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
