"""Systematic differences between neighboring states' equations near their shared border.

Two states' regression equations rarely agree on a basin at their border: each was fitted to
its own gages. ``tools/border_qa.py`` measured the difference on every border. At each gage
within :data:`BORDER_KM` of a neighboring state it delineated the gage's basin, evaluated both
states' equations on it, and compared them at the 1 percent AEP. A pair is **systematic** when
3 or more clean comparisons (neither side extrapolated, no urban set picked by default) have a
median ratio beyond 1.26 (0.1 in log10). ``docs/BORDER_QA.md`` has every pair and
``docs/VERIFICATION_PLAN.md`` the reading of them.

This module ships the systematic pairs (``flowfreq/data/border_offsets.csv``, written by
``python tools/border_qa.py offsets``) so that an estimate near such a border can say so. It
holds no geometry: how far a point is from a state line is the caller's to measure.

Nothing here adjusts an estimate. A difference is reported, not smoothed.

A pair is directional. ``home`` is the state the gages are in, and the basins compared are
the home state's. An offset measured at Arkansas gages near Missouri is one row; Missouri
gages near Arkansas are another, which may not be systematic or may not have been measured.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from flowfreq.regression.jurisdictions import jurisdiction

#: The distance from the neighboring state's outline, in km, within which gages were compared.
BORDER_KM = 20.0

OFFSETS_TABLE = Path(__file__).parent.parent / "data" / "border_offsets.csv"


@dataclass(frozen=True)
class BorderOffset:
    """One systematic difference: the neighbor's equations on the home state's border basins.

    Attributes
    ----------
    home : str
        Postal code of the state the compared gages are in.
    neighbor : str
        Postal code of the adjacent state whose equations were evaluated on the same basins.
    gages : int
        Number of clean comparisons.
    neighbor_over_home : float
        Median of the neighbor's 1 percent AEP flow over the home state's.
    with_fit : int
        How many of those gages have an at-site Bulletin 17C fit (station skew) to compare with.
    home_over_at_site, neighbor_over_at_site : float or None
        Median of each state's 1 percent AEP flow over the at-site fit; ``None`` when no gage
        has a fit.
    home_nearer : int
        At how many of the ``with_fit`` gages the home state's value was nearer the fit.
    """

    home: str
    neighbor: str
    gages: int
    neighbor_over_home: float
    with_fit: int
    home_over_at_site: Optional[float]
    neighbor_over_at_site: Optional[float]
    home_nearer: int

    def describe(self) -> str:
        """The difference in a sentence or two, naming both states.

        Returns
        -------
        str
            For example "At 11 Arkansas gages within 20 km of Missouri, the Missouri
            equations gave a median 0.73 times the Arkansas equations at the 1 percent AEP.
            The Arkansas equations were nearer the at-site Bulletin 17C fit at 9 of 11 gages."
        """
        home, neighbor = jurisdiction(self.home).name, jurisdiction(self.neighbor).name
        text = (
            f"At {self.gages} {home} gages within {BORDER_KM:g} km of {neighbor}, "
            f"the {neighbor} equations gave a median {self.neighbor_over_home:g} times "
            f"the {home} equations at the 1 percent AEP."
        )
        if self.with_fit:
            of = f"{self.with_fit} gages" + ("" if self.with_fit == self.gages else " with a fit")
            text += (
                f" The {home} equations were nearer the at-site Bulletin 17C fit at "
                f"{self.home_nearer} of {of}."
            )
        return text


def _optional(value: Any) -> Optional[float]:
    return None if pd.isna(value) else float(value)


@lru_cache(maxsize=None)
def _load(path: Path) -> Tuple[BorderOffset, ...]:
    # Only an empty cell is missing: the at-site medians of a pair no gage of which has a fit.
    table = pd.read_csv(
        path, dtype={"home": str, "neighbor": str}, keep_default_na=False, na_values=[""]
    )
    rows: List[Dict[str, Any]] = [
        {str(k): v for k, v in r.items()} for r in table.to_dict("records")
    ]
    return tuple(
        BorderOffset(
            home=str(row["home"]),
            neighbor=str(row["neighbor"]),
            gages=int(row["gages"]),
            neighbor_over_home=float(row["neighbor_over_home"]),
            with_fit=int(row["with_fit"]),
            home_over_at_site=_optional(row["home_over_at_site"]),
            neighbor_over_at_site=_optional(row["neighbor_over_at_site"]),
            home_nearer=int(row["home_nearer"]),
        )
        for row in rows
    )


def border_offsets(path: Optional[Path] = None) -> Tuple[BorderOffset, ...]:
    """Every systematic border difference, ordered by home state and neighbor.

    Parameters
    ----------
    path : Path, optional
        A table in the layout of ``flowfreq/data/border_offsets.csv``. Defaults to the
        shipped table.

    Returns
    -------
    tuple of BorderOffset

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    """
    return _load(Path(path) if path is not None else OFFSETS_TABLE)


def border_offsets_between(
    state: str, other: str, path: Optional[Path] = None
) -> Tuple[BorderOffset, ...]:
    """The systematic differences on the border between two states, in either direction.

    Parameters
    ----------
    state, other : str
        Two-letter postal codes, case-insensitive.
    path : Path, optional
        As for :func:`border_offsets`.

    Returns
    -------
    tuple of BorderOffset
        Up to two: the one measured at ``state``'s gages first, then the one measured at
        ``other``'s. Empty when neither direction is systematic, which also covers a pair that
        was never compared (states that are not adjacent, or without verified equations).

    Raises
    ------
    KeyError
        If either code is not a jurisdiction.
    """
    a, b = jurisdiction(state).code, jurisdiction(other).code
    found = {(o.home, o.neighbor): o for o in border_offsets(path)}
    return tuple(found[key] for key in ((a, b), (b, a)) if key in found)


def border_offsets_for(state: str, path: Optional[Path] = None) -> Tuple[BorderOffset, ...]:
    """The systematic differences on any of a state's borders, in either direction.

    Parameters
    ----------
    state : str
        Two-letter postal code, case-insensitive.
    path : Path, optional
        As for :func:`border_offsets`.

    Returns
    -------
    tuple of BorderOffset
        Those measured at ``state``'s gages first, then those measured at a neighbor's.

    Raises
    ------
    KeyError
        If the code is not a jurisdiction.
    """
    code = jurisdiction(state).code
    offsets = border_offsets(path)
    return tuple(o for o in offsets if o.home == code) + tuple(
        o for o in offsets if o.neighbor == code
    )
