"""The 56 jurisdictions and the rollout wave each belongs to.

This is the single source of truth for ``docs/MASTER_ROADMAP.md`` §3.2's wave
table. ``tools/gen_regression_coverage.py`` reads it, so the table in the
roadmap and the generated coverage report cannot drift apart.

Waves follow river basins because regional skew studies, regression
calibration datasets, and donor gages cross state lines along them.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple


@dataclass(frozen=True)
class Jurisdiction:
    """A state, district, or territory."""

    code: str
    name: str
    wave: str
    is_territory: bool = False


WAVE_NAMES: Dict[str, str] = {
    "1": "Columbia River basin",
    "2": "Colorado River basin",
    "3": "West coast and Pacific",
    "4a": "Upper Mississippi and Great Lakes",
    "4b": "Missouri basin",
    "4c": "Ohio basin",
    "4d": "Lower Mississippi, Arkansas-Red, and Gulf",
    "5": "East coast",
    "6": "Territories",
}

# Order within a wave is the planned order of work.
_TABLE: Tuple[Tuple[str, str, str], ...] = (
    ("WA", "Washington", "1"),
    ("OR", "Oregon", "1"),
    ("ID", "Idaho", "1"),
    ("MT", "Montana", "1"),
    ("CO", "Colorado", "2"),
    ("UT", "Utah", "2"),
    ("WY", "Wyoming", "2"),
    ("NM", "New Mexico", "2"),
    ("AZ", "Arizona", "2"),
    ("NV", "Nevada", "2"),
    ("CA", "California", "3"),
    ("AK", "Alaska", "3"),
    ("HI", "Hawaii", "3"),
    ("MN", "Minnesota", "4a"),
    ("WI", "Wisconsin", "4a"),
    ("IA", "Iowa", "4a"),
    ("IL", "Illinois", "4a"),
    ("MO", "Missouri", "4a"),
    ("MI", "Michigan", "4a"),
    ("ND", "North Dakota", "4b"),
    ("SD", "South Dakota", "4b"),
    ("NE", "Nebraska", "4b"),
    ("KS", "Kansas", "4b"),
    ("OH", "Ohio", "4c"),
    ("IN", "Indiana", "4c"),
    ("KY", "Kentucky", "4c"),
    ("TN", "Tennessee", "4c"),
    ("WV", "West Virginia", "4c"),
    ("AR", "Arkansas", "4d"),
    ("LA", "Louisiana", "4d"),
    ("MS", "Mississippi", "4d"),
    ("OK", "Oklahoma", "4d"),
    ("TX", "Texas", "4d"),
    ("ME", "Maine", "5"),
    ("NH", "New Hampshire", "5"),
    ("VT", "Vermont", "5"),
    ("MA", "Massachusetts", "5"),
    ("RI", "Rhode Island", "5"),
    ("CT", "Connecticut", "5"),
    ("NY", "New York", "5"),
    ("NJ", "New Jersey", "5"),
    ("PA", "Pennsylvania", "5"),
    ("DE", "Delaware", "5"),
    ("MD", "Maryland", "5"),
    ("DC", "District of Columbia", "5"),
    ("VA", "Virginia", "5"),
    ("NC", "North Carolina", "5"),
    ("SC", "South Carolina", "5"),
    ("GA", "Georgia", "5"),
    ("FL", "Florida", "5"),
    ("AL", "Alabama", "5"),
    ("PR", "Puerto Rico", "6"),
    ("VI", "U.S. Virgin Islands", "6"),
    ("GU", "Guam", "6"),
    ("MP", "Northern Mariana Islands", "6"),
    ("AS", "American Samoa", "6"),
)

JURISDICTIONS: Tuple[Jurisdiction, ...] = tuple(
    Jurisdiction(code, name, wave, is_territory=(wave == "6")) for code, name, wave in _TABLE
)
BY_CODE: Dict[str, Jurisdiction] = {j.code: j for j in JURISDICTIONS}


def jurisdiction(code: str) -> Jurisdiction:
    """Look up a jurisdiction by postal code.

    Parameters
    ----------
    code : str
        Two-letter postal code, case-insensitive.

    Returns
    -------
    Jurisdiction

    Raises
    ------
    KeyError
        For an unknown code.
    """
    try:
        return BY_CODE[code.upper()]
    except KeyError:
        raise KeyError(f"Unknown jurisdiction code {code!r}") from None


def wave(wave_id: str) -> List[Jurisdiction]:
    """Jurisdictions in a wave, in planned order of work.

    Parameters
    ----------
    wave_id : str
        A key of :data:`WAVE_NAMES`.

    Returns
    -------
    list of Jurisdiction

    Raises
    ------
    KeyError
        For an unknown wave.
    """
    if wave_id not in WAVE_NAMES:
        raise KeyError(f"Unknown wave {wave_id!r}; choose from {list(WAVE_NAMES)}")
    return [j for j in JURISDICTIONS if j.wave == wave_id]
