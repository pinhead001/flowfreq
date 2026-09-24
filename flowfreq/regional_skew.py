"""Regional (generalized) skew lookup for Bulletin 17C.

Bulletin 17C weights the at-site skew with a regional skew and its MSE. Until
now flowfreq took that pair as two bare floats from the caller. This module
gives them a versioned table, with a citation for every value, and a lookup
that raises rather than guessing.

Two rules the table enforces:

- A row is usable only when ``status == "verified"``, meaning it has been
  transcribed from its cited report and second-checked. ``pending`` rows exist
  so that the gaps are visible, and they are never returned as values.
- There is no silent fallback to the Bulletin 17B Plate I map. A caller who
  wants a fallback must pass one explicitly, and it is labelled as such.

Roadmap: ``docs/MASTER_ROADMAP.md`` §1.3, issue #34.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union

import pandas as pd

logger = logging.getLogger(__name__)

DEFAULT_TABLE = Path(__file__).parent / "data" / "regional_skew.csv"
COLUMNS = (
    "state",
    "skew_region",
    "skew",
    "skew_mse",
    "effective_record_years",
    "citation",
    "table",
    "status",
    "notes",
)
STATUSES = frozenset({"pending", "verified"})


class RegionalSkewUnavailable(LookupError):
    """No verified regional skew for the requested state or region."""


@dataclass(frozen=True)
class RegionalSkew:
    """A regional skew with its uncertainty and source."""

    state: str
    skew_region: str
    skew: float
    skew_mse: float
    citation: str
    table: str = ""
    effective_record_years: Optional[float] = None

    def __post_init__(self) -> None:
        if not self.citation.strip():
            raise ValueError("RegionalSkew requires a citation")
        if not math.isfinite(self.skew):
            raise ValueError(f"Non-finite skew {self.skew}")
        if not (self.skew_mse > 0 and math.isfinite(self.skew_mse)):
            raise ValueError(f"skew_mse must be positive and finite, got {self.skew_mse}")


def load_table(path: Union[str, Path, None] = None) -> pd.DataFrame:
    """Load and validate the regional skew table.

    Parameters
    ----------
    path : str or Path, optional
        Defaults to the packaged ``flowfreq/data/regional_skew.csv``.

    Returns
    -------
    pandas.DataFrame

    Raises
    ------
    ValueError
        On missing columns, an unknown status, a duplicate (state, region)
        pair, or a verified row missing a value or citation.
    """
    df = pd.read_csv(path or DEFAULT_TABLE, dtype={"state": str, "skew_region": str, "table": str})
    missing = [c for c in COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"regional skew table missing columns: {missing}")
    bad = set(df["status"]) - STATUSES
    if bad:
        raise ValueError(f"unknown status value(s): {sorted(bad)}")
    dup = df.duplicated(subset=["state", "skew_region"])
    if dup.any():
        raise ValueError(
            f"duplicate (state, skew_region): {df.loc[dup, ['state', 'skew_region']].values.tolist()}"
        )
    verified = df[df["status"] == "verified"]
    for col in ("skew", "skew_mse", "citation"):
        if verified[col].isna().any():
            raise ValueError(f"verified row(s) missing {col}")
    return df


def available(state: str, path: Union[str, Path, None] = None) -> List[str]:
    """List the verified skew regions for a state.

    Parameters
    ----------
    state : str
        Two-letter postal code.
    path : str or Path, optional
        Table override.

    Returns
    -------
    list of str
    """
    df = load_table(path)
    rows = df[(df["state"] == state.upper()) & (df["status"] == "verified")]
    return sorted(rows["skew_region"].tolist())


def regional_skew_for(
    state: str,
    skew_region: Optional[str] = None,
    *,
    path: Union[str, Path, None] = None,
) -> RegionalSkew:
    """Return the verified regional skew for a state and skew region.

    Parameters
    ----------
    state : str
        Two-letter postal code.
    skew_region : str, optional
        Needed when the state has more than one verified skew region.
    path : str or Path, optional
        Table override.

    Returns
    -------
    RegionalSkew

    Raises
    ------
    RegionalSkewUnavailable
        If there is no verified row, or the region is ambiguous.
    """
    df = load_table(path)
    st = state.upper()
    rows = df[df["state"] == st]
    if skew_region is not None:
        rows = rows[rows["skew_region"] == skew_region]
    verified = rows[rows["status"] == "verified"]
    if verified.empty:
        pending = rows[rows["status"] == "pending"]
        why = "pending transcription" if not pending.empty else "no entry"
        raise RegionalSkewUnavailable(
            f"No verified regional skew for {st}"
            + (f" region {skew_region!r}" if skew_region else "")
            + f" ({why}). Supply regional_skew/regional_skew_mse explicitly."
        )
    if len(verified) > 1:
        raise RegionalSkewUnavailable(
            f"{st} has {len(verified)} skew regions {sorted(verified['skew_region'])}; name one"
        )
    r = verified.iloc[0]
    erl = r["effective_record_years"]
    return RegionalSkew(
        state=st,
        skew_region=str(r["skew_region"]),
        skew=float(r["skew"]),
        skew_mse=float(r["skew_mse"]),
        citation=str(r["citation"]),
        table="" if pd.isna(r["table"]) else str(r["table"]),
        effective_record_years=None if pd.isna(erl) else float(erl),
    )
