"""Offline library of published regression equations, one JSON file per state.

Files live in ``flowfreq/data/regression/<CODE>.json``. A file with
``status: "pending"`` records that the state is on the plan but nothing has
been transcribed yet. It loads, and it carries no equations. Asking it for an
estimate raises :class:`EquationsUnavailable`, never a default value.

Populate the files from the cited report or from ``tools/snapshot_nss.py``,
never by hand from memory.

Roadmap: ``docs/MASTER_ROADMAP.md`` §3.1-3.2, issues #35 and #37-#40.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from flowfreq.regression.equations import RegressionEquation
from flowfreq.regression.jurisdictions import jurisdiction

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "regression"
SCHEMA_VERSION = 1
STATUSES = frozenset({"pending", "partial", "verified"})


class EquationsUnavailable(LookupError):
    """No equations are available for the requested state, region, or AEP."""


@dataclass
class StateLibrary:
    """The regression equations published for one jurisdiction."""

    state: str
    status: str
    equations: List[RegressionEquation] = field(default_factory=list)
    source_reports: List[Dict[str, Any]] = field(default_factory=list)
    issue: Optional[int] = None
    notes: str = ""

    @property
    def regions(self) -> List[str]:
        """Region codes with at least one equation, sorted."""
        return sorted({e.region_code for e in self.equations})

    def aeps(self, region_code: str) -> List[float]:
        """AEPs published for a region, largest first."""
        return sorted({e.aep for e in self.equations if e.region_code == region_code}, reverse=True)

    def equation(self, region_code: str, aep: float) -> RegressionEquation:
        """Return the equation for one region and AEP.

        Raises
        ------
        EquationsUnavailable
            If the state is pending or the pair is not published.
        """
        if not self.equations:
            raise EquationsUnavailable(
                f"{self.state}: no equations transcribed (status {self.status!r}"
                + (f", issue #{self.issue}" if self.issue else "")
                + ")"
            )
        for e in self.equations:
            if e.region_code == region_code and abs(e.aep - aep) < 1e-12:
                return e
        raise EquationsUnavailable(
            f"{self.state}: no equation for region {region_code!r} at AEP {aep}"
        )


def parse_state(data: Dict[str, Any]) -> StateLibrary:
    """Validate and build a :class:`StateLibrary` from its JSON form.

    Raises
    ------
    ValueError
        On a schema-version mismatch, unknown status or jurisdiction, a
        ``verified`` file with no equations, or duplicate (region, AEP) pairs.
    """
    if data.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"schema_version must be {SCHEMA_VERSION}, got {data.get('schema_version')}"
        )
    state = jurisdiction(str(data["state"])).code
    status = str(data["status"])
    if status not in STATUSES:
        raise ValueError(f"{state}: unknown status {status!r}")
    eqs = [RegressionEquation.from_dict(e) for e in data.get("equations", [])]
    if status == "verified" and not eqs:
        raise ValueError(f"{state}: status 'verified' with no equations")
    keys = [(e.region_code, e.aep) for e in eqs]
    if len(set(keys)) != len(keys):
        raise ValueError(f"{state}: duplicate (region, AEP) equations")
    return StateLibrary(
        state=state,
        status=status,
        equations=eqs,
        source_reports=list(data.get("source_reports", [])),
        issue=data.get("issue"),
        notes=str(data.get("notes", "")),
    )


def load_state(code: str, data_dir: Union[str, Path, None] = None) -> StateLibrary:
    """Load one jurisdiction's equation library.

    Parameters
    ----------
    code : str
        Two-letter postal code.
    data_dir : str or Path, optional
        Override the packaged data directory.

    Returns
    -------
    StateLibrary

    Raises
    ------
    EquationsUnavailable
        If the jurisdiction has no file yet (its wave has not started).
    """
    code = jurisdiction(code).code
    path = Path(data_dir or DATA_DIR) / f"{code}.json"
    if not path.exists():
        raise EquationsUnavailable(f"{code}: no regression library file; its wave has not started")
    return parse_state(json.loads(path.read_text(encoding="utf-8")))


def available_states(data_dir: Union[str, Path, None] = None) -> Dict[str, str]:
    """Map every jurisdiction with a library file to its status."""
    out: Dict[str, str] = {}
    for p in sorted(Path(data_dir or DATA_DIR).glob("*.json")):
        out[p.stem] = parse_state(json.loads(p.read_text(encoding="utf-8"))).status
    return out
