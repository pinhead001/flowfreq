"""Offline library of published regression equations, one JSON file per state.

Files live in ``flowfreq/data/regression/<CODE>.json``. A file with
``status: "pending"`` records that the state is on the plan but nothing has
been transcribed yet. It loads, and it carries no equations. Asking it for an
estimate raises :class:`EquationsUnavailable`, never a default value.

Populate the files from the cited report or from ``tools/snapshot_nss.py``,
never by hand from memory.

Equation sets flowfreq develops itself (for example a Bulletin 17C-era GLS
update of a state's B17B-era equations) are *not* published USGS equations.
They live in ``flowfreq/data/regression/provisional/<CODE>.json`` with
``status: "provisional"``, and :func:`load_state` returns one only when the
caller passes ``allow_provisional=True``. The default lookup is never changed by
a provisional file, mirroring :mod:`flowfreq.regional_skew`.

Roadmap: ``docs/MASTER_ROADMAP.md`` §3.1-3.2, issues #35 and #37-#40.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from flowfreq.regression.equations import RegressionEquation
from flowfreq.regression.jurisdictions import jurisdiction

logger = logging.getLogger(__name__)

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "regression"
#: Subdirectory of the data directory holding flowfreq-developed equation sets.
PROVISIONAL_SUBDIR = "provisional"
SCHEMA_VERSION = 1
#: How faithfully a state file transcribes its source reports -- a statement about
#: the transcription, not about how much the reports themselves publish:
#:
#: - ``pending``: no equations transcribed yet.
#: - ``partial``: something the cited reports *do* publish is missing or unresolved
#:   -- a region or AEP not stored (e.g. a form the schema cannot express), a
#:   published covariance or other statistic deliberately left out, or a
#:   cross-check that disagrees with the report or NSS beyond rounding.
#: - ``verified``: everything the cited reports publish for the stored regions is
#:   stored, double-entered against the source, and cross-checked -- against live
#:   NSS where NSS serves the region, and otherwise against the report's own
#:   worked examples and per-gage estimates. A report that publishes no covariance
#:   (or no AVP) can still be ``verified``; intervals then fall back to ``sep_log``.
#: - ``provisional``: developed by flowfreq following published USGS procedure, not
#:   a USGS publication. Allowed only under ``provisional/``, and returned only on
#:   opt-in (``load_state(code, allow_provisional=True)``).
STATUSES = frozenset({"pending", "partial", "verified", "provisional"})


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
    advisory: str = ""

    @property
    def provisional(self) -> bool:
        """True for a flowfreq-developed (non-USGS) equation set."""
        return self.status == "provisional"

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


def _provisional_advisory(code: str) -> str:
    return (
        f"PROVISIONAL (not USGS): the {code} equations were developed by flowfreq following "
        "published USGS procedure and have not been reviewed or published by the USGS. The "
        f"published set is load_state({code!r}). A USGS report supersedes this one."
    )


def load_state(
    code: str,
    data_dir: Union[str, Path, None] = None,
    *,
    allow_provisional: bool = False,
) -> StateLibrary:
    """Load one jurisdiction's equation library.

    Parameters
    ----------
    code : str
        Two-letter postal code.
    data_dir : str or Path, optional
        Override the packaged data directory.
    allow_provisional : bool
        Return the jurisdiction's provisional (flowfreq-developed, non-USGS) set
        from ``provisional/<CODE>.json`` when there is one, with ``advisory`` set
        and a logged warning. Without one, the published set is returned as usual.

    Returns
    -------
    StateLibrary

    Raises
    ------
    EquationsUnavailable
        If the jurisdiction has no file yet (its wave has not started).
    ValueError
        If a published file claims ``provisional`` status, or a file under
        ``provisional/`` claims any other status or another jurisdiction.
    """
    code = jurisdiction(code).code
    root = Path(data_dir or DATA_DIR)
    if allow_provisional:
        prov = root / PROVISIONAL_SUBDIR / f"{code}.json"
        if prov.exists():
            lib = parse_state(json.loads(prov.read_text(encoding="utf-8")))
            if lib.status != "provisional" or lib.state != code:
                raise ValueError(
                    f"{prov}: a provisional set must have status 'provisional' and state "
                    f"{code!r}, got {lib.status!r} / {lib.state!r}"
                )
            lib.advisory = _provisional_advisory(code)
            logger.warning("%s", lib.advisory)
            return lib
    path = root / f"{code}.json"
    if not path.exists():
        raise EquationsUnavailable(f"{code}: no regression library file; its wave has not started")
    lib = parse_state(json.loads(path.read_text(encoding="utf-8")))
    if lib.provisional:
        raise ValueError(
            f"{path}: status 'provisional' belongs under {PROVISIONAL_SUBDIR}/, not in place of "
            "the published set"
        )
    return lib


def available_states(data_dir: Union[str, Path, None] = None) -> Dict[str, str]:
    """Map every jurisdiction with a published library file to its status.

    Provisional sets are not listed; see :func:`provisional_states`.
    """
    out: Dict[str, str] = {}
    for p in sorted(Path(data_dir or DATA_DIR).glob("*.json")):
        out[p.stem] = parse_state(json.loads(p.read_text(encoding="utf-8"))).status
    return out


def provisional_states(data_dir: Union[str, Path, None] = None) -> List[str]:
    """Jurisdictions with a provisional (flowfreq-developed) equation set, sorted."""
    return sorted(p.stem for p in (Path(data_dir or DATA_DIR) / PROVISIONAL_SUBDIR).glob("*.json"))
