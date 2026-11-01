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

USGS-published equations that NSS/StreamStats does not serve yet (a newer report
than the one NSS uses, or one NSS holds only as a staging draft) live in
``flowfreq/data/regression/awaiting_nss/<CODE>.json`` with
``status: "awaiting_nss"``, and :func:`load_state` returns one only when the caller
passes ``allow_awaiting_nss=True``. They get their own folder rather than sharing
``provisional/``: the loader keys each opt-in folder to exactly one status, the two
are different claims (USGS-published versus flowfreq-developed), and a state can
have both at once (one ``<CODE>.json`` per folder).

Equations harvested from the USGS NSS web service and parsed by flowfreq, but not
checked against the source report, live in
``flowfreq/data/regression/nss_snapshot/<CODE>.json`` with ``status:
"nss_snapshot"``, and :func:`load_state` returns one only when the caller passes
``allow_nss_snapshot=True``. :func:`flowfreq.regression.nss_library.snapshot_to_library`
generates them from ``tools/snapshot_nss.py`` output.

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
#: Subdirectory holding USGS-published sets that NSS/StreamStats does not serve yet.
AWAITING_NSS_SUBDIR = "awaiting_nss"
#: Subdirectory holding sets harvested from NSS and not checked against the report.
NSS_SNAPSHOT_SUBDIR = "nss_snapshot"
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
#: - ``awaiting_nss``: published by the USGS and transcribed to the ``verified``
#:   standard (double-entered against the report and cross-checked), but not yet
#:   served by NSS/StreamStats, which still serves the set in ``<CODE>.json``.
#:   Allowed only under ``awaiting_nss/``, and returned only on opt-in
#:   (``load_state(code, allow_awaiting_nss=True)``).
#: - ``nss_snapshot``: harvested from the USGS NSS web service and parsed by
#:   flowfreq, but not checked against the source report ("available with
#:   caveat", never ``verified``). Equations NSS serves that flowfreq cannot parse,
#:   and regions NSS would not resolve, are listed in ``unparsed``. Allowed only
#:   under ``nss_snapshot/``, and returned only on opt-in
#:   (``load_state(code, allow_nss_snapshot=True)``).
STATUSES = frozenset(
    {"pending", "partial", "verified", "provisional", "awaiting_nss", "nss_snapshot"}
)
#: Opt-in statuses, each with the one subdirectory allowed to hold it.
OPT_IN_SUBDIRS = {
    "provisional": PROVISIONAL_SUBDIR,
    "awaiting_nss": AWAITING_NSS_SUBDIR,
    "nss_snapshot": NSS_SNAPSHOT_SUBDIR,
}


class EquationsUnavailable(LookupError):
    """No equations are available for the requested state, region, or AEP."""


@dataclass
class RegionGate:
    """A rule that selects a region by a basin characteristic that is not a regression
    variable, such as Arizona's High Elevation region (mean basin elevation 7,500 ft or
    more, SIR 2014-5211).

    NSS applies these through "gating" parameters it does not publish thresholds for; the
    threshold comes from the report. ``minimum`` is inclusive and ``maximum`` exclusive,
    so complementary gates (``minimum=7500`` and ``maximum=7500``) partition the values.

    Attributes
    ----------
    variable : str
        StreamStats characteristic code (e.g. ``"ELEV"``).
    minimum, maximum : float, optional
        At least one is given.
    units : str
    source : str
        Where the threshold is printed.
    """

    variable: str
    minimum: Optional[float] = None
    maximum: Optional[float] = None
    units: str = ""
    source: str = ""

    def admits(self, value: float) -> bool:
        """Whether a basin with this characteristic value belongs to the region."""
        if self.minimum is not None and value < self.minimum:
            return False
        if self.maximum is not None and value >= self.maximum:
            return False
        return True

    def describe(self) -> str:
        """The rule in words, e.g. ``"ELEV >= 7500 ft"``."""
        parts = []
        if self.minimum is not None:
            parts.append(f"{self.variable} >= {self.minimum:g}")
        if self.maximum is not None:
            parts.append(f"{self.variable} < {self.maximum:g}")
        return " and ".join(parts) + (f" {self.units}" if self.units else "")

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RegionGate":
        gate = cls(
            variable=str(d["variable"]),
            minimum=None if d.get("minimum") is None else float(d["minimum"]),
            maximum=None if d.get("maximum") is None else float(d["maximum"]),
            units=str(d.get("units") or ""),
            source=str(d.get("source") or ""),
        )
        if gate.minimum is None and gate.maximum is None:
            raise ValueError(f"region gate on {gate.variable} has no minimum or maximum")
        return gate


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
    #: Published equations or regions that could not be stored (``nss_snapshot`` sets):
    #: dicts with ``region_code``, ``region_name``, ``statistic_code``, ``aep``,
    #: ``equation`` (empty for an unresolved region) and ``reason``.
    unparsed: List[Dict[str, Any]] = field(default_factory=list)
    #: When the NSS snapshot behind an ``nss_snapshot`` set was retrieved (ISO 8601 UTC).
    retrieved_utc: str = ""
    #: Region-selection rules on characteristics that are not regression variables
    #: (:class:`RegionGate`), by region code. ``assemble`` uses them to choose between
    #: overlapping regions.
    region_gates: Dict[str, RegionGate] = field(default_factory=dict)

    @property
    def provisional(self) -> bool:
        """True for a flowfreq-developed (non-USGS) equation set."""
        return self.status == "provisional"

    @property
    def awaiting_nss(self) -> bool:
        """True for a USGS-published set that NSS/StreamStats does not serve yet."""
        return self.status == "awaiting_nss"

    @property
    def nss_snapshot(self) -> bool:
        """True for a set harvested from NSS and not checked against the source report."""
        return self.status == "nss_snapshot"

    @property
    def regions(self) -> List[str]:
        """Region codes with at least one equation, sorted."""
        return sorted({e.region_code for e in self.equations})

    def aeps(self, region_code: str) -> List[float]:
        """AEPs published for a region, largest first.

        A printed constant (``form="constant"``, such as Nevada Region 6's 2-year
        ``Q=0``) is published, so its AEP is included.
        """
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
    gates = {
        str(code): RegionGate.from_dict(g) for code, g in (data.get("region_gates") or {}).items()
    }
    unknown = sorted(set(gates) - {e.region_code for e in eqs})
    if unknown:
        raise ValueError(f"{state}: region_gates for regions with no equations: {unknown}")
    return StateLibrary(
        state=state,
        status=status,
        equations=eqs,
        source_reports=list(data.get("source_reports", [])),
        issue=data.get("issue"),
        notes=str(data.get("notes") or ""),
        unparsed=[dict(u) for u in data.get("unparsed") or []],
        retrieved_utc=str(data.get("retrieved_utc") or ""),
        region_gates=gates,
    )


def _provisional_advisory(code: str) -> str:
    return (
        f"PROVISIONAL (not USGS): the {code} equations were developed by flowfreq following "
        "published USGS procedure and have not been reviewed or published by the USGS. The "
        f"published set is load_state({code!r}). A USGS report supersedes this one."
    )


def _awaiting_nss_advisory(code: str, lib: StateLibrary) -> str:
    reports = "; ".join(
        str(r.get("citation", "")) if isinstance(r, dict) else str(r) for r in lib.source_reports
    )
    return (
        f"NOT YET IN STREAMSTATS (USGS-published): the {code} equations from {reports} are "
        "published by the USGS but not yet served by NSS/StreamStats, which still uses the "
        f"set in load_state({code!r}). StreamStats estimates will not match these until NSS "
        "adopts them."
    )


def _nss_snapshot_advisory(code: str, lib: StateLibrary) -> str:
    date = lib.retrieved_utc[:10] or "an unrecorded date"
    return (
        f"NSS SNAPSHOT (not report-verified): the {code} equations were harvested from the "
        f"USGS NSS web service on {date} and parsed by flowfreq. They have not been checked "
        "against the source report and carry no covariance or model-error variance; "
        "intervals use NSS's average standard error of prediction where NSS reports one. "
        f"{len(lib.unparsed)} published equations or regions could not be parsed and are "
        "listed in StateLibrary.unparsed."
    )


def _load_opt_in(root: Path, code: str, status: str, warn: bool = True) -> Optional[StateLibrary]:
    """Load ``<subdir>/<CODE>.json`` for an opt-in status, or None if there is none."""
    path = root / OPT_IN_SUBDIRS[status] / f"{code}.json"
    if not path.exists():
        return None
    lib = parse_state(json.loads(path.read_text(encoding="utf-8")))
    if lib.status != status or lib.state != code:
        raise ValueError(
            f"{path}: a set under {OPT_IN_SUBDIRS[status]}/ must have status {status!r} and "
            f"state {code!r}, got {lib.status!r} / {lib.state!r}"
        )
    if status == "provisional":
        lib.advisory = _provisional_advisory(code)
    elif status == "awaiting_nss":
        lib.advisory = _awaiting_nss_advisory(code, lib)
    else:
        lib.advisory = _nss_snapshot_advisory(code, lib)
    log = logger.warning if warn else logger.debug
    log("%s", lib.advisory)
    return lib


def load_state(
    code: str,
    data_dir: Union[str, Path, None] = None,
    *,
    allow_provisional: bool = False,
    allow_awaiting_nss: bool = False,
    allow_nss_snapshot: bool = False,
    warn: bool = True,
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
    allow_awaiting_nss : bool
        Return the jurisdiction's USGS-published set that NSS/StreamStats does not
        serve yet, from ``awaiting_nss/<CODE>.json``, when there is one, with
        ``advisory`` set and a logged warning. Without one, the published set is
        returned as usual. If both flags are set and both files exist, this set
        wins, because a USGS report supersedes a flowfreq-developed one.
    allow_nss_snapshot : bool
        Return the jurisdiction's set harvested from the NSS web service, from
        ``nss_snapshot/<CODE>.json``, when there is one, with ``advisory`` set and
        a logged warning. Without one, the published set is returned as usual.
        It ranks last among the opt-in sets: ``awaiting_nss``, then
        ``provisional``, then ``nss_snapshot``.
    warn : bool
        Log the advisory of a returned opt-in set with ``logger.warning`` (the
        default). When False it is logged with ``logger.debug`` instead; ``advisory``
        is set either way.

    Returns
    -------
    StateLibrary

    Raises
    ------
    EquationsUnavailable
        If the jurisdiction has no file yet (its wave has not started).
    ValueError
        If a published file claims an opt-in status (``provisional``,
        ``awaiting_nss`` or ``nss_snapshot``), or a file in an opt-in folder claims
        any other status or another jurisdiction.
    """
    code = jurisdiction(code).code
    root = Path(data_dir or DATA_DIR)
    for status, wanted in (
        ("awaiting_nss", allow_awaiting_nss),
        ("provisional", allow_provisional),
        ("nss_snapshot", allow_nss_snapshot),
    ):
        if wanted:
            opt_in = _load_opt_in(root, code, status, warn)
            if opt_in is not None:
                return opt_in
    path = root / f"{code}.json"
    if not path.exists():
        raise EquationsUnavailable(f"{code}: no regression library file; its wave has not started")
    lib = parse_state(json.loads(path.read_text(encoding="utf-8")))
    if lib.status in OPT_IN_SUBDIRS:
        raise ValueError(
            f"{path}: status {lib.status!r} belongs under {OPT_IN_SUBDIRS[lib.status]}/, not in "
            "place of the published set"
        )
    return lib


def available_states(data_dir: Union[str, Path, None] = None) -> Dict[str, str]:
    """Map every jurisdiction with a published library file to its status.

    Opt-in sets are not listed; see :func:`provisional_states`,
    :func:`awaiting_nss_states` and :func:`nss_snapshot_states`.
    """
    out: Dict[str, str] = {}
    for p in sorted(Path(data_dir or DATA_DIR).glob("*.json")):
        out[p.stem] = parse_state(json.loads(p.read_text(encoding="utf-8"))).status
    return out


def provisional_states(data_dir: Union[str, Path, None] = None) -> List[str]:
    """Jurisdictions with a provisional (flowfreq-developed) equation set, sorted."""
    return sorted(p.stem for p in (Path(data_dir or DATA_DIR) / PROVISIONAL_SUBDIR).glob("*.json"))


def awaiting_nss_states(data_dir: Union[str, Path, None] = None) -> List[str]:
    """Jurisdictions with a USGS-published set NSS/StreamStats does not serve yet, sorted."""
    root = Path(data_dir or DATA_DIR) / AWAITING_NSS_SUBDIR
    return sorted(p.stem for p in root.glob("*.json"))


def nss_snapshot_states(data_dir: Union[str, Path, None] = None) -> List[str]:
    """Jurisdictions with a set harvested from NSS (not report-verified), sorted."""
    root = Path(data_dir or DATA_DIR) / NSS_SNAPSHOT_SUBDIR
    return sorted(p.stem for p in root.glob("*.json"))
