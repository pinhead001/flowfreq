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
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Union

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
class RegionLimit:
    """Where an equation set can apply, for telling it apart from an overlapping set.

    A set excluded by a limit is no alternative for this basin, so it drops out of the overlap
    without the caller choosing. It stays reachable with ``region=``. Unlike a
    :class:`RegionGate`, which assigns a basin among sets that all have gates, a limit needs no
    gate on the other sets, and a basin it does not exclude remains the caller's to decide:
    an urban set's limit sends a rural basin to the rural set, and leaves an urban one alone.

    Attributes
    ----------
    variable : str
        StreamStats characteristic code. A value that is unknown excludes nothing.
    minimum, maximum : float, optional
        Both inclusive, as an equation's range is.
    below : str, optional
        Another characteristic that ``variable`` must be less than, as Oklahoma's regulated
        sets need ``DAUNREG`` under ``CONTDA``.
    units : str
    source : str
        Where the limit is printed, or that it is a flowfreq method choice.
    """

    variable: str
    minimum: Optional[float] = None
    maximum: Optional[float] = None
    below: Optional[str] = None
    units: str = ""
    source: str = ""

    def excludes(self, characteristics: Mapping[str, float]) -> bool:
        """Whether this basin is known to be outside the limit."""
        value = characteristics.get(self.variable)
        if value is None or not math.isfinite(value):
            return False
        if self.minimum is not None and value < self.minimum:
            return True
        if self.maximum is not None and value > self.maximum:
            return True
        if self.below is not None:
            other = characteristics.get(self.below)
            if other is not None and math.isfinite(other) and value >= other:
                return True
        return False

    def describe(self) -> str:
        """The rule in words, e.g. ``"DRNAREA <= 1 mi2"``."""
        parts = []
        if self.minimum is not None:
            parts.append(f"{self.variable} >= {self.minimum:g}")
        if self.maximum is not None:
            parts.append(f"{self.variable} <= {self.maximum:g}")
        if self.below is not None:
            parts.append(f"{self.variable} < {self.below}")
        return " and ".join(parts) + (f" {self.units}" if self.units else "")

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RegionLimit":
        limit = cls(
            variable=str(d["variable"]),
            minimum=None if d.get("minimum") is None else float(d["minimum"]),
            maximum=None if d.get("maximum") is None else float(d["maximum"]),
            below=None if d.get("below") is None else str(d["below"]),
            units=str(d.get("units") or ""),
            source=str(d.get("source") or ""),
        )
        if limit.minimum is None and limit.maximum is None and limit.below is None:
            raise ValueError(f"region limit on {limit.variable} has no bound")
        return limit


@dataclass(frozen=True)
class KnownMiss:
    """A gage where the regression is known to miss the record by far, with no fault found.

    Estimates near it carry a caution. It changes no value: the stored equations are what the
    report prints. Each entry records one investigated case in ``docs/VERIFICATION_PLAN.md``.

    Attributes
    ----------
    site_no, name : str
    lat, lon : float
        The gage, decimal degrees.
    radius_km : float
        How far from the gage an estimate is cautioned.
    regression_q1pct, at_site_q1pct : float
        Flood flows at the 1 percent AEP, ft3/s. The at-site value is Bulletin 17C with
        station skew.
    note : str
        What else is known, in a sentence.
    checked : str
        The date of the check.
    """

    site_no: str
    name: str
    lat: float
    lon: float
    radius_km: float
    regression_q1pct: float
    at_site_q1pct: float
    note: str = ""
    checked: str = ""

    def distance_km(self, lat: float, lon: float) -> float:
        """Great-circle distance from the gage to a point."""
        p1, p2 = math.radians(self.lat), math.radians(lat)
        dp, dl = p2 - p1, math.radians(lon - self.lon)
        a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
        return 2 * 6371.0 * math.asin(math.sqrt(a))

    def caution(self) -> str:
        """The caution in words."""
        ratio = self.regression_q1pct / self.at_site_q1pct
        size = (
            f"about {ratio:.1f} times the record's"
            if ratio >= 1
            else f"about 1/{1 / ratio:.1f} of the record's"
        )
        return (
            f"Known miss nearby: at USGS {self.site_no} {self.name}, the regression's 1 percent flow "
            f"of {self.regression_q1pct:,.0f} ft3/s is {size} {self.at_site_q1pct:,.0f}. "
            + (f"{self.note} " if self.note else "")
            + "Compare with gage records before relying on this estimate."
        )

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "KnownMiss":
        return cls(
            site_no=str(d["site_no"]),
            name=str(d["name"]),
            lat=float(d["lat"]),
            lon=float(d["lon"]),
            radius_km=float(d["radius_km"]),
            regression_q1pct=float(d["regression_q1pct"]),
            at_site_q1pct=float(d["at_site_q1pct"]),
            note=str(d.get("note") or ""),
            checked=str(d.get("checked") or ""),
        )


@dataclass(frozen=True)
class RegionSelector:
    """A basin characteristic that names the region at the point of interest.

    Indiana's report selects the equation by where the point lies, not by area weighting, and
    StreamStats returns that region as ``HIGHREG``. Where the value is known and listed, the region
    is used alone for the whole basin, whatever NSS's ``percent_weight`` says. ``region=`` overrides.

    Attributes
    ----------
    variable : str
        StreamStats characteristic code.
    values : dict of str to str
        The characteristic's value, as a whole number in text, to the library ``region_code``.
    source : str
        Where the rule is printed.
    """

    variable: str
    values: Dict[str, str]
    source: str = ""

    def region(self, characteristics: Mapping[str, float]) -> Optional[str]:
        """The library region this basin's value names, or None if unknown or unlisted."""
        value = characteristics.get(self.variable)
        if value is None or not math.isfinite(value) or value != int(value):
            return None
        return self.values.get(str(int(value)))

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RegionSelector":
        return cls(
            variable=str(d["variable"]),
            values={str(k): str(v) for k, v in d["values"].items()},
            source=str(d.get("source") or ""),
        )


@dataclass
class RegionPreference:
    """An order of preference among overlapping equation sets.

    A method choice, not a published rule, which is why it is not a :class:`RegionGate`.
    Iowa publishes a full set (SIR 2013-5086) and a drainage-area-only set (SIR 2015-5055)
    for each region, and nothing printed says which to use where both apply.

    Attributes
    ----------
    tiers : list of list of str
        Region codes, most preferred tier first. ``assemble`` uses the first tier whose
        regions, among those the basin has the characteristics for, cover the basin.
    basis : str
        Why this order, recorded with the estimate.
    """

    tiers: List[List[str]]
    basis: str = ""

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "RegionPreference":
        return cls(tiers=[[str(c) for c in t] for t in d["tiers"]], basis=str(d.get("basis") or ""))


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
    #: The order to try overlapping sets in when no published rule decides it.
    region_preference: Optional[RegionPreference] = None
    #: Where a set can apply (:class:`RegionLimit`), by region code. ``assemble`` drops a set
    #: a limit excludes from an overlap, before any other rule.
    region_limits: Dict[str, List[RegionLimit]] = field(default_factory=dict)
    #: The characteristic that names the region at the point, where a report selects that way.
    region_selector: Optional[RegionSelector] = None
    #: Investigated gages where the regression misses the record by far, with no fault found.
    known_misses: List[KnownMiss] = field(default_factory=list)
    #: A caution added to every estimate that uses the region, keyed by library region code.
    region_cautions: Dict[str, str] = field(default_factory=dict)
    #: The region used when NSS locates none of the library's regions and the caller names none,
    #: with why: ``{"region": code, "basis": text}``.
    fallback_region: Optional[Dict[str, str]] = None
    #: How a basin spanning regions is averaged, where the report says: "linear" (in flow, as
    #: South Dakota's procedure 2 does) or "log". None leaves it to the caller's default.
    weighting_space: Optional[str] = None
    #: A region's companion set, used beside it wherever NSS locates the region but not the
    #: companion: New Hampshire's and New York's drainage-area-only sets, which NSS never
    #: locates. Keyed by library region code.
    companion_regions: Dict[str, str] = field(default_factory=dict)

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
    pref = data.get("region_preference")
    preference = RegionPreference.from_dict(pref) if pref else None
    if preference is not None:
        unknown = sorted({c for t in preference.tiers for c in t} - {e.region_code for e in eqs})
        if unknown:
            raise ValueError(f"{state}: region_preference for regions with no equations: {unknown}")
    limits = {
        str(code): [RegionLimit.from_dict(x) for x in entries]
        for code, entries in (data.get("region_limits") or {}).items()
    }
    unknown = sorted(set(limits) - {e.region_code for e in eqs})
    if unknown:
        raise ValueError(f"{state}: region_limits for regions with no equations: {unknown}")
    regions_known = {e.region_code for e in eqs}
    cautions = {str(k): str(v) for k, v in (data.get("region_cautions") or {}).items()}
    unknown = sorted(set(cautions) - regions_known)
    if unknown:
        raise ValueError(f"{state}: region_cautions for regions with no equations: {unknown}")
    fallback = data.get("fallback_region")
    if fallback is not None:
        fallback = {"region": str(fallback["region"]), "basis": str(fallback.get("basis") or "")}
        if fallback["region"] not in regions_known:
            raise ValueError(f"{state}: fallback_region {fallback['region']!r} has no equations")
    weighting_space = data.get("weighting_space")
    if weighting_space not in (None, "log", "linear"):
        raise ValueError(
            f"{state}: weighting_space must be 'log' or 'linear', got {weighting_space!r}"
        )
    companions = {str(k): str(v) for k, v in (data.get("companion_regions") or {}).items()}
    unknown = sorted((set(companions) | set(companions.values())) - regions_known)
    if unknown:
        raise ValueError(f"{state}: companion_regions for regions with no equations: {unknown}")
    selector_data = data.get("region_selector")
    selector = RegionSelector.from_dict(selector_data) if selector_data else None
    if selector is not None:
        unknown = sorted(set(selector.values.values()) - {e.region_code for e in eqs})
        if unknown:
            raise ValueError(f"{state}: region_selector for regions with no equations: {unknown}")
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
        region_preference=preference,
        region_limits=limits,
        region_selector=selector,
        known_misses=[KnownMiss.from_dict(x) for x in data.get("known_misses") or []],
        region_cautions=cautions,
        fallback_region=fallback,
        weighting_space=weighting_space,
        companion_regions=companions,
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
        ``nss_snapshot/<CODE>.json``, when there is one and the published set is not
        ``verified``, with ``advisory`` set and a logged warning. Otherwise the published
        set is returned as usual: the flag permits a snapshot where nothing better exists,
        and a state that gains a verified set keeps its snapshot file as the evidence that
        the two agree. It ranks last among the opt-in sets: ``awaiting_nss``, then
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
    ):
        if wanted:
            opt_in = _load_opt_in(root, code, status, warn)
            if opt_in is not None:
                return opt_in
    path = root / f"{code}.json"
    lib = parse_state(json.loads(path.read_text(encoding="utf-8"))) if path.exists() else None
    if lib is not None and lib.status in OPT_IN_SUBDIRS:
        raise ValueError(
            f"{path}: status {lib.status!r} belongs under {OPT_IN_SUBDIRS[lib.status]}/, not in "
            "place of the published set"
        )
    if allow_nss_snapshot and (lib is None or lib.status != "verified"):
        opt_in = _load_opt_in(root, code, "nss_snapshot", warn)
        if opt_in is not None:
            return opt_in
    if lib is None:
        raise EquationsUnavailable(f"{code}: no regression library file; its wave has not started")
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
