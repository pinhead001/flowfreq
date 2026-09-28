"""Legacy NWIS vs. Water Data OGC API: do both instantaneous-value backends agree?

Issue #29's gate for switching the default of
``USGSgage.download_instantaneous_flow`` / ``download_instantaneous_stage`` to
``"waterdata-ogc"``, modelled on ``tests/test_peak_backend_parity.py``. Both
backends read the same USGS database, so over the same window they must return
row for row the same frame:

- the same UTC index (same instants, none missing on either side) -- including
  across the autumn fall-back night, where legacy's local clock repeats an hour
  and the OGC backend has to derive that repeated hour itself;
- the same values, exactly (both are the database's published decimal);
- the same ``datetime_local`` and ``tz_cd``, which the OGC backend derives from
  the monitoring location's time zone while legacy reports them verbatim;
- the same ``qualification_code`` after one documented mapping, below.

**The one difference found (live, 2026-09-27): qualifier tokens legacy hides.**
On rows that *have* a value, the OGC API's ``qualifier`` list can carry tokens
the legacy RDB ``_cd`` column leaves out -- legacy only surfaces them as text in
place of a missing value (``Ice``, ``Eqp``). Seen: ``ICE`` on Methow at Pateros
(12449950) in January-February 2023 (legacy ``A:e``, OGC ``A:e:ICE``; legacy
``A``, OGC ``A:ICE``), and ``EQUIP`` on Big Sandy (03606500) on 2024-03-25
(legacy ``A:e``, OGC ``A:EQUIP:e``). The approval code and every token legacy
does report agree; values, instants and local times are identical. The library
keeps the extra tokens (``flowfreq.waterdata.map_qualification_code`` passes
unmapped tokens through verbatim -- more information, not different data), and
the comparison removes exactly the tokens in :data:`LEGACY_SUPPRESSED_TOKENS`
from the OGC side. Any *other* extra token is reported, so a new kind of
difference surfaces here rather than being absorbed.

Rows with no value (ice-affected with no estimate, equipment outage) are
dropped by both backends, and that is part of what is compared: 06214500's
provisional ice-null rows of January 2026 are absent from both.

**Multi-sensor sites** need an explicit ``ts_id`` on each backend, and the two
name a series differently -- legacy by DD number, OGC by 32-hex
``time_series_id``. :func:`map_dd_to_time_series_id` maps one to the other the
way a user would: legacy prints each DD number with a description
(``Gage height, feet, HEADWATER, [Headwater]``) whose sublocation token matches
the OGC series' ``sublocation_identifier``, among the series active in the
window. Olmsted (03612600) has two live 00065 series plus two retired
TAILWATER ones, so the window filter is load-bearing.

Needs *both* hosts. The legacy fetch skips -- rather than fails -- when NWIS is
unreachable (it is blocked from Claude Code *web* sessions; reachable from a CLI
session on a developer machine on 2026-09-27, where every case below agreed).
Run with::

    pytest tests/test_iv_backend_parity.py -m requires_network -v
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd
import pytest
import requests

from flowfreq.usgs import USGSgage
from flowfreq.waterdata import (
    TimeSeriesInfo,
    fetch_monitoring_location,
    list_instantaneous_series,
    local_day_bounds,
    resolve_local_zone,
)

#: OGC qualifier tokens seen on rows *with* a value that the legacy RDB ``_cd``
#: column omits. See the module docstring for where each was observed.
LEGACY_SUPPRESSED_TOKENS = frozenset({"ICE", "EQUIP"})

#: (site, method, start, end, why). Single-series sites: no ``ts_id`` needed.
SINGLE_SERIES_CASES: List[Tuple[str, str, str, str, str]] = [
    ("03606500", "flow", "2024-11-01", "2024-11-05", "Big Sandy, 2024-11-03 fall-back"),
    ("03606500", "stage", "2024-11-01", "2024-11-05", "Big Sandy stage, fall-back"),
    ("03606500", "flow", "2024-03-08", "2024-03-26", "spring-forward; EQUIP 2024-03-25"),
    ("12449950", "flow", "2024-05-20", "2024-05-24", "Methow at Pateros, snowmelt"),
    ("12449950", "stage", "2024-05-20", "2024-05-24", "Methow stage, snowmelt"),
    ("12449950", "flow", "2023-01-01", "2023-03-31", "Methow winter: ICE tokens"),
    ("06191500", "flow", "2024-01-10", "2024-01-20", "Yellowstone, MST no DST: A:e"),
    ("06214500", "flow", "2026-01-15", "2026-01-25", "Yellowstone at Billings: ice nulls"),
]

#: Olmsted, Ohio River: headwater and tailwater stage, both flagged primary.
MULTI_SENSOR_SITE = "03612600"
MULTI_SENSOR_WINDOW = ("2024-11-02", "2024-11-04")
#: What :func:`map_dd_to_time_series_id` must find (legacy header 2026-09-27).
EXPECTED_MAPPING = {
    "60629": "15beb94252164ce0b9a77a90edce7528",  # HEADWATER
    "323512": "7ca465077c26415c93fd0549c697a270",  # TAILWATER
}

_PARAM = {"flow": "00060", "stage": "00065"}


# ----------------------------------------------------------------------------
# Comparison helpers
# ----------------------------------------------------------------------------


def _split_code(code: object) -> Tuple[str, frozenset]:
    """``"A:e:ICE"`` -> (``"A"``, {``"e"``, ``"ICE"``}). The approval code leads."""
    parts = [p for p in str(code or "").split(":") if p]
    if not parts:
        return "", frozenset()
    return parts[0], frozenset(parts[1:])


def qualification_mismatch(legacy_code: object, ogc_code: object) -> Optional[str]:
    """Why two qualification codes disagree, or ``None`` if they agree.

    Approval codes must match. Qualifier tokens must match once the OGC side's
    :data:`LEGACY_SUPPRESSED_TOKENS` are set aside; token order is not
    meaningful.
    """
    la, lt = _split_code(legacy_code)
    oa, ot = _split_code(ogc_code)
    if la != oa:
        return f"approval {la!r} vs {oa!r}"
    missing = lt - ot
    extra = ot - lt - LEGACY_SUPPRESSED_TOKENS
    if missing or extra:
        return f"tokens only in legacy {sorted(missing)}, only in OGC {sorted(extra)}"
    return None


def iv_mismatches(legacy: pd.DataFrame, ogc: pd.DataFrame, value_col: str) -> List[str]:
    """Every disagreement between two instantaneous frames, as readable lines."""
    problems: List[str] = []
    for name, frame in (("legacy", legacy), ("OGC", ogc)):
        if list(frame.columns) != [value_col, "datetime_local", "tz_cd", "qualification_code"]:
            problems.append(f"{name} columns {list(frame.columns)}")
            return problems
        if str(frame.index.tz) != "UTC":
            problems.append(f"{name} index tz is {frame.index.tz}, not UTC")
        if not frame.index.is_unique or not frame.index.is_monotonic_increasing:
            problems.append(f"{name} index is not unique and increasing")

    only_legacy = legacy.index.difference(ogc.index)
    only_ogc = ogc.index.difference(legacy.index)
    if len(only_legacy):
        problems.append(
            f"{len(only_legacy)} instant(s) only in legacy, e.g. "
            f"{[t.isoformat() for t in only_legacy[:3]]}"
        )
    if len(only_ogc):
        problems.append(
            f"{len(only_ogc)} instant(s) only in OGC, e.g. "
            f"{[t.isoformat() for t in only_ogc[:3]]}"
        )

    common = legacy.index.intersection(ogc.index)
    a, b = legacy.loc[common], ogc.loc[common]
    for col in (value_col, "datetime_local", "tz_cd"):
        differ = a[col].to_numpy() != b[col].to_numpy()
        if differ.any():
            first = common[differ.argmax()]
            problems.append(
                f"{col}: {int(differ.sum())} row(s) differ, first at {first.isoformat()}: "
                f"legacy {a.at[first, col]!r} vs OGC {b.at[first, col]!r}"
            )

    reasons = [
        qualification_mismatch(x, y)
        for x, y in zip(a["qualification_code"], b["qualification_code"])
    ]
    bad = [(t, r) for t, r in zip(common, reasons) if r is not None]
    if bad:
        t, r = bad[0]
        problems.append(
            f"qualification_code: {len(bad)} row(s) differ, first at {t.isoformat()}: {r} "
            f"(legacy {a.at[t, 'qualification_code']!r}, OGC {b.at[t, 'qualification_code']!r})"
        )
    return problems


# ----------------------------------------------------------------------------
# Legacy DD number -> OGC time_series_id
# ----------------------------------------------------------------------------

_TS_HEADER = re.compile(r"^#\s+(\d+)\s+(\d{5})\s+(.+?)\s*$")


def legacy_series_descriptions(rdb_text: str, param_cd: str) -> Dict[str, str]:
    """``{DD number: description}`` from a legacy IV RDB header.

    The header lists every series in the response under ``TS_ID Parameter
    Description``, e.g. ``#    60629       00065     Gage height, feet,
    HEADWATER, [Headwater]``.
    """
    found: Dict[str, str] = {}
    for line in rdb_text.splitlines():
        if not line.startswith("#"):
            continue
        match = _TS_HEADER.match(line)
        if match and match.group(2) == param_cd:
            found[match.group(1)] = match.group(3)
    return found


def map_dd_to_time_series_id(
    description: str,
    series: Sequence[TimeSeriesInfo],
    window: Tuple[pd.Timestamp, pd.Timestamp],
) -> str:
    """The OGC ``time_series_id`` of the series a legacy description names.

    A candidate's ``sublocation_identifier`` must appear as one of the
    description's comma-separated fields (``HEADWATER``); a description with no
    such field matches only series without a sublocation. Series whose known
    ``begin``/``end`` lie wholly outside `window` are not candidates, since a
    site can carry retired series for the same sublocation.

    Raises
    ------
    LookupError
        Zero or several candidates -- the mapping is never guessed.
    """
    fields = {f.strip().upper() for f in description.split(",")}
    lo, hi = window

    def active(s: TimeSeriesInfo) -> bool:
        return not ((s.end is not None and s.end < lo) or (s.begin is not None and s.begin > hi))

    def subloc(s: TimeSeriesInfo) -> str:
        return (s.sublocation_identifier or "").strip().upper()

    named = [s for s in series if subloc(s) and subloc(s) in fields]
    pool = named if named else [s for s in series if not subloc(s)]
    candidates = [s for s in pool if active(s)]
    if len(candidates) != 1:
        raise LookupError(
            f"{len(candidates)} OGC series match legacy description {description!r} in "
            f"the window: {[s.describe() for s in candidates]}"
        )
    return candidates[0].time_series_id


# ----------------------------------------------------------------------------
# Live
# ----------------------------------------------------------------------------


def _skip_if_legacy_unreachable(exc: requests.RequestException) -> None:
    """Skip on a transport failure, a proxy block or a legacy 5xx; else fail.

    The legacy service is being retired and answered a transient 503 during
    this work, which says nothing about parity. The legacy IV path re-raises a
    chunk failure as a plain ``RequestException`` naming the window, so the
    original is its ``__cause__``.
    """
    for err in (exc, exc.__cause__):
        if isinstance(err, (requests.ConnectionError, requests.Timeout)):
            pytest.skip(f"legacy NWIS unreachable from here: {exc}")
        if isinstance(err, requests.HTTPError) and err.response is not None:
            status = err.response.status_code
            if status in (403, 407) or status >= 500:
                pytest.skip(f"legacy NWIS unavailable from here (HTTP {status}): {exc}")


def _legacy(site_no: str, method: str, start: str, end: str, **kw: object) -> pd.DataFrame:
    fetch = getattr(USGSgage(site_no), f"download_instantaneous_{method}")
    try:
        return fetch(start, end, backend="nwis-legacy", **kw)
    except requests.RequestException as exc:
        _skip_if_legacy_unreachable(exc)
        raise


def _legacy_rdb(site_no: str, param_cd: str, day: str) -> str:
    try:
        response = requests.get(
            USGSgage.BASE_URL_IV,
            params={
                "format": "rdb",
                "sites": site_no,
                "parameterCd": param_cd,
                "startDT": day,
                "endDT": day,
            },
            timeout=60,
        )
        response.raise_for_status()
    except requests.RequestException as exc:
        _skip_if_legacy_unreachable(exc)
        raise
    return response.text


@pytest.mark.requires_network
@pytest.mark.parametrize(
    ("site_no", "method", "start", "end", "why"),
    SINGLE_SERIES_CASES,
    ids=[f"{c[0]}-{c[1]}-{c[2]}" for c in SINGLE_SERIES_CASES],
)
def test_backends_agree(site_no: str, method: str, start: str, end: str, why: str) -> None:
    legacy = _legacy(site_no, method, start, end)
    ogc = getattr(USGSgage(site_no), f"download_instantaneous_{method}")(
        start, end, backend="waterdata-ogc"
    )
    problems = iv_mismatches(legacy, ogc, legacy.columns[0])
    assert not problems, f"{site_no} {method} {start}..{end} ({why}):\n" + "\n".join(problems)


@pytest.mark.requires_network
def test_fall_back_night_is_on_both_backends() -> None:
    """The case the UTC axis exists for: 01:00-01:59 local happens twice."""
    legacy = _legacy("03606500", "flow", "2024-11-03", "2024-11-03")
    ogc = USGSgage("03606500").download_instantaneous_flow("2024-11-03", "2024-11-03")
    for frame in (legacy, ogc):
        repeated = frame[frame["datetime_local"].dt.hour == 1]
        assert set(repeated["tz_cd"]) == {"CDT", "CST"}
        assert repeated["datetime_local"].duplicated().any()
    assert not iv_mismatches(legacy, ogc, "flow_cfs")


@pytest.mark.requires_network
@pytest.mark.parametrize("method", ["flow", "stage"])
def test_default_backend_agrees_with_legacy(method: str) -> None:
    """What every existing caller now gets must match what it got."""
    start, end = "2024-11-01", "2024-11-05"
    legacy = _legacy("03606500", method, start, end)
    new = getattr(USGSgage("03606500"), f"download_instantaneous_{method}")(start, end)
    problems = iv_mismatches(legacy, new, legacy.columns[0])
    assert not problems, "\n".join(problems)


@pytest.mark.requires_network
def test_multi_sensor_site_agrees_series_by_series() -> None:
    """Olmsted stage: map each legacy DD number to its OGC series, then compare."""
    start, end = MULTI_SENSOR_WINDOW
    descriptions = legacy_series_descriptions(
        _legacy_rdb(MULTI_SENSOR_SITE, "00065", start), "00065"
    )
    assert len(descriptions) == 2, descriptions

    gage = USGSgage(MULTI_SENSOR_SITE)
    series = list_instantaneous_series(MULTI_SENSOR_SITE, "00065")
    loc = fetch_monitoring_location(MULTI_SENSOR_SITE)
    zone = resolve_local_zone(loc.get("time_zone_abbreviation"), loc.get("uses_daylight_savings"))
    window = local_day_bounds(start, end, zone)

    mapping = {
        dd: map_dd_to_time_series_id(desc, series, window) for dd, desc in descriptions.items()
    }
    assert mapping == EXPECTED_MAPPING

    for dd, uuid in mapping.items():
        legacy = _legacy(MULTI_SENSOR_SITE, "stage", start, end, ts_id=dd)
        ogc = gage.download_instantaneous_stage(start, end, ts_id=uuid)
        problems = iv_mismatches(legacy, ogc, "gage_height_ft")
        assert not problems, f"DD {dd} / {uuid} ({descriptions[dd]}):\n" + "\n".join(problems)

    # Neither backend picks a sensor for the caller.
    with pytest.raises(ValueError, match="ts_id"):
        gage.download_instantaneous_stage(start, end)


# ----------------------------------------------------------------------------
# Offline: the comparison itself, so a live pass means what it says
# ----------------------------------------------------------------------------


def _frame(rows: list, value_col: str = "flow_cfs") -> pd.DataFrame:
    """Rows of (utc, value, local, tz_cd, code)."""
    index = pd.DatetimeIndex([pd.Timestamp(r[0], tz="UTC") for r in rows], name="datetime")
    return pd.DataFrame(
        {
            value_col: [float(r[1]) for r in rows],
            "datetime_local": pd.to_datetime([r[2] for r in rows]),
            "tz_cd": [r[3] for r in rows],
            "qualification_code": [r[4] for r in rows],
        },
        index=index,
    )


_FALL_BACK = [
    ("2024-11-03 05:30", 100, "2024-11-03 00:30", "CDT", "A"),
    ("2024-11-03 06:30", 101, "2024-11-03 01:30", "CDT", "A"),
    ("2024-11-03 07:30", 102, "2024-11-03 01:30", "CST", "A:e"),
]


class TestIvMismatches:
    def test_identical_frames_agree(self) -> None:
        f = _frame(_FALL_BACK)
        assert iv_mismatches(f, f.copy(), "flow_cfs") == []

    def test_legacy_suppressed_tokens_are_tolerated(self) -> None:
        a = _frame(_FALL_BACK)
        b = a.copy()
        b["qualification_code"] = ["A:ICE", "A", "A:EQUIP:e"]
        assert iv_mismatches(a, b, "flow_cfs") == []

    def test_token_order_does_not_matter(self) -> None:
        assert qualification_mismatch("A:e:ICE", "A:ICE:e") is None

    def test_an_unknown_extra_token_is_reported(self) -> None:
        a = _frame(_FALL_BACK)
        b = a.copy()
        b["qualification_code"] = ["A:Bkw", "A", "A:e"]
        report = "\n".join(iv_mismatches(a, b, "flow_cfs"))
        assert "only in OGC ['Bkw']" in report

    def test_a_token_missing_from_ogc_is_reported(self) -> None:
        assert "only in legacy ['e']" in str(qualification_mismatch("A:e", "A:ICE"))

    def test_approval_difference_is_reported_even_with_suppressed_tokens(self) -> None:
        assert "approval 'A' vs 'P'" in str(qualification_mismatch("A", "P:ICE"))

    def test_every_kind_of_disagreement_is_reported(self) -> None:
        a = _frame(_FALL_BACK)
        b = _frame(
            [
                ("2024-11-03 05:30", 100.5, "2024-11-03 00:30", "CDT", "A"),
                ("2024-11-03 06:30", 101, "2024-11-03 00:30", "CST", "A"),
                ("2024-11-03 08:30", 102, "2024-11-03 02:30", "CST", "A"),
            ]
        )
        report = "\n".join(iv_mismatches(a, b, "flow_cfs"))
        assert "1 instant(s) only in legacy" in report
        assert "1 instant(s) only in OGC" in report
        assert "flow_cfs: 1 row(s) differ" in report
        assert "datetime_local: 1 row(s) differ" in report
        assert "tz_cd: 1 row(s) differ" in report

    def test_a_non_utc_or_disordered_index_is_reported(self) -> None:
        a = _frame(_FALL_BACK)
        b = a.copy()
        b.index = b.index.tz_convert("America/Chicago")
        assert any("not UTC" in p for p in iv_mismatches(a, b, "flow_cfs"))
        c = a.iloc[::-1]
        assert any("not unique and increasing" in p for p in iv_mismatches(a, c, "flow_cfs"))

    def test_wrong_columns_are_reported(self) -> None:
        a = _frame(_FALL_BACK)
        b = a.rename(columns={"flow_cfs": "gage_height_ft"})
        assert "OGC columns" in iv_mismatches(a, b, "flow_cfs")[0]


class TestLegacySkip:
    """Which legacy failures skip the live gate, and which fail it."""

    @staticmethod
    def _wrapped(status: int) -> requests.RequestException:
        response = requests.Response()
        response.status_code = status
        wrapped = requests.RequestException("Instantaneous-value request failed ...")
        wrapped.__cause__ = requests.HTTPError(f"{status}", response=response)
        return wrapped

    @pytest.mark.parametrize("status", [403, 407, 502, 503])
    def test_blocked_or_down_skips(self, status: int) -> None:
        with pytest.raises(pytest.skip.Exception):
            _skip_if_legacy_unreachable(self._wrapped(status))

    def test_connection_error_skips(self) -> None:
        with pytest.raises(pytest.skip.Exception):
            _skip_if_legacy_unreachable(requests.ConnectionError("refused"))

    @pytest.mark.parametrize("status", [400, 404])
    def test_client_errors_do_not_skip(self, status: int) -> None:
        _skip_if_legacy_unreachable(self._wrapped(status))  # returns; caller re-raises


_OLMSTED_HEADER = """\
# Data provided for site 03612600
#    TS_ID       Parameter Description
#    60629       00065     Gage height, feet, HEADWATER, [Headwater]
#    323512      00065     Gage height, feet, TAILWATER
#    151767      00060     Discharge, cubic feet per second
#
agency_cd\tsite_no\tdatetime
"""


def _series(
    ts: str, subloc: Optional[str], begin: Optional[str], end: Optional[str]
) -> TimeSeriesInfo:
    return TimeSeriesInfo(
        time_series_id=ts,
        parameter_code="00065",
        sublocation_identifier=subloc,
        web_description=None,
        primary="Primary",
        begin=pd.Timestamp(begin, tz="UTC") if begin else None,
        end=pd.Timestamp(end, tz="UTC") if end else None,
    )


#: Shaped on Olmsted's live 00065 instantaneous metadata (2026-09-27): two
#: current series and a TAILWATER one retired in 2023.
_OLMSTED_SERIES = [
    _series("15beb94252164ce0b9a77a90edce7528", "HEADWATER", "2013-03-02", "2026-09-27"),
    _series("5407ef71596646da90836fefa909757a", "TAILWATER", "2019-02-24", "2023-04-01"),
    _series("7ca465077c26415c93fd0549c697a270", "TAILWATER", "2023-04-13", "2026-09-27"),
]
_WINDOW_2024 = (pd.Timestamp("2024-11-02 05:00", tz="UTC"), pd.Timestamp("2024-11-05", tz="UTC"))


class TestDdMapping:
    def test_header_is_parsed_per_parameter(self) -> None:
        assert legacy_series_descriptions(_OLMSTED_HEADER, "00065") == {
            "60629": "Gage height, feet, HEADWATER, [Headwater]",
            "323512": "Gage height, feet, TAILWATER",
        }
        assert list(legacy_series_descriptions(_OLMSTED_HEADER, "00060")) == ["151767"]

    def test_sublocation_and_window_pick_one_series(self) -> None:
        descriptions = legacy_series_descriptions(_OLMSTED_HEADER, "00065")
        mapping = {
            dd: map_dd_to_time_series_id(d, _OLMSTED_SERIES, _WINDOW_2024)
            for dd, d in descriptions.items()
        }
        assert mapping == EXPECTED_MAPPING

    def test_retired_series_is_chosen_for_an_old_window(self) -> None:
        window = (pd.Timestamp("2020-06-01", tz="UTC"), pd.Timestamp("2020-06-02", tz="UTC"))
        uuid = map_dd_to_time_series_id("Gage height, feet, TAILWATER", _OLMSTED_SERIES, window)
        assert uuid == "5407ef71596646da90836fefa909757a"

    def test_ambiguity_raises_rather_than_guessing(self) -> None:
        both = [
            _series("a" * 32, "TAILWATER", None, None),
            _series("b" * 32, "TAILWATER", None, None),
        ]
        with pytest.raises(LookupError, match="2 OGC series"):
            map_dd_to_time_series_id("Gage height, feet, TAILWATER", both, _WINDOW_2024)

    def test_plain_description_matches_the_series_without_sublocation(self) -> None:
        series = [_series("c" * 32, None, None, None), *_OLMSTED_SERIES]
        assert map_dd_to_time_series_id("Gage height, feet", series, _WINDOW_2024) == "c" * 32

    def test_no_match_raises(self) -> None:
        with pytest.raises(LookupError, match="0 OGC series"):
            map_dd_to_time_series_id("Gage height, feet, BACKUP", _OLMSTED_SERIES, _WINDOW_2024)
