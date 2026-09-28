"""Legacy NWIS vs. Water Data OGC API: do both daily-value backends agree?

Issue #29's gate for switching ``USGSgage.download_daily_flow`` to
``backend="waterdata-ogc"``, modelled on ``tests/test_peak_backend_parity.py``.
For the same site and window both must return:

- the same days (index ``date``, the gage's local calendar day);
- the same daily mean on each, exactly (both publish the same decimal);
- the same qualification, comparing legacy's RDB ``_cd`` code with
  ``flowfreq.waterdata.download_daily``'s ``qualification_code``. The approval
  code and every legacy token must match; OGC may add the tokens legacy's
  daily RDB never shows (:data:`LEGACY_SUPPRESSED_TOKENS`: ``ICE`` on an
  ice-affected estimate, ``EQUIP`` on an equipment-affected one -- legacy
  writes both as plain ``A:e``).

**Ice.** An ice-affected day is either an estimate or nothing, and the two
backends agree on which (live, 2026-09-27):

- *estimate* -- numeric on both; OGC ``["ESTIMATED", "ICE"]`` (``A:e:ICE``),
  legacy ``A:e``. Methow at Twisp (12449500), 2022-12-20..24.
- *no value* -- legacy writes the text ``Ice`` in the value column, which its
  parser coerces to NaN and drops; the API returns ``value: null`` with
  ``["ICE"]``, which ``parse_daily_features`` drops. 06214500 has 26 such
  provisional days in WY2026 (e.g. 2026-01-18..25), absent from both frames.

A day that legacy writes as ``Ice`` but the API publishes as a number would be
the one *deliberate* difference: the API's value is USGS's estimate, and legacy
simply never exposed it. :func:`daily_mismatches` classifies such days
separately (``ice_estimates_only_in_ogc``) rather than as a failure. None was
found on any record compared.

Also compared, over whole records through ``download_daily_flow``'s defaults
(``DEFAULT_START_DATE`` to today, UTC): Big Sandy 03606500 (1929-) and Methow
at Twisp 12449500 (1919-) -- every day identical.

Needs *both* hosts. The legacy fetch skips -- rather than fails -- when NWIS is
unreachable, blocked or answering 5xx (it returned a transient 503 during this
work). Run with::

    pytest tests/test_daily_backend_parity.py -m requires_network -v
"""

from __future__ import annotations

from dataclasses import dataclass, field
from io import StringIO
from typing import List, Optional, Tuple
from unittest.mock import Mock, patch

import pandas as pd
import pytest
import requests

from flowfreq.usgs import USGSgage
from flowfreq.waterdata import download_daily, parse_daily_features
from tests.fixtures.waterdata_ogc_responses import load, load_rdb

#: OGC qualifier tokens the legacy daily RDB's ``_cd`` column never shows.
LEGACY_SUPPRESSED_TOKENS = frozenset({"ICE", "EQUIP"})

#: (site, start, end, why) compared live, one water year each.
WATER_YEARS: List[Tuple[str, str, str, str]] = [
    ("03606500", "2023-10-01", "2024-09-30", "Big Sandy WY2024: A:e and EQUIP days"),
    ("12449500", "2022-10-01", "2023-09-30", "Methow at Twisp WY2023: ice estimates"),
    ("06214500", "2025-10-01", "2026-09-26", "Yellowstone at Billings WY2026: ice nulls"),
]

#: Whole records through download_daily_flow's own default range.
FULL_RECORD_SITES = ["03606500", "12449500"]

#: (fixture stem, legacy RDB stem) pairs captured 2026-09-27.
CAPTURED_PAIRS = [
    ("daily_03606500_00060_equip", "dv_03606500_2024-03-20_2024-03-31"),
    ("daily_12449500_00060_ice_estimated", "dv_12449500_2022-12-15_2022-12-28"),
    ("daily_06214500_00060_ice_null", "dv_06214500_2026-01-10_2026-01-25"),
]


# ----------------------------------------------------------------------------
# Comparison helpers
# ----------------------------------------------------------------------------


@dataclass
class LegacyDaily:
    """A legacy daily-values RDB, parsed two ways."""

    frame: pd.DataFrame  # exactly what download_daily_flow(backend="nwis-legacy") returns
    codes: pd.Series  # the ``_cd`` column, by date, for every row
    text_values: pd.Series  # non-numeric value text (``Ice``, ``Eqp``), by date


def parse_legacy_daily(rdb_text: str, site_no: str = "00000000") -> LegacyDaily:
    """Run the real legacy parser on `rdb_text`, and read its codes alongside."""
    response = Mock(text=rdb_text, status_code=200)
    response.raise_for_status.return_value = None
    with patch("flowfreq.usgs.requests.get", return_value=response):
        frame = USGSgage(site_no).download_daily_flow(
            "1850-01-01", "2100-01-01", backend="nwis-legacy"
        )

    lines = [l for l in rdb_text.splitlines() if l.strip() and not l.startswith("#")]
    raw = pd.read_csv(StringIO("\n".join(lines)), sep="\t", skiprows=[1], dtype=str)
    value_col = next(c for c in raw.columns if "00060" in c and not c.endswith("_cd"))
    index = pd.DatetimeIndex(pd.to_datetime(raw["datetime"]), name="date")
    values = pd.Series(raw[value_col].to_numpy(), index=index)
    codes = pd.Series(raw[f"{value_col}_cd"].fillna("").to_numpy(), index=index)
    numeric = pd.to_numeric(values, errors="coerce")
    return LegacyDaily(frame, codes, values[numeric.isna()])


def _split_code(code: object) -> Tuple[str, frozenset]:
    parts = [p for p in str(code or "").split(":") if p]
    return (parts[0], frozenset(parts[1:])) if parts else ("", frozenset())


def qualification_mismatch(legacy_code: object, ogc_code: object) -> Optional[str]:
    """Why two codes disagree (``None`` if they agree); see the module docstring."""
    la, lt = _split_code(legacy_code)
    oa, ot = _split_code(ogc_code)
    if la != oa:
        return f"approval {la!r} vs {oa!r}"
    missing, extra = lt - ot, ot - lt - LEGACY_SUPPRESSED_TOKENS
    if missing or extra:
        return f"tokens only in legacy {sorted(missing)}, only in OGC {sorted(extra)}"
    return None


@dataclass
class DailyComparison:
    problems: List[str] = field(default_factory=list)
    #: Days legacy wrote as ``Ice`` text but the API publishes a number for --
    #: the deliberate, documented difference. Not a problem.
    ice_estimates_only_in_ogc: List[pd.Timestamp] = field(default_factory=list)


def daily_mismatches(legacy: LegacyDaily, ogc: pd.DataFrame) -> DailyComparison:
    """Every disagreement between a legacy daily record and an OGC one.

    `ogc` is :func:`flowfreq.waterdata.download_daily`'s frame
    (``flow_cfs`` and ``qualification_code``).
    """
    result = DailyComparison()
    lg = legacy.frame
    for name, frame in (("legacy", lg), ("OGC", ogc)):
        if frame.index.name != "date" or frame.index.tz is not None:
            result.problems.append(f"{name} index is {frame.index.name!r}/{frame.index.tz}")
        if not frame.index.is_unique or not frame.index.is_monotonic_increasing:
            result.problems.append(f"{name} index is not unique and increasing")

    only_legacy = lg.index.difference(ogc.index)
    only_ogc = ogc.index.difference(lg.index)
    ice_text = {d for d, v in legacy.text_values.items() if str(v).strip().lower() == "ice"}
    result.ice_estimates_only_in_ogc = [d for d in only_ogc if d in ice_text]
    unexplained = [d for d in only_ogc if d not in ice_text]
    if len(only_legacy):
        result.problems.append(
            f"{len(only_legacy)} day(s) only in legacy, e.g. {[str(d.date()) for d in only_legacy[:3]]}"
        )
    if unexplained:
        result.problems.append(
            f"{len(unexplained)} day(s) only in OGC, e.g. {[str(d.date()) for d in unexplained[:3]]}"
        )

    common = lg.index.intersection(ogc.index)
    a = lg.loc[common, "flow_cfs"].astype(float).to_numpy()
    b = ogc.loc[common, "flow_cfs"].astype(float).to_numpy()
    differ = a != b
    if differ.any():
        first = common[differ.argmax()]
        result.problems.append(
            f"flow_cfs: {int(differ.sum())} day(s) differ, first {first.date()}: legacy "
            f"{lg.at[first, 'flow_cfs']} vs OGC {ogc.at[first, 'flow_cfs']}"
        )

    bad = [
        (d, r)
        for d in common
        if (r := qualification_mismatch(legacy.codes.get(d), ogc.at[d, "qualification_code"]))
    ]
    if bad:
        d, r = bad[0]
        result.problems.append(
            f"qualification: {len(bad)} day(s) differ, first {d.date()}: {r} (legacy "
            f"{legacy.codes.get(d)!r}, OGC {ogc.at[d, 'qualification_code']!r})"
        )
    return result


# ----------------------------------------------------------------------------
# Live
# ----------------------------------------------------------------------------


def _legacy_rdb(site_no: str, start: str, end: str) -> str:
    """The legacy daily RDB, fetched as ``download_daily_flow`` fetches it."""
    try:
        response = requests.get(
            USGSgage.BASE_URL_DAILY,
            params={
                "format": "rdb",
                "sites": site_no,
                "parameterCd": "00060",
                "statCd": "00003",
                "startDT": start,
                "endDT": end,
            },
            timeout=180,
        )
        if response.status_code in (403, 407) or response.status_code >= 500:
            pytest.skip(f"legacy NWIS unavailable from here: HTTP {response.status_code}")
        response.raise_for_status()
    except (requests.ConnectionError, requests.Timeout) as exc:
        pytest.skip(f"legacy NWIS unreachable from here: {exc}")
    return response.text


@pytest.mark.requires_network
@pytest.mark.parametrize(
    ("site_no", "start", "end", "why"), WATER_YEARS, ids=[w[0] for w in WATER_YEARS]
)
def test_backends_agree_over_a_water_year(site_no: str, start: str, end: str, why: str) -> None:
    legacy = parse_legacy_daily(_legacy_rdb(site_no, start, end), site_no)
    ogc = download_daily(site_no, start_date=start, end_date=end)
    result = daily_mismatches(legacy, ogc)
    assert not result.problems, f"{site_no} {start}..{end} ({why}):\n" + "\n".join(result.problems)


@pytest.mark.requires_network
def test_ice_days_without_a_value_are_absent_from_both() -> None:
    """06214500 WY2026: provisional ``Ice`` days are dropped by both backends."""
    start, end = "2026-01-10", "2026-01-25"
    legacy = parse_legacy_daily(_legacy_rdb("06214500", start, end), "06214500")
    ogc = download_daily("06214500", start_date=start, end_date=end)
    ice = pd.date_range("2026-01-18", "2026-01-25", name="date")
    assert set(legacy.text_values.index) >= set(ice)
    assert not ogc.index.isin(ice).any() and not legacy.frame.index.isin(ice).any()
    assert not daily_mismatches(legacy, ogc).problems


@pytest.mark.requires_network
@pytest.mark.parametrize("site_no", FULL_RECORD_SITES)
def test_download_daily_flow_default_agrees_with_legacy(site_no: str) -> None:
    """The whole record, through the method's own defaults, now on the OGC API."""
    try:
        legacy = USGSgage(site_no).download_daily_flow(backend="nwis-legacy")
    except requests.HTTPError as exc:
        status = exc.response.status_code if exc.response is not None else 0
        if status in (403, 407) or status >= 500:
            pytest.skip(f"legacy NWIS unavailable from here: {exc}")
        raise
    except (requests.ConnectionError, requests.Timeout) as exc:
        pytest.skip(f"legacy NWIS unreachable from here: {exc}")
    new = USGSgage(site_no).download_daily_flow()

    assert list(new.columns) == ["flow_cfs"] and new.index.name == "date"
    assert new.index.equals(legacy.index), (
        f"{len(legacy.index.difference(new.index))} day(s) only in legacy, "
        f"{len(new.index.difference(legacy.index))} only in OGC"
    )
    assert (new["flow_cfs"].to_numpy() == legacy["flow_cfs"].astype(float).to_numpy()).all()


# ----------------------------------------------------------------------------
# Offline: the captured pairs, and the comparison itself
# ----------------------------------------------------------------------------


@pytest.mark.parametrize(("ogc_name", "rdb_name"), CAPTURED_PAIRS)
def test_captured_pairs_agree(ogc_name: str, rdb_name: str) -> None:
    """Both parsers, on the same day's real payloads, return the same record."""
    legacy = parse_legacy_daily(load_rdb(rdb_name))
    ogc = parse_daily_features(load(ogc_name)["features"])
    result = daily_mismatches(legacy, ogc)
    assert not result.problems, "\n".join(result.problems)
    assert not result.ice_estimates_only_in_ogc
    assert len(ogc) > 0


def test_captured_ice_null_days_are_dropped_by_both() -> None:
    legacy = parse_legacy_daily(load_rdb("dv_06214500_2026-01-10_2026-01-25"))
    ogc = parse_daily_features(load("daily_06214500_00060_ice_null")["features"])
    ice = pd.date_range("2026-01-18", "2026-01-25", name="date")
    assert set(legacy.text_values.str.strip()) == {"Ice"}
    assert set(legacy.text_values.index) == set(ice)
    assert not ogc.index.isin(ice).any() and not legacy.frame.index.isin(ice).any()


def test_captured_ice_estimates_are_numeric_on_both() -> None:
    legacy = parse_legacy_daily(load_rdb("dv_12449500_2022-12-15_2022-12-28"))
    ogc = parse_daily_features(load("daily_12449500_00060_ice_estimated")["features"])
    days = pd.date_range("2022-12-20", "2022-12-24", name="date")
    assert (ogc.loc[days, "qualification_code"] == "A:e:ICE").all()
    assert (legacy.codes.loc[days] == "A:e").all()
    assert ogc.loc[days, "flow_cfs"].tolist() == legacy.frame.loc[days, "flow_cfs"].tolist()


def _legacy_from_rows(rows: List[Tuple[str, str, str]]) -> LegacyDaily:
    body = "\n".join(f"USGS\t00000000\t{d}\t{v}\t{c}" for d, v, c in rows)
    text = (
        "# synthetic\n"
        "agency_cd\tsite_no\tdatetime\t1_00060_00003\t1_00060_00003_cd\n"
        "5s\t15s\t20d\t14n\t10s\n" + body + "\n"
    )
    return parse_legacy_daily(text)


def _ogc(rows: List[Tuple[str, float, str]]) -> pd.DataFrame:
    return pd.DataFrame(
        {"flow_cfs": [r[1] for r in rows], "qualification_code": [r[2] for r in rows]},
        index=pd.DatetimeIndex(pd.to_datetime([r[0] for r in rows]), name="date"),
    )


class TestDailyMismatches:
    def test_identical_records_agree(self) -> None:
        legacy = _legacy_from_rows([("2024-01-01", "10.5", "A"), ("2024-01-02", "11", "A:e")])
        ogc = _ogc([("2024-01-01", 10.5, "A"), ("2024-01-02", 11.0, "A:e:ICE")])
        assert daily_mismatches(legacy, ogc).problems == []

    def test_ice_text_day_with_an_ogc_estimate_is_classified_not_failed(self) -> None:
        legacy = _legacy_from_rows([("2024-01-01", "10", "P"), ("2024-01-02", "Ice", "P")])
        ogc = _ogc([("2024-01-01", 10.0, "P"), ("2024-01-02", 9.0, "P:e:ICE")])
        result = daily_mismatches(legacy, ogc)
        assert result.problems == []
        assert result.ice_estimates_only_in_ogc == [pd.Timestamp("2024-01-02")]

    def test_every_kind_of_disagreement_is_reported(self) -> None:
        legacy = _legacy_from_rows(
            [
                ("2024-01-01", "10", "A"),
                ("2024-01-02", "11", "A"),
                ("2024-01-03", "12", "A:e"),
                ("2024-01-04", "Eqp", "P"),
            ]
        )
        ogc = _ogc(
            [
                ("2024-01-01", 10.1, "A"),
                ("2024-01-03", 12.0, "P:e"),
                ("2024-01-04", 5.0, "P"),
                ("2024-01-05", 13.0, "A:Bkw"),
            ]
        )
        report = "\n".join(daily_mismatches(legacy, ogc).problems)
        assert "1 day(s) only in legacy" in report  # 01-02
        assert "2 day(s) only in OGC" in report  # 01-04 (Eqp is not ice), 01-05
        assert "flow_cfs: 1 day(s) differ" in report
        assert "approval 'A' vs 'P'" in report

    def test_unknown_extra_token_is_reported(self) -> None:
        assert "only in OGC ['Bkw']" in str(qualification_mismatch("A", "A:Bkw"))
        assert qualification_mismatch("A:e", "A:EQUIP:e") is None

    def test_legacy_frame_is_the_real_parser_output(self) -> None:
        legacy = _legacy_from_rows([("2024-01-01", "10", "A"), ("2024-01-02", "Ice", "P")])
        assert list(legacy.frame.columns) == ["flow_cfs"]
        assert len(legacy.frame) == 1
        assert legacy.text_values.tolist() == ["Ice"]
