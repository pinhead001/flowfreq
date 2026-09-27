"""Legacy NWIS vs. Water Data OGC API: do both peak backends return the same record?

Issue #29's gate for switching ``peak_sources.DEFAULT_BACKEND``. Both backends
read the same USGS database, so for a given site they must agree on:

- which water years have a peak;
- the peak discharge in each;
- the qualification codes, compared as code *sets* (order is not meaningful),
  after removing legacy's date-precision codes ``Bd``/``Bm`` (day/month
  unknown). The OGC API carries those as ``DAYUNKNOWN``/``MONTHUNKNOWN``
  qualifier tokens, which ``qualifiers_to_codes`` deliberately leaves out of
  the code string -- the placeholder date already encodes them;
- the date, to within one day. The OGC API's ``time`` is the UTC date when the
  time of day is known, while legacy ``peak_dt`` is the local date, so an
  evening peak differs by a day (live-verified 2026-09-25: Big Sandy's WY2023
  peak is 2023-01-04 on the API, 2023-01-03 20:45 CST). Placeholder dates for
  an unknown day or month match exactly on both.

Needs *both* hosts. Legacy NWIS is blocked from Claude Code *web* sessions
(TODO.md, environment constraints) but was reachable from a CLI session on a
developer machine on 2026-09-26, where all three sites below agreed. The legacy
fetch skips -- rather than fails -- when NWIS is unreachable. Run with::

    pytest tests/test_peak_backend_parity.py -m requires_network -v
"""

from __future__ import annotations

from typing import List

import pandas as pd
import pytest
import requests

from flowfreq.peak_codes import parse_codes
from flowfreq.peak_sources import LegacyNwisBackend, WaterDataApiBackend
from flowfreq.usgs import USGSgage

#: Big Sandy: unknown-day historic peaks. Orestimba: zero flows with unknown
#: month (B17C Appendix 10). Potomac at Point of Rocks: a long record with
#: many code types.
PARITY_SITES = ["03606500", "11274500", "01638500"]


def _legacy(site_no: str, via_gage: bool = False) -> pd.DataFrame:
    try:
        if via_gage:
            return USGSgage(site_no).download_peak_flow(backend="nwis-legacy")
        return LegacyNwisBackend().fetch_peaks(site_no)
    except (requests.ConnectionError, requests.Timeout) as exc:
        pytest.skip(f"legacy NWIS unreachable from here: {exc}")
    except requests.HTTPError as exc:
        if exc.response is not None and exc.response.status_code in (403, 407):
            pytest.skip(f"legacy NWIS blocked from here: {exc}")
        raise


#: Legacy date-precision peak codes (NWIS RDB header: "Bd ... Day of occurrence
#: is unknown or not exact", "Bm ... Month ...").
DATE_PRECISION_CODES = frozenset({"Bd", "Bm"})


def _discharge_codes(raw: object) -> frozenset:
    """Codes that describe the discharge, not the date's precision."""
    tokens = [t.strip() for t in str(raw or "").split(",")]
    kept = ",".join(t for t in tokens if t and t not in DATE_PRECISION_CODES)
    return frozenset(parse_codes(kept))


def _mismatches(legacy: pd.DataFrame, ogc: pd.DataFrame) -> List[str]:
    """Every disagreement between the two frames, as readable lines."""
    problems: List[str] = []
    lg = legacy.set_index("water_year")
    og = ogc.set_index("water_year")

    only_legacy = sorted(set(lg.index) - set(og.index))
    only_ogc = sorted(set(og.index) - set(lg.index))
    if only_legacy:
        problems.append(f"water years only in legacy: {only_legacy}")
    if only_ogc:
        problems.append(f"water years only in OGC: {only_ogc}")

    for wy in sorted(set(lg.index) & set(og.index)):
        a, b = lg.loc[wy], og.loc[wy]
        if a["peak_flow_cfs"] != pytest.approx(b["peak_flow_cfs"]):
            problems.append(f"WY{wy} flow: legacy {a['peak_flow_cfs']} vs OGC {b['peak_flow_cfs']}")
        ca = _discharge_codes(a["qualification_code"])
        cb = _discharge_codes(b["qualification_code"])
        if ca != cb:
            problems.append(f"WY{wy} codes: legacy {sorted(ca)} vs OGC {sorted(cb)}")
        da, db = pd.Timestamp(a["peak_date"]), pd.Timestamp(b["peak_date"])
        if pd.isna(da) or pd.isna(db) or abs((da - db).days) > 1:
            problems.append(f"WY{wy} date: legacy {da} vs OGC {db}")
    return problems


@pytest.mark.requires_network
@pytest.mark.parametrize("site_no", PARITY_SITES)
def test_backends_agree(site_no: str) -> None:
    legacy = _legacy(site_no)
    ogc = WaterDataApiBackend().fetch_peaks(site_no)
    problems = _mismatches(legacy, ogc)
    assert not problems, f"site {site_no}:\n" + "\n".join(problems)


@pytest.mark.requires_network
@pytest.mark.parametrize("site_no", PARITY_SITES)
def test_download_peak_flow_default_agrees_with_legacy(site_no: str) -> None:
    """The same gate through ``USGSgage.download_peak_flow``, whose default is now
    the OGC backend: what every existing caller gets must match what it got."""
    legacy = _legacy(site_no, via_gage=True)
    new = USGSgage(site_no).download_peak_flow()
    problems = _mismatches(legacy, new)
    assert not problems, f"site {site_no}:\n" + "\n".join(problems)


class TestMismatchReport:
    """The comparison itself, offline -- so a live pass means what it says."""

    @staticmethod
    def _frame(rows: list) -> pd.DataFrame:
        return pd.DataFrame(
            rows, columns=["water_year", "peak_date", "peak_flow_cfs", "qualification_code"]
        ).assign(peak_date=lambda d: pd.to_datetime(d["peak_date"]))

    def test_identical_frames_agree(self) -> None:
        f = self._frame([[1990, "1990-03-01", 100.0, "7"], [1991, "1991-01-02", 50.0, ""]])
        assert _mismatches(f, f.copy()) == []

    def test_a_one_day_utc_shift_is_tolerated(self) -> None:
        a = self._frame([[2023, "2023-01-03", 2100.0, ""]])
        b = self._frame([[2023, "2023-01-04", 2100.0, ""]])
        assert _mismatches(a, b) == []

    def test_code_order_does_not_matter(self) -> None:
        a = self._frame([[2000, "2000-05-01", 10.0, "2,7"]])
        b = self._frame([[2000, "2000-05-01", 10.0, "7,2"]])
        assert _mismatches(a, b) == []

    def test_date_precision_codes_are_not_discharge_codes(self) -> None:
        a = self._frame([[1897, "1897-03-01", 25000.0, "7,Bd"], [1947, "1947-01-01", 0.0, "Bm"]])
        b = self._frame([[1897, "1897-03-01", 25000.0, "7"], [1947, "1947-01-01", 0.0, ""]])
        assert _mismatches(a, b) == []

    def test_every_kind_of_disagreement_is_reported(self) -> None:
        a = self._frame(
            [
                [1990, "1990-03-01", 100.0, "7"],
                [1991, "1991-01-02", 50.0, ""],
                [1992, "1992-01-02", 60.0, ""],
            ]
        )
        b = self._frame(
            [
                [1990, "1990-03-01", 101.0, "2"],
                [1991, "1991-01-09", 50.0, ""],
                [1993, "1993-01-02", 70.0, ""],
            ]
        )
        report = "\n".join(_mismatches(a, b))
        assert "only in legacy: [1992]" in report
        assert "only in OGC: [1993]" in report
        assert "WY1990 flow" in report
        assert "WY1990 codes" in report
        assert "WY1991 date" in report
