"""Trimmed live captures of the USGS Water Data OGC API, for ``flowfreq.waterdata`` tests.

Captured 2026-09-25 (UTC timeStamp 2026-09-26) by ``tools/capture_waterdata_fixtures.py``
from ``https://api.waterdata.usgs.gov/ogcapi/v1/collections``. Trimming only removes
content the parser never reads: ``links`` other than ``self``/``next``, the
``thresholds`` block of time-series metadata, and all monitoring-location properties
but a handful (the time-zone fields among them). Feature rows are verbatim.

Files (``tests/fixtures/waterdata_ogc/*.json``):

- ``loc_<site>`` -- monitoring-locations record. 03606500 and 03612600 are
  ``CST``/``Y``; 06214500 and 06191500 are ``MST``/``Y``.
- ``meta_<site>_<param>`` -- time-series-metadata filtered to one parameter.
  03606500 has one instantaneous 00060 series (``6fcce55a...``, begin
  2002-01-01T12:00Z) alongside daily-mean and annual-max series; 03612600 has two
  instantaneous 00065 series, HEADWATER (``15beb942...``, 2013-) and TAILWATER
  (``7ca46507...``, 2023-04-13-), both ``primary``.
- ``iv_03606500_00060_dst`` -- Big Sandy, 2024-11-03 04:00-10:00Z, the CDT->CST
  fall-back night (07:00Z). 12 irregularly spaced rows.
- ``iv_03606500_00060_page{1,2,3}`` -- the same window at ``limit=5``, following
  ``links[rel=next]``.
- ``iv_03606500_00060_empty`` -- 1990-01-01, before the record: 200, zero features.
- ``iv_03612600_00065_unfiltered`` -- both Olmsted stage series interleaved at
  identical timestamps; ``..._tailwater`` the same window with ``time_series_id``.
- ``iv_06214500_00060_ice`` -- Yellowstone at Billings, ``qualifier: ["ICE"]`` with
  null values.
- ``iv_06191500_00060_estimated`` -- one ``["ESTIMATED"]`` approved row.

Daily values, captured 2026-09-27 by ``tools/capture_daily_fixtures.py``. Each
``daily_*.json`` is the ``daily`` collection (00060, statistic 00003) and has a
twin ``dv_<site>_<start>_<end>.rdb``, the legacy NWIS daily-values RDB for the
same window, so the offline tests compare the two parsers on real payloads:

- ``daily_03606500_00060_equip`` / ``dv_03606500_2024-03-20_2024-03-31`` --
  Big Sandy, approved; 2024-03-25..27 are ``["EQUIP", "ESTIMATED"]`` (legacy
  ``A:e``). Rows arrive unsorted.
- ``daily_03606500_00060_page{1,2,3}`` -- the same window at ``limit=5``,
  following ``links[rel=next]``.
- ``daily_12449500_00060_ice_estimated`` / ``dv_12449500_2022-12-15_2022-12-28``
  -- Methow at Twisp; 2022-12-20..24 are ice-affected *estimates*, numeric on
  both (``["ESTIMATED", "ICE"]``, legacy ``A:e``).
- ``daily_06214500_00060_ice_null`` / ``dv_06214500_2026-01-10_2026-01-25`` --
  Yellowstone at Billings, provisional; 2026-01-18..25 have no value: null with
  ``["ICE"]`` on the API, the text ``Ice`` on legacy.
- ``daily_03606500_00060_empty`` -- 1900, before the record: zero features.

Site information, captured live 2026-10-03 for ``tests/test_site_info.py``, all
Methow River at Twisp (12449500):

- ``loc_12449500`` -- monitoring-locations feature, *with* its ``Point``
  geometry, properties trimmed to the ones ``fetch_site_info`` reads plus the
  datums.
- ``meta_12449500_00060`` -- discharge time-series-metadata (annual max,
  instantaneous, daily mean), ``thresholds`` removed.
- ``site_catalog_12449500.rdb`` -- legacy ``seriesCatalogOutput=true``, every
  row except the water-quality (``qw``) ones after the first two. The first
  ``dv`` row is water temperature (00010, 2002-04-18..2002-09-29).
- ``site_expanded_12449500.rdb`` -- legacy ``siteOutput=expanded``, verbatim.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict

FIXTURE_DIR = Path(__file__).resolve().parent / "waterdata_ogc"

BIG_SANDY = "03606500"
BIG_SANDY_IV_00060 = "6fcce55acfa94ec49bf59f175c3b6764"
OLMSTED = "03612600"
OLMSTED_HEADWATER = "15beb94252164ce0b9a77a90edce7528"
OLMSTED_TAILWATER = "7ca465077c26415c93fd0549c697a270"


def load(name: str) -> Dict[str, Any]:
    """Return a fresh copy of the named capture (file stem, no ``.json``)."""
    with open(FIXTURE_DIR / f"{name}.json", encoding="utf-8") as handle:
        data: Dict[str, Any] = json.load(handle)
    return data


def load_rdb(name: str) -> str:
    """Return the named legacy RDB capture (file stem, no ``.rdb``) as text."""
    return (FIXTURE_DIR / f"{name}.rdb").read_text(encoding="utf-8")
