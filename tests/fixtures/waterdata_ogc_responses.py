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
