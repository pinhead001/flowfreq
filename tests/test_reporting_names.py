"""A station or site name is drawn as written: matplotlib must not read it as mathematics."""

from __future__ import annotations

from typing import Any, Dict

import pytest

from flowfreq.reporting import build_section
from tests.test_reporting_daily import STATION
from tests.test_reporting_sites import GAGE
from tests.test_reporting_summary import LOCATIONS, SUMMARY

NAME = "Cash $^$ Creek"

CASES: Dict[str, Dict[str, Any]] = {
    "gaged_ffa": {**GAGE, "site_name": NAME, "stationarity": True},
    "b17b_ffa": {
        **{k: GAGE[k] for k in ("peaks", "regional_skew", "regional_skew_se")},
        "site_name": NAME,
    },
    "low_flow": {**STATION, "site_name": NAME},
    "flow_duration": {**STATION, "site_name": NAME},
    "monthly_statistics": {**STATION, "site_name": NAME},
    "summary_hydrograph": {**STATION, "site_name": NAME},
    "baseflow": {**STATION, "site_name": NAME},
    "site_summary": {**SUMMARY, "site_name": NAME},
}


@pytest.mark.parametrize("kind", CASES)
def test_a_name_with_a_dollar_sign_builds(kind: str) -> None:
    built = build_section(kind, dict(CASES[kind]))
    assert NAME in built.title


def test_so_does_a_location_on_the_summary_map() -> None:
    odd = [{**LOCATIONS[0], "name": NAME}, {**LOCATIONS[1], "name": "_hidden"}]
    built = build_section("site_summary", {**SUMMARY, "locations": odd})
    assert any(NAME in str(getattr(b, "header", "")) for b in built.blocks)
