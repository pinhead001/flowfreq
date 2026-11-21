"""Four data-retrieval behaviours from the review of 2026-11-19.

Each was read in the code and is pinned here with a stubbed service:

- the legacy daily parser took the first of several discharge series;
- a peak row with no value was dropped at info level, with no years named;
- ``fetch_nwis_batch`` ran six workers whether or not there was an API key;
- every NSS estimate asked for the national statistic-group list again.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from unittest.mock import patch

import pytest

from flowfreq import streamstats, usgs, waterdata
from flowfreq.usgs import USGSgage, fetch_nwis_batch
from tests.fixtures.nwis_rdb import DV_BASIC
from tests.test_streamstats import (
    _GOOD_CHARACTERISTICS,
    NSS_ESTIMATE_RESPONSE_PFS,
    _mock_response,
    _nss_get_dispatch,
)
from tests.test_usgs import _PEAKS_OGC
from tests.test_usgs import _mock_response as _rdb_response
from tests.test_usgs import _RoutedGet

_FIXTURES = Path(__file__).parent / "fixtures"


def test_legacy_daily_refuses_two_discharge_series() -> None:
    two = DV_BASIC.replace(
        "152418_00060_00003\t152418_00060_00003_cd",
        "152418_00060_00003\t152418_00060_00003_cd\t99_00060_00003\t99_00060_00003_cd",
    )
    lines = []
    for line in two.split("\n"):
        if line.startswith("USGS\t"):
            line += "\t7\tA"
        elif line.startswith("5s\t"):
            line += "\t14n\t10s"
        lines.append(line)
    with patch("flowfreq.usgs.requests.get", return_value=_rdb_response("\n".join(lines))):
        with pytest.raises(ValueError, match="2 daily discharge series.*99_00060_00003"):
            USGSgage("12449500").download_daily_flow(backend="nwis-legacy")
    with patch("flowfreq.usgs.requests.get", return_value=_rdb_response(DV_BASIC)):
        assert len(USGSgage("12449500").download_daily_flow(backend="nwis-legacy")) >= 5


def test_a_peak_with_no_value_is_dropped_out_loud(caplog) -> None:
    data = json.loads((_FIXTURES / "waterdata_peaks_03606500.json").read_text(encoding="utf-8"))
    pages = [copy.deepcopy(data["page1"]), copy.deepcopy(data["page2"])]
    year = pages[0]["features"][1]["properties"]["water_year"]
    pages[0]["features"][1]["properties"]["value"] = None
    with (
        patch("flowfreq.usgs.requests.get", _RoutedGet({_PEAKS_OGC: pages})),
        caplog.at_level("WARNING", logger="flowfreq.peak_sources"),
    ):
        frame = USGSgage("03606500").download_peak_flow(site_metadata=False)
    assert year not in set(frame["water_year"]) and len(frame) == 6
    assert f"water years {year}" in caplog.text and "treated as a gap" in caplog.text
    assert "As in peakfq" in caplog.text


@pytest.mark.parametrize("key, expected", [(None, 2), ("abc", 6)])
def test_batch_workers_follow_the_api_key(monkeypatch, key, expected) -> None:
    monkeypatch.delenv(waterdata.API_KEY_ENV, raising=False)
    waterdata.set_api_key(key)
    seen = {}

    class Pool:
        def __init__(self, max_workers):
            seen["workers"] = max_workers

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

        def submit(self, fn, *args):
            raise AssertionError("no site was given")

    try:
        with patch("concurrent.futures.ThreadPoolExecutor", Pool):
            assert fetch_nwis_batch([]) == ({}, {})
            assert seen["workers"] == expected
            fetch_nwis_batch([], workers=3)
            assert seen["workers"] == 3
        with pytest.raises(ValueError, match="at least 1"):
            fetch_nwis_batch([], workers=0)
    finally:
        waterdata.set_api_key(None)
    assert usgs.fetch_nwis_batch is fetch_nwis_batch


def test_the_national_statistic_groups_are_read_once() -> None:
    with (
        patch("flowfreq.streamstats.requests.get", side_effect=_nss_get_dispatch) as get,
        patch(
            "flowfreq.streamstats.requests.post",
            return_value=_mock_response(NSS_ESTIMATE_RESPONSE_PFS),
        ),
    ):
        for _ in range(3):
            streamstats.estimate_flow_statistics(
                "WA", _GOOD_CHARACTERISTICS, statistic_group_codes=["PFS"]
            )
    national = [c for c in get.call_args_list if c.args[0].endswith("/nssservices/statisticgroups")]
    assert len(national) == 1
