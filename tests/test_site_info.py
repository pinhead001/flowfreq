"""Offline tests for ``USGSgage.fetch_site_info`` on both backends.

Every payload is a trimmed real capture of Methow River at Twisp, WA
(12449500), taken 2026-10-03 (see ``tests/fixtures/waterdata_ogc_responses.py``):

- ``site_catalog_12449500.rdb`` -- the legacy ``seriesCatalogOutput=true``
  listing. Its first ``dv`` row is water temperature (00010) for 2002, which
  is what the period-of-record lookup used to read: a caller bounding a daily
  download with those dates got 165 days of discharge instead of 1919-2026.
- ``site_expanded_12449500.rdb`` -- the legacy ``siteOutput=expanded`` row.
- ``loc_12449500`` / ``meta_12449500_00060`` -- the Water Data API's
  ``monitoring-locations`` feature and discharge ``time-series-metadata``.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Mapping, Optional
from unittest.mock import Mock, patch

import pandas as pd
import pytest
import requests

from flowfreq.peak_sources import DEFAULT_BACKEND
from flowfreq.usgs import USGSgage, _catalog_period, _read_site_rdb
from flowfreq.waterdata import (
    DAILY_MEAN_STATISTIC,
    INSTANTANEOUS_STATISTIC,
    resolve_local_zone,
    series_period_of_record,
    site_attributes,
    site_zone,
)
from tests.fixtures.waterdata_ogc_responses import load, load_rdb

METHOW = "12449500"

#: What both backends report for 12449500 on 2026-10-03, where they agree.
METHOW_COMMON: Dict[str, Any] = {
    "site_name": "METHOW RIVER AT TWISP, WA",
    "drainage_area": 1301.0,
    "state_code": "53",
    "daily_por_start": "1919-06-01",
    "iv_por_start": "1991-04-10",
}


def _response(*, text: str = "", payload: Optional[Any] = None, status: int = 200) -> Mock:
    response = Mock()
    response.status_code = status
    response.text = text
    response.headers = {}
    response.json.return_value = payload
    if status >= 400:
        response.raise_for_status.side_effect = requests.HTTPError(
            f"{status} error", response=response
        )
    else:
        response.raise_for_status.return_value = None
    return response


Handler = Callable[[str, Mapping[str, Any]], Mock]


def _legacy_service(expanded: str, catalog: str) -> Handler:
    """A fake legacy site service: ``siteOutput`` gets ``expanded``, the catalog ``catalog``."""

    def handler(url: str, params: Mapping[str, Any]) -> Mock:
        assert url == USGSgage.BASE_URL_SITE, url
        if params.get("seriesCatalogOutput") == "true":
            return _response(text=catalog)
        assert params.get("siteOutput") == "expanded"
        return _response(text=expanded)

    return handler


def _waterdata(location: Optional[Dict[str, Any]], meta: Optional[Dict[str, Any]]) -> Handler:
    """A fake Water Data API. ``None`` for either collection makes it answer 500."""

    def handler(url: str, params: Mapping[str, Any]) -> Mock:
        assert "api.waterdata.usgs.gov" in url, url
        if "/monitoring-locations/items/" in url:
            assert url.endswith(f"USGS-{METHOW}")
            return _response(payload=location) if location else _response(status=500)
        if "/time-series-metadata/items" in url:
            assert params["parameter_code"] == "00060"
            assert params["monitoring_location_id"] == f"USGS-{METHOW}"
            return _response(payload=meta) if meta else _response(status=500)
        raise AssertionError(f"unexpected request {url}")

    return handler


def _fetch(handler: Handler, **kwargs: Any) -> USGSgage:
    gage = USGSgage(METHOW)

    def get(url: str, params: Optional[Mapping[str, Any]] = None, **_: Any) -> Mock:
        return handler(url, params or {})

    with patch("requests.get", side_effect=get):
        gage.fetch_site_info(use_local_first=False, **kwargs)
    return gage


def _legacy_methow() -> USGSgage:
    return _fetch(
        _legacy_service(load_rdb(f"site_expanded_{METHOW}"), load_rdb(f"site_catalog_{METHOW}")),
        backend="nwis-legacy",
    )


def _waterdata_methow() -> USGSgage:
    return _fetch(
        _waterdata(load(f"loc_{METHOW}"), load(f"meta_{METHOW}_00060")),
        backend="waterdata-ogc",
    )


# ----------------------------------------------------------------------------
# Legacy: the period of record is discharge's, not the first row's
# ----------------------------------------------------------------------------


class TestLegacyCatalogPeriod:
    def test_daily_period_is_discharge_not_water_temperature(self):
        """The regression: 12449500's first dv row is 00010 for 2002-04-18..2002-09-29."""
        gage = _legacy_methow()
        assert gage.daily_por_start == "1919-06-01"
        assert gage.daily_por_end == "2026-10-02"

    def test_iv_period_is_discharge_not_gage_height(self):
        gage = _legacy_methow()
        assert gage.iv_por_start == "1991-04-10"
        assert gage.iv_por_end == "2026-10-03"

    def test_capture_really_leads_with_another_parameter(self):
        """Guard the fixture: without the filter the first dv row is temperature."""
        df = _read_site_rdb(load_rdb(f"site_catalog_{METHOW}"))
        assert df is not None
        first_dv = df[df["data_type_cd"] == "dv"].iloc[0]
        assert (first_dv["parm_cd"], first_dv["begin_date"]) == ("00010", "2002-04-18")

    def test_codes_are_read_as_strings(self):
        """A numeric read would make parm_cd 60 and drop the filter's match."""
        df = _read_site_rdb(load_rdb(f"site_catalog_{METHOW}"))
        assert df is not None
        assert "00060" in set(df["parm_cd"].dropna())

    def test_daily_needs_the_mean_statistic(self):
        """A dv discharge row with another statistic (max, 00001) is not the daily mean."""
        df = _read_site_rdb(load_rdb(f"site_catalog_{METHOW}"))
        assert df is not None
        assert _catalog_period(df, "dv", "00001") == (None, None)
        assert _catalog_period(df, "dv", "00003") == ("1919-06-01", "2026-10-02")

    def test_no_discharge_series_gives_none(self):
        df = _read_site_rdb(load_rdb(f"site_catalog_{METHOW}"))
        assert df is not None
        no_q = df[df["parm_cd"] != "00060"]
        assert _catalog_period(no_q, "dv", "00003") == (None, None)
        assert _catalog_period(no_q, "uv", None) == (None, None)

    def test_several_sensors_give_the_envelope(self):
        df = _read_site_rdb(load_rdb(f"site_catalog_{METHOW}"))
        assert df is not None
        extra = df[(df["data_type_cd"] == "uv") & (df["parm_cd"] == "00060")].copy()
        extra["begin_date"], extra["end_date"] = "1985-01-01", "2000-01-01"
        both = pd.concat([df, extra], ignore_index=True)
        assert _catalog_period(both, "uv", None) == ("1985-01-01", "2026-10-03")

    def test_empty_response_is_none(self):
        assert _read_site_rdb("# nothing\n") is None


class TestLegacyMetadata:
    def test_every_attribute(self):
        gage = _legacy_methow()
        for name, value in METHOW_COMMON.items():
            assert getattr(gage, name) == value, name
        assert gage.latitude == pytest.approx(48.36514507)
        assert gage.longitude == pytest.approx(-120.1161917)
        assert gage.huc == "17020008"
        assert gage._last_api_error is None

    def test_huc_keeps_leading_zero(self):
        """03606500's HUC is 06040005; a numeric read drops the zero."""
        rdb = load_rdb(f"site_expanded_{METHOW}").replace("\t17020008\t", "\t06040005\t")
        gage = _fetch(
            _legacy_service(rdb, load_rdb(f"site_catalog_{METHOW}")), backend="nwis-legacy"
        )
        assert gage.huc == "06040005"

    def test_failure_is_logged_not_raised(self, caplog):
        def down(url: str, params: Mapping[str, Any]) -> Mock:
            return _response(status=503)

        gage = _fetch(down, backend="nwis-legacy")
        assert gage.site_name is None and gage.daily_por_start is None
        assert "503" in (gage._last_api_error or "")
        assert "Period-of-record request failed" in caplog.text


# ----------------------------------------------------------------------------
# Water Data API
# ----------------------------------------------------------------------------


class TestWaterDataSiteInfo:
    def test_is_the_default_backend(self):
        assert DEFAULT_BACKEND == "waterdata-ogc"
        gage = _fetch(_waterdata(load(f"loc_{METHOW}"), load(f"meta_{METHOW}_00060")))
        assert gage.daily_por_start == "1919-06-01"

    def test_every_attribute(self):
        gage = _waterdata_methow()
        for name, value in METHOW_COMMON.items():
            assert getattr(gage, name) == value, name
        assert gage.latitude == pytest.approx(48.3651450711904, abs=1e-12)
        assert gage.longitude == pytest.approx(-120.116191668067, abs=1e-12)
        assert gage.huc == "170200080610"
        # 1919-06-01T07:00Z .. 2026-10-01T07:00Z, local midnights (PDT).
        assert gage.daily_por_end == "2026-10-01"
        # Last reading 2026-10-03T06:15Z is 23:15 PDT on 2026-10-02.
        assert gage.iv_por_end == "2026-10-02"
        assert gage._last_api_error is None

    def test_types_match_the_legacy_backend(self):
        ogc, legacy = _waterdata_methow(), _legacy_methow()
        for name in (
            "site_name",
            "drainage_area",
            "latitude",
            "longitude",
            "huc",
            "state_code",
            "daily_por_start",
            "daily_por_end",
            "iv_por_start",
            "iv_por_end",
        ):
            assert type(getattr(ogc, name)) is type(getattr(legacy, name)), name

    def test_the_peak_series_is_not_a_period_of_record(self):
        """The annual-max series (1920-2025) shares 00060 but is neither period."""
        gage = _waterdata_methow()
        assert "1920-06-16" not in (gage.daily_por_start, gage.iv_por_start)

    def test_local_attributes_are_kept(self):
        gage = USGSgage(METHOW)
        gage.site_name = "Methow at Twisp"
        gage.drainage_area = 1300.0
        handler = _waterdata(load(f"loc_{METHOW}"), load(f"meta_{METHOW}_00060"))

        def get(url: str, params: Optional[Mapping[str, Any]] = None, **_: Any) -> Mock:
            return handler(url, params or {})

        with patch("requests.get", side_effect=get):
            gage.fetch_site_info(use_local_first=False)
        assert gage.site_name == "Methow at Twisp"
        assert gage.drainage_area == 1300.0
        assert gage.latitude is not None

    def test_location_failure_still_sets_periods(self, caplog):
        gage = _fetch(_waterdata(None, load(f"meta_{METHOW}_00060")))
        assert gage.site_name is None and gage.latitude is None
        # No zone: UTC dates, which for a US site are the same days.
        assert gage.daily_por_start == "1919-06-01"
        assert gage._last_api_error
        assert "Monitoring-location request failed" in caplog.text

    def test_metadata_failure_still_sets_location(self):
        gage = _fetch(_waterdata(load(f"loc_{METHOW}"), None))
        assert gage.site_name == METHOW_COMMON["site_name"]
        assert gage.daily_por_start is None and gage.iv_por_start is None
        assert gage._last_api_error

    def test_unknown_backend_raises_before_any_request(self):
        with patch("requests.get", side_effect=AssertionError("no request")):
            with pytest.raises(ValueError, match="Unknown site-information backend"):
                USGSgage(METHOW).fetch_site_info(backend="nwis")


class TestWaterDataHelpers:
    def test_site_attributes_blank_fields_are_none(self):
        feature = {"properties": {"drainage_area": None, "monitoring_location_name": " "}}
        attrs = site_attributes(feature)
        assert attrs == {
            "site_name": None,
            "drainage_area": None,
            "latitude": None,
            "longitude": None,
            "huc": None,
            "state_code": None,
        }

    def test_site_zone_unknown_is_none(self):
        assert site_zone({"properties": {"time_zone_abbreviation": "XYZ"}}) is None
        assert site_zone(load(f"loc_{METHOW}")) == "America/Los_Angeles"

    def test_daily_dates_round_to_the_local_day(self):
        """A local midnight an hour off the zone rule still lands on its own day."""
        meta = {
            "features": [
                {
                    "properties": {
                        "statistic_id": "00003",
                        "primary": "Primary",
                        "begin": "1940-01-01T08:00:00+00:00",  # 00:00 PST
                        # 01:00 PDT: a standard-time midnight on a DST day.
                        "end": "2020-06-01T08:00:00+00:00",
                    }
                }
            ]
        }
        zone = resolve_local_zone("PST", "Y")
        assert series_period_of_record(meta["features"], DAILY_MEAN_STATISTIC, zone) == (
            "1940-01-01",
            "2020-06-01",
        )

    def test_statistic_falls_back_to_computation_identifier(self):
        features = [
            {
                "properties": {
                    "statistic_id": None,
                    "computation_identifier": "Instantaneous",
                    "begin": "2007-10-01T07:00:00+00:00",
                    "end": "2024-09-30T07:00:00+00:00",
                }
            }
        ]
        assert series_period_of_record(features, INSTANTANEOUS_STATISTIC) == (
            "2007-10-01",
            "2024-09-30",
        )
        assert series_period_of_record(features, DAILY_MEAN_STATISTIC) == (None, None)

    def test_primary_series_win(self):
        def series(begin: str, primary: Optional[str]) -> Dict[str, Any]:
            return {
                "properties": {
                    "statistic_id": "00003",
                    "primary": primary,
                    "begin": begin,
                    "end": "2020-01-01T08:00:00+00:00",
                }
            }

        features = [
            series("1950-01-01T08:00:00+00:00", "Primary"),
            series("1900-01-01T08:00:00+00:00", None),
        ]
        assert series_period_of_record(features, DAILY_MEAN_STATISTIC)[0] == "1950-01-01"

    def test_unknown_statistic_raises(self):
        with pytest.raises(ValueError, match="Unsupported statistic_id"):
            series_period_of_record([], "00001")
