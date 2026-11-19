"""Faults a review of 2026-11-19 reproduced, each with the input that showed it."""

from __future__ import annotations

from typing import Any, Dict, List
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest
import requests

from flowfreq import waterdata
from flowfreq.lowflow import LowFlowFrequency
from flowfreq.regime import flow_duration_curve
from flowfreq.ungaged import assemble
from tests.test_region_limits import KW, _located


def _daily(values: np.ndarray, start: str = "1990-04-01") -> pd.DataFrame:
    return pd.DataFrame(
        {"flow_cfs": values}, index=pd.date_range(start, periods=len(values), freq="D")
    )


class TestWaterDataPaging:
    """A page that is not a feature collection is an error, not the end of the record."""

    def _pages(self, pages: List[Dict[str, Any]]):
        return patch.object(waterdata, "_get_json", side_effect=pages)

    def test_an_error_body_with_http_200_is_not_the_last_page(self):
        first = {
            "features": [{"id": 1}, {"id": 2}],
            "links": [{"rel": "next", "href": "https://example.invalid/next"}],
        }
        with self._pages([first, {"code": "ServerError", "description": "upstream timeout"}]):
            with pytest.raises(requests.RequestException, match="page 2 is not a feature"):
                waterdata._get_all_features("https://example.invalid/items", {}, 30, "daily")

    def test_fewer_rows_than_the_service_counted_is_an_error(self):
        short = {"features": [{"id": 1}], "numberMatched": 366, "links": []}
        with self._pages([short]):
            with pytest.raises(requests.RequestException, match="1 of 366 rows"):
                waterdata._get_all_features("https://example.invalid/items", {}, 30, "daily")

    def test_a_whole_record_and_an_empty_one_still_pass(self):
        whole = {"features": [{"id": 1}, {"id": 2}], "numberMatched": 2, "links": []}
        with self._pages([whole]):
            assert len(waterdata._get_all_features("u", {}, 30, "daily")) == 2
        with self._pages([{"features": [], "links": []}]):
            assert waterdata._get_all_features("u", {}, 30, "daily") == []


def test_flow_duration_treats_negative_flows_as_missing(caplog):
    values = np.linspace(10.0, 1000.0, 1000)
    values[:30] = -999999.0
    curve = flow_duration_curve(_daily(values))
    assert curve["flow_cfs"].min() > 0
    clean = flow_duration_curve(_daily(values[30:]))
    assert list(curve["flow_cfs"]) == pytest.approx(list(clean["flow_cfs"]))
    assert "30 negative flow(s) treated as missing" in caplog.text
    with pytest.raises(ValueError, match="no non-null flows"):
        flow_duration_curve(_daily(np.full(10, -1.0)))


@pytest.mark.parametrize("confidence,pct", [(0.58, 58), (0.975, 98), (0.90, 90), (0.95, 95)])
def test_low_flow_limit_columns_are_named_for_the_level_given(confidence, pct):
    rng = np.random.default_rng(5)
    days = pd.date_range("1990-04-01", "2015-03-31", freq="D")
    flows = (
        100.0 + 50.0 * np.sin(2 * np.pi * days.dayofyear / 365.25) + rng.gamma(2.0, 5.0, len(days))
    )
    fit = LowFlowFrequency(pd.DataFrame({"flow_cfs": flows}, index=days))
    fit.run_analysis()
    limits = fit.compute_confidence_limits(confidence=confidence)
    assert f"lower_{pct}pct" in limits.columns and f"upper_{pct}pct" in limits.columns


class TestCoverage:
    """A basin only partly inside the regions used says so."""

    def test_partial_coverage_is_recorded(self):
        from flowfreq.regression.library import load_state

        lib = load_state("TN")
        values = {"DRNAREA": 205.0}
        part = assemble(lib, _located(("GC2022", 40.0)), values, **KW)
        whole = assemble(lib, _located(("GC2022", 100.0)), values, **KW)
        assert part.covered_percent == pytest.approx(40.0)
        assert whole.covered_percent == pytest.approx(100.0)
        # The flow is the same: the weight is rescaled. The record of it is what was missing.
        assert list(part.quantiles["flow_cfs"]) == pytest.approx(list(whole.quantiles["flow_cfs"]))

    def test_an_estimate_cautions_when_the_regions_cover_part_of_the_basin(self, monkeypatch):
        from dataclasses import replace

        from flowfreq.ungaged import estimate_at_site
        from tests.test_ungaged import WA_LAT, WA_LON, Fakes, _wa_site

        fakes = Fakes(_wa_site())
        monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", fakes.delineate)
        monkeypatch.setattr("flowfreq.ungaged.locate_regression_regions", fakes.locate)
        whole = estimate_at_site(WA_LAT, WA_LON, state="WA")
        assert whole.provenance["covered_percent"] == pytest.approx(100.0)
        assert "cover only" not in whole.advisory
        fakes.regions = [replace(r, percent_weight=0.6 * r.percent_weight) for r in fakes.regions]
        part = estimate_at_site(WA_LAT, WA_LON, state="WA")
        assert part.provenance["covered_percent"] == pytest.approx(60.0)
        assert "cover only 60 percent of the basin" in part.advisory
        assert "cover only 60 percent" in part.to_markdown()
