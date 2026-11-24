"""Chat-tool arguments the section kinds already took but the tools did not pass on.

From the review of 2026-11-19 (TODO.md, item 5): perception thresholds, peaks typed by hand,
a regional skew by state, a basin outline for the summary's map, and the ungaged options
``aeps``, ``weighting_space`` and ``min_coverage_percent``. Each test calls the tool as a
chat client does and checks what reached the stored section options.
"""

from __future__ import annotations

from typing import Any, Dict

import pytest

pytest.importorskip("mcp", reason='the MCP server needs: pip install "flowfreq[mcp]"')

from flowfreq import mcp_server  # noqa: E402
from tests.test_mcp_server import REGRESSION, SKEW, Chat, chat, peaks, session  # noqa: E402,F401
from tests.test_reporting_sites import GAGE  # noqa: E402
from tests.test_ungaged import WA_LAT, WA_LON, Fakes, _wa_site  # noqa: E402

HAND_PEAKS = {
    "water_year": GAGE["peaks"]["water_year"],
    "peak_flow_cfs": GAGE["peaks"]["peak_flow_cfs"],
}
SQUARE = {
    "type": "Polygon",
    "coordinates": [
        [[-88.30, 36.00], [-88.20, 36.00], [-88.20, 36.10], [-88.30, 36.10], [-88.30, 36.00]]
    ],
}


def _description(chat: Chat, name: str) -> str:
    (tool,) = [t for t in chat.tools() if t.name == name]
    return str(tool.description)


class TestRunFfaPeaksByHand:
    def test_peaks_without_a_station_number_download_nothing(self, chat, session, peaks):
        added = chat.call(
            "run_ffa",
            peaks=HAND_PEAKS,
            site_name="Big Sandy, typed in",
            stationarity=False,
            **SKEW,
        )
        assert peaks == []
        assert added["kind"] == "gaged_ffa"
        assert "| 1 | 100 | 18,069 |" in added["result"]
        options = session.spec.sections[0].options
        assert options["peaks"]["water_year"] == HAND_PEAKS["water_year"]
        assert options["site_name"] == "Big Sandy, typed in"
        assert "site_no" not in options

    def test_peaks_with_a_station_number_are_used_as_given(self, chat, session, peaks):
        chat.call("run_ffa", site_no="03606500", peaks=HAND_PEAKS, stationarity=False, **SKEW)
        assert peaks == []
        assert session.spec.sections[0].options["site_no"] == "03606500"

    def test_neither_station_nor_peaks_is_refused(self, chat, peaks):
        said = chat.refusal("run_ffa", skew="station")
        assert "site_no" in said and "peaks" in said
        assert peaks == [] and chat.call("report_outline")["sections"] == []

    def test_regional_skew_from_the_site_needs_a_station(self, chat, peaks):
        said = chat.refusal("run_ffa", peaks=HAND_PEAKS)
        assert "regional_skew_state" in said
        assert chat.call("report_outline")["sections"] == []

    def test_the_description_says_how(self, chat):
        said = _description(chat, "run_ffa")
        assert "peaks" in said and "water_year" in said and "peak_flow_cfs" in said


class TestRunFfaThresholds:
    def test_thresholds_reach_the_section(self, chat, session, peaks):
        thresholds = [{"start_year": 1897, "end_year": 1929, "threshold_cfs": 30000.0}]
        chat.call(
            "run_ffa",
            site_no="03606500",
            perception_thresholds=thresholds,
            stationarity=False,
            **SKEW,
        )
        assert session.spec.sections[0].options["perception_thresholds"] == thresholds

    def test_a_malformed_threshold_is_refused(self, chat, peaks):
        said = chat.refusal(
            "run_ffa",
            site_no="03606500",
            perception_thresholds=[{"start_year": 1897, "threshold_cfs": 30000.0}],
            stationarity=False,
            **SKEW,
        )
        assert "start_year, end_year, threshold_cfs" in said
        assert chat.call("report_outline")["sections"] == []

    def test_the_description_says_how(self, chat):
        said = _description(chat, "run_ffa")
        assert "perception_thresholds" in said and "threshold_cfs" in said


class TestRunFfaSkewByState:
    def test_state_lookup_options(self):
        options = mcp_server._gaged_options(
            "03606500",
            "regional",
            None,
            None,
            None,
            False,
            regional_skew_state="TN",
        )
        assert options == {"site_no": "03606500", "regional_skew_state": "TN"}

    def test_state_and_region_options(self):
        options = mcp_server._gaged_options(
            "03606500",
            "regional",
            None,
            None,
            None,
            True,
            regional_skew_state="ID",
            regional_skew_region="Region 3",
        )
        assert options == {
            "site_no": "03606500",
            "regional_skew_state": "ID",
            "regional_skew_region": "Region 3",
            "allow_provisional_skew": True,
        }

    def test_skew_from_the_state_is_used(self, chat, session, peaks):
        chat.call("run_ffa", site_no="03606500", regional_skew_state="TN", stationarity=False)
        options = session.spec.sections[0].options
        # The lookup is resolved into the values it found, which is what is stored.
        assert options["regional_skew"] == 0.048
        assert "regional_skew_state" not in options and "regional_skew_from_site" not in options
        assert "2025-5110" in options["regional_skew_citation"]
        assert chat.call("get_section", section_id="s1").count("0.048") >= 1

    @pytest.mark.parametrize(
        "arguments, said",
        [
            ({"skew": "station", "regional_skew_state": "TN"}, "skew='regional'"),
            ({"regional_skew_state": "TN", **SKEW}, "not both"),
            (
                {"regional_skew_region": "Region 3"},
                "regional_skew_region needs regional_skew_state",
            ),
        ],
    )
    def test_refusals(self, chat, peaks, arguments, said):
        assert said in chat.refusal("run_ffa", site_no="03606500", **arguments)
        assert peaks == [] and chat.call("report_outline")["sections"] == []


class TestUngagedOptions:
    def test_aeps_and_weighting_space_from_characteristics(self, chat, session):
        chat.call(
            "estimate_from_characteristics",
            aeps=[0.01, 0.1],
            weighting_space="linear",
            **REGRESSION,
        )
        options = session.spec.sections[0].options
        assert options["aeps"] == [0.01, 0.1]
        assert options["weighting_space"] == "linear"

    def test_a_bad_weighting_space_is_refused(self, chat):
        chat.refusal("estimate_from_characteristics", weighting_space="cubic", **REGRESSION)
        assert chat.call("report_outline")["sections"] == []

    def test_through_streamstats(self, chat, session, monkeypatch):
        fakes = Fakes(_wa_site())
        monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", fakes.delineate)
        monkeypatch.setattr("flowfreq.ungaged.locate_regression_regions", fakes.locate)
        import flowfreq.ungaged as ungaged

        seen: Dict[str, Any] = {}
        real = ungaged.estimate_at_site

        def spy(*args: Any, **kwargs: Any) -> Any:
            seen.update(kwargs)
            return real(*args, **kwargs)

        monkeypatch.setattr(ungaged, "estimate_at_site", spy)
        chat.call(
            "estimate_ungaged",
            state="WA",
            lat=WA_LAT,
            lon=WA_LON,
            aeps=[0.01, 0.1],
            weighting_space="log",
            min_coverage_percent=50.0,
        )
        assert seen["aeps"] == [0.01, 0.1]
        assert seen["weighting_space"] == "log"
        assert seen["min_coverage_percent"] == 50.0
        options = session.spec.sections[0].options
        assert options["min_coverage_percent"] == 50.0
        assert options["aeps"] == [0.01, 0.1] and options["weighting_space"] == "log"

    def test_coverage_below_the_minimum_is_refused(self, chat, monkeypatch):
        fakes = Fakes(_wa_site())
        monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", fakes.delineate)
        monkeypatch.setattr("flowfreq.ungaged.locate_regression_regions", fakes.locate)
        said = chat.refusal(
            "estimate_ungaged", state="WA", lat=WA_LAT, lon=WA_LON, min_coverage_percent=100.1
        )
        assert "min_coverage_percent" in said
        assert chat.call("report_outline")["sections"] == []

    def test_the_descriptions_say_how(self, chat):
        for name in ("estimate_ungaged", "estimate_from_characteristics"):
            said = _description(chat, name)
            assert "aeps" in said and "weighting_space" in said
        assert "min_coverage_percent" in _description(chat, "estimate_ungaged")

    def test_supplied_characteristics_and_the_cautions_reach_the_chat(
        self, chat, session, monkeypatch
    ):
        from tests.test_ungaged_named_region import SUPPLIED, _other

        site = {
            **_wa_site(),
            "streamstats_region": "TX",
            "characteristics": {"DRNAREA": 760.0, "PRECIPCENT": 24.5, "OmegaEM": 0.33},
        }
        fakes = Fakes(site)
        fakes.regions = [_other("TX")]
        monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", fakes.delineate)
        monkeypatch.setattr("flowfreq.ungaged.locate_regression_regions", fakes.locate)
        added = chat.call(
            "estimate_ungaged",
            state="TX",
            lat=29.43,
            lon=-100.0,
            region="Statewide OmegaEM",
            supplied_characteristics=SUPPLIED,
        )
        assert any("contributing drainage area" in text for text in added["cautions"])
        options = session.spec.sections[0].options
        assert options["characteristics"]["MCSLOPE"] == SUPPLIED["MCSLOPE"]
        assert options["streamstats"]["supplied"] == ["CONTDA", "MCSLOPE"]


class TestSummaryBasin:
    def test_a_basin_outline_reaches_the_map(self, chat, session, peaks):
        chat.call("run_ffa", site_no="03606500", stationarity=False, **SKEW)
        chat.call("estimate_from_characteristics", **REGRESSION)
        chat.call(
            "add_site_summary",
            section_ids=["s1", "s2"],
            labels=["Gage", "Regression"],
            site_name="the gage",
            locations=[
                {"name": "Gage", "lat": 36.04, "lon": -88.23, "basin": SQUARE},
            ],
        )
        summary = session.spec.sections[0]
        assert summary.kind == "site_summary"
        assert summary.options["locations"][0]["basin"] == SQUARE

    def test_the_description_says_how(self, chat):
        said = _description(chat, "add_site_summary")
        assert "basin" in said and "GeoJSON" in said
