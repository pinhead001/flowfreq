"""The MCP server: flowfreq's report tools, called as a chat client calls them.

Each test talks to the server through an in-process MCP client, so what is checked is what a
model would see: the tool list, the results, and the text of a refusal. Nothing makes a
request: peaks are stubbed where a tool would download them, and StreamStats is faked.
"""

from __future__ import annotations

import asyncio
import json
import time
from typing import Any, Dict, List

import pytest

pytest.importorskip("mcp", reason='the MCP server needs: pip install "flowfreq[mcp]"')

from mcp import Client  # noqa: E402

from flowfreq import __version__  # noqa: E402
from flowfreq import mcp_server  # noqa: E402
from flowfreq.bulletin17c import Bulletin17C  # noqa: E402
from flowfreq.mcp_server import (  # noqa: E402
    INSTRUCTIONS,
    REPORT_DIR_ENV,
    build_server,
    default_report_dir,
    main,
)
from flowfreq.reporting import gaged  # noqa: E402
from flowfreq.reporting import DRAFT_LABEL, sections  # noqa: E402
from flowfreq.reporting.session import ReportSession  # noqa: E402
from tests.test_reporting_engine_check import (  # noqa: E402
    AGREE,
    DISAGREE,
    GAGE_HEADING,
    REFERENCE,
    _as_emafitpr,
)
from tests.test_reporting_sites import GAGE, TENNESSEE, UNUSED  # noqa: E402
from tests.test_reporting_transposed import CITATION, EXPONENTS  # noqa: E402
from tests.test_ungaged import WA_LAT, WA_LON, Fakes, _wa_site  # noqa: E402

TOOLS = [
    "find_gage",
    "report_outline",
    "set_report_details",
    "run_ffa",
    "estimate_ungaged",
    "estimate_from_characteristics",
    "weight_gage_with_regression",
    "transpose_gage_to_site",
    "check_gage_against_reference",
    "low_flow_frequency",
    "flow_duration",
    "monthly_statistics",
    "summary_hydrograph",
    "baseflow",
    "apply_future_change_factors",
    "extend_gage_record",
    "compare_estimates",
    "add_site_summary",
    "add_text_section",
    "edit_narrative",
    "mark_narrative_reviewed",
    "get_section",
    "remove_section",
    "move_section",
    "render_report",
    "save_report",
    "load_report",
]
REGRESSION = {
    "state": "TN",
    "lat": 36.04,
    "lon": -88.23,
    "characteristics": {"DRNAREA": 205.0},
    "region": "HA4 DA only",
}
SKEW = {"regional_skew": GAGE["regional_skew"], "regional_skew_se": GAGE["regional_skew_se"]}


class Chat:
    """A client's side of the conversation with one server."""

    def __init__(self, server: Any) -> None:
        self.server = server

    def _run(self, work: Any) -> Any:
        async def go() -> Any:
            async with Client(self.server) as client:
                return await work(client)

        return asyncio.run(go())

    def result(self, name: str, **arguments: Any) -> Any:
        """The raw result of a tool call."""
        return self._run(lambda client: client.call_tool(name, arguments))

    def call(self, name: str, **arguments: Any) -> Any:
        """What a tool returned; fails the test if the tool refused."""
        result = self.result(name, **arguments)
        assert not result.is_error, result.content[0].text
        return result.structured_content["result"]

    def refusal(self, name: str, **arguments: Any) -> str:
        """The text of a tool's refusal; fails the test if it did not refuse."""
        result = self.result(name, **arguments)
        assert result.is_error, result.content
        return str(result.content[0].text)

    def tools(self) -> List[Any]:
        return list(self._run(lambda client: client.list_tools()).tools)

    def instructions(self) -> str:
        async def read(client: Any) -> Any:
            return client.instructions

        return str(self._run(read))


@pytest.fixture
def session(tmp_path) -> ReportSession:
    return ReportSession("Big Sandy River", output_dir=tmp_path)


@pytest.fixture
def chat(session) -> Chat:
    return Chat(build_server(session))


@pytest.fixture
def peaks(monkeypatch) -> List[str]:
    """Stub the peak download with Big Sandy's record; returns the stations asked for."""
    asked: List[str] = []

    def fake(site_no: str) -> Dict[str, Any]:
        asked.append(site_no)
        return {
            "peaks": {**GAGE["peaks"], "source": "test double", "retrieved": "2026-11-15"},
            "site_name": GAGE["site_name"],
            "drainage_area_sqmi": GAGE["drainage_area_sqmi"],
        }

    monkeypatch.setattr(gaged, "download_peaks", fake)
    return asked


class TestServer:
    def test_tools_and_instructions(self, chat):
        tools = chat.tools()
        assert [tool.name for tool in tools] == TOOLS
        for tool in tools:
            # The description is all a model knows of a tool: none may be empty or indented.
            assert tool.description and len(tool.description) > 40, tool.name
            assert "\n    " not in tool.description, tool.name
        assert chat.instructions() == INSTRUCTIONS
        assert "Every number in the report comes from a tool" in INSTRUCTIONS

    def test_arguments_are_typed(self, chat):
        by_name = {tool.name: tool for tool in chat.tools()}
        schema = by_name["run_ffa"].input_schema
        # site_no is optional since peaks can be given by hand instead.
        assert "required" not in schema
        assert schema["properties"]["skew"]["default"] == "regional"
        assert schema["properties"]["allow_regulated"] == {
            "default": False,
            "title": "Allow Regulated",
            "type": "boolean",
        }
        assert by_name["estimate_from_characteristics"].input_schema["required"] == [
            "state",
            "lat",
            "lon",
            "characteristics",
            "region",
        ]

    def test_version_and_name(self, chat):
        assert chat.server.name == "flowfreq"
        assert chat.server.version == __version__

    def test_a_new_server_has_its_own_empty_report(self, tmp_path):
        chat = Chat(build_server(output_dir=tmp_path))
        assert chat.call("report_outline") == {"title": "Untitled report", "sections": []}
        assert chat.call("add_text_section", title="One", narrative="Text.")["section_id"] == "s1"
        assert Chat(build_server(output_dir=tmp_path)).call("report_outline")["sections"] == []

    def test_figures_are_drawn_without_a_window(self, chat):
        import matplotlib

        # Tools run on worker threads, where an interactive backend's windows are fatal.
        assert matplotlib.get_backend().lower() == "agg"

    def test_calls_sent_together_take_turns(self, chat, session, monkeypatch):
        running: List[int] = []
        most: List[int] = []
        outline = session.outline

        def slow() -> Dict[str, Any]:
            running.append(1)
            most.append(len(running))
            time.sleep(0.05)
            running.pop()
            return outline()

        monkeypatch.setattr(session, "outline", slow)

        async def together(client: Any) -> Any:
            return await asyncio.gather(*(client.call_tool("report_outline", {}) for _ in range(4)))

        results = chat._run(together)
        assert [r.is_error for r in results] == [False] * 4
        assert most == [1, 1, 1, 1]

    def test_default_report_dir(self, monkeypatch, tmp_path):
        monkeypatch.setenv(REPORT_DIR_ENV, str(tmp_path / "reports"))
        assert default_report_dir() == tmp_path / "reports"
        monkeypatch.delenv(REPORT_DIR_ENV)
        monkeypatch.setattr(mcp_server.Path, "home", classmethod(lambda cls: tmp_path))
        assert default_report_dir() == tmp_path / "flowfreq-reports"


class TestFindGage:
    def test_by_number_and_by_name(self, chat):
        found = chat.call("find_gage", query="03606500")
        assert found == {
            "count": 1,
            "gages": [
                {
                    "site_no": "03606500",
                    "site_name": "BIG SANDY RIVER AT BRUCETON, TN",
                    "state": "TN",
                    "latitude": pytest.approx(36.037, abs=1e-3),
                    "longitude": pytest.approx(-88.229, abs=1e-3),
                    "drainage_area_sqmi": 205,
                    "n_peaks": found["gages"][0]["n_peaks"],
                    "first_water_year": 1897,
                    "last_water_year": found["gages"][0]["last_water_year"],
                    "regulation_class": "unknown",
                }
            ],
        }
        assert isinstance(found["gages"][0]["n_peaks"], int)
        named = chat.call("find_gage", query="big sandy bruceton", state="TN", limit=3)
        assert [gage["site_no"] for gage in named["gages"]] == ["03606500"]

    def test_no_match_and_a_refusal(self, chat):
        assert chat.call("find_gage", query="no such river anywhere") == {"count": 0, "gages": []}
        assert chat.refusal("find_gage", query="  ") == (
            "Error executing tool find_gage: ValueError: search needs a station number or "
            "words from a station name"
        )


class TestRunFfa:
    def test_adds_the_section_with_its_peaks(self, chat, session, peaks):
        added = chat.call("run_ffa", site_no=" 03606500 ", stationarity=False, **SKEW)
        assert peaks == ["03606500"]
        assert added["section_id"] == "s1" and added["kind"] == "gaged_ffa"
        assert added["title"] == (
            "Flood frequency at USGS 03606500 Big Sandy River at Bruceton, TN"
        )
        assert "| 1 | 100 | 18,069 |" in added["result"]
        options = session.spec.sections[0].options
        assert options["peaks"]["water_year"] == GAGE["peaks"]["water_year"]
        assert options["stationarity"] is False and options["peaks_table"] is True

    def test_get_section_is_markdown(self, chat, peaks):
        chat.call("run_ffa", site_no="03606500", stationarity=False, peaks_table=False, **SKEW)
        printed = chat.call("get_section", section_id="s1")
        assert printed.startswith("## Flood frequency at USGS 03606500")
        assert "| 1 | 100 | 18,069 |" in printed and "(Figure, in the report: " in printed

    @pytest.mark.parametrize(
        "skew, given, expected",
        [
            ("regional", {}, {"regional_skew_from_site": True}),
            (
                "regional",
                {"allow_provisional_skew": True},
                {"regional_skew_from_site": True, "allow_provisional_skew": True},
            ),
            (
                "regional",
                {"regional_skew": -0.2, "regional_skew_se": 0.4, "regional_skew_citation": "X"},
                {"regional_skew": -0.2, "regional_skew_se": 0.4, "regional_skew_citation": "X"},
            ),
            ("station", {}, {"station_skew_only": True}),
            ("default", {}, {"use_default_skew": True}),
        ],
    )
    def test_skew_choices(self, skew, given, expected):
        arguments = {
            "regional_skew": None,
            "regional_skew_se": None,
            "regional_skew_citation": None,
            "allow_provisional_skew": False,
            **given,
        }
        options = mcp_server._gaged_options("03606500", skew, **arguments)
        assert options == {"site_no": "03606500", **expected}

    def test_station_skew_and_a_threshold(self, chat, session, peaks):
        chat.call(
            "run_ffa",
            site_no="03606500",
            skew="station",
            low_outlier_threshold=1500.0,
            stationarity=False,
            title="At the gage",
            allow_regulated=True,
        )
        section = session.spec.sections[0]
        assert section.title == "At the gage"
        assert section.options["station_skew_only"] is True
        assert section.options["low_outlier_threshold"] == 1500.0
        assert section.options["allow_regulated"] is True

    @pytest.mark.parametrize(
        "arguments, said",
        [
            ({"site_no": "Big Sandy"}, "ReportError: site_no is a USGS station number"),
            ({"site_no": "03606500", "skew": "weighted"}, "skew must be one of regional, station"),
            (
                {"site_no": "03606500", "regional_skew": -0.2},
                "give regional_skew and regional_skew_se together",
            ),
            (
                {"site_no": "03606500", "skew": "station", **SKEW},
                "a regional_skew value goes with skew='regional', not 'station'",
            ),
        ],
    )
    def test_refuses_before_any_download(self, chat, peaks, arguments, said):
        assert said in chat.refusal("run_ffa", **arguments)
        assert peaks == [] and chat.call("report_outline")["sections"] == []

    def test_a_refusal_from_the_analysis_reaches_the_client(self, chat, monkeypatch):
        def regulated(site_no: str) -> Dict[str, Any]:
            from flowfreq.regulation import RegulatedRecordError

            raise RegulatedRecordError(f"USGS {site_no} is regulated")

        monkeypatch.setattr(gaged, "download_peaks", regulated)
        assert chat.refusal("run_ffa", site_no="01019000", skew="station") == (
            "Error executing tool run_ffa: RegulatedRecordError: USGS 01019000 is regulated"
        )

    def test_a_fault_is_reported_with_its_type(self, chat, monkeypatch, caplog):
        def broken(site_no: str) -> Dict[str, Any]:
            raise RuntimeError("the disk is on fire")

        monkeypatch.setattr(gaged, "download_peaks", broken)
        said = chat.refusal("run_ffa", site_no="03606500", skew="station")
        assert said == "Error executing tool run_ffa: RuntimeError: the disk is on fire"
        # A fault, unlike a refusal, leaves its traceback in the server's log.
        (record,) = [r for r in caplog.records if r.name == "flowfreq.mcp_server"]
        assert record.exc_info is not None


class TestEstimates:
    def test_from_characteristics(self, chat, session):
        added = chat.call("estimate_from_characteristics", site_name="Big Sandy", **REGRESSION)
        assert added["kind"] == "ungaged" and "| 1 | 24,552 | 12,155 | 49,589 |" in added["result"]
        assert session.spec.sections[0].options["site_name"] == "Big Sandy"

    def test_from_characteristics_refuses(self, chat):
        said = chat.refusal(
            "estimate_from_characteristics", **{**REGRESSION, "characteristics": {"ELEV": 400.0}}
        )
        assert "DRNAREA" in said
        said = chat.refusal("estimate_from_characteristics", **{**REGRESSION, "region": "nope"})
        assert said.startswith("Error executing tool estimate_from_characteristics: ")
        assert chat.call("report_outline")["sections"] == []

    def test_extrapolation_is_opt_in(self, chat, session):
        big = {**REGRESSION, "characteristics": {"DRNAREA": 90000.0}}
        assert "DRNAREA" in chat.refusal("estimate_from_characteristics", **big)
        chat.call("estimate_from_characteristics", allow_extrapolation=True, **big)
        assert session.spec.sections[0].options["allow_extrapolation"] is True

    def test_ungaged_through_streamstats(self, chat, session, monkeypatch):
        site = _wa_site()
        site["characteristics"][UNUSED] = 5000.0
        fakes = Fakes(site)
        monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", fakes.delineate)
        monkeypatch.setattr("flowfreq.ungaged.locate_regression_regions", fakes.locate)
        added = chat.call(
            "estimate_ungaged", state="WA", lat=WA_LAT, lon=WA_LON, site_name="Above the bridge"
        )
        assert added["kind"] == "ungaged" and added["result"].startswith("Regression flood")
        assert len(fakes.delineate_calls) == 1
        options = session.spec.sections[0].options
        assert options["streamstats"]["polygon_area_sq_mi"] == 82.2
        assert options["characteristics"] == site["characteristics"]

    def test_streamstats_refusing_is_reported(self, chat, monkeypatch):
        import requests

        def refused(*args: Any, **kwargs: Any) -> Any:
            raise requests.HTTPError("403 Client Error: Forbidden for url: streamstats")

        monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", refused)
        said = chat.refusal("estimate_ungaged", state="WA", lat=WA_LAT, lon=WA_LON)
        assert "403 Client Error" in said
        assert chat.call("report_outline")["sections"] == []


class TestTranspose:
    @pytest.fixture
    def gage(self, chat, peaks) -> str:
        added = chat.call("run_ffa", site_no="03606500", stationarity=False, **SKEW)
        return str(added["section_id"])

    def test_exponents_by_aep(self, chat, session, gage):
        added = chat.call(
            "transpose_gage_to_site",
            gage_section_id=gage,
            target_area_sqmi=150.0,
            exponent_citation=CITATION,
            exponent_aeps=EXPONENTS["aeps"],
            exponents=EXPONENTS["exponents"],
            exponent_region="Region 1",
            target_name="Big Sandy River above Bruceton",
        )
        assert added["section_id"] == "s2" and added["number"] == 2
        assert added["kind"] == "transposed"
        assert added["title"] == "Transposed flood frequency at Big Sandy River above Bruceton"
        assert (
            "| 1 | 18,069 | 0.620 | 0.8239 | 14,888 | 11,539 | 22,339 | published |"
            in added["result"]
        )
        # The gage's own section stays, and the options hold what was asked.
        outline = chat.call("report_outline")
        assert [s["kind"] for s in outline["sections"]] == ["gaged_ffa", "transposed"]
        options = session.spec.sections[1].options
        assert options["exponents"] == EXPONENTS and options["detail"] is False
        assert set(options) == {
            "donor",
            "donor_area_sqmi",
            "target_area_sqmi",
            "target_name",
            "exponents",
            "detail",
        }
        text = chat.call("get_section", section_id="s2")
        assert f"- Exponent source: {CITATION}" in text and "- Area ratio: 0.7317" in text
        assert "Regulation was not checked" in text

    def test_one_exponent_and_every_option(self, chat, session, gage):
        added = chat.call(
            "transpose_gage_to_site",
            gage_section_id=gage,
            target_area_sqmi=50.0,
            exponent_citation=CITATION,
            exponent=0.8,
            donor_area_sqmi=180.0,
            max_aep=0.5,
            allow_out_of_range=True,
            donor_site_class="regulated",
            target_site_class="reference",
            allow_regulation_mismatch=True,
            detail=True,
            title="Upstream site",
            position=1,
        )
        assert added["title"] == "Upstream site" and added["number"] == 1
        options = session.spec.sections[0].options
        assert options["exponents"] == {"citation": CITATION, "exponent": 0.8}
        assert options["donor_area_sqmi"] == 180.0 and options["max_aep"] == 0.5
        assert options["allow_out_of_range"] is True
        assert options["allow_regulation_mismatch"] is True
        text = chat.call("get_section", section_id="s2")
        assert "- Area ratio: 0.2778" in text and "| 99 |" not in text.split("### Analysis")[0]
        assert "The area ratio 0.2778 is outside the applicable range" in text
        assert "Regulation mismatch, transferred by override" in text
        assert "### Analysis at the donor gage" in text

    @pytest.mark.parametrize(
        "change, said",
        [
            ({"gage_section_id": "s7"}, "ReportError: no section 's7'"),
            ({}, "ReportError: transposed: exponents takes either exponent"),
            (
                {"exponent": 0.8, "exponents": [0.7]},
                "ReportError: transposed: exponents takes either exponent",
            ),
            (
                {"exponent": 0.8, "exponent_citation": " "},
                "ReportError: transposed: exponents needs a citation",
            ),
            (
                {"exponent": 0.8, "target_area_sqmi": 50.0},
                "ValueError: area ratio 0.2439 is outside the applicable range (0.5, 1.5)",
            ),
            (
                {"exponent": 0.8, "donor_site_class": "regulated", "target_site_class": "urban"},
                "RegulationMismatchError: regulation mismatch: donor is regulated",
            ),
            (
                {"exponent": 0.8, "target_site_class": "dammed"},
                "ReportError: transposed: target_site_class must be one of reference, regulated",
            ),
        ],
    )
    def test_refusals(self, chat, gage, change, said):
        arguments = {
            "gage_section_id": gage,
            "target_area_sqmi": 150.0,
            "exponent_citation": CITATION,
            **change,
        }
        refusal = chat.refusal("transpose_gage_to_site", **arguments)
        assert refusal.startswith(f"Error executing tool transpose_gage_to_site: {said}"), refusal
        assert [s["kind"] for s in chat.call("report_outline")["sections"]] == ["gaged_ffa"]

    def test_the_citation_and_the_areas_are_required(self, chat):
        schema = {t.name: t for t in chat.tools()}["transpose_gage_to_site"].input_schema
        assert schema["required"] == ["gage_section_id", "target_area_sqmi", "exponent_citation"]
        assert schema["properties"]["detail"]["default"] is False


class TestEngineCheck:
    PRINTOUT = "A second program's printout for this gage, table of computed quantiles"

    @pytest.fixture
    def gage(self, chat, peaks) -> str:
        added = chat.call("run_ffa", site_no="03606500", stationarity=False, **SKEW)
        return str(added["section_id"])

    @pytest.fixture
    def emafitpr(self, monkeypatch) -> List[int]:
        """The live peakfq call, answered from stored numbers; the list counts the calls."""
        calls: List[int] = []

        def live(self, *args, **kwargs):
            calls.append(1)
            return _as_emafitpr(REFERENCE)

        monkeypatch.setattr(Bulletin17C, "fortran_reference", live)
        return calls

    @pytest.fixture
    def no_extension(self, monkeypatch) -> None:
        def absent(self, *args, **kwargs):
            raise ImportError("flowfreq.peakfqr is not built")

        monkeypatch.setattr(Bulletin17C, "fortran_reference", absent)

    def test_against_peakfq(self, chat, session, gage, emafitpr):
        added = chat.call("check_gage_against_reference", gage_section_id=gage)
        assert emafitpr == [1]
        assert added["section_id"] == "s2" and added["kind"] == "engine_check"
        assert added["title"] == (
            "Computation check at USGS 03606500 Big Sandy River at Bruceton, TN"
        )
        assert added["result"].startswith("The Bulletin 17C fit at USGS 03606500 Big Sandy")
        assert AGREE in added["result"]
        assert "| 1 | 18,069 | 18,069 | <0.001 | yes |" in added["result"]
        outline = chat.call("report_outline")
        assert [s["kind"] for s in outline["sections"]] == ["gaged_ffa", "engine_check"]
        options = session.spec.sections[1].options
        assert set(options) == {"gage", "reference", "detail"} and options["detail"] is False
        assert options["reference"]["source"].startswith("USGS peakfq 8.1.0")
        assert options["reference"]["quantiles"] == REFERENCE["quantiles"]
        text = chat.call("get_section", section_id="s2")
        assert "90-percent confidence limits from flowfreq and from the reference." in text
        assert GAGE_HEADING not in text

    def test_against_numbers_the_user_gives(self, chat, session, gage, no_extension):
        added = chat.call(
            "check_gage_against_reference",
            gage_section_id=gage,
            reference_source=self.PRINTOUT,
            reference_aeps=[0.1, 0.01],
            reference_flows_cfs=[10600.0, 18500.0],
            reference_lower_cfs=[8900.0, 14000.0],
            reference_upper_cfs=[13100.0, 27100.0],
        )
        assert DISAGREE in added["result"]
        assert "| 1 | 18,069 | 18,500 | 2.330 | no |" in added["result"]
        assert session.spec.sections[1].options["reference"] == {
            "source": self.PRINTOUT,
            "quantiles": [
                {"aep": 0.1, "flow_cfs": 10600.0, "lower_cfs": 8900.0, "upper_cfs": 13100.0},
                {"aep": 0.01, "flow_cfs": 18500.0, "lower_cfs": 14000.0, "upper_cfs": 27100.0},
            ],
        }
        # The same numbers, with the tolerance the user accepts, the gage repeated, and first.
        loose = chat.call(
            "check_gage_against_reference",
            gage_section_id=gage,
            reference_source=self.PRINTOUT,
            reference_aeps=[0.1, 0.01],
            reference_flows_cfs=[10600.0, 18500.0],
            tolerance_pct=3,
            ci_tolerance_pct=5,
            detail=True,
            title="Second computation",
            position=1,
        )
        assert AGREE in loose["result"] and "(tolerance 3 percent)" in loose["result"]
        assert loose["title"] == "Second computation" and loose["number"] == 1
        options = session.spec.sections[0].options
        assert options["tolerance_pct"] == 3 and options["ci_tolerance_pct"] == 5
        assert "lower_cfs" not in options["reference"]["quantiles"][0]
        text = chat.call("get_section", section_id=loose["section_id"])
        assert f"### {GAGE_HEADING}" in text
        assert "confidence limits from flowfreq and from the reference" not in text

    @pytest.mark.parametrize(
        "change, said",
        [
            (
                {},
                "ReportError: engine_check: the reference fit is made by flowfreq's Fortran "
                "extension (USGS peakfq 8.1.0), which is not built in this environment.",
            ),
            ({"gage_section_id": "s7", "reference_aeps": [0.01]}, "ValueError: a reference needs"),
            (
                {"reference_source": "A printout"},
                "ValueError: a reference needs reference_source, reference_aeps and "
                "reference_flows_cfs; with none of the reference_* arguments the record is run "
                "through peakfq",
            ),
            (
                {"reference_aeps": [0.01], "reference_flows_cfs": [18069.0]},
                "ReportError: engine_check: reference needs a source",
            ),
            (
                {
                    "reference_source": "A printout",
                    "reference_aeps": [0.01],
                    "reference_flows_cfs": [18069.0],
                    "reference_lower_cfs": [14000.0],
                },
                "ValueError: give reference_lower_cfs and reference_upper_cfs together",
            ),
            (
                {
                    "reference_source": "A printout",
                    "reference_aeps": [0.1, 0.01],
                    "reference_flows_cfs": [18069.0],
                },
                "ValueError: reference_aeps, reference_flows_cfs must have one value per AEP: 2",
            ),
            (
                {
                    "reference_source": "A printout",
                    "reference_aeps": [1.0],
                    "reference_flows_cfs": [18069.0],
                },
                "ReportError: engine_check: a reference aep is a fraction between 0 and 1",
            ),
            (
                {
                    "gage_section_id": "s7",
                    "reference_source": "A printout",
                    "reference_aeps": [0.01],
                    "reference_flows_cfs": [18069.0],
                },
                "ReportError: no section 's7'",
            ),
        ],
    )
    def test_refusals(self, chat, gage, no_extension, change, said):
        arguments = {"gage_section_id": gage, **change}
        refusal = chat.refusal("check_gage_against_reference", **arguments)
        assert refusal.startswith(
            f"Error executing tool check_gage_against_reference: {said}"
        ), refusal
        assert [s["kind"] for s in chat.call("report_outline")["sections"]] == ["gaged_ffa"]

    def test_only_the_gage_is_required(self, chat):
        schema = {t.name: t for t in chat.tools()}["check_gage_against_reference"].input_schema
        assert schema["required"] == ["gage_section_id"]
        assert schema["properties"]["detail"]["default"] is False
        assert "check_gage_against_reference" in INSTRUCTIONS


class TestWholeReport:
    """A conversation from an empty report to a Word document."""

    def test_gage_regression_weight_narrative_render(self, chat, session, peaks, tmp_path):
        docx = pytest.importorskip("docx")
        chat.call("set_report_details", author="A. Hydrologist", date="2026-11-16")
        chat.call("add_text_section", title="Introduction", narrative="A bridge is planned.")
        gage = chat.call("run_ffa", site_no="03606500", stationarity=False, **SKEW)
        regression = chat.call("estimate_from_characteristics", **REGRESSION)
        weighted = chat.call(
            "weight_gage_with_regression",
            gage_section_id=gage["section_id"],
            regression_section_id=regression["section_id"],
        )
        assert "| 1 | 18,069 | 24,552 | 0.83 | 19,019 | 14,333 | 25,237 |" in weighted["result"]
        assert weighted["section_id"] == "s4"

        edited = chat.call(
            "edit_narrative",
            section_id="s4",
            narrative="The weighted 1-percent flow is 19,019 cfs.",
        )
        assert edited["narrative"] == "draft"
        assert DRAFT_LABEL in chat.call("get_section", section_id="s4")
        outline = chat.call("move_section", section_id="s1", position=2)
        assert [s["section_id"] for s in outline["sections"]] == ["s4", "s1"]
        outline = chat.call("move_section", section_id="s1", position=1)
        assert [(s["number"], s["kind"], s["narrative"]) for s in outline["sections"]] == [
            (1, "text", "draft"),
            (2, "weighted", "draft"),
        ]

        written = chat.call("render_report", filename="big_sandy.docx")
        assert written["sections"] == 2 and written["draft_narratives"] == 2
        text = "\n".join(p.text for p in docx.Document(written["report"]).paragraphs)
        assert "The weighted 1-percent flow is 19,019 cfs." in text and DRAFT_LABEL in text

        # The user reads one narrative; the next rendering labels only the other.
        assert chat.call("mark_narrative_reviewed", section_id="s4")["narrative"] == "reviewed"
        assert chat.call("render_report", filename="big_sandy.md")["draft_narratives"] == 1
        assert (tmp_path / "big_sandy.md").read_text(encoding="utf-8").count(DRAFT_LABEL) == 1

        # The specification beside the report rebuilds it in a later conversation, offline.
        peaks.clear()
        later = Chat(build_server(output_dir=tmp_path))
        outline = later.call("load_report", filename="big_sandy.json")
        assert outline["author"] == "A. Hydrologist" and peaks == []
        assert [(s["kind"], s["narrative"]) for s in outline["sections"]] == [
            ("text", "draft"),
            ("weighted", "reviewed"),
        ]
        assert later.call("get_section", section_id="s2") == chat.call(
            "get_section", section_id="s4"
        )

    def test_save_remove_and_the_refusals(self, chat, tmp_path):
        chat.call("add_text_section", title="One", narrative="Text.")
        chat.call("add_text_section", title="Two", narrative="Text.", position=1)
        saved = chat.call("save_report", filename="draft.json")
        assert saved["sections"] == 2
        assert [
            s["title"] for s in json.loads((tmp_path / "draft.json").read_text())["sections"]
        ] == [
            "Two",
            "One",
        ]
        assert [s["title"] for s in chat.call("remove_section", section_id="s2")["sections"]] == [
            "One"
        ]
        for name, arguments, said in [
            ("remove_section", {"section_id": "s2"}, "no section 's2' in the report"),
            ("get_section", {"section_id": "s9"}, "no section 's9'"),
            ("move_section", {"section_id": "s1", "position": 4}, "position must be from 1 to 1"),
            ("mark_narrative_reviewed", {"section_id": "s9"}, "no section 's9'"),
            (
                "edit_narrative",
                {"section_id": "s1", "narrative": "| a |\n|---|\n| 1 |"},
                "a narrative cannot contain a table",
            ),
            (
                "add_text_section",
                {"title": "Three", "narrative": "# Heading"},
                "a narrative cannot contain a heading",
            ),
            ("render_report", {"filename": "../out.docx"}, "give a file name ending in"),
            ("load_report", {"filename": "missing.json"}, "no specification missing.json"),
            ("set_report_details", {"title": ""}, "a report needs a title"),
            (
                "weight_gage_with_regression",
                {"gage_section_id": "s1", "regression_section_id": "s1"},
                "section s1 is a text section, not a gaged_ffa section",
            ),
        ]:
            refusal = chat.refusal(name, **arguments)
            assert refusal.startswith(f"Error executing tool {name}: ReportError: "), refusal
            assert said in refusal
        assert sorted(p.name for p in tmp_path.iterdir()) == ["draft.json"]

    def test_a_missing_argument_is_an_error_not_a_crash(self, chat):
        result = chat.result("get_section")
        assert result.is_error and "section_id" in result.content[0].text


class TestDailyAndFutureTools:
    """low_flow_frequency, flow_duration and apply_future_change_factors."""

    FACTORS = {
        "factor_set_name": "made-up 2050",
        "scenario": "SSP5-8.5",
        "horizon": "2050",
        "citation": "Example, A.B., 2024, A made-up climate guidance for the tests, table 2.",
        "factor_aeps": [0.5, 0.1, 0.01],
        "factors": [1.05, 1.10, 1.20],
    }

    @pytest.fixture
    def daily_flows(self, monkeypatch) -> List[Any]:
        """Stub the daily download with a synthetic record; returns what was asked for."""
        from flowfreq.reporting import daily
        from tests.test_reporting_daily import STATION

        asked: List[Any] = []

        def fake(site_no: str, start_date: Any = None, end_date: Any = None) -> Dict[str, Any]:
            asked.append((site_no, start_date, end_date))
            return {
                "daily": dict(STATION["daily"]),
                "site_name": STATION["site_name"],
                "drainage_area_sqmi": STATION["drainage_area_sqmi"],
            }

        monkeypatch.setattr(daily, "download_daily", fake)
        return asked

    def test_low_flow_frequency(self, chat, session, daily_flows):
        added = chat.call("low_flow_frequency", site_no=" 01234567 ", non_exceedance=[0.5, 0.1])
        assert daily_flows == [("01234567", None, None)]
        assert added["kind"] == "low_flow"
        assert added["title"] == "Low-flow frequency at USGS 01234567 Test Creek near Nowhere"
        assert "| 7Q10 |" in added["result"] and "| 7Q2 |" in added["result"]
        options = session.spec.sections[0].options
        assert options["n_day"] == 7 and options["daily"]["start"] == "1990-04-01"

    def test_flow_duration(self, chat, daily_flows):
        added = chat.call(
            "flow_duration",
            site_no="01234567",
            exceedance_pct=[10, 50, 90],
            start_date="2000-01-01",
        )
        assert daily_flows == [("01234567", "2000-01-01", None)]
        assert added["kind"] == "flow_duration"
        assert added["result"].count("| q") == 3

    def test_monthly_statistics(self, chat, session, daily_flows):
        added = chat.call("monthly_statistics", site_no="01234567", end_date="2010-09-30")
        assert daily_flows == [("01234567", None, "2010-09-30")]
        assert added["kind"] == "monthly_statistics"
        assert "| October |" in added["result"] and "| September |" in added["result"]
        assert session.spec.sections[0].options["annual_table"] is False

    def test_summary_hydrograph(self, chat, daily_flows):
        added = chat.call("summary_hydrograph", site_no="01234567", start_date="2000-10-01")
        assert daily_flows == [("01234567", "2000-10-01", None)]
        assert added["kind"] == "summary_hydrograph"
        assert "| Oct 15 | 20 |" in added["result"]

    def test_baseflow(self, chat, session, daily_flows):
        added = chat.call(
            "baseflow", site_no="01234567", method="hysep_sliding", drainage_area_sqmi=59.1
        )
        assert daily_flows == [("01234567", None, None)]
        assert added["kind"] == "baseflow"
        assert "| Separation method | hysep_sliding |" in added["result"]
        options = session.spec.sections[0].options
        assert options["drainage_area_sqmi"] == 59.1 and options["annual_table"] is True

    def test_baseflow_takes_the_stations_own_drainage_area(self, chat, session, daily_flows):
        chat.call("baseflow", site_no="01234567")
        options = session.spec.sections[0].options
        assert options["drainage_area_sqmi"] == 42.0
        assert options["method"] == "ih_smoothed_minima"

    def test_baseflow_refuses_an_unknown_method_in_words(self, chat, daily_flows):
        refusal = chat.refusal("baseflow", site_no="01234567", method="eckhardt")
        assert "method must be one of ih_smoothed_minima" in refusal
        assert daily_flows == []

    def test_a_station_number_is_checked(self, chat, daily_flows):
        for name in (
            "low_flow_frequency",
            "flow_duration",
            "monthly_statistics",
            "summary_hydrograph",
            "baseflow",
        ):
            refusal = chat.refusal(name, site_no="Test Creek")
            assert "site_no is a USGS station number, all digits" in refusal
        assert daily_flows == []

    def test_apply_future_change_factors(self, chat, session, peaks):
        chat.call("run_ffa", site_no="03606500", stationarity=False, peaks_table=False, **SKEW)
        added = chat.call("apply_future_change_factors", section_id="s1", **self.FACTORS)
        assert added["section_id"] == "s2" and added["kind"] == "future_flow"
        assert added["title"].startswith("Future-condition flood frequency at USGS 03606500")
        assert "| 1 | 18,069 | 1.200 | 21,683 | published |" in added["result"]
        stored = session.spec.sections[1].options["factor_set"]
        assert stored["level"] == "site" and stored["legal_status"] == "informational"

    def test_factors_are_refused_in_words(self, chat, peaks):
        chat.call("run_ffa", site_no="03606500", stationarity=False, peaks_table=False, **SKEW)
        refusal = chat.refusal(
            "apply_future_change_factors",
            section_id="s1",
            **{**self.FACTORS, "factors": [1.05, 1.10]},
        )
        assert "aeps and factors must be non-empty and equal in length" in refusal
        assert "apply_future_change_factors" in INSTRUCTIONS


class TestMain:
    def test_runs_the_server_on_stdio(self, monkeypatch, tmp_path):
        ran: List[Any] = []

        class Server:
            def run(self, transport: str) -> None:
                ran.append(transport)

        def build(**kwargs: Any) -> Server:
            ran.append(kwargs)
            return Server()

        monkeypatch.setattr(mcp_server, "build_server", build)
        main(["--dir", str(tmp_path)])
        assert ran == [{"output_dir": tmp_path}, "stdio"]

    def test_version(self, capsys):
        with pytest.raises(SystemExit) as exit_:
            main(["--version"])
        assert exit_.value.code == 0
        assert capsys.readouterr().out.strip() == f"flowfreq {__version__}"
