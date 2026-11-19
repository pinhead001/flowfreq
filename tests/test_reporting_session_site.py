"""ReportSession.compare, .summarize and .extend_record, and their flowfreq-mcp tools.

Nothing here makes a request: the gages carry their peaks (Bulletin 17C Appendix 8's
Suwanee Creek and Etowah River, and Big Sandy) and the regression its characteristics.
"""

from __future__ import annotations

import asyncio
from typing import Any, Dict, List

import pytest

pytest.importorskip("mcp")

from mcp import Client  # noqa: E402

from flowfreq.mcp_server import INSTRUCTIONS, build_server  # noqa: E402
from flowfreq.reporting import ReportError  # noqa: E402
from flowfreq.reporting.session import MAIN_TABLES, ReportSession  # noqa: E402
from tests.test_reporting_extension import ETOWAH, SUWANEE  # noqa: E402
from tests.test_reporting_sites import GAGE, TENNESSEE  # noqa: E402
from tests.test_reporting_transposed import EXPONENTS  # noqa: E402

INDEX: Dict[str, Any] = {
    "site_no": "02392000",
    "site_name": "Etowah River at Canton, GA",
    "peaks": ETOWAH,
    "station_skew_only": True,
    "apply_peak_codes": False,
    "stationarity": False,
    "peaks_table": False,
}
LOCATIONS: List[Dict[str, Any]] = [
    {"name": "Gage", "lat": 36.04, "lon": -88.23, "characteristics": {"DRNAREA": 205.0}},
    {"name": "Site", "lat": 36.10, "lon": -88.30, "characteristics": {"DRNAREA": 150.0}},
]


def _ids(session: ReportSession) -> List[str]:
    return [s["section_id"] for s in session.outline()["sections"]]


@pytest.fixture
def session(tmp_path) -> ReportSession:
    made = ReportSession("Big Sandy River", output_dir=tmp_path)
    made.add_section("gaged_ffa", GAGE)  # s1
    made.add_section("ungaged", TENNESSEE)  # s2
    made.transpose("s1", 150.0, EXPONENTS, target_name="the site")  # s3, placed after s1
    return made


class TestCompare:
    def test_sets_the_sections_side_by_side_and_keeps_them(self, session):
        added = session.compare(
            ["s1", "s3", "s2"], labels=["Gage", "Transposed", "Regression"], site_name="the site"
        )
        assert added["kind"] == "comparison" and added["section_id"] == "s4"
        assert added["title"] == "Flood quantiles at the site by 3 methods"
        assert added["result"].startswith(MAIN_TABLES["comparison"])
        assert "| Gage (cfs) | Transposed (cfs) | Regression (cfs) |" in added["result"]
        assert _ids(session) == ["s1", "s3", "s2", "s4"]
        stored = session.spec.sections[3].options
        assert stored["detail"] is False and stored["estimates"][1]["options"]["detail"] is False

    def test_the_sections_options_are_copied_not_shared(self, session):
        session.compare(["s1", "s2"])
        gage = session.spec.sections[0].options
        nested = session.spec.sections[3].options["estimates"][0]["options"]
        assert nested["peaks"] == gage["peaks"] and nested["peaks"] is not gage["peaks"]

    def test_labels_default_to_titles(self, session):
        added = session.compare(["s1", "s2"])
        assert "Flood frequency at USGS 03606500" in added["result"]

    @pytest.mark.parametrize(
        "ids,labels,said",
        [
            (["s1"], None, "two or more estimates"),
            (["s1", "s9"], None, "no section 's9'"),
            (["s1", "s2"], ["only one"], "one label for each section"),
        ],
    )
    def test_refusals_leave_the_report_alone(self, session, ids, labels, said):
        with pytest.raises(ReportError, match=said):
            session.compare(ids, labels=labels)
        assert _ids(session) == ["s1", "s3", "s2"]

    def test_a_text_section_is_not_an_estimate(self, session):
        session.add_text("Notes", "A paragraph.")
        with pytest.raises(ReportError, match="is a text section; an estimate is"):
            session.compare(["s1", "s4"])


class TestSummarize:
    def test_goes_first_with_the_flows_table(self, session):
        added = session.summarize(
            ["s3", "s2"],
            LOCATIONS,
            labels=["Transposed", "Regression"],
            site_name="the site",
            characteristics_source="made up for the tests.",
            title="Executive summary",
        )
        assert added["number"] == 1 and added["title"] == "Executive summary"
        assert added["result"].startswith(MAIN_TABLES["site_summary"])
        assert _ids(session) == ["s4", "s1", "s3", "s2"]
        text = session.section_markdown("s4")
        assert "| Characteristic | Code | Gage | Site |" in text
        assert "Source: made up for the tests." in text

    def test_bad_locations_are_refused(self, session):
        with pytest.raises(ReportError, match="needs lat and lon"):
            session.summarize(["s3", "s2"], [{"name": "Site"}])
        assert len(session.outline()["sections"]) == 3


class TestExtendRecord:
    @pytest.fixture
    def gages(self, tmp_path) -> ReportSession:
        made = ReportSession("Suwanee Creek", output_dir=tmp_path)
        made.add_section("gaged_ffa", {**SUWANEE, "stationarity": False})  # s1
        made.add_section("gaged_ffa", INDEX)  # s2
        return made

    def test_adds_the_extended_fit_after_the_gage(self, gages):
        added = gages.extend_record("s1", "s2")
        assert added["section_id"] == "s3" and added["number"] == 2
        assert added["title"].endswith("record extended")
        assert _ids(gages) == ["s1", "s3", "s2"]
        text = gages.section_markdown("s3")
        assert "| Peaks added | 13 |" in text and "| 0.8519 |" in text
        extension = gages.spec.sections[1].options["extension"]
        assert extension["index"]["site_no"] == "02392000"
        assert extension["index"]["peaks"]["water_year"][0] == 1892

    def test_options_reach_the_extension(self, gages):
        gages.extend_record("s1", "s2", n_extend=8, split_year=1990, title="Extended")
        text = gages.section_markdown("s3")
        assert "| Peaks added | 8 |" in text and "| 1985 to 1989 |" in text

    def test_it_can_be_compared_like_any_gage(self, gages):
        gages.extend_record("s1", "s2")
        added = gages.compare(["s1", "s3"], labels=["Observed", "Extended"])
        assert "| Observed (cfs) | Extended (cfs) |" in added["result"]

    @pytest.mark.parametrize(
        "short,long,said",
        [
            ("s1", "s1", "cannot extend its own record"),
            ("s2", "s1", "could not be extended"),
        ],
    )
    def test_refusals(self, gages, short, long, said):
        with pytest.raises(ReportError, match=said):
            gages.extend_record(short, long)
        assert _ids(gages) == ["s1", "s2"]

    def test_an_extended_record_is_not_extended_again(self, gages):
        gages.extend_record("s1", "s2")
        with pytest.raises(ReportError, match="already an extended record"):
            gages.extend_record("s3", "s2")


class TestTools:
    def _call(self, server: Any, name: str, **arguments: Any) -> Any:
        async def go() -> Any:
            async with Client(server) as client:
                return await client.call_tool(name, arguments)

        return asyncio.run(go())

    def test_the_three_tools_through_a_client(self, session):
        server = build_server(session)
        compared = self._call(
            server,
            "compare_estimates",
            section_ids=["s3", "s2"],
            labels=["Transposed", "Regression"],
        )
        assert not compared.is_error, compared.content
        assert compared.structured_content["result"]["kind"] == "comparison"
        summary = self._call(
            server,
            "add_site_summary",
            section_ids=["s3", "s2"],
            locations=LOCATIONS,
            labels=["Transposed", "Regression"],
            site_name="the site",
        )
        assert not summary.is_error, summary.content
        assert summary.structured_content["result"]["number"] == 1
        refused = self._call(server, "compare_estimates", section_ids=["s1"])
        assert refused.is_error and "two or more estimates" in refused.content[0].text

    def test_extend_gage_record_through_a_client(self, tmp_path):
        made = ReportSession("Suwanee Creek", output_dir=tmp_path)
        made.add_section("gaged_ffa", {**SUWANEE, "stationarity": False})
        made.add_section("gaged_ffa", INDEX)
        result = self._call(
            build_server(made),
            "extend_gage_record",
            gage_section_id="s1",
            index_gage_section_id="s2",
        )
        assert not result.is_error, result.content
        assert result.structured_content["result"]["title"].endswith("record extended")

    def test_the_instructions_name_them(self):
        for name in ("compare_estimates", "add_site_summary", "extend_gage_record"):
            assert name in INSTRUCTIONS

    def test_a_states_factors_are_labelled_as_the_states(self, session):
        result = self._call(
            build_server(session),
            "apply_future_change_factors",
            section_id="s1",
            factor_set_name="made-up",
            scenario="S",
            horizon="2050",
            citation="Example, A.B., 2024, made up for the tests.",
            factor_aeps=[0.5, 0.01],
            factors=[1.05, 1.2],
            level="state",
            geography=["TN"],
        )
        assert not result.is_error, result.content
        section_id = result.structured_content["result"]["section_id"]
        assert "state-level" in session.section_markdown(section_id)
