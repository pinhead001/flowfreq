"""Tests for flowfreq.reporting: the section kinds and the built report.

Everything here runs offline. The gaged sections carry their peaks inline (Big Sandy's
systematic record), and the one test of a download replaces ``sections.download_peaks``.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List

import numpy as np
import pandas as pd
import pytest
import requests

from flowfreq.regional_skew import regional_skew_for
from flowfreq.regression.equations import OutOfRangeError
from flowfreq.regulation import RegulatedRecordError
from flowfreq.reporting import (
    DRAFT_LABEL,
    SECTION_KINDS,
    Bullets,
    Figure,
    Heading,
    Narrative,
    Paragraph,
    Report,
    ReportError,
    ReportSpec,
    SectionSpec,
    Table,
    build_report,
    build_section,
    render_markdown,
    sections,
)
from flowfreq.reporting.sections import BULLETIN_17C, cfs, download_peaks, short_citation
from flowfreq.workflow import run_ffa
from tests.fixtures import big_sandy as bs

YEARS = sorted(bs.SYSTEMATIC_PEAKS)
FLOWS = [bs.SYSTEMATIC_PEAKS[y] for y in YEARS]
PEAKS: Dict[str, Any] = {"water_year": YEARS, "peak_flow_cfs": FLOWS}
SKEW = {"regional_skew": bs.REGIONAL_SKEW, "regional_skew_se": bs.REGIONAL_SKEW_SD}
#: A gaged section without the slow or long parts, for tests that are about something else.
FAST: Dict[str, Any] = {"peaks": PEAKS, "stationarity": False, "peaks_table": False}
REGULATED_SITE = "01019000"  # GAGES-II: ~415 days of dam storage
FLORIDA = {
    "state": "FL",
    "characteristics": {"DRNAREA": 105.0, "LC92STOR": 43.3},
    "region": "1",
    "lat": 30.4,
    "lon": -84.3,
}
PNG = b"\x89PNG\r\n\x1a\n"


def _tables(blocks: List[Any]) -> Dict[str, Table]:
    """Captioned tables by caption (number removed)."""
    out = {}
    for block in blocks:
        if isinstance(block, Table) and block.caption:
            out[
                block.caption.split(". ", 1)[-1] if block.caption[:5] == "Table" else block.caption
            ] = block
    return out


def _text(blocks: List[Any]) -> str:
    return "\n".join(b.text for b in blocks if isinstance(b, Paragraph))


@pytest.fixture(scope="module")
def big_sandy_spec() -> ReportSpec:
    return ReportSpec(
        title="Big Sandy River flood frequency",
        subtitle="Example",
        author="A. Analyst",
        project="Demo",
        sections=[
            SectionSpec(
                "text",
                title="Purpose",
                narrative="This report gives **flood quantiles**.\n\n- one\n- two",
                narrative_status="reviewed",
            ),
            SectionSpec(
                "gaged_ffa",
                options={
                    "site_no": "03606500",
                    "site_name": "Big Sandy River at Bruceton, TN",
                    "drainage_area_sqmi": 205.0,
                    "peaks": PEAKS,
                    **SKEW,
                },
                narrative="The 1-percent flow is in the table above.",
            ),
            SectionSpec("ungaged", options=dict(FLORIDA)),
        ],
    )


@pytest.fixture(scope="module")
def big_sandy(big_sandy_spec: ReportSpec) -> Report:
    return build_report(big_sandy_spec)


# -- small helpers ------------------------------------------------------------------------


class TestCfs:
    @pytest.mark.parametrize(
        "value, text",
        [(18069.4, "18,069"), (100, "100"), (99.96, "100"), (12.345, "12.3"), (0.4567, "0.457")],
    )
    def test_formats(self, value, text):
        assert cfs(value) == text

    @pytest.mark.parametrize("value", [None, float("nan"), float("inf")])
    def test_missing_is_a_dash(self, value):
        assert cfs(value) == "-"


class TestShortCitation:
    @pytest.mark.parametrize(
        "citation, short",
        [
            ("Smith, J.A., 2001, A title: U.S. Geological Survey", "Smith (2001)"),
            ("Southard, R.E., and Veilleux, A.G., 2014, Methods", "Southard and Veilleux (2014)"),
            (
                "Eash, D.A., Barnes, K.K., and Veilleux, A.G., 2013, Methods",
                "Eash and others (2013)",
            ),
            ("England, J.F., Jr., 2019a, Guidelines", "England (2019)"),
        ],
    )
    def test_authors_and_year(self, citation, short):
        assert short_citation(citation) == short

    def test_bulletin_17c(self):
        assert short_citation(BULLETIN_17C) == "England and others (2019)"

    def test_another_form_is_returned_whole(self):
        assert (
            short_citation(" Regional skew study, state DOT ") == "Regional skew study, state DOT"
        )


# -- section kinds ------------------------------------------------------------------------


class TestBuildSection:
    def test_kinds_registered(self):
        assert set(SECTION_KINDS) == {
            "text",
            "gaged_ffa",
            "ungaged",
            "weighted",
            "transposed",
            "engine_check",
            "future_flow",
            "comparison",
            "site_summary",
            "low_flow",
            "flow_duration",
        }
        assert all(kind.summary for kind in SECTION_KINDS.values())

    def test_unknown_kind(self):
        with pytest.raises(ReportError, match="unknown section kind 'lowflow'.*gaged_ffa"):
            build_section("lowflow", {})

    def test_unknown_option(self):
        with pytest.raises(ReportError, match="gaged_ffa: unknown option.*regional_skw"):
            build_section("gaged_ffa", {**FAST, "regional_skw": -0.5})

    def test_missing_option(self):
        options = {k: v for k, v in FLORIDA.items() if k != "lat"}
        with pytest.raises(ReportError, match="ungaged: missing option.*lat"):
            build_section("ungaged", options)

    def test_the_callers_options_are_not_changed(self):
        options = {**FAST, "regional_skew_state": "TN"}
        build_section("gaged_ffa", options)
        assert "regional_skew_state" in options and "regional_skew" not in options


class TestGagedFfa:
    def test_quantile_table_is_run_ffa(self, big_sandy):
        direct = run_ffa(
            np.array(FLOWS), np.array(YEARS), site_no="03606500", stationarity=False, **SKEW
        )
        table = _tables(big_sandy.blocks)["Flood quantiles with 90-percent confidence limits."]
        frame = direct["quantile_df"]
        assert len(table.rows) == len(frame) > 0
        for row, (_, expected) in zip(table.rows, frame.iterrows()):
            assert row[2] == cfs(expected["Flow (cfs)"])
            assert row[3] == cfs(expected["Lower 90% CI"])
            assert row[4] == cfs(expected["Upper 90% CI"])
        by_aep = {row[0]: row for row in table.rows}
        assert by_aep["1"][1:3] == ("100", "18,069")

    def test_parameters_table(self, big_sandy):
        rows = dict(
            _tables(big_sandy.blocks)[
                "Fitted log-Pearson Type III parameters and low-outlier screening."
            ].rows
        )
        assert rows["Regional skew"] == "-0.500"
        assert rows["Weighted skew"] == rows["Skew used"] == "-0.281"
        assert rows["Low-outlier threshold (cfs)"] == "none"

    def test_record_table(self, big_sandy):
        rows = dict(_tables(big_sandy.blocks)["Station and peak-flow record."].rows)
        assert rows["USGS station number"] == "03606500"
        assert rows["Drainage area (square miles)"] == "205"
        assert rows["Annual peaks"] == "44"
        assert rows["Water years"] == "1930 to 1973"
        assert rows["Largest peak (cfs)"] == "17,000 (water year 1935)"
        assert rows["Regulation screen"].startswith("unknown")

    def test_figures_and_stationarity(self, big_sandy):
        figures = [b for b in big_sandy.blocks if isinstance(b, Figure)]
        assert [f.caption.split(".")[0] for f in figures[:3]] == [
            "Figure 1",
            "Figure 2",
            "Figure 3",
        ]
        assert all(f.png.startswith(PNG) for f in figures)
        assert Heading(2, "Stationarity screen") in big_sandy.blocks
        assert any("Mann-Kendall" in str(b.rows) for b in big_sandy.blocks if isinstance(b, Table))

    def test_peaks_table_has_only_the_columns_the_record_carries(self, big_sandy):
        table = _tables(big_sandy.blocks)["Annual peak flows used in the analysis."]
        assert table.header == ("Water year", "Peak flow (cfs)")
        assert table.rows[0] == ("1930", "9,100") and len(table.rows) == 44

    def test_dates_and_codes_get_columns(self):
        peaks = {
            **PEAKS,
            "peak_date": [f"{y}-03-01" for y in YEARS],
            "qualification_code": ["2" if y == 1941 else None for y in YEARS],
        }
        content = build_section(
            "gaged_ffa", {"peaks": peaks, "stationarity": False, "apply_peak_codes": False, **SKEW}
        )
        table = _tables(content.blocks)["Annual peak flows used in the analysis."]
        assert table.header == ("Water year", "Date", "Peak flow (cfs)", "Qualification codes")
        assert table.rows[11] == ("1941", "1941-03-01", "1,200", "2")
        assert "qualification codes were not applied" in _text(content.blocks)

    def test_an_uncited_skew_is_said_to_be_uncited(self, big_sandy):
        text = _text(big_sandy.blocks)
        assert "regional skew of -0.500 (standard error 0.550), supplied by the analyst" in text
        assert "without a citation" in text

    def test_station_skew_only(self):
        content = build_section("gaged_ffa", {**FAST, "station_skew_only": True})
        assert "used without regional weighting" in _text(content.blocks)
        assert content.references == [BULLETIN_17C]
        rows = dict(
            _tables(content.blocks)[
                "Fitted log-Pearson Type III parameters and low-outlier screening."
            ].rows
        )
        assert "Regional skew" not in rows and rows["Skew used"] == rows["Station skew"]

    def test_default_skew_is_called_unsourced(self):
        content = build_section("gaged_ffa", {**FAST, "use_default_skew": True})
        assert "no located source" in _text(content.blocks)

    def test_no_skew_choice_is_refused(self):
        with pytest.raises(ValueError, match="No regional skew chosen"):
            build_section("gaged_ffa", dict(FAST))

    def test_skew_looked_up_by_state_is_frozen_and_cited(self):
        found = regional_skew_for("TN")
        content = build_section("gaged_ffa", {**FAST, "regional_skew_state": "TN"})
        options = content.options
        assert "regional_skew_state" not in options
        assert options["regional_skew"] == pytest.approx(0.048)
        assert options["regional_skew_se"] == pytest.approx(math.sqrt(0.16))
        assert options["regional_skew_citation"] == found.citation
        assert options["regional_skew_note"].startswith("Skew Region 2")
        assert (
            "regional skew of 0.048 (standard error 0.400) from Messinger and others (2026): "
            "Skew Region 2" in _text(content.blocks)
        )
        assert content.references == [BULLETIN_17C, found.citation]
        # The frozen options build the same section with no lookup.
        again = build_section("gaged_ffa", options)
        assert again.options == options
        assert [b for b in again.blocks if not isinstance(b, Figure)] == [
            b for b in content.blocks if not isinstance(b, Figure)
        ]

    def test_ambiguous_state_skew_is_a_report_error(self):
        with pytest.raises(ReportError, match="no regional skew: CA has 2 skew regions"):
            build_section("gaged_ffa", {**FAST, "regional_skew_state": "CA"})

    @pytest.mark.parametrize(
        "extra, message",
        [
            ({"regional_skew_state": "TN", "regional_skew": 0.1}, "not both"),
            ({"regional_skew_state": "TN", "regional_skew_from_site": True}, "not both"),
            ({"regional_skew_region": "Skew Region 2"}, "needs regional_skew_state"),
            ({"regional_skew_from_site": True}, "needs site_no"),
        ],
    )
    def test_conflicting_skew_options(self, extra, message):
        with pytest.raises(ReportError, match=message):
            build_section("gaged_ffa", {**FAST, **extra})

    def test_regulated_record_is_refused(self):
        with pytest.raises(RegulatedRecordError, match="allow_regulated=True"):
            build_section(
                "gaged_ffa", {**FAST, "site_no": REGULATED_SITE, "station_skew_only": True}
            )

    def test_regulated_override_is_stated(self):
        content = build_section(
            "gaged_ffa",
            {**FAST, "site_no": REGULATED_SITE, "station_skew_only": True, "allow_regulated": True},
        )
        assert "**This record is classified as regulated.**" in _text(content.blocks)
        rows = dict(_tables(content.blocks)["Station and peak-flow record."].rows)
        assert rows["Regulation screen"].startswith("regulated")

    def test_perception_thresholds_are_stated_and_used(self):
        thresholds = [{"start_year": 1890, "end_year": 1929, "threshold_cfs": 18000}]
        content = build_section("gaged_ffa", {**FAST, **SKEW, "perception_thresholds": thresholds})
        assert "Perception thresholds: 1890 to 1929, 18,000 cfs." in _text(content.blocks)
        plain = build_section("gaged_ffa", {**FAST, **SKEW})
        caption = "Flood quantiles with 90-percent confidence limits."
        assert _tables(content.blocks)[caption].rows != _tables(plain.blocks)[caption].rows

    def test_low_outlier_override_is_stated(self):
        content = build_section("gaged_ffa", {**FAST, **SKEW, "low_outlier_threshold": 2000.0})
        assert "A low-outlier threshold of 2,000 cfs was set by the analyst" in _text(
            content.blocks
        )

    @pytest.mark.parametrize(
        "peaks, message",
        [
            ([1, 2, 3], "peaks must be an object"),
            ({"water_year": [2001], "flow": [10.0]}, "unknown peaks key.*flow"),
            ({"water_year": [], "peak_flow_cfs": []}, "needs water_year and peak_flow_cfs"),
            ({"water_year": [2001, 2002], "peak_flow_cfs": [10.0]}, "must be a list of 2 values"),
            ({"water_year": [2001], "peak_flow_cfs": ["high"]}, "not numeric"),
            ({"water_year": [2001], "peak_flow_cfs": [float("inf")]}, "non-finite"),
        ],
    )
    def test_bad_peaks(self, peaks, message):
        with pytest.raises(ReportError, match=message):
            build_section("gaged_ffa", {"peaks": peaks, "station_skew_only": True})

    def test_bad_perception_threshold(self):
        with pytest.raises(ReportError, match="start_year, end_year, threshold_cfs"):
            build_section(
                "gaged_ffa", {**FAST, **SKEW, "perception_thresholds": [{"start_year": 1890}]}
            )

    def test_neither_site_nor_peaks(self):
        with pytest.raises(ReportError, match="give site_no .* or peaks"):
            build_section("gaged_ffa", {"station_skew_only": True})

    def test_title_without_a_site(self):
        assert build_section("gaged_ffa", {**FAST, **SKEW}).title == "Flood frequency"


class TestDownload:
    def test_downloaded_peaks_are_written_into_the_options(self, monkeypatch):
        calls = []

        def fake(site_no: str) -> Dict[str, Any]:
            calls.append(site_no)
            return {
                "peaks": {**PEAKS, "source": "test double", "retrieved": "2026-11-14"},
                "site_name": "Big Sandy River at Bruceton, TN",
                "drainage_area_sqmi": 205.0,
            }

        monkeypatch.setattr(sections, "download_peaks", fake)
        asked = {"site_no": "03606500", "stationarity": False, "peaks_table": False, **SKEW}
        content = build_section("gaged_ffa", asked)
        assert calls == ["03606500"]
        assert content.options["peaks"]["water_year"] == YEARS
        assert content.options["site_name"] == "Big Sandy River at Bruceton, TN"
        assert content.title == "Flood frequency at USGS 03606500 Big Sandy River at Bruceton, TN"
        rows = dict(_tables(content.blocks)["Station and peak-flow record."].rows)
        assert rows["Data source"] == "test double, retrieved 2026-11-14"

        # Built again from the options as resolved, nothing is downloaded.
        monkeypatch.setattr(sections, "download_peaks", lambda site_no: pytest.fail("downloaded"))
        again = build_section("gaged_ffa", content.options)
        assert again.options == content.options

    def test_a_given_site_name_is_kept(self, monkeypatch):
        monkeypatch.setattr(
            sections,
            "download_peaks",
            lambda site_no: {"peaks": PEAKS, "site_name": "FROM NWIS", "drainage_area_sqmi": None},
        )
        content = build_section(
            "gaged_ffa", {**FAST, "peaks": None, "site_no": "03606500", "site_name": "Mine", **SKEW}
        )
        assert content.options["site_name"] == "Mine"
        assert "drainage_area_sqmi" not in content.options

    def test_download_peaks_shapes_the_record(self, monkeypatch):
        class FakeGage:
            def __init__(self, site_no: str) -> None:
                self.site_no = site_no
                self.site_name = "SOMEWHERE CREEK"
                self.drainage_area = 12.5

            def download_peak_flow(self) -> pd.DataFrame:
                return pd.DataFrame(
                    {
                        "water_year": [2002, 2001, 2003],
                        "peak_flow_cfs": [20.0, 10.0, np.nan],
                        "qualification_code": ["", "6", None],
                        "peak_date": pd.to_datetime(["2002-05-01", "2001-04-02", None]),
                    }
                )

        monkeypatch.setattr("flowfreq.usgs.USGSgage", FakeGage)
        got = download_peaks("01234567")
        assert got["site_name"] == "SOMEWHERE CREEK" and got["drainage_area_sqmi"] == 12.5
        peaks = got["peaks"]
        assert peaks["water_year"] == [2001, 2002]  # sorted, the missing flow dropped
        assert peaks["peak_flow_cfs"] == [10.0, 20.0]
        assert peaks["qualification_code"] == ["6", None]
        assert peaks["peak_date"] == ["2001-04-02", "2002-05-01"]
        assert "01234567" in peaks["source"] and peaks["retrieved"]

    def test_download_failure_is_a_report_error(self, monkeypatch):
        class Down:
            def __init__(self, site_no: str) -> None:
                self.site_no = site_no

            def download_peak_flow(self) -> pd.DataFrame:
                raise requests.ConnectionError("no route")

        monkeypatch.setattr("flowfreq.usgs.USGSgage", Down)
        with pytest.raises(ReportError, match="01234567 could not be downloaded: no route"):
            download_peaks("01234567")

    def test_no_peaks_is_a_report_error(self, monkeypatch):
        class Empty:
            def __init__(self, site_no: str) -> None:
                self.site_no = site_no

            def download_peak_flow(self) -> pd.DataFrame:
                return pd.DataFrame(
                    {
                        "water_year": [],
                        "peak_flow_cfs": [],
                        "qualification_code": [],
                        "peak_date": [],
                    }
                )

        monkeypatch.setattr("flowfreq.usgs.USGSgage", Empty)
        with pytest.raises(ReportError, match="has no annual peaks"):
            download_peaks("01234567")


class TestUngaged:
    def test_florida_estimate(self, big_sandy):
        table = _tables(big_sandy.blocks)["Regression flood quantiles with prediction limits."]
        by_aep = {row[0]: row for row in table.rows}
        assert by_aep["1"][1:4] == ("5,163", "2,130", "12,516")
        assert any(
            isinstance(b, Bullets) and "Method: regression" in b.items for b in big_sandy.blocks
        )
        assert "were not measured by it" in _text(big_sandy.blocks)

    def test_reference_is_the_equations_report(self):
        content = build_section("ungaged", dict(FLORIDA))
        (reference,) = content.references
        assert reference.startswith("Verdi, R.J., and Dixon, J.F., 2011,")
        assert "(2011)" not in reference  # the publication already gives its year
        assert reference.endswith("https://pubs.usgs.gov/publication/sir20115034")

    def test_site_name_titles_the_section(self):
        content = build_section("ungaged", {**FLORIDA, "site_name": "Proposed culvert"})
        assert content.title == "Regression estimate at Proposed culvert"

    def test_out_of_range_passes_through(self):
        options = {**FLORIDA, "characteristics": {"DRNAREA": 1e6, "LC92STOR": 43.3}}
        with pytest.raises(OutOfRangeError, match="DRNAREA"):
            build_section("ungaged", options)

    def test_unknown_state_is_a_lookup_error(self):
        with pytest.raises(LookupError, match="ZZ"):
            build_section("ungaged", {**FLORIDA, "state": "ZZ"})

    @pytest.mark.parametrize(
        "change, message",
        [
            ({"characteristics": {}}, "characteristics must be an object"),
            ({"characteristics": {"DRNAREA": "big"}}, "must be numbers"),
            ({"region": 1}, "region must be a region code"),
            ({"region": None}, "with characteristics, region must be a region code"),
            ({"region": ["1", "2"]}, "with characteristics, region must be a region code"),
            ({"lat": "north"}, "lat and lon must be numbers"),
        ],
    )
    def test_bad_inputs(self, change, message):
        with pytest.raises(ReportError, match=message):
            build_section("ungaged", {**FLORIDA, **change})


# -- the report ---------------------------------------------------------------------------


class TestBuildReport:
    def test_outline(self, big_sandy):
        headings = [b.text for b in big_sandy.blocks if isinstance(b, Heading) and b.level <= 1]
        assert headings == [
            "Big Sandy River flood frequency",
            "1. Purpose",
            "2. Flood frequency at USGS 03606500 Big Sandy River at Bruceton, TN",
            "3. Regression estimate at 30.4000, -84.3000 (FL)",
            "References",
            "About this report",
        ]
        assert big_sandy.blocks[1] == Paragraph("*Example*")
        assert big_sandy.blocks[2].text.startswith("Project: Demo. Prepared by: A. Analyst. Date: ")

    def test_captions_are_numbered_through_the_report(self, big_sandy):
        captions = [b.caption for b in big_sandy.blocks if isinstance(b, Table) and b.caption]
        assert [c.split(".")[0] for c in captions] == [f"Table {i}" for i in range(1, 6)]

    def test_references_are_sorted_and_single(self, big_sandy):
        start = big_sandy.blocks.index(Heading(1, "References"))
        end = big_sandy.blocks.index(Heading(1, "About this report"))
        listed = [b.text for b in big_sandy.blocks[start + 1 : end]]
        assert listed[0] == BULLETIN_17C and listed[1].startswith("Verdi")
        twice = build_report(
            ReportSpec("T", [SectionSpec("ungaged", options=dict(FLORIDA))] * 2, date="2026-01-01")
        )
        assert sum("Verdi" in b.text for b in twice.blocks if isinstance(b, Paragraph)) == 1

    def test_the_date_is_fixed_in_the_spec_as_built(self, big_sandy, big_sandy_spec):
        assert big_sandy_spec.date is None  # the input is not changed
        assert big_sandy.spec.date and big_sandy.spec.date in big_sandy.blocks[2].text
        given = build_report(
            ReportSpec("T", [SectionSpec("text", "A", narrative="x")], date="1 May")
        )
        assert given.blocks[1] == Paragraph("Date: 1 May.")

    def test_rebuilding_the_saved_spec_gives_the_same_report(self, big_sandy, tmp_path):
        path = big_sandy.spec.save(tmp_path / "spec.json")
        again = build_report(ReportSpec.load(path))
        assert again.to_markdown() == big_sandy.to_markdown()
        assert again.figures() == big_sandy.figures()
        assert again.spec == big_sandy.spec

    def test_draft_narrative_is_labeled_and_reviewed_is_not(self, big_sandy):
        narratives = [b for b in big_sandy.blocks if isinstance(b, Narrative)]
        assert [n.status for n in narratives] == ["reviewed", "draft"]
        markdown = big_sandy.to_markdown()
        assert markdown.count(DRAFT_LABEL) == 1
        assert f"> **{DRAFT_LABEL}**\n>\n> The 1-percent flow is in the table above." in markdown
        assert "\nThis report gives **flood quantiles**.\n\n- one\n- two\n" in markdown
        assert "written by a language model and has not been reviewed" in big_sandy.blocks[-1].text

    def test_no_draft_no_disclaimer(self):
        report = build_report(
            ReportSpec("T", [SectionSpec("text", "A", narrative="x", narrative_status="reviewed")])
        )
        assert DRAFT_LABEL not in report.to_markdown()
        assert "language model" not in report.blocks[-1].text
        assert "building the saved specification again reproduces them" in report.blocks[-1].text

    def test_text_section_needs_a_narrative(self):
        with pytest.raises(
            ReportError, match=r"section 1 \(text\): a text section needs a narrative"
        ):
            build_report(ReportSpec("T", [SectionSpec("text", title="Empty", narrative="  ")]))

    def test_text_section_needs_a_title(self):
        with pytest.raises(ReportError, match=r"section 1 \(text\): needs a title"):
            build_report(ReportSpec("T", [SectionSpec("text", narrative="x")]))

    def test_a_table_in_a_narrative_is_refused(self):
        narrative = "The quantiles:\n\n| AEP | Flow |\n|---|---|\n| 1% | 99,999 |\n"
        spec = ReportSpec(
            "T",
            [
                SectionSpec("text", "A", narrative="fine"),
                SectionSpec("ungaged", options=dict(FLORIDA), narrative=narrative),
            ],
        )
        with pytest.raises(
            ReportError, match=r"section 2 \(ungaged\): a narrative cannot contain a table"
        ):
            build_report(spec)

    def test_section_errors_name_the_section(self):
        spec = ReportSpec(
            "T", [SectionSpec("text", "A", narrative="x"), SectionSpec("hydraulics", "B")]
        )
        with pytest.raises(ReportError, match=r"section 2 \(hydraulics\): unknown section kind"):
            build_report(spec)

    def test_a_spec_title_replaces_the_sections_own(self):
        report = build_report(
            ReportSpec("T", [SectionSpec("ungaged", title="Culvert 4", options=dict(FLORIDA))])
        )
        assert Heading(1, "1. Culvert 4") in report.blocks


class TestOutput:
    def test_markdown_file_and_figures(self, big_sandy, tmp_path):
        path = big_sandy.save(tmp_path / "report.md")
        text = path.read_text(encoding="utf-8")
        assert text == big_sandy.to_markdown("report_files")
        assert text.startswith("# Big Sandy River flood frequency\n")
        assert "![Figure 1. Fitted frequency curve" in text
        assert "(report_files/figure_1.png)" in text
        written = sorted(p.name for p in (tmp_path / "report_files").iterdir())
        assert written == ["figure_1.png", "figure_2.png", "figure_3.png"]
        assert (tmp_path / "report_files" / "figure_1.png").read_bytes() == big_sandy.figures()[
            "figure_1.png"
        ]

    def test_markdown_without_figures_makes_no_folder(self, tmp_path):
        report = build_report(ReportSpec("T", [SectionSpec("text", "A", narrative="x")]))
        report.save(tmp_path / "plain.md")
        assert not (tmp_path / "plain_files").exists()

    @pytest.mark.parametrize("name", ["report.pdf", "report"])
    def test_other_suffix_is_refused(self, big_sandy, tmp_path, name):
        with pytest.raises(ReportError, match="use .docx or .md"):
            big_sandy.save(tmp_path / name)

    def test_render_markdown_blocks(self):
        blocks = [
            Heading(0, "T"),
            Heading(2, "Part"),
            Table(("A", "B"), (("x|y", "1"),), ("left", "right"), "Table 1. Cap."),
            Figure(PNG, "Figure 1. Fig."),
            Narrative("Line one.\n\nLine two."),
        ]
        assert render_markdown(blocks, image_dir="") == (
            "# T\n\n### Part\n\nTable 1. Cap.\n\n| A | B |\n|---|---:|\n| x\\|y | 1 |\n\n"
            "![Figure 1. Fig.](figure_1.png)\n\nFigure 1. Fig.\n\n"
            f"> **{DRAFT_LABEL}**\n>\n> Line one.\n>\n> Line two.\n"
        )
