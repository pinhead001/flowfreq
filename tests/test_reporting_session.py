"""`ReportSession`: a report built a step at a time, as a chat builds one.

Nothing here makes a request: the gage section carries its peaks and the regression section
its characteristics.
"""

from __future__ import annotations

import json
from typing import Any, Dict

import pytest

from flowfreq.reporting import DRAFT_LABEL, ReportError, ReportSpec, build_section
from flowfreq.reporting.session import LONG_TABLE_ROWS, MAIN_TABLES, ReportSession
from tests.test_reporting_sites import GAGE, TENNESSEE

#: The same equations at a smaller basin, for a second section that is not the first.
UPSTREAM: Dict[str, Any] = {**TENNESSEE, "characteristics": {"DRNAREA": 120.0}}
TEXT = "The basin is rural.\n\n- Mostly forest\n- No large dams"


@pytest.fixture
def session(tmp_path) -> ReportSession:
    return ReportSession("Big Sandy River", output_dir=tmp_path)


def _ids(session: ReportSession):
    return [section["section_id"] for section in session.outline()["sections"]]


class TestDetails:
    def test_new_session(self, tmp_path):
        session = ReportSession(output_dir=tmp_path)
        assert session.outline() == {"title": "Untitled report", "sections": []}
        assert session.output_dir == tmp_path

    def test_set_details_keeps_what_is_not_given(self, session):
        session.set_details(subtitle="Flood frequency", author="A. Hydrologist")
        outline = session.set_details(project="Bridge 12", date="2026-11-16")
        assert outline == {
            "title": "Big Sandy River",
            "subtitle": "Flood frequency",
            "author": "A. Hydrologist",
            "project": "Bridge 12",
            "date": "2026-11-16",
            "sections": [],
        }
        # An empty string clears a detail; the title cannot be cleared.
        assert "subtitle" not in session.set_details(subtitle=" ")
        with pytest.raises(ReportError, match="needs a title"):
            session.set_details(title="  ")
        assert session.outline()["title"] == "Big Sandy River"

    def test_spec_is_a_copy(self, session):
        session.add_section("ungaged", TENNESSEE)
        spec = session.spec
        assert isinstance(spec, ReportSpec) and spec.title == "Big Sandy River"
        spec.sections[0].options["characteristics"]["DRNAREA"] = 1.0
        spec.sections.clear()
        assert session.spec.sections[0].options["characteristics"] == {"DRNAREA": 205.0}


class TestAddSection:
    def test_returns_the_id_and_the_main_table(self, session):
        added = session.add_section("ungaged", TENNESSEE)
        assert added["section_id"] == "s1" and added["number"] == 1
        assert added["kind"] == "ungaged" and added["narrative"] == "none"
        assert added["title"] == "Regression estimate at 36.0400, -88.2300 (TN)"
        assert added["result"].startswith(MAIN_TABLES["ungaged"])
        assert "| 1 | 24,552 | 12,155 | 49,589 |" in added["result"]

    def test_gage_section(self, session):
        added = session.add_section("gaged_ffa", GAGE, title="At the gage")
        assert added["title"] == "At the gage"
        assert "| 1 | 100 | 18,069 |" in added["result"]

    def test_options_are_stored_as_resolved(self, session):
        session.add_section("ungaged", TENNESSEE)
        assert session.spec.sections[0].options == build_section("ungaged", TENNESSEE).options

    def test_position_and_ids(self, session):
        session.add_section("ungaged", TENNESSEE)
        session.add_text("Summary", TEXT)
        first = session.add_section("ungaged", UPSTREAM, position=1)
        assert first["section_id"] == "s3" and first["number"] == 1
        assert _ids(session) == ["s3", "s1", "s2"]
        # An id is not issued twice, so an old one cannot come to mean a new section.
        session.remove_section("s3")
        assert session.add_text("Later", TEXT, position=3)["section_id"] == "s4"
        assert _ids(session) == ["s1", "s2", "s4"]

    @pytest.mark.parametrize(
        "kind, options, position, match",
        [
            ("ungaged", TENNESSEE, 3, "position must be from 1 to 2"),
            ("ungaged", TENNESSEE, 0, "position must be from 1 to 2"),
            ("flood", {}, None, "unknown section kind"),
            ("ungaged", {"state": "TN"}, None, "missing option"),
            ("text", {}, None, "a text section needs a narrative"),
        ],
    )
    def test_a_refused_section_changes_nothing(self, session, kind, options, position, match):
        session.add_section("ungaged", TENNESSEE)
        with pytest.raises(ReportError, match=match):
            session.add_section(kind, options, position=position)
        assert _ids(session) == ["s1"]
        # The refused section took no id either.
        assert session.add_text("Next", TEXT)["section_id"] == "s2"


class TestText:
    def test_add_text(self, session):
        added = session.add_text(" Basin ", TEXT)
        assert added == {
            "section_id": "s1",
            "number": 1,
            "kind": "text",
            "title": "Basin",
            "narrative": "draft",
        }

    @pytest.mark.parametrize(
        "title, narrative, match",
        [
            (" ", TEXT, "needs a title"),
            ("Basin", "  ", "needs a narrative"),
            ("Basin", "| a | b |\n|---|---|\n| 1 | 2 |", "cannot contain a table"),
            ("Basin", "## Setting\n\nRural.", "cannot contain a heading"),
        ],
    )
    def test_refuses(self, session, title, narrative, match):
        with pytest.raises(ReportError, match=match):
            session.add_text(title, narrative)
        assert _ids(session) == []


class TestNarrative:
    def test_draft_then_reviewed_then_draft_again(self, session):
        session.add_section("ungaged", TENNESSEE)
        assert session.set_narrative("s1", "The 1-percent flow is large.")["narrative"] == "draft"
        assert DRAFT_LABEL in session.section_markdown("s1")
        assert session.mark_reviewed("s1")["narrative"] == "reviewed"
        printed = session.section_markdown("s1")
        assert DRAFT_LABEL not in printed and "The 1-percent flow is large." in printed
        # New text has not been read, whatever the old text's status.
        assert session.set_narrative("s1", "Rewritten.")["narrative"] == "draft"
        assert session.spec.sections[0].narrative == "Rewritten."

    def test_empty_text_removes_it(self, session):
        session.add_section("ungaged", TENNESSEE)
        session.set_narrative("s1", "Some text.")
        assert session.set_narrative("s1", "  ")["narrative"] == "none"
        assert session.spec.sections[0].narrative is None

    def test_refuses(self, session):
        session.add_section("ungaged", TENNESSEE)
        session.add_text("Basin", TEXT)
        with pytest.raises(ReportError, match="cannot contain a table"):
            session.set_narrative("s1", "| a |\n|---|\n| 1 |")
        with pytest.raises(ReportError, match="remove the section instead"):
            session.set_narrative("s2", "")
        with pytest.raises(ReportError, match="no narrative to review"):
            session.mark_reviewed("s1")
        with pytest.raises(ReportError, match=r"no section 's9'.*its sections: s1, s2"):
            session.set_narrative("s9", "Text.")
        assert session.spec.sections[0].narrative is None
        assert session.spec.sections[1].narrative == TEXT

    def test_unknown_section_in_an_empty_report(self, session):
        with pytest.raises(ReportError, match="its sections: none yet"):
            session.mark_reviewed("s1")


class TestSectionMarkdown:
    def test_figures_and_long_tables_are_named(self, session):
        session.add_section("gaged_ffa", {**GAGE, "peaks_table": True})
        printed = session.section_markdown("s1")
        assert printed.startswith(
            "## Flood frequency at USGS 03606500 Big Sandy River at Bruceton, TN\n"
        )
        assert "| 1 | 100 | 18,069 |" in printed
        assert "(Figure, in the report: " in printed and "![" not in printed
        n_peaks = len(GAGE["peaks"]["water_year"])
        assert n_peaks > LONG_TABLE_ROWS
        assert f"(Table of {n_peaks} rows, in the report: " in printed
        assert f"| {GAGE['peaks']['water_year'][0]} |" not in printed

    def test_a_short_peaks_table_is_printed(self, session):
        short = {key: values[:20] for key, values in GAGE["peaks"].items()}
        session.add_section("gaged_ffa", {**GAGE, "peaks": short, "peaks_table": True})
        printed = session.section_markdown("s1")
        assert "(Table of" not in printed and f"| {short['water_year'][0]} |" in printed

    def test_unknown_section(self, session):
        with pytest.raises(ReportError, match="no section 's1'"):
            session.section_markdown("s1")


class TestWeight:
    @pytest.fixture
    def parts(self, session) -> ReportSession:
        session.add_text("Introduction", TEXT)
        session.add_section("gaged_ffa", GAGE)
        session.add_section("ungaged", TENNESSEE)
        return session

    def test_takes_the_place_of_its_parts(self, parts):
        weighted = parts.weight("s2", "s3")
        assert weighted["section_id"] == "s4" and weighted["number"] == 2
        assert weighted["kind"] == "weighted"
        assert weighted["result"].startswith(MAIN_TABLES["weighted"])
        assert "| 1 | 18,069 | 24,552 | 0.83 | 19,019 | 14,333 | 25,237 |" in weighted["result"]
        assert _ids(parts) == ["s1", "s4"]
        options = parts.spec.sections[1].options
        assert options["gage"]["site_no"] == "03606500" and options["detail"] is True
        assert "### At-site analysis" in parts.section_markdown("s4")

    def test_keep_parts_and_no_detail(self, parts):
        weighted = parts.weight("s2", "s3", detail=False, keep_parts=True)
        assert weighted["number"] == 2 and _ids(parts) == ["s1", "s4", "s2", "s3"]
        assert "At-site analysis" not in parts.section_markdown("s4")

    def test_a_part_with_a_narrative_is_not_dropped(self, parts):
        parts.set_narrative("s3", "The regression is the Tennessee equation.")
        with pytest.raises(ReportError, match="section s3 has a narrative"):
            parts.weight("s2", "s3")
        assert _ids(parts) == ["s1", "s2", "s3"]
        assert parts.weight("s2", "s3", keep_parts=True)["section_id"] == "s4"

    @pytest.mark.parametrize(
        "gage, regression, match",
        [
            ("s3", "s2", "section s3 is a ungaged section, not a gaged_ffa section"),
            ("s2", "s1", "section s1 is a text section, not a ungaged section"),
            ("s2", "s7", "no section 's7'"),
        ],
    )
    def test_refuses(self, parts, gage, regression, match):
        with pytest.raises(ReportError, match=match):
            parts.weight(gage, regression)
        assert _ids(parts) == ["s1", "s2", "s3"]

    def test_a_refused_weighting_changes_nothing(self, parts):
        # Nevada's region 10 publishes no variance of prediction, so nothing can be weighted.
        parts.add_section(
            "ungaged",
            {
                "state": "NV",
                "characteristics": {"DRNAREA": 205.0},
                "region": "10",
                "lat": 36.2,
                "lon": -115.1,
            },
        )
        with pytest.raises(ReportError, match="no AEP could be weighted"):
            parts.weight("s2", "s4")
        assert _ids(parts) == ["s1", "s2", "s3", "s4"]


class TestOrder:
    def test_move_and_remove(self, session):
        for title in ("One", "Two", "Three"):
            session.add_text(title, TEXT)
        outline = session.move_section("s3", 1)
        assert [s["title"] for s in outline["sections"]] == ["Three", "One", "Two"]
        assert [s["number"] for s in outline["sections"]] == [1, 2, 3]
        assert _ids(session) == ["s3", "s1", "s2"]
        session.move_section("s3", 3)
        assert _ids(session) == ["s1", "s2", "s3"]
        assert [s["section_id"] for s in session.remove_section("s2")["sections"]] == ["s1", "s3"]

    def test_refuses(self, session):
        session.add_text("One", TEXT)
        with pytest.raises(ReportError, match="position must be from 1 to 1, got 2"):
            session.move_section("s1", 2)
        with pytest.raises(ReportError, match="no section 's2'"):
            session.move_section("s2", 1)
        with pytest.raises(ReportError, match="no section 's2'"):
            session.remove_section("s2")
        assert _ids(session) == ["s1"]


class TestFiles:
    @pytest.fixture
    def report(self, session) -> ReportSession:
        session.set_details(author="A. Hydrologist", date="2026-11-16")
        session.add_text("Introduction", TEXT)
        session.add_section("ungaged", TENNESSEE)
        session.set_narrative("s2", "The regression is the Tennessee equation.")
        session.mark_reviewed("s2")
        return session

    def test_render_markdown_and_the_specification(self, report, tmp_path):
        written = report.render("big_sandy.md")
        assert written == {
            "report": str(tmp_path / "big_sandy.md"),
            "specification": str(tmp_path / "big_sandy.json"),
            "sections": 2,
            "draft_narratives": 1,
        }
        text = (tmp_path / "big_sandy.md").read_text(encoding="utf-8")
        assert text.startswith("# Big Sandy River\n")
        assert "## 1. Introduction" in text and "| 1 | 24,552 |" in text
        assert text.count(DRAFT_LABEL) == 1
        stored = json.loads((tmp_path / "big_sandy.json").read_text(encoding="utf-8"))
        assert [s["kind"] for s in stored["sections"]] == ["text", "ungaged"]
        assert stored["sections"][1]["narrative_status"] == "reviewed"
        assert stored["date"] == "2026-11-16"

    def test_render_docx(self, report, tmp_path):
        docx = pytest.importorskip("docx")
        written = report.render("big_sandy.docx")
        document = docx.Document(written["report"])
        assert document.paragraphs[0].text == "Big Sandy River"
        assert (tmp_path / "big_sandy.json").is_file()

    def test_render_makes_the_directory(self, tmp_path):
        session = ReportSession("New", output_dir=tmp_path / "reports" / "2026")
        session.add_text("One", TEXT)
        assert (tmp_path / "reports" / "2026" / "new.md").samefile(
            session.render("new.md")["report"]
        )

    @pytest.mark.parametrize(
        "filename",
        ["../big_sandy.md", "sub/big_sandy.md", "big_sandy.txt", "big_sandy", ".md", ""],
    )
    def test_render_refuses_a_bad_name(self, report, tmp_path, filename):
        with pytest.raises(ReportError, match="give a file name ending in .docx or .md"):
            report.render(filename)
        assert list(tmp_path.iterdir()) == []

    def test_render_refuses_an_empty_report(self, session, tmp_path):
        with pytest.raises(ReportError, match="no sections yet"):
            session.render("empty.md")
        assert list(tmp_path.iterdir()) == []

    def test_save_and_load(self, report, tmp_path):
        assert report.save("big_sandy.json") == {
            "specification": str(tmp_path / "big_sandy.json"),
            "sections": 2,
        }
        other = ReportSession("Something else", output_dir=tmp_path)
        other.add_text("Old", TEXT)
        outline = other.load("big_sandy.json")
        # The sections are new to this session, so they take its next ids.
        assert [s["section_id"] for s in outline["sections"]] == ["s2", "s3"]
        assert [s["narrative"] for s in outline["sections"]] == ["draft", "reviewed"]
        assert {**outline, "sections": []} == {**report.outline(), "sections": []}
        assert other.spec.to_dict() == report.spec.to_dict()
        assert other.section_markdown("s3") == report.section_markdown("s2")

    def test_load_refuses_and_keeps_the_report(self, report, tmp_path):
        before = report.spec.to_dict()
        with pytest.raises(ReportError, match="no specification missing.json"):
            report.load("missing.json")
        with pytest.raises(ReportError, match="give a file name ending in .json"):
            report.load("big_sandy.md")
        with pytest.raises(ReportError, match="give a file name ending in .json"):
            report.save("big_sandy.docx")
        bad = ReportSpec("Bad", sections=session_sections()).to_dict()
        (tmp_path / "bad.json").write_text(json.dumps(bad), encoding="utf-8")
        with pytest.raises(ReportError, match=r"section 2 \(ungaged\): .*missing option"):
            report.load("bad.json")
        (tmp_path / "not_a_spec.json").write_text("[1, 2]", encoding="utf-8")
        with pytest.raises(ReportError, match="a report spec is a JSON object"):
            report.load("not_a_spec.json")
        assert report.spec.to_dict() == before
        # The failed loads issued no ids.
        assert report.add_text("Next", TEXT)["section_id"] == "s3"


def session_sections():
    """A text section that builds, then an ungaged section that cannot."""
    from flowfreq.reporting import SectionSpec

    return [
        SectionSpec("text", title="Fine", narrative=TEXT),
        SectionSpec("ungaged", options={"state": "TN"}),
    ]
