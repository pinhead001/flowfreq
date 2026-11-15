"""Tests for flowfreq.reporting.word: the DOCX a report writes, read back with python-docx."""

from __future__ import annotations

import builtins
from typing import Any, Dict

import pytest

from flowfreq.reporting import DRAFT_LABEL, Report, ReportSpec, SectionSpec, Table, build_report
from flowfreq.reporting.word import _column_widths, _text_in
from tests.fixtures import big_sandy as bs

YEARS = sorted(bs.SYSTEMATIC_PEAKS)
PEAKS: Dict[str, Any] = {
    "water_year": YEARS,
    "peak_flow_cfs": [bs.SYSTEMATIC_PEAKS[y] for y in YEARS],
}


@pytest.fixture(scope="module")
def report() -> Report:
    return build_report(
        ReportSpec(
            title="Big Sandy River flood frequency",
            author="A. Analyst",
            date="2026-11-14",
            sections=[
                SectionSpec(
                    "text",
                    title="Purpose",
                    narrative="This report gives **flood quantiles** for *design*.\n\n- one\n- two",
                    narrative_status="reviewed",
                ),
                SectionSpec(
                    "gaged_ffa",
                    options={
                        "site_no": "03606500",
                        "peaks": PEAKS,
                        "regional_skew": bs.REGIONAL_SKEW,
                        "regional_skew_se": bs.REGIONAL_SKEW_SD,
                        "stationarity": False,
                    },
                    narrative="The 1-percent flow is in the table above.",
                ),
            ],
        )
    )


@pytest.fixture(scope="module")
def document(report: Report, tmp_path_factory: pytest.TempPathFactory) -> Any:
    docx = pytest.importorskip("docx")
    path = report.save(tmp_path_factory.mktemp("docx") / "report.docx")
    assert path.stat().st_size > 10_000
    return docx.Document(str(path))


class TestDocx:
    def test_title_and_headings(self, document):
        styled = [(p.style.name, p.text) for p in document.paragraphs if p.text]
        assert styled[0] == ("Title", "Big Sandy River flood frequency")
        assert ("Heading 1", "1. Purpose") in styled
        assert ("Heading 1", "2. Flood frequency at USGS 03606500") in styled
        assert ("Heading 2", "Annual peaks") in styled
        assert ("Heading 1", "About this report") in styled
        assert document.core_properties.title == "Big Sandy River flood frequency"
        assert document.core_properties.author == "A. Analyst"

    def test_tables_carry_the_computed_cells(self, document, report):
        tables = [b for b in report.blocks if isinstance(b, Table)]
        assert len(document.tables) == len(tables) == 4
        for written, block in zip(document.tables, tables):
            cells = [tuple(c.text for c in row.cells) for row in written.rows]
            assert cells == [block.header, *block.rows]
        quantiles = document.tables[2]
        assert [c.text for c in quantiles.rows[7].cells][:3] == ["1", "100", "18,069"]

    def test_table_header_repeats_and_is_bold(self, document):
        header = document.tables[2].rows[0]
        assert header._tr.xpath("./w:trPr/w:tblHeader")
        assert all(run.bold for cell in header.cells for run in cell.paragraphs[0].runs)
        body = document.tables[2].rows[1]
        assert not body._tr.xpath("./w:trPr/w:tblHeader")

    def test_captions(self, document):
        captions = [p.text for p in document.paragraphs if p.style.name == "Caption"]
        assert captions[0] == "Table 1. Station and peak-flow record."
        assert "Figure 1. Fitted frequency curve with 90-percent confidence limits" in captions[3]
        assert captions[-1] == "Table 4. Annual peak flows used in the analysis."

    def test_figures_are_embedded(self, document, report):
        assert len(document.inline_shapes) == len(report.figures()) == 2
        assert all(abs(shape.width.inches - 6.5) < 0.01 for shape in document.inline_shapes)

    def test_draft_narrative_is_labeled_and_indented(self, document):
        paragraphs = list(document.paragraphs)
        (index,) = [i for i, p in enumerate(paragraphs) if p.text == DRAFT_LABEL]
        label, draft = paragraphs[index], paragraphs[index + 1]
        assert label.runs[0].bold and label.runs[0].italic
        assert draft.text == "The 1-percent flow is in the table above."
        assert draft.paragraph_format.left_indent == label.paragraph_format.left_indent
        assert label.paragraph_format.left_indent.inches == pytest.approx(0.3, abs=0.01)

    def test_reviewed_narrative_is_ordinary_text_with_inline_styles(self, document):
        (purpose,) = [p for p in document.paragraphs if p.text.startswith("This report gives")]
        assert purpose.paragraph_format.left_indent is None
        assert [(r.text, r.bold, r.italic) for r in purpose.runs] == [
            ("This report gives ", None, None),
            ("flood quantiles", True, None),
            (" for ", None, None),
            ("design", None, True),
            (".", None, None),
        ]
        bullets = [p.text for p in document.paragraphs if p.style.name == "List Bullet"]
        assert bullets == ["one", "two"]

    def test_page_number_field_in_the_footer(self, document):
        footer = document.sections[0].footer.paragraphs[0]
        from docx.oxml.ns import qn

        (field,) = footer._p.xpath("./w:fldSimple")
        assert field.get(qn("w:instr")) == "PAGE"

    def test_quantile_columns_hold_their_numbers(self, document):
        for row in document.tables[2].rows:
            widths = [cell.width.inches for cell in row.cells]
            assert sum(widths) <= 6.5 + 1e-3
            assert min(widths) >= 0.17 + _text_in("18,069") - 1e-3

    def test_rows_do_not_split_and_a_short_table_stays_on_one_page(self, document):
        quantiles, peaks = document.tables[2], document.tables[3]
        assert len(quantiles.rows) == 10 and len(peaks.rows) == 45
        for table in (quantiles, peaks):
            assert all(row._tr.xpath("./w:trPr/w:cantSplit") for row in table.rows)

        def kept(table):
            return [
                all(c.paragraphs[0].paragraph_format.keep_with_next for c in row.cells)
                for row in table.rows
            ]

        assert kept(quantiles) == [True] * 9 + [False]  # the last row ends the block
        assert not any(kept(peaks))  # a long table breaks, repeating its header

    def test_margins_leave_the_text_width(self, document):
        section = document.sections[0]
        width = section.page_width.inches - section.left_margin.inches - section.right_margin.inches
        assert width == pytest.approx(6.5, abs=0.01)


class TestColumnWidths:
    def test_text_width_counts_capitals_digits_and_bold(self):
        assert _text_in("") == 0
        assert _text_in("LC92STOR") > _text_in("12345678") > _text_in("abcdefgh")
        assert _text_in("mw") > _text_in("ae") > _text_in("il")
        assert _text_in("Flow", bold=True) == pytest.approx(1.1 * _text_in("Flow"))

    def test_a_table_that_fits_takes_its_cells_widths(self):
        table = Table(("Water year", "Item"), (("1930", "a fairly long description"),))
        year, item = _column_widths(table)
        assert year == pytest.approx(0.17 + _text_in("Water year", bold=True))  # the header
        assert item == pytest.approx(0.17 + _text_in("a fairly long description"))  # the cell
        assert year + item < 6.5  # not stretched across the page

    def test_a_wide_table_fills_the_text_width_and_keeps_numbers_whole(self):
        table = Table(
            ("Annual exceedance probability (percent)", "Flow (cfs)", "Notes", "Source"),
            (("0.2", "1,234,567", "some remarks " * 20, "a long citation " * 20),),
        )
        widths = _column_widths(table)
        assert sum(widths) == pytest.approx(6.5)
        # The flow column holds its longest word, the number, on one line.
        assert widths[1] >= 0.17 + _text_in("1,234,567")
        assert min(widths[2:]) > widths[0] > widths[1]  # running text gets the most room
        assert max(widths) < 3.0  # and no column more than the cap

    def test_one_long_column_takes_what_the_others_leave(self):
        table = Table(("Item", "Value"), (("Station", "x " * 250),))
        item, value = _column_widths(table)
        assert item == pytest.approx(0.17 + _text_in("Station"))
        assert item + value == pytest.approx(6.5)

    def test_too_many_columns_share_the_width(self):
        table = Table(("Extrapolated",) * 12, (("no",) * 12,))
        widths = _column_widths(table)
        assert sum(widths) == pytest.approx(6.5) and len(set(round(w, 9) for w in widths)) == 1

    def test_markup_does_not_count_as_width(self):
        bold = _column_widths(Table(("A",), (("**12,345**",),)))
        assert bold == _column_widths(Table(("A",), (("12,345",),)))


def test_missing_python_docx_says_which_extra(report, tmp_path, monkeypatch):
    real_import = builtins.__import__

    def no_docx(name: str, *args: Any, **kwargs: Any) -> Any:
        if name == "docx" or name.startswith("docx."):
            raise ImportError("No module named 'docx'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_docx)
    with pytest.raises(ImportError, match=r"pip install 'flowfreq\[report\]'"):
        report.save(tmp_path / "report.docx")
    assert not (tmp_path / "report.docx").exists()
