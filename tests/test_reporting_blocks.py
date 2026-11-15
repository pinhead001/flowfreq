"""Tests for flowfreq.reporting.blocks: the block types and the Markdown subset."""

from __future__ import annotations

import pytest

from flowfreq.reporting.blocks import (
    Bullets,
    Figure,
    Heading,
    Narrative,
    Paragraph,
    Table,
    narrative_blocks,
    number_captions,
    parse_markdown,
    plain,
    split_inline,
    without_title,
)


class TestParseMarkdown:
    def test_headings_paragraphs_and_bullets(self):
        blocks = parse_markdown(
            "# Title\n\nFirst line\nsecond line.\n\n## Part ##\n\n- one\n- two\n  wrapped\n* three\n"
        )
        assert blocks == [
            Heading(1, "Title"),
            Paragraph("First line second line."),
            Heading(2, "Part"),
            Bullets(("one", "two wrapped", "three")),
        ]

    def test_table_with_alignment(self):
        (table,) = parse_markdown("| A | B | C |\n|---|---:|:---:|\n| 1 | 2 | 3 |\n| x | y | z |\n")
        assert table == Table(
            ("A", "B", "C"), (("1", "2", "3"), ("x", "y", "z")), ("left", "right", "center")
        )

    def test_escaped_pipe_stays_in_its_cell(self):
        (table,) = parse_markdown("| A | B |\n|---|---|\n| a \\| b | c |\n")
        assert table.rows == (("a | b", "c"),)

    def test_short_row_is_padded_and_long_row_folded(self):
        (table,) = parse_markdown("| A | B |\n|---|---|\n| only |\n| 1 | 2 | 3 |\n")
        assert table.rows == (("only", ""), ("1", "2 | 3"))

    def test_rule_is_dropped_and_numbered_line_starts_a_paragraph(self):
        blocks = parse_markdown("Before\n\n---\n\nIntro\n1. first\n2. second\n")
        assert blocks == [
            Paragraph("Before"),
            Paragraph("Intro"),
            Paragraph("1. first"),
            Paragraph("2. second"),
        ]

    def test_a_pipe_in_prose_is_not_a_table(self):
        assert parse_markdown("a | b\n") == [Paragraph("a | b")]

    def test_empty_text(self):
        assert parse_markdown("") == []


class TestInline:
    def test_runs(self):
        assert split_inline("a **b** and *c* then `d`.") == [
            ("a ", ""),
            ("b", "bold"),
            (" and ", ""),
            ("c", "italic"),
            (" then ", ""),
            ("d", "code"),
            (".", ""),
        ]

    def test_a_lone_asterisk_is_text(self):
        assert split_inline("2 * 3 * 4") == [("2 * 3 * 4", "")]

    def test_plain(self):
        assert plain("**Q1%** is *large*") == "Q1% is large"


class TestNarrativeBlocks:
    def test_paragraphs_and_bullets_pass(self):
        assert narrative_blocks("One.\n\n- a\n- b\n") == [Paragraph("One."), Bullets(("a", "b"))]

    def test_a_table_is_refused(self):
        with pytest.raises(ValueError, match="cannot contain a table"):
            narrative_blocks("| Q | cfs |\n|---|---|\n| 1% | 9,999 |\n")

    def test_a_heading_is_refused(self):
        with pytest.raises(ValueError, match="cannot contain a heading"):
            narrative_blocks("## Findings\n\nText.")


class TestNarrative:
    def test_default_is_draft(self):
        assert Narrative("text").status == "draft"

    def test_unknown_status_is_refused(self):
        with pytest.raises(ValueError, match="draft, reviewed"):
            Narrative("text", "approved")


def test_without_title_splits_a_leading_section_heading():
    blocks = [Heading(1, "Site"), Paragraph("x")]
    assert without_title(blocks) == ("Site", [Paragraph("x")])
    assert without_title([Paragraph("x")]) == (None, [Paragraph("x")])
    assert without_title([Heading(2, "Part")]) == (None, [Heading(2, "Part")])


def test_number_captions_counts_tables_and_figures_apart():
    blocks = [
        Table(("A",), (("1",),), caption="First."),
        Table(("A",), (("1",),)),  # uncaptioned: not numbered
        Figure(b"png", "Curve."),
        Table(("A",), (("1",),), caption="Second."),
        Figure(b"png", "Peaks."),
    ]
    numbered = number_captions(blocks)
    assert [getattr(b, "caption") for b in numbered] == [
        "Table 1. First.",
        None,
        "Figure 1. Curve.",
        "Table 2. Second.",
        "Figure 2. Peaks.",
    ]
    assert blocks[0].caption == "First."  # the input is not changed
