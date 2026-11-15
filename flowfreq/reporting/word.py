"""DOCX output, through python-docx (the ``report`` extra)."""

from __future__ import annotations

import io
from pathlib import Path
from typing import TYPE_CHECKING, Any, List, Tuple, Union

from .blocks import (
    DRAFT_LABEL,
    Block,
    Bullets,
    Figure,
    Heading,
    Paragraph,
    Table,
    narrative_blocks,
    plain,
    split_inline,
)

if TYPE_CHECKING:
    from .build import Report

_TABLE_POINTS = 9
_CODE_FONT = "Consolas"
_DRAFT_INDENT_IN = 0.3
_TEXT_WIDTH_IN = 6.5
#: Rough widths, in inches, of the characters of table text by class, of a cell's two side
#: margins, and how much wider bold is. They only have to be good enough to size columns:
#: too narrow breaks a number across lines, too wide takes the room running text needs.
_NARROW, _NARROW_IN = frozenset("iljtfrI .,;:'!()[]|-/"), 0.04
_WIDE, _WIDE_IN = frozenset("mwMW%@"), 0.106
_UPPER_IN, _DIGIT_IN, _OTHER_IN = 0.081, 0.069, 0.065
_CELL_PAD_IN = 0.17
_BOLD = 1.1
#: A table that has to be squeezed does not give one column more than this before its
#: cells wrap.
_MAX_COLUMN_IN = 3.0
#: A table of no more rows than this is kept on one page.
_KEEP_TOGETHER_ROWS = 12


def _text_in(text: str, bold: bool = False) -> float:
    """Rough width of a run of table text, in inches."""

    def one(ch: str) -> float:
        if ch in _NARROW:
            return _NARROW_IN
        if ch in _WIDE:
            return _WIDE_IN
        return _UPPER_IN if ch.isupper() else _DIGIT_IN if ch.isdigit() else _OTHER_IN

    return sum(map(one, text)) * (_BOLD if bold else 1.0)


def _widths_in(text: str, bold: bool = False) -> Tuple[float, float]:
    """Widths of the longest word of ``text`` and of the whole of it."""
    return max((_text_in(w, bold) for w in text.split()), default=0.0), _text_in(text, bold)


def _column_widths(table: Table) -> List[float]:
    """Column widths in inches, together no wider than the text.

    A column is as wide as its longest cell when the table fits. When it does not, every
    column keeps room for its longest word, so a number is never broken across two lines,
    and the columns of running text give up the rest.
    """
    least: List[float] = []
    full: List[float] = []
    for i, head in enumerate(table.header):
        sizes = [_widths_in(plain(head), bold=True)]
        sizes += [_widths_in(plain(row[i])) for row in table.rows]
        least.append(_CELL_PAD_IN + max(word for word, _ in sizes))
        full.append(_CELL_PAD_IN + max(whole for _, whole in sizes))
    if sum(full) <= _TEXT_WIDTH_IN:
        return full
    if sum(least) >= _TEXT_WIDTH_IN:
        return [_TEXT_WIDTH_IN * w / sum(least) for w in least]
    wanted = [max(low, min(high, _MAX_COLUMN_IN)) for low, high in zip(least, full)]
    if sum(wanted) <= _TEXT_WIDTH_IN:
        # Only the capped columns were too wide: they take what the others leave.
        spare = _TEXT_WIDTH_IN - sum(wanted)
        capped = [high > want for high, want in zip(full, wanted)]
        return [w + spare / sum(capped) if c else w for w, c in zip(wanted, capped)]
    share = (_TEXT_WIDTH_IN - sum(least)) / (sum(wanted) - sum(least))
    return [low + share * (high - low) for low, high in zip(least, wanted)]


def _runs(paragraph: Any, text: str, *, size: Any = None, bold: bool = False) -> None:
    """Add ``text`` to a python-docx paragraph, one run per inline style."""
    for content, style in split_inline(text):
        run = paragraph.add_run(content)
        run.bold = bold or style == "bold" or None
        run.italic = style == "italic" or None
        if style == "code":
            run.font.name = _CODE_FONT
        if size is not None:
            run.font.size = size


def render_docx(report: "Report", path: Union[str, Path]) -> Path:
    """Write a report as a Word document.

    Parameters
    ----------
    report : Report
    path : str or Path

    Returns
    -------
    Path

    Raises
    ------
    ImportError
        python-docx is not installed.
    """
    try:
        import docx
    except ImportError as exc:
        raise ImportError(
            "writing a .docx report needs python-docx: pip install 'flowfreq[report]'"
        ) from exc
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Inches, Pt

    alignment = {
        "left": WD_ALIGN_PARAGRAPH.LEFT,
        "right": WD_ALIGN_PARAGRAPH.RIGHT,
        "center": WD_ALIGN_PARAGRAPH.CENTER,
    }
    document = docx.Document()
    document.core_properties.title = report.spec.title
    if report.spec.author:
        document.core_properties.author = report.spec.author
    for section in document.sections:
        section.left_margin = section.right_margin = Inches((8.5 - _TEXT_WIDTH_IN) / 2)
        section.top_margin = section.bottom_margin = Inches(1)

    # Page number, centered in the footer: a PAGE field, which python-docx has no call for.
    footer = document.sections[0].footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    field = OxmlElement("w:fldSimple")
    field.set(qn("w:instr"), "PAGE")
    footer._p.append(field)

    def add(blocks: List[Block], indent: Any = None) -> None:
        for block in blocks:
            if isinstance(block, Heading):
                document.add_heading(block.text, level=min(block.level, 9))
            elif isinstance(block, Paragraph):
                paragraph = document.add_paragraph()
                paragraph.paragraph_format.left_indent = indent
                _runs(paragraph, block.text)
            elif isinstance(block, Bullets):
                for item in block.items:
                    paragraph = document.add_paragraph(style="List Bullet")
                    if indent is not None:
                        paragraph.paragraph_format.left_indent = indent + Inches(0.25)
                    _runs(paragraph, item)
            elif isinstance(block, Table):
                if block.caption:
                    caption = document.add_paragraph(block.caption, style="Caption")
                    caption.paragraph_format.keep_with_next = True
                table = document.add_table(rows=1, cols=len(block.header))
                table.style = "Table Grid"
                align = block.align or ("left",) * len(block.header)
                header = table.rows[0]
                # Repeat the header row when the table breaks across pages.
                repeat = OxmlElement("w:tblHeader")
                repeat.set(qn("w:val"), "true")
                header._tr.get_or_add_trPr().append(repeat)
                # Word reads a column's width from its cells, so every cell carries it.
                widths = [Inches(w) for w in _column_widths(block)]
                table.autofit = False
                for column, width in zip(table.columns, widths):
                    column.width = width
                for cell, text, side, width in zip(header.cells, block.header, align, widths):
                    cell.width = width
                    cell.paragraphs[0].alignment = alignment[side]
                    _runs(cell.paragraphs[0], text, size=Pt(_TABLE_POINTS), bold=True)
                for row in block.rows:
                    for cell, text, side, width in zip(table.add_row().cells, row, align, widths):
                        cell.width = width
                        cell.paragraphs[0].alignment = alignment[side]
                        _runs(cell.paragraphs[0], text, size=Pt(_TABLE_POINTS))
                # No row breaks across a page, and a short table stays whole: Word keeps a
                # row with the next one when its paragraphs are kept with the next.
                written = list(table.rows)
                for row_ in written:
                    row_._tr.get_or_add_trPr().append(OxmlElement("w:cantSplit"))
                if len(block.rows) <= _KEEP_TOGETHER_ROWS:
                    for row_ in written[:-1]:
                        for cell in row_.cells:
                            cell.paragraphs[0].paragraph_format.keep_with_next = True
                document.add_paragraph()  # space below the table
            elif isinstance(block, Figure):
                document.add_picture(io.BytesIO(block.png), width=Inches(block.width_in))
                picture = document.paragraphs[-1]
                picture.alignment = WD_ALIGN_PARAGRAPH.CENTER
                picture.paragraph_format.keep_with_next = True
                caption = document.add_paragraph(block.caption, style="Caption")
                caption.alignment = WD_ALIGN_PARAGRAPH.CENTER
            elif block.status == "draft":
                label = document.add_paragraph()
                label.paragraph_format.left_indent = Inches(_DRAFT_INDENT_IN)
                label.paragraph_format.keep_with_next = True
                run = label.add_run(DRAFT_LABEL)
                run.bold = run.italic = True
                add(narrative_blocks(block.text), Inches(_DRAFT_INDENT_IN))
            else:
                add(narrative_blocks(block.text))

    add(report.blocks)
    out = Path(path)
    document.save(str(out))
    return out
