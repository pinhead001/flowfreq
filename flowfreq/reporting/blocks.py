"""Report content blocks, and the Markdown subset flowfreq's own outputs are written in.

A report is a flat list of blocks. Every renderer (Markdown, DOCX) reads the same list, so a
table or figure is produced once, by flowfreq, and laid out per format. :func:`parse_markdown`
turns the ``to_markdown`` text of an existing result (a site estimate, a stationarity page)
into blocks, which is how those results enter a report without a second implementation.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from typing import List, Optional, Sequence, Tuple, Union

#: Narrative states. A draft renders under a visible label; a reviewed one as ordinary text.
NARRATIVE_STATUSES: Tuple[str, ...] = ("draft", "reviewed")

#: The label a draft narrative carries in every output format.
DRAFT_LABEL = "Draft narrative (written by a language model; not reviewed)"


@dataclass(frozen=True)
class Heading:
    """A heading. Level 0 is the report title, 1 a section, 2 and deeper its parts."""

    level: int
    text: str


@dataclass(frozen=True)
class Paragraph:
    """A paragraph; ``text`` may carry ``**bold**``, ``*italic*`` and ```code```."""

    text: str


@dataclass(frozen=True)
class Bullets:
    """A bullet list; each item is inline text as :class:`Paragraph` takes it."""

    items: Tuple[str, ...]


@dataclass(frozen=True)
class Table:
    """A table of already-formatted cells.

    Attributes
    ----------
    header : tuple of str
    rows : tuple of tuple of str
        Each row as long as ``header``.
    align : tuple of str
        ``"left"``, ``"right"`` or ``"center"`` per column; empty means all left.
    caption : str, optional
        Without its number; :func:`number_captions` adds ``Table N.``.
    """

    header: Tuple[str, ...]
    rows: Tuple[Tuple[str, ...], ...]
    align: Tuple[str, ...] = ()
    caption: Optional[str] = None


@dataclass(frozen=True)
class Figure:
    """A figure as PNG bytes, with its caption (unnumbered) and printed width in inches."""

    png: bytes
    caption: str
    width_in: float = 6.5


@dataclass(frozen=True)
class Narrative:
    """Prose that flowfreq did not compute: an author's, or a language model's draft.

    ``text`` holds paragraphs and bullet lists only (:func:`narrative_blocks` enforces it), so
    a table in a report always comes from flowfreq.
    """

    text: str
    status: str = "draft"

    def __post_init__(self) -> None:
        if self.status not in NARRATIVE_STATUSES:
            raise ValueError(
                f"narrative status must be one of {', '.join(NARRATIVE_STATUSES)}, "
                f"got {self.status!r}"
            )


Block = Union[Heading, Paragraph, Bullets, Table, Figure, Narrative]

_HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")
_BULLET = re.compile(r"^\s*[-*+]\s+(.*)$")
_NUMBERED = re.compile(r"^\s*\d+[.)]\s+")
_RULE = re.compile(r"^\s*([-*_])(\s*\1){2,}\s*$")
_SEPARATOR_CELL = re.compile(r"^:?-+:?$")
_UNESCAPED_PIPE = re.compile(r"(?<!\\)\|")
#: Inline markup, longest first: bold, code, then italic that hugs its text.
_INLINE = re.compile(r"(\*\*.+?\*\*|`[^`]+`|(?<![\w*])\*[^\s*](?:[^*]*[^\s*])?\*(?![\w*]))")


def _cells(line: str) -> List[str]:
    body = line.strip()
    if body.startswith("|"):
        body = body[1:]
    if body.endswith("|") and not body.endswith("\\|"):
        body = body[:-1]
    return [c.strip().replace("\\|", "|") for c in _UNESCAPED_PIPE.split(body)]


def _is_separator(line: str) -> bool:
    cells = _cells(line)
    return bool(cells) and "|" in line and all(_SEPARATOR_CELL.match(c) for c in cells)


def _alignment(cell: str) -> str:
    if cell.startswith(":") and cell.endswith(":"):
        return "center"
    return "right" if cell.endswith(":") else "left"


def _starts_block(line: str, following: Optional[str]) -> bool:
    """Whether ``line`` opens something other than a paragraph continuation."""
    if _HEADING.match(line) or _BULLET.match(line) or _NUMBERED.match(line) or _RULE.match(line):
        return True
    return "|" in line and following is not None and _is_separator(following)


def parse_markdown(text: str) -> List[Block]:
    """Blocks for the Markdown flowfreq writes: headings, paragraphs, bullets, pipe tables.

    Parameters
    ----------
    text : str
        Markdown. ``#`` is read as heading level 1 (a section), ``##`` as 2, and so on.
        Horizontal rules are dropped; a numbered line starts its own paragraph.

    Returns
    -------
    list of Block
        :class:`Heading`, :class:`Paragraph`, :class:`Bullets` and :class:`Table` only.
    """
    lines = text.splitlines()
    blocks: List[Block] = []
    i, n = 0, len(lines)
    while i < n:
        line = lines[i].rstrip()
        following = lines[i + 1] if i + 1 < n else None
        if not line.strip() or _RULE.match(line):
            i += 1
            continue
        heading = _HEADING.match(line)
        if heading:
            blocks.append(Heading(len(heading.group(1)), heading.group(2)))
            i += 1
            continue
        if "|" in line and following is not None and _is_separator(following):
            header = _cells(line)
            align = tuple(_alignment(c) for c in _cells(following))
            i += 2
            rows: List[Tuple[str, ...]] = []
            while i < n and lines[i].strip().startswith("|"):
                cells = _cells(lines[i])
                # A short row is padded and a long one folded into its last cell, so a
                # stray pipe in a cell cannot shift the columns of the rows after it.
                if len(cells) > len(header):
                    cells = cells[: len(header) - 1] + [" | ".join(cells[len(header) - 1 :])]
                rows.append(tuple(cells + [""] * (len(header) - len(cells))))
                i += 1
            blocks.append(Table(tuple(header), tuple(rows), align[: len(header)]))
            continue
        bullet = _BULLET.match(line)
        if bullet:
            items: List[str] = []
            while i < n:
                current = lines[i].rstrip()
                item = _BULLET.match(current)
                if item:
                    items.append(item.group(1).strip())
                elif current.startswith((" ", "\t")) and current.strip() and items:
                    items[-1] += " " + current.strip()  # a wrapped item
                else:
                    break
                i += 1
            blocks.append(Bullets(tuple(items)))
            continue
        parts = [line.strip()]
        i += 1
        while i < n and lines[i].strip():
            nxt = lines[i + 1] if i + 1 < n else None
            if _starts_block(lines[i], nxt):
                break
            parts.append(lines[i].strip())
            i += 1
        blocks.append(Paragraph(" ".join(parts)))
    return blocks


def split_inline(text: str) -> List[Tuple[str, str]]:
    """Inline text as ``(text, style)`` runs; style is ``""``, ``"bold"``, ``"italic"``, ``"code"``."""
    runs: List[Tuple[str, str]] = []
    for part in _INLINE.split(text):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**") and len(part) > 4:
            runs.append((part[2:-2], "bold"))
        elif part.startswith("`") and part.endswith("`") and len(part) > 2:
            runs.append((part[1:-1], "code"))
        elif _INLINE.fullmatch(part):
            runs.append((part[1:-1], "italic"))
        else:
            runs.append((part, ""))
    return runs


def plain(text: str) -> str:
    """Inline text with its markup removed."""
    return "".join(run for run, _ in split_inline(text))


def narrative_blocks(text: str) -> List[Block]:
    """A narrative's paragraphs and bullet lists.

    Raises
    ------
    ValueError
        The text holds a heading or a table. Headings would change the report's outline,
        and a table in a report must come from flowfreq, not from prose.
    """
    blocks = parse_markdown(text)
    for block in blocks:
        if isinstance(block, Table):
            raise ValueError(
                "a narrative cannot contain a table: tables in a report are computed by "
                "flowfreq. Describe the result in words, or add the section that computes it."
            )
        if isinstance(block, Heading):
            raise ValueError(
                "a narrative cannot contain a heading; put the text under its own section"
            )
    return blocks


def without_title(blocks: Sequence[Block]) -> Tuple[Optional[str], List[Block]]:
    """Split a parsed document into its leading level-1 heading (if any) and the rest."""
    if blocks and isinstance(blocks[0], Heading) and blocks[0].level == 1:
        return blocks[0].text, list(blocks[1:])
    return None, list(blocks)


def number_captions(blocks: Sequence[Block]) -> List[Block]:
    """Number the captioned tables and the figures in order: ``Table 1.``, ``Figure 1.``."""
    out: List[Block] = []
    tables = figures = 0
    for block in blocks:
        if isinstance(block, Table) and block.caption:
            tables += 1
            out.append(replace(block, caption=f"Table {tables}. {block.caption}"))
        elif isinstance(block, Figure):
            figures += 1
            out.append(replace(block, caption=f"Figure {figures}. {block.caption}"))
        else:
            out.append(block)
    return out
