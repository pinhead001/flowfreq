"""Build a report from its specification, and write it as Markdown or DOCX."""

from __future__ import annotations

import datetime as _dt
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, List, Tuple, Union

from .blocks import (
    DRAFT_LABEL,
    Block,
    Bullets,
    Figure,
    Heading,
    Narrative,
    Paragraph,
    Table,
    narrative_blocks,
    number_captions,
)
from .sections import _version, build_section, section_title
from .spec import ReportError, ReportSpec, SectionSpec

_ALIGN_RULE = {"left": "---", "right": "---:", "center": ":---:"}


def _cell(text: str) -> str:
    return text.replace("|", "\\|").replace("\n", " ")


def _figure_names(blocks: List[Block]) -> List[Tuple[str, Figure]]:
    figures = [b for b in blocks if isinstance(b, Figure)]
    return [(f"figure_{i}.png", figure) for i, figure in enumerate(figures, 1)]


def render_markdown(blocks: List[Block], image_dir: str = "figures") -> str:
    """Blocks as Markdown.

    Parameters
    ----------
    blocks : list of Block
    image_dir : str
        The directory figure links point into, relative to the Markdown file. Figures are
        named ``figure_1.png``, ``figure_2.png`` in order.

    Returns
    -------
    str
        A table's caption sits above it and a figure's below. A draft narrative is a block
        quote under the draft label.
    """
    out: List[str] = []
    figure_no = 0
    for block in blocks:
        if isinstance(block, Heading):
            out.append("#" * (block.level + 1) + " " + block.text)
        elif isinstance(block, Paragraph):
            out.append(block.text)
        elif isinstance(block, Bullets):
            out.append("\n".join(f"- {item}" for item in block.items))
        elif isinstance(block, Table):
            align = block.align or ("left",) * len(block.header)
            lines = [block.caption, ""] if block.caption else []
            lines.append("| " + " | ".join(_cell(c) for c in block.header) + " |")
            lines.append("|" + "|".join(_ALIGN_RULE[a] for a in align) + "|")
            lines.extend("| " + " | ".join(_cell(c) for c in row) + " |" for row in block.rows)
            out.append("\n".join(lines))
        elif isinstance(block, Figure):
            figure_no += 1  # the order _figure_names counts in
            name = f"figure_{figure_no}.png"
            link = f"{image_dir}/{name}" if image_dir else name
            out.append(f"![{block.caption}]({link})\n\n{block.caption}")
        elif block.status == "draft":
            quoted = [f"> **{DRAFT_LABEL}**", ">"]
            for line in block.text.strip().splitlines():
                quoted.append(f"> {line}" if line.strip() else ">")
            out.append("\n".join(quoted))
        else:
            out.append(block.text.strip())
    return "\n\n".join(out) + "\n"


@dataclass(frozen=True)
class Report:
    """A built report.

    Attributes
    ----------
    spec : ReportSpec
        The specification as built: the date fixed and everything a section fetched written
        into its options. Save it (``report.spec.save(path)``) and :func:`build_report` gives
        the same report again without the network.
    blocks : list of Block
        The whole report in order, captions numbered.
    """

    spec: ReportSpec
    blocks: List[Block]

    def figures(self) -> Dict[str, bytes]:
        """The figures as ``{file name: PNG bytes}``, named as the Markdown links them."""
        return {name: figure.png for name, figure in _figure_names(self.blocks)}

    def to_markdown(self, image_dir: str = "figures") -> str:
        """The report as Markdown; see :func:`render_markdown`."""
        return render_markdown(self.blocks, image_dir)

    def save(self, path: Union[str, Path]) -> Path:
        """Write the report, in the format the suffix names.

        Parameters
        ----------
        path : str or Path
            ``.docx`` (needs the ``report`` extra) or ``.md``. Markdown writes its figures
            beside it, in ``<name>_files/``.

        Returns
        -------
        Path
            The file written.

        Raises
        ------
        ReportError
            Any other suffix.
        ImportError
            ``.docx`` without python-docx installed.
        """
        out = Path(path)
        suffix = out.suffix.lower()
        if suffix == ".docx":
            from .word import render_docx

            render_docx(self, out)
        elif suffix == ".md":
            folder = out.with_name(out.stem + "_files")
            figures = self.figures()
            if figures:
                folder.mkdir(parents=True, exist_ok=True)
                for name, png in figures.items():
                    (folder / name).write_bytes(png)
            out.write_text(self.to_markdown(folder.name), encoding="utf-8")
        else:
            raise ReportError(
                f"cannot write a report as {out.suffix or 'a file with no suffix'}: "
                "use .docx or .md (Word exports a .docx to PDF)"
            )
        return out


def build_report(spec: ReportSpec) -> Report:
    """Compute every section of a report.

    Parameters
    ----------
    spec : ReportSpec

    Returns
    -------
    Report
        Its ``spec`` is the specification as built.

    Raises
    ------
    ReportError
        A section cannot be built as specified; the message names the section. A text
        section without a narrative, or a narrative holding a table or a heading.
    ValueError, LookupError
        Refusals the analyses raise themselves pass through: a regulated record
        (:class:`~flowfreq.regulation.RegulatedRecordError`), a characteristic outside an
        equation's range, a state with no equations.
    """
    date = spec.date or _dt.date.today().isoformat()
    blocks: List[Block] = [Heading(0, spec.title)]
    if spec.subtitle:
        blocks.append(Paragraph(f"*{spec.subtitle}*"))
    byline = [
        f"{label}: {value}"
        for label, value in (("Project", spec.project), ("Prepared by", spec.author))
        if value
    ]
    blocks.append(Paragraph(". ".join([*byline, f"Date: {date}"]) + "."))

    references: List[str] = []
    sections: List[SectionSpec] = []
    drafts = 0
    for number, section in enumerate(spec.sections, 1):
        where = f"section {number} ({section.kind})"
        try:
            content = build_section(section.kind, section.options)
        except ReportError as exc:
            raise ReportError(f"{where}: {exc}") from exc
        title = section_title(section.kind, section.title, content.title)
        if not title:
            raise ReportError(f"{where}: needs a title")
        blocks.append(Heading(1, f"{number}. {title}"))
        blocks.extend(content.blocks)
        if section.narrative and section.narrative.strip():
            try:
                narrative_blocks(section.narrative)
            except ValueError as exc:
                raise ReportError(f"{where}: {exc}") from exc
            blocks.append(Narrative(section.narrative, section.narrative_status))
            drafts += section.narrative_status == "draft"
        elif section.kind == "text":
            raise ReportError(f"{where}: a text section needs a narrative")
        references.extend(r for r in content.references if r not in references)
        sections.append(replace(section, options=content.options))

    if references:
        blocks.append(Heading(1, "References"))
        blocks.extend(Paragraph(r) for r in sorted(references, key=str.lower))
    about = (
        f"The tables and figures in this report were computed by flowfreq {_version()} from "
        "the inputs recorded in the report specification; building the saved specification "
        "again reproduces them."
    )
    if drafts:
        about += (
            " Text marked as draft narrative was written by a language model and has not "
            "been reviewed: it is not a finding of the analysis until an analyst has checked "
            "it against the tables."
        )
    blocks.append(Heading(1, "About this report"))
    blocks.append(Paragraph(about))
    return Report(spec=replace(spec, date=date, sections=sections), blocks=number_captions(blocks))
