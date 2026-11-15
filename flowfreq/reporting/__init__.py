"""Reports assembled from flowfreq results: a JSON specification in, Markdown or DOCX out.

A report is described by a :class:`ReportSpec`: a title and a list of sections, each a kind
(``gaged_ffa``, ``ungaged``, ``weighted``, ``text``) with its inputs. :func:`build_report` computes every
table and figure from those inputs, so no number in a report is typed by hand or written by
a language model. Prose that flowfreq did not compute rides along as a section's
``narrative``, and renders under a draft label until someone marks it reviewed.

>>> from flowfreq.reporting import ReportSpec, SectionSpec, build_report
>>> spec = ReportSpec(
...     title="Big Sandy River flood frequency",
...     sections=[SectionSpec("gaged_ffa", options={"site_no": "03606500",
...                                                 "regional_skew_from_site": True})],
... )
>>> report = build_report(spec)          # downloads the peaks, fits, draws
>>> report.save("big_sandy.docx")        # needs: pip install "flowfreq[report]"
>>> report.spec.save("big_sandy.json")   # the peaks are in it: rebuilds offline
"""

from .blocks import (
    DRAFT_LABEL,
    NARRATIVE_STATUSES,
    Block,
    Bullets,
    Figure,
    Heading,
    Narrative,
    Paragraph,
    Table,
    parse_markdown,
)
from .build import Report, build_report, render_markdown
from .sections import SECTION_KINDS, SectionContent, SectionKind, build_section
from .spec import SCHEMA, ReportError, ReportSpec, SectionSpec

__all__ = [
    "DRAFT_LABEL",
    "NARRATIVE_STATUSES",
    "SCHEMA",
    "SECTION_KINDS",
    "Block",
    "Bullets",
    "Figure",
    "Heading",
    "Narrative",
    "Paragraph",
    "Report",
    "ReportError",
    "ReportSpec",
    "SectionContent",
    "SectionKind",
    "SectionSpec",
    "Table",
    "build_report",
    "build_section",
    "parse_markdown",
    "render_markdown",
]
