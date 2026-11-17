"""A report put together one step at a time: the state a chat holds between tool calls.

A :class:`ReportSession` wraps a :class:`~flowfreq.reporting.ReportSpec`. Each section is
built when it is added, so an input the analysis refuses is refused then and not when the
report is rendered, and the section is kept with its options as resolved: the downloaded
peaks, the looked-up skew, the values StreamStats returned. What a method returns is small
and JSON-serializable (the section's id, its title, its main table), which is what a chat
needs to go on with.

The division of labor is the reporting package's: every table and figure is computed by
flowfreq from the stored inputs, and the caller contributes prose only, through
:meth:`ReportSession.set_narrative`. That prose prints under a draft label until
:meth:`ReportSession.mark_reviewed` says a person has read it.

>>> session = ReportSession("Big Sandy River flood frequency")
>>> added = session.add_section("gaged_ffa", {"site_no": "03606500",
...                                           "regional_skew_from_site": True})
>>> session.set_narrative(added["section_id"], "The 1-percent flow is ...")
>>> session.render("big_sandy.docx")     # and big_sandy.json, the specification as built
"""

from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Union

from .blocks import Block, Figure, Heading, Narrative, Paragraph, Table, narrative_blocks
from .build import build_report, render_markdown
from .daily import FLOW_DURATION_CAPTION, LOW_FLOW_CAPTION
from .sections import FUTURE_FLOW_CAPTION, SectionContent, build_section
from .spec import ReportError, ReportSpec, SectionSpec

#: The table of each kind that carries its result, by caption.
MAIN_TABLES: Dict[str, str] = {
    "gaged_ffa": "Flood quantiles with 90-percent confidence limits.",
    "ungaged": "Regression flood quantiles with prediction limits.",
    "weighted": "At-site and regression flood quantiles and their weighted estimate.",
    "transposed": "Flood quantiles transposed from the gage by the drainage-area ratio.",
    "engine_check": "Flood quantiles from flowfreq and from the reference.",
    "low_flow": LOW_FLOW_CAPTION,
    "flow_duration": FLOW_DURATION_CAPTION,
    "future_flow": FUTURE_FLOW_CAPTION,
}

#: Kinds whose opening paragraph is part of the result: a check's verdict is its result.
_VERDICT_KINDS = frozenset({"engine_check"})

#: A table longer than this is named, not printed, when a section is read back.
LONG_TABLE_ROWS = 40

_REPORT_SUFFIXES = (".docx", ".md")
_DETAILS = ("title", "subtitle", "author", "project", "date")


@dataclass
class _Entry:
    """A section in the session: its id, its spec (options as resolved) and what was built."""

    id: str
    spec: SectionSpec
    content: SectionContent

    @property
    def title(self) -> str:
        return self.spec.title or self.content.title


class ReportSession:
    """A report under construction.

    Parameters
    ----------
    title : str
        The report's title; :meth:`set_details` changes it.
    output_dir : str or Path, optional
        Where :meth:`render` and :meth:`save` write and :meth:`load` reads. Defaults to the
        working directory. Files are named, never given by path, so nothing is written or
        read outside it.

    Attributes
    ----------
    output_dir : Path
    """

    def __init__(
        self, title: str = "Untitled report", output_dir: Union[str, Path, None] = None
    ) -> None:
        self.output_dir = Path(output_dir) if output_dir is not None else Path.cwd()
        self._details: Dict[str, Optional[str]] = {key: None for key in _DETAILS}
        self._details["title"] = title
        self._entries: List[_Entry] = []
        self._issued = 0

    # -- state ----------------------------------------------------------------------------

    @property
    def spec(self) -> ReportSpec:
        """The report as it stands, as a specification (a copy)."""
        details = self._details
        return ReportSpec(
            title=details["title"] or "Untitled report",
            sections=[copy.deepcopy(entry.spec) for entry in self._entries],
            subtitle=details["subtitle"],
            author=details["author"],
            project=details["project"],
            date=details["date"],
        )

    def _entry(self, section_id: str) -> _Entry:
        for entry in self._entries:
            if entry.id == section_id:
                return entry
        known = ", ".join(entry.id for entry in self._entries) or "none yet"
        raise ReportError(f"no section {section_id!r} in the report; its sections: {known}")

    def _describe(self, entry: _Entry) -> Dict[str, Any]:
        spec = entry.spec
        if not (spec.narrative and spec.narrative.strip()):
            narrative = "none"
        else:
            narrative = spec.narrative_status
        return {
            "section_id": entry.id,
            "number": self._entries.index(entry) + 1,
            "kind": spec.kind,
            "title": entry.title,
            "narrative": narrative,
        }

    def _index(self, position: Optional[int]) -> int:
        """Where a new section goes; checked before the section is built."""
        if position is None:
            return len(self._entries)
        if not 1 <= position <= len(self._entries) + 1:
            raise ReportError(
                f"position must be from 1 to {len(self._entries) + 1}, got {position}"
            )
        return position - 1

    def _new(self, spec: SectionSpec) -> _Entry:
        """Build a section; nothing in the session changes if it is refused."""
        if spec.narrative and spec.narrative.strip():
            _check_narrative(spec.narrative)
        elif spec.kind == "text":
            raise ReportError("a text section needs a narrative")
        content = build_section(spec.kind, spec.options)
        if not (spec.title or content.title):
            raise ReportError(f"a {spec.kind} section needs a title")
        self._issued += 1
        built = SectionSpec(
            spec.kind,
            title=spec.title,
            options=content.options,
            narrative=spec.narrative,
            narrative_status=spec.narrative_status,
        )
        return _Entry(f"s{self._issued}", built, content)

    def outline(self) -> Dict[str, Any]:
        """The report's title block and its sections in order.

        Returns
        -------
        dict
            ``title``, ``subtitle``, ``author``, ``project``, ``date`` (those set) and
            ``sections``: for each, ``section_id``, ``number``, ``kind``, ``title`` and
            ``narrative`` (``"none"``, ``"draft"`` or ``"reviewed"``).
        """
        return {
            **{key: value for key, value in self._details.items() if value is not None},
            "sections": [self._describe(entry) for entry in self._entries],
        }

    def set_details(
        self,
        *,
        title: Optional[str] = None,
        subtitle: Optional[str] = None,
        author: Optional[str] = None,
        project: Optional[str] = None,
        date: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Change the title block; an argument left ``None`` keeps its value.

        Parameters
        ----------
        title, subtitle, author, project : str, optional
        date : str, optional
            The report date as text. Unset, a report carries the date it is rendered.

        Returns
        -------
        dict
            :meth:`outline`.

        Raises
        ------
        ReportError
            An empty title.
        """
        if title is not None and not title.strip():
            raise ReportError("a report needs a title")
        given = {
            "title": title,
            "subtitle": subtitle,
            "author": author,
            "project": project,
            "date": date,
        }
        for key, value in given.items():
            if value is not None:
                self._details[key] = value.strip() or None
        return self.outline()

    # -- sections -------------------------------------------------------------------------

    def add_section(
        self,
        kind: str,
        options: Mapping[str, Any],
        *,
        title: Optional[str] = None,
        position: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Build a section and add it to the report.

        Parameters
        ----------
        kind : str
            A key of :data:`flowfreq.reporting.SECTION_KINDS`.
        options : mapping
            The kind's inputs.
        title : str, optional
            The heading, in place of the kind's own.
        position : int, optional
            Where it goes, counting from 1. Defaults to the end.

        Returns
        -------
        dict
            ``section_id``, ``number``, ``kind``, ``title``, ``narrative`` and ``result``:
            the section's main table as Markdown (``None`` for a kind without one).

        Raises
        ------
        ReportError
            The section cannot be built as given, or a bad position. Refusals the analyses
            raise themselves pass through, as from
            :func:`~flowfreq.reporting.build_section`. The report is unchanged.
        """
        index = self._index(position)
        entry = self._new(SectionSpec(kind, title=title, options=dict(options)))
        self._entries.insert(index, entry)
        return {**self._describe(entry), "result": _main_table(entry)}

    def add_text(
        self, title: str, narrative: str, *, position: Optional[int] = None
    ) -> Dict[str, Any]:
        """Add a section of prose: paragraphs and bullet lists, printed as a draft.

        Parameters
        ----------
        title : str
        narrative : str
        position : int, optional
            Where it goes, counting from 1. Defaults to the end.

        Returns
        -------
        dict
            ``section_id``, ``number``, ``kind``, ``title``, ``narrative``.

        Raises
        ------
        ReportError
            No title, no text, or text holding a table or a heading.
        """
        if not title.strip():
            raise ReportError("a text section needs a title")
        index = self._index(position)
        entry = self._new(SectionSpec("text", title=title.strip(), narrative=narrative))
        self._entries.insert(index, entry)
        return self._describe(entry)

    def weight(
        self,
        gage_id: str,
        regression_id: str,
        *,
        detail: bool = True,
        keep_parts: bool = False,
    ) -> Dict[str, Any]:
        """Weight a gage section with a regression section already in the report.

        Parameters
        ----------
        gage_id, regression_id : str
            The ``gaged_ffa`` section and the ``ungaged`` section, which must be the
            regression at the gage.
        detail : bool, default True
            Carry both analyses below the weighted table.
        keep_parts : bool, default False
            Keep the two sections as well. By default the weighted section takes the
            gage section's place and both are removed, since it holds what they held.

        Returns
        -------
        dict
            As :meth:`add_section`.

        Raises
        ------
        ReportError
            Sections of the wrong kind, a part with a narrative that removing it would
            drop, or a weighting :func:`~flowfreq.reporting.build_section` refuses.
        """
        gage, regression = self._entry(gage_id), self._entry(regression_id)
        for entry, kind in ((gage, "gaged_ffa"), (regression, "ungaged")):
            if entry.spec.kind != kind:
                raise ReportError(
                    f"section {entry.id} is a {entry.spec.kind} section, not a {kind} section"
                )
        if not keep_parts:
            for entry in (gage, regression):
                if entry.spec.narrative and entry.spec.narrative.strip():
                    raise ReportError(
                        f"section {entry.id} has a narrative, which replacing it would drop; "
                        "keep the parts or move the text first"
                    )
        built = self._new(
            SectionSpec(
                "weighted",
                options={
                    "gage": gage.spec.options,
                    "regression": regression.spec.options,
                    "detail": detail,
                },
            )
        )
        self._entries.insert(self._entries.index(gage), built)
        if not keep_parts:
            self._entries.remove(gage)
            self._entries.remove(regression)
        return {**self._describe(built), "result": _main_table(built)}

    def transpose(
        self,
        gage_id: str,
        target_area_sqmi: float,
        exponents: Mapping[str, Any],
        *,
        target_name: Optional[str] = None,
        options: Optional[Mapping[str, Any]] = None,
        detail: bool = False,
        title: Optional[str] = None,
        position: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Transpose a gage section's quantiles to an ungaged site on the same stream.

        Parameters
        ----------
        gage_id : str
            The ``gaged_ffa`` section that is the donor. It stays in the report.
        target_area_sqmi : float
            The site's drainage area, in square miles.
        exponents : mapping
            The drainage-area exponents and where they came from: ``citation`` and either
            ``exponent`` (one value for every AEP) or ``aeps`` and ``exponents``, with
            ``region`` and ``valid_area_range_sqmi`` optional.
        target_name : str, optional
            A name for the site.
        options : mapping, optional
            The ``transposed`` kind's other options (``donor_area_sqmi``, ``max_aep``,
            ``allow_out_of_range``, ``donor_site_class``, ``target_site_class``,
            ``allow_regulation_mismatch``).
        detail : bool, default False
            Repeat the gage's analysis below the transposed table. Off by default, since
            the gage's own section is still in the report.
        title : str, optional
            The heading, in place of the kind's own.
        position : int, optional
            Where it goes, counting from 1. Defaults to just after the gage's section.

        Returns
        -------
        dict
            As :meth:`add_section`.

        Raises
        ------
        ReportError
            A section of the wrong kind, a bad position, or a transposition
            :func:`~flowfreq.reporting.build_section` refuses. Refusals the analysis
            raises itself (an area ratio outside the method's range, a regulation
            mismatch) pass through. The report is unchanged.
        """
        gage = self._entry(gage_id)
        if gage.spec.kind != "gaged_ffa":
            raise ReportError(
                f"section {gage.id} is a {gage.spec.kind} section, not a gaged_ffa section"
            )
        index = self._entries.index(gage) + 1 if position is None else self._index(position)
        given: Dict[str, Any] = dict(options or {})
        for key in ("donor", "target_area_sqmi", "exponents", "target_name", "detail"):
            if key in given:
                raise ReportError(f"transposed: {key} is an argument, not one of the options")
        given.update(
            donor=gage.spec.options,
            target_area_sqmi=target_area_sqmi,
            exponents=dict(exponents),
            detail=detail,
        )
        if target_name:
            given["target_name"] = target_name
        built = self._new(SectionSpec("transposed", title=title, options=given))
        self._entries.insert(index, built)
        return {**self._describe(built), "result": _main_table(built)}

    def engine_check(
        self,
        gage_id: str,
        *,
        reference: Optional[Mapping[str, Any]] = None,
        tolerance_pct: Optional[float] = None,
        parameter_tolerance_pct: Optional[float] = None,
        ci_tolerance_pct: Optional[float] = None,
        detail: bool = False,
        title: Optional[str] = None,
        position: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Check a gage section's fit against a reference computation of the same record.

        Parameters
        ----------
        gage_id : str
            The ``gaged_ffa`` section to check. It stays in the report.
        reference : mapping, optional
            The reference's numbers: ``source``, ``quantiles`` (objects with ``aep`` and
            ``flow_cfs``, and ``lower_cfs`` and ``upper_cfs`` for 90-percent limits), and
            optionally ``parameters``, ``low_outlier_count`` and
            ``low_outlier_threshold_cfs``. Without it the record is run through USGS
            peakfq 8.1.0 by flowfreq's Fortran extension, and what that returned is stored
            in the section's options.
        tolerance_pct, parameter_tolerance_pct, ci_tolerance_pct : float, optional
            Tolerances in percent for the quantiles (default 1), the mean and standard
            deviation of log10 flow (0.5) and the confidence limits (2).
        detail : bool, default False
            Repeat the gage's analysis below the comparison. Off by default, since the
            gage's own section is still in the report.
        title : str, optional
            The heading, in place of the kind's own.
        position : int, optional
            Where it goes, counting from 1. Defaults to just after the gage's section.

        Returns
        -------
        dict
            As :meth:`add_section`, with the verdict paragraph ahead of the quantile table
            in ``result``.

        Raises
        ------
        ReportError
            A section of the wrong kind, a bad position, a reference
            :func:`~flowfreq.reporting.build_section` refuses, or no reference given and
            the Fortran extension not built. The report is unchanged.
        """
        gage = self._entry(gage_id)
        if gage.spec.kind != "gaged_ffa":
            raise ReportError(
                f"section {gage.id} is a {gage.spec.kind} section, not a gaged_ffa section"
            )
        index = self._entries.index(gage) + 1 if position is None else self._index(position)
        given: Dict[str, Any] = {"gage": gage.spec.options, "detail": detail}
        if reference is not None:
            given["reference"] = dict(reference)
        for name, value in (
            ("tolerance_pct", tolerance_pct),
            ("parameter_tolerance_pct", parameter_tolerance_pct),
            ("ci_tolerance_pct", ci_tolerance_pct),
        ):
            if value is not None:
                given[name] = value
        built = self._new(SectionSpec("engine_check", title=title, options=given))
        self._entries.insert(index, built)
        return {**self._describe(built), "result": _main_table(built)}

    def future_flow(
        self,
        base_id: str,
        factor_set: Mapping[str, Any],
        *,
        aeps: Optional[Sequence[float]] = None,
        detail: bool = False,
        title: Optional[str] = None,
        position: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Scale a section's flood quantiles by change factors for a future condition.

        Parameters
        ----------
        base_id : str
            The ``gaged_ffa`` or ``ungaged`` section whose quantiles are scaled. It stays in
            the report.
        factor_set : mapping
            The change factors and where they come from: ``name``, ``aeps``, ``factors``,
            ``scenario``, ``horizon``, ``citation``, and optionally ``level``,
            ``geography`` and ``legal_status``
            (:class:`~flowfreq.future_flow.ChangeFactorSet`).
        aeps : sequence of float, optional
            The AEPs to scale. Defaults to the section's AEPs within the factors' range.
        detail : bool, default False
            Repeat the present-condition analysis below the table.
        title : str, optional
            The heading, in place of the kind's own.
        position : int, optional
            Where it goes, counting from 1. Defaults to just after the base section.

        Returns
        -------
        dict
            As :meth:`add_section`.

        Raises
        ------
        ReportError
            A section of the wrong kind, a bad position, or factors
            :func:`~flowfreq.reporting.build_section` refuses. The report is unchanged.
        """
        base = self._entry(base_id)
        names = {"gaged_ffa": "gage", "ungaged": "regression"}
        if base.spec.kind not in names:
            raise ReportError(
                f"section {base.id} is a {base.spec.kind} section, not a gaged_ffa or an "
                "ungaged section"
            )
        index = self._entries.index(base) + 1 if position is None else self._index(position)
        given: Dict[str, Any] = {
            names[base.spec.kind]: base.spec.options,
            "factor_set": dict(factor_set),
            "detail": detail,
        }
        if aeps is not None:
            given["aeps"] = [float(a) for a in aeps]
        built = self._new(SectionSpec("future_flow", title=title, options=given))
        self._entries.insert(index, built)
        return {**self._describe(built), "result": _main_table(built)}

    def section_markdown(self, section_id: str) -> str:
        """A section as it will print, as Markdown.

        Parameters
        ----------
        section_id : str

        Returns
        -------
        str
            The heading, the computed content and the narrative. A figure is its caption
            and a table of more than :data:`LONG_TABLE_ROWS` rows is its caption and its
            size; both are whole in the rendered report.

        Raises
        ------
        ReportError
            An unknown section.
        """
        entry = self._entry(section_id)
        blocks: List[Block] = [Heading(1, entry.title)]
        for block in entry.content.blocks:
            if isinstance(block, Table) and len(block.rows) > LONG_TABLE_ROWS:
                what = block.caption or ", ".join(block.header)
                blocks.append(
                    Paragraph(f"(Table of {len(block.rows)} rows, in the report: {what})")
                )
            elif isinstance(block, Figure):
                blocks.append(Paragraph(f"(Figure, in the report: {block.caption})"))
            else:
                blocks.append(block)
        spec = entry.spec
        if spec.narrative and spec.narrative.strip():
            blocks.append(Narrative(spec.narrative, spec.narrative_status))
        return render_markdown(blocks)

    def set_narrative(self, section_id: str, narrative: str) -> Dict[str, Any]:
        """Set a section's prose, which prints after its computed content as a draft.

        Parameters
        ----------
        section_id : str
        narrative : str
            Paragraphs and bullet lists. Empty text removes the narrative.

        Returns
        -------
        dict
            ``section_id``, ``number``, ``kind``, ``title``, ``narrative``.

        Raises
        ------
        ReportError
            An unknown section, text holding a table or a heading, or empty text for a
            text section, which is nothing else.
        """
        entry = self._entry(section_id)
        if narrative.strip():
            _check_narrative(narrative)
        elif entry.spec.kind == "text":
            raise ReportError("a text section needs a narrative; remove the section instead")
        entry.spec.narrative = narrative if narrative.strip() else None
        # New text is unread text, whatever the old text's status was.
        entry.spec.narrative_status = "draft"
        return self._describe(entry)

    def mark_reviewed(self, section_id: str) -> Dict[str, Any]:
        """Record that a person has read a section's narrative, so it prints unlabeled.

        Parameters
        ----------
        section_id : str

        Returns
        -------
        dict
            ``section_id``, ``number``, ``kind``, ``title``, ``narrative``.

        Raises
        ------
        ReportError
            An unknown section, or one with no narrative.
        """
        entry = self._entry(section_id)
        if not (entry.spec.narrative and entry.spec.narrative.strip()):
            raise ReportError(f"section {section_id} has no narrative to review")
        entry.spec.narrative_status = "reviewed"
        return self._describe(entry)

    def remove_section(self, section_id: str) -> Dict[str, Any]:
        """Take a section out of the report.

        Returns
        -------
        dict
            :meth:`outline`.

        Raises
        ------
        ReportError
            An unknown section.
        """
        self._entries.remove(self._entry(section_id))
        return self.outline()

    def move_section(self, section_id: str, position: int) -> Dict[str, Any]:
        """Move a section to a position, counting from 1.

        Returns
        -------
        dict
            :meth:`outline`.

        Raises
        ------
        ReportError
            An unknown section or a position outside the report.
        """
        entry = self._entry(section_id)
        if not 1 <= position <= len(self._entries):
            raise ReportError(f"position must be from 1 to {len(self._entries)}, got {position}")
        self._entries.remove(entry)
        self._entries.insert(position - 1, entry)
        return self.outline()

    # -- files ----------------------------------------------------------------------------

    def _path(self, filename: str, suffixes: Any) -> Path:
        name = Path(filename)
        if name.name != filename or not name.stem or name.suffix.lower() not in suffixes:
            raise ReportError(
                f"give a file name ending in {' or '.join(suffixes)}, with no directory "
                f"(files go in {self.output_dir}); got {filename!r}"
            )
        return self.output_dir / name.name

    def render(self, filename: str) -> Dict[str, Any]:
        """Write the report, and beside it the specification that rebuilds it.

        Parameters
        ----------
        filename : str
            A name ending in ``.docx`` (needs the ``report`` extra) or ``.md``, written in
            :attr:`output_dir`. The specification takes the same name with ``.json``.

        Returns
        -------
        dict
            ``report`` and ``specification`` (the paths written), ``sections`` and
            ``draft_narratives`` (how many print under the draft label).

        Raises
        ------
        ReportError
            A report with no sections, or a bad file name.
        ImportError
            ``.docx`` without python-docx installed.
        """
        path = self._path(filename, _REPORT_SUFFIXES)
        if not self._entries:
            raise ReportError("the report has no sections yet")
        report = build_report(self.spec)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        written = report.save(path)
        specification = report.spec.save(path.with_suffix(".json"))
        return {
            "report": str(written),
            "specification": str(specification),
            "sections": len(self._entries),
            "draft_narratives": sum(
                self._describe(entry)["narrative"] == "draft" for entry in self._entries
            ),
        }

    def save(self, filename: str) -> Dict[str, Any]:
        """Write the report's specification, to pick the report up again with :meth:`load`.

        Parameters
        ----------
        filename : str
            A name ending in ``.json``, written in :attr:`output_dir`.

        Returns
        -------
        dict
            ``specification`` (the path written) and ``sections``.

        Raises
        ------
        ReportError
            A bad file name.
        """
        path = self._path(filename, (".json",))
        self.output_dir.mkdir(parents=True, exist_ok=True)
        return {"specification": str(self.spec.save(path)), "sections": len(self._entries)}

    def load(self, filename: str) -> Dict[str, Any]:
        """Replace the report with a saved specification, building every section of it.

        Parameters
        ----------
        filename : str
            A ``.json`` name in :attr:`output_dir`.

        Returns
        -------
        dict
            :meth:`outline`.

        Raises
        ------
        ReportError
            A bad file name, no such file, a file that is not a report specification, or
            a section that cannot be built. The report in the session is unchanged.
        """
        path = self._path(filename, (".json",))
        if not path.is_file():
            raise ReportError(f"no specification {path.name} in {self.output_dir}")
        loaded = ReportSpec.load(path)
        issued = self._issued
        entries: List[_Entry] = []
        try:
            for number, section in enumerate(loaded.sections, 1):
                try:
                    entries.append(self._new(section))
                except ReportError as exc:
                    raise ReportError(f"section {number} ({section.kind}): {exc}") from exc
        except Exception:
            self._issued = issued
            raise
        self._entries = entries
        self._details = {key: getattr(loaded, key) for key in _DETAILS}
        return self.outline()


def _check_narrative(text: str) -> None:
    try:
        narrative_blocks(text)
    except ValueError as exc:
        raise ReportError(str(exc)) from exc


def _main_table(entry: _Entry) -> Optional[str]:
    caption = MAIN_TABLES.get(entry.spec.kind)
    lead = entry.content.blocks[:1] if entry.spec.kind in _VERDICT_KINDS else []
    for block in entry.content.blocks:
        if isinstance(block, Table) and block.caption == caption:
            return render_markdown([*lead, block])
    return None


__all__ = ["LONG_TABLE_ROWS", "MAIN_TABLES", "ReportSession"]
