"""The report specification: what a report contains, as plain data that saves to JSON.

A spec records inputs, never results. Building it (:func:`flowfreq.reporting.build_report`)
computes every section, and returns the spec *as built*: anything a section downloaded (a
gage's annual peaks) is written into its options, so saving that spec and building it again
gives the same report without the network.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from .blocks import NARRATIVE_STATUSES

#: Version of the JSON layout. Bumped when a saved spec would be read differently.
SCHEMA = 1

_SECTION_KEYS = ("kind", "title", "options", "narrative", "narrative_status")
_REPORT_KEYS = ("schema", "title", "subtitle", "author", "project", "date", "sections")


class ReportError(ValueError):
    """A report cannot be built as specified: the message says what to change."""


def _unknown(given: Any, allowed: Any, what: str) -> None:
    extra = sorted(set(given) - set(allowed))
    if extra:
        raise ReportError(
            f"{what}: unknown key(s) {', '.join(extra)}; allowed: {', '.join(allowed)}"
        )


@dataclass
class SectionSpec:
    """One section of a report.

    Attributes
    ----------
    kind : str
        A key of :data:`flowfreq.reporting.SECTION_KINDS`.
    title : str, optional
        The section heading. ``None`` takes the kind's own (``text`` sections need one).
    options : dict
        The kind's inputs, JSON-serializable.
    narrative : str, optional
        Prose placed after the section's computed content: paragraphs and bullet lists.
    narrative_status : {"draft", "reviewed"}
        ``draft`` (the default) renders under a label saying a language model wrote it and
        nobody has reviewed it. Set ``reviewed`` once a person has.
    """

    kind: str
    title: Optional[str] = None
    options: Dict[str, Any] = field(default_factory=dict)
    narrative: Optional[str] = None
    narrative_status: str = "draft"

    def __post_init__(self) -> None:
        if self.narrative_status not in NARRATIVE_STATUSES:
            raise ReportError(
                f"narrative_status must be one of {', '.join(NARRATIVE_STATUSES)}, "
                f"got {self.narrative_status!r}"
            )

    def to_dict(self) -> Dict[str, Any]:
        """Plain-dict form, omitting what is unset."""
        out: Dict[str, Any] = {"kind": self.kind}
        if self.title is not None:
            out["title"] = self.title
        if self.options:
            out["options"] = self.options
        if self.narrative:
            out["narrative"] = self.narrative
            out["narrative_status"] = self.narrative_status
        return out

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SectionSpec":
        """Read :meth:`to_dict`'s form. Raises :class:`ReportError` on an unknown key."""
        if not isinstance(data, dict) or "kind" not in data:
            raise ReportError(f"a section needs a 'kind'; got {data!r}")
        _unknown(data, _SECTION_KEYS, f"section {data['kind']!r}")
        options = data.get("options") or {}
        if not isinstance(options, dict):
            raise ReportError(f"section {data['kind']!r}: options must be an object")
        return cls(
            kind=str(data["kind"]),
            title=data.get("title"),
            options=dict(options),
            narrative=data.get("narrative"),
            narrative_status=data.get("narrative_status", "draft"),
        )


@dataclass
class ReportSpec:
    """A report: its title block and its sections, in order.

    Attributes
    ----------
    title : str
    sections : list of SectionSpec
    subtitle, author, project : str, optional
        Shown under the title.
    date : str, optional
        The report date as text (ISO ``YYYY-MM-DD`` suggested). ``None`` takes the build date,
        which the spec as built then records.
    """

    title: str
    sections: List[SectionSpec] = field(default_factory=list)
    subtitle: Optional[str] = None
    author: Optional[str] = None
    project: Optional[str] = None
    date: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        """Plain-dict form with the schema version, omitting what is unset."""
        out: Dict[str, Any] = {"schema": SCHEMA, "title": self.title}
        for key in ("subtitle", "author", "project", "date"):
            value = getattr(self, key)
            if value is not None:
                out[key] = value
        out["sections"] = [s.to_dict() for s in self.sections]
        return out

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ReportSpec":
        """Read :meth:`to_dict`'s form.

        Raises
        ------
        ReportError
            No title, an unknown key, or a schema newer than this flowfreq reads.
        """
        if not isinstance(data, dict):
            raise ReportError("a report spec is a JSON object")
        _unknown(data, _REPORT_KEYS, "report spec")
        schema = data.get("schema", SCHEMA)
        if schema != SCHEMA:
            raise ReportError(
                f"report spec schema {schema!r} is not the schema this flowfreq reads ({SCHEMA})"
            )
        if not data.get("title"):
            raise ReportError("a report spec needs a title")
        sections = data.get("sections") or []
        if not isinstance(sections, list):
            raise ReportError("report spec: sections must be a list")
        return cls(
            title=str(data["title"]),
            sections=[SectionSpec.from_dict(s) for s in sections],
            subtitle=data.get("subtitle"),
            author=data.get("author"),
            project=data.get("project"),
            date=data.get("date"),
        )

    def to_json(self) -> str:
        """JSON text. Raises ``ValueError`` if an option holds NaN or a non-JSON value."""
        return json.dumps(self.to_dict(), indent=1, ensure_ascii=False, allow_nan=False) + "\n"

    @classmethod
    def from_json(cls, text: str) -> "ReportSpec":
        """Read :meth:`to_json`'s text."""
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ReportError(f"report spec is not valid JSON: {exc}") from exc
        return cls.from_dict(data)

    def save(self, path: Union[str, Path]) -> Path:
        """Write the spec as UTF-8 JSON and return the path."""
        out = Path(path)
        out.write_text(self.to_json(), encoding="utf-8")
        return out

    @classmethod
    def load(cls, path: Union[str, Path]) -> "ReportSpec":
        """Read a spec saved by :meth:`save`."""
        return cls.from_json(Path(path).read_text(encoding="utf-8"))
