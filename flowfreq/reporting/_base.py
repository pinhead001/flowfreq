"""What every section kind shares: the content a builder returns, formatting, figures."""

from __future__ import annotations

import io
import math
import re
from dataclasses import dataclass, field, replace
from typing import Any, Callable, Dict, List, Mapping, Sequence, Tuple

from .blocks import Block, Heading
from .spec import ReportError

#: The reference every Bulletin 17C section cites.
BULLETIN_17C = (
    "England, J.F., Jr., Cohn, T.A., Faber, B.A., Stedinger, J.R., Thomas, W.O., Jr., "
    "Veilleux, A.G., Kiang, J.E., and Mason, R.R., Jr., 2019, Guidelines for determining "
    "flood flow frequency, Bulletin 17C (ver. 1.1, May 2019): U.S. Geological Survey "
    "Techniques and Methods, book 4, chap. B5, 148 p., https://doi.org/10.3133/tm4B5."
)


_PEAK_KEYS = (
    "water_year",
    "peak_flow_cfs",
    "qualification_code",
    "peak_date",
    "source",
    "retrieved",
)


_FIGURE_DPI = 200


#: The authors and year a USGS-style citation opens with: "Eash, D.A., and Veilleux, A.G., 2013,".
_AUTHORS_YEAR = re.compile(r"^(?P<authors>.+?),\s*(?P<year>(?:18|19|20)\d{2})[a-z]?,")


_INITIALS = re.compile(r"(?:[A-Z]\.[ -]?)+|Jr\.?|Sr\.?|II|III")


@dataclass(frozen=True)
class SectionContent:
    """What a section builder returns.

    Attributes
    ----------
    title : str
        The section's own heading, used when the spec gives none.
    blocks : list of Block
        Tables, figures and boilerplate text, all computed.
    options : dict
        The options as resolved: enough to build the same section offline.
    references : list of str
        Full citations for the report's reference list.
    cautions : list of str
        What the regression behind the flows warns of (``SiteEstimate.advisory``), for a
        reply that shows only the section's table. Already in ``blocks``.
    """

    title: str
    blocks: List[Block]
    options: Dict[str, Any]
    references: List[str]
    cautions: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class SectionKind:
    """A registered section kind: its name, what it takes and the function that builds it."""

    name: str
    summary: str
    required: Tuple[str, ...]
    optional: Tuple[str, ...]
    build: Callable[[Dict[str, Any]], SectionContent]

    def check(self, options: Mapping[str, Any]) -> None:
        """Raise :class:`ReportError` for an unknown or a missing option."""
        allowed = (*self.required, *self.optional)
        extra = sorted(set(options) - set(allowed))
        if extra:
            raise ReportError(
                f"{self.name}: unknown option(s) {', '.join(extra)}; "
                f"allowed: {', '.join(allowed) or 'none'}"
            )
        missing = [key for key in self.required if options.get(key) is None]
        if missing:
            raise ReportError(f"{self.name}: missing option(s) {', '.join(missing)}")


#: The section kinds a spec may name. Filled by :mod:`flowfreq.reporting.sections`, and by
#: the modules that register further kinds (``daily``, ``summary``).
SECTION_KINDS: Dict[str, SectionKind] = {}


def _version() -> str:
    from .. import __version__

    return str(__version__)


def cfs(value: Any) -> str:
    """A flow for a table cell: whole cfs at 100 and above, three significant figures below."""
    if value is None:
        return "-"
    number = float(value)
    if not math.isfinite(number):
        return "-"
    return f"{number:,.0f}" if abs(number) >= 100 else f"{number:.3g}"


def short_citation(citation: str) -> str:
    """``Eash and others (2013)`` for a citation that opens with authors and a year.

    One author gives ``Eash (2013)`` and two ``Southard and Veilleux (2014)``. A citation in
    any other form is returned whole, so nothing is guessed.
    """
    match = _AUTHORS_YEAR.match(citation.strip())
    if match is None:
        return citation.strip()
    # "Eash, D.A., Barnes, K.K., and Veilleux, A.G.": the surnames are the parts that are
    # not initials.
    parts = [p.strip().removeprefix("and ").strip() for p in match["authors"].split(",")]
    surnames = [p for p in parts if p and not _INITIALS.fullmatch(p)]
    if not surnames:
        return citation.strip()
    if len(surnames) == 1:
        who = surnames[0]
    elif len(surnames) == 2:
        who = f"{surnames[0]} and {surnames[1]}"
    else:
        who = f"{surnames[0]} and others"
    return f"{who} ({match['year']})"


def _figure_png(fig: Any) -> bytes:
    """A pyplot figure as PNG bytes; the figure is closed."""
    import matplotlib.pyplot as plt

    buffer = io.BytesIO()
    # No Software entry: it carries the matplotlib version, which is not report content.
    fig.savefig(
        buffer, format="png", dpi=_FIGURE_DPI, bbox_inches="tight", metadata={"Software": None}
    )
    plt.close(fig)
    return buffer.getvalue()


def _demote(blocks: Sequence[Block], by: int = 1) -> List[Block]:
    """Push every heading ``by`` levels down, to nest a parsed page under a subsection."""
    return [replace(b, level=b.level + by) if isinstance(b, Heading) else b for b in blocks]


# -- text ---------------------------------------------------------------------------------


def _build_text(options: Dict[str, Any]) -> SectionContent:
    return SectionContent(title="", blocks=[], options=dict(options), references=[])
