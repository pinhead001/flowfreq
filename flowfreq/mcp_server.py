"""flowfreq's report tools as a Model Context Protocol server: a report built through a chat.

A chat client (Claude Code, Claude Desktop) starts this server and calls its tools. The
server holds one :class:`~flowfreq.reporting.session.ReportSession`, so the report takes
shape over the conversation: find a gage, fit it, estimate at a site, weight the two, write
a few paragraphs, render the DOCX. The division of labor is the reporting package's. Every
table and figure is computed here from stored inputs; the model contributes prose, which
prints under a draft label until the user says they have read it.

Run it with ``flowfreq-mcp --dir <output directory>`` (``pip install "flowfreq[mcp]"``);
``docs/CHAT_REPORTS.md`` has the client configuration.

>>> from flowfreq.mcp_server import build_server
>>> server = build_server()            # the tools over a new, empty report
>>> server.run("stdio")                # what flowfreq-mcp does
"""

from __future__ import annotations

import argparse
import functools
import inspect
import logging
import os
import sys
import threading
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Union

try:
    from mcp.server.mcpserver import MCPServer
    from mcp.server.mcpserver.exceptions import ToolError
except ImportError as exc:  # pragma: no cover - exercised only without the extra
    raise ImportError(
        'flowfreq.mcp_server needs the mcp package (2.3 or later): pip install "flowfreq[mcp]"'
    ) from exc

import pandas as pd

from . import __version__
from .reporting.session import ReportSession
from .reporting.spec import ReportError

logger = logging.getLogger(__name__)

#: The environment variable naming the directory reports are written to.
REPORT_DIR_ENV = "FLOWFREQ_REPORT_DIR"

#: What the server tells a client's model about using it.
INSTRUCTIONS = """\
flowfreq builds a hydrology report section by section and writes it as a Word document.

Rules:
- Every number in the report comes from a tool. Do not state a flow, a probability or a \
statistic that a tool did not return, and do not compute one yourself.
- Tables and figures are made by the tools. A narrative holds paragraphs and bullet lists \
only; a narrative with a table or a heading is refused.
- Narrative you write is saved as a draft and prints under a label saying a language model \
wrote it and nobody has reviewed it. Call mark_narrative_reviewed only when the user says \
they have read that section's text.
- When a tool refuses (a regulated record, a characteristic outside an equation's range, \
no regional skew), tell the user what it said. The allow_* arguments override a refusal; \
set one only when the user asks for it.

A usual order: find_gage, run_ffa, then estimate_ungaged or estimate_from_characteristics \
at the gage and weight_gage_with_regression, then get_section to read what was computed, \
edit_narrative, and render_report. The report is kept for the life of the server; \
save_report and load_report carry it between conversations.
"""

_GAGE_FIELDS = (
    "site_no",
    "site_name",
    "state",
    "latitude",
    "longitude",
    "drainage_area_sqmi",
    "n_peaks",
    "first_water_year",
    "last_water_year",
    "regulation_class",
)
_SKEW_CHOICES = ("regional", "station", "default")


def default_report_dir() -> Path:
    """Where reports go when the server is given no directory.

    Returns
    -------
    Path
        ``$FLOWFREQ_REPORT_DIR`` if set, otherwise ``flowfreq-reports`` in the home directory.
    """
    named = os.environ.get(REPORT_DIR_ENV)
    return Path(named).expanduser() if named else Path.home() / "flowfreq-reports"


def _reported(fn: Callable[..., Any], lock: threading.Lock) -> Callable[..., Any]:
    """A tool as the server runs it: one at a time, and refusing in words a client sees.

    The server runs each call on a worker thread, and a model may send several calls at
    once. The report and matplotlib's pyplot are both shared state, so ``lock`` makes the
    calls take turns. Whatever a tool raises becomes a ``ToolError``, the one error whose
    text the server passes on.
    """

    @functools.wraps(fn)
    def call(*args: Any, **kwargs: Any) -> Any:
        try:
            with lock:
                return fn(*args, **kwargs)
        except ToolError:
            raise
        except Exception as exc:
            # A refusal is an answer; anything else is a fault worth a traceback in the log.
            refusal = isinstance(exc, (ValueError, LookupError, ImportError, OSError))
            logger.warning("%s: %s", fn.__name__, exc, exc_info=not refusal)
            text = str(exc.args[0]) if isinstance(exc, KeyError) and exc.args else str(exc)
            raise ToolError(f"{type(exc).__name__}: {text}") from exc

    return call


def _plain(value: Any) -> Any:
    """A catalog cell as JSON: missing is ``None``, a whole number an ``int``."""
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return None
    if isinstance(value, str):
        return value
    number = float(value)
    return int(number) if number.is_integer() else number


def _gaged_options(
    site_no: str,
    skew: str,
    regional_skew: Optional[float],
    regional_skew_se: Optional[float],
    regional_skew_citation: Optional[str],
    allow_provisional_skew: bool,
) -> Dict[str, Any]:
    """The ``gaged_ffa`` options for a station and a choice of skew."""
    if not site_no.strip().isdigit():
        raise ReportError(f"site_no is a USGS station number, all digits; got {site_no!r}")
    if skew not in _SKEW_CHOICES:
        raise ReportError(f"skew must be one of {', '.join(_SKEW_CHOICES)}; got {skew!r}")
    options: Dict[str, Any] = {"site_no": site_no.strip()}
    given = regional_skew is not None or regional_skew_se is not None
    if skew != "regional":
        if given:
            raise ReportError(f"a regional_skew value goes with skew='regional', not {skew!r}")
        options["station_skew_only" if skew == "station" else "use_default_skew"] = True
    elif given:
        if regional_skew is None or regional_skew_se is None:
            raise ReportError("give regional_skew and regional_skew_se together")
        options["regional_skew"] = float(regional_skew)
        options["regional_skew_se"] = float(regional_skew_se)
        if regional_skew_citation:
            options["regional_skew_citation"] = regional_skew_citation
    else:
        options["regional_skew_from_site"] = True
        if allow_provisional_skew:
            options["allow_provisional_skew"] = True
    return options


def build_server(
    session: Optional[ReportSession] = None,
    *,
    output_dir: Union[str, Path, None] = None,
) -> MCPServer:
    """The MCP server: flowfreq's report tools over one report.

    Parameters
    ----------
    session : ReportSession, optional
        The report the tools work on. Defaults to a new, empty one.
    output_dir : str or Path, optional
        Where a new session writes; ignored when ``session`` is given. Defaults to
        :func:`default_report_dir`.

    Returns
    -------
    mcp.server.mcpserver.MCPServer
        Run it with ``server.run("stdio")``.

    Raises
    ------
    ImportError
        On import of this module, when the ``mcp`` package is not installed.
    """
    report = (
        session
        if session is not None
        else ReportSession(
            output_dir=output_dir if output_dir is not None else default_report_dir()
        )
    )
    # Figures are drawn on the server's worker threads and only ever saved. An interactive
    # backend (TkAgg, the default on Windows) makes a window for each figure, and a Tk window
    # made off the main thread takes the process down when it is collected.
    import matplotlib

    matplotlib.use("Agg")
    server = MCPServer("flowfreq", instructions=INSTRUCTIONS, version=__version__)
    lock = threading.Lock()

    def find_gage(
        query: str,
        state: Optional[str] = None,
        min_years: Optional[int] = None,
        limit: int = 10,
    ) -> Dict[str, Any]:
        """Find USGS streamgages with annual peak records, by station number or name.

        query: digits match the start of the station number; anything else is words that
        must all appear in the station name, in any order (for example "big sandy bruceton").
        state: two-letter postal code. min_years: fewest annual peaks. limit: most rows.

        Returns the matches, longest record first: station number, name, state, latitude,
        longitude, drainage area (square miles), number of peaks, first and last water
        year, and regulation class (reference, regulated, urban or unknown). The catalog
        is a file that ships with flowfreq; nothing is downloaded.
        """
        from .catalog import load_catalog, search

        found = search(load_catalog(), query, state=state, min_years=min_years, limit=limit)
        return {
            "count": len(found),
            "gages": [
                {key: _plain(row[key]) for key in _GAGE_FIELDS} for _, row in found.iterrows()
            ],
        }

    def report_outline() -> Dict[str, Any]:
        """The report so far: its title block, and each section's id, number, kind, title
        and narrative status (none, draft or reviewed). Section ids (s1, s2, ...) are what
        the other tools take; they do not change when sections move.
        """
        return report.outline()

    def set_report_details(
        title: Optional[str] = None,
        subtitle: Optional[str] = None,
        author: Optional[str] = None,
        project: Optional[str] = None,
        date: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Set the report's title block. An argument left out keeps its value; an empty
        string clears it (the title cannot be empty). date is text, ISO YYYY-MM-DD
        suggested; unset, the report carries the day it is rendered. Returns the outline.
        """
        return report.set_details(
            title=title, subtitle=subtitle, author=author, project=project, date=date
        )

    def run_ffa(
        site_no: str,
        skew: str = "regional",
        regional_skew: Optional[float] = None,
        regional_skew_se: Optional[float] = None,
        regional_skew_citation: Optional[str] = None,
        allow_provisional_skew: bool = False,
        allow_regulated: bool = False,
        low_outlier_threshold: Optional[float] = None,
        stationarity: bool = True,
        peaks_table: bool = True,
        title: Optional[str] = None,
        position: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Add a Bulletin 17C flood-frequency section for a USGS gage.

        Downloads the gage's annual peaks from USGS, fits the log-Pearson Type III
        distribution by the Expected Moments Algorithm with the Multiple Grubbs-Beck
        low-outlier test, and adds a section with the record, the method, the fitted
        parameters, the flood quantiles with 90-percent confidence limits, the frequency
        curve and a stationarity screen. The peaks are stored in the report, so it
        rebuilds later without a download.

        site_no: the USGS station number (find_gage looks it up).
        skew: "regional" (the default) weights the station skew with the regional skew
        flowfreq has on file for the gage's location, or with regional_skew and
        regional_skew_se (its standard error) when you give them, citing
        regional_skew_citation; "station" uses the station skew alone; "default" uses
        an unsourced generalized skew of -0.302 and is a last resort.
        allow_provisional_skew: accept a regional skew flowfreq has not verified.
        allow_regulated: fit a record flowfreq screens as regulated, which it otherwise
        refuses because Bulletin 17C assumes unregulated peaks.
        low_outlier_threshold: a low-outlier threshold in cfs, in place of the test.
        stationarity, peaks_table: include the trend screen and the table of peaks.
        title: the section heading, in place of the station's name. position: where the
        section goes, counting from 1 (default: the end).

        Returns the section's id and the quantile table. Call get_section for the rest.
        """
        options = _gaged_options(
            site_no,
            skew,
            regional_skew,
            regional_skew_se,
            regional_skew_citation,
            allow_provisional_skew,
        )
        if allow_regulated:
            options["allow_regulated"] = True
        if low_outlier_threshold is not None:
            options["low_outlier_threshold"] = float(low_outlier_threshold)
        options["stationarity"] = stationarity
        options["peaks_table"] = peaks_table
        return report.add_section("gaged_ffa", options, title=title, position=position)

    def estimate_ungaged(
        state: str,
        lat: float,
        lon: float,
        site_name: Optional[str] = None,
        region: Optional[str] = None,
        allow_extrapolation: bool = False,
        allow_provisional: bool = False,
        title: Optional[str] = None,
        position: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Add a regional regression estimate at a point, with the basin from StreamStats.

        Delineates the basin above the point with USGS StreamStats, takes the basin
        characteristics and the regression region it reports, and evaluates the state's
        published regression equations: flood quantiles with prediction limits. What
        StreamStats returned is stored in the report, so it rebuilds later without a
        request. This takes a minute or more, and StreamStats is sometimes unavailable;
        estimate_from_characteristics needs no service when the characteristics are known.

        state: two-letter postal code. lat, lon: decimal degrees (lon negative in the
        Americas), on the stream. site_name: a name for the site. region: a regression
        region code, when the basin lies in several and one is to be used.
        allow_extrapolation: evaluate equations outside their calibrated range, flagged
        in the table. allow_provisional: use equations flowfreq has not verified.

        Returns the section's id and the quantile table.
        """
        options: Dict[str, Any] = {"state": state, "lat": lat, "lon": lon}
        if site_name:
            options["site_name"] = site_name
        if region:
            options["region"] = region
        if allow_extrapolation:
            options["allow_extrapolation"] = True
        if allow_provisional:
            options["allow_provisional"] = True
        return report.add_section("ungaged", options, title=title, position=position)

    def estimate_from_characteristics(
        state: str,
        lat: float,
        lon: float,
        characteristics: Dict[str, float],
        region: str,
        site_name: Optional[str] = None,
        allow_extrapolation: bool = False,
        allow_provisional: bool = False,
        title: Optional[str] = None,
        position: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Add a regional regression estimate from basin characteristics the user supplies.

        Evaluates the state's published regression equations for one region: flood
        quantiles with prediction limits. No service is called.

        state: two-letter postal code. lat, lon: the site, decimal degrees.
        characteristics: values by StreamStats code, for example {"DRNAREA": 205.0}
        (drainage area in square miles); a refusal names a code that is missing.
        region: the regression region's code in flowfreq's library; a refusal for an
        unknown code lists the state's regions.
        allow_extrapolation, allow_provisional, title, position: as estimate_ungaged.

        Returns the section's id and the quantile table.
        """
        options: Dict[str, Any] = {
            "state": state,
            "lat": lat,
            "lon": lon,
            "characteristics": dict(characteristics),
            "region": region,
        }
        if site_name:
            options["site_name"] = site_name
        if allow_extrapolation:
            options["allow_extrapolation"] = True
        if allow_provisional:
            options["allow_provisional"] = True
        return report.add_section("ungaged", options, title=title, position=position)

    def weight_gage_with_regression(
        gage_section_id: str,
        regression_section_id: str,
        detail: bool = True,
        keep_parts: bool = False,
    ) -> Dict[str, Any]:
        """Weight a gage's flood quantiles with the regression estimate at the gage.

        Bulletin 17C Appendix 9: the two estimates are averaged in log space by the
        inverse of their variances, which narrows the confidence limits. Both sections
        must already be in the report: a run_ffa section, and a regression estimate made
        at the gage's own location and drainage area.

        detail: carry both analyses below the weighted table. keep_parts: keep the two
        sections as well; by default the weighted section takes their place, since it
        holds what they held (a part with a narrative is not replaced).

        Returns the new section's id and the weighted table.
        """
        return report.weight(
            gage_section_id, regression_section_id, detail=detail, keep_parts=keep_parts
        )

    def add_text_section(
        title: str, narrative: str, position: Optional[int] = None
    ) -> Dict[str, Any]:
        """Add a section of prose only (an introduction, a summary, a description of the
        basin). narrative is paragraphs and bullet lists, with **bold** and *italic*; no
        tables and no headings. It prints as a draft until marked reviewed. position
        counts from 1 (default: the end).
        """
        return report.add_text(title, narrative, position=position)

    def edit_narrative(section_id: str, narrative: str) -> Dict[str, Any]:
        """Set the prose that follows a section's computed tables and figures, replacing
        any there. Paragraphs and bullet lists only. Read the section with get_section
        first and take every number from it. New text is a draft again, whatever the old
        text was. An empty narrative removes it.
        """
        return report.set_narrative(section_id, narrative)

    def mark_narrative_reviewed(section_id: str) -> Dict[str, Any]:
        """Record that the user has read a section's narrative, so it prints without the
        draft label. Call this only when the user says so for that section.
        """
        return report.mark_reviewed(section_id)

    def get_section(section_id: str) -> str:
        """A section as it will print, in Markdown: every table flowfreq computed, and the
        narrative. Figures and very long tables are named, not shown.
        """
        return report.section_markdown(section_id)

    def remove_section(section_id: str) -> Dict[str, Any]:
        """Take a section out of the report. Returns the outline."""
        return report.remove_section(section_id)

    def move_section(section_id: str, position: int) -> Dict[str, Any]:
        """Move a section to a position, counting from 1. Returns the outline."""
        return report.move_section(section_id, position)

    def render_report(filename: str) -> Dict[str, Any]:
        """Write the report to the output directory. filename is a plain file name ending
        in .docx (Word) or .md (Markdown), with no directory. The specification that
        rebuilds the report is written beside it as .json. Returns both paths and how
        many narratives still print as drafts; tell the user where the file is.
        """
        return report.render(filename)

    def save_report(filename: str) -> Dict[str, Any]:
        """Save the report's specification (a .json file name, no directory) so a later
        conversation can pick it up with load_report.
        """
        return report.save(filename)

    def load_report(filename: str) -> Dict[str, Any]:
        """Replace the report with one saved by save_report or render_report (a .json file
        name in the output directory), rebuilding every section. Returns the outline.
        """
        return report.load(filename)

    tools: Sequence[Callable[..., Any]] = (
        find_gage,
        report_outline,
        set_report_details,
        run_ffa,
        estimate_ungaged,
        estimate_from_characteristics,
        weight_gage_with_regression,
        add_text_section,
        edit_narrative,
        mark_narrative_reviewed,
        get_section,
        remove_section,
        move_section,
        render_report,
        save_report,
        load_report,
    )
    for tool in tools:
        server.add_tool(_reported(tool, lock), description=inspect.cleandoc(tool.__doc__ or ""))
    return server


def main(argv: Optional[List[str]] = None) -> None:
    """Run the server on standard input and output (the ``flowfreq-mcp`` command).

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments; defaults to ``sys.argv[1:]``.
    """
    parser = argparse.ArgumentParser(
        prog="flowfreq-mcp",
        description=(
            "flowfreq's report tools as a Model Context Protocol server, for a chat client "
            "to start. It speaks MCP on standard input and output."
        ),
    )
    parser.add_argument(
        "--dir",
        type=Path,
        default=None,
        help=(
            "the directory reports are written to and saved reports are read from "
            f"(default: ${REPORT_DIR_ENV}, or flowfreq-reports in the home directory)"
        ),
    )
    parser.add_argument("--version", action="version", version=f"flowfreq {__version__}")
    args = parser.parse_args(argv)
    # Standard output carries the protocol, so the log goes to standard error.
    logging.basicConfig(
        stream=sys.stderr, level=logging.INFO, format="%(levelname)s %(name)s: %(message)s"
    )
    build_server(output_dir=args.dir).run("stdio")


__all__ = ["INSTRUCTIONS", "REPORT_DIR_ENV", "build_server", "default_report_dir", "main"]


if __name__ == "__main__":
    main()
