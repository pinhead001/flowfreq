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
at the gage and weight_gage_with_regression (or transpose_gage_to_site for a site on \
the gage's stream), then get_section to read what was computed, \
edit_narrative, and render_report. check_gage_against_reference adds a comparison of a \
gage's computation with USGS peakfq or with numbers the user gives. low_flow_frequency and \
flow_duration add sections from a gage's daily flows, and apply_future_change_factors \
scales a section's flood quantiles by factors the user supplies. The report is kept for \
the life of the server; \
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


def _station(site_no: str) -> str:
    """A USGS station number, checked."""
    if not site_no.strip().isdigit():
        raise ReportError(f"site_no is a USGS station number, all digits; got {site_no!r}")
    return site_no.strip()


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

    def transpose_gage_to_site(
        gage_section_id: str,
        target_area_sqmi: float,
        exponent_citation: str,
        exponent: Optional[float] = None,
        exponent_aeps: Optional[List[float]] = None,
        exponents: Optional[List[float]] = None,
        exponent_region: Optional[str] = None,
        target_name: Optional[str] = None,
        donor_area_sqmi: Optional[float] = None,
        max_aep: Optional[float] = None,
        allow_out_of_range: bool = False,
        donor_site_class: Optional[str] = None,
        target_site_class: Optional[str] = None,
        allow_regulation_mismatch: bool = False,
        detail: bool = False,
        title: Optional[str] = None,
        position: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Transpose a gage's flood quantiles to an ungaged site on the same stream.

        The drainage-area ratio method: each of the gage's quantiles is multiplied by
        (site area / gage area) raised to the drainage-area exponent of a regional
        regression. The gage's run_ffa section must already be in the report, and stays.

        gage_section_id: the run_ffa section. target_area_sqmi: the site's drainage area.
        exponent_citation: the report and table the exponents were published in, or why a
        value was assumed. The exponents are the user's to supply from that report: ask
        for them, and do not give values from memory. Give either exponent (one value for
        every AEP) or exponent_aeps and exponents (the published values, AEP as a
        fraction, for example [0.5, 0.1, 0.01] and [0.70, 0.66, 0.62]); between published
        AEPs the exponent is interpolated, and beyond them held, which the table flags.
        exponent_region: the regression region, for the record. target_name: a name for
        the site. donor_area_sqmi: the gage's drainage area, when run_ffa did not find it.
        max_aep: leave out quantiles more frequent than this AEP (a fraction).
        allow_out_of_range: transpose although the area ratio is outside 0.5 to 1.5,
        where the regression itself is the better estimate.
        donor_site_class, target_site_class: reference, regulated, urban or unknown; a
        regulated site on one side only is refused unless allow_regulation_mismatch.
        detail: repeat the gage's analysis below the table. title, position: as run_ffa
        (default position: just after the gage's section).

        Returns the new section's id and the transposed table.
        """
        spec: Dict[str, Any] = {"citation": exponent_citation}
        if exponent is not None:
            spec["exponent"] = exponent
        if exponent_aeps is not None:
            spec["aeps"] = list(exponent_aeps)
        if exponents is not None:
            spec["exponents"] = list(exponents)
        if exponent_region:
            spec["region"] = exponent_region
        options: Dict[str, Any] = {}
        if donor_area_sqmi is not None:
            options["donor_area_sqmi"] = donor_area_sqmi
        if max_aep is not None:
            options["max_aep"] = max_aep
        if allow_out_of_range:
            options["allow_out_of_range"] = True
        if donor_site_class:
            options["donor_site_class"] = donor_site_class
        if target_site_class:
            options["target_site_class"] = target_site_class
        if allow_regulation_mismatch:
            options["allow_regulation_mismatch"] = True
        return report.transpose(
            gage_section_id,
            target_area_sqmi,
            spec,
            target_name=target_name,
            options=options,
            detail=detail,
            title=title,
            position=position,
        )

    def check_gage_against_reference(
        gage_section_id: str,
        reference_source: Optional[str] = None,
        reference_aeps: Optional[List[float]] = None,
        reference_flows_cfs: Optional[List[float]] = None,
        reference_lower_cfs: Optional[List[float]] = None,
        reference_upper_cfs: Optional[List[float]] = None,
        tolerance_pct: Optional[float] = None,
        ci_tolerance_pct: Optional[float] = None,
        detail: bool = False,
        title: Optional[str] = None,
        position: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Check a gage's flood frequency computation against a second computation of
        the same record, and add the comparison as a section. The gage's run_ffa section
        must already be in the report, and stays.

        With no reference_* argument the record is run through USGS peakfq 8.1.0, the
        reference implementation of Bulletin 17C, by flowfreq's Fortran extension. Where
        that extension is not built the tool refuses and says how to build it.

        Otherwise the reference is the user's own run of another program (peakfq,
        HEC-SSP): ask for its numbers, and do not give values from memory.
        reference_source: the program, its version and where the numbers were taken from.
        reference_aeps (fractions, for example [0.1, 0.01]) and reference_flows_cfs: its
        flood quantiles, in the same order. reference_lower_cfs, reference_upper_cfs: its
        90-percent confidence limits in that order, both or neither.

        tolerance_pct: the largest quantile difference that counts as agreement, in
        percent (default 1). ci_tolerance_pct: the same for the confidence limits
        (default 2). detail: repeat the gage's analysis below the comparison. title,
        position: as run_ffa (default position: just after the gage's section).

        Returns the new section's id, and whether the two agree with the quantile table.
        """
        columns = {
            "reference_aeps": reference_aeps,
            "reference_flows_cfs": reference_flows_cfs,
            "reference_lower_cfs": reference_lower_cfs,
            "reference_upper_cfs": reference_upper_cfs,
        }
        reference: Optional[Dict[str, Any]] = None
        if reference_source is not None or any(v is not None for v in columns.values()):
            if reference_aeps is None or reference_flows_cfs is None:
                raise ValueError(
                    "a reference needs reference_source, reference_aeps and "
                    "reference_flows_cfs; with none of the reference_* arguments the record "
                    "is run through peakfq"
                )
            if (reference_lower_cfs is None) != (reference_upper_cfs is None):
                raise ValueError("give reference_lower_cfs and reference_upper_cfs together")
            given = {name: list(v) for name, v in columns.items() if v is not None}
            if len({len(v) for v in given.values()}) != 1:
                raise ValueError(
                    f"{', '.join(given)} must have one value per AEP: {len(reference_aeps)}"
                )
            rows: List[Dict[str, Any]] = []
            for i, aep in enumerate(reference_aeps):
                row: Dict[str, Any] = {"aep": aep, "flow_cfs": reference_flows_cfs[i]}
                if reference_lower_cfs is not None and reference_upper_cfs is not None:
                    row["lower_cfs"] = reference_lower_cfs[i]
                    row["upper_cfs"] = reference_upper_cfs[i]
                rows.append(row)
            reference = {"source": reference_source, "quantiles": rows}
        return report.engine_check(
            gage_section_id,
            reference=reference,
            tolerance_pct=tolerance_pct,
            ci_tolerance_pct=ci_tolerance_pct,
            detail=detail,
            title=title,
            position=position,
        )

    def low_flow_frequency(
        site_no: str,
        n_day: int = 7,
        non_exceedance: Optional[List[float]] = None,
        year_type: str = "climatic",
        distribution: str = "lp3",
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        annual_table: bool = True,
        title: Optional[str] = None,
        position: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Add a low-flow frequency section for a USGS gage (7Q10, 7Q2 and the like).

        Downloads the gage's daily mean flows from USGS, takes the lowest n_day-day mean
        flow of each year, fits a log-Pearson Type III distribution to those minima and
        adds a section with the record, the method, the fitted parameters, the low-flow
        quantiles with 90-percent confidence limits, the fitted curve and the annual
        minima. The daily flows are stored in the report, so it rebuilds later without a
        download; a long record makes the saved report several hundred kilobytes.

        site_no: the USGS station number. n_day: the averaging period in days (7 for
        7Q10). non_exceedance: annual non-exceedance probabilities as fractions, for
        example [0.5, 0.1] for the 2-year and 10-year low flows; the default is 0.5, 0.2,
        0.1, 0.05, 0.02 and 0.01. year_type: "climatic" (April 1 to March 31, the
        default for low flows), "water" or "calendar". distribution: "lp3" (station
        skew) or "lognormal" (zero skew, steadier on a short record). start_date,
        end_date: limit the record, YYYY-MM-DD. annual_table: include the table of
        annual minima. title, position: as run_ffa.

        The tool refuses a record with fewer than 10 complete years. Returns the
        section's id and the quantile table.
        """
        options: Dict[str, Any] = {
            "site_no": _station(site_no),
            "n_day": n_day,
            "year_type": year_type,
            "distribution": distribution,
            "annual_table": annual_table,
        }
        if non_exceedance is not None:
            options["non_exceedance"] = list(non_exceedance)
        for name, value in (("start_date", start_date), ("end_date", end_date)):
            if value:
                options[name] = value
        return report.add_section("low_flow", options, title=title, position=position)

    def flow_duration(
        site_no: str,
        exceedance_pct: Optional[List[float]] = None,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        title: Optional[str] = None,
        position: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Add a flow-duration section for a USGS gage.

        Downloads the gage's daily mean flows from USGS and adds a section with the
        record, the flow equaled or exceeded on a given percent of days, and the
        flow-duration curve. Write these statistics with a lowercase q: q10 is the flow
        exceeded on 10 percent of days (a high flow) and q90 a low one. They are not
        flood quantiles and carry no annual probability. The daily flows are stored in
        the report, so it rebuilds later without a download.

        site_no: the USGS station number. exceedance_pct: percents of days, each
        between 0 and 100 (default 1, 5, 10, 20, 50, 80, 90, 95 and 99). start_date,
        end_date: limit the record, YYYY-MM-DD. title, position: as run_ffa.

        Returns the section's id and the table.
        """
        options: Dict[str, Any] = {"site_no": _station(site_no)}
        if exceedance_pct is not None:
            options["exceedance_pct"] = list(exceedance_pct)
        for name, value in (("start_date", start_date), ("end_date", end_date)):
            if value:
                options[name] = value
        return report.add_section("flow_duration", options, title=title, position=position)

    def apply_future_change_factors(
        section_id: str,
        factor_set_name: str,
        scenario: str,
        horizon: str,
        citation: str,
        factor_aeps: List[float],
        factors: List[float],
        legal_status: str = "informational",
        aeps: Optional[List[float]] = None,
        detail: bool = False,
        title: Optional[str] = None,
        position: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Scale a section's flood quantiles by change factors for a future condition,
        and add the result as a section. The section scaled is a run_ffa section or a
        regression estimate already in the report, and stays.

        flowfreq ships no change factors. They come from a source the user names (a
        state's climate guidance, an agency report): ask for the numbers and the
        citation, and do not give factors from memory.
        factor_set_name: a short name for the set. scenario: the emissions or warming
        scenario, for example "SSP5-8.5". horizon: the planning horizon, for example
        "2050". citation: the source document and table.
        factor_aeps (fractions, for example [0.1, 0.01]) and factors: the multipliers
        as published, in the same order. legal_status: "required", "recommended" or
        "informational".
        aeps: the AEPs to scale; by default those of the section that lie within
        factor_aeps. An AEP outside them takes the nearest published factor and is
        marked. detail: repeat the present-condition analysis below the table. title,
        position: as run_ffa (default position: just after the section scaled).

        Returns the new section's id and the table of present and future flows.
        """
        return report.future_flow(
            section_id,
            {
                "name": factor_set_name,
                "scenario": scenario,
                "horizon": horizon,
                "citation": citation,
                "aeps": list(factor_aeps),
                "factors": list(factors),
                "level": "site",
                "legal_status": legal_status,
            },
            aeps=aeps,
            detail=detail,
            title=title,
            position=position,
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
        transpose_gage_to_site,
        check_gage_against_reference,
        low_flow_frequency,
        flow_duration,
        apply_future_change_factors,
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
