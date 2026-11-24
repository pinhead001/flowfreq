"""The ``site_summary`` report section: what a study of one site found, on a page or two.

It sets out the basin characteristics of the site and of the gages used, a map of their
watersheds, the flood quantiles each method gives at the site, and the methods' frequency
curves on one plot. The estimates are the ones a ``comparison`` section takes, and are
computed the same way; the characteristics and watershed outlines are the analyst's inputs
(from StreamStats or a GIS), stored in the specification and printed as given. Discussion
goes in the section's ``narrative``.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..aep import AEP_MATCH
from .blocks import Block, Figure, Paragraph, Table
from .regional import _caution_paragraphs
from .sections import (
    SECTION_KINDS,
    SectionContent,
    SectionKind,
    _compared,
    _figure_png,
    _flow_table,
    _matched_aeps,
    _row_at,
    _version,
    cfs,
)
from .spec import ReportError

#: The AEPs summarized by default: the 1.5-year flood to the 500-year flood.
SUMMARY_AEPS: Tuple[float, ...] = (0.667, 0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)
#: Captions of the section's two tables.
SUMMARY_FLOWS_CAPTION = "Flood quantiles at the site by each method."
SUMMARY_CHARACTERISTICS_CAPTION = "Basin characteristics at each location."

_LOCATION_KEYS = ("name", "lat", "lon", "characteristics", "basin", "marker")
#: The point markers a location may take (matplotlib's).
_MARKERS = ("o", "s", "^", "v", "<", ">", "*", "D", "P", "X")
_OWNER = "site_summary"


def _rings(geometry: Any, where: str) -> List[np.ndarray]:
    """The outer ring of each polygon of a GeoJSON Polygon or MultiPolygon, as lon/lat."""
    if not isinstance(geometry, dict):
        raise ReportError(f"{where}: basin must be a GeoJSON Polygon or MultiPolygon geometry")
    kind, coordinates = geometry.get("type"), geometry.get("coordinates")
    polygons: Any
    if kind == "Polygon":
        polygons = [coordinates]
    elif kind == "MultiPolygon":
        polygons = coordinates
    else:
        raise ReportError(f"{where}: basin must be a GeoJSON Polygon or MultiPolygon geometry")
    try:
        rings = [np.array(polygon[0], dtype=np.float64)[:, :2] for polygon in polygons]
    except (TypeError, ValueError, IndexError) as exc:
        raise ReportError(f"{where}: basin coordinates are not rings of lon, lat pairs") from exc
    if not rings or any(ring.ndim != 2 or len(ring) < 4 for ring in rings):
        raise ReportError(f"{where}: basin coordinates are not rings of lon, lat pairs")
    return rings


def _locations(given: Any) -> List[Dict[str, Any]]:
    if not isinstance(given, list) or not given:
        raise ReportError(f"{_OWNER}: locations must be a non-empty list")
    out = []
    for index, each in enumerate(given):
        where = f"{_OWNER}: location {index + 1}"
        if not isinstance(each, dict):
            raise ReportError(f"{where} must be an object")
        extra = sorted(set(each) - set(_LOCATION_KEYS))
        if extra:
            raise ReportError(
                f"{where}: unknown key(s) {', '.join(extra)}; allowed: {', '.join(_LOCATION_KEYS)}"
            )
        name = str(each.get("name") or "").strip()
        if not name:
            raise ReportError(f"{where} needs a name")
        try:
            lat, lon = float(each["lat"]), float(each["lon"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ReportError(f"{where} needs lat and lon in decimal degrees") from exc
        if not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0):
            raise ReportError(f"{where}: lat and lon are not decimal degrees")
        characteristics = each.get("characteristics") or {}
        if not isinstance(characteristics, dict):
            raise ReportError(f"{where}: characteristics must be an object of code to value")
        try:
            values = {str(code): float(value) for code, value in characteristics.items()}
        except (TypeError, ValueError) as exc:
            raise ReportError(f"{where}: a characteristic is not a number") from exc
        if any(not math.isfinite(v) for v in values.values()):
            raise ReportError(f"{where}: a characteristic is not a finite number")
        marker = str(each.get("marker") or "o")
        if marker not in _MARKERS:
            raise ReportError(
                f"{where}: marker must be one of {' '.join(_MARKERS)}; got {marker!r}"
            )
        rings = None if each.get("basin") is None else _rings(each["basin"], where)
        out.append(
            {
                "name": name,
                "lat": lat,
                "lon": lon,
                "characteristics": values,
                "rings": rings,
                "marker": marker,
            }
        )
    if len({loc["name"] for loc in out}) != len(out):
        raise ReportError(f"{_OWNER}: two locations have the same name")
    return out


def _value(value: Optional[float]) -> str:
    if value is None:
        return "-"
    return f"{value:,.0f}" if abs(value) >= 100 else f"{value:.4g}"


def _map_figure(locations: Sequence[Dict[str, Any]]) -> bytes:
    """Watershed outlines, largest first, with each location marked and named."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    fig, ax = plt.subplots(figsize=(7, 8))
    outlined = [loc for loc in locations if loc["rings"]]
    order = sorted(outlined, key=lambda loc: -sum(len(r) for r in loc["rings"]))
    for index, loc in enumerate(order):
        color = f"C{locations.index(loc)}"
        for ring in loc["rings"]:
            ax.fill(ring[:, 0], ring[:, 1], color=color, alpha=0.10, zorder=1 + index)
            ax.plot(ring[:, 0], ring[:, 1], color=color, linewidth=1.2, zorder=2 + index)
    for index, loc in enumerate(locations):
        ax.plot(
            loc["lon"],
            loc["lat"],
            loc["marker"],
            color=f"C{index}",
            markeredgecolor="black",
            markersize=9,
            zorder=20,
            label=loc["name"] + (" (watershed outlined)" if loc["rings"] else ""),
        )
    middle = float(np.mean([loc["lat"] for loc in locations]))
    ax.set_aspect(1.0 / max(math.cos(math.radians(middle)), 0.1))
    ax.set_xlabel("Longitude (degrees)")
    ax.set_ylabel("Latitude (degrees)")
    ax.ticklabel_format(useOffset=False)
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.09), frameon=False)
    return _figure_png(fig)


def _curve_figure(
    curves: Sequence[Tuple[str, pd.DataFrame]], aeps: Sequence[float], where: str
) -> bytes:
    """Each method's flood quantiles against AEP on a normal probability axis, with limits."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    from scipy.special import ndtri

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for index, (label, flows) in enumerate(curves):
        rows = flows[flows["flow_cfs"].notna()].sort_values("aep", ascending=False)
        top, bottom = max(aeps) * (1 + AEP_MATCH), min(aeps) * (1 - AEP_MATCH)
        rows = rows[(rows["aep"] <= top) & (rows["aep"] >= bottom)]
        x = -ndtri(rows["aep"].to_numpy(dtype=float))
        color = f"C{index}"
        ax.plot(x, rows["flow_cfs"], "-o", color=color, markersize=4, label=label)
        limits = rows.dropna(subset=["lower", "upper"])
        if not limits.empty:
            xl = -ndtri(limits["aep"].to_numpy(dtype=float))
            ax.fill_between(xl, limits["lower"], limits["upper"], color=color, alpha=0.12)
            ax.plot(xl, limits["lower"], "--", color=color, linewidth=0.8)
            ax.plot(xl, limits["upper"], "--", color=color, linewidth=0.8)
    ticks = sorted(aeps, reverse=True)
    ax.set_xticks(-ndtri(np.array(ticks)))
    ax.set_xticklabels([f"{100 * a:.3g}" for a in ticks])
    ax.set_yscale("log")
    plain = FuncFormatter(lambda value, _: f"{value:,.0f}")
    ax.yaxis.set_major_formatter(plain)
    ax.yaxis.set_minor_formatter(plain)
    ax.tick_params(axis="y", which="minor", labelsize=8)
    ax.set_xlabel("Annual exceedance probability (percent)")
    ax.set_ylabel("Flow (cfs)")
    ax.set_title(f"Flood frequency at {where}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right")
    return _figure_png(fig)


def _build_site_summary(given: Dict[str, Any]) -> SectionContent:
    options = dict(given)
    locations = _locations(options["locations"])
    listed = options["estimates"]
    if not isinstance(listed, list) or not listed:
        raise ReportError(f"{_OWNER}: estimates must be a non-empty list")
    built = [_compared(i, each, _OWNER) for i, each in enumerate(listed)]
    labels = [label for label, _, _, _ in built]
    if len(set(labels)) != len(labels):
        raise ReportError(f"{_OWNER}: two estimates have the same label")
    tables = [_flow_table(flows) for _, _, flows, _ in built]
    aeps, apart = _matched_aeps(_OWNER, tables, options.get("aeps"), SUMMARY_AEPS)
    plotted = options.get("curves")
    if plotted is None:
        plotted = labels
    elif not isinstance(plotted, list) or not plotted or any(p not in labels for p in plotted):
        raise ReportError(
            f"{_OWNER}: curves must be a list of estimate labels; the labels are "
            + ", ".join(labels)
        )

    site = str(options["site_name"]).strip() if options.get("site_name") else ""
    where = site or "the site"
    blocks: List[Block] = []

    codes: List[str] = []
    for loc in locations:
        codes.extend(c for c in loc["characteristics"] if c not in codes)
    names = options.get("characteristic_names") or {}
    if not isinstance(names, dict):
        raise ReportError(f"{_OWNER}: characteristic_names must be an object of code to text")
    if codes:
        source = str(options.get("characteristics_source") or "").strip()
        blocks.append(
            Paragraph(
                "The basin characteristics below are inputs to this study, printed as "
                "supplied; flowfreq did not compute them."
                + (f" Source: {source}" if source else "")
            )
        )
        blocks.append(
            Table(
                ("Characteristic", "Code", *(loc["name"] for loc in locations)),
                tuple(
                    (
                        str(names.get(code, code)),
                        code,
                        *(_value(loc["characteristics"].get(code)) for loc in locations),
                    )
                    for code in codes
                ),
                ("left", "left", *("right",) * len(locations)),
                SUMMARY_CHARACTERISTICS_CAPTION,
            )
        )
    outlined = [loc["name"] for loc in locations if loc["rings"]]
    blocks.append(
        Figure(
            _map_figure(locations),
            "Locations used in this study"
            + (", with the watershed above " + ", ".join(outlined) if outlined else "")
            + ". Outlines are as supplied; there is no base map.",
        )
    )
    blocks.append(
        Paragraph(
            f"The flood quantiles at {where} by each method are summarized below; the "
            "sections that follow give each method's inputs, arithmetic and limits. The "
            "methods are compared, not combined. The computation was made with flowfreq "
            f"{_version()}."
        )
    )
    blocks.append(
        Table(
            (
                "Annual exceedance probability (percent)",
                "Recurrence interval (years)",
                *(f"{label} (cfs)" for label in labels),
            ),
            tuple(
                (
                    f"{aep * 100:.3g}",
                    f"{1.0 / aep:.3g}",
                    *(cfs(float(_row_at(table, aep)["flow_cfs"])) for table in tables),
                )
                for aep in aeps
            ),
            ("right",) * (len(labels) + 2),
            SUMMARY_FLOWS_CAPTION,
        )
    )
    if apart:
        blocks.append(
            Paragraph(
                "The methods tabulate some rows at slightly different probabilities, which "
                "are shown together: " + "; ".join(apart) + "."
            )
        )
    for label, content, _, _ in built:
        blocks.extend(_caution_paragraphs(label, content))
    curves = [(label, flows) for label, _, flows, _ in built if label in plotted]
    said = "; ".join(f"{label}, {limits}" for label, _, _, limits in built if label in plotted)
    blocks.append(
        Figure(
            _curve_figure(curves, aeps, where),
            f"Flood quantiles at {where} by method. Shaded bands are each method's limits "
            f"({said}).",
        )
    )
    references: List[str] = []
    resolved = []
    for (label, content, _, _), each in zip(built, listed):
        references.extend(r for r in content.references if r not in references)
        resolved.append({"label": label, "kind": each["kind"], "options": content.options})
    options["estimates"] = resolved
    return SectionContent(
        title=f"Summary: {site}" if site else "Summary",
        blocks=blocks,
        options=options,
        references=references,
        cautions=[f"{label}: {text}" for label, c, _, _ in built for text in c.cautions],
    )


SECTION_KINDS["site_summary"] = SectionKind(
    name="site_summary",
    summary=(
        "A summary of a study of one site: the basin characteristics of the site and the "
        "gages used (as supplied), a map of their watersheds, the flood quantiles each "
        "method gives at the site from the 1.5-year to the 500-year flood, and the methods' "
        "frequency curves on one plot."
    ),
    required=("locations", "estimates"),
    optional=(
        "site_name",
        "aeps",
        "curves",
        "characteristic_names",
        "characteristics_source",
    ),
    build=_build_site_summary,
)

__all__ = ["SUMMARY_AEPS", "SUMMARY_CHARACTERISTICS_CAPTION", "SUMMARY_FLOWS_CAPTION"]
