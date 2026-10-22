"""
Which gages are flow-connected: one upstream of the other on the same river.

Transfers between gages, and the pair fit in :mod:`flowfreq.transfer_variance`,
need to know which gages are nested. Two independent sources answer it:

- **Basin polygons.** A gage is upstream of another when its snapped outlet
  lies inside the other's watershed polygon and its area is not larger
  (:func:`flow_relation`). The polygons come from StreamStats
  (``WatershedCharacteristics.polygon_geojson``, via :func:`basin_from_watershed`).
- **NLDI navigation.** The Network-Linked Data Index lists every NWIS site
  upstream of a site on the NHDPlus network, tributaries included
  (:func:`nldi_upstream_sites`). Two sites on the same flowline are ordered by
  their NHDPlus measure, which runs from 0 at the flowline's downstream end to
  100 at its upstream end.

:func:`nested_gage_pairs` runs both where it can and reports each verdict, so a
pair the two sources disagree on is flagged rather than silently decided.

What this module cannot do
--------------------------
- **It cannot fix a bad snap.** A polygon test inherits StreamStats' snapping:
  an outlet snapped onto a side gully has a basin that contains nothing. NLDI
  indexes sites to flowlines by its own snapping, which can put a gage on a
  tributary next to the main stem. Disagreements between the two are the
  signal; resolving them is the analyst's job.
- **It does not judge regulation, diversions or losing reaches.** Nested is a
  statement about the network, not about whether the flows are comparable.

References
----------
U.S. Geological Survey, Network-Linked Data Index (NLDI) API,
https://api.water.usgs.gov/nldi/ (verified live 2026-10-21).
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from itertools import combinations
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import pandas as pd
import requests

from flowfreq.streamstats import (
    WatershedCharacteristics,
    _haversine_distance_m,
    _point_in_geometry,
    _request_with_backoff,
)

logger = logging.getLogger(__name__)

NLDI_BASE_URL = "https://api.water.usgs.gov/nldi/linked-data"

#: Two outlets this close, with areas within :data:`COLOCATED_AREA_TOLERANCE`,
#: are the same point on the river. Needed because a basin's own outlet sits
#: on its boundary, where a point-in-polygon test can fall either way.
COLOCATED_DISTANCE_M: float = 100.0

#: Fractional area agreement for two outlets to count as co-located.
COLOCATED_AREA_TOLERANCE: float = 0.05

#: Navigation distance, km: far enough to reach any headwater.
_NLDI_DISTANCE_KM = 9999

NESTED = "nested"
COLOCATED = "colocated"
NOT_NESTED = "not nested"


@dataclass(frozen=True)
class Basin:
    """One delineated basin: an outlet, an area and the polygon.

    Attributes
    ----------
    site_id : str
        Gage number or other label.
    area_sqmi : float
        Drainage area.
    outlet_lat, outlet_lon : float
        The snapped outlet -- the basin's real pour point, not the coordinate
        that was sent.
    geometry : dict
        GeoJSON Polygon or MultiPolygon geometry, WGS84.
    """

    site_id: str
    area_sqmi: float
    outlet_lat: float
    outlet_lon: float
    geometry: Dict[str, Any]


def basin_from_watershed(site_id: str, watershed: WatershedCharacteristics) -> Basin:
    """Build a :class:`Basin` from a StreamStats result.

    Raises
    ------
    ValueError
        The result carries no polygon or no positive ``DRNAREA``. A cache entry
        written without polygons cannot be tested for nesting.
    """
    polygon = watershed.polygon_geojson
    if polygon is None:
        raise ValueError(f"{site_id}: StreamStats result has no watershed polygon")
    geometry = polygon.get("geometry", polygon) if polygon.get("type") == "Feature" else polygon
    area = watershed.characteristics.get("DRNAREA")
    if area is None or area.value is None or not area.value > 0:
        raise ValueError(f"{site_id}: StreamStats result has no positive DRNAREA")
    return Basin(
        site_id=site_id,
        area_sqmi=float(area.value),
        outlet_lat=float(watershed.snap.snapped_lat),
        outlet_lon=float(watershed.snap.snapped_lon),
        geometry=geometry,
    )


def flow_relation(site: Basin, other: Basin) -> Optional[str]:
    """Where ``other`` sits on the river relative to ``site``, if at all.

    Returns
    -------
    str or None
        ``"colocated"`` when the outlets are within
        :data:`COLOCATED_DISTANCE_M` and the areas within
        :data:`COLOCATED_AREA_TOLERANCE`; ``"upstream"`` when ``other``'s
        outlet lies inside ``site``'s basin and its area is not larger;
        ``"downstream"`` the reverse; ``None`` on different streams.
        The area condition stops a point that lands inside a much smaller
        basin, by a ragged DEM boundary, from being called upstream of it.
    """
    distance = _haversine_distance_m(
        site.outlet_lat, site.outlet_lon, other.outlet_lat, other.outlet_lon
    )
    larger = max(site.area_sqmi, other.area_sqmi)
    if (
        distance <= COLOCATED_DISTANCE_M
        and abs(site.area_sqmi - other.area_sqmi) <= COLOCATED_AREA_TOLERANCE * larger
    ):
        return "colocated"
    if other.area_sqmi <= site.area_sqmi and _point_in_geometry(
        other.outlet_lon, other.outlet_lat, site.geometry
    ):
        return "upstream"
    if site.area_sqmi <= other.area_sqmi and _point_in_geometry(
        site.outlet_lon, site.outlet_lat, other.geometry
    ):
        return "downstream"
    return None


@dataclass(frozen=True)
class NLDISite:
    """An NWIS site as NLDI indexes it to the NHDPlus network."""

    site_no: str
    name: str
    comid: Optional[int]
    reachcode: str
    measure: Optional[float]


def _strip_agency(identifier: str) -> str:
    return identifier.split("-", 1)[1] if "-" in identifier else identifier


def nldi_upstream_sites(site_no: str, *, timeout: float = 60.0) -> List[NLDISite]:
    """Every NWIS site upstream of ``site_no``, tributaries included.

    The list includes ``site_no`` itself, and every NWIS site type NLDI
    indexes (lakes, diversions, wells on the network), not only stream gages.

    Parameters
    ----------
    site_no : str
        USGS site number, without the ``USGS-`` prefix.
    timeout : float
        Seconds per request.

    Returns
    -------
    list of NLDISite

    Raises
    ------
    requests.HTTPError
        NLDI does not know the site (404) or fails.
    """
    url = f"{NLDI_BASE_URL}/nwissite/USGS-{site_no}/navigation/UT/nwissite"
    response = _request_with_backoff(
        requests.get, url, params={"distance": _NLDI_DISTANCE_KM}, timeout=timeout
    )
    response.raise_for_status()
    out = []
    for feature in response.json().get("features", []):
        p = feature.get("properties", {})
        comid = p.get("comid")
        measure = p.get("measure")
        out.append(
            NLDISite(
                site_no=_strip_agency(str(p.get("identifier", ""))),
                name=str(p.get("name", "")),
                comid=None if comid is None else int(comid),
                reachcode=str(p.get("reachcode") or ""),
                measure=None if measure is None else float(measure),
            )
        )
    return out


Navigator = Callable[[str], Sequence[NLDISite]]


def _nldi_verdict(
    a: str, b: str, upstream_of: Mapping[str, Optional[Dict[str, NLDISite]]]
) -> Tuple[Optional[str], Optional[Tuple[str, str]]]:
    """NLDI's verdict on the pair, and ``(upstream, downstream)`` when nested."""
    ua, ub = upstream_of.get(a), upstream_of.get(b)
    for start, other, listing in ((a, b, ua), (b, a, ub)):
        if listing is None or other not in listing:
            continue
        here, there = listing.get(start), listing[other]
        same_line = here is not None and here.comid is not None and here.comid == there.comid
        if (
            same_line
            and here is not None
            and here.measure is not None
            and there.measure is not None
        ):
            if math.isclose(here.measure, there.measure, abs_tol=1e-6):
                return COLOCATED, (other, start)
            if there.measure < here.measure:
                return NESTED, (start, other)  # ``other`` is downstream on the flowline
        return NESTED, (other, start)
    if ua is not None and ub is not None:
        return NOT_NESTED, None
    return None, None


def _polygon_verdict(
    a: str, b: str, basins: Mapping[str, Basin]
) -> Tuple[Optional[str], Optional[Tuple[str, str]]]:
    if a not in basins or b not in basins:
        return None, None
    rel = flow_relation(basins[a], basins[b])
    if rel == "upstream":
        return NESTED, (b, a)
    if rel == "downstream":
        return NESTED, (a, b)
    if rel == "colocated":
        return COLOCATED, (a, b) if basins[a].area_sqmi <= basins[b].area_sqmi else (b, a)
    return NOT_NESTED, None


def nested_gage_pairs(
    site_nos: Sequence[str],
    *,
    basins: Optional[Mapping[str, Basin]] = None,
    navigator: Optional[Navigator] = nldi_upstream_sites,
) -> pd.DataFrame:
    """Every flow-connected pair among ``site_nos``, by polygons and by NLDI.

    Parameters
    ----------
    site_nos : sequence of str
        The gages to pair. Sites outside this list are ignored, even when
        NLDI lists them.
    basins : mapping of str to Basin, optional
        Polygons by site number. Pairs with both present get a polygon verdict.
    navigator : callable, optional
        Site number to its upstream sites, itself included. Defaults to
        :func:`nldi_upstream_sites`; ``None`` skips NLDI. A site whose
        navigation fails is logged and left without an NLDI verdict.

    Returns
    -------
    pd.DataFrame
        One row per pair that either source calls nested or co-located:
        ``upstream``, ``downstream``, ``polygon`` and ``nldi`` (each
        ``"nested"``, ``"colocated"``, ``"not nested"`` or ``None`` when that
        source could not say), ``relation`` (the polygon verdict where there
        is one, else NLDI's) and ``conflict`` (both sources answered and
        differ, in kind or in direction).

    Raises
    ------
    ValueError
        Neither ``basins`` nor ``navigator`` is given.
    """
    if basins is None and navigator is None:
        raise ValueError("nested_gage_pairs needs basins or a navigator")
    sites = list(dict.fromkeys(str(s) for s in site_nos))
    upstream_of: Dict[str, Optional[Dict[str, NLDISite]]] = {}
    if navigator is not None:
        for s in sites:
            try:
                upstream_of[s] = {x.site_no: x for x in navigator(s)}
            except Exception as exc:  # network, 404, malformed payload
                logger.warning("NLDI navigation failed for %s: %s", s, exc)
                upstream_of[s] = None

    rows = []
    for a, b in combinations(sites, 2):
        pv, pdir = _polygon_verdict(a, b, basins or {})
        nv, ndir = _nldi_verdict(a, b, upstream_of) if navigator is not None else (None, None)
        if NESTED not in (pv, nv) and COLOCATED not in (pv, nv):
            continue
        direction = pdir or ndir
        assert direction is not None
        conflict = (
            pv is not None and nv is not None and (pv != nv or (pv == NESTED and pdir != ndir))
        )
        rows.append(
            {
                "upstream": direction[0],
                "downstream": direction[1],
                "polygon": pv,
                "nldi": nv,
                "relation": pv if pv is not None and pv != NOT_NESTED else nv,
                "conflict": bool(conflict),
            }
        )
    columns = ["upstream", "downstream", "polygon", "nldi", "relation", "conflict"]
    frame = pd.DataFrame(rows, columns=columns)
    if frame["conflict"].any():
        logger.warning("polygons and NLDI disagree on %d pair(s)", int(frame["conflict"].sum()))
    return frame


def connections(
    pairs: pd.DataFrame,
    *,
    drop_conflicts: bool = False,
    include_colocated: bool = False,
) -> List[Tuple[str, str]]:
    """``(upstream, downstream)`` tuples for :func:`flowfreq.transfer_variance.build_pairs`.

    Parameters
    ----------
    pairs : pd.DataFrame
        From :func:`nested_gage_pairs`.
    drop_conflicts : bool, default False
        Leave out pairs the two sources disagree on.
    include_colocated : bool, default False
        Keep co-located pairs. They measure the same flow, so they say nothing
        about transfer error and are left out by default.
    """
    keep = pairs["relation"].isin([NESTED, COLOCATED] if include_colocated else [NESTED])
    if drop_conflicts:
        keep &= ~pairs["conflict"].astype(bool)
    return [
        (str(u), str(d)) for u, d in zip(pairs.loc[keep, "upstream"], pairs.loc[keep, "downstream"])
    ]
