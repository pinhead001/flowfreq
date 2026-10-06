"""Wave-level QA for the regression rollout: cross-border consistency and LOOCV.

``docs/MASTER_ROADMAP.md`` §3.2 asks for three things before a wave is done:

1. a **cross-border consistency check** -- evaluate each state's equations on
   border gages and on basins that straddle the line, and report the
   discontinuity at the border by AEP;
2. **transposition LOOCV** (§5.2) on the wave's gage network, as bias and RMSE in
   log space by AEP;
3. a wave vignette (``docs/vignettes/``).

This module is the offline half of 1 and 2: geometry on GeoJSON with no GIS
dependency, the mapping from NSS/StreamStats region codes to the offline
library's region codes, evaluation of one state's equations on a basin, the
discrepancy test, and the LOOCV summaries. ``tools/wave_qa.py`` does the network
half (StreamStats delineation, NSS region location, peaks, Bulletin 17C) with a
cache, and writes ``docs/WAVE_QA_W1_W2.md``.

**The discrepancy threshold.** Two equations evaluated on the same basin are two
estimates of one quantile. If each is unbiased with prediction standard error
``SEP_a`` and ``SEP_b`` (log10 units), and their errors are independent, their
log difference has standard deviation ``sqrt(SEP_a**2 + SEP_b**2)``. A pair is
flagged when ``|log10 Q_a - log10 Q_b| > z * sqrt(SEP_a**2 + SEP_b**2)`` with
``z = 1.645``, the two-sided 90 % level the library's own prediction intervals
use (:func:`flowfreq.regression.evaluate`). Independence is the conservative
reading: two regressions fitted partly on the same border gages have positively
correlated errors, which makes the true spread of their difference *smaller*, so
a flag under this test is a real inconsistency rather than noise. The same test
compares an equation with the gage's own Bulletin 17C estimate, using
``SEP_eq**2 + var(log Q_B17C)``, the latter from the B17C 90 % confidence limits.

Nothing here changes an equation. It reports.
"""

from __future__ import annotations

import logging
import math
import re
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from flowfreq.regression.equations import OutOfRangeError, evaluate
from flowfreq.regression.library import EquationsUnavailable, StateLibrary
from flowfreq.streamstats import _point_in_geometry, _polygons_of, geojson_area_sq_mi

logger = logging.getLogger(__name__)

#: Wave membership (``docs/MASTER_ROADMAP.md`` §3.2).
WAVE1: Tuple[str, ...] = ("WA", "OR", "ID", "MT")
WAVE2: Tuple[str, ...] = ("CO", "UT", "WY", "NM", "AZ", "NV")

#: Shared state lines checked: inside Wave 1, inside Wave 2, and between the two.
#: The Four Corners point pairs (CO-AZ, UT-NM) share a single point, not a line,
#: and are left out; borders with Wave 3+ states (CA, SD, NE, ...) wait for them.
BORDER_PAIRS: Tuple[Tuple[str, str], ...] = (
    ("WA", "OR"),
    ("WA", "ID"),
    ("OR", "ID"),
    ("ID", "MT"),
    ("CO", "UT"),
    ("CO", "WY"),
    ("CO", "NM"),
    ("UT", "WY"),
    ("UT", "NV"),
    ("UT", "AZ"),
    ("NM", "AZ"),
    ("AZ", "NV"),
    ("ID", "UT"),
    ("ID", "NV"),
    ("ID", "WY"),
    ("MT", "WY"),
    ("OR", "NV"),
)

#: Two-sided 90 % standard-normal critical value.
Z90: float = 1.6448536269514722

#: Eastern Oregon's six flood regions (Cooper, 2006, OWRD OFR SW 06-001) are not in
#: NSS. StreamStats' ``nss/regions`` MapServer layer 31 (``grid_name='regeastor'``)
#: carries them as gc1754-gc1759, named ``Reg1_East_Slope_Cascade_2006_OROFR_06001``
#: through ``Reg6_Southeast_Eastern_2006_OROFR_06001`` (read live 2026-10-04).
OR_EASTERN_GRIDCODES: Dict[str, str] = {
    "GC1754": "E1",
    "GC1755": "E2",
    "GC1756": "E3",
    "GC1757": "E4",
    "GC1758": "E5",
    "GC1759": "E6",
}

#: StreamStats ``nss/regions`` MapServer peak-flow region layers: state -> (layer id,
#: ``grid_name`` values), read live 2026-10-06. NSS ``bylocation`` cannot stand in for
#: these everywhere: New Mexico's nine regions share one state-wide polygon
#: (``pkdummy``, every region at 100 %), and the selector is the ``HIGHREG``
#: characteristic StreamStats computes on an NM delineation. Wyoming and Nevada have
#: no layer, so they are located with NSS ``bylocation``.
REGION_LAYERS: Dict[str, Tuple[int, Tuple[str, ...]]] = {
    "WA": (40, ("nss_peak",)),
    "OR": (31, ("or_regions", "regeastor")),
    "ID": (11, ("reg6n3",)),
    "MT": (22, ("pkregionid",)),
    "CO": (5, ("newpkregid",)),
    "UT": (37, ("ut_pkregs",)),
    "AZ": (3, ("az_regions_c",)),
}

#: Overlay regions chosen by mean basin elevation rather than by polygon: Arizona's
#: region 1 (High Elevation, gc1618) is drawn over every other region's polygon, and
#: NSS takes it when ``ELEV >= 7,500 ft`` (AZ.json notes; SIR 2014-5211).
ELEVATION_OVERLAYS: Dict[str, Tuple[str, float]] = {"AZ": ("GC1618", 7500.0)}

#: Stand-ins for a characteristic a StreamStats region does not compute, in order of
#: preference: the same physical quantity, from another source vintage. Drainage area
#: and contributing drainage area are equal on the basins checked here (no closed
#: basins at gages that pass the screen); PRISM mean annual precipitation normals for
#: 1961-90 (``PRECIP`` in UT/NM/NV), 1971-2000 (``PRECIP`` in CO/MT/WY), 1981-2010
#: (``PRECPRIS10``) and 1991-2020 (``PRECPRIS20``) differ by a few percent; NLCD
#: forest cover for 2001/2011/2016/2021 likewise. Nothing else is substituted.
PROXIES: Dict[str, Tuple[str, ...]] = {
    "CONTDA": ("DRNAREA",),
    "DRNAREA": ("CONTDA",),
    "PRECPRIS10": ("PRECPRIS20", "PRECIP"),
    "PRECIP": ("PRECPRIS10", "PRECPRIS20"),
    "LC11FOREST": ("LC16FOREST", "LC21FOREST", "FOREST"),
    "FOREST": ("LC11FOREST", "LC16FOREST", "LC21FOREST"),
}

#: Region codes left out of a state-vs-state comparison, with the reason. The
#: Navajo Nation equations (Waltemeyer, SIR 2006-5306) apply on Navajo lands in
#: four states and are not any one state's border equations; Montana's
#: channel-width equations need a field-measured width.
_EXCLUDED_REGIONS: Tuple[Tuple[str, str], ...] = (
    (r"^Navajo", "Navajo Nation equations (SIR 2006-5306), not a state equation"),
    (r"-(AC|BF|RS)$", "channel-width equation (needs a measured width)"),
)

_GC_PATTERN = re.compile(r"\bGC\d+\b")


# ----------------------------------------------------------------------------
# Geometry (GeoJSON, lon/lat, no GIS dependency)
# ----------------------------------------------------------------------------


def _geometry(obj: Mapping[str, Any]) -> Dict[str, Any]:
    if obj.get("type") == "Feature":
        geom = obj.get("geometry")
        if not isinstance(geom, dict):
            raise ValueError("Feature has no geometry")
        return geom
    return dict(obj)


def point_in(lat: float, lon: float, geometry: Mapping[str, Any]) -> bool:
    """Whether a point lies inside a GeoJSON Polygon/MultiPolygon (or a Feature of one)."""
    return _point_in_geometry(lon, lat, _geometry(geometry))


def distance_km(lat: float, lon: float, geometry: Mapping[str, Any]) -> float:
    """Distance from a point to a polygon, km: 0 inside, else to the nearest edge.

    Parameters
    ----------
    lat, lon : float
        WGS84 decimal degrees.
    geometry : dict
        GeoJSON Polygon/MultiPolygon, or a Feature carrying one.

    Returns
    -------
    float
    """
    geom = _geometry(geometry)
    if bool(contains(geom, np.array([lat]), np.array([lon]))[0]):
        return 0.0
    # Vectorised form of flowfreq.streamstats._distance_to_boundary_m (local
    # equirectangular about the point), fast on state- and region-sized rings.
    r = 6378137.0
    kx = math.cos(math.radians(lat)) * math.pi * r / 180.0
    ky = math.pi * r / 180.0
    best = math.inf
    for polygon in _polygons_of(geom):
        for ring in polygon:
            a = np.asarray([p[:2] for p in ring], dtype=float)
            if len(a) < 2:
                continue
            x1, y1 = (a[:-1, 0] - lon) * kx, (a[:-1, 1] - lat) * ky
            dx, dy = (a[1:, 0] - a[:-1, 0]) * kx, (a[1:, 1] - a[:-1, 1]) * ky
            seg2 = dx * dx + dy * dy
            with np.errstate(invalid="ignore", divide="ignore"):
                t = np.where(seg2 > 0, -(x1 * dx + y1 * dy) / seg2, 0.0)
            t = np.clip(t, 0.0, 1.0)
            best = min(best, float(np.min(np.hypot(x1 + t * dx, y1 + t * dy))))
    return best / 1000.0


def circle(lat: float, lon: float, radius_km: float, n: int = 64) -> Dict[str, Any]:
    """A GeoJSON Polygon approximating a circle (local equirectangular).

    Raises
    ------
    ValueError
        For a non-positive radius or fewer than 8 vertices.
    """
    if not radius_km > 0:
        raise ValueError(f"radius_km must be positive, got {radius_km}")
    if n < 8:
        raise ValueError("a circle needs at least 8 vertices")
    dlat = radius_km / 111.32
    dlon = dlat / max(math.cos(math.radians(lat)), 1e-6)
    ring = [
        [lon + dlon * math.cos(2 * math.pi * k / n), lat + dlat * math.sin(2 * math.pi * k / n)]
        for k in range(n)
    ]
    ring.append(ring[0])
    return {"type": "Polygon", "coordinates": [ring]}


def sample_points(geometry: Mapping[str, Any], n: int = 400) -> List[Tuple[float, float]]:
    """About ``n`` points on a regular lon/lat grid inside a polygon, as (lat, lon).

    A grid sample of a basin is how the fraction of it inside another polygon (a
    state, a regression region) is estimated without a polygon-intersection
    library. With 400 points the fraction is good to about 1-3 % of the basin.

    Raises
    ------
    ValueError
        If ``n`` is less than 10.
    """
    if n < 10:
        raise ValueError("n must be at least 10")
    geom = _geometry(geometry)
    xs: List[float] = []
    ys: List[float] = []
    for polygon in _polygons_of(geom):
        if polygon:
            xs.extend(p[0] for p in polygon[0])
            ys.extend(p[1] for p in polygon[0])
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    kx = math.cos(math.radians((y0 + y1) / 2))
    # Grid spacing from the bounding box; refine until enough points fall inside.
    out: List[Tuple[float, float]] = []
    m = max(int(math.sqrt(n)), 4)
    for _ in range(8):
        step = max((x1 - x0) * kx, (y1 - y0)) / m
        sx, sy = step / max(kx, 1e-6), step
        gx, gy = np.meshgrid(np.arange(x0 + sx / 2, x1, sx), np.arange(y0 + sy / 2, y1, sy))
        lon_g, lat_g = gx.ravel(), gy.ravel()
        inside = contains(geom, lat_g, lon_g)
        out = [(float(a), float(b)) for a, b in zip(lat_g[inside], lon_g[inside])]
        if len(out) >= n * 0.6:
            break
        m = int(m * 1.6) + 1
    if not out:
        # A sliver smaller than one cell: fall back to its vertices' centroid.
        out = [(float(np.mean(ys)), float(np.mean(xs)))]
    return out


def contains(geometry: Mapping[str, Any], lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Vectorised point-in-polygon: a boolean array, holes excluded.

    Uses :class:`matplotlib.path.Path` (a flowfreq dependency already), which is
    C-speed on the 10 m DEM basin polygons StreamStats returns (tens of thousands
    of vertices); the scalar ray cast in :mod:`flowfreq.streamstats` is not.
    """
    from matplotlib.path import Path as MplPath

    geom = _geometry(geometry)
    pts = np.column_stack([np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)])
    result = np.zeros(len(pts), dtype=bool)
    for polygon in _polygons_of(geom):
        if not polygon or len(polygon[0]) < 3:
            continue
        ring0 = np.asarray([p[:2] for p in polygon[0]], dtype=float)
        hit = MplPath(ring0).contains_points(pts)
        for hole in polygon[1:]:
            if len(hole) >= 3:
                hit &= ~MplPath(np.asarray([p[:2] for p in hole], dtype=float)).contains_points(pts)
        result |= hit
    return result


def fraction_inside(points: Sequence[Tuple[float, float]], geometry: Mapping[str, Any]) -> float:
    """Share of sample points ``(lat, lon)`` inside a polygon, 0-1."""
    if not points:
        raise ValueError("no sample points")
    arr = np.asarray(points, dtype=float)
    return float(contains(geometry, arr[:, 0], arr[:, 1]).mean())


def area_sq_mi(geometry: Mapping[str, Any]) -> float:
    """Geodesic area of a GeoJSON polygon, square miles."""
    return geojson_area_sq_mi(_geometry(geometry))


def esri_rings_to_geojson(rings: Sequence[Sequence[Sequence[float]]]) -> Dict[str, Any]:
    """An ArcGIS ``rings`` polygon as a GeoJSON MultiPolygon, holes kept.

    ArcGIS orders exterior rings clockwise and holes counter-clockwise (in lon/lat,
    a clockwise ring has negative shoelace area). Each hole is attached to the
    first exterior that contains its first vertex.
    """
    exteriors: List[List[List[List[float]]]] = []
    holes: List[List[List[float]]] = []
    for ring in rings:
        pts = [[float(p[0]), float(p[1])] for p in ring]
        if len(pts) < 4:
            continue
        signed = sum(
            (x2 - x1) * (y2 + y1) for (x1, y1), (x2, y2) in zip(pts[:-1], pts[1:])
        )  # > 0 is clockwise
        if signed > 0:
            exteriors.append([pts])
        else:
            holes.append(pts)
    for hole in holes:
        x, y = hole[0]
        for poly in exteriors:
            if _point_in_geometry(x, y, {"type": "Polygon", "coordinates": [poly[0]]}):
                poly.append(hole)
                break
    return {"type": "MultiPolygon", "coordinates": exteriors}


# ----------------------------------------------------------------------------
# Region codes: NSS/StreamStats -> offline library
# ----------------------------------------------------------------------------


def excluded_reason(region_code: str) -> Optional[str]:
    """Why a library region is left out of the border comparison, or None."""
    for pattern, reason in _EXCLUDED_REGIONS:
        if re.search(pattern, region_code):
            return reason
    return None


def nss_region_map(lib: StateLibrary) -> Dict[str, str]:
    """Map NSS/StreamStats region codes (``GC1750``) to the library's region codes.

    The library records each region's NSS code in its ``region_name``
    ("Regression Region 1 (NSS GC1750, ...)"). Oregon's eastern regions, which
    NSS does not carry, are added from :data:`OR_EASTERN_GRIDCODES`. Excluded
    regions (:func:`excluded_reason`) are left out.

    Returns
    -------
    dict of str to str
    """
    out: Dict[str, str] = {}
    for eq in lib.equations:
        if excluded_reason(eq.region_code):
            continue
        for code in _GC_PATTERN.findall(eq.region_name):
            out[code.upper()] = eq.region_code
    if lib.state == "OR":
        out.update(OR_EASTERN_GRIDCODES)
    return out


def resolve_overlay(
    state: str, codes: Sequence[str], characteristics: Mapping[str, float]
) -> List[str]:
    """The region codes that apply in a polygon that names several (an overlay).

    Arizona's polygons each carry their own region and region 1 (``gc1619,gc1618``);
    region 1 applies when mean basin elevation ``ELEV`` is at least 7,500 ft, the
    other otherwise (:data:`ELEVATION_OVERLAYS`). Oregon's western-interior polygon
    lists Regions 1, 2A and 2B; it is Region 2 (2A and 2B, blended on elevation by
    :func:`estimate_state`). Other states' codes pass through.

    Raises
    ------
    KeyError
        An overlay polygon and no ``ELEV`` to decide it.
    """
    codes = [c.strip().upper() for c in codes if c.strip()]
    if state == "OR" and ("GC730" in codes or "GC731" in codes):
        # Layer 31's western-interior polygon is 'gc10001,gc730,gc731,gc729': the
        # Region 2A/2B polygon, with Region 1 listed as well. Inside it Region 2
        # (the elevation blend) applies.
        return ["GC730", "GC731"]
    overlay = ELEVATION_OVERLAYS.get(state)
    if overlay is None or overlay[0] not in codes or len(codes) == 1:
        return codes
    code, threshold = overlay
    if "ELEV" not in characteristics:
        raise KeyError(f"{state}: ELEV needed to choose region {code} or {codes}")
    if float(characteristics["ELEV"]) >= threshold:
        return [code]
    return [c for c in codes if c != code]


def highreg_region(value: float) -> str:
    """New Mexico's ``HIGHREG`` characteristic as an NSS code (1096 -> ``GC1096``).

    Raises
    ------
    ValueError
        For a value that is not a positive whole number.
    """
    if not (value > 0 and float(value).is_integer()):
        raise ValueError(f"HIGHREG must be a positive integer GC number, got {value}")
    return f"GC{int(value)}"


def region_weights(
    lib: StateLibrary, located: Iterable[Tuple[str, float]]
) -> Tuple[Dict[str, float], float]:
    """Library region weights from located NSS regions.

    Parameters
    ----------
    lib : StateLibrary
    located : iterable of (nss_code, percent)
        NSS ``bylocation`` answers (or grid-sample fractions as percents). Codes
        that are not one of the state's peak-flow regions (low-flow regions, urban
        and national studies) are ignored.

    Returns
    -------
    (weights, coverage)
        ``weights`` maps library region codes to fractions summing to 1;
        ``coverage`` is the share (0-1) of the area the regions covered before
        normalising. Oregon's Regions 2A and 2B share one polygon and are an
        elevation split, so they are returned together as region ``"2"``.
        Empty weights when no peak-flow region was located.
    """
    mapping = nss_region_map(lib)
    raw: Dict[str, float] = {}
    for code, pct in located:
        region = mapping.get(str(code).upper())
        if region is None:
            continue
        if lib.state == "OR" and region in ("2A", "2B"):
            region = "2"
        # 2A and 2B both answer for the same polygon: count it once.
        raw[region] = max(raw.get(region, 0.0), float(pct))
    total = sum(raw.values())
    if total <= 0:
        return {}, 0.0
    return {r: w / total for r, w in raw.items()}, min(total / 100.0, 1.0)


# ----------------------------------------------------------------------------
# Characteristics and one state's estimate
# ----------------------------------------------------------------------------


def prepare_characteristics(
    characteristics: Mapping[str, float], outlet_lat: float, outlet_lon: float
) -> Tuple[Dict[str, float], Tuple[str, ...]]:
    """Fill the characteristics an equation can take from the outlet, and note proxies.

    ``LAT_OUT``/``LAT_GAGE`` (decimal degrees) and ``LONG_OUT`` (positive degrees
    west, as Wyoming's equations define it) are exact from the outlet coordinate.
    ``CONTDA`` and ``DRNAREA`` stand in for each other when a region computes only
    one: equal for a basin with no noncontributing area, which is the case at
    every gage this check uses (Arizona and Montana calibrate on contributing
    area). The other :data:`PROXIES` are the same quantity from a different data
    vintage (PRISM normals, NLCD year), which is how a neighbour's equation gets
    evaluated on a basin delineated in a StreamStats region that does not compute
    its exact code. Each proxy is listed, and :func:`estimate_state` reports the
    ones an estimate actually used, so a flagged pair can be traced to them.

    Returns
    -------
    (values, substitutions)
        ``substitutions`` entries read ``"TARGET=SOURCE"``.
    """
    out = {str(k): float(v) for k, v in characteristics.items()}
    subs: List[str] = []
    for code in ("LAT_OUT", "LAT_GAGE"):
        out.setdefault(code, float(outlet_lat))
    out.setdefault("LONG_OUT", abs(float(outlet_lon)))
    have = dict(out)
    for target, sources in PROXIES.items():
        if target in have:
            continue
        src = next((s for s in sources if s in have), None)
        if src is not None:
            out[target] = have[src]
            subs.append(f"{target}={src}")
    return out, tuple(subs)


@dataclass(frozen=True)
class StateEstimate:
    """One state's equations evaluated on one basin at one AEP.

    ``flow_cfs`` is the area-weighted mean of the located regions' flows (linear
    weighting, as NSS's own ``areaave`` does). ``sep_log`` is the same weighting of
    the regions' SEPs, ``None`` if any contributing region publishes none.
    """

    state: str
    aep: float
    flow_cfs: float
    sep_log: Optional[float]
    regions: Dict[str, float]
    out_of_range: Tuple[str, ...] = ()
    substitutions: Tuple[str, ...] = ()

    @property
    def log_flow(self) -> float:
        """log10 of :attr:`flow_cfs`."""
        return math.log10(self.flow_cfs)


def estimate_state(
    lib: StateLibrary,
    weights: Mapping[str, float],
    characteristics: Mapping[str, float],
    aep: float,
    *,
    substitutions: Sequence[str] = (),
) -> StateEstimate:
    """Evaluate one state's equations on a basin, area-weighting across regions.

    Inputs outside a region's calibrated range are evaluated anyway (this is a
    QA comparison, and the range violation is the finding) and listed in
    ``out_of_range``. Oregon's region ``"2"`` is the Region 2A/2B elevation blend
    (:func:`flowfreq.regression.oregon.estimate_region2`).

    Raises
    ------
    EquationsUnavailable
        No weights, or a region publishes no equation at this AEP.
    KeyError
        A characteristic an equation needs is missing.
    ValueError
        A transform's argument is outside its domain.
    """
    if not weights:
        raise EquationsUnavailable(f"{lib.state}: no peak-flow region located")
    flows: List[Tuple[float, float]] = []
    seps: List[Tuple[float, Optional[float]]] = []
    oor: List[str] = []
    codes_used: set = set()
    for region, w in weights.items():
        if lib.state == "OR" and region == "2":
            from flowfreq.regression.oregon import estimate_region2

            r2 = estimate_region2(aep, characteristics, lib=lib, allow_extrapolation=True)
            flows.append((w, r2.flow_cfs))
            comp = [(r2.weights[k], e.equation.sep_log) for k, e in r2.estimates.items()]
            sep = (
                None
                if any(s is None for _, s in comp)
                else sum(cw * float(s) for cw, s in comp if s is not None)
            )
            seps.append((w, sep))
            for k, e in r2.estimates.items():
                oor.extend(f"{k}: {p}" for p in e.out_of_range)
                codes_used.update(e.inputs)
            continue
        est = evaluate(lib.equation(region, aep), characteristics, allow_extrapolation=True)
        flows.append((w, est.flow_cfs))
        seps.append((w, est.equation.sep_log))
        oor.extend(f"{region}: {p}" for p in est.out_of_range)
        codes_used.update(est.inputs)
    flow = sum(w * q for w, q in flows)
    sep_log = (
        None
        if any(s is None for _, s in seps)
        else sum(w * float(s) for w, s in seps if s is not None)
    )
    return StateEstimate(
        state=lib.state,
        aep=aep,
        flow_cfs=flow,
        sep_log=sep_log,
        regions=dict(weights),
        out_of_range=tuple(oor),
        substitutions=tuple(s for s in substitutions if s.split("=")[0] in codes_used),
    )


def try_estimate_state(
    lib: StateLibrary,
    weights: Mapping[str, float],
    characteristics: Mapping[str, float],
    aep: float,
    *,
    substitutions: Sequence[str] = (),
) -> Tuple[Optional[StateEstimate], str]:
    """:func:`estimate_state`, returning ``(None, reason)`` instead of raising."""
    try:
        return (
            estimate_state(lib, weights, characteristics, aep, substitutions=substitutions),
            "",
        )
    except EquationsUnavailable as exc:
        return None, f"no equation: {exc}"
    except KeyError as exc:
        return None, f"missing characteristic: {exc.args[0] if exc.args else exc}"
    except (OutOfRangeError, ValueError) as exc:
        return None, f"not evaluable: {exc}"


# ----------------------------------------------------------------------------
# The discrepancy test
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class Discrepancy:
    """``log10 a - log10 b`` against its threshold ``z * sqrt(var_a + var_b)``."""

    diff_log: float
    threshold_log: Optional[float]
    note: str = ""

    @property
    def flagged(self) -> bool:
        """Whether the difference exceeds the threshold (False when there is none)."""
        return self.threshold_log is not None and abs(self.diff_log) > self.threshold_log

    @property
    def ratio(self) -> float:
        """``a / b``."""
        return float(10.0**self.diff_log)


def discrepancy(
    log_a: float,
    sd_a: Optional[float],
    log_b: float,
    sd_b: Optional[float],
    *,
    z: float = Z90,
) -> Discrepancy:
    """Compare two log10 estimates against ``z`` times their combined standard error.

    When one standard error is unknown (Nevada's hybrid regions publish none) the
    other alone sets the threshold, which makes the test stricter, and the note
    says so. With neither, no threshold is set and nothing is flagged.

    Raises
    ------
    ValueError
        For a non-positive ``z`` or a negative standard error.
    """
    if not z > 0:
        raise ValueError(f"z must be positive, got {z}")
    for sd in (sd_a, sd_b):
        if sd is not None and sd < 0:
            raise ValueError(f"standard errors must be non-negative, got {sd}")
    known = [s for s in (sd_a, sd_b) if s is not None]
    note = "" if len(known) == 2 else ("one SE unknown" if known else "no SE")
    threshold = z * math.sqrt(sum(s * s for s in known)) if known else None
    return Discrepancy(diff_log=log_a - log_b, threshold_log=threshold, note=note)


#: A B17C fit is not used as the gage's "observed" quantile when its 2-year flow is
#: below 1 cfs or its 1 %-to-2-year ratio exceeds 1,000. Both happen only at short,
#: mostly zero-flow records (20-30 years, most of them zero or censored by MGBT),
#: where the fitted LP3 is extrapolated far past any data: one 21-year Columbia
#: Plateau record gives Q2 = 1e-11 cfs and Q1% = 1e15 cfs.
DEGENERATE_Q2_CFS: float = 1.0
DEGENERATE_RATIO: float = 1000.0


def b17c_degenerate(q2: float, q100: float) -> bool:
    """Whether a B17C fit is unusable as a reference (see :data:`DEGENERATE_RATIO`).

    Non-finite or non-positive quantiles count as degenerate.
    """
    if not (math.isfinite(q2) and math.isfinite(q100)) or q2 <= 0 or q100 <= 0:
        return True
    return q2 < DEGENERATE_Q2_CFS or q100 / q2 > DEGENERATE_RATIO


def b17c_log_sd(lower: float, upper: float, *, z: float = Z90) -> float:
    """Standard error of a B17C log quantile from its two-sided confidence limits.

    ``(log10 upper - log10 lower) / (2 z)``; with the default ``z`` the limits are
    the 90 % (5 %/95 %) limits :func:`flowfreq.workflow.run_ffa` reports.

    Raises
    ------
    ValueError
        For non-positive or inverted limits.
    """
    if not (0 < lower <= upper):
        raise ValueError(f"need 0 < lower <= upper, got {lower}, {upper}")
    return (math.log10(upper) - math.log10(lower)) / (2.0 * z)


def compare_basin(
    libs: Mapping[str, StateLibrary],
    home: str,
    neighbor: str,
    characteristics: Mapping[str, float],
    outlet: Tuple[float, float],
    located_home: Iterable[Tuple[str, float]],
    located_neighbor: Iterable[Tuple[str, float]],
    b17c: Optional[Mapping[float, Tuple[float, float, float]]],
    aeps: Sequence[float],
    *,
    min_b17c_cfs: float = 1.0,
) -> List[Dict[str, Any]]:
    """Both states' equations on one basin, against each other and the gage, by AEP.

    Parameters
    ----------
    libs : mapping of state to StateLibrary
    home, neighbor : str
        The gage's state and the state across the line.
    characteristics : mapping
        Basin characteristics from the home state's StreamStats region.
    outlet : (lat, lon)
        The snapped outlet, for the coordinate characteristics.
    located_home, located_neighbor : iterable of (nss_code, percent)
        Each state's regions for the basin (:func:`region_weights` input).
    b17c : mapping of AEP to (quantile, lower 90 %, upper 90 %), optional
        The gage's Bulletin 17C estimate.
    aeps : sequence of float
    min_b17c_cfs : float, default 1.0
        B17C quantiles below this are not compared. At an ephemeral gage whose record
        is mostly zero-flow years the frequent quantiles collapse toward zero (5e-17
        cfs at one Columbia Plateau tributary), and a log comparison with them means
        nothing.

    Returns
    -------
    list of dict
        One row per AEP: ``home_q``/``neighbor_q``, their ``_sep`` and
        ``_oor`` (inputs outside the calibrated range), ``_why`` when a side is
        not evaluable, ``diff_log``/``thr_log``/``flag_hn`` for home vs.
        neighbour, and ``b17c_q``, ``<side>_vs_b17c_log``/``_thr`` and
        ``flag_<side>_b17c`` for each side against the gage.
    """
    chars, subs = prepare_characteristics(characteristics, outlet[0], outlet[1])
    wh, cov_h = region_weights(libs[home], located_home)
    wn, cov_n = region_weights(libs[neighbor], located_neighbor)
    weights = {"home": wh, "neighbor": wn}
    states = {"home": home, "neighbor": neighbor}
    rows: List[Dict[str, Any]] = []
    for aep in aeps:
        row: Dict[str, Any] = {
            "aep": aep,
            "home_regions": ";".join(f"{r}:{v:.2f}" for r, v in weights["home"].items()),
            "neighbor_regions": ";".join(f"{r}:{v:.2f}" for r, v in weights["neighbor"].items()),
            "home_coverage": round(cov_h, 3),
            "neighbor_coverage": round(cov_n, 3),
        }
        ests: Dict[str, Optional[StateEstimate]] = {}
        for role in ("home", "neighbor"):
            est, why = try_estimate_state(
                libs[states[role]], weights[role], chars, aep, substitutions=subs
            )
            ests[role] = est
            row[f"{role}_q"] = None if est is None else est.flow_cfs
            row[f"{role}_sep"] = None if est is None else est.sep_log
            row[f"{role}_oor"] = "" if est is None else "; ".join(est.out_of_range)
            row[f"{role}_proxies"] = "" if est is None else ";".join(est.substitutions)
            row[f"{role}_why"] = why
        h, n = ests["home"], ests["neighbor"]
        if h is not None and n is not None:
            d = discrepancy(h.log_flow, h.sep_log, n.log_flow, n.sep_log)
            row.update(
                diff_log=d.diff_log, thr_log=d.threshold_log, flag_hn=d.flagged, note_hn=d.note
            )
        fit = (b17c or {}).get(aep)
        if fit is not None and fit[0] >= min_b17c_cfs:
            q, lo, hi = fit
            sd = b17c_log_sd(lo, hi) if (0 < lo <= hi) else None
            row.update(b17c_q=q, b17c_sd=sd)
            for role, e in ests.items():
                if e is None:
                    continue
                d = discrepancy(e.log_flow, e.sep_log, math.log10(q), sd)
                row[f"{role}_vs_b17c_log"] = d.diff_log
                row[f"{role}_vs_b17c_thr"] = d.threshold_log
                row[f"flag_{role}_b17c"] = d.flagged
        rows.append(row)
    return rows


# ----------------------------------------------------------------------------
# LOOCV summaries
# ----------------------------------------------------------------------------


def summarize_log_errors(
    frame: pd.DataFrame,
    by: Sequence[str],
    *,
    estimate: str = "log_est",
    observed: str = "log_obs",
) -> pd.DataFrame:
    """Bias and RMSE in log10 space, grouped.

    Parameters
    ----------
    frame : DataFrame
        One row per (site, AEP) with ``estimate`` and ``observed`` log10 columns.
    by : sequence of str
        Grouping columns, e.g. ``["state", "region", "aep"]``.

    Returns
    -------
    DataFrame
        ``by`` columns plus ``n``, ``bias_log`` (mean of estimate - observed),
        ``rmse_log``, ``median_abs_log`` (robust to the few wild sites RMSE is
        dominated by) and ``bias_pct`` (``10**bias - 1``, in percent).

    Raises
    ------
    KeyError
        A named column is missing.
    """
    missing = [c for c in [*by, estimate, observed] if c not in frame.columns]
    if missing:
        raise KeyError(f"missing columns: {missing}")
    e = frame[estimate].astype(float) - frame[observed].astype(float)
    work = frame.loc[:, list(by)].copy()
    work["_e"] = e
    work["_e2"] = e * e
    g = work.groupby(list(by), dropna=False)
    out = g["_e"].agg(n="count", bias_log="mean").reset_index()
    out["rmse_log"] = np.sqrt(g["_e2"].mean().to_numpy())
    out["median_abs_log"] = g["_e"].apply(lambda s: float(np.median(np.abs(s)))).to_numpy()
    out["bias_pct"] = (10.0 ** out["bias_log"] - 1.0) * 100.0
    return out


@dataclass(frozen=True)
class DarSettings:
    """Donor eligibility for the drainage-area-ratio LOOCV.

    The ratio band is §5.1's ``0.5 <= A_u/A_g <= 1.5``; a donor must share the
    target's HUC8 (a cheap stand-in for the same-stream check §5.2 has not built),
    and the nearest eligible donor is used.
    """

    ratio_band: Tuple[float, float] = (0.5, 1.5)
    same_huc: int = 8
    max_distance_km: float = 100.0
    extra: Dict[str, Any] = field(default_factory=dict)


def _haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dp, dl = p2 - p1, math.radians(lon2 - lon1)
    a = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * 6371.0088 * math.asin(math.sqrt(a))


def dar_loocv(
    network: pd.DataFrame,
    exponents: Mapping[Tuple[str, float], float],
    aeps: Sequence[float],
    settings: DarSettings = DarSettings(),
) -> pd.DataFrame:
    """Leave-one-out drainage-area-ratio transposition over a gage network.

    Each gage in turn is the ungaged target: its nearest eligible donor's B17C
    quantile is moved to it by ``log Q_t = log Q_d + b(p) log10(A_t / A_d)``, and
    compared with the target's own B17C quantile. ``b(p)`` is the drainage-area
    exponent of the target's own regional equation at that AEP (§5.1).

    Parameters
    ----------
    network : DataFrame
        One row per gage: ``site_no``, ``state``, ``huc8``, ``latitude``,
        ``longitude``, ``drainage_area_sqmi`` and ``log_q_<aep>`` columns.
    exponents : mapping of (site_no, aep) to float
        ``b(p)`` for each target; a target without one is skipped.
    aeps : sequence of float
    settings : DarSettings

    Returns
    -------
    DataFrame
        One row per (target, AEP): ``site_no``, ``state``, ``donor``,
        ``area_ratio``, ``distance_km``, ``aep``, ``b``, ``log_est``, ``log_obs``.
    """
    rows: List[Dict[str, Any]] = []
    net = network.reset_index(drop=True)
    lo, hi = settings.ratio_band
    huc = net["huc8"].astype(str).str[: settings.same_huc].to_numpy()
    area = net["drainage_area_sqmi"].to_numpy(dtype=float)
    lat = net["latitude"].to_numpy(dtype=float)
    lon = net["longitude"].to_numpy(dtype=float)
    for i in range(len(net)):
        cands = []
        for j in np.flatnonzero(huc == huc[i]):
            if j == i:
                continue
            ratio = float(area[i] / area[j])
            if not lo <= ratio <= hi:
                continue
            dist = _haversine_km(lat[i], lon[i], lat[j], lon[j])
            if dist <= settings.max_distance_km:
                cands.append((dist, ratio, int(j)))
        if not cands:
            continue
        dist, ratio, jbest = min(cands, key=lambda c: c[0])
        t, d = net.iloc[i], net.iloc[jbest]
        for aep in aeps:
            col = f"log_q_{aep:g}"
            b = exponents.get((str(t["site_no"]), aep))
            if b is None or col not in net.columns:
                continue
            obs, don = t[col], d[col]
            if pd.isna(obs) or pd.isna(don):
                continue
            rows.append(
                {
                    "site_no": t["site_no"],
                    "state": t["state"],
                    "donor": d["site_no"],
                    "area_ratio": ratio,
                    "distance_km": dist,
                    "aep": aep,
                    "b": b,
                    "log_est": float(don) + b * math.log10(ratio),
                    "log_obs": float(obs),
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            "site_no",
            "state",
            "donor",
            "area_ratio",
            "distance_km",
            "aep",
            "b",
            "log_est",
            "log_obs",
        ],
    )


def area_exponent(lib: StateLibrary, region: str, aep: float) -> Optional[float]:
    """The drainage-area exponent of a region's equation (DRNAREA or CONTDA, log10).

    ``None`` when the region has no equation at that AEP or its area term is not a
    plain ``log10`` of area.
    """
    try:
        eq = lib.equation(region, aep)
    except EquationsUnavailable:
        return None
    for v, c in zip(eq.variables, eq.coefficients):
        if (
            v.code in ("DRNAREA", "CONTDA")
            and v.transform == "log10"
            and v.scale == 1.0
            and v.offset == 0.0
        ):
            return float(c)
    return None
