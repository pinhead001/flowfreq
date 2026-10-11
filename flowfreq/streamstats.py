"""
flowfreq.streamstats - USGS StreamStats watershed delineation and basin characteristics

Phase 1 (``docs/STREAMSTATS_MODULE_DESIGN.md``): given a pour point, snap it to the
stream network, delineate its watershed, and return basin characteristics (drainage
area, precipitation, canopy, ...) as regression predictors.

Phase 2 (``docs/STREAMSTATS_NSS_ADDENDUM.md``): feed those characteristics into NSS
(National Streamflow Statistics) to compute actual flow-statistic estimates (peak-flow,
low-flow, ...), each carrying its own regression equation and a resolved citation.
**Region selection across NSS's regressionRegions** is by the watershed polygon Phase 1
returns: :func:`locate_regression_regions` asks NSS which regions it falls in, with
area percentages (addendum S5), and :func:`estimate_flow_statistics` estimates only
those when given the polygon. Without one it returns every in-range region, unlabelled.

**Governing principle (design doc S4): a 200 is not an answer.** StreamStats can
delineate a hillslope sliver for an unsnappable point and return HTTP 200 with a
structurally valid polygon, signalled only by a ``WarningMsg`` string. Every response
here is validated against what it is supposed to contain; a result that fails
validation raises rather than being returned.

**The ``pourpoint`` snap response's ``output``, confirmed live 2026-09-11.** It is a
GeoJSON Point: ``{"type": "Point", "coordinates": [lon, lat]}`` -- GeoJSON coordinate
order is always longitude-then-latitude (RFC 7946 S3.1.1), the reverse of the
``(lat, lon)`` convention this module's own public API uses everywhere else.
:func:`_extract_point_lat_lon` handles this, trying the GeoJSON ``coordinates`` shape
first since that is what the live service actually returns, and falling back to a
handful of flat ``lat``/``lon``-style keys only in case some other region or a future
service version answers differently.

**The watershed polygon, verified live 2026-09-27 (design doc S10).** It was in the
``delineate/sshydro`` response all along, at ``bcrequest.wsresp.featurecollection[0]``
as the feature named ``globalwatershed`` -- the same place USGS's own published
workflow notebook reads it from. (An earlier version of this module called
``ss-delineate/v1/delineate/features/{region}`` for it and, called with the snapped
point, got back only a zero-area Point, so it was removed and ``polygon_geojson`` was
``None`` until this.) :func:`delineate_and_get_characteristics` now validates that
polygon before returning it -- rings closed, snapped point inside, geodesic area within
2% of ``DRNAREA`` -- and FR-3's ``WarningMsg`` check (:func:`_find_warning_msg`, a
recursive scan) is confirmed live to find the warning in the ``globalwatershedpoint``
and ``globalwatershed`` feature properties.
"""

from __future__ import annotations

import json
import logging
import math
import os
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import requests

logger = logging.getLogger(__name__)

#: Service versions this module was written and verified against (design doc S3,
#: exercised live 2026-09-09/10). The cache key includes these, so bumping them after
#: a service change invalidates stale cached results instead of silently reusing them
#: (NFR-1) -- update them here when ``/ss-delineate/openapi.json`` or
#: ``/ss-hydro/openapi.json`` report a new ``info.version``.
SS_DELINEATE_VERSION: str = "1.2.0"
SS_HYDRO_VERSION: str = "1.4.0"

#: Host with no server affinity, safe for the initial snap call and for the
#: region-list lookup, neither of which depends on SS-Delineate's temporary files.
GENERIC_HOST: str = "streamstats.usgs.gov"

#: Server used for the first SS-Delineate call of a point, before the service's own
#: ``usgswim-hostname`` response header is known. Whichever server actually answers is
#: then used for every subsequent call for that point (server stickiness, design doc
#: S3): SS-Hydro's basin-characteristics computation reads temporary files SS-Delineate
#: left on that specific node.
DEFAULT_SERVER: str = "prodweba"

#: Identifies this client to a free public service with no API key, per NFR-2.
USER_AGENT: str = "flowfreq/streamstats (https://github.com/pinhead001/flowfreq)"

#: StreamStats documents a hard limit of four simultaneous requests per user. This is
#: enforced as a ceiling on caller-requested concurrency, not a target -- NFR-2 asks for
#: serial-by-default besides.
MAX_CONCURRENCY: int = 4

#: Square metres per international square mile (1609.344 m squared).
SQ_M_PER_SQ_MI: float = 1609.344**2

#: Watershed-polygon validation tolerances (design doc S10, verified live 2026-09-27).
#: The polygon's geodesic area must match ss-hydro's own ``DRNAREA`` to within
#: ``POLYGON_AREA_RTOL`` relative, or ``POLYGON_AREA_ATOL_SQ_MI`` absolute, whichever is
#: larger. Observed live: 0.06-0.07% on Goat Creek, Lost River and the Methow at Pateros,
#: 0.0% on a 4.96 mi^2 WI basin -- the gap being ss-hydro's projected-CRS area and its
#: three-significant-figure rounding. 2% is ~30x that, loose enough never to trip on
#: rounding and tight enough to catch a polygon of the wrong basin.
POLYGON_AREA_RTOL: float = 0.02
POLYGON_AREA_ATOL_SQ_MI: float = 0.01

#: The snapped pour point is the outlet, so it sits on the basin boundary by
#: construction. Observed live it falls 15-21 m *inside* the polygon (one 10 m DEM cell
#: or two); it is accepted inside the polygon, or outside it by no more than this.
POUR_POINT_BOUNDARY_TOL_M: float = 100.0

#: Sphere radius for the polygon area: the WGS84 semi-major axis, the value the
#: Chamberlain-Duquette spherical approximation is normally paired with (as in
#: ``geojson-area``/Turf). Observed within 0.1% of ss-delineate's own projected
#: ``Shape_Area`` on every basin checked live.
_EARTH_RADIUS_M: float = 6378137.0


class UnsnappablePointError(ValueError):
    """A coordinate would not snap to the StreamStats stream network (FR-1).

    A data problem with the point itself, not a transient failure: retrying the same
    coordinate will not help. The caller must move the node.
    """


class UnsupportedRegionError(ValueError):
    """A region code is not one StreamStats recognizes, or rejected the request (FR-8).

    Never inferred from a coordinate by this module -- a wrong region would otherwise
    silently delineate against the wrong data layers.
    """


class DegenerateDelineationError(ValueError):
    """A delineation response failed validation: a warning, or a structurally invalid
    or degenerate polygon (FR-3). The accompanying geometry must not be used.
    """


class StreamStatsResponseError(ValueError):
    """A response did not have the shape this module expects.

    Distinct from a data problem with the request: this means a bug in this module, or
    that the service's response shape changed underneath it, and must fail loudly
    rather than return something that merely resembles a result.
    """


class StreamStatsTransportError(requests.RequestException):
    """A transport failure, timeout, or 5xx persisted after retrying. Retryable."""


@dataclass
class SnapResult:
    """Result of snapping a coordinate to the stream network (FR-1/FR-2).

    Attributes
    ----------
    requested_lat, requested_lon : float
        The coordinate as given by the caller.
    snapped_lat, snapped_lon : float
        The coordinate StreamStats actually delineates. The delineation belongs to
        this point, not the requested one.
    could_snap : bool
        Always True on a returned ``SnapResult`` -- construction raises
        :class:`UnsnappablePointError` otherwise, so a caller never has to remember to
        check this before trusting the snapped coordinates.
    distance_m : float
        Great-circle distance between requested and snapped coordinates, in meters.
        A data-quality signal the caller can threshold; this module does not impose a
        cutoff of its own (design doc S9.3 -- a snap tolerance is a caller decision).
    """

    requested_lat: float
    requested_lon: float
    snapped_lat: float
    snapped_lon: float
    could_snap: bool
    distance_m: float

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SnapResult":
        return cls(**data)


@dataclass
class Characteristic:
    """One StreamStats basin characteristic (FR-5/FR-6).

    Never flattened to a bare float: ``unit`` and ``msg`` are how a caller notices a
    characteristic came back in an unexpected unit or with a service-attached caveat
    (for example, a note about local versus total drainage area).
    """

    code: str
    name: str
    description: str
    value: float
    unit: str
    msg: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Characteristic":
        return cls(**data)


@dataclass
class Provenance:
    """Records how a :class:`WatershedCharacteristics` result was obtained (NFR-4).

    StreamStats updates its underlying data layers over time, so a characteristic
    obtained today is not guaranteed reproducible from the same coordinate next year.
    Without this, a published number obtained through this module cannot be defended.
    """

    service_versions: Dict[str, str]
    request_urls: List[str]
    requested_at_utc: str
    server_used: str

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Provenance":
        return cls(**data)


@dataclass
class WatershedCharacteristics:
    """A delineated watershed's basin characteristics, indexable by StreamStats code.

    Attributes
    ----------
    region : str
        The region code the caller supplied (FR-8 -- never inferred).
    snap : SnapResult
        The snap that preceded delineation.
    polygon_geojson : dict, optional
        The delineated watershed as a GeoJSON ``Feature`` (RFC 7946: WGS84,
        ``[lon, lat]`` order) whose geometry is a ``Polygon`` or ``MultiPolygon``,
        exactly as ``ss-delineate`` returned it as its ``globalwatershed`` feature, and
        validated before being returned: every ring closed, the snapped pour point
        inside it, and its area consistent with ``DRNAREA`` (see
        :func:`delineate_and_get_characteristics`). A plain dict, not a shapely
        geometry: NFR-6 keeps this module free of a hard GIS dependency. ``None`` when
        the caller passed ``include_polygon=False`` (the polygon is still validated),
        or for a result deserialized from a dict saved before polygons were returned.
    characteristics : dict of str to Characteristic
        Keyed by StreamStats parameter code (e.g. ``DRNAREA``, ``PRECPRIS10``).
    provenance : Provenance
        How and when this result was obtained.
    polygon_area_sq_mi : float, optional
        Geodesic area of the watershed polygon, square miles (:func:`geojson_area_sq_mi`).
        Set whenever a polygon was validated, including with ``include_polygon=False``.
    unavailable : dict of str to str
        Requested characteristics ss-hydro answered with its -999 "not found"
        sentinel, code to the service's message (design doc S11). Never in
        ``characteristics``.
    """

    region: str
    snap: SnapResult
    polygon_geojson: Optional[Dict[str, Any]]
    characteristics: Dict[str, Characteristic] = field(default_factory=dict)
    provenance: Optional[Provenance] = None
    polygon_area_sq_mi: Optional[float] = None
    unavailable: Dict[str, str] = field(default_factory=dict)

    def __getitem__(self, code: str) -> Characteristic:
        return self.characteristics[code]

    def __contains__(self, code: str) -> bool:
        return code in self.characteristics

    def to_dict(self) -> Dict[str, Any]:
        return {
            "region": self.region,
            "snap": self.snap.to_dict(),
            "polygon_geojson": self.polygon_geojson,
            "characteristics": {k: v.to_dict() for k, v in self.characteristics.items()},
            "provenance": self.provenance.to_dict() if self.provenance else None,
            "polygon_area_sq_mi": self.polygon_area_sq_mi,
            "unavailable": dict(self.unavailable),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "WatershedCharacteristics":
        provenance = data.get("provenance")
        return cls(
            region=data["region"],
            snap=SnapResult.from_dict(data["snap"]),
            polygon_geojson=data["polygon_geojson"],
            characteristics={
                code: Characteristic.from_dict(c) for code, c in data["characteristics"].items()
            },
            provenance=Provenance.from_dict(provenance) if provenance else None,
            polygon_area_sq_mi=data.get("polygon_area_sq_mi"),
            unavailable=dict(data.get("unavailable") or {}),
        )


def _default_cache_path() -> Path:
    return Path.home() / ".flowfreq" / "streamstats_cache.json"


class StreamStatsCache:
    """Disk-backed cache for delineation/characteristics results (NFR-1/NFR-5).

    Keyed on (region, requested lat/lon, service versions, requested characteristic
    codes) -- see :func:`_cache_key` for why the *requested* rather than snapped
    coordinate is used -- so a service-version bump or a different characteristic
    filter never silently reuses a stale entry. Loads on construction and writes
    through on every ``set``, so a populated cache works with no network at all,
    including no snap call.

    Safe to share between the worker threads of :func:`batch_get_characteristics`:
    ``set`` serializes under a lock and replaces the file atomically, so a crash
    mid-write leaves the previous file rather than a truncated one. It is not
    safe to share one cache *file* between processes -- the last writer wins.
    Every ``set`` rewrites the whole file, so with polygons included (tens of MB
    for a hundred basins) writes dominate a batch of cache misses; concurrent
    workers wait on each other's writes, not on StreamStats.
    """

    def __init__(self, path: Optional[Path] = None) -> None:
        self._path = path or _default_cache_path()
        self._data: Dict[str, Dict[str, Any]] = {}
        self._lock = threading.Lock()
        self._load()

    def _load(self) -> None:
        if self._path.exists():
            try:
                self._data = json.loads(self._path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                logger.warning(
                    "Could not read StreamStats cache at %s; starting empty",
                    self._path,
                    exc_info=True,
                )
                self._data = {}

    def get(self, key: str) -> Optional[Dict[str, Any]]:
        with self._lock:
            return self._data.get(key)

    def set(self, key: str, value: Dict[str, Any]) -> None:
        with self._lock:
            # Committed to memory only once it is on disk, so an entry that
            # cannot be serialized or written never stays behind to fail every
            # later write too. The copy is shallow: one dict of references.
            updated = {**self._data, key: value}
            # Compact: polygons make the file tens of MB, and indenting more
            # than doubled it and every write that rewrites it.
            text = json.dumps(updated, separators=(",", ":"))
            self._path.parent.mkdir(parents=True, exist_ok=True)
            fd, tmp_name = tempfile.mkstemp(
                dir=self._path.parent, prefix=f".{self._path.name}.", suffix=".tmp"
            )
            try:
                with os.fdopen(fd, "w", encoding="utf-8") as handle:
                    handle.write(text)
                os.replace(tmp_name, self._path)
            except BaseException:
                Path(tmp_name).unlink(missing_ok=True)
                raise
            self._data = updated


#: Bumped when a cached result's *content* changes shape, so an entry written by an
#: older version of this module is re-fetched rather than returned missing a field.
#: ``poly1``: results carry a validated watershed polygon (2026-09-27); entries
#: written before that have ``polygon_geojson=None`` and must not satisfy a lookup.
_CACHE_SCHEMA: str = "poly1"


def _cache_key(
    region: str, lat: float, lon: float, bc_labels: str, include_polygon: bool = True
) -> str:
    """Cache key from the *requested* coordinate, not the snapped one.

    The design doc (NFR-1) describes the key as snapped lat/lon, but keying on the
    snapped coordinate would force a snap call -- a network round trip -- before a
    cache hit could even be recognized, defeating NFR-5's stronger requirement that a
    populated cache complete a run with no network at all. Keying on what the caller
    already has in hand achieves both: service-version and characteristic-filter
    changes still invalidate stale entries, and a hit never touches the network.
    """
    return "|".join(
        [
            region,
            f"{round(lat, 6):.6f}",
            f"{round(lon, 6):.6f}",
            SS_DELINEATE_VERSION,
            SS_HYDRO_VERSION,
            bc_labels,
            _CACHE_SCHEMA,
            "polygon" if include_polygon else "nopolygon",
        ]
    )


def _haversine_distance_m(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance in meters, no geospatial dependency (NFR-6)."""
    r = 6371000.0
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


def _first_present(d: Dict[str, Any], keys: Sequence[str]) -> Optional[float]:
    """First of *keys* present in *d* with a numeric value, or None."""
    for key in keys:
        if key in d and d[key] is not None:
            try:
                return float(d[key])
            except (TypeError, ValueError):
                continue
    return None


def _extract_point_lat_lon(point: Dict[str, Any]) -> Tuple[Optional[float], Optional[float]]:
    """(lat, lon) from a ``pourpoint`` snap response's ``output`` object.

    Confirmed live 2026-09-11: it is a GeoJSON Point, ``{"type": "Point", "coordinates":
    [lon, lat]}`` -- GeoJSON coordinate order is always longitude-then-latitude (RFC
    7946 S3.1.1), the reverse of every other coordinate in this module's public API.
    Falls back to flat ``lat``/``lon``-style keys only if ``coordinates`` is absent, in
    case some other region or a future service version answers differently.
    """
    coords = point.get("coordinates")
    if isinstance(coords, (list, tuple)) and len(coords) >= 2:
        try:
            return float(coords[1]), float(coords[0])
        except (TypeError, ValueError):
            pass
    lat = _first_present(point, ("lat", "y", "snappedLat"))
    lon = _first_present(point, ("lon", "x", "snappedLon"))
    return lat, lon


def _request_with_backoff(
    method: Callable[..., requests.Response],
    url: str,
    *,
    params: Optional[Dict[str, Any]] = None,
    json_body: Optional[Union[Dict[str, Any], List[Any]]] = None,
    timeout: float,
    max_retries: int = 3,
    backoff_base: float = 2.0,
) -> requests.Response:
    """Issue one polite request with exponential backoff on 5xx or timeout (NFR-2).

    A 4xx is a data problem, not a transient one, and is never retried here --
    retrying a request that will fail the same way again just adds load to a free
    public service with no key.
    """
    last_exc: Optional[BaseException] = None
    for attempt in range(max_retries):
        try:
            kwargs: Dict[str, Any] = {"headers": {"User-Agent": USER_AGENT}, "timeout": timeout}
            if params is not None:
                kwargs["params"] = params
            if json_body is not None:
                kwargs["json"] = json_body
            response = method(url, **kwargs)
        except (requests.Timeout, requests.ConnectionError) as exc:
            last_exc = exc
            logger.warning(
                "StreamStats request to %s failed (attempt %d/%d): %s",
                url,
                attempt + 1,
                max_retries,
                exc,
            )
        else:
            if response.status_code >= 500:
                last_exc = requests.HTTPError(f"{response.status_code} error from {url}")
                logger.warning(
                    "StreamStats %s returned %d (attempt %d/%d)",
                    url,
                    response.status_code,
                    attempt + 1,
                    max_retries,
                )
            else:
                return response
        if attempt < max_retries - 1:
            time.sleep(backoff_base * (2**attempt))
    raise StreamStatsTransportError(
        f"StreamStats request to {url} failed after {max_retries} attempts"
    ) from last_exc


def _resolve_server(response: requests.Response, default: str) -> str:
    """Read the ``usgswim-hostname`` server-stickiness header (design doc S3).

    Falls back to *default* if the header is absent, rather than raising -- an
    undocumented header omission should not itself be treated as a validation
    failure the way a warning or a degenerate polygon is.
    """
    hostname = response.headers.get("usgswim-hostname")
    return hostname.lower() if hostname else default


def snap_point(region: str, lat: float, lon: float, *, timeout: float = 45.0) -> SnapResult:
    """Snap a coordinate to the StreamStats stream network (FR-1).

    Parameters
    ----------
    region : str
        StreamStats region code (e.g. ``"WA"``). Required, never inferred (FR-8).
    lat, lon : float
        Decimal-degree coordinate, WGS84.
    timeout : float
        Per-request timeout in seconds.

    Returns
    -------
    SnapResult

    Raises
    ------
    UnsnappablePointError
        The point will not snap. Move it; retrying will not help.
    UnsupportedRegionError
        StreamStats rejected the region code for this coordinate.
    StreamStatsResponseError
        The response was not valid JSON, or reported success with no usable output
        coordinates.
    StreamStatsTransportError
        The request failed after retrying.
    """
    url = f"https://{GENERIC_HOST}/pourpoint/v1/snap/str900"
    params = {"region": region, "lat": lat, "lon": lon}
    response = _request_with_backoff(requests.get, url, params=params, timeout=timeout)

    if response.status_code == 422:
        raise UnsupportedRegionError(
            f"StreamStats rejected region {region!r} while snapping ({lat}, {lon}): "
            f"{response.text}"
        )
    response.raise_for_status()

    try:
        data = response.json()
    except ValueError as exc:
        raise StreamStatsResponseError(
            f"Snap response for region {region!r} at ({lat}, {lon}) was not valid JSON"
        ) from exc

    could_snap = bool(data.get("couldSnap", False))
    if not could_snap:
        raise UnsnappablePointError(
            f"StreamStats could not snap ({lat}, {lon}) to the stream network in "
            f"region {region!r}; move the point rather than retrying"
        )

    output = data.get("output") or {}
    snapped_lat, snapped_lon = _extract_point_lat_lon(output)
    if snapped_lat is None or snapped_lon is None:
        raise StreamStatsResponseError(
            f"Snap response for region {region!r} reported couldSnap=true but no "
            f"usable output coordinates: {data!r}"
        )

    return SnapResult(
        requested_lat=lat,
        requested_lon=lon,
        snapped_lat=snapped_lat,
        snapped_lon=snapped_lon,
        could_snap=could_snap,
        distance_m=_haversine_distance_m(lat, lon, snapped_lat, snapped_lon),
    )


def _find_warning_msg(obj: Any) -> str:
    """Recursively search a JSON-like structure for a non-empty ``WarningMsg`` string.

    Confirmed live 2026-09-27 on the design doc's off-network point: the ``sshydro``
    response carries it in the ``properties`` of both the ``globalwatershedpoint`` and
    the ``globalwatershed`` features under ``bcrequest.wsresp.featurecollection[0]``.
    The whole structure is still scanned rather than those two paths, so a warning
    that moves elsewhere is caught rather than silently missed.
    """
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key.lower() == "warningmsg" and isinstance(value, str) and value.strip():
                return value.strip()
        for value in obj.values():
            found = _find_warning_msg(value)
            if found:
                return found
    elif isinstance(obj, list):
        for item in obj:
            found = _find_warning_msg(item)
            if found:
                return found
    return ""


def _polygons_of(geometry: Dict[str, Any]) -> List[List[List[List[float]]]]:
    """The polygons of a GeoJSON Polygon/MultiPolygon, each a list of rings."""
    gtype = geometry.get("type")
    coords = geometry.get("coordinates")
    if gtype == "Polygon" and isinstance(coords, list):
        return [coords]
    if gtype == "MultiPolygon" and isinstance(coords, list):
        return list(coords)
    raise DegenerateDelineationError(
        f"Watershed geometry is {gtype!r}, not a Polygon or MultiPolygon with coordinates"
    )


def _ring_area_m2(ring: Sequence[Sequence[float]]) -> float:
    """Unsigned area of one closed ``[lon, lat]`` ring on a sphere, square metres.

    The Chamberlain-Duquette approximation (JPL Publication 07-03, 2007) that
    ``geojson-area`` and Turf use: exact for a sphere up to the great-circle-vs-rhumb
    treatment of each edge, which on a basin's metre-scale DEM edges is negligible.
    """
    total = 0.0
    for (lon1, lat1, *_), (lon2, lat2, *_) in zip(ring[:-1], ring[1:]):
        total += math.radians(lon2 - lon1) * (
            2.0 + math.sin(math.radians(lat1)) + math.sin(math.radians(lat2))
        )
    return abs(total) * _EARTH_RADIUS_M**2 / 2.0


def geojson_area_sq_mi(geometry: Dict[str, Any]) -> float:
    """Geodesic area of a GeoJSON Polygon or MultiPolygon, in square miles.

    Coordinates are ``[lon, lat]`` in WGS84 (RFC 7946). Each polygon's first ring is its
    exterior and the rest are holes, which are subtracted. No GIS dependency (NFR-6).

    Parameters
    ----------
    geometry : dict
        A GeoJSON geometry object (not a Feature).

    Returns
    -------
    float
        Area in square miles.

    Raises
    ------
    DegenerateDelineationError
        *geometry* is not a Polygon or MultiPolygon.
    """
    area_m2 = 0.0
    for polygon in _polygons_of(geometry):
        if not polygon:
            continue
        area_m2 += _ring_area_m2(polygon[0]) - sum(_ring_area_m2(h) for h in polygon[1:])
    return area_m2 / SQ_M_PER_SQ_MI


def _point_in_ring(lon: float, lat: float, ring: Sequence[Sequence[float]]) -> bool:
    """Even-odd ray-casting test in lon/lat, adequate at basin scale."""
    inside = False
    for (x1, y1, *_), (x2, y2, *_) in zip(ring[:-1], ring[1:]):
        if (y1 > lat) != (y2 > lat) and lon < (x2 - x1) * (lat - y1) / (y2 - y1) + x1:
            inside = not inside
    return inside


def _point_in_geometry(lon: float, lat: float, geometry: Dict[str, Any]) -> bool:
    for polygon in _polygons_of(geometry):
        if polygon and _point_in_ring(lon, lat, polygon[0]):
            if not any(_point_in_ring(lon, lat, hole) for hole in polygon[1:]):
                return True
    return False


def _distance_to_boundary_m(lon: float, lat: float, geometry: Dict[str, Any]) -> float:
    """Shortest distance from a point to any ring edge, metres (local equirectangular)."""
    kx = math.cos(math.radians(lat)) * math.pi * _EARTH_RADIUS_M / 180.0
    ky = math.pi * _EARTH_RADIUS_M / 180.0
    best = math.inf
    for polygon in _polygons_of(geometry):
        for ring in polygon:
            for (x1, y1, *_), (x2, y2, *_) in zip(ring[:-1], ring[1:]):
                ax, ay = (x1 - lon) * kx, (y1 - lat) * ky
                dx, dy = (x2 - x1) * kx, (y2 - y1) * ky
                seg2 = dx * dx + dy * dy
                t = 0.0 if seg2 == 0 else max(0.0, min(1.0, -(ax * dx + ay * dy) / seg2))
                best = min(best, math.hypot(ax + t * dx, ay + t * dy))
    return best


def _extract_global_watershed(sshydro_data: Dict[str, Any]) -> Dict[str, Any]:
    """The ``globalwatershed`` GeoJSON Feature from a ``delineate/sshydro`` response.

    Verified live 2026-09-27 (design doc S10), and matching USGS's own published
    workflow notebook: the polygon is at
    ``bcrequest.wsresp.featurecollection[0]``, a list of ``{"name", "feature"}``
    entries, where the entry named ``globalwatershed`` holds a FeatureCollection whose
    Feature with ``properties.GlobalWshd == 1`` is the whole watershed.
    """
    try:
        collections = sshydro_data["bcrequest"]["wsresp"]["featurecollection"]
        entries = collections[0]
    except (KeyError, IndexError, TypeError) as exc:
        raise StreamStatsResponseError(
            "ss-delineate sshydro response carried no bcrequest.wsresp.featurecollection "
            "to take the watershed polygon from"
        ) from exc

    for entry in entries if isinstance(entries, list) else []:
        if not isinstance(entry, dict) or entry.get("name") != "globalwatershed":
            continue
        features = (entry.get("feature") or {}).get("features") or []
        whole = [
            f
            for f in features
            if isinstance(f, dict) and (f.get("properties") or {}).get("GlobalWshd") == 1
        ]
        if len(whole) == 1:
            return whole[0]
        if not whole and len(features) == 1 and isinstance(features[0], dict):
            return features[0]
        raise DegenerateDelineationError(
            f"globalwatershed carried {len(features)} features, {len(whole)} marked "
            "GlobalWshd=1; expected exactly one whole-watershed feature"
        )
    raise DegenerateDelineationError(
        "ss-delineate sshydro response carried no globalwatershed feature (FR-3)"
    )


def _validate_watershed_polygon(
    feature: Dict[str, Any], snapped_lat: float, snapped_lon: float
) -> float:
    """FR-3 geometry checks on the ``globalwatershed`` feature; returns its area (mi^2).

    Raises :class:`DegenerateDelineationError` unless every ring is closed with at least
    four positions, the area is positive, and the snapped pour point is inside the
    polygon or within :data:`POUR_POINT_BOUNDARY_TOL_M` of its boundary. The area-vs-
    ``DRNAREA`` check needs ss-hydro's answer and is done separately
    (:func:`_check_polygon_area`).
    """
    geometry = feature.get("geometry")
    if not isinstance(geometry, dict):
        raise DegenerateDelineationError("globalwatershed feature has no geometry")
    polygons = _polygons_of(geometry)
    if not polygons:
        raise DegenerateDelineationError("globalwatershed geometry has no polygons")
    for polygon in polygons:
        if not polygon:
            raise DegenerateDelineationError("globalwatershed polygon has no rings")
        for ring in polygon:
            if not isinstance(ring, list) or len(ring) < 4:
                raise DegenerateDelineationError(
                    "globalwatershed ring has fewer than 4 positions (RFC 7946 S3.1.6)"
                )
            if list(ring[0][:2]) != list(ring[-1][:2]):
                raise DegenerateDelineationError(
                    f"globalwatershed ring is not closed: first {ring[0]} != last {ring[-1]}"
                )

    area_sq_mi = geojson_area_sq_mi(geometry)
    if not area_sq_mi > 0.0:
        raise DegenerateDelineationError("globalwatershed polygon has zero area")

    if not _point_in_geometry(snapped_lon, snapped_lat, geometry):
        distance = _distance_to_boundary_m(snapped_lon, snapped_lat, geometry)
        if distance > POUR_POINT_BOUNDARY_TOL_M:
            raise DegenerateDelineationError(
                f"Snapped pour point ({snapped_lat}, {snapped_lon}) lies {distance:.0f} m "
                f"outside the delineated watershed (tolerance "
                f"{POUR_POINT_BOUNDARY_TOL_M:.0f} m); the polygon is not this point's basin"
            )
    return area_sq_mi


def _check_polygon_area(area_sq_mi: float, characteristics: Dict[str, Characteristic]) -> None:
    """Raise unless the polygon's area agrees with ss-hydro's own ``DRNAREA``.

    Skipped, with a log message, when ``DRNAREA`` was not among the characteristics
    returned: there is then nothing independent to compare against.
    """
    drnarea = characteristics.get("DRNAREA")
    if drnarea is None:
        logger.info("DRNAREA not returned; watershed polygon area check skipped")
        return
    unit = drnarea.unit.strip().lower()
    if unit in ("square miles", "mi^2", "sq mi", "square mile"):
        expected = drnarea.value
    elif unit in ("square kilometers", "square kilometres", "km^2", "sq km"):
        expected = drnarea.value * 1.0e6 / SQ_M_PER_SQ_MI
    else:
        raise StreamStatsResponseError(
            f"DRNAREA came back in unrecognised unit {drnarea.unit!r}; cannot check the "
            "watershed polygon's area against it"
        )
    tolerance = max(POLYGON_AREA_RTOL * expected, POLYGON_AREA_ATOL_SQ_MI)
    if abs(area_sq_mi - expected) > tolerance:
        raise DegenerateDelineationError(
            f"Watershed polygon area {area_sq_mi:.4g} mi^2 disagrees with DRNAREA "
            f"{expected:.4g} mi^2 by more than {tolerance:.3g} mi^2"
        )


#: ss-hydro's value for a requested characteristic it cannot compute (design doc S11):
#: HTTP 200, ``value: -999.0``, ``msg: "Basin Characteristic not found in database"``.
#: NSS uses -999.99 as its own "no value" placeholder; both are refused as data.
_MISSING_VALUE_SENTINELS: Tuple[float, ...] = (-999.0, -999.99)


def _parse_characteristics(hydro_data: Any) -> Tuple[Dict[str, Characteristic], Dict[str, str]]:
    """Parse the ss-hydro basin-characteristics response (design doc S3): a list of
    ``{name, description, code, unit, value, msg}`` objects, optionally wrapped in a
    ``{"parameters": [...]}`` envelope.

    Returns ``(characteristics, unavailable)``. An entry carrying the service's
    missing-value sentinel is never a characteristic: it goes into ``unavailable``,
    code to the service's message, so -999 can never reach a regression equation.
    """
    if isinstance(hydro_data, dict) and "parameters" in hydro_data:
        items = hydro_data["parameters"]
    elif isinstance(hydro_data, list):
        items = hydro_data
    else:
        raise StreamStatsResponseError(
            f"ss-hydro response was neither a list nor a dict with 'parameters': "
            f"{type(hydro_data).__name__}"
        )

    characteristics: Dict[str, Characteristic] = {}
    unavailable: Dict[str, str] = {}
    for item in items:
        try:
            code = str(item["code"])
            value = float(item["value"])
        except (KeyError, TypeError, ValueError) as exc:
            raise StreamStatsResponseError(
                f"ss-hydro characteristic entry missing code/value: {item!r}"
            ) from exc
        if value in _MISSING_VALUE_SENTINELS:
            unavailable[code] = str(item.get("msg") or f"value {value}")
            logger.warning("ss-hydro could not compute %s: %s", code, unavailable[code])
            continue
        characteristics[code] = Characteristic(
            code=code,
            name=str(item.get("name", "")),
            description=str(item.get("description", "")),
            value=value,
            unit=str(item.get("unit", "")),
            msg=str(item.get("msg") or ""),
        )
    return characteristics, unavailable


def delineate_and_get_characteristics(
    region: str,
    lat: float,
    lon: float,
    *,
    characteristic_codes: Optional[Sequence[str]] = None,
    cache: Optional[StreamStatsCache] = None,
    default_server: str = DEFAULT_SERVER,
    timeout: float = 60.0,
    include_polygon: bool = True,
) -> WatershedCharacteristics:
    """Snap, delineate, and compute basin characteristics for one pour point.

    Implements exactly the protocol design doc S3 verified live: a ``pourpoint`` snap
    (FR-1), ``ss-delineate``'s ``delineate/sshydro`` call to obtain the request body
    ``ss-hydro`` needs (validated for a ``WarningMsg`` per FR-3), and the ``ss-hydro``
    POST itself (FR-5/FR-6). Both ``ss-delineate``/``ss-hydro`` calls are pinned to the
    same server (design doc S3 -- ``ss-hydro`` reads temporary files ``ss-delineate``
    left on that specific node).

    The watershed polygon comes from the same ``delineate/sshydro`` response (design
    doc S10, verified live 2026-09-27) and is validated before anything is returned:
    every ring closed, the snapped pour point inside it (or within
    :data:`POUR_POINT_BOUNDARY_TOL_M` of the boundary it sits on), and its geodesic
    area within :data:`POLYGON_AREA_RTOL` of ss-hydro's own ``DRNAREA``. A polygon
    failing any of these raises :class:`DegenerateDelineationError`.

    Parameters
    ----------
    region : str
        StreamStats region code. Required, never inferred from the coordinate (FR-8).
    lat, lon : float
        Decimal-degree pour point, WGS84.
    characteristic_codes : sequence of str, optional
        StreamStats characteristic codes to request (resolves design doc S9.4).
        Sent as ss-hydro's ``BCs`` query parameter, semicolon-delimited; defaults to
        all (``BCs=*``). Every request states it explicitly rather than relying on a
        service default (FR-4). A code the region does not compute is left out of
        ``characteristics`` and listed in ``unavailable`` (design doc S11).
    cache : StreamStatsCache, optional
        When given, a cache hit skips all network calls; a miss populates it.
    default_server : str
        Server for the first ``ss-delineate`` call, before the service's own
        ``usgswim-hostname`` header is known.
    timeout : float
        Per-request timeout in seconds. StreamStats delineation and characteristics
        calls have been observed to take 6-11 s each under good conditions (NFR-3).
    include_polygon : bool
        Keep the validated polygon on the result (default). ``False`` still validates
        it and records its area, but drops the geometry itself -- a large basin's
        polygon is several hundred kB, which a cache of many points rewrites on every
        entry.

    Returns
    -------
    WatershedCharacteristics

    Raises
    ------
    UnsnappablePointError, UnsupportedRegionError, DegenerateDelineationError,
    StreamStatsResponseError, StreamStatsTransportError
        See each exception's docstring; see also S7 of the design doc.
    """
    bc_labels = ";".join(characteristic_codes) if characteristic_codes else "*"
    cache_key = _cache_key(region, lat, lon, bc_labels, include_polygon)

    if cache is not None:
        cached = cache.get(cache_key)
        if cached is not None:
            logger.info("StreamStats cache hit for %s", cache_key)
            return WatershedCharacteristics.from_dict(cached)

    snap = snap_point(region, lat, lon, timeout=timeout)

    request_urls: List[str] = []

    sshydro_url = (
        f"https://{default_server}.{GENERIC_HOST}/ss-delineate/v1/delineate/sshydro/{region}"
    )
    sshydro_params = {"lat": snap.snapped_lat, "lon": snap.snapped_lon}
    sshydro_response = _request_with_backoff(
        requests.get, sshydro_url, params=sshydro_params, timeout=timeout
    )
    request_urls.append(getattr(sshydro_response, "url", sshydro_url))
    if sshydro_response.status_code == 422:
        raise StreamStatsResponseError(
            f"ss-delineate sshydro request for region {region!r} at "
            f"({snap.snapped_lat}, {snap.snapped_lon}) was rejected: {sshydro_response.text}"
        )
    sshydro_response.raise_for_status()
    server_used = _resolve_server(sshydro_response, default_server)

    try:
        sshydro_data = sshydro_response.json()
    except ValueError as exc:
        raise StreamStatsResponseError(
            f"ss-delineate sshydro response for region {region!r} was not valid JSON"
        ) from exc

    warning_msg = _find_warning_msg(sshydro_data)
    if warning_msg:
        raise DegenerateDelineationError(
            f"StreamStats delineation warning at ({snap.snapped_lat}, "
            f"{snap.snapped_lon}) in region {region!r}: {warning_msg}"
        )

    bcrequest = sshydro_data.get("bcrequest")
    if not bcrequest:
        raise StreamStatsResponseError(
            f"ss-delineate sshydro response for region {region!r} carried no "
            f"bcrequest payload to pass to ss-hydro: {sshydro_data!r}"
        )

    # The watershed polygon rides in the same response (design doc S10). Geometry is
    # checked here, before the ss-hydro call is paid for; its area can only be checked
    # against DRNAREA once ss-hydro has answered.
    watershed_feature = _extract_global_watershed(sshydro_data)
    polygon_area_sq_mi = _validate_watershed_polygon(
        watershed_feature, snap.snapped_lat, snap.snapped_lon
    )

    hydro_url = (
        f"https://{server_used}.{GENERIC_HOST}"
        "/ss-hydro/v1/basin-characteristics/calculate-using-ssdelineate/"
    )
    hydro_params = {
        "region": region,
        "lat": snap.snapped_lat,
        "lon": snap.snapped_lon,
        # ss-hydro 1.4.0's OpenAPI names this `BCs`. The `bcLabels` this module sent until
        # 2026-09-28 was silently ignored, in the query string and in the body alike
        # (design doc S11).
        "BCs": bc_labels,
    }
    hydro_response = _request_with_backoff(
        requests.post, hydro_url, params=hydro_params, json_body=bcrequest, timeout=timeout
    )
    request_urls.append(getattr(hydro_response, "url", hydro_url))
    if hydro_response.status_code == 422:
        raise StreamStatsResponseError(
            f"ss-hydro rejected the basin-characteristics request for region "
            f"{region!r}: {hydro_response.text}"
        )
    hydro_response.raise_for_status()

    try:
        hydro_data = hydro_response.json()
    except ValueError as exc:
        raise StreamStatsResponseError(
            f"ss-hydro response for region {region!r} was not valid JSON"
        ) from exc

    characteristics, unavailable = _parse_characteristics(hydro_data)
    _check_polygon_area(polygon_area_sq_mi, characteristics)

    polygon_geojson: Optional[Dict[str, Any]] = None
    if include_polygon:
        polygon_geojson = {
            "type": "Feature",
            "geometry": watershed_feature["geometry"],
            "properties": dict(watershed_feature.get("properties") or {}),
        }

    provenance = Provenance(
        service_versions={"ss-delineate": SS_DELINEATE_VERSION, "ss-hydro": SS_HYDRO_VERSION},
        request_urls=request_urls,
        requested_at_utc=datetime.now(timezone.utc).isoformat(),
        server_used=server_used,
    )

    result = WatershedCharacteristics(
        region=region,
        snap=snap,
        polygon_geojson=polygon_geojson,
        characteristics=characteristics,
        provenance=provenance,
        unavailable=unavailable,
        polygon_area_sq_mi=polygon_area_sq_mi,
    )

    if cache is not None:
        cache.set(cache_key, result.to_dict())

    return result


def list_regions(*, timeout: float = 30.0) -> List[Dict[str, Any]]:
    """List StreamStats region codes (resolves design doc S9.2).

    Used to validate a caller's region codes up front rather than silently
    delineating against the wrong data layers for a mistyped one (FR-8). See
    :func:`batch_get_characteristics`'s ``validate_regions`` parameter.

    Returns
    -------
    list of dict
        Each carries at least an ``id``, ``name``, and ``code``.
    """
    url = f"https://{GENERIC_HOST}/nssservices/regions"
    response = _request_with_backoff(requests.get, url, timeout=timeout)
    response.raise_for_status()
    try:
        data = response.json()
    except ValueError as exc:
        raise StreamStatsResponseError("Regions list response was not valid JSON") from exc
    if not isinstance(data, list):
        raise StreamStatsResponseError(
            f"Expected a list of regions from {url}, got {type(data).__name__}"
        )
    return data


def batch_get_characteristics(
    points: Sequence[Tuple[str, str, float, float]],
    *,
    concurrency: int = 1,
    cache: Optional[StreamStatsCache] = None,
    characteristic_codes: Optional[Sequence[str]] = None,
    validate_regions: bool = True,
    timeout: float = 60.0,
    include_polygon: bool = True,
) -> Tuple[Dict[str, WatershedCharacteristics], Dict[str, str]]:
    """Delineate and fetch characteristics for many pour points (FR-7).

    A single bad point (unsnappable, wrong region, a service hiccup) never aborts the
    batch -- its reason is collected in ``errors`` instead, mirroring
    :func:`flowfreq.usgs.fetch_nwis_batch`'s two-dict return.

    Parameters
    ----------
    points : sequence of (point_id, region, lat, lon)
        ``point_id`` is any caller-chosen key (e.g. a reach node ID) used to key both
        returned dicts.
    concurrency : int
        Parallel workers. Serial (1) by default per NFR-2; capped at
        :data:`MAX_CONCURRENCY` regardless of what is requested, since StreamStats
        documents a hard limit of four simultaneous requests per user.
    cache : StreamStatsCache, optional
        Shared across all points in the batch.
    characteristic_codes : sequence of str, optional
        Passed through to :func:`delineate_and_get_characteristics` for every point.
    validate_regions : bool
        When True (default), fetch :func:`list_regions` once up front and raise
        :class:`UnsupportedRegionError` immediately for any point using an unknown
        region code, before issuing any delineation calls.
    timeout : float
        Per-request timeout in seconds, passed through to every point.
    include_polygon : bool
        Passed through to :func:`delineate_and_get_characteristics`.

    Returns
    -------
    tuple
        ``(results, errors)`` -- ``results`` maps ``point_id`` to
        :class:`WatershedCharacteristics` for points that succeeded; ``errors`` maps
        ``point_id`` to the failure reason for points that did not.
    """
    if concurrency > MAX_CONCURRENCY:
        logger.warning(
            "Requested concurrency %d exceeds StreamStats' documented "
            "4-simultaneous-request limit; capping at %d",
            concurrency,
            MAX_CONCURRENCY,
        )
    workers = max(1, min(concurrency, MAX_CONCURRENCY))

    if validate_regions:
        valid_codes = {str(r.get("code")) for r in list_regions(timeout=timeout)}
        for point_id, region, _, _ in points:
            if region not in valid_codes:
                raise UnsupportedRegionError(
                    f"Point {point_id!r} uses region {region!r}, not one of "
                    f"StreamStats' known region codes"
                )

    def _one(region: str, lat: float, lon: float) -> WatershedCharacteristics:
        return delineate_and_get_characteristics(
            region,
            lat,
            lon,
            characteristic_codes=characteristic_codes,
            cache=cache,
            timeout=timeout,
            include_polygon=include_polygon,
        )

    results: Dict[str, WatershedCharacteristics] = {}
    errors: Dict[str, str] = {}

    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_id = {
            executor.submit(_one, region, lat, lon): point_id
            for point_id, region, lat, lon in points
        }
        for future in as_completed(future_to_id):
            point_id = future_to_id[future]
            try:
                results[point_id] = future.result()
            except Exception as e:
                errors[point_id] = str(e)

    return results, errors


# =============================================================================
# Phase 2 -- NSS (National Streamflow Statistics) flow-statistics estimation.
# See docs/STREAMSTATS_NSS_ADDENDUM.md for the live-verification pass this
# implements. The real protocol differs from every NSS-related detail the
# ff-idea02 PDF/py transcript guessed: GET /nssservices/regions/{region} carries
# no statistic groups, POST /nssservices/estimate does not exist at all, and the
# real estimate call (POST /nssservices/Scenarios/Estimate) takes a bare JSON
# array body, not one wrapped in {"scenarioList": [...]} -- the addendum's S2
# records how each of those was found.
# =============================================================================


@dataclass
class RegressionCitation:
    """A published citation for a regression-equation set (addendum S2/S4 item 3).

    Resolves the "uncited exponent" concern the design doc's FR-9 raised for a
    later NSS module: every regression estimate can be traced to a real,
    DOI-linked reference via this, not just an opaque ``citationID``.
    """

    citation_id: int
    title: str
    author: str
    citation_url: str
    last_year_of_data: Optional[int] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RegressionCitation":
        return cls(**data)


@dataclass
class FlowStatisticEstimate:
    """One computed flow statistic from one NSS regression region.

    Attributes
    ----------
    equation : str
        The literal regression equation NSS used, e.g.
        ``"3.846*DRNAREA^0.745*10^(0.032*PRECPRIS10)/10^(0.0078*CANOPY_PCT)"`` --
        returned directly by the service, not reconstructed.
    standard_error_pct : float, optional
        NSS's own standard error, in percent, from the confusingly named ``errors``
        list of the raw response (a prediction-error statistic, not a failure
        indicator). Peak-flow equations report it as ``ASEp`` (average standard error
        of prediction); WA's low-flow equations as ``SE`` ("average standard error
        (of either estimate or prediction)"), confirmed live 2026-09-28 (addendum
        S6). ``None`` when NSS returned neither -- observed live for a wildly
        out-of-range input (addendum S3), so its absence is a soft signal worth
        noticing even though it is not a hard failure marker.
    standard_error_code : str, optional
        Which of those two NSS codes ``standard_error_pct`` came from.
    """

    code: str
    name: str
    description: str
    value: float
    unit: str
    equation: str
    standard_error_pct: Optional[float] = None
    interval_lower: Optional[float] = None
    interval_upper: Optional[float] = None
    standard_error_code: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "FlowStatisticEstimate":
        return cls(**data)


@dataclass
class RegionFlowEstimates:
    """Every flow statistic NSS computed for one regressionRegion, one statistic group.

    A state commonly defines several independently-calibrated regressionRegions per
    statistic group (WA Peak-Flow: four) -- this is one of them, not "the" answer
    for a point. ``located`` says whether NSS placed the watershed in it (``None`` when
    no polygon was given); see :func:`estimate_flow_statistics` on region selection.
    """

    region_code: str
    region_name: str
    statistic_group_code: str
    statistic_group_name: str
    estimates: Dict[str, FlowStatisticEstimate] = field(default_factory=dict)
    citation: Optional[RegressionCitation] = None
    #: ``True`` when NSS placed the watershed polygon in this region, ``False`` when it
    #: did not, ``None`` when no polygon was supplied (location was never checked).
    located: Optional[bool] = None
    #: NSS's own area percentage of the watershed in this region (``percentWeight``
    #: from ``bylocation``); ``None`` unless located. For the area-averaged result,
    #: ``None`` too -- it is the whole basin.
    percent_weight: Optional[float] = None
    #: ``True`` only for NSS's own ``areaave`` result: the percent-weighted mean of the
    #: located regions' estimates, returned when the basin spans several regions.
    area_averaged: bool = False

    def __getitem__(self, code: str) -> FlowStatisticEstimate:
        return self.estimates[code]

    def __contains__(self, code: str) -> bool:
        return code in self.estimates

    def to_dict(self) -> Dict[str, Any]:
        return {
            "region_code": self.region_code,
            "region_name": self.region_name,
            "statistic_group_code": self.statistic_group_code,
            "statistic_group_name": self.statistic_group_name,
            "estimates": {k: v.to_dict() for k, v in self.estimates.items()},
            "citation": self.citation.to_dict() if self.citation else None,
            "located": self.located,
            "percent_weight": self.percent_weight,
            "area_averaged": self.area_averaged,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RegionFlowEstimates":
        citation = data.get("citation")
        return cls(
            region_code=data["region_code"],
            region_name=data["region_name"],
            statistic_group_code=data["statistic_group_code"],
            statistic_group_name=data["statistic_group_name"],
            estimates={
                code: FlowStatisticEstimate.from_dict(e) for code, e in data["estimates"].items()
            },
            citation=RegressionCitation.from_dict(citation) if citation else None,
            located=data.get("located"),
            percent_weight=data.get("percent_weight"),
            area_averaged=bool(data.get("area_averaged", False)),
        )


@dataclass
class RegressionRegionLocation:
    """One NSS regression region a watershed polygon overlaps (addendum S5).

    Attributes
    ----------
    percent_weight : float
        NSS's own percentage of the watershed's area in this region, as NSS reports it
        (rounded to a whole percent). It is what NSS weights by when area-averaging.
    area_sq_mi : float
        The overlap area, square miles, unrounded.
    """

    region_id: int
    code: str
    name: str
    percent_weight: float
    area_sq_mi: float
    citation_id: Optional[int] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "RegressionRegionLocation":
        return cls(**data)


#: How far NSS's rounded ``percentWeight`` may sit from ``100 * area / basin area``
#: before a ``bylocation`` answer is refused: 0.5 for NSS's whole-percent rounding plus
#: 1.0 for the gap between NSS's projected area and this module's spherical one
#: (0.3% at most live). Observed live: 89/11 against 89.46/10.54 (Goat Creek, national
#: regions), 44/32/24 against 44.12/31.86/24.01 (Ogeechee, GA).
_PERCENT_WEIGHT_TOL: float = 1.5

#: How close the ``percentWeight`` of the located regions of one statistic group must
#: sum to 100 for NSS to area-average them. Confirmed live: 60+40 returns an
#: ``areaave`` result; 60+30 returns none, and no error either.
_WEIGHT_SUM_TOL: float = 0.5


def _geometry_of(watershed_polygon: Dict[str, Any]) -> Dict[str, Any]:
    """The geometry of a GeoJSON Feature, or the object itself if it is a geometry."""
    if watershed_polygon.get("type") == "Feature":
        geometry = watershed_polygon.get("geometry")
        if not isinstance(geometry, dict):
            raise ValueError("watershed_polygon is a Feature with no geometry")
        return geometry
    return watershed_polygon


def locate_regression_regions(
    region: str,
    watershed_polygon: Dict[str, Any],
    *,
    timeout: float = 60.0,
) -> List[RegressionRegionLocation]:
    """The NSS regression regions a watershed polygon falls in, with area percentages.

    ``POST /nssservices/regions/{region}/regressionregions/bylocation`` with the bare
    GeoJSON geometry as the body, verified live 2026-09-27 (addendum S5). Scoped to
    *region*: the unscoped ``/regressionregions/bylocation`` also returns national
    studies and neighbouring states' regions (Oregon's low-flow regions for a WA basin
    on the Columbia), none of which a WA scenario can use.

    Every answer is checked against the polygon it was asked about: each region's
    overlap area may not exceed the basin's by more than :data:`POLYGON_AREA_RTOL`, and
    each ``percentWeight`` must agree with ``100 * area / basin area`` to within
    1.5 percentage points.

    Parameters
    ----------
    region : str
        StreamStats/NSS region code, e.g. ``"WA"``. Never inferred.
    watershed_polygon : dict
        A GeoJSON Polygon/MultiPolygon, or a Feature carrying one -- for instance
        :attr:`WatershedCharacteristics.polygon_geojson`. NSS refuses a bare Point
        (``400 "Geometry is not of type: Polygon,MultiPolygon"``, confirmed live).

    Returns
    -------
    list of RegressionRegionLocation
        Possibly empty: a basin outside every region that has geometry in NSS. Some
        regions have none (WA's low-flow ``GC1434``, ``GC1556``, ``GC1558``,
        ``statusID`` 3) and so can never be located.

    Raises
    ------
    StreamStatsResponseError
        A rejected request, or an answer that fails the checks above.
    StreamStatsTransportError
        Transport failure after retrying.
    """
    geometry = _geometry_of(watershed_polygon)
    basin_sq_mi = geojson_area_sq_mi(geometry)
    url = f"https://{GENERIC_HOST}/nssservices/regions/{region}/regressionregions/bylocation"
    response = _request_with_backoff(requests.post, url, json_body=geometry, timeout=timeout)
    if response.status_code >= 400:
        raise StreamStatsResponseError(
            f"NSS bylocation for region {region!r} was rejected "
            f"({response.status_code}): {response.text}"
        )
    try:
        data = response.json()
    except ValueError as exc:
        raise StreamStatsResponseError(
            f"NSS bylocation response for region {region!r} was not valid JSON"
        ) from exc
    if not isinstance(data, list):
        raise StreamStatsResponseError(
            f"Expected a list of regression regions from {url}, got {type(data).__name__}"
        )

    located: List[RegressionRegionLocation] = []
    for item in data:
        try:
            loc = RegressionRegionLocation(
                region_id=int(item["id"]),
                code=str(item["code"]),
                name=str(item.get("name", "")),
                percent_weight=float(item["percentWeight"]),
                area_sq_mi=float(item["area"]),
                citation_id=(
                    int(item["citationID"]) if item.get("citationID") is not None else None
                ),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise StreamStatsResponseError(
                f"NSS bylocation entry missing id/code/percentWeight/area: {item!r}"
            ) from exc
        if loc.area_sq_mi < 0 or loc.area_sq_mi > basin_sq_mi * (1.0 + POLYGON_AREA_RTOL):
            raise StreamStatsResponseError(
                f"NSS placed {loc.area_sq_mi:.4g} mi^2 of a {basin_sq_mi:.4g} mi^2 basin in "
                f"region {loc.code}; the answer is not about this polygon"
            )
        implied = 100.0 * loc.area_sq_mi / basin_sq_mi
        if abs(loc.percent_weight - implied) > _PERCENT_WEIGHT_TOL:
            raise StreamStatsResponseError(
                f"NSS gave region {loc.code} percentWeight {loc.percent_weight} but an "
                f"area of {loc.area_sq_mi:.4g} of {basin_sq_mi:.4g} mi^2 ({implied:.1f}%)"
            )
        located.append(loc)
    return located


#: NSS ``errors`` codes read as the standard error, in order of preference: ``ASEp``
#: on peak-flow equations, ``SE`` on WA's low-flow ones (both confirmed live).
_STANDARD_ERROR_CODES: Tuple[str, ...] = ("ASEp", "SE")


def list_statistic_groups(
    region: Optional[str] = None, *, timeout: float = 30.0
) -> List[Dict[str, Any]]:
    """List NSS statistic groups (addendum S2).

    With ``region=None``, lists every statistic group NSS defines (``PC``, ``PFS``,
    ``LFS``, ``FDS``, ...). With ``region`` given, lists only the groups actually
    valid for that region (confirmed live: WA supports only ``PFS`` and ``LFS``) --
    the real replacement for the ff-idea02 transcript's wrong assumption that
    ``GET /nssservices/regions/{region}`` itself carries this list.

    Returns
    -------
    list of dict
        Each carries at least ``id``, ``name``, ``code``, ``defType``.
    """
    if region:
        url = f"https://{GENERIC_HOST}/nssservices/regions/{region}/statisticgroups"
    else:
        url = f"https://{GENERIC_HOST}/nssservices/statisticgroups"
    response = _request_with_backoff(requests.get, url, timeout=timeout)
    response.raise_for_status()
    try:
        data = response.json()
    except ValueError as exc:
        raise StreamStatsResponseError("Statistic groups response was not valid JSON") from exc
    if not isinstance(data, list):
        raise StreamStatsResponseError(
            f"Expected a list of statistic groups from {url}, got {type(data).__name__}"
        )
    return data


def _fill_and_validate_regions(
    scenario: Dict[str, Any], characteristics: Dict[str, Characteristic]
) -> Tuple[Dict[str, Any], Dict[str, str]]:
    """Fill a scenario template's parameters from real characteristics.

    Keeps only the regressionRegions whose every required parameter is present in
    *characteristics* and within that region's own declared ``[min, max]``
    (addendum S3): NSS will not check this for you, and returns a plausible,
    silently-wrong extrapolated value for an out-of-range input with no warning at
    all -- worse than Phase 1's ``WarningMsg``-bearing failures, which at least
    signal something. Validating before submission, using limits the service
    itself just echoed back, means an invalid region is never even asked for --
    stronger than catching a bad answer after the fact.

    Returns
    -------
    tuple
        ``(scenario_with_only_valid_regions, {region_code: skip_reason})``.
    """
    valid_regions = []
    skipped: Dict[str, str] = {}
    for region in scenario.get("regressionRegions", []):
        region_code = str(region.get("code", "?"))
        ok = True
        for param in region.get("parameters", []):
            code = param.get("code")
            char = characteristics.get(code)
            if char is None:
                skipped[region_code] = f"no basin characteristic for required parameter {code!r}"
                ok = False
                break
            limits = param.get("limits") or {}
            lo, hi = limits.get("min"), limits.get("max")
            if (lo is not None and char.value < lo) or (hi is not None and char.value > hi):
                skipped[region_code] = (
                    f"{code}={char.value} is outside this region's valid range [{lo}, {hi}]"
                )
                ok = False
                break
            param["value"] = char.value
        if ok:
            valid_regions.append(region)
    filled = dict(scenario)
    filled["regressionRegions"] = valid_regions
    return filled, skipped


def _fetch_citations(region_ids: Sequence[int], *, timeout: float) -> Dict[int, RegressionCitation]:
    """Resolve citations for a set of regressionRegion IDs, keyed by citation ID.

    Confirmed live: a returned citation's own ``id`` is the same value each
    regressionRegion's own ``citationID`` field references, and several regions
    commonly share one citation (WA's four Peak-Flow regions all resolve to a
    single Mastin et al. 2016 citation) -- the response lists each distinct
    citation once, matched back to a region by ``citationID``, not by request
    order or count.
    """
    url = f"https://{GENERIC_HOST}/nssservices/citations"
    params = {"regressionregions": ",".join(str(r) for r in region_ids)}
    response = _request_with_backoff(requests.get, url, params=params, timeout=timeout)
    response.raise_for_status()
    try:
        data = response.json()
    except ValueError as exc:
        raise StreamStatsResponseError("NSS citations response was not valid JSON") from exc
    if not isinstance(data, list):
        raise StreamStatsResponseError(f"Expected a list of citations, got {type(data).__name__}")
    return {
        int(c["id"]): RegressionCitation(
            citation_id=int(c["id"]),
            title=str(c.get("title", "")),
            author=str(c.get("author", "")),
            citation_url=str(c.get("citationURL", "")),
            last_year_of_data=c.get("lastYearOfData"),
        )
        for c in data
    }


def _parse_region_results(rr: Dict[str, Any]) -> Dict[str, FlowStatisticEstimate]:
    """Parse one regressionRegion's ``results`` from a ``Scenarios/Estimate`` answer."""
    estimates: Dict[str, FlowStatisticEstimate] = {}
    for stat in rr.get("results", []) or []:
        errors = {e.get("code"): e.get("value") for e in stat.get("errors") or []}
        se_code = next((c for c in _STANDARD_ERROR_CODES if c in errors), None)
        sep = errors[se_code] if se_code is not None else None
        bounds = stat.get("intervalBounds") or {}
        try:
            code = str(stat["code"])
            value = float(stat["value"])
        except (KeyError, TypeError, ValueError) as exc:
            raise StreamStatsResponseError(
                f"NSS result entry missing code/value: {stat!r}"
            ) from exc
        estimates[code] = FlowStatisticEstimate(
            code=code,
            name=str(stat.get("name", "")),
            description=str(stat.get("description", "")),
            value=value,
            unit=str((stat.get("unit") or {}).get("abbr", "")),
            equation=str(stat.get("equation", "")),
            standard_error_pct=sep,
            standard_error_code=se_code,
            interval_lower=bounds.get("lower"),
            interval_upper=bounds.get("upper"),
        )
    return estimates


#: NSS's own area-averaged pseudo-region in a ``Scenarios/Estimate`` answer (live).
_AREA_AVERAGED_CODE: str = "areaave"


def _check_area_average(
    average: RegionFlowEstimates,
    weighted: Sequence[RegionFlowEstimates],
    weights: Dict[str, float],
) -> None:
    """Raise unless NSS's ``areaave`` values are the percent-weighted mean it claims.

    Confirmed live (addendum S5) that NSS computes ``sum(w_i * Q_i) / 100`` over the
    values it returns for each region -- 0.6*116 + 0.4*33.3 = 82.92 exactly. Checked
    here to 0.5% rather than trusted.
    """
    for code, est in average.estimates.items():
        parts = [(weights[r.region_code], r.estimates.get(code)) for r in weighted]
        if any(p is None for _, p in parts):
            raise StreamStatsResponseError(
                f"NSS area-averaged {code} has no counterpart in every weighted region"
            )
        expected = sum(w * p.value for w, p in parts if p is not None) / 100.0
        if not math.isclose(est.value, expected, rel_tol=0.005, abs_tol=1e-9):
            raise StreamStatsResponseError(
                f"NSS area-averaged {code} = {est.value} but the percent-weighted mean of "
                f"its regions is {expected:.6g}"
            )


def estimate_flow_statistics(
    region: str,
    characteristics: Dict[str, Characteristic],
    statistic_group_codes: Optional[Sequence[str]] = None,
    *,
    unit_system: int = 2,
    timeout: float = 60.0,
    watershed_polygon: Optional[Dict[str, Any]] = None,
    include_unlocated: bool = False,
) -> Tuple[List[RegionFlowEstimates], Dict[str, str]]:
    """Estimate NSS regression flow statistics from real basin characteristics.

    Phase 2 (``docs/STREAMSTATS_NSS_ADDENDUM.md``). Given characteristics Phase 1's
    :func:`delineate_and_get_characteristics` already produced, computes the
    regression equations NSS defines for the region, skipping (never silently
    computing) any regressionRegion whose required parameters are missing or fall
    outside that region's own declared valid range (addendum S3).

    **Region selection.** NSS defines several independently calibrated
    regressionRegions per statistic group within a state (WA Peak-Flow has four) and
    computes an answer for every one it is handed, whether or not the basin is in it.
    Pass ``watershed_polygon`` (e.g. :attr:`WatershedCharacteristics.polygon_geojson`)
    and the regions are resolved by :func:`locate_regression_regions` (addendum S5):

    - only the regions the polygon falls in are estimated, each labelled
      ``located=True`` with NSS's own area ``percent_weight``;
    - when the located regions of one statistic group are all in range and their
      weights sum to 100, NSS's own area-weighted mean comes back as an extra
      :class:`RegionFlowEstimates` with ``area_averaged=True`` (code ``areaave``),
      checked against the weighted mean of the regions' values;
    - ``include_unlocated=True`` keeps the other in-range regions too, labelled
      ``located=False``, and then requests no area average (NSS computes none for a
      mix of weighted and unweighted regions -- confirmed live).

    Without a polygon the behaviour is unchanged from before: every in-range region is
    returned with ``located=None``, and choosing among them is the caller's job.

    Parameters
    ----------
    region : str
        StreamStats region code (e.g. ``"WA"``). Required, never inferred.
    characteristics : dict of str to Characteristic
        Real basin characteristics, e.g. ``WatershedCharacteristics.characteristics``
        from :func:`delineate_and_get_characteristics`.
    statistic_group_codes : sequence of str, optional
        NSS statistic group codes to estimate (e.g. ``["PFS"]``). Defaults to every
        group :func:`list_statistic_groups` reports as valid for this region -- a
        live lookup, not a guess.
    unit_system : int
        1=Metric, 2=US Customary (default), 3=Universal.
    timeout : float
        Per-request timeout in seconds.
    watershed_polygon : dict, optional
        The basin as a GeoJSON Polygon/MultiPolygon or a Feature carrying one.
    include_unlocated : bool
        With a polygon, also return in-range regions the basin is not in.

    Returns
    -------
    tuple
        ``(region_estimates, skipped)`` -- ``region_estimates`` is a list of
        :class:`RegionFlowEstimates`, one per (statistic group, regressionRegion)
        estimated, plus any area-averaged result; ``skipped`` maps
        ``"{statistic_group_code}:{region_code}"`` to why that region was excluded,
        and ``"{statistic_group_code}:areaave"`` to why no area average was requested
        for a basin that spans several regions. One bad regression region never
        excludes its siblings (FR-7's batch shape, one level deeper).

    Raises
    ------
    StreamStatsResponseError, StreamStatsTransportError
        A malformed response, or a transport failure after retrying.
    """
    if statistic_group_codes is None:
        statistic_group_codes = [g["code"] for g in list_statistic_groups(region, timeout=timeout)]

    id_to_code = {g["id"]: g["code"] for g in list_statistic_groups(timeout=timeout)}

    locations: Optional[Dict[str, RegressionRegionLocation]] = None
    if watershed_polygon is not None:
        locations = {
            loc.code: loc
            for loc in locate_regression_regions(region, watershed_polygon, timeout=timeout)
        }

    templates_url = f"https://{GENERIC_HOST}/nssservices/regions/{region}/Scenarios"
    templates_params = {
        "statisticgroups": ",".join(statistic_group_codes),
        "unitsystem": unit_system,
    }
    templates_response = _request_with_backoff(
        requests.get, templates_url, params=templates_params, timeout=timeout
    )
    templates_response.raise_for_status()
    try:
        templates = templates_response.json()
    except ValueError as exc:
        raise StreamStatsResponseError(
            f"NSS scenario-template response for region {region!r} was not valid JSON"
        ) from exc
    if not isinstance(templates, list):
        raise StreamStatsResponseError(
            f"Expected a list of scenario templates for region {region!r}, got "
            f"{type(templates).__name__}"
        )

    skipped: Dict[str, str] = {}
    to_submit: List[Dict[str, Any]] = []
    #: statistic group code -> {region code: percent weight} submitted for averaging
    averaged_weights: Dict[str, Dict[str, float]] = {}
    for scenario in templates:
        group_code = id_to_code.get(scenario.get("statisticGroupID"), "")
        template_regions = scenario.get("regressionRegions", [])
        if locations is not None and not include_unlocated:
            for rr in template_regions:
                code = str(rr.get("code", "?"))
                if code not in locations:
                    skipped[f"{group_code}:{code}"] = (
                        "the watershed is not in this regression region (NSS bylocation)"
                    )
            scenario = dict(scenario)
            scenario["regressionRegions"] = [
                rr for rr in template_regions if str(rr.get("code")) in locations
            ]
        filled, region_skips = _fill_and_validate_regions(scenario, characteristics)
        for region_code, reason in region_skips.items():
            skipped[f"{group_code}:{region_code}"] = reason
        if not filled["regressionRegions"]:
            continue

        if locations is not None and not include_unlocated:
            located_here = [str(rr.get("code")) for rr in scenario["regressionRegions"]]
            if len(located_here) > 1:
                weights = {c: locations[c].percent_weight for c in located_here}
                total = sum(weights.values())
                if region_skips:
                    skipped[f"{group_code}:areaave"] = (
                        f"basin spans {', '.join(located_here)} but "
                        f"{', '.join(sorted(region_skips))} could not be estimated, so no "
                        "area-weighted average"
                    )
                elif abs(total - 100.0) > _WEIGHT_SUM_TOL:
                    skipped[f"{group_code}:areaave"] = (
                        f"located regions' weights sum to {total:g}, not 100, so NSS "
                        "would not area-average them"
                    )
                else:
                    for rr in filled["regressionRegions"]:
                        rr["percentWeight"] = weights[str(rr.get("code"))]
                    averaged_weights[group_code] = weights
        to_submit.append(filled)

    if not to_submit:
        return [], skipped

    estimate_url = f"https://{GENERIC_HOST}/nssservices/Scenarios/Estimate"
    estimate_response = _request_with_backoff(
        requests.post, estimate_url, json_body=to_submit, timeout=timeout
    )
    if estimate_response.status_code >= 400:
        raise StreamStatsResponseError(
            f"NSS estimate request for region {region!r} was rejected "
            f"({estimate_response.status_code}): {estimate_response.text}"
        )
    try:
        results = estimate_response.json()
    except ValueError as exc:
        raise StreamStatsResponseError(
            f"NSS estimate response for region {region!r} was not valid JSON"
        ) from exc
    if not isinstance(results, list):
        raise StreamStatsResponseError(
            f"Expected a list of scenario results for region {region!r}, got "
            f"{type(results).__name__}"
        )

    parsed: List[Tuple[Optional[int], RegionFlowEstimates]] = []
    region_ids: List[int] = []
    for scenario in results:
        group_code = id_to_code.get(scenario.get("statisticGroupID"), "")
        group_name = str(scenario.get("statisticGroupName", ""))
        group_items: List[RegionFlowEstimates] = []
        average: Optional[RegionFlowEstimates] = None
        for rr in scenario.get("regressionRegions", []):
            code = str(rr.get("code", ""))
            item = RegionFlowEstimates(
                region_code=code,
                region_name=str(rr.get("name", "")),
                statistic_group_code=group_code,
                statistic_group_name=group_name,
                estimates=_parse_region_results(rr),
            )
            if code == _AREA_AVERAGED_CODE:
                if group_code not in averaged_weights:
                    raise StreamStatsResponseError(
                        f"NSS returned an area average for {group_code} that was not asked for"
                    )
                item.area_averaged = True
                item.located = True
                average = item
                parsed.append((None, item))
                continue
            if locations is not None:
                loc = locations.get(code)
                item.located = loc is not None
                item.percent_weight = loc.percent_weight if loc is not None else None
            region_id = rr.get("id")
            if isinstance(region_id, int) and region_id > 0:
                region_ids.append(region_id)
            group_items.append(item)
            parsed.append((rr.get("citationID"), item))

        if group_code in averaged_weights:
            if average is None:
                raise StreamStatsResponseError(
                    f"NSS returned no area average for {group_code} although the located "
                    "regions were submitted with weights summing to 100"
                )
            _check_area_average(average, group_items, averaged_weights[group_code])

    citations_by_id: Dict[int, RegressionCitation] = {}
    if region_ids:
        try:
            citations_by_id = _fetch_citations(region_ids, timeout=timeout)
        except (StreamStatsResponseError, StreamStatsTransportError, requests.RequestException):
            logger.warning(
                "Could not resolve NSS citations for regions %s", region_ids, exc_info=True
            )

    region_estimates: List[RegionFlowEstimates] = []
    for citation_id, item in parsed:
        if citation_id is not None:
            item.citation = citations_by_id.get(citation_id)
        region_estimates.append(item)

    return region_estimates, skipped


def batch_estimate_flow_statistics(
    points: Sequence[Tuple[str, str, float, float]],
    *,
    statistic_group_codes: Optional[Sequence[str]] = None,
    unit_system: int = 2,
    concurrency: int = 1,
    cache: Optional[StreamStatsCache] = None,
    validate_regions: bool = True,
    timeout: float = 60.0,
    select_by_location: bool = True,
    include_unlocated: bool = False,
) -> Tuple[Dict[str, List[RegionFlowEstimates]], Dict[str, Dict[str, str]], Dict[str, str]]:
    """Delineate, then estimate flow statistics, for many pour points.

    Composes :func:`batch_get_characteristics` (Phase 1) with
    :func:`estimate_flow_statistics` for each point that successfully delineates.
    Mirrors FR-7's batch shape at two levels: one point's delineation failure never
    blocks another point's estimate, and within an otherwise-successful point, one
    bad regressionRegion never blocks its siblings.

    Parameters
    ----------
    points : sequence of (point_id, region, lat, lon)
    statistic_group_codes, unit_system
        Passed through to :func:`estimate_flow_statistics` for every point.
    concurrency, cache, validate_regions, timeout
        Passed through to :func:`batch_get_characteristics`.
    select_by_location : bool
        Pass each point's own watershed polygon to :func:`estimate_flow_statistics`
        so NSS regions are chosen by location (default). ``False`` restores the old
        every-in-range-region behaviour.
    include_unlocated : bool
        Passed through to :func:`estimate_flow_statistics`.

    Returns
    -------
    tuple
        ``(estimates, skipped, errors)`` -- ``estimates`` maps ``point_id`` to a
        list of :class:`RegionFlowEstimates` for points that both delineated and
        estimated successfully; ``skipped`` maps ``point_id`` to
        :func:`estimate_flow_statistics`'s own per-region skip dict, for an
        otherwise-successful point with some regions excluded; ``errors`` maps
        ``point_id`` to a failure reason, for points that never delineated
        (:func:`batch_get_characteristics`'s own errors) or whose estimate call
        itself failed outright.
    """
    characteristics_by_point, errors = batch_get_characteristics(
        points,
        concurrency=concurrency,
        cache=cache,
        validate_regions=validate_regions,
        timeout=timeout,
    )

    region_by_point = {point_id: region for point_id, region, _, _ in points}
    estimates: Dict[str, List[RegionFlowEstimates]] = {}
    skipped: Dict[str, Dict[str, str]] = {}

    for point_id, watershed in characteristics_by_point.items():
        try:
            region_estimates, point_skipped = estimate_flow_statistics(
                region_by_point[point_id],
                watershed.characteristics,
                statistic_group_codes=statistic_group_codes,
                unit_system=unit_system,
                timeout=timeout,
                watershed_polygon=watershed.polygon_geojson if select_by_location else None,
                include_unlocated=include_unlocated,
            )
            estimates[point_id] = region_estimates
            if point_skipped:
                skipped[point_id] = point_skipped
        except Exception as e:
            errors[point_id] = str(e)

    return estimates, skipped, errors
