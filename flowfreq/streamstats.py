"""
flowfreq.streamstats - USGS StreamStats watershed delineation and basin characteristics

Phase 1 (``docs/STREAMSTATS_MODULE_DESIGN.md``): given a pour point, snap it to the
stream network, delineate its watershed, and return basin characteristics (drainage
area, precipitation, canopy, ...) as regression predictors.

Phase 2 (``docs/STREAMSTATS_NSS_ADDENDUM.md``): feed those characteristics into NSS
(National Streamflow Statistics) to compute actual flow-statistic estimates (peak-flow,
low-flow, ...), each carrying its own regression equation and a resolved citation.
**Region selection across NSS's regressionRegions is not automatic** -- see
:func:`estimate_flow_statistics`'s docstring and the addendum S4 before assuming a
single "the" answer for a point; StreamStats itself does not resolve this without a
watershed polygon, which Phase 1 does not produce.

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

**No watershed polygon, confirmed live 2026-09-11.** An earlier version of this module
also called ``ss-delineate/v1/delineate/features/{region}`` (borrowed from the
unverified ff-idea02 PDF/py transcript, the same source whose ``ss-hydro`` endpoint
guess was already known wrong) to obtain and FR-3-validate the ``globalwatershed``
polygon. Called live with the region's real, snapped ``lat``/``lon``, it returns HTTP
200 with a single ``Point`` feature of zero ``Shape_Area``/``Shape_Leng`` -- an echo of
the snapped pour point, not a delineated basin -- so it was never the right call and has
been removed from the pipeline. This module now follows exactly the three-call protocol
docs/STREAMSTATS_MODULE_DESIGN.md S3 verified end to end (snap, ``delineate/sshydro``,
``ss-hydro``), none of which returns a polygon; :attr:`WatershedCharacteristics.
polygon_geojson` is always ``None``. FR-3's ``WarningMsg`` check is done by
:func:`_find_warning_msg`, a recursive scan of the ``sshydro`` response, rather than a
polygon-degeneracy check, since no polygon is available to check. Obtaining the actual
watershed geometry remains an open question -- worth its own live-verification pass
before attempting again, per TODO.md.
"""

from __future__ import annotations

import json
import logging
import math
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
        Reserved for the watershed's GeoJSON geometry; always ``None`` in this version.
        The verified delineation protocol (:func:`delineate_and_get_characteristics`)
        has no call that returns the polygon -- see the module docstring. Would be a
        plain dict, not a shapely geometry, if populated: NFR-6 keeps this module free
        of a hard GIS dependency.
    characteristics : dict of str to Characteristic
        Keyed by StreamStats parameter code (e.g. ``DRNAREA``, ``PRECPRIS10``).
    provenance : Provenance
        How and when this result was obtained.
    """

    region: str
    snap: SnapResult
    polygon_geojson: Optional[Dict[str, Any]]
    characteristics: Dict[str, Characteristic] = field(default_factory=dict)
    provenance: Optional[Provenance] = None

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
    """

    def __init__(self, path: Optional[Path] = None) -> None:
        self._path = path or _default_cache_path()
        self._data: Dict[str, Dict[str, Any]] = {}
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
        return self._data.get(key)

    def set(self, key: str, value: Dict[str, Any]) -> None:
        self._data[key] = value
        self._path.parent.mkdir(parents=True, exist_ok=True)
        self._path.write_text(json.dumps(self._data, indent=2), encoding="utf-8")


def _cache_key(region: str, lat: float, lon: float, bc_labels: str) -> str:
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

    The design doc's own live testing (S4) found this key on an unsnappable point's
    delineation response; the exact nesting within the ``sshydro`` chaining variant's
    response was not independently re-confirmed (see the module docstring), so this
    scans the whole structure defensively rather than assuming one exact path.
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


def _parse_characteristics(hydro_data: Any) -> Dict[str, Characteristic]:
    """Parse the ss-hydro basin-characteristics response (design doc S3): a list of
    ``{name, description, code, unit, value, msg}`` objects, optionally wrapped in a
    ``{"parameters": [...]}`` envelope.
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
    for item in items:
        try:
            code = str(item["code"])
            value = float(item["value"])
        except (KeyError, TypeError, ValueError) as exc:
            raise StreamStatsResponseError(
                f"ss-hydro characteristic entry missing code/value: {item!r}"
            ) from exc
        characteristics[code] = Characteristic(
            code=code,
            name=str(item.get("name", "")),
            description=str(item.get("description", "")),
            value=value,
            unit=str(item.get("unit", "")),
            msg=str(item.get("msg") or ""),
        )
    return characteristics


def delineate_and_get_characteristics(
    region: str,
    lat: float,
    lon: float,
    *,
    characteristic_codes: Optional[Sequence[str]] = None,
    cache: Optional[StreamStatsCache] = None,
    default_server: str = DEFAULT_SERVER,
    timeout: float = 60.0,
) -> WatershedCharacteristics:
    """Snap, delineate, and compute basin characteristics for one pour point.

    Implements exactly the protocol design doc S3 verified live: a ``pourpoint`` snap
    (FR-1), ``ss-delineate``'s ``delineate/sshydro`` call to obtain the request body
    ``ss-hydro`` needs (validated for a ``WarningMsg`` per FR-3), and the ``ss-hydro``
    POST itself (FR-5/FR-6). Both ``ss-delineate``/``ss-hydro`` calls are pinned to the
    same server (design doc S3 -- ``ss-hydro`` reads temporary files ``ss-delineate``
    left on that specific node). No watershed polygon is returned -- see the module
    docstring for why.

    Parameters
    ----------
    region : str
        StreamStats region code. Required, never inferred from the coordinate (FR-8).
    lat, lon : float
        Decimal-degree pour point, WGS84.
    characteristic_codes : sequence of str, optional
        StreamStats characteristic codes to request (resolves design doc S9.4).
        Defaults to all (``bcLabels=*``). Every request explicitly states this
        parameter rather than relying on a service default (FR-4).
    cache : StreamStatsCache, optional
        When given, a cache hit skips all network calls; a miss populates it.
    default_server : str
        Server for the first ``ss-delineate`` call, before the service's own
        ``usgswim-hostname`` header is known.
    timeout : float
        Per-request timeout in seconds. StreamStats delineation and characteristics
        calls have been observed to take 6-11 s each under good conditions (NFR-3).

    Returns
    -------
    WatershedCharacteristics

    Raises
    ------
    UnsnappablePointError, UnsupportedRegionError, DegenerateDelineationError,
    StreamStatsResponseError, StreamStatsTransportError
        See each exception's docstring; see also S7 of the design doc.
    """
    bc_labels = ",".join(characteristic_codes) if characteristic_codes else "*"
    cache_key = _cache_key(region, lat, lon, bc_labels)

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

    # No watershed polygon is available through this protocol: the one endpoint that
    # looked plausible for it (`ss-delineate/v1/delineate/features/{region}`, borrowed
    # from the unverified PDF/py script) was tried live and returns an unrelated,
    # zero-area Point feature echoing the snapped pour point, not a basin polygon --
    # not the design doc's own literally-verified protocol, which never called it. See
    # the module docstring.
    polygon_geojson: Optional[Dict[str, Any]] = None

    hydro_url = (
        f"https://{server_used}.{GENERIC_HOST}"
        "/ss-hydro/v1/basin-characteristics/calculate-using-ssdelineate/"
    )
    hydro_params = {
        "region": region,
        "lat": snap.snapped_lat,
        "lon": snap.snapped_lon,
        "bcLabels": bc_labels,
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

    characteristics = _parse_characteristics(hydro_data)

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
        NSS's own "ASEp" (average standard error of prediction), confusingly
        labelled ``errors`` in the raw response (a prediction-error statistic, not
        a failure indicator). ``None`` when NSS did not return one -- observed
        live for a wildly out-of-range input (addendum S3), so its absence is a
        soft signal worth noticing even though it is not a hard failure marker.
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
    for a point. See :func:`estimate_flow_statistics`'s docstring on region
    selection before treating one of these as authoritative for a given location.
    """

    region_code: str
    region_name: str
    statistic_group_code: str
    statistic_group_name: str
    estimates: Dict[str, FlowStatisticEstimate] = field(default_factory=dict)
    citation: Optional[RegressionCitation] = None

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
        )


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


def estimate_flow_statistics(
    region: str,
    characteristics: Dict[str, Characteristic],
    statistic_group_codes: Optional[Sequence[str]] = None,
    *,
    unit_system: int = 2,
    timeout: float = 60.0,
) -> Tuple[List[RegionFlowEstimates], Dict[str, str]]:
    """Estimate NSS regression flow statistics from real basin characteristics.

    Phase 2 (``docs/STREAMSTATS_NSS_ADDENDUM.md``). Given characteristics Phase 1's
    :func:`delineate_and_get_characteristics` already produced, computes every
    regression equation NSS defines for the region -- across every regressionRegion
    within every requested statistic group -- skipping (never silently computing)
    any region whose required parameters are missing or fall outside that region's
    own declared valid range (addendum S3).

    **Region selection is not automatic.** NSS defines multiple, independently
    calibrated regressionRegions per statistic group within a state (WA Peak-Flow
    has four), and no available call filters them by location: Phase 1 provides no
    watershed polygon, and confirmed live, NSS's own ``ByLocation`` call does not
    filter on a bare point either. Every geographically-plausible region computes a
    result with no error as long as its parameters are merely in numeric range --
    picking the geographically-correct one among the returned regions is the
    caller's responsibility, the same way FR-8 makes the StreamStats region itself
    the caller's responsibility one level up.

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

    Returns
    -------
    tuple
        ``(region_estimates, skipped)`` -- ``region_estimates`` is a list of
        :class:`RegionFlowEstimates`, one per (statistic group, regressionRegion)
        pair that had every required parameter in range; ``skipped`` maps
        ``"{statistic_group_code}:{region_code}"`` to why it was excluded. Mirrors
        FR-7's batch shape one level deeper: one bad regression region never
        excludes its siblings.

    Raises
    ------
    StreamStatsResponseError, StreamStatsTransportError
        A malformed response, or a transport failure after retrying.
    """
    if statistic_group_codes is None:
        statistic_group_codes = [g["code"] for g in list_statistic_groups(region, timeout=timeout)]

    id_to_code = {g["id"]: g["code"] for g in list_statistic_groups(timeout=timeout)}

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
    for scenario in templates:
        group_code = id_to_code.get(scenario.get("statisticGroupID"), "")
        filled, region_skips = _fill_and_validate_regions(scenario, characteristics)
        for region_code, reason in region_skips.items():
            skipped[f"{group_code}:{region_code}"] = reason
        if filled["regressionRegions"]:
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
        for rr in scenario.get("regressionRegions", []):
            estimates: Dict[str, FlowStatisticEstimate] = {}
            for stat in rr.get("results", []) or []:
                errors = stat.get("errors") or []
                sep = next((e.get("value") for e in errors if e.get("code") == "ASEp"), None)
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
                    interval_lower=bounds.get("lower"),
                    interval_upper=bounds.get("upper"),
                )
            region_id = rr.get("id")
            citation_id = rr.get("citationID")
            if region_id is not None:
                region_ids.append(region_id)
            parsed.append(
                (
                    citation_id,
                    RegionFlowEstimates(
                        region_code=str(rr.get("code", "")),
                        region_name=str(rr.get("name", "")),
                        statistic_group_code=group_code,
                        statistic_group_name=group_name,
                        estimates=estimates,
                    ),
                )
            )

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
            )
            estimates[point_id] = region_estimates
            if point_skipped:
                skipped[point_id] = point_skipped
        except Exception as e:
            errors[point_id] = str(e)

    return estimates, skipped, errors
