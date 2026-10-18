"""Oregon peak-flow regression update, task 3: basin characteristics.

NOT a USGS study. Oregon regression update plan
(``docs/superpowers/plans/2026-10-12-oregon-regression-update.md``), task 3.
Oregon StreamStats was re-checked on 2026-10-13 and still computes only 9-13
characteristics: none of PRECIP, ELEV, I24H2Y or JANMAXT2K in the west, and a
PRECIP of 0.013 in. for the John Day at Service Creek. So only its drainage
area and basin polygon are used, and every other characteristic is computed
here from one national source each, over that polygon.

The polygon is StreamStats' when its area is within 10 % of the NWIS drainage
(or contributing) area, else NLDI's basin for the gage's indexed NHDPlus
flowline under the same test, else none (``basin_error``). StreamStats snaps the
gage to its own stream grid and lands on another stream at roughly a quarter
of these gages.

====================  ==========================================================
``DRNAREA``           the accepted polygon's area: StreamStats' DRNAREA, or the
                      NLDI polygon's geodesic area (mi2)
``PRECIP``            PRISM 1991-2020 normals, annual precipitation, 800 m (in.)
``JANMAXT``           PRISM 1991-2020 normals, January mean maximum temperature,
                      800 m (deg F)
``I24H2Y``            NOAA Atlas 2, 2-year 24-hour precipitation, western U.S.
                      digitized grid (NOAA HDSC, 2016) (in.)
``ELEV``              USGS 3DEP 1/3 arc-second DEM, basin mean, from the 3DEP
                      ImageServer (ft)
``BSLOPD``            the same, through its ``Slope Degrees`` function (degrees)
``SOILKSAT``          STATSGO2 saturated hydraulic conductivity (in./hr): each
                      component's horizon ``ksat_r`` averaged by thickness over
                      its profile, then the components averaged by
                      ``comppct_r`` (those with no horizon data left out), as
                      Wolock's (1997) STATSGO PERMAVE; rasterized at 15
                      arc-seconds from the polygons and tables NRCS Soil Data
                      Access serves (``gsmmupolygon``, areasymbol ``US``)
``VOLCQP``            percent of the basin on volcanic rock no older than
                      Pliocene: USGS SGMC (Horton and others, 2017) units with
                      ``generalize`` "Igneous, volcanic" and an ``age_max`` in the
                      Pliocene or Quaternary; rasterized at 15 arc-seconds
``VOLCQP_OG``         as ``VOLCQP``, but classified in Oregon from the DOGAMI
                      Oregon Geologic Data Compilation, release 8 (OGDC-8;
                      Darin and others, 2025): units whose ``ThematicRockType``
                      is volcanic, vent/pyroclastic or volcaniclastic and whose
                      ``ThematicAge`` begins no earlier than the Pliocene. Outside
                      Oregon, which OGDC does not map, the SGMC mask is used.
``LAKEAREA``          percent of the basin under lakes, ponds and reservoirs
                      of at least 0.05 km2 (NHD high-resolution waterbodies,
                      FTYPE 390 and 436); rasterized at 5 arc-seconds
====================  ==========================================================

Each zonal mean is the mean of the raster cells whose centres fall inside the
basin polygon (holes excluded), or the cell under the polygon's centroid for a
basin smaller than one cell. Raster reading needs ``rasterio``
(``pip install -e ".[study]"``); nothing under ``flowfreq/`` imports it.

Stages, each cached under ``--cache`` and resumable:

- ``basins``: StreamStats delineation for every gage
  ``oregon_stations.csv`` keeps, in its own state's StreamStats region.
- ``characteristics``: the zonal statistics; writes
  ``data/regression_study/oregon_basin_characteristics.csv``.
- ``soils`` (before ``characteristics``): STATSGO2 polygons and
  component/horizon tables from Soil Data Access over the cached basins'
  extent, rasterized to ``<rasters>/statsgo2_ksat_15s.tif``.
- ``geology`` (before ``characteristics``): young volcanic SGMC units and NHD
  lakes over the same extent, from their ArcGIS feature services, rasterized
  to ``<rasters>/sgmc_volcqp_15s.tif`` and ``<rasters>/nhd_lakes_5s.tif``.
- ``ogdc`` (after ``geology``, before ``characteristics``): young volcanic
  OGDC-8 units from DOGAMI's OGDC map service, on the SGMC mask's grid, with
  the SGMC mask outside Oregon (TIGERweb state boundary), rasterized to
  ``<rasters>/ogdc_volcqp_15s.tif``.

Usage::

    python tools/build_oregon_basin_characteristics.py --cache <dir> basins
    python tools/build_oregon_basin_characteristics.py --cache <dir> --rasters <dir> soils
    python tools/build_oregon_basin_characteristics.py --cache <dir> --rasters <dir> geology
    python tools/build_oregon_basin_characteristics.py --cache <dir> --rasters <dir> ogdc
    python tools/build_oregon_basin_characteristics.py --cache <dir> --rasters <dir> characteristics
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
STATIONS = ROOT / "data" / "regression_study" / "oregon_stations.csv"
OUT = ROOT / "data" / "regression_study" / "oregon_basin_characteristics.csv"

#: name -> (path under --rasters, scale to the stored unit, offset after scaling).
RASTERS: Dict[str, Tuple[str, float, float]] = {
    # mm -> in.
    "PRECIP": ("prism_ppt_us_30s_2020_avg_30y/prism_ppt_us_30s_2020_avg_30y.tif", 1 / 25.4, 0.0),
    # deg C -> deg F
    "JANMAXT": (
        "prism_tmax_us_30s_202001_avg_30y/prism_tmax_us_30s_202001_avg_30y.tif",
        1.8,
        32.0,
    ),
    # in. * 100000 -> in.
    "I24H2Y": ("na2_westus_2yr24hr/na2_westus_2yr24hr.asc", 1e-5, 0.0),
    # um/s -> in./hr (the raster the ``soils`` stage writes)
    "SOILKSAT": ("statsgo2_ksat_15s.tif", 3.6 / 25.4, 0.0),
    # 0/1 cells -> percent (the ``geology`` stage)
    "VOLCQP": ("sgmc_volcqp_15s.tif", 100.0, 0.0),
    "LAKEAREA": ("nhd_lakes_5s.tif", 100.0, 0.0),
    # 0/1 cells -> percent (the ``ogdc`` stage)
    "VOLCQP_OG": ("ogdc_volcqp_15s.tif", 100.0, 0.0),
}
SOURCES = {
    "DRNAREA": "StreamStats delineation, or NLDI basin where StreamStats mis-snapped "
    "(basin_source; area within 10 % of NWIS either way), 2026-10-13",
    "PRECIP": "PRISM Climate Group, 1991-2020 normals, annual ppt, 30 arc-second",
    "JANMAXT": "PRISM Climate Group, 1991-2020 normals, January tmax, 30 arc-second",
    "I24H2Y": "NOAA Atlas 2 (Miller and others, 1973), western U.S. 2-yr 24-hr grid, "
    "NOAA HDSC digitization (15 arc-second, NAD27)",
    "ELEV": "USGS 3DEP 1/3 arc-second DEM via the 3DEPElevation ImageServer, "
    "computeStatisticsHistograms at 30 m",
    "BSLOPD": "USGS 3DEP via the 3DEPElevation ImageServer, Slope Degrees raster "
    "function, basin mean at 30 m",
    "SOILKSAT": "NRCS STATSGO2 (U.S. General Soil Map) via Soil Data Access, retrieved "
    "2026-10-13: horizon ksat_r thickness-weighted per component, components weighted "
    "by comppct_r, rasterized at 15 arc-seconds (in./hr)",
    "VOLCQP": "USGS State Geologic Map Compilation (Horton and others, 2017, DS 1052, "
    "doi:10.5066/F7WH2N65) via the energy.usgs.gov State_Geological_Map_Compilation "
    "feature service, retrieved 2026-10-14: percent of basin on units with generalize "
    "'Igneous, volcanic' and age_max Pliocene or Quaternary, 15 arc-second raster",
    "VOLCQP_OG": "Oregon: DOGAMI Oregon Geologic Data Compilation release 8 (OGDC-8; "
    "Darin and others, 2025) via the gis.dogami.oregon.gov Public/OGDC MapServer layer 2, "
    "retrieved 2026-10-18: percent of basin on units with ThematicRockType volcanic rocks, "
    "vent and pyroclastic rocks or volcaniclastic rocks and a ThematicAge beginning no "
    "earlier than the Pliocene; outside Oregon (TIGERweb boundary) the VOLCQP SGMC mask; "
    "15 arc-second raster",
    "LAKEAREA": "USGS NHD high-resolution Waterbody (hydro.nationalmap.gov nhd MapServer "
    "layer 12), retrieved 2026-10-14: percent of basin under FTYPE 390 (LakePond) or "
    "436 (Reservoir) with AREASQKM >= 0.05, 5 arc-second raster",
}
#: A basin is accepted when its area is within this fraction of the NWIS area.
AREA_TOL = 0.10
NLDI = "https://api.water.usgs.gov/nldi/linked-data/nwissite/USGS-{site}/basin"
DEP = (
    "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/"
    "computeStatisticsHistograms"
)

SDA = "https://sdmdataaccess.sc.egov.usda.gov/tabular/post.rest"
#: Douglas-Peucker tolerance (degrees) for the STATSGO2 polygons: ~200 m, against
#: a 1:250,000 source and a 15 arc-second (~460 m) raster.
SOIL_REDUCE = 0.002
SOIL_RES = 1.0 / 240.0
SOIL_NODATA = -9999.0

SGMC = (
    "https://energy.usgs.gov/arcgis/rest/services/Hosted/State_Geological_Map_Compilation/"
    "FeatureServer/1/query"
)
#: Volcanic units no older than Pliocene (their oldest age, ``age_max``).
SGMC_WHERE = (
    "generalize = 'Igneous, volcanic' AND "
    "(age_max LIKE '%Quaternary%' OR age_max LIKE '%Pliocene%')"
)
OGDC = "https://gis.dogami.oregon.gov/arcgis/rest/services/Public/OGDC/MapServer/2/query"
#: OGDC-8 ``ThematicRockType`` values counted as volcanic. Volcaniclastic rocks
#: include the Mazama pumice and ash; intrusive rocks (plugs, dikes) are left out.
OGDC_VOLCANIC = ("volcanic rocks", "vent and pyroclastic rocks", "volcaniclastic rocks")
#: OGDC-8 ``ThematicAge`` values that begin no earlier than the Pliocene. The
#: Tertiary/Quaternary volcanic units are late High Cascade lavas ("QT" labels);
#: Miocene/Pliocene and other ranges reaching into the Miocene (early High
#: Cascade, Deschutes Formation) are left out.
OGDC_YOUNG = (
    "Holocene",
    "Pleistocene",
    "Quaternary",
    "Pliocene",
    "Pliocene/Pleistocene",
    "Pliocene/Quaternary",
    "Tertiary/Quaternary",
)


def _sql_in(field: str, values: Sequence[str]) -> str:
    return f"{field} IN (" + ", ".join(f"'{v}'" for v in values) + ")"


OGDC_WHERE = (
    f"{_sql_in('ThematicRockType', OGDC_VOLCANIC)} AND {_sql_in('ThematicAge', OGDC_YOUNG)}"
)
STATE_BOUNDARY = (
    "https://tigerweb.geo.census.gov/arcgis/rest/services/TIGERweb/State_County/"
    "MapServer/0/query"
)
NHD_WATERBODY = "https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer/12/query"
NHD_WHERE = "FTYPE IN (390, 436) AND AREASQKM >= 0.05"
LAKE_RES = 1.0 / 720.0

log = logging.getLogger("or_basins")


def used_gages() -> pd.DataFrame:
    df = pd.read_csv(STATIONS, dtype={"site_no": str, "huc": str})
    return df[df["excluded"].isna()].reset_index(drop=True)


def delineate(cache: Path, concurrency: int) -> Dict[str, Dict[str, Any]]:
    """StreamStats basin per gage: ``{"DRNAREA": .., "geometry": ..}`` or ``{"error": ..}``."""
    from flowfreq.streamstats import StreamStatsCache, batch_get_characteristics

    d = cache / "basins"
    d.mkdir(parents=True, exist_ok=True)
    gages = used_gages()
    out: Dict[str, Dict[str, Any]] = {}
    todo = []
    for r in gages.itertuples():
        f = d / f"{r.site_no}.json"
        if f.exists():
            out[r.site_no] = json.loads(f.read_text(encoding="utf-8"))
        else:
            todo.append((r.site_no, r.state, float(r.lat), float(r.lon)))
    log.info("basins: %d cached, %d to delineate", len(out), len(todo))
    ss_cache = StreamStatsCache(cache / "streamstats_raw")
    for i in range(0, len(todo), 12):
        chunk = todo[i : i + 12]
        res, errs = batch_get_characteristics(
            chunk,
            concurrency=concurrency,
            cache=ss_cache,
            characteristic_codes=["DRNAREA"],
            validate_regions=False,
            timeout=180.0,
        )
        for site, _, _, _ in chunk:
            rec: Dict[str, Any]
            if site in res:
                w = res[site]
                c = w.characteristics.get("DRNAREA")
                rec = {
                    "DRNAREA": None if c is None else float(c.value),
                    "geometry": w.polygon_geojson["geometry"],
                }
            else:
                rec = {"error": errs.get(site, "unknown")[:300]}
            (d / f"{site}.json").write_text(json.dumps(rec), encoding="utf-8")
            out[site] = rec
        log.info("basins %d/%d", min(i + 12, len(todo)), len(todo))
    return out


def _polygons(geometry: Dict[str, Any]) -> List[List[List[List[float]]]]:
    if geometry["type"] == "Polygon":
        return [geometry["coordinates"]]
    return list(geometry["coordinates"])


def centroid(geometry: Dict[str, Any]) -> Tuple[float, float]:
    """Area-weighted (lon, lat) centroid of the outer rings."""
    a_sum = cx = cy = 0.0
    for poly in _polygons(geometry):
        ring = np.asarray(poly[0], dtype=float)
        x, y = ring[:, 0], ring[:, 1]
        cross = x[:-1] * y[1:] - x[1:] * y[:-1]
        cx += float(((x[:-1] + x[1:]) * cross).sum()) / 6.0
        cy += float(((y[:-1] + y[1:]) * cross).sum()) / 6.0
        a_sum += float(cross.sum()) / 2.0
    return cx / a_sum, cy / a_sum


def zonal_mean(src: Any, geometry: Dict[str, Any]) -> Tuple[float, int]:
    """Mean of the cells whose centres fall in ``geometry``, and how many there were.

    A basin smaller than one cell gets the cell under its centroid (count 0).
    ``src`` is an open rasterio dataset in geographic coordinates.
    """
    from rasterio.features import geometry_mask
    from rasterio.windows import from_bounds

    xs = np.concatenate([np.asarray(p[0], dtype=float)[:, 0] for p in _polygons(geometry)])
    ys = np.concatenate([np.asarray(p[0], dtype=float)[:, 1] for p in _polygons(geometry)])
    rx, ry = src.res
    win = (
        from_bounds(
            xs.min() - rx, ys.min() - ry, xs.max() + rx, ys.max() + ry, transform=src.transform
        )
        .round_offsets()
        .round_lengths()
    )
    data = src.read(1, window=win, masked=True).astype(float)
    inside = ~geometry_mask(
        [geometry], out_shape=data.shape, transform=src.window_transform(win), invert=False
    )
    vals = data[inside & ~np.ma.getmaskarray(data)]
    if vals.size:
        return float(vals.mean()), int(vals.size)
    lon, lat = centroid(geometry)
    v = next(src.sample([(lon, lat)], masked=True))[0]
    return (float("nan") if np.ma.is_masked(v) else float(v)), 0


def dem_stats(geometry: Dict[str, Any], timeout: float = 300.0) -> Dict[str, Optional[float]]:
    """Basin-mean 3DEP elevation (ft) and slope (degrees) from the ImageServer.

    ``computeStatisticsHistograms`` over the polygon at a 30 m request pixel,
    coarsened for large basins to keep under ~4 million pixels (the service
    returns no statistics above that; basins over ~7,000 mi2 at 30 m). The slope
    uses the service's own ``Slope Degrees`` raster function.
    """
    from flowfreq.streamstats import geojson_area_sq_mi

    rings = [ring for poly in _polygons(geometry) for ring in poly]
    area_m2 = geojson_area_sq_mi(geometry) * 2.59e6
    px = max(30.0, float(np.sqrt(area_m2 / 4e6)))
    out: Dict[str, Optional[float]] = {}
    for name, rule, scale in (
        ("ELEV", None, 3.28084),
        ("BSLOPD", {"rasterFunction": "Slope Degrees"}, 1.0),
    ):
        data = {
            "geometry": json.dumps({"rings": rings, "spatialReference": {"wkid": 4326}}),
            "geometryType": "esriGeometryPolygon",
            "pixelSize": json.dumps({"x": px, "y": px, "spatialReference": {"wkid": 3857}}),
            "f": "json",
        }
        if rule:
            data["renderingRule"] = json.dumps(rule)
        r = requests.post(DEP, data=data, timeout=timeout)
        r.raise_for_status()
        stats = r.json().get("statistics") or []
        mean = stats[0].get("mean") if stats else None
        out[name] = None if mean is None else float(mean) * scale
    return out


def nldi_basin(cache: Path, site: str) -> Dict[str, Any]:
    """NLDI's basin for the gage's indexed NHDPlus flowline (cached), or ``{"error": ..}``."""
    d = cache / "nldi"
    d.mkdir(exist_ok=True)
    f = d / f"{site}.json"
    if f.exists():
        return dict(json.loads(f.read_text(encoding="utf-8")))
    try:
        r = requests.get(NLDI.format(site=site), params={"simplified": "false"}, timeout=180)
        rec = (
            {"geometry": r.json()["features"][0]["geometry"]}
            if r.ok
            else {"error": f"HTTP {r.status_code}"}
        )
    except (requests.RequestException, ValueError, KeyError, IndexError) as exc:
        return {"error": f"{type(exc).__name__}: {exc}"[:200]}  # not cached: transient
    f.write_text(json.dumps(rec), encoding="utf-8")
    return rec


def _ratio(area: Optional[float], nwis: Sequence[Optional[float]]) -> float:
    """``area`` over the NWIS drainage (or contributing) area it is closest to."""
    refs = [float(x) for x in nwis if x is not None and pd.notna(x) and float(x) > 0]
    if area is None or not refs or not area > 0:
        return float("nan")
    return min((area / x for x in refs), key=lambda r: abs(np.log(r)))


def resolve_basin(
    cache: Path, gage: Any, streamstats: Dict[str, Any]
) -> Tuple[Optional[Dict[str, Any]], Optional[float], str, float]:
    """(geometry, DRNAREA, source, area ratio) for one gage; geometry None if neither fits.

    StreamStats snaps the gage location to its own stream grid and lands on the
    wrong stream at about a third of these gages (a tributary gage snapped onto
    the Klamath mainstem, 10,638 mi2 for 13). So a basin is accepted only when
    its area is within :data:`AREA_TOL` of the NWIS drainage or contributing
    area: StreamStats first, then NLDI's basin for the gage's indexed flowline.
    """
    from flowfreq.streamstats import geojson_area_sq_mi

    nwis = (gage.drainage_area, gage.contributing_area)
    if "geometry" in streamstats:
        ratio = _ratio(streamstats.get("DRNAREA"), nwis)
        if abs(ratio - 1) <= AREA_TOL:
            return streamstats["geometry"], streamstats.get("DRNAREA"), "StreamStats", ratio
    nl = nldi_basin(cache, gage.site_no)
    if "geometry" in nl:
        area = geojson_area_sq_mi(nl["geometry"])
        ratio = _ratio(area, nwis)
        if abs(ratio - 1) <= AREA_TOL:
            return nl["geometry"], area, "NLDI", ratio
    ss_ratio = _ratio(streamstats.get("DRNAREA"), nwis)
    why = f"no basin within {AREA_TOL:.0%} of NWIS area (StreamStats ratio {ss_ratio:.2f}"
    why += f"; NLDI {nl.get('error') or 'ratio %.2f' % _ratio(geojson_area_sq_mi(nl['geometry']), nwis)})"
    return None, None, why, ss_ratio


def sda_query(query: str, timeout: float = 600.0) -> List[Dict[str, Any]]:
    """Rows of a Soil Data Access query, as dicts keyed by column name."""
    r = requests.post(SDA, json={"query": query, "format": "JSON+COLUMNNAME"}, timeout=timeout)
    r.raise_for_status()
    table = r.json().get("Table") or []
    if not table:
        return []
    cols = table[0]
    return [dict(zip(cols, row)) for row in table[1:]]


def _wkt_polygon(body: str) -> List[List[List[float]]]:
    """Rings of one WKT polygon body, ``(x y, ...), (x y, ...)`` with or without outer parens."""
    rings = []
    for ring in re.split(r"\)\s*,\s*\(", body):
        ring = ring.replace("(", "").replace(")", "").strip()
        if ring:
            rings.append([[float(v) for v in xy.split()[:2]] for xy in ring.split(",")])
    return rings


def wkt_to_geojson(wkt: str) -> Dict[str, Any]:
    """GeoJSON for a WKT ``POLYGON`` or ``MULTIPOLYGON`` (SDA's ``STAsText``).

    Raises ``ValueError`` for any other geometry type.
    """
    m = re.match(r"\s*([A-Za-z]+)\s*(.*)$", wkt, re.S)
    if m is None:
        raise ValueError("not WKT")
    kind, body = m.group(1).upper(), m.group(2).strip()
    if kind == "POLYGON":
        return {"type": "Polygon", "coordinates": _wkt_polygon(body[1:-1])}
    if kind == "MULTIPOLYGON":
        parts = re.split(r"\)\)\s*,\s*\(\(", body[1:-1])
        return {"type": "MultiPolygon", "coordinates": [_wkt_polygon(p) for p in parts]}
    raise ValueError(f"unsupported WKT geometry {kind!r}")


def mapunit_ksat(rows: pd.DataFrame) -> pd.Series:
    """Map-unit Ksat (um/s) from component/horizon rows, indexed by ``mukey``.

    ``rows`` has ``mukey, cokey, comppct_r, hzdept_r, hzdepb_r, ksat_r``, one per
    horizon (a component without horizons has NaN depths). A component's Ksat
    is its horizons' ``ksat_r`` weighted by thickness; a map unit's is its
    components' weighted by ``comppct_r``, over the components that have one.
    A map unit with none (water, rock outcrop only) is absent from the result.
    """
    h = rows.dropna(subset=["hzdept_r", "hzdepb_r", "ksat_r"]).copy()
    h["thick"] = (h["hzdepb_r"] - h["hzdept_r"]).clip(lower=0.0)
    h = h[h["thick"] > 0]
    h["kt"] = h["ksat_r"] * h["thick"]
    comp = h.groupby(["mukey", "cokey"]).agg(kt=("kt", "sum"), t=("thick", "sum")).reset_index()
    comp["ksat"] = comp["kt"] / comp["t"]
    pct = rows.drop_duplicates("cokey").set_index("cokey")["comppct_r"]
    comp["w"] = comp["cokey"].map(pct).fillna(0.0).astype(float)
    comp = comp[comp["w"] > 0]
    comp["kw"] = comp["ksat"] * comp["w"]
    g = comp.groupby("mukey")
    return (g["kw"].sum() / g["w"].sum()).rename("ksat_um_s")


def _basin_bounds(cache: Path) -> Tuple[float, float, float, float]:
    """Bounding box (lon0, lat0, lon1, lat1) of every cached basin polygon, padded 0.05 deg."""
    xs: List[float] = []
    ys: List[float] = []
    for sub in ("basins", "nldi"):
        for f in (cache / sub).glob("*.json"):
            g = json.loads(f.read_text(encoding="utf-8")).get("geometry")
            for poly in _polygons(g) if g else []:
                ring = np.asarray(poly[0], dtype=float)
                xs += [float(ring[:, 0].min()), float(ring[:, 0].max())]
                ys += [float(ring[:, 1].min()), float(ring[:, 1].max())]
    return min(xs) - 0.05, min(ys) - 0.05, max(xs) + 0.05, max(ys) + 0.05


def soils(cache: Path, rasters: Path) -> Path:
    """Stage ``soils``: the STATSGO2 Ksat raster (um/s) over the basins' extent.

    Soil Data Access responses are cached under ``<cache>/statsgo2``, so a rerun
    only rasterizes.
    """
    import rasterio
    from rasterio.features import rasterize
    from rasterio.transform import from_origin

    d = cache / "statsgo2"
    d.mkdir(parents=True, exist_ok=True)
    x0, y0, x1, y1 = _basin_bounds(cache)
    bbox = f"POLYGON(({x0} {y0},{x1} {y0},{x1} {y1},{x0} {y1},{x0} {y0}))"
    where = f"mupolygongeo.STIntersects(geometry::STGeomFromText('{bbox}',4326))=1"
    f_keys = d / "mukeys.json"
    if not f_keys.exists():
        rows = sda_query(f"SELECT DISTINCT mukey FROM gsmmupolygon WHERE {where}")
        f_keys.write_text(json.dumps(sorted(int(r["mukey"]) for r in rows)))
    mukeys: List[int] = json.loads(f_keys.read_text())
    log.info("soils: %d STATSGO2 map units in %s", len(mukeys), (x0, y0, x1, y1))

    f_tab = d / "horizons.csv"
    if not f_tab.exists():
        parts = []
        for i in range(0, len(mukeys), 300):
            keys = ",".join(str(k) for k in mukeys[i : i + 300])
            q = (
                "SELECT c.mukey, c.cokey, c.comppct_r, h.hzdept_r, h.hzdepb_r, h.ksat_r "
                "FROM component c LEFT JOIN chorizon h ON h.cokey = c.cokey "
                f"WHERE c.mukey IN ({keys})"
            )
            parts.append(pd.DataFrame(sda_query(q)))
        pd.concat(parts).to_csv(f_tab, index=False)
    tab = pd.read_csv(f_tab)
    for c in ("comppct_r", "hzdept_r", "hzdepb_r", "ksat_r"):
        tab[c] = pd.to_numeric(tab[c], errors="coerce")
    ksat = mapunit_ksat(tab)
    log.info("soils: Ksat for %d of %d map units", len(ksat), len(mukeys))

    shapes: List[Tuple[Dict[str, Any], int]] = []
    for i in range(0, len(mukeys), 40):
        f = d / f"poly_{i:05d}.json"
        if not f.exists():
            keys = ",".join(str(k) for k in mukeys[i : i + 40])
            rows = sda_query(
                f"SELECT mukey, mupolygongeo.Reduce({SOIL_REDUCE}).STAsText() AS wkt "
                f"FROM gsmmupolygon WHERE mukey IN ({keys}) AND {where}"
            )
            f.write_text(json.dumps(rows))
            log.info("soils: polygons %d/%d", min(i + 40, len(mukeys)), len(mukeys))
        for r in json.loads(f.read_text()):
            try:
                shapes.append((wkt_to_geojson(r["wkt"]), int(r["mukey"])))
            except ValueError as exc:
                log.warning("mukey %s: %s", r["mukey"], exc)

    width = int(np.ceil((x1 - x0) / SOIL_RES))
    height = int(np.ceil((y1 - y0) / SOIL_RES))
    transform = from_origin(x0, y1, SOIL_RES, SOIL_RES)
    mu = rasterize(shapes, out_shape=(height, width), transform=transform, fill=0, dtype="int32")
    uniq, inv = np.unique(mu, return_inverse=True)
    lut = ksat.to_dict()
    vals = np.array([lut.get(int(k), SOIL_NODATA) for k in uniq], dtype="float32")
    out = vals[inv].reshape(mu.shape)
    path = rasters / RASTERS["SOILKSAT"][0]
    rasters.mkdir(parents=True, exist_ok=True)
    profile = dict(driver="GTiff", width=width, height=height, count=1, dtype="float32",
                   crs="EPSG:4326", transform=transform, nodata=SOIL_NODATA, compress="deflate")  # fmt: skip
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(out, 1)
    covered = 100.0 * float((out != SOIL_NODATA).mean())
    log.info("soils: wrote %s (%d x %d), %.1f %% of cells with Ksat", path, width, height, covered)
    return path


def arcgis_features(
    url: str,
    where: str,
    bounds: Tuple[float, float, float, float],
    order_by: str,
    cache: Path,
    reduce_deg: float,
    page: int = 2000,
) -> List[Dict[str, Any]]:
    """GeoJSON features of an ArcGIS feature query over ``bounds``, paged and cached.

    Each page is cached as ``<cache>/page_<offset>.json``; geometry comes back in
    EPSG:4326, generalized by ``reduce_deg`` degrees.
    """
    cache.mkdir(parents=True, exist_ok=True)
    x0, y0, x1, y1 = bounds
    env = {"xmin": x0, "ymin": y0, "xmax": x1, "ymax": y1, "spatialReference": {"wkid": 4326}}
    feats: List[Dict[str, Any]] = []
    offset = 0
    while True:
        f = cache / f"page_{offset:06d}.json"
        if not f.exists():
            params = {
                "where": where,
                "geometry": json.dumps(env),
                "geometryType": "esriGeometryEnvelope",
                "inSR": 4326,
                "spatialRel": "esriSpatialRelIntersects",
                "outFields": order_by,
                "outSR": 4326,
                "maxAllowableOffset": reduce_deg,
                "orderByFields": order_by,
                "resultOffset": offset,
                "resultRecordCount": page,
                "f": "geojson",
            }
            r = requests.get(url, params=params, timeout=600)
            r.raise_for_status()
            d = r.json()
            if "features" not in d:
                raise RuntimeError(f"{url}: {str(d)[:300]}")
            f.write_text(json.dumps(d["features"]), encoding="utf-8")
            log.info("%s: page at %d, %d features", cache.name, offset, len(d["features"]))
        got = json.loads(f.read_text(encoding="utf-8"))
        feats += got
        if len(got) < page:
            return feats
        offset += page


def _write_mask(
    path: Path, feats: List[Dict[str, Any]], bounds: Tuple[float, float, float, float], res: float
) -> float:
    """Rasterize ``feats`` to a 0/1 uint8 GeoTIFF over ``bounds``; returns the fraction set."""
    import rasterio
    from rasterio.features import rasterize
    from rasterio.transform import from_origin

    x0, y0, x1, y1 = bounds
    width = int(np.ceil((x1 - x0) / res))
    height = int(np.ceil((y1 - y0) / res))
    transform = from_origin(x0, y1, res, res)
    shapes = [(ft["geometry"], 1) for ft in feats if ft.get("geometry")]
    mask = rasterize(shapes, out_shape=(height, width), transform=transform, fill=0, dtype="uint8")
    profile = dict(driver="GTiff", width=width, height=height, count=1, dtype="uint8",
                   crs="EPSG:4326", transform=transform, nodata=255, compress="deflate")  # fmt: skip
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(mask, 1)
    return float(mask.mean())


def geology(cache: Path, rasters: Path) -> None:
    """Stage ``geology``: SGMC young-volcanic and NHD lake masks over the basins' extent.

    Both are presence masks (1 inside a selected polygon, 0 elsewhere), so a
    basin's zonal mean is the fraction of its area covered. SGMC's units are not
    reconciled across state lines, so a unit's age and class follow the source
    state map.
    """
    bounds = _basin_bounds(cache)
    rasters.mkdir(parents=True, exist_ok=True)
    volc = arcgis_features(SGMC, SGMC_WHERE, bounds, "fid", cache / "sgmc", SOIL_REDUCE / 2)
    frac = _write_mask(rasters / RASTERS["VOLCQP"][0], volc, bounds, SOIL_RES)
    log.info(
        "geology: %d SGMC young volcanic polygons, %.1f %% of the extent", len(volc), 100 * frac
    )
    lakes = arcgis_features(NHD_WATERBODY, NHD_WHERE, bounds, "OBJECTID", cache / "nhd", 0.0002)
    frac = _write_mask(rasters / RASTERS["LAKEAREA"][0], lakes, bounds, LAKE_RES)
    log.info("geology: %d NHD lakes and reservoirs, %.2f %% of the extent", len(lakes), 100 * frac)


def oregon_polygon(cache: Path) -> Dict[str, Any]:
    """Oregon's boundary (GeoJSON, EPSG:4326, generalized to 0.001 degree), cached."""
    f = cache / "oregon_state_0p001.json"
    if not f.exists():
        r = requests.get(
            STATE_BOUNDARY,
            params={
                "where": "STUSAB='OR'",
                "outFields": "STUSAB",
                "outSR": 4326,
                "maxAllowableOffset": 0.001,
                "f": "geojson",
            },
            timeout=300,
        )
        r.raise_for_status()
        f.write_text(json.dumps(r.json()["features"][0]["geometry"]), encoding="utf-8")
    return dict(json.loads(f.read_text(encoding="utf-8")))


def merge_state_mask(inside: np.ndarray, state: np.ndarray, fallback: np.ndarray) -> np.ndarray:
    """``state`` where ``inside`` is set, ``fallback`` elsewhere (all 0/1 arrays, one grid).

    OGDC maps Oregon only; a basin reaching into a neighbouring state takes
    that state's SGMC classification there instead of zero.
    """
    if not inside.shape == state.shape == fallback.shape:
        raise ValueError(f"grids differ: {inside.shape}, {state.shape}, {fallback.shape}")
    return np.where(inside.astype(bool), state, fallback).astype("uint8")


def ogdc(cache: Path, rasters: Path) -> None:
    """Stage ``ogdc``: young volcanic OGDC-8 units in Oregon, SGMC's outside it.

    Written on the grid of the ``geology`` stage's SGMC mask, which must exist.
    """
    import rasterio
    from rasterio.features import rasterize

    sgmc_path = rasters / RASTERS["VOLCQP"][0]
    if not sgmc_path.exists():
        raise SystemExit(f"{sgmc_path} is missing: run the geology stage first")
    bounds = _basin_bounds(cache)
    feats = arcgis_features(OGDC, OGDC_WHERE, bounds, "OBJECTID", cache / "ogdc8", SOIL_REDUCE / 2)
    with rasterio.open(sgmc_path) as src:
        fallback = src.read(1)
        profile = {
            k: v for k, v in src.profile.items() if k not in ("blockxsize", "blockysize", "tiled")
        }
        transform = src.transform
    shapes = [(ft["geometry"], 1) for ft in feats if ft.get("geometry")]
    young = rasterize(shapes, out_shape=fallback.shape, transform=transform, fill=0, dtype="uint8")
    inside = rasterize(
        [(oregon_polygon(cache), 1)],
        out_shape=fallback.shape,
        transform=transform,
        fill=0,
        dtype="uint8",
    )
    mask = merge_state_mask(inside, young, fallback)
    with rasterio.open(rasters / RASTERS["VOLCQP_OG"][0], "w", **profile) as dst:
        dst.write(mask, 1)
    log.info(
        "ogdc: %d OGDC-8 young volcanic polygons; %.1f %% of the extent (SGMC %.1f %%)",
        len(feats),
        100 * float(mask.mean()),
        100 * float(fallback.mean()),
    )


def characteristics(cache: Path, rasters: Path) -> pd.DataFrame:
    import rasterio

    gages = used_gages()
    basins = {}
    for site in gages["site_no"]:
        f = cache / "basins" / f"{site}.json"
        basins[site] = json.loads(f.read_text(encoding="utf-8")) if f.exists() else {}
    rows: List[Dict[str, Any]] = []
    srcs = {k: rasterio.open(rasters / v[0]) for k, v in RASTERS.items()}
    dem_dir = cache / "dem"
    dem_dir.mkdir(exist_ok=True)
    try:
        for i, gage in enumerate(gages.itertuples()):
            site = gage.site_no
            row: Dict[str, Any] = {"site_no": site}
            g, area, source, ratio = resolve_basin(cache, gage, basins[site])
            row["area_ratio"] = ratio
            if g is None:
                row["basin_error"] = source
                rows.append(row)
                continue
            lon, lat = centroid(g)
            row.update(DRNAREA=area, basin_source=source, lat_cent=lat, lon_cent=lon)
            for name, (_, scale, offset) in RASTERS.items():
                v, n = zonal_mean(srcs[name], g)
                row[name] = v * scale + offset
                row[f"{name}_cells"] = n
            f = dem_dir / f"{site}_{source}.json"
            if not f.exists():
                try:
                    f.write_text(json.dumps(dem_stats(g)))
                except requests.RequestException as exc:
                    log.warning("%s: 3DEP failed: %s", site, exc)
            if f.exists():
                row.update(json.loads(f.read_text()))
            rows.append(row)
            if (i + 1) % 50 == 0:
                log.info("characteristics %d/%d", i + 1, len(gages))
    finally:
        for s in srcs.values():
            s.close()
    df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False, lineterminator="\n", float_format="%.6g")
    OUT.with_suffix(".sources.json").write_text(json.dumps(SOURCES, indent=1) + "\n")
    n_basin = int(df["DRNAREA"].notna().sum()) if "DRNAREA" in df else 0
    log.info("wrote %s: %d gages, %d with a basin", OUT.name, len(df), n_basin)
    return df


def main(argv: Sequence[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("stage", choices=["basins", "soils", "geology", "ogdc", "characteristics"])
    ap.add_argument("--cache", type=Path, default=ROOT / ".cache" / "or_basins")
    ap.add_argument("--rasters", type=Path, default=None)
    ap.add_argument("--concurrency", type=int, default=3)
    args = ap.parse_args(list(argv))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    logging.getLogger("flowfreq").setLevel(logging.ERROR)
    args.cache.mkdir(parents=True, exist_ok=True)
    if args.stage == "basins":
        out = delineate(args.cache, args.concurrency)
        bad: List[str] = [s for s, r in out.items() if "error" in r]
        log.info("basins: %d delineated, %d failed", len(out) - len(bad), len(bad))
        return 0
    if args.rasters is None:
        raise SystemExit(f"{args.stage} needs --rasters")
    if args.stage == "geology":
        geology(args.cache, args.rasters)
        return 0
    if args.stage == "ogdc":
        ogdc(args.cache, args.rasters)
        return 0
    if args.stage == "soils":
        soils(args.cache, args.rasters)
        return 0
    characteristics(args.cache, args.rasters)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
