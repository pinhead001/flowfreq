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

Usage::

    python tools/build_oregon_basin_characteristics.py --cache <dir> basins
    python tools/build_oregon_basin_characteristics.py --cache <dir> --rasters <dir> characteristics
"""

from __future__ import annotations

import argparse
import json
import logging
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
}
#: A basin is accepted when its area is within this fraction of the NWIS area.
AREA_TOL = 0.10
NLDI = "https://api.water.usgs.gov/nldi/linked-data/nwissite/USGS-{site}/basin"
DEP = (
    "https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/"
    "computeStatisticsHistograms"
)

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
    ap.add_argument("stage", choices=["basins", "characteristics"])
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
        raise SystemExit("characteristics needs --rasters")
    characteristics(args.cache, args.rasters)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
