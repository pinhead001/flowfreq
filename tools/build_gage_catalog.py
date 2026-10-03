"""Build ``flowfreq/data/gage_catalog.csv.gz`` from the USGS Water Data OGC API.

Every USGS peak-flow site, active or inactive, with at least ``--min-years``
(default 10) water years of discharge peaks. Issue #33, roadmap section 1.1.

How (live-verified 2026-09-30):

1. ``peaks`` collection, one bulk query per state: ``state_code=<FIPS>``,
   ``parameter_code=00060``, ``properties=monitoring_location_id,water_year,
   value,qualifier``, 50,000 features a page. Washington is 27,162 peaks in
   one page.
2. ``monitoring-locations``, one bulk query per state with
   ``site_type_code=ST``, for name, point, drainage area and HUC. The few
   peak sites of other types (canals, lakes, tidal streams) are fetched one by
   one.
3. Regulation class: :func:`flowfreq.regulation.classify_site` with the
   site's peak codes and the packaged GAGES-II screen (#32).
4. Regression region (WA, OR, ID, MT only): point-in-polygon of the gage
   location against the StreamStats peak-flow region polygons
   (``gis.streamstats.usgs.gov/.../nss/regions/MapServer``, layers 40, 31, 11,
   22). Blank elsewhere; see :mod:`flowfreq.catalog` for what the code means.

Every Water Data request goes through :func:`flowfreq.waterdata.request`,
which sends ``USGS_API_KEY`` if set and backs off on 429/503. Without a key
the API allows 1000 requests per hour per IP; a national build is about 120
bulk requests plus a few hundred single-site ones, paced by ``--pace``.
Responses are cached per state under ``--cache``, so an interrupted build
resumes where it stopped (``--refresh`` ignores the cache).

Usage::

    python tools/build_gage_catalog.py --all                  # national build
    python tools/build_gage_catalog.py --states WA OR ID MT    # merge these states
    python tools/build_gage_catalog.py --state WA --sites wa.txt   # listed sites only

``--sites`` keeps the old per-site mode (one peak and one location request per
site) with a selectable ``--backend`` (default ``waterdata-ogc``, previously
``nwis-legacy``). The output is never edited by hand; rerun the tool.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
import requests
from matplotlib.path import Path as MplPath

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from flowfreq import waterdata  # noqa: E402
from flowfreq.catalog import (  # noqa: E402
    CATALOG_PATH,
    COLUMNS,
    META_PATH,
    MIN_YEARS,
    build_rows,
    catalog_row,
    load_catalog,
    location_fields,
    peaks_by_site,
    sites_needing_location,
    write_catalog,
)
from flowfreq.peak_sources import get_backend  # noqa: E402
from flowfreq.regional_skew import STATE_FIPS  # noqa: E402

logger = logging.getLogger("build_gage_catalog")

BASE = waterdata.WATERDATA_BASE_URL
POSTAL_FIPS = {v: k for k, v in STATE_FIPS.items()}

SS_REGIONS = "https://gis.streamstats.usgs.gov/arcgis/rest/services/nss/regions/MapServer"
#: Wave 1 peak-flow region polygons: state -> (layer id, grid_name values).
#: WA SIR 2016-5118 (nss_peak; gc0 "National Urban" excluded), OR SIR 2005-5116
#: west (or_regions) and OWRD OFR SW 06-001 east (regeastor), ID SIR 2016-5083
#: (reg6n3), MT SIR 2015-5019F (pkregionid).
REGION_LAYERS: Dict[str, tuple] = {
    "WA": (40, ("nss_peak",)),
    "OR": (31, ("or_regions", "regeastor")),
    "ID": (11, ("reg6n3",)),
    "MT": (22, ("pkregionid",)),
}
#: GRIDCODE values that are not regression regions (helper/urban codes).
NON_REGION_CODES = frozenset({"GC0", "GC10001", "GC10003"})


class Paced:
    """Space Water Data API calls at least ``pace`` seconds apart and count them."""

    def __init__(self, pace: float) -> None:
        self.pace = pace
        self.calls = 0
        self._last = 0.0

    def get(self, url: str, params: Optional[Mapping[str, Any]]) -> Dict[str, Any]:
        wait = self._last + self.pace - time.monotonic()
        if wait > 0:
            time.sleep(wait)
        self._last = time.monotonic()
        self.calls += 1
        return dict(waterdata.request(url, params, 300).json())

    def all_features(self, url: str, params: Mapping[str, Any]) -> List[Dict[str, Any]]:
        feats: List[Dict[str, Any]] = []
        payload = self.get(url, params)
        seen = set()
        while True:
            feats.extend(payload.get("features") or [])
            nxt = next(
                (lk["href"] for lk in payload.get("links") or [] if lk.get("rel") == "next"),
                None,
            )
            if nxt is None:
                return feats
            if nxt in seen:
                raise RuntimeError(f"paging loop at {nxt}")
            seen.add(nxt)
            payload = self.get(nxt, None)


def _cached(path: Path, refresh: bool, fetch) -> Any:  # type: ignore[no-untyped-def]
    if path.exists() and not refresh:
        return json.loads(path.read_text(encoding="utf-8"))
    data = fetch()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data), encoding="utf-8")
    tmp.replace(path)  # atomic: a half-written cache file never looks complete
    return data


def state_peak_features(api: Paced, fips: str, cache: Path, refresh: bool) -> List[dict]:
    params = {
        "f": "json",
        "state_code": fips,
        "parameter_code": "00060",
        "properties": "monitoring_location_id,water_year,value,qualifier",
        "skipGeometry": "true",
        "limit": 50000,
    }
    return list(
        _cached(
            cache / f"peaks_{fips}.json",
            refresh,
            lambda: [
                {"properties": f["properties"]}
                for f in api.all_features(f"{BASE}/peaks/items", params)
            ],
        )
    )


def state_locations(api: Paced, fips: str, cache: Path, refresh: bool) -> List[dict]:
    params = {
        "f": "json",
        "state_code": fips,
        "site_type_code": "ST",
        "properties": "monitoring_location_name,state_code,hydrologic_unit_code,drainage_area",
        "limit": 50000,
    }
    return list(
        _cached(
            cache / f"locations_{fips}.json",
            refresh,
            lambda: api.all_features(f"{BASE}/monitoring-locations/items", params),
        )
    )


def single_location(api: Paced, site: str, cache: Path, refresh: bool) -> Optional[dict]:
    def fetch() -> Optional[dict]:
        try:
            return api.get(f"{BASE}/monitoring-locations/items/USGS-{site}", {"f": "json"})
        except requests.HTTPError as exc:
            if exc.response is not None and exc.response.status_code == 404:
                return None
            raise

    return _cached(cache / "sites" / f"{site}.json", refresh, fetch)  # type: ignore[no-any-return]


def region_polygons(state: str, cache: Path, refresh: bool) -> List[dict]:
    layer, grids = REGION_LAYERS[state]

    def fetch() -> List[dict]:
        out: List[dict] = []
        offset = 0
        where = " OR ".join(f"grid_name='{g}'" for g in grids)
        while True:
            r = requests.get(
                f"{SS_REGIONS}/{layer}/query",
                params={
                    "where": where,
                    "outFields": "grid_name,GRIDCODE",
                    "returnGeometry": "true",
                    "outSR": "4326",
                    # ~50 m generalization, as tools/build_regional_skew_hucs.py
                    # uses; without it Idaho's full-resolution rings fail.
                    "maxAllowableOffset": "0.0005",
                    "orderByFields": "OBJECTID",
                    "resultOffset": offset,
                    "resultRecordCount": 50,
                    "f": "json",
                },
                timeout=300,
            )
            r.raise_for_status()
            payload = r.json()
            if "error" in payload:  # ArcGIS reports failures with HTTP 200
                raise RuntimeError(f"StreamStats regions layer {layer}: {payload['error']}")
            feats = payload.get("features") or []
            out.extend(
                {"gridcode": f["attributes"]["GRIDCODE"], "rings": f["geometry"]["rings"]}
                for f in feats
            )
            if not feats or not payload.get("exceededTransferLimit"):
                return out
            offset += len(feats)

    return list(_cached(cache / f"regions_{state}.json", refresh, fetch))


def assign_regions(
    sites: Mapping[str, Mapping[str, Any]], polygons: Sequence[Mapping[str, Any]]
) -> Dict[str, str]:
    """Region code(s) of the polygon(s) containing each site's point (even-odd rule)."""
    ids = [s for s, loc in sites.items() if loc.get("latitude") is not None]
    if not ids:
        return {}
    pts = np.array([[sites[s]["longitude"], sites[s]["latitude"]] for s in ids], dtype=float)
    hits: Dict[str, List[str]] = {s: [] for s in ids}
    for poly in polygons:
        count = np.zeros(len(pts), dtype=int)
        for ring in poly["rings"]:
            if len(ring) >= 3:
                count += MplPath(np.asarray(ring, dtype=float)).contains_points(pts)
        inside = count % 2 == 1
        codes = [
            c.strip().upper()
            for c in str(poly["gridcode"]).split(",")
            if c.strip() and c.strip().upper() not in NON_REGION_CODES
        ]
        for s, hit in zip(ids, inside):
            if hit:
                hits[s].extend(c for c in codes if c not in hits[s])
    return {s: ";".join(sorted(c)) for s, c in hits.items() if c}


def build_state(api: Paced, state: str, cache: Path, refresh: bool, min_years: int) -> pd.DataFrame:
    fips = POSTAL_FIPS[state]
    peaks = peaks_by_site(state_peak_features(api, fips, cache, refresh))
    locs = {}
    for feat in state_locations(api, fips, cache, refresh):
        f = location_fields(feat)
        locs[f["site_no"]] = f
    for site in sites_needing_location(peaks, locs, min_years):
        feat = single_location(api, site, cache, refresh)
        if feat:
            locs[site] = location_fields(feat)
    regions: Dict[str, str] = {}
    if state in REGION_LAYERS:
        regions = assign_regions(
            {s: locs[s] for s in peaks if s in locs}, region_polygons(state, cache, refresh)
        )
    df = build_rows(peaks, locs, min_years=min_years, regions=regions, default_state=state)
    # A site is listed under the state its location record names; drop rows the
    # state's peak query returned for a site located elsewhere (it is built there).
    return df[df["state"] == state]


def bulk(args: argparse.Namespace, states: Sequence[str]) -> pd.DataFrame:
    api = Paced(args.pace)
    frames = []
    for st in states:
        t0 = time.monotonic()
        df = build_state(api, st, args.cache, args.refresh, args.min_years)
        logger.info(
            "%s: %d sites (%d API calls so far, %.0f s)",
            st,
            len(df),
            api.calls,
            time.monotonic() - t0,
        )
        frames.append(df)
    logger.info("Water Data API calls this run: %d", api.calls)
    frames = [f for f in frames if len(f)]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=list(COLUMNS))


def per_site(args: argparse.Namespace) -> pd.DataFrame:
    """The original per-site mode, for a hand-picked list."""
    backend = get_backend(args.backend)
    sites = [
        s.strip()
        for s in args.sites.read_text().splitlines()
        if s.strip() and not s.startswith("#")
    ]
    api = Paced(args.pace)
    rows, errors = [], {}
    for site in sites:
        try:
            peaks = backend.fetch_peaks(site)
            loc = location_fields(single_location(api, site, args.cache, args.refresh) or {})
            row = catalog_row(
                site,
                args.state,
                peaks,
                site_name=loc.get("site_name") or "",
                latitude=loc.get("latitude"),
                longitude=loc.get("longitude"),
                drainage_area_sqmi=loc.get("drainage_area_sqmi"),
                huc8=loc.get("huc8"),
            )
        except Exception as exc:  # one bad site never aborts the batch
            errors[site] = f"{type(exc).__name__}: {exc}"
            continue
        if row["n_peaks"] >= args.min_years:
            rows.append(row)
    for site, msg in errors.items():
        logger.warning("%s: %s", site, msg)
    return pd.DataFrame(rows, columns=list(COLUMNS))


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--all", action="store_true", help="every state and territory")
    mode.add_argument("--states", nargs="+", help="postal codes to (re)build and merge")
    mode.add_argument("--sites", type=Path, help="per-site mode: file of site numbers")
    ap.add_argument("--state", help="postal code for --sites mode")
    ap.add_argument("--min-years", type=int, default=MIN_YEARS)
    ap.add_argument("--backend", default="waterdata-ogc", help="--sites mode peak backend")
    ap.add_argument("--out", type=Path, default=CATALOG_PATH)
    ap.add_argument("--meta", type=Path, default=META_PATH)
    ap.add_argument("--cache", type=Path, default=REPO_ROOT / ".cache" / "gage_catalog")
    ap.add_argument("--refresh", action="store_true", help="ignore cached responses")
    ap.add_argument("--pace", type=float, default=1.0, help="seconds between API calls")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    if args.sites:
        if not args.state:
            ap.error("--sites needs --state")
        new = per_site(args)
        built = sorted(set(new["state"]))
    else:
        built = sorted(STATE_FIPS.values()) if args.all else [s.upper() for s in args.states]
        unknown = [s for s in built if s not in POSTAL_FIPS]
        if unknown:
            ap.error(f"unknown state code(s): {unknown}")
        try:
            new = bulk(args, built)
        except requests.HTTPError as exc:
            if exc.response is None or exc.response.status_code != 429:
                raise
            # Every completed response is cached, so nothing is lost: rerun later
            # (or with USGS_API_KEY set) and the build resumes where it stopped.
            logger.error("rate limited; nothing written. Rerun to resume. (%s)", exc)
            return 2

    if not args.all and args.out.exists():
        old = load_catalog(args.out)
        keep = old[~old["site_no"].isin(new["site_no"])]
        if not args.sites:
            keep = keep[~keep["state"].isin(built)]
        new = pd.concat([keep, new], ignore_index=True)
    new = new.drop_duplicates("site_no", keep="last")
    new = new.sort_values(["state", "site_no"]).reset_index(drop=True)
    write_catalog(new, args.out)
    load_catalog(args.out)  # validate what was written

    meta = {
        "built": dt.date.today().isoformat(),
        "rows": int(len(new)),
        "states": int(new["state"].nunique()),
        "min_years": args.min_years,
        "peaks_source": "USGS Water Data OGC API, collections/peaks (parameter_code 00060)",
        "sites_source": "USGS Water Data OGC API, collections/monitoring-locations",
        "regulation_source": "flowfreq.regulation (GAGES-II screen + peak code 6), issue #32",
        "regression_region_source": (
            "StreamStats nss/regions MapServer peak-flow polygons, gage point-in-polygon; "
            "WA, OR, ID, MT only"
        ),
        "regulation_class_counts": {
            str(k): int(v) for k, v in new["regulation_class"].value_counts().items()
        },
        "regression_region_filled": int(new["regression_region"].notna().sum()),
        "tool": "tools/build_gage_catalog.py",
    }
    args.meta.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    logger.info(
        "wrote %d rows (%.0f KiB) to %s", len(new), args.out.stat().st_size / 1024, args.out
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
