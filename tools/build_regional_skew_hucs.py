"""Regenerate the hydrologic-unit definitions of the regional skew studies.

Writes two packaged tables that :mod:`flowfreq.regional_skew` reads offline:

``flowfreq/data/regional_skew_hucs.csv``
    Which hydrologic units each skew study covers (``member``), which it
    explicitly excludes (``excluded``), and which its report leaves
    ambiguous (``unresolved``). The most specific (longest) matching prefix
    decides.

``flowfreq/data/regional_skew_huc8_states.csv``
    The states each HUC8 touches, from the Watershed Boundary Dataset, so a
    lookup by HUC alone can tell a Wyoming or Canadian unit of HUC2 17 from an
    Idaho one.

Sources, all live-verified when this script was written (2026-09-27):

- Pacific Northwest (PNW) study region: SIR 2016-5118 p. 1 (abstract) and
  p. 23 -- "the Columbia River Basin south of the U.S.-Canada border, the
  Oregon and Washington coastal areas, and the remainder of Washington",
  i.e. the areas in Oregon, Washington, Idaho and western Montana within the
  Columbia drainage, the Oregon/western-Washington coast and Puget Sound.
  That is HUC4 1701-1709 (Columbia, including the Snake and the Willamette),
  1710 (Oregon-Washington coastal) and 1711 (Puget Sound), restricted to
  ID, MT, OR and WA. All 290 study gages (SIR 2016-5083 Table B1) are in
  those states; 285 are in those HUC4s.
- The five that are not -- four in HUC4 1712 (Oregon closed basins) and one
  in 1801 (Klamath, 11502500) -- were used by the study, but the stated region
  does not include those basins. They are recorded ``unresolved``: the lookup
  refuses them rather than decide for the report.
- Snake River Plain exclusion: SIR 2016-5083 p. 52 says the PNW model is not
  valid in its "region 0". Region 0 is the part of Idaho outside every 2016-5083
  peak-flow region (p. 17; fig. 3, p. 9). Those region polygons are served by
  StreamStats (``nss/regions`` layer 11, ``Peak_Flow_Region_*_2016_5083``). The
  gap does not follow HU boundaries, so a 12-digit HU is recorded ``excluded``
  only when at least 95% of its area lies in it, and ``boundary`` when 2-95%
  does. A boundary HU cannot be decided from its code alone -- study gage
  13116500 (Medicine Lodge Creek) sits in region 6_8 inside an HU that is 79%
  region 0 -- so ``regional_skew_at(lat, lon)`` resolves it with a point query
  against the same layer, and a HUC-only lookup raises. Area fractions are
  estimated by grid sampling (~250 m).
- HU boundaries and states: USGS Watershed Boundary Dataset map service
  (``hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer``).

Only ``requests``, numpy and matplotlib (for point-in-polygon) are used; all
are flowfreq dependencies.

Usage::

    python tools/build_regional_skew_hucs.py            # rewrite both tables
    python tools/build_regional_skew_hucs.py --check    # fail if they differ
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import requests
from matplotlib.path import Path as MplPath

ROOT = Path(__file__).resolve().parents[1]
HUCS_CSV = ROOT / "flowfreq" / "data" / "regional_skew_hucs.csv"
STATES_CSV = ROOT / "flowfreq" / "data" / "regional_skew_huc8_states.csv"

WBD = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer"
SS_REGIONS = "https://gis.streamstats.usgs.gov/arcgis/rest/services/nss/regions/MapServer/11"

PNW = "Pacific Northwest"
SRP = "Snake River Plain"
PNW_STATES = "ID;MT;OR;WA"
PNW_SOURCE = (
    "SIR 2016-5118 p. 1 (abstract) and p. 23 (https://doi.org/10.3133/sir20165118); "
    "SIR 2016-5083 app. B (https://doi.org/10.3133/sir20165083)"
)
SRP_SOURCE = (
    "SIR 2016-5083 p. 52 (PNW skew invalid in region 0), p. 17 and p. 9 fig. 3 (region 0 defined); "
    "region polygons from StreamStats nss/regions layer 11 (Peak_Flow_Region_*_2016_5083)"
)

#: HUC4s of the PNW study region, with the part of the report's wording each
#: one satisfies.
PNW_HUC4: Dict[str, str] = {
    "1701": "Kootenai-Pend Oreille-Spokane (Columbia basin)",
    "1702": "Upper Columbia (Columbia basin)",
    "1703": "Yakima (Columbia basin)",
    "1704": "Upper Snake (Columbia basin)",
    "1705": "Middle Snake (Columbia basin)",
    "1706": "Lower Snake (Columbia basin)",
    "1707": "Middle Columbia (Columbia basin)",
    "1708": "Lower Columbia (Columbia basin)",
    "1709": "Willamette (Columbia basin)",
    "1710": "Oregon-Washington coastal",
    "1711": "Puget Sound",
}

UNRESOLVED: Dict[str, str] = {
    "1712": (
        "Oregon closed basins. Not in the study region as SIR 2016-5118 pp. 1, 23 define it "
        "(Columbia basin, OR/WA coast, remainder of WA), yet 4 of the 290 study gages are "
        "here (SIR 2016-5083 Table B1, e.g. 10396000). Unresolved: use "
        "regional_skew_for('OR', 'Pacific Northwest') deliberately, or consult the USGS."
    ),
    "1801": (
        "Klamath-Northern California coastal. Not in the study region as SIR 2016-5118 "
        "pp. 1, 23 define it, yet 1 of the 290 study gages (11502500, Williamson River, OR) "
        "is here (SIR 2016-5083 Table B1). Unresolved for the same reason as 1712; the "
        "California part is outside the study entirely."
    ),
}

#: A 12-digit HU is excluded when at least this fraction of its area is in region 0,
EXCLUDE_FRACTION = 0.95
#: and on the boundary (decided per point) when more than this fraction is.
BOUNDARY_FRACTION = 0.02
#: Grid spacing, degrees, for area-fraction sampling.
SAMPLE_STEP = 0.0025


def _get(url: str, params: Dict[str, Any], what: str) -> Dict[str, Any]:
    for attempt in range(4):
        try:
            r = requests.get(url, params=params, timeout=300)
            r.raise_for_status()
            payload = r.json()
            if "error" in payload:
                raise RuntimeError(f"{what}: {payload['error']}")
            return dict(payload)
        except (requests.RequestException, ValueError, RuntimeError):
            if attempt == 3:
                raise
            time.sleep(5 * (attempt + 1))
    raise AssertionError("unreachable")


def _query_all(layer_url: str, where: str, fields: str, geometry: bool) -> List[Dict[str, Any]]:
    """Every feature of an ArcGIS layer query, paging with resultOffset."""
    out: List[Dict[str, Any]] = []
    offset = 0
    while True:
        params = {
            "where": where,
            "outFields": fields,
            "returnGeometry": "true" if geometry else "false",
            "outSR": "4326",
            "maxAllowableOffset": "0.0005",
            "orderByFields": "OBJECTID",
            "resultOffset": offset,
            "resultRecordCount": 500,
            "f": "json",
        }
        payload = _get(f"{layer_url}/query", params, f"{layer_url} {where}")
        feats = payload.get("features") or []
        out.extend(feats)
        if not feats or not payload.get("exceededTransferLimit"):
            return out
        offset += len(feats)


def _paths(rings: Sequence[Sequence[Sequence[float]]]) -> List[MplPath]:
    return [MplPath(np.asarray(r, dtype=float)) for r in rings if len(r) >= 3]


def _inside(paths: Sequence[MplPath], pts: np.ndarray) -> np.ndarray:
    """Even-odd point-in-polygon over every ring (handles holes and multipart)."""
    count = np.zeros(len(pts), dtype=int)
    for p in paths:
        count += p.contains_points(pts)
    return count % 2 == 1


def huc8_states() -> List[Tuple[str, str, str]]:
    feats = _query_all(
        f"{WBD}/4", "huc8 LIKE '17%' OR huc8 LIKE '1801%'", "huc8,name,states", geometry=False
    )
    rows = {
        (f["attributes"]["huc8"], f["attributes"]["name"], f["attributes"]["states"]) for f in feats
    }
    return sorted(rows)


def snake_river_plain() -> List[Tuple[str, str, str, float]]:
    """(huc12, rule, name, fraction) for HUs of HUC4 1704/1705 touching 2016-5083 region 0."""
    regions = _query_all(
        SS_REGIONS, "NAME LIKE 'Peak_Flow_Region%2016_5083'", "NAME", geometry=True
    )
    names = sorted(f["attributes"]["NAME"] for f in regions)
    if len(names) != 6:
        raise RuntimeError(f"expected the six 2016-5083 peak-flow regions, got {names}")
    region_paths = [p for f in regions for p in _paths(f["geometry"]["rings"])]

    out: List[Tuple[str, str, str, float]] = []
    for huc4 in ("1704", "1705"):
        hus = _query_all(f"{WBD}/6", f"huc12 LIKE '{huc4}%'", "huc12,name,states", geometry=True)
        for f in hus:
            a = f["attributes"]
            if "ID" not in (a["states"] or "").split(","):
                continue  # region 0 is inside Idaho
            hu_paths = _paths(f["geometry"]["rings"])
            allpts = np.concatenate([p.vertices for p in hu_paths])
            (x0, y0), (x1, y1) = allpts.min(0), allpts.max(0)
            gx, gy = np.meshgrid(
                np.arange(x0, x1 + SAMPLE_STEP, SAMPLE_STEP),
                np.arange(y0, y1 + SAMPLE_STEP, SAMPLE_STEP),
            )
            pts = np.c_[gx.ravel(), gy.ravel()]
            pts = pts[_inside(hu_paths, pts)]
            if len(pts) == 0:
                continue
            frac = float((~_inside(region_paths, pts)).mean())
            if frac >= EXCLUDE_FRACTION:
                if a["states"] != "ID":
                    raise RuntimeError(f"{a['huc12']} ({a['states']}) straddles Idaho's border")
                out.append((a["huc12"], "excluded", a["name"], frac))
            elif frac > BOUNDARY_FRACTION:
                out.append((a["huc12"], "boundary", a["name"], frac))
    n_ex = sum(1 for r in out if r[1] == "excluded")
    print(f"Snake River Plain: {n_ex} HU12 excluded, {len(out) - n_ex} on the boundary")
    return sorted(out)


def build() -> Tuple[str, str]:
    hucs = io.StringIO()
    w = csv.writer(hucs, lineterminator="\n")
    w.writerow(["huc", "skew_region", "rule", "valid_states", "source", "notes"])
    for huc4, what in PNW_HUC4.items():
        w.writerow([huc4, PNW, "member", PNW_STATES, PNW_SOURCE, what])
    for huc4, why in UNRESOLVED.items():
        w.writerow([huc4, PNW, "unresolved", "", PNW_SOURCE, why])
    for huc12, rule, name, frac in snake_river_plain():
        w.writerow([huc12, SRP, rule, "", SRP_SOURCE, f"{name}; {frac:.0%} of area in region 0"])

    states = io.StringIO()
    w = csv.writer(states, lineterminator="\n")
    w.writerow(["huc8", "name", "states"])
    for row in huc8_states():
        w.writerow(row)
    return hucs.getvalue(), states.getvalue()


def main(argv: Sequence[str] = ()) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="fail if the tables would change")
    args = parser.parse_args(list(argv))
    hucs, states = build()
    if args.check:
        stale = [
            p.name
            for p, text in ((HUCS_CSV, hucs), (STATES_CSV, states))
            if not p.exists() or p.read_text(encoding="utf-8") != text
        ]
        if stale:
            print(f"out of date: {', '.join(stale)}")
            return 1
        print("up to date")
        return 0
    HUCS_CSV.write_text(hucs, encoding="utf-8")
    STATES_CSV.write_text(states, encoding="utf-8")
    print(f"wrote {HUCS_CSV.relative_to(ROOT)} and {STATES_CSV.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
