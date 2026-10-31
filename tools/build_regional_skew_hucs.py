"""Regenerate the hydrologic-unit definitions of the regional skew studies.

Writes the packaged table :mod:`flowfreq.regional_skew` reads offline,
``flowfreq/data/regional_skew_hucs.csv``: which hydrologic units each skew study
covers (``member``), which it explicitly excludes (``excluded``), and which its
report leaves ambiguous (``unresolved``). For the site's state, the most specific
(longest) matching prefix decides; a HUC may carry one row per study when the rows'
``valid_states`` are disjoint (overlapping studies, #204).

The states each HUC8 touches (``regional_skew_huc8_states.csv``) are built by
``tools/build_huc8_states.py`` (#205); this script reads WBD only to select the
Montana and Arizona HUC8s, and no longer writes that table.

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
- Arizona study region (added 2026-10-11): SIR 2014-5211 p. 1 (abstract) and
  p. 54 -- a "statewide regional skew" (B-GLS constant model), so the study
  covers Arizona and stops at the state line, although its 307 gages include
  some in neighboring states (app. 1, p. 49). Every HUC8 the WBD lists as
  touching AZ is a ``member`` valid in AZ only; one shared with NM, UT, NV,
  CA, CO or Mexico therefore needs the site's state.
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

AZ = "Arizona"
AZ_SOURCE = (
    "SIR 2014-5211 p. 1 (abstract), p. 54 and Table 1-2 p. 55 "
    "(https://doi.org/10.3133/sir20145211)"
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

#: Provisional (flowfreq-tested, not USGS) extension of the PNW study.
OR_CLOSED = "Pacific Northwest (Oregon closed basins, provisional)"
OR_CLOSED_SOURCE = (
    "SIR 2016-5083 app. B (https://doi.org/10.3133/sir20165083) value; applicability tested "
    "in docs/OREGON_UNRESOLVED_SKEW.md (tools/pnw_edge_skew_test.py)"
)
OR_CLOSED_HUC4: Dict[str, str] = {
    "1712": (
        "Oregon closed basins. Not in the PNW study region as SIR 2016-5118 pp. 1, 23 define "
        "it, though 4 study gages are here (SIR 2016-5083 Table B1). A B-GLS test of 8 "
        "screened Oregon gages here is consistent with G = -0.07 (t = -0.59); provisional, "
        "Oregon only."
    ),
}

#: Provisional flowfreq study (not USGS): Montana east of the Continental Divide.
MT_EAST = "Montana east of the Divide (provisional)"
MT_EAST_SOURCE = (
    "flowfreq Montana skew study v2 (docs/MONTANA_REGIONAL_SKEW_V2.md; "
    "tools/build_montana_skew_study_v2.py); NOT a USGS study"
)

UNRESOLVED: Dict[str, str] = {
    "1801": (
        "Klamath-Northern California coastal. Not in the study region as SIR 2016-5118 "
        "pp. 1, 23 define it, yet 1 of the 290 study gages (11502500, Williamson River, OR) "
        "is here (SIR 2016-5083 Table B1). Unresolved: only 2 screened Oregon gages, too few "
        "to test (docs/OREGON_UNRESOLVED_SKEW.md); the California part is outside the study."
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
        f"{WBD}/4",
        "huc8 LIKE '17%' OR huc8 LIKE '1801%' OR states LIKE '%MT%' OR states LIKE '%AZ%'",
        "huc8,name,states",
        geometry=False,
    )
    rows = {
        (f["attributes"]["huc8"], f["attributes"]["name"], f["attributes"]["states"]) for f in feats
    }
    return sorted(rows)


def montana_east_huc8() -> List[Tuple[str, str, str]]:
    """(huc8, name, states) for HUC8s in regions 10 and 09 that touch Montana."""
    rows = [r for r in huc8_states() if r[0][:2] in ("10", "09") and "MT" in r[2].split(",")]
    if not rows:
        raise RuntimeError("WBD returned no HUC8 in regions 09/10 touching MT")
    print(f"Montana east of the Divide: {len(rows)} HUC8")
    return rows


def arizona_huc8() -> List[Tuple[str, str, str]]:
    """(huc8, name, states) for every HUC8 the WBD lists as touching Arizona."""
    rows = [r for r in huc8_states() if "AZ" in r[2].split(",")]
    if not rows:
        raise RuntimeError("WBD returned no HUC8 touching AZ")
    shared = sum(1 for r in rows if r[2] != "AZ")
    print(f"Arizona: {len(rows)} HUC8 ({shared} shared with another state or Mexico)")
    return rows


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


#: Studies whose regions the reports define by hydrologic unit, transcribed without a
#: network query (#204). ``valid_states`` is the states whose verified row in
#: ``regional_skew.csv`` adopts the study, so overlapping studies (HU 05 is Great Lakes
#: and Ohio in Ohio, Skew Region 2 in West Virginia) never give a state two rows.
GLO = "Great Lakes and Ohio River Basins"
GLO_STATES = "IN;MI;NY;OH;PA"
GLO_SOURCE = (
    "SIR 2019-5105 p. 2 (https://doi.org/10.3133/sir20195105): most of the Great Lakes "
    "Basin (HU 04) and part of the Ohio River Basin (HU 05); HUs 0511 and 0513 and Canada "
    "were not considered"
)
ENYPA = "Eastern New York and Pennsylvania"
ENYPA_STATES = "NY;PA"
ENYPA_SOURCE = (
    "SIR 2021-5015 pp. 1-2 (https://doi.org/10.3133/sir20215015): hydrologic units 0202, "
    "0204, 0205, 0206, 0207 and 0110; the skew does not apply in HU 0203"
)
SR1 = "Skew Region 1 (Atlantic drainage)"
SR2 = "Skew Region 2 (Gulf of America drainage)"
SR_SOURCE = (
    "SIR 2025-5110 p. 25 and fig. 2.1 (https://doi.org/10.3133/sir20255110): Skew Region 1 "
    "is the study area draining to the Atlantic, Skew Region 2 the study area draining to "
    "the Gulf of America (VA, WV, KY, TN and part of western MD)"
)
NE = "New England"
NE_SOURCE = (
    "SIR 2017-5037 p. 2 (https://doi.org/10.3133/sir20175037): CT, ME, NH, RI, VT, MA "
    "excluding Cape Cod and the islands, and easternmost NY"
)
ARLA = "Arkansas-Louisiana"
ARLA_SOURCE = (
    "SIR 2016-5081 app. 1 pp. 120-121 (https://doi.org/10.3133/sir20165081): Arkansas, "
    "Louisiana, southern Missouri and eastern Oklahoma; not valid in 13 southern Louisiana "
    "HUC8s and part of 12040201"
)
#: SIR 2016-5081 app. 1 pp. 120-121, verbatim.
ARLA_EXCLUDED = (
    "08070100", "08070204", "08070300", "08080101", "08080102", "08080103", "08080202",
    "08080206", "08090100", "08090202", "08090203", "08090302", "08090301",
)  # fmt: skip


AK_RSA1 = "Regional Skew Area 1"
AK_RSA2 = "Regional Skew Area 2"
AK_EXCLUSION = "Regional skew exclusion area"
AK_SOURCE = (
    "SIR 2016-5024 Table 5 p. 17 (https://doi.org/10.3133/sir20165024): HUC8s, or HUC10s "
    "where a HUC8 is split; outside both areas only station skew applies"
)
#: SIR 2016-5024 Table 5, transcribed from the page image: HUC8 -> listed HUC10s, or
#: ``()`` for "All". Area 2's '19020104218' (11 digits) is read as 1902010418, the only
#: WBD unit that fits (Hinchinbrook Island). 1902010413 is printed but is not in the
#: 2026 WBD; its row is kept as printed and matches nothing.
AK_TABLE_5: Dict[str, Dict[str, Tuple[str, ...]]] = {
    AK_RSA1: {
        **{h: () for h in ("19020101", "19020102", "19020103")},
        "19020104": tuple(f"19020104{n:02d}" for n in range(1, 5)),
        "19020301": tuple(f"19020301{n:02d}" for n in range(1, 11)),
        **{h: () for h in ("19020302", "19020401", "19020402", "19020501", "19020502",
                           "19020503", "19020505")},
        **{f"190803{n:02d}": () for n in range(1, 10)},
        **{f"190804{n:02d}": () for n in range(1, 5)},
        "19090101": (),
        "19090102": (),
    },
    AK_RSA2: {
        **{f"190101{n:02d}": () for n in range(2, 6)},
        "19010106": ("1901010601", "1901010602", "1901010603"),
        "19010107": ("1901010701",),
        "19010204": (),
        "19010206": (),
        "19010207": ("1901020701", "1901020702", "1901020703"),
        **{f"190102{n:02d}": () for n in range(8, 13)},
        **{f"190103{n:02d}": () for n in range(1, 4)},
        "19010304": ("1901030404",),
        "19010402": (),
        "19010404": tuple(f"19010404{n:02d}" for n in range(2, 6)),
        "19010405": (),
        "19010406": (),
        "19020104": tuple(f"19020104{n:02d}" for n in range(5, 19)),
        **{f"190202{n:02d}": () for n in range(1, 4)},
        "19020301": ("1902030111", "1902030112"),
        "19020602": tuple(f"19020602{n:02d}" for n in range(8, 14)),
        "19020701": (),
        "19020702": (),
        "19030101": (),
        "19030102": ("1903010201", "1903010202", "1903010203", "1903010204", "1903010206"),
    },
}  # fmt: skip


def alaska_rows() -> List[Tuple[str, str, str, str, str, str]]:
    """Alaska's two skew areas (SIR 2016-5024 Table 5); the rest of HU 19 is excluded."""
    rows: List[Tuple[str, str, str, str, str, str]] = [
        ("19", AK_EXCLUSION, "excluded", "AK", AK_SOURCE, "outside Regional Skew Areas 1 and 2"),
    ]
    for region, units in AK_TABLE_5.items():
        for huc8, huc10s in units.items():
            if not huc10s:
                rows.append((huc8, region, "member", "AK", AK_SOURCE, "whole HUC8"))
            for huc10 in huc10s:
                note = "HUC10 of a split HUC8"
                if huc10 == "1902010418":
                    note += "; printed '19020104218' in Table 5"
                elif huc10 == "1902010413":
                    note += "; not in the 2026 WBD"
                rows.append((huc10, region, "member", "AK", AK_SOURCE, note))
    return rows


def report_defined_rows() -> List[Tuple[str, str, str, str, str, str]]:
    """``(huc, skew_region, rule, valid_states, source, notes)`` rows from report text."""
    rows: List[Tuple[str, str, str, str, str, str]] = [
        ("04", GLO, "member", GLO_STATES, GLO_SOURCE, "Great Lakes Basin"),
    ]
    rows += [
        (f"05{n:02d}", GLO, "member", GLO_STATES, GLO_SOURCE, "Ohio River Basin")
        for n in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14)
    ]
    rows += [
        (h, ENYPA, "member", ENYPA_STATES, ENYPA_SOURCE, what)
        for h, what in (
            ("0110", "Connecticut Coastal"),
            ("0202", "Upper Hudson"),
            ("0204", "Delaware"),
            ("0205", "Susquehanna"),
            ("0206", "Upper Chesapeake"),
            ("0207", "Potomac"),
        )
    ]
    rows.append(
        ("0203", ENYPA, "excluded", "NJ;NY", ENYPA_SOURCE, "Lower Hudson-Long Island: "
         "unregulated, non-urban data were insufficient")
    )  # fmt: skip
    rows += [
        ("02", SR1, "member", "VA;WV", SR_SOURCE, "Atlantic drainage (Mid-Atlantic)"),
        ("03", SR1, "member", "VA;WV", SR_SOURCE, "Atlantic drainage (Chowan-Roanoke)"),
        ("05", SR2, "member", "KY;MD;TN;VA;WV", SR_SOURCE, "Gulf drainage (Ohio)"),
        ("06", SR2, "member", "KY;TN;VA", SR_SOURCE, "Gulf drainage (Tennessee)"),
        ("08", SR2, "member", "KY;TN", SR_SOURCE, "Gulf drainage (Lower Mississippi)"),
    ]
    rows += [
        ("01", NE, "member", "CT;MA;NH;RI", NE_SOURCE, "New England region"),
        ("01090002", NE, "excluded", "MA", NE_SOURCE, "Cape Cod, the islands and "
         "southeasternmost Massachusetts (Atlantic Coastal Pine Barrens)"),
    ]  # fmt: skip
    rows += [
        ("08", ARLA, "member", "AR;LA", ARLA_SOURCE, "Lower Mississippi"),
        ("11", ARLA, "member", "AR;LA", ARLA_SOURCE, "Arkansas-White-Red"),
        ("1201", ARLA, "member", "LA", ARLA_SOURCE, "Sabine"),
        ("12040201", ARLA, "unresolved", "LA", ARLA_SOURCE, "Sabine Lake: the report "
         "excludes part of this HUC8"),
    ]  # fmt: skip
    rows += [
        (h, ARLA, "excluded", "LA", ARLA_SOURCE, "tidal, regulated, diverted or backwater")
        for h in ARLA_EXCLUDED
    ]
    return rows


def build() -> str:
    hucs = io.StringIO()
    w = csv.writer(hucs, lineterminator="\n")
    w.writerow(["huc", "skew_region", "rule", "valid_states", "source", "notes"])
    for huc4, what in PNW_HUC4.items():
        w.writerow([huc4, PNW, "member", PNW_STATES, PNW_SOURCE, what])
    for huc4, why in OR_CLOSED_HUC4.items():
        w.writerow([huc4, OR_CLOSED, "member", "OR", OR_CLOSED_SOURCE, why])
    for huc8, name, touches in montana_east_huc8():
        w.writerow([huc8, MT_EAST, "member", "MT", MT_EAST_SOURCE, f"{name}; WBD states {touches}"])
    for huc4, why in UNRESOLVED.items():
        w.writerow([huc4, PNW, "unresolved", "", PNW_SOURCE, why])
    for huc12, rule, name, frac in snake_river_plain():
        w.writerow([huc12, SRP, rule, "", SRP_SOURCE, f"{name}; {frac:.0%} of area in region 0"])
    for huc8, name, touches in arizona_huc8():
        w.writerow([huc8, AZ, "member", "AZ", AZ_SOURCE, f"{name}; WBD states {touches}"])
    for row in report_defined_rows():
        w.writerow(row)
    for row in alaska_rows():
        w.writerow(row)
    return hucs.getvalue()


def main(argv: Sequence[str] = ()) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="fail if the tables would change")
    args = parser.parse_args(list(argv))
    hucs = build()
    if args.check:
        if not HUCS_CSV.exists() or HUCS_CSV.read_text(encoding="utf-8") != hucs:
            print(f"out of date: {HUCS_CSV.name}")
            return 1
        print("up to date")
        return 0
    HUCS_CSV.write_text(hucs, encoding="utf-8")
    print(f"wrote {HUCS_CSV.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
