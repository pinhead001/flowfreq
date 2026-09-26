"""Capture and trim live Water Data OGC API responses for tests/test_waterdata.py.

Needs network access. Rewrites every file under tests/fixtures/waterdata_ogc/; see
tests/fixtures/waterdata_ogc_responses.py for what each capture is and why.
Run: python tools/capture_waterdata_fixtures.py
"""

import json
from pathlib import Path

import requests

B = "https://api.waterdata.usgs.gov/ogcapi/v1/collections"
OUT = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "waterdata_ogc"
OUT.mkdir(parents=True, exist_ok=True)

LOC_KEEP = (
    "monitoring_location_name",
    "state_name",
    "drainage_area",
    "time_zone_abbreviation",
    "uses_daylight_savings",
)


def get(url, params=None):
    r = requests.get(url, params=params, timeout=120)
    r.raise_for_status()
    return r.json()


def trim_links(doc):
    doc["links"] = [l for l in doc.get("links", []) if l["rel"] in ("next", "self")]
    return doc


def save(name, doc):
    (OUT / f"{name}.json").write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print(name, doc.get("numberReturned", "-"))


def loc(site):
    d = get(f"{B}/monitoring-locations/items/USGS-{site}", {"f": "json"})
    d["properties"] = {k: d["properties"].get(k) for k in LOC_KEEP}
    d["geometry"] = None
    d.pop("links", None)
    return d


def meta(site, pc):
    d = get(
        f"{B}/time-series-metadata/items",
        {
            "f": "json",
            "monitoring_location_id": f"USGS-{site}",
            "parameter_code": pc,
            "skipGeometry": "true",
            "limit": 1000,
        },
    )
    for f in d["features"]:
        f["properties"].pop("thresholds", None)
    return trim_links(d)


def iv(site, pc, time, limit=50000, ts=None):
    p = {
        "f": "json",
        "monitoring_location_id": f"USGS-{site}",
        "parameter_code": pc,
        "time": time,
        "skipGeometry": "true",
        "limit": limit,
    }
    if ts:
        p["time_series_id"] = ts
    return get(f"{B}/continuous/items", p)


save("loc_03606500", loc("03606500"))
save("meta_03606500_00060", meta("03606500", "00060"))
save("meta_03606500_00065", meta("03606500", "00065"))
save("loc_03612600", loc("03612600"))
save("meta_03612600_00065", meta("03612600", "00065"))
save("loc_06214500", loc("06214500"))
save("meta_06214500_00060", meta("06214500", "00060"))
save("loc_06191500", loc("06191500"))
save("meta_06191500_00060", meta("06191500", "00060"))

BS = "6fcce55acfa94ec49bf59f175c3b6764"
DST = "2024-11-03T04:00:00Z/2024-11-03T10:00:00Z"
save("iv_03606500_00060_dst", trim_links(iv("03606500", "00060", DST, ts=BS)))

# Paging: the same window at limit=5 -> three pages.
page = iv("03606500", "00060", DST, limit=5, ts=BS)
n = 1
while True:
    nxt = [l["href"] for l in page["links"] if l["rel"] == "next"]
    save(f"iv_03606500_00060_page{n}", trim_links(page))
    if not nxt:
        break
    page = get(nxt[0])
    n += 1

save(
    "iv_03606500_00060_empty",
    trim_links(iv("03606500", "00060", "1990-01-01T00:00:00Z/1990-01-02T00:00:00Z", ts=BS)),
)

OH = "2024-06-01T00:00:00Z/2024-06-01T00:30:00Z"
save("iv_03612600_00065_unfiltered", trim_links(iv("03612600", "00065", OH)))
save(
    "iv_03612600_00065_tailwater",
    trim_links(iv("03612600", "00065", OH, ts="7ca465077c26415c93fd0549c697a270")),
)

save(
    "iv_06214500_00060_ice",
    trim_links(iv("06214500", "00060", "2026-01-20T00:00:00Z/2026-01-20T03:00:00Z")),
)
save(
    "iv_06191500_00060_estimated",
    trim_links(iv("06191500", "00060", "2024-01-15T00:00:00Z/2024-01-15T02:00:00Z")),
)
