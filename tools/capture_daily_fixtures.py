"""Capture trimmed live daily-value responses for tests/test_waterdata_daily.py.

Needs network access to both the Water Data OGC API and legacy NWIS. For each
window it saves the OGC ``daily`` collection response (``daily_*.json``) and
the legacy NWIS daily-values RDB for the same site and window
(``dv_*.rdb``), so the offline tests can check that the two parsers agree on
real payloads. Rewrites only the ``daily_*`` / ``dv_*`` files under
tests/fixtures/waterdata_ogc/; see tests/fixtures/waterdata_ogc_responses.py for
what each capture is. Trimming keeps only ``self``/``next`` links; feature rows
and RDB lines are verbatim.

Run: python tools/capture_daily_fixtures.py
"""

import json
from pathlib import Path

import requests

B = "https://api.waterdata.usgs.gov/ogcapi/v1/collections"
LEGACY_DV = "https://waterservices.usgs.gov/nwis/dv/"
OUT = Path(__file__).resolve().parents[1] / "tests" / "fixtures" / "waterdata_ogc"
OUT.mkdir(parents=True, exist_ok=True)


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


def daily(site, time, limit=50000):
    return get(
        f"{B}/daily/items",
        {
            "f": "json",
            "monitoring_location_id": f"USGS-{site}",
            "parameter_code": "00060",
            "statistic_id": "00003",
            "time": time,
            "skipGeometry": "true",
            "limit": limit,
        },
    )


def legacy(site, start, end):
    r = requests.get(
        LEGACY_DV,
        params={
            "format": "rdb",
            "sites": site,
            "parameterCd": "00060",
            "statCd": "00003",
            "startDT": start,
            "endDT": end,
        },
        timeout=120,
    )
    r.raise_for_status()
    # Drop the "retrieved:" stamp so a re-capture diffs only on data.
    lines = [l for l in r.text.splitlines() if not l.startswith("# retrieved:")]
    name = f"dv_{site}_{start}_{end}"
    (OUT / f"{name}.rdb").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(name, sum(1 for l in lines if l.startswith("USGS")))


# (site, start, end, fixture name, what it shows)
WINDOWS = [
    # Approved, with ESTIMATED and EQUIP days (legacy A:e).
    ("03606500", "2024-03-20", "2024-03-31", "daily_03606500_00060_equip"),
    # Approved ice-affected estimates: numeric, ["ESTIMATED", "ICE"] (legacy A:e).
    ("12449500", "2022-12-15", "2022-12-28", "daily_12449500_00060_ice_estimated"),
    # Provisional ice days with no value: null + ["ICE"] (legacy text "Ice").
    ("06214500", "2026-01-10", "2026-01-25", "daily_06214500_00060_ice_null"),
]

for site, start, end, name in WINDOWS:
    save(name, trim_links(daily(site, f"{start}/{end}")))
    legacy(site, start, end)

# Paging: the Big Sandy window at limit=5 -> three pages, rows unsorted.
site, start, end, _ = WINDOWS[0]
page = daily(site, f"{start}/{end}", limit=5)
n = 1
while True:
    nxt = [l["href"] for l in page["links"] if l["rel"] == "next"]
    save(f"daily_03606500_00060_page{n}", trim_links(page))
    if not nxt:
        break
    page = get(nxt[0])
    n += 1

# Before the record: 200 with zero features.
save("daily_03606500_00060_empty", trim_links(daily("03606500", "1900-01-01/1900-01-05")))
