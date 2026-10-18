"""Extract each gage's published regression region from WY, NV and NM station tables.

``tools/build_gage_catalog.py`` fills ``regression_region`` by point-in-polygon
wherever region polygons exist. Wyoming, Nevada and New Mexico have none: NSS
returns no geometry for WY's and NV's regions, and NM's nine regions all carry
one statewide polygon (see ``NSS_REGION_IDS`` there). Each report does print
the region of every gage it used, though, and this tool transcribes those tables
into ``data/catalog_regions/``:

- **WY**, Miller (2003), WRIR 03-4107, Table 10 (pp. 52-61): map number, station,
  latitude, longitude and region for 364 gages.
- **NV**, Thomas, Hjalmarson and Waltemeyer (1997), WSP 2433, data section
  "Basin, climatic, and flood characteristics for streamflow-gaging stations"
  (pp. 109-195): station, name and flood region for 1,320 gages in 16 regions.
  The PDF is a scanned image with an OCR text layer, so the region is read by
  word position, from the cell under the "Flood region" header.
- **NM**, Waltemeyer (2008), SIR 2008-5119, Appendix 1: station, name and map
  identifier for 293 gages, grouped under "... flood region (N)" headers.

Each table is checked against the report before it is written: WY's and NM's
per-region counts must equal ``n_sites`` in ``flowfreq/data/regression/{WY,NM}.json``
(the counts the reports' tables give). NV's hybrid regions 6 and 10 use every
listed station, so those counts must equal 80 and 104 exactly (WSP 2433 Tables 10
and 14); its GLS regions drop undefined and outlier stations, so each listed count
must be at least the table's station count.

The tables keep the report's own region numbers. The catalog builder maps them
to NSS codes and accepts a row for a gage only when its name (and, where the table
has one, its latitude) agrees with the gage's Water Data location record, so a
misread or misprinted station number is never matched to some other gage.

Needs poppler's ``pdftotext`` (``-bbox`` for NV, ``-layout`` for WY and NM; the
xpdf build has no ``-bbox``).

Usage::

    python tools/extract_station_regions.py --pdftotext /path/to/poppler/pdftotext
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "data" / "catalog_regions"
LIBRARY = REPO_ROOT / "flowfreq" / "data" / "regression"

COLUMNS = ("site_no", "report_region", "report_name", "latitude", "where")

#: state -> (output file, report PDF).
REPORTS: Dict[str, tuple] = {
    "WY": ("WY_wri034107.csv", "https://pubs.usgs.gov/wri/wri034107/pdf/wri034107.pdf"),
    "NV": ("NV_wsp2433.csv", "https://pubs.usgs.gov/wsp/2433/report.pdf"),
    "NM": ("NM_sir20085119.csv", "https://pubs.usgs.gov/sir/2008/5119/pdf/sir2008-5119.pdf"),
}

#: WSP 2433's data section, PDF pages (printed pp. 109-195 are PDF pp. 118-204).
NV_PAGES = (118, 205)
#: WSP 2433 Tables 5, 6, 7, 9, 10 and 14: stations each Nevada region's equations use.
NV_TABLE_COUNTS = {1: 165, 2: 108, 3: 35, 5: 37, 6: 80, 10: 104}
#: The hybrid regions use every listed station (p. 16), so their counts are exact.
NV_HYBRID = (6, 10)

_STATION = re.compile(r"\d{8}|\d{15}")


def library_counts(state: str) -> Dict[int, int]:
    """``n_sites`` per numeric region code in the state's packaged equations."""
    data = json.loads((LIBRARY / f"{state}.json").read_text(encoding="utf-8"))
    return {
        int(e["region_code"]): int(e["n_sites"])
        for e in data["equations"]
        if str(e["region_code"]).isdigit()
    }


def parse_wy(layout_text: str) -> List[Dict[str, str]]:
    """WRIR 03-4107 Table 10 rows: map number, station, latitude, longitude, region."""
    lines = layout_text.split("\n")
    start = next(i for i, ln in enumerate(lines) if "Table 10. Basin characteristics" in ln)
    row = re.compile(
        r"^\s+(\d{1,3})\s+(\d{8}|\d{15})\s+(\d{2}\.\d{4})\s+(\d{2,3}\.\d{4})\s+(\d)\s+[\d,.]+\s"
    )
    rows: Dict[int, Dict[str, str]] = {}
    for i, ln in enumerate(lines[start:], start):
        m = row.match(ln)
        if m and int(m.group(1)) not in rows:
            rows[int(m.group(1))] = {
                "site_no": m.group(2),
                "report_region": m.group(5),
                "report_name": "",
                "latitude": m.group(3),
                "where": f"Table 10, map number {m.group(1)}",
            }
    return [rows[k] for k in sorted(rows)]


def parse_nm(layout_text: str) -> List[Dict[str, str]]:
    """SIR 2008-5119 Appendix 1 rows: station, name and map id under region headers."""
    lines = layout_text.split("\n")
    start = next(i for i, ln in enumerate(lines) if ln.startswith("Appendix 1. Selected basin"))
    header = re.compile(r"flood region \((\d)\)", re.I)
    # Station, a name (it has letters), the characteristics, and the map id last. The
    # appendix's second table (peak discharges) has no name and does not match.
    row = re.compile(r"^\s*(\d{8}|\d{15})\s+([A-Za-z].*?)\s{2,}.*?\s(\d{1,3})\s*$")
    region: Optional[str] = None
    rows: Dict[int, Dict[str, str]] = {}
    for ln in lines[start:]:
        m = header.search(ln)
        if m:
            region = m.group(1)
            continue
        m = row.match(ln)
        if m and region and int(m.group(3)) not in rows:
            rows[int(m.group(3))] = {
                "site_no": m.group(1),
                "report_region": region,
                "report_name": m.group(2).strip(),
                "latitude": "",
                "where": f"Appendix 1, map identifier {m.group(3)}",
            }
    return [rows[k] for k in sorted(rows)]


_WORD = re.compile(
    r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([^<]*)</word>'
)


def parse_nv(bbox_html: str, first_pdf_page: int) -> List[Dict[str, str]]:
    """WSP 2433 data-section rows from ``pdftotext -bbox`` word boxes.

    A listing page has a "Flood region" column; the station's region is the one- or
    two-digit word centered under that header's "region" word on the station's line.
    """
    rows: List[Dict[str, str]] = []
    for k, page in enumerate(bbox_html.split("<page ")[1:]):
        words = [
            (float(a), float(b), float(c), float(d), e) for a, b, c, d, e in _WORD.findall(page)
        ]
        hdr = [w for w in words if w[4] == "region" and w[1] < 200]
        if not hdr:
            continue  # a peak-discharge page, or the section's introduction
        lat_hdr = [w for w in words if w[4].startswith("Lati") and w[1] < 200]
        # OCR misses the "Lati-" header on one page; the column starts ~13 pt after "region".
        lat_x = lat_hdr[0][0] if lat_hdr else hdr[0][2] + 13
        x0, x1 = hdr[0][0] - 2, hdr[0][2] + 4
        for s in words:
            # OCR puts a comma in three station numbers (10261,530).
            site = s[4].replace(",", "")
            if not (_STATION.fullmatch(site) and s[0] < 100 and s[1] > hdr[0][3]):
                continue
            line = sorted(
                (w for w in words if abs(w[1] - s[1]) < 3 and w is not s), key=lambda w: w[0]
            )
            region = [
                w for w in line if x0 <= (w[0] + w[2]) / 2 <= x1 and re.fullmatch(r"\d{1,2}", w[4])
            ]
            if len(region) != 1:
                raise ValueError(f"PDF p. {first_pdf_page + k}: no single region cell for {site}")
            name = " ".join(w[4] for w in line if w[2] < x0).replace("&apos;", "'")
            lat_text = re.sub(
                r"\s*\.\s*", ".", " ".join(w[4] for w in line if lat_x - 6 <= w[0] < lat_x + 30)
            )
            lat = re.fullmatch(r"\d{2}\.\d{3}", lat_text)
            rows.append(
                {
                    "site_no": site,
                    "report_region": region[0][4],
                    "report_name": name,
                    "latitude": lat.group(0) if lat else "",
                    "where": f"data section, PDF p. {first_pdf_page + k}",
                }
            )
    return rows


def check_counts(state: str, rows: Sequence[Mapping[str, str]]) -> Dict[int, int]:
    """Raise unless the per-region counts match what the report says (module docstring)."""
    got = collections.Counter(int(r["report_region"]) for r in rows)
    if state == "NV":
        for region, n in NV_TABLE_COUNTS.items():
            ok = got[region] == n if region in NV_HYBRID else got[region] >= n
            if not ok:
                raise ValueError(f"NV region {region}: {got[region]} listed, table uses {n}")
    else:
        want = library_counts(state)
        if dict(got) != want:
            raise ValueError(f"{state}: parsed {dict(sorted(got.items()))}, report {want}")
    return dict(sorted(got.items()))


def write_table(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator="\n")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["site_no"], r["where"])):
            w.writerow({k: r[k] for k in COLUMNS})


def _pdf(url: str, cache: Path) -> Path:
    path = cache / url.rsplit("/", 1)[-1]
    if not path.exists():
        cache.mkdir(parents=True, exist_ok=True)
        r = requests.get(url, timeout=300)
        r.raise_for_status()
        path.write_bytes(r.content)
    return path


def _pdftotext(exe: str, pdf: Path, *args: str) -> str:
    out = subprocess.run(
        [exe, *args, str(pdf), "-"], check=True, capture_output=True, encoding="utf-8"
    )
    return out.stdout


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pdftotext", default="pdftotext", help="poppler pdftotext executable")
    ap.add_argument("--states", nargs="+", default=sorted(REPORTS), choices=sorted(REPORTS))
    ap.add_argument("--cache", type=Path, default=REPO_ROOT / ".cache" / "station_regions")
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    args = ap.parse_args(argv)
    for state in args.states:
        name, url = REPORTS[state]
        pdf = _pdf(url, args.cache)
        if state == "NV":
            first, last = NV_PAGES
            html = _pdftotext(args.pdftotext, pdf, "-bbox", "-f", str(first), "-l", str(last))
            rows = parse_nv(html, first)
        else:
            text = _pdftotext(args.pdftotext, pdf, "-layout")
            rows = parse_wy(text) if state == "WY" else parse_nm(text)
        counts = check_counts(state, rows)
        write_table(args.out / name, rows)
        print(f"{state}: {len(rows)} rows, by region {counts} -> {args.out / name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
