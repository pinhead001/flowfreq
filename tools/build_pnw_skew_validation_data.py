"""Build the Pacific Northwest regional-skew validation inputs.

Writes ``tests/fixtures/skew_study/pnw_table_b1.csv``: SIR 2016-5083 appendix B
Table B1 (pp. 39-47; identical to SIR 2016-5118 Table A1), one row per
streamgage, plus the historical period of each gage's NWIS peak record through
water year 2012, which the report's pseudo concurrent record length (eq. B10)
needs and Table B1 does not print.

Table B1 columns (transcribed from the report PDF with ``pdftotext -layout``;
the input to this script is that transcription, ``--table-b1``, as a JSON list
of ``[index, site_no, PRL, lat, lon, G, MSE_G, used]`` rows):

- ``prl``: skew pseudo record length, years;
- ``lat``, ``lon``: basin centroid, decimal degrees;
- ``skew``: station skew (EMA/MGB, PeakFQ 7.1, data through WY2012);
- ``skew_mse``: its MSE;
- ``used``: ``Yes``, ``no-P`` (PRL < 35 yr) or ``no-R`` (redundant).

Added from NWIS (legacy peak service, the one the study itself used):

- ``yb``, ``ye``: first and last water year of any peak, historic included,
  capped at 2012; ``ph = ye - yb + 1`` is the historical-period length.

This is an approximation of the study's own historical periods, which the
report does not tabulate: a PeakFQ historic period can begin before the
earliest recorded peak. It only enters the cross-correlation of skews
(B8-B10), never the skews or their variances.

Usage::

    python tools/build_pnw_skew_validation_data.py --table-b1 b1rows.json
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "tests" / "fixtures" / "skew_study" / "pnw_table_b1.csv"
PEAK_URL = "https://nwis.waterdata.usgs.gov/nwis/peak"
LAST_WY = 2012


def _water_year(peak_dt: str) -> int:
    y, m = int(peak_dt[:4]), int(peak_dt[5:7] or 0)
    return y + 1 if m >= 10 else y


def periods(sites: Sequence[str]) -> Dict[str, Tuple[int, int]]:
    out: Dict[str, Tuple[int, int]] = {}
    for i in range(0, len(sites), 50):
        chunk = list(sites[i : i + 50])
        for attempt in range(4):
            try:
                r = requests.get(
                    PEAK_URL,
                    params={
                        "multiple_site_no": ",".join(chunk),
                        "agency_cd": "USGS",
                        "format": "rdb",
                    },
                    timeout=300,
                )
                r.raise_for_status()
                break
            except requests.RequestException:
                if attempt == 3:
                    raise
                time.sleep(10 * (attempt + 1))
        lines = [ln for ln in r.text.splitlines() if not ln.startswith("#")]
        df = pd.read_csv(io.StringIO("\n".join([lines[0]] + lines[2:])), sep="\t", dtype=str)
        df = df.dropna(subset=["peak_dt"])
        df["wy"] = df["peak_dt"].map(_water_year)
        df = df[df["wy"] <= LAST_WY]
        for site, g in df.groupby("site_no"):
            out[str(site)] = (int(g["wy"].min()), int(g["wy"].max()))
    return out


def main(argv: Sequence[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--table-b1", required=True, help="JSON transcription of Table B1")
    args = ap.parse_args(list(argv))
    rows: List[List[str]] = json.loads(Path(args.table_b1).read_text(encoding="utf-8"))
    df = pd.DataFrame(
        rows, columns=["index", "site_no", "prl", "lat", "lon", "skew", "skew_mse", "used"]
    )
    if len(df) != 461 or df["site_no"].duplicated().any():
        raise SystemExit("Table B1 has 461 distinct streamgages")
    counts = df["used"].value_counts().to_dict()
    if counts != {"Yes": 290, "no-P": 135, "no-R": 36}:
        raise SystemExit(f"Table B1 flags do not match the report (290/135/36): {counts}")
    per = periods(df["site_no"].tolist())
    df["yb"] = df["site_no"].map(lambda s: per.get(s, (None, None))[0])
    df["ye"] = df["site_no"].map(lambda s: per.get(s, (None, None))[1])
    df["ph"] = df["ye"] - df["yb"] + 1
    missing = df.loc[df["yb"].isna() & (df["used"] == "Yes"), "site_no"].tolist()
    if missing:
        raise SystemExit(f"no NWIS peaks for used gages: {missing}")
    df["index"] = df["index"].astype(int)
    df = df.sort_values("index")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False, lineterminator="\n")
    print(f"wrote {OUT.relative_to(ROOT)} ({len(df)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
