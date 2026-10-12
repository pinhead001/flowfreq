"""Check flowfreq's station skews against USGS PeakFQ for the SIR 2025-5019 examples.

Northwest studies plan, task 6. The Montana skew study (v1 and v2) uses
flowfreq's own EMA station skews; this compares them with PeakFQ 7.5.1's, on
identical inputs, for the example analyses USGS published with SIR 2025-5019
(Siefken and others, 2025; data release https://doi.org/10.5066/P1WHRK8H):

- ``ffa_2022examples.psf``: the PeakFQ specifications (perception thresholds,
  historic periods, PILF settings) for 63 analyses in MT/ND/SD/WY.
- ``ffa_2022examples_WATSTORE.txt``: the peaks those specifications read.
- ``ffa_2022examples.xlsx`` Table 1-5: PeakFQ's mean, standard deviation and
  station skew for each analysis.

Each analysis is converted with :func:`flowfreq.psf_convert.convert_psf` and
fitted natively (:meth:`~flowfreq.psf_convert.StationInputs.run`); the station
skew, log mean and SD are compared with Table 1-5's. Writes
``data/skew_study/sir2025_5019_station_skews.csv``.

Usage::

    python tools/check_montana_station_skews.py --cache <dir>
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Sequence

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _xlsx import read_sheet  # noqa: E402

from flowfreq.psf_convert import convert_psf  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "skew_study" / "sir2025_5019_station_skews.csv"
SB = "https://www.sciencebase.gov/catalog/item/{item}?format=json"
FILES = {
    "ffa_2022examples.psf": "6750da38d34e5c4500ced636",
    "ffa_2022examples_WATSTORE.txt": "6750c827d34e5c4500ced5e1",
    "ffa_2022examples.xlsx": "6751b9f8d34e5c4500cf1dd9",
}


def _download(cache: Path) -> None:
    for name, item in FILES.items():
        path = cache / name
        if path.exists():
            continue
        meta = requests.get(SB.format(item=item), timeout=120).json()
        url = next(f["url"] for f in meta["files"] if f["name"] == name)
        r = requests.get(url, timeout=600)
        r.raise_for_status()
        path.write_bytes(r.content)


def _num(s: str) -> float:
    s = s.strip()
    return float(s) if s not in ("", "--") else float("nan")


ANALYSIS_ID = re.compile(r"(\d+)(?:\.(\d+))?")


def normalize_id(aid: str) -> str:
    """Canonical ``SSSSSSSS.NN`` analysis ID.

    Table 1-5 stores many IDs as numbers, so they read back as floats: a dropped
    leading zero (``6411900`` for 06411900.00), a missing ``.00``, or binary
    noise (``5059000.0999999996`` for 05059000.10, ``12362500.1`` for .10).
    """
    m = ANALYSIS_ID.fullmatch(aid.strip())
    if m is None:
        return ""
    site, frac = m.group(1), m.group(2) or "0"
    suffix = round(float("0." + frac) * 100)
    return f"{site.zfill(8) if len(site) < 8 else site}.{suffix:02d}"


def table_1_5(path: Path) -> Dict[str, Dict[str, object]]:
    """analysis ID -> PeakFQ's results (Table 1-5)."""
    rows = read_sheet(path, 4)
    out: Dict[str, Dict[str, object]] = {}
    for r in rows[3:]:
        aid = normalize_id(r[0])
        if not aid:
            continue
        out[aid] = {
            "usgs_type": r[4].strip(),
            "usgs_n": _num(r[5]),
            "usgs_mean": _num(r[8]),
            "usgs_sd": _num(r[9]),
            "usgs_skew_type": r[10].strip(),
            "usgs_station_skew": _num(r[11]),
            "usgs_regional_skew": _num(r[12]),
        }
    return out


def basin(site: str) -> str:
    """Major basin by station-number prefix: the study's MT basins vs the rest."""
    if site.startswith("12"):
        return "Columbia (17)"
    if site.startswith("06"):
        return "Missouri (10)"
    if site.startswith("13"):
        return "Snake (17)"
    if site.startswith("0501") or site.startswith("0502"):
        return "Hudson Bay (09, St. Mary/Belly)"
    if len(site) > 8:
        return "other (15-digit site)"
    return "other (Red River, ND)"


def check(cache: Path) -> pd.DataFrame:
    usgs = table_1_5(cache / "ffa_2022examples.xlsx")
    stations = convert_psf(cache / "ffa_2022examples.psf")
    rows: List[Dict[str, object]] = []
    for aid, inputs in stations.items():
        site = aid.split(".")[0]
        row: Dict[str, object] = {"analysis_id": aid, "site_no": site, "basin": basin(site)}
        key = normalize_id(aid)
        if key not in usgs:
            # The .psf may name an analysis without its suffix; match on the site
            # when Table 1-5 has exactly one analysis for it.
            same = [k for k in usgs if k.split(".")[0] == site]
            key = same[0] if len(same) == 1 else key
        row.update(usgs.get(key, {}))
        try:
            res = inputs.run()
            row.update(
                flowfreq_mean=res.mean_log,
                flowfreq_sd=res.std_log,
                flowfreq_station_skew=res.skew_station,
                flowfreq_n=res.n_peaks,
            )
        except Exception as exc:  # noqa: BLE001 -- recorded per analysis
            row["error"] = f"{type(exc).__name__}: {exc}"[:200]
        rows.append(row)
    df = pd.DataFrame(rows)
    df["skew_diff"] = df["flowfreq_station_skew"] - df["usgs_station_skew"]
    return df


def main(argv: Sequence[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cache", type=Path, required=True)
    args = ap.parse_args(list(argv))
    args.cache.mkdir(parents=True, exist_ok=True)
    _download(args.cache)
    df = check(args.cache)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False, lineterminator="\n", float_format="%.6g")
    ok = df.dropna(subset=["skew_diff"])
    summary = {
        "analyses": len(df),
        "compared": len(ok),
        "median_abs_diff": round(float(ok["skew_diff"].abs().median()), 4),
        "over_0.2": ok.loc[ok["skew_diff"].abs() > 0.2, "analysis_id"].tolist(),
        "errors": (
            df.loc[df.get("error", pd.Series(dtype=str)).notna(), "analysis_id"].tolist()
            if "error" in df
            else []
        ),
    }
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
