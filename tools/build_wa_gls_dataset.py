"""Build the Washington SIR 2016-5118 regression dataset for the GLS refit.

Northwest studies plan, task 3. Writes
``data/regression_study/wa_2016_5118_sites.csv``: one row per Table 5 gage.

Sources (Mastin and others, 2016, SIR 2016-5118 ver. 1.2,
https://doi.org/10.3133/sir20165118), the report's own spreadsheets:

- ``sir20165118_table5.xlsx``: basin characteristics, regression region and
  the "Used in Regional Regression?" flag.
- ``sir20165118_table2.xlsx``: EMA inputs. The station row's number of peaks,
  historical period length and record span.
- ``sir20165118_table8.xlsx``: flood quantiles at AEP 0.5 to 0.002. Line 1 is
  the station EMA estimate, line 2 the regression estimate, line 3 the weighted
  (WIE) estimate, and line 4 the 95-percent interval of line 3.

Line 1 is a log-Pearson III fit with the weighted skew, so its log10 moments
(``mu``, ``sd``, ``skew_w``) are recovered from its eight quantiles by least
squares. The ``fit_max_rel_err`` column records how well they reproduce them.

Usage::

    python tools/build_wa_gls_dataset.py --cache <dir>
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import requests
from scipy.optimize import minimize_scalar
from scipy.stats import pearson3

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _xlsx import read_sheet  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "regression_study" / "wa_2016_5118_sites.csv"
URL = "https://pubs.usgs.gov/sir/2016/5118/sir20165118_table{n}.xlsx"
AEPS = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)
STATION = re.compile(r"\d{8,15}")
INTERVAL = re.compile(r"\(\s*([\d.,]+)\s*-\s*([\d.,]+)\s*\)")


def _download(n: int, cache: Path) -> Path:
    path = cache / f"sir20165118_table{n}.xlsx"
    if not path.exists():
        r = requests.get(URL.format(n=n), timeout=300)
        r.raise_for_status()
        path.write_bytes(r.content)
    return path


def _num(s: str) -> float:
    s = s.strip().replace(",", "")
    return float(s) if s else float("nan")


def table5(path: Path) -> List[Dict[str, str]]:
    rows = read_sheet(path)
    header = [h.strip() for h in rows[2]]
    out = []
    for r in rows[3:]:
        if not STATION.fullmatch(r[1].strip()):
            continue
        out.append(dict(zip(header, (v.strip() for v in r))))
    return out


def table2(path: Path) -> Dict[str, Tuple[float, float, float, float]]:
    """station -> (number of peaks, historical period length, first year, last year)."""
    out = {}
    for r in read_sheet(path):
        if len(r) > 9 and STATION.fullmatch(r[1].strip()):
            out[r[1].strip()] = (_num(r[4]), _num(r[5]), _num(r[8]), _num(r[9]))
    return out


def table8(path: Path) -> Dict[str, Dict[str, List[float]]]:
    """station -> {station, regression, weighted, ci_lo, ci_hi} quantiles at AEPS.

    Section A gages have 4 lines (station, regression, weighted, interval of
    the weighted). Sections B (basins extending into Canada) and C (> 5 percent
    impervious) have 2: station, and its interval; their regression and
    weighted entries are NaN.
    """
    rows = read_sheet(path)
    out: Dict[str, Dict[str, List[float]]] = {}
    section = "A"
    nan8 = [float("nan")] * len(AEPS)
    for i, r in enumerate(rows):
        if r[0].strip().rstrip(".") in ("B", "C"):
            section = r[0].strip().rstrip(".")
        if not STATION.fullmatch(r[1].strip()):
            continue
        n_lines = 4 if section == "A" else 2
        lines = [rows[i + k][7:15] for k in range(n_lines)]
        lo, hi = [], []
        for cell in lines[-1]:
            m = INTERVAL.search(cell)
            lo.append(_num(m.group(1)) if m else float("nan"))
            hi.append(_num(m.group(2)) if m else float("nan"))
        out[r[1].strip()] = {
            "station": [_num(v) for v in lines[0]],
            "regression": [_num(v) for v in lines[1]] if n_lines == 4 else nan8,
            "weighted": [_num(v) for v in lines[2]] if n_lines == 4 else nan8,
            "ci_lo": lo,
            "ci_hi": hi,
        }
    return out


def lp3_moments(q: Sequence[float]) -> Tuple[float, float, float, float]:
    """(mu, sd, skew, max relative error) of the LP3 that best fits the quantiles."""
    y = np.log10(np.asarray(q, dtype=float))
    p = 1.0 - np.asarray(AEPS)

    def fit(g: float) -> Tuple[float, float, float]:
        k = pearson3.ppf(p, g)
        A = np.column_stack([np.ones_like(k), k])
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        return float(np.sum((A @ coef - y) ** 2)), float(coef[0]), float(coef[1])

    grid = np.linspace(-3.0, 3.0, 121)
    g0 = float(grid[int(np.argmin([fit(g)[0] for g in grid]))])
    g = float(
        minimize_scalar(lambda x: fit(x)[0], bounds=(g0 - 0.06, g0 + 0.06), method="bounded").x
    )
    _, mu, sd = fit(g)
    est = 10 ** (mu + sd * pearson3.ppf(p, g))
    rel = float(np.max(np.abs(est / 10**y - 1.0)))
    return mu, sd, g, rel


def build(cache: Path) -> List[Dict[str, object]]:
    t5 = table5(_download(5, cache))
    t2 = table2(_download(2, cache))
    t8 = table8(_download(8, cache))
    out: List[Dict[str, object]] = []
    for r in t5:
        site = r["STANUM"]
        row: Dict[str, object] = {
            "site_no": site,
            "state": r["STATE"],
            "region": r["Regression Region"],
            "used": r["Used in Regional Regression?"].lower() == "yes",
            "lat_cent": r["LAT_CENT"],
            "lon_cent": r["LONG_CENT"],
            "DRNAREA": r["DRNAREA"],
            "PRECIP": r["PRECIP"],
            "CANOPY_PCT": r["CANOPY_PCT"],
        }
        n_peaks, hist, y0, y1 = t2.get(site, (float("nan"),) * 4)
        row.update(n_peaks=n_peaks, hist_len=hist, first_year=y0, last_year=y1)
        q = t8.get(site)
        if q is not None:
            for key, prefix in (
                ("station", "q"),
                ("regression", "qreg"),
                ("weighted", "qw"),
                ("ci_lo", "ci_lo"),
                ("ci_hi", "ci_hi"),
            ):
                for a, v in zip(AEPS, q[key]):
                    row[f"{prefix}_{a}"] = v
            if all(np.isfinite(q["station"])) and min(q["station"]) > 0:
                mu, sd, g, rel = lp3_moments(q["station"])
                row.update(mu=round(mu, 6), sd=round(sd, 6), skew_w=round(g, 5))
                row["fit_max_rel_err"] = round(rel, 6)
        out.append(row)
    return out


def main(argv: Sequence[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cache", type=Path, required=True, help="download cache directory")
    args = ap.parse_args(list(argv))
    args.cache.mkdir(parents=True, exist_ok=True)
    rows = build(args.cache)
    fields: List[str] = []
    for r in rows:
        fields += [k for k in r if k not in fields]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n_used = sum(1 for r in rows if r["used"])
    print(f"wrote {OUT.relative_to(ROOT)}: {len(rows)} gages, {n_used} used in regression")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
