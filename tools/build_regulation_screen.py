"""Regenerate ``flowfreq/data/regulation_screen.csv.gz`` from GAGES-II.

The table holds, per GAGES-II gage, the regulation/urbanization screen that
:mod:`flowfreq.regulation` reads offline (issue #32, roadmap section 1.1).

Source (live-verified 2026-09-30):

- Falcone, J., 2011, GAGES-II: Geospatial Attributes of Gages for Evaluating
  Streamflow: U.S. Geological Survey data release,
  https://doi.org/10.5066/P96CPHOT (ScienceBase item 631405bbd34e36012efa304a,
  file ``basinchar_and_report_sept_2011.zip``, 55 MB). 9,322 gages.
  Fields used: ``CLASS`` (Ref/Non-ref, Bas_Classif), ``STOR_NID_2009``
  (NID 2009 dam storage, megaliters per km2, HydroMod_Dams), ``NDAMS_2009``,
  ``RUNAVE7100`` (1971-2000 mean annual runoff, mm/yr, Hydro; conterminous
  US only) and ``IMPNLCD06`` (NLCD 2006 percent impervious, Pop_Infrastr;
  conterminous US only).

Classification rules, each from a published source (first match wins):

1. ``reference`` -- GAGES-II ``CLASS == "Ref"``: Falcone's own least-disturbed
   class. Dudley and others (2018) likewise let their HCDN (reference) class
   stand whatever the storage.
2. ``regulated`` -- dam storage normalized by mean annual runoff greater than
   127.8 days. Dudley, R.W., Archfield, S.A., Hodgkins, G.A., Renard, B., and
   Ryberg, K.R., 2018, Peak-streamflow trends and change-points and basin
   characteristics for 2,683 U.S. Geological Survey streamgages in the
   conterminous U.S. (ver. 3.0, April 2019): U.S. Geological Survey data
   release, https://doi.org/10.5066/P9AEGXY0. Their "Regulated" class is
   storage "in the highest 25-percentile among all study basins"; 127.8 days
   is that 75th percentile of ``dam.storage.2009.norm.mean.Q`` in their
   ``gage_characteristics.csv`` (their smallest Regulated value is 128.4).
   They normalized ``STOR_NID_2009`` by observed 1966-2015 mean flow; this
   tool uses GAGES-II ``RUNAVE7100`` instead (1 ML/km2 = 1 mm, so days =
   STOR_NID_2009 / RUNAVE7100 * 365) so that every conterminous GAGES-II gage
   gets a value. On the 2,262 gages where both exist, the two agree to a
   median ratio of 1.02 (interquartile range 0.93-1.13); ``--check-dudley``
   reports the class agreement.
3. ``urban`` -- NLCD 2006 impervious cover greater than 5 percent.
   Mastin, M.C., Konrad, C.P., Veilleux, A.G., and Tecca, A.E., 2016,
   Magnitude, frequency, and trends of floods at gaged and ungaged sites in
   Washington, based on data through water year 2014: U.S. Geological Survey
   Scientific Investigations Report 2016-5118, p. 23 ("urbanized basins
   (defined as basins with more than 5 percent land cover classified as
   impervious)"). The HCDN-2009 criterion in GAGES-II is the same 5 percent.
4. ``unknown`` -- non-reference with neither threshold met, or (Alaska,
   Hawaii, Puerto Rico) without the runoff and impervious fields.

Peak-code evidence (code 6) is added at run time by
:func:`flowfreq.regulation.classify_site`, not stored here.

``runave7100_mm`` (GAGES-II ``RUNAVE7100``) is stored so that
:func:`flowfreq.regulation.classify_site`'s opt-in current-NID refinement can
normalize today's storage by the same runoff the 2009 value used.

The current National Inventory of Dams is not summed here: that needs a
watershed polygon per gage, which only a StreamStats delineation provides, so
it is an opt-in, per-site refinement at run time
(``classify_site(..., use_current_nid=True)``), not a table column.

Usage::

    python tools/build_regulation_screen.py                 # rewrite the table
    python tools/build_regulation_screen.py --check         # fail if it differs
    python tools/build_regulation_screen.py --check-dudley  # agreement report
"""

from __future__ import annotations

import argparse
import io
import sys
import zipfile
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parent.parent
OUT = REPO_ROOT / "flowfreq" / "data" / "regulation_screen.csv.gz"
#: Deterministic gzip (no timestamp) so a rebuild of unchanged data is byte-identical.
GZIP = {"method": "gzip", "mtime": 0, "compresslevel": 9}

GAGES2_ITEM = "631405bbd34e36012efa304a"
GAGES2_ZIP = "basinchar_and_report_sept_2011.zip"
DUDLEY_ITEM = "5b183960e4b092d965219d62"
SB = "https://www.sciencebase.gov/catalog/item/"

SOURCE = "GAGES-II"
SOURCE_DATE = "2011-09-30"
#: Dudley and others (2018), 75th percentile of normalized dam storage, days.
REGULATED_STORAGE_DAYS = 127.8
#: Mastin and others (2016), SIR 2016-5118 p. 23, percent impervious.
URBAN_IMPERVIOUS_PCT = 5.0

COLUMNS = (
    "site_no",
    "site_class",
    "gagesii_class",
    "ndams_2009",
    "stor_nid_2009_ml_km2",
    "norm_storage_days",
    "imperv_pct_2006",
    "runave7100_mm",
    "basis",
    "source",
    "source_date",
)


def _sciencebase_file(item: str, name: str, cache: Path) -> Path:
    """Download one named file of a ScienceBase item into ``cache`` (once)."""
    dest = cache / name
    if dest.exists():
        return dest
    meta = requests.get(SB + item, params={"format": "json"}, timeout=60)
    meta.raise_for_status()
    url = next(f["url"] for f in meta.json()["files"] if f["name"] == name)
    r = requests.get(url, timeout=600)
    r.raise_for_status()
    cache.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(r.content)
    return dest


def _gages2_tables(cache: Path) -> Dict[str, pd.DataFrame]:
    outer = zipfile.ZipFile(_sciencebase_file(GAGES2_ITEM, GAGES2_ZIP, cache))
    inner = zipfile.ZipFile(io.BytesIO(outer.read("spreadsheets-in-csv-format.zip")))

    def read(name: str) -> pd.DataFrame:
        return pd.read_csv(io.BytesIO(inner.read(name)), dtype={"STAID": str}, encoding="latin-1")

    return {
        n: read(f"{n}.txt")
        for n in (
            "conterm_bas_classif",
            "conterm_hydromod_dams",
            "conterm_hydro",
            "conterm_pop_infrastr",
            "AKHIPR_bas_classif",
            "AKHIPR_hydromod_dams",
        )
    }


def classify(
    gagesii_class: str, norm_days: Optional[float], imperv: Optional[float]
) -> tuple[str, str]:
    """Apply the rules in the module docstring; return (class, basis)."""
    if str(gagesii_class).strip().lower() == "ref":
        return "reference", "GAGES-II CLASS=Ref"
    if norm_days is not None and np.isfinite(norm_days) and norm_days > REGULATED_STORAGE_DAYS:
        return "regulated", f"dam storage {norm_days:.0f} d > {REGULATED_STORAGE_DAYS} d"
    if imperv is not None and np.isfinite(imperv) and imperv > URBAN_IMPERVIOUS_PCT:
        return "urban", f"impervious {imperv:.1f}% > {URBAN_IMPERVIOUS_PCT:g}%"
    if norm_days is None or not np.isfinite(norm_days):
        return "unknown", "Non-ref; no runoff/impervious data"
    return "unknown", "Non-ref; below storage and impervious thresholds"


def build(cache: Path) -> pd.DataFrame:
    t = _gages2_tables(cache)
    con = (
        t["conterm_bas_classif"][["STAID", "CLASS"]]
        .merge(t["conterm_hydromod_dams"][["STAID", "NDAMS_2009", "STOR_NID_2009"]], on="STAID")
        .merge(t["conterm_hydro"][["STAID", "RUNAVE7100"]], on="STAID", how="left")
        .merge(t["conterm_pop_infrastr"][["STAID", "IMPNLCD06"]], on="STAID", how="left")
    )
    runoff = pd.to_numeric(con["RUNAVE7100"], errors="coerce")
    con["RUNAVE7100"] = runoff
    con["norm"] = np.where(runoff > 0, con["STOR_NID_2009"] / runoff * 365.0, np.nan)
    ak = t["AKHIPR_bas_classif"][["STAID", "CLASS"]].merge(
        t["AKHIPR_hydromod_dams"][["STAID", "NDAMS_2009", "STOR_NID_2009"]], on="STAID"
    )
    ak["norm"] = np.nan
    ak["IMPNLCD06"] = np.nan
    ak["RUNAVE7100"] = np.nan
    both = pd.concat([con, ak], ignore_index=True)
    rows = []
    for r in both.itertuples(index=False):
        # Classify on the rounded values that are stored, so the table is self-consistent.
        norm = None if pd.isna(r.norm) else round(float(r.norm), 1)
        imp = None if pd.isna(r.IMPNLCD06) else round(float(r.IMPNLCD06), 2)
        cls, basis = classify(r.CLASS, norm, imp)
        rows.append(
            {
                "site_no": str(r.STAID).strip(),
                "site_class": cls,
                "gagesii_class": "Ref" if str(r.CLASS).strip().lower() == "ref" else "Non-ref",
                "ndams_2009": int(r.NDAMS_2009),
                "stor_nid_2009_ml_km2": round(float(r.STOR_NID_2009), 2),
                "norm_storage_days": None if norm is None else round(norm, 1),
                "imperv_pct_2006": None if imp is None else round(imp, 2),
                "runave7100_mm": None if pd.isna(r.RUNAVE7100) else round(float(r.RUNAVE7100), 1),
                "basis": basis,
                "source": SOURCE,
                "source_date": SOURCE_DATE,
            }
        )
    out = pd.DataFrame(rows, columns=list(COLUMNS)).sort_values("site_no")
    if out["site_no"].duplicated().any():
        raise SystemExit("duplicate STAID in GAGES-II")
    return out.reset_index(drop=True)


def dudley_agreement(table: pd.DataFrame, cache: Path) -> pd.DataFrame:
    """Cross-tabulate this table's class against Dudley and others (2018)."""
    path = _sciencebase_file(DUDLEY_ITEM, "gage_characteristics.csv", cache)
    d = pd.read_csv(path)
    d["site_no"] = d["site_id"].str[2:]
    d["dudley"] = d["gage.class"].fillna("(none)")
    j = d.merge(table, on="site_no")
    return pd.crosstab(j["dudley"], j["site_class"])


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--cache", type=Path, default=REPO_ROOT / ".cache" / "gagesii")
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--check", action="store_true", help="fail if the table would change")
    ap.add_argument("--check-dudley", action="store_true", help="print class agreement")
    args = ap.parse_args(argv)

    table = build(args.cache)
    print(table["site_class"].value_counts().to_string())
    if args.check_dudley:
        print(dudley_agreement(table, args.cache).to_string())
    if args.check:
        old = pd.read_csv(args.out, dtype={"site_no": str})
        same = old.astype(str).equals(table.astype(str).replace({"None": "nan"}))
        print("unchanged" if same else "DIFFERS")
        return 0 if same else 1
    # ~890 KiB as plain CSV, ~120 KiB gzipped; pandas reads either directly.
    table.to_csv(args.out, index=False, compression=GZIP)
    print(f"wrote {len(table)} rows, {args.out.stat().st_size / 1024:.0f} KiB to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
