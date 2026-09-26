"""Build ``flowfreq/data/gage_catalog.csv`` from the USGS peak and site services.

This needs network access to USGS hosts, which Claude Code web sessions do not
have (see ``TODO.md``, "Environment constraints"). Run it from a CLI session.

Usage::

    python tools/build_gage_catalog.py --state WA --sites wa_sites.txt
    python tools/build_gage_catalog.py --state WA --sites wa_sites.txt --append

``--sites`` is a text file with one USGS site number per line. Discovering
every peak site in a state automatically is left to the Water Data API
backend (issue #29), once that has been verified live. Rows with fewer than
``--min-years`` peaks are dropped. The output is never edited by hand;
rerun the tool instead.

Roadmap: ``docs/MASTER_ROADMAP.md`` §1.1, issue #33.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from flowfreq.catalog import CATALOG_PATH, COLUMNS, catalog_row, load_catalog  # noqa: E402
from flowfreq.peak_sources import get_backend  # noqa: E402

logger = logging.getLogger("build_gage_catalog")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--state", required=True)
    ap.add_argument("--sites", required=True, type=Path)
    ap.add_argument("--min-years", type=int, default=10)
    ap.add_argument("--backend", default="nwis-legacy")
    ap.add_argument("--out", type=Path, default=CATALOG_PATH)
    ap.add_argument("--append", action="store_true", help="merge into an existing catalog")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    from flowfreq.usgs import USGSgage

    backend = get_backend(args.backend)
    sites = [
        s.strip()
        for s in args.sites.read_text().splitlines()
        if s.strip() and not s.startswith("#")
    ]
    rows, errors = [], {}
    for site in sites:
        try:
            peaks = backend.fetch_peaks(site)
            gage = USGSgage(site)
            gage.fetch_site_info(use_local_first=False)
            row = catalog_row(
                site,
                args.state,
                peaks,
                site_name=gage.site_name or "",
                latitude=gage.latitude,
                longitude=gage.longitude,
                drainage_area_sqmi=gage.drainage_area,
            )
        except Exception as exc:  # one bad site never aborts the batch
            errors[site] = f"{type(exc).__name__}: {exc}"
            continue
        if row["n_peaks"] >= args.min_years:
            rows.append(row)
    new = pd.DataFrame(rows, columns=list(COLUMNS))
    if args.append and args.out.exists():
        old = load_catalog(args.out)
        new = pd.concat([old[~old["site_no"].isin(new["site_no"])], new], ignore_index=True)
    new = new.sort_values(["state", "site_no"]).reset_index(drop=True)
    new.to_csv(args.out, index=False)
    load_catalog(args.out)  # validate what was written
    logger.info("wrote %d rows to %s; %d site(s) failed", len(new), args.out, len(errors))
    for site, msg in errors.items():
        logger.warning("%s: %s", site, msg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
