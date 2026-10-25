"""Nationwide vignette: ``estimate_at_site`` at one Washington and one Tennessee site.

Unlike ``wave1_wave2_ungaged.py``, which is offline from captured StreamStats values,
this calls the live services: StreamStats delineates the basin from a pour point and
returns its characteristics and regression regions, then flowfreq evaluates the state's
equations and prints the report. It needs the network and takes a minute or two.

* Washington (verified library): Little Klickitat River near Goldendale, treated as
  ungaged; its gage's Bulletin 17C 1 % AEP flow is printed for comparison.
* Tennessee (NSS snapshot library, so ``allow_nss_snapshot=True``): Sycamore Creek near
  Ashland City, a Cumberland tributary northwest of Nashville. Downtown Nashville is not
  on StreamStats' stream network and the Cumberland River there is outside the
  equations' range, so a smaller stream nearby stands in.

Run::

    python docs/vignettes/ungaged_nationwide.py
"""

from __future__ import annotations

import json
from pathlib import Path

from flowfreq import estimate_at_site

SITES = Path(__file__).with_name("wave1_wave2_sites.json")


def main() -> None:
    wa = json.loads(SITES.read_text(encoding="utf-8"))["WA"]
    lat, lon = wa["outlet"]
    print(f"Washington: {wa['site_name']} (gage B17C Q1% = {wa['b17c_q1pct']:,.0f} cfs)\n")
    print(estimate_at_site(lat, lon, state="WA").to_markdown())

    print("\n" + "=" * 72 + "\n")
    print("Tennessee: Sycamore Creek near Ashland City (USGS 03431800), NSS snapshot\n")
    est = estimate_at_site(36.320022, -87.051194, state="TN", allow_nss_snapshot=True)
    print(est.to_markdown())


if __name__ == "__main__":
    main()
