"""Live check of :func:`flowfreq.estimate_at_site` against StreamStats and NSS.

Everything else in ``test_ungaged.py`` runs offline against canned responses; this is
the one place the whole chain (snap, delineate, basin characteristics, NSS region
location, equation library, evaluation) meets the real services. Run with::

    PYTHONSAFEPATH=1 PYTHONUTF8=1 pytest tests/test_ungaged_live.py -m requires_network -v

Values observed 2026-10-25:

* Little Klickitat River near Goldendale, WA (14112000; outlet 45.844124, -120.796392):
  snapped 0.02 m, 82.2 sq mi. Located regions Peak_Region_1 (GC1750, 1 %) and
  Peak_Region_4 (GC1753, 99 %). StreamStats could not compute CANOPY_PCT for the
  watershed that day ("Cannot find a 'wshforest' raster"), so region 1 was skipped and
  the estimate rests on region 4 alone: 1 % AEP 6,439 cfs (90 % PI 2,720 - 15,240),
  against the gage's Bulletin 17C 7,183 cfs. Library status ``verified``, no advisory.
* Tennessee: the brief's downtown Nashville point (36.1627, -86.7816) is not on
  StreamStats' stream network (``UnsnappablePointError``), and the Cumberland River at
  Nashville (36.159, -86.772) snaps but is far outside the equations' range
  (``OutOfRangeError``: CONTDA 0.05 below its 0.173 minimum and CSL10_85 208.5
  above its 132 maximum). The Tennessee set mixes urban and rural equation ranges
  (GC346 DA <= 30 sq mi, GC348 DA >= 30 sq mi), and none covers a 40-200 sq mi urban
  Nashville stream, so the point used is Sycamore Creek near Ashland City (USGS
  03431800; 36.320022, -87.051194), about 25 miles northwest of Nashville, in the same
  Cumberland basin: snapped 9 m, 97.7 sq mi. Located regions GC344, GC2016 and GC2020
  at 100 % each (GC1419 and GC2069 had no library region). Library status
  ``nss_snapshot`` with the "NSS SNAPSHOT (not report-verified)" advisory. 1 % AEP
  19,624 cfs, with no prediction interval (the snapshot carries no covariance).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from flowfreq import estimate_at_site

SITES = Path(__file__).parents[1] / "docs" / "vignettes" / "wave1_wave2_sites.json"


@pytest.mark.requires_network
def test_washington_little_klickitat_matches_gage():
    site = json.loads(SITES.read_text(encoding="utf-8"))["WA"]
    lat, lon = site["outlet"]
    est = estimate_at_site(lat, lon, state="WA")

    located = {r["code"]: r["percent_weight"] / 100.0 for r in est.provenance["regions_located"]}
    by_nss = {est.region_codes[c]: w for c, w in located.items() if c in est.region_codes}
    assert set(by_nss) == {"1", "4"}
    assert by_nss["4"] > 0.95

    q100 = float(est.quantiles.loc[est.quantiles["aep"] == 0.01, "flow_cfs"].iloc[0])
    gage = site["b17c_q1pct"]
    assert gage / 2 < q100 < gage * 2
    assert est.library_status == "verified"


@pytest.mark.requires_network
def test_tennessee_snapshot_estimate_carries_advisory():
    # Sycamore Creek near Ashland City, TN (03431800); see the module docstring for why
    # not downtown Nashville.
    est = estimate_at_site(36.320022, -87.051194, state="TN", allow_nss_snapshot=True)

    q100 = float(est.quantiles.loc[est.quantiles["aep"] == 0.01, "flow_cfs"].iloc[0])
    assert q100 > 0
    assert est.library_status == "nss_snapshot"
    assert est.advisory
