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
  ``nss_snapshot`` with the "NSS SNAPSHOT (not report-verified)" advisory. Averaged
  1/3 each (before overlap detection), the 1 % AEP was 19,624 cfs with no interval.

  The three are not a partition but alternative equation sets for the same area:
  GC344 ``MultiVariable_Area_2`` (WRIR 03-4176, 2003; CONTDA, CSL10_85), superseded by
  SIR 2024-5130 (2025), which NSS serves as two sets, GC2016
  ``TN_Peakflow_Region_2_SIR_2024-5130`` (DRNAREA, CSL10_85, LC11IMPALT) and GC2020
  ``TN_DAO_Peakflow_Region_2_SIR_2024_5130`` (drainage area only). Both 2024 sets are
  rural; they differ as multivariable versus drainage-area-only, not rural versus
  urban, so the lower-numbered one, GC2016 (also the report's full multivariable set),
  is chosen with ``region="GC2016"``: 1 % AEP 19,572 cfs (90 % PI 10,540 - 36,344).
  Without ``region`` the estimate is refused with
  :class:`~flowfreq.ungaged.OverlappingRegionsError` listing the three.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from flowfreq import estimate_at_site
from flowfreq.ungaged import OverlappingRegionsError

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


TN_LAT, TN_LON = 36.320022, -87.051194


@pytest.mark.requires_network
def test_tennessee_overlapping_regions_are_refused_without_a_choice():
    # Sycamore Creek near Ashland City, TN (03431800); see the module docstring for why
    # not downtown Nashville, and for the three overlapping equation sets.
    with pytest.raises(OverlappingRegionsError) as info:
        estimate_at_site(TN_LAT, TN_LON, state="TN", allow_nss_snapshot=True)
    years = {code: year for code, _, year in info.value.candidates}
    assert years == {"GC344": 2003, "GC2016": 2025, "GC2020": 2025}
    assert "--region" in str(info.value)


@pytest.mark.requires_network
def test_tennessee_snapshot_estimate_with_chosen_region():
    # GC2016, SIR 2024-5130's multivariable set: see the module docstring.
    est = estimate_at_site(TN_LAT, TN_LON, state="TN", allow_nss_snapshot=True, region="GC2016")

    assert [r.code for r in est.regions] == ["GC2016"]
    assert {a["code"] for a in est.provenance["alternative_regions"]} == {"GC344", "GC2020"}
    q100 = float(est.quantiles.loc[est.quantiles["aep"] == 0.01, "flow_cfs"].iloc[0])
    assert q100 > 0
    assert est.library_status == "nss_snapshot"
    assert est.advisory
