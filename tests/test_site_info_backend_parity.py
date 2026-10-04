"""Legacy NWIS vs. Water Data OGC API: does ``fetch_site_info`` agree?

The gate for moving ``USGSgage.fetch_site_info`` to ``backend="waterdata-ogc"``,
after the legacy site service returned 503s on 2026-10-03. For each site both
backends must give the same:

- ``site_name`` and ``drainage_area`` and ``state_code``, exactly;
- ``daily_por_start`` and ``iv_por_start``, exactly (the API's UTC
  ``begin`` converted to the gage's local day);
- ``latitude``/``longitude`` to 1e-6 degree (~0.1 m).

Where the services genuinely differ, the comparison is looser, on purpose:

- **HUC.** The API reports the 12-digit ``hydrologic_unit_code``, legacy the
  8-digit ``huc_cd``; the 8 digits must be the 12's prefix.
- **Coordinates.** Same NAD83 position (the API's ``geometry``, legacy's
  ``dec_lat_va``/``dec_long_va`` with ``dec_coord_datum_cd`` NAD83), but legacy
  rounds to 7-8 decimal places and the API carries ~15 significant digits.
- **End dates.** Each service refreshes its catalog on its own schedule. On
  2026-10-03 legacy listed daily means to 2026-10-02 and the API to
  2026-10-01, and legacy's IV end was "today" where the API's last reading
  converted to local time fell late the previous day. The ends must agree
  within :data:`END_DATE_SLACK_DAYS`.

Live, 2026-10-03, every site below passed:

=========  ==========  ==========  ==========  ==========
site       daily start  (legacy/API ends)       iv start
=========  ==========  ======================  ==========
03606500   1929-08-01  2026-10-02 / 2026-10-01  2002-01-01
12449500   1919-06-01  2026-10-02 / 2026-10-01  1991-04-10
12358500   1939-10-01  2026-10-02 / 2026-10-01  1995-10-01
=========  ==========  ======================  ==========

Needs *both* hosts. The legacy side skips -- rather than fails -- when NWIS is
unreachable, blocked or answering 5xx, as the other backend-parity tests do.
Run with::

    pytest tests/test_site_info_backend_parity.py -m requires_network -v
"""

from __future__ import annotations

import re
from typing import Optional

import pandas as pd
import pytest

from flowfreq.usgs import USGSgage

SITES = ["03606500", "12449500", "12358500"]

#: Largest gap allowed between the two backends' end dates.
END_DATE_SLACK_DAYS = 3

#: Coordinate agreement, decimal degrees.
COORD_TOL = 1e-6


def _fetch(site: str, backend: str) -> USGSgage:
    gage = USGSgage(site)
    gage.fetch_site_info(use_local_first=False, backend=backend)
    return gage


def _legacy(site: str) -> USGSgage:
    """The legacy backend's attributes, or a skip when the service is down.

    ``fetch_site_info`` logs a failed request rather than raising, so the skip
    reads the error it recorded.
    """
    gage = _fetch(site, "nwis-legacy")
    error = getattr(gage, "_last_api_error", None)
    if error:
        if re.search(r"\b5\d\d (Server Error|error)", error):
            pytest.skip(f"legacy NWIS unavailable from here: {error}")
        if "403" in error or "407" in error or "Connection" in error or "timed out" in error:
            pytest.skip(f"legacy NWIS unreachable from here: {error}")
        pytest.fail(f"legacy site service failed for {site}: {error}")
    return gage


def _days_apart(a: Optional[str], b: Optional[str]) -> int:
    assert a is not None and b is not None
    return abs((pd.Timestamp(a) - pd.Timestamp(b)).days)


@pytest.mark.requires_network
@pytest.mark.parametrize("site", SITES)
def test_site_info_backends_agree(site: str) -> None:
    ogc = _fetch(site, "waterdata-ogc")
    assert ogc._last_api_error is None, ogc._last_api_error
    legacy = _legacy(site)

    for name in ("site_name", "drainage_area", "state_code", "daily_por_start", "iv_por_start"):
        assert getattr(ogc, name) == getattr(legacy, name), name
        assert getattr(ogc, name) is not None, name

    assert ogc.latitude == pytest.approx(legacy.latitude, abs=COORD_TOL)
    assert ogc.longitude == pytest.approx(legacy.longitude, abs=COORD_TOL)

    assert ogc.huc is not None and legacy.huc is not None
    assert len(ogc.huc) == 12 and len(legacy.huc) == 8
    assert ogc.huc.startswith(legacy.huc)

    assert _days_apart(ogc.daily_por_end, legacy.daily_por_end) <= END_DATE_SLACK_DAYS
    assert _days_apart(ogc.iv_por_end, legacy.iv_por_end) <= END_DATE_SLACK_DAYS


@pytest.mark.requires_network
def test_methow_daily_period_is_discharge_on_both_backends() -> None:
    """12449500's catalog lists 2002 water temperature first; neither backend may use it."""
    ogc = _fetch("12449500", "waterdata-ogc")
    assert ogc.daily_por_start == "1919-06-01"
    legacy = _legacy("12449500")
    assert legacy.daily_por_start == "1919-06-01"
