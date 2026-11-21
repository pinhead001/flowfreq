"""NCEI daily air temperature (GHCN-Daily summaries).

Reads one station's daily minimum and maximum air temperature from the NCEI Data
Service. The request goes through :func:`flowfreq.waterdata.request` for its retries,
without the USGS API key, which belongs to another host.
"""

from __future__ import annotations

import io
import logging

import pandas as pd

from flowfreq import waterdata

logger = logging.getLogger(__name__)

NCEI_URL = "https://www.ncei.noaa.gov/access/services/data/v1"


def fetch_ghcn_daily(station_id: str, start: str, end: str, timeout: int = 60) -> pd.DataFrame:
    """Daily TMIN and TMAX for one NCEI station, in degrees Celsius.

    Parameters
    ----------
    station_id : str
        GHCN station identifier, such as ``USC00123456``.
    start, end : str
        First and last day, ``YYYY-MM-DD``, inclusive.
    timeout : int, optional
        Per-attempt timeout, seconds.

    Returns
    -------
    pd.DataFrame
        Indexed by a daily DatetimeIndex named ``date`` that covers every day
        from ``start`` to ``end``, with float columns ``tmin_c`` and ``tmax_c``.
        A day the service does not report, or an element the station does not
        measure, is NaN. A response with no rows gives an all-NaN frame and a
        warning.

    Raises
    ------
    ValueError
        If ``start`` is after ``end``.
    requests.RequestException
        If the request fails after its retries.
    """
    first, last = pd.Timestamp(start), pd.Timestamp(end)
    if first > last:
        raise ValueError(f"start {start} is after end {end}")
    full = pd.date_range(first, last, freq="D", name="date")
    params = {
        "dataset": "daily-summaries",
        "stations": station_id,
        "startDate": start,
        "endDate": end,
        "dataTypes": "TMIN,TMAX",
        "units": "metric",
        "format": "csv",
    }
    response = waterdata.request(NCEI_URL, params, timeout, send_api_key=False)
    text = response.text or ""
    raw = pd.read_csv(io.StringIO(text)) if text.strip() else pd.DataFrame()
    out = pd.DataFrame(index=full, columns=["tmin_c", "tmax_c"], dtype="float64")
    if raw.empty or "DATE" not in raw.columns:
        logger.warning("NCEI returned no daily rows for %s, %s to %s", station_id, start, end)
        return out
    day = pd.to_datetime(raw["DATE"]).dt.normalize()
    for src, dst in (("TMIN", "tmin_c"), ("TMAX", "tmax_c")):
        if src in raw.columns:
            col = pd.Series(pd.to_numeric(raw[src], errors="coerce").to_numpy(), index=day)
            col = col[~col.index.duplicated(keep="last")]
            out[dst] = col.reindex(full).astype("float64")
    return out
