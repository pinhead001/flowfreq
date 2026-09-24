"""Pluggable sources of annual peak-flow data.

USGS is retiring the legacy NWIS web services (``nwis.waterdata.usgs.gov``,
``waterservices.usgs.gov``) in favour of the Water Data OGC APIs
(``api.waterdata.usgs.gov``). This module puts one interface in front of both,
so callers pick a backend by name and every backend returns the same frame.

Only the legacy backend works today. The Water Data API backend raises
:class:`NotImplementedError` on purpose: its endpoints have not been verified
against the live service, and this repository's rule (``TODO.md``,
``docs/STREAMSTATS_NSS_ADDENDUM.md``) is that no client is written against an
endpoint until it has been exercised live. Guessed endpoints have been wrong
three times here already.

Roadmap: ``docs/MASTER_ROADMAP.md`` §1.1, issue #29.
"""

from __future__ import annotations

import logging
from typing import Dict, Protocol, Type, runtime_checkable

import pandas as pd

logger = logging.getLogger(__name__)

#: Columns every backend returns, in this order.
PEAK_COLUMNS = ("water_year", "peak_date", "peak_flow_cfs", "qualification_code")


@runtime_checkable
class PeakDataBackend(Protocol):
    """A source of annual peak-flow records for one site."""

    name: str

    def fetch_peaks(self, site_no: str) -> pd.DataFrame:
        """Return the site's annual peaks with columns :data:`PEAK_COLUMNS`."""
        ...  # pragma: no cover


def validate_peak_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Check a backend's output against the shared contract.

    Parameters
    ----------
    df : pandas.DataFrame
        Backend output.

    Returns
    -------
    pandas.DataFrame
        ``df`` restricted to :data:`PEAK_COLUMNS`, in order.

    Raises
    ------
    ValueError
        On missing columns, duplicate water years, or negative flows.
    """
    missing = [c for c in PEAK_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Peak frame missing columns: {missing}")
    if df["water_year"].duplicated().any():
        dups = sorted(df.loc[df["water_year"].duplicated(), "water_year"].unique())
        raise ValueError(f"Duplicate water years in peak frame: {dups}")
    if (df["peak_flow_cfs"] < 0).any():
        raise ValueError("Negative peak discharge in peak frame")
    return df.loc[:, list(PEAK_COLUMNS)]


class LegacyNwisBackend:
    """The existing NWIS RDB peak service, via :class:`flowfreq.usgs.USGSgage`."""

    name = "nwis-legacy"

    def fetch_peaks(self, site_no: str) -> pd.DataFrame:
        """Download peaks from the legacy NWIS peak service.

        Parameters
        ----------
        site_no : str
            USGS site number.

        Returns
        -------
        pandas.DataFrame
        """
        from flowfreq.usgs import USGSgage  # deferred: usgs imports plotting deps

        return validate_peak_frame(USGSgage(site_no).download_peak_flow())


class WaterDataApiBackend:
    """USGS Water Data OGC API. Not implemented until live-verified (issue #29)."""

    name = "waterdata-ogc"

    def fetch_peaks(self, site_no: str) -> pd.DataFrame:
        """Not implemented; see the class docstring.

        Raises
        ------
        NotImplementedError
            Always, until the endpoint is verified live.
        """
        raise NotImplementedError(
            "The Water Data OGC API peak endpoint has not been verified live. "
            "Use backend='nwis-legacy', or complete issue #29 first."
        )


_BACKENDS: Dict[str, Type[PeakDataBackend]] = {
    LegacyNwisBackend.name: LegacyNwisBackend,
    WaterDataApiBackend.name: WaterDataApiBackend,
}

#: Backend used when none is named. Switch only after #29's parity test passes.
DEFAULT_BACKEND = LegacyNwisBackend.name


def get_backend(name: str = DEFAULT_BACKEND) -> PeakDataBackend:
    """Return a backend instance by name.

    Parameters
    ----------
    name : str
        One of ``"nwis-legacy"`` or ``"waterdata-ogc"``.

    Returns
    -------
    PeakDataBackend

    Raises
    ------
    KeyError
        For an unknown backend name.
    """
    try:
        return _BACKENDS[name]()
    except KeyError:
        raise KeyError(f"Unknown peak backend {name!r}; choose from {sorted(_BACKENDS)}") from None
