"""Peak records and analysis settings for the Wyoming/Montana test sites.

Reads the two CSVs the peakfq R package ships in its own test data:

``wymt_ffa_2022A_EMPdata_7_4.csv``
    One row per annual peak -- ``site_no``, ``peak_WY``, ``peak_va`` -- which
    is the *input* record. The station fixtures in ``wyoming_montana.py`` carry
    only expected outputs, so this is the one place the peaks themselves live.

``wymt_ffa_2022A_EXPinfo_7_4.csv``
    One row per site: analysis period, skew option, regional skew and its MSE,
    and peakfq's own results.

Those results are from peakfq **7.4**, not the 8.1.0 the repository vendors, so
they are a sanity cross-check and never a parity target. Parity is against
golden files generated from ``vendor/peakfqr`` by
``tools/gen_fortran_golden.py``, exactly as for Big Sandy.

Everything here resolves through :mod:`tests.fixtures.paths`, so it degrades to
a clean skip if the reference tree is ever absent rather than failing.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from tests.fixtures.paths import TESTDATA_AVAILABLE, testdata_path

_PEAKS_CSV = "wymt_ffa_2022A_EMPdata_7_4.csv"
_INFO_CSV = "wymt_ffa_2022A_EXPinfo_7_4.csv"


@dataclass(frozen=True)
class WymtSite:
    """One site's input record and the settings peakfq analysed it under.

    Attributes
    ----------
    site_no : str
        USGS station number as it appears in the CSVs, e.g. ``"06326500.00"``.
    name : str
        Station name.
    peaks : dict of int to float
        Water year to annual peak discharge (cfs).
    begyear, endyear : int
        Analysis period.
    regional_skew, regional_skew_mse : float
        ``RegSkew`` and ``RegMSEG``. A skew option of ``"Station"`` is reported
        as ``RegSkew = -999``; those sites carry ``None`` here instead, because
        -999 is a sentinel and not a skew.
    skew_option : str
        ``"Weighted"`` or ``"Station"``.
    n_historical : int
        ``HistPeaks``. Zero for every site used as a parity case, which is why
        their inputs need no perception thresholds.
    expected_74 : dict of str to float
        peakfq 7.4's own results -- cross-check only, never a parity target.
    """

    site_no: str
    name: str
    peaks: Dict[int, float]
    begyear: int
    endyear: int
    regional_skew: Optional[float]
    regional_skew_mse: float
    skew_option: str
    n_historical: int
    expected_74: Dict[str, float] = field(default_factory=dict)

    @property
    def is_contiguous(self) -> bool:
        """True when every water year in the record has a peak.

        A parity case built from this site assumes so: with no gaps and no
        historic peaks there is nothing to censor, so the EMA inputs are just
        the peaks themselves.
        """
        years = sorted(self.peaks)
        return years == list(range(years[0], years[-1] + 1))


def _rows(name: str) -> List[dict]:
    with open(testdata_path(name), newline="") as handle:
        return list(csv.DictReader(handle))


def _site_key(raw: str) -> str:
    """Normalise a site number for lookup; the CSVs quote it, e.g. '06326500.00'."""
    return raw.strip().strip('"')


def load_site(site_no: str) -> WymtSite:
    """Load one site's peaks and settings.

    Parameters
    ----------
    site_no : str
        Station number as it appears in the CSVs, e.g. ``"06326500.00"``.

    Returns
    -------
    WymtSite

    Raises
    ------
    FileNotFoundError
        If the peakfqr reference test data is not present.
    KeyError
        If no site matches, with the available site numbers in the message.
    """
    if not TESTDATA_AVAILABLE:
        raise FileNotFoundError(
            "peakfqr reference test data not present; see tests/fixtures/paths.py"
        )

    info_rows = {_site_key(r["site_no"]): r for r in _rows(_INFO_CSV)}
    key = _site_key(site_no)
    if key not in info_rows:
        raise KeyError(f"no site {site_no!r}; available: {sorted(info_rows)}")
    info = info_rows[key]

    peaks: Dict[int, float] = {}
    for row in _rows(_PEAKS_CSV):
        if _site_key(row["site_no"]) != key:
            continue
        # A negative water year is peakfq's marker for a historic peak. No site
        # used as a parity case has any, and silently folding one in as if it
        # were systematic would misstate the record.
        year = int(row["peak_WY"])
        if year < 0:
            raise ValueError(
                f"site {site_no} carries a historic peak (water year {year}); "
                "build its perception thresholds explicitly rather than "
                "treating it as a systematic record"
            )
        peaks[year] = float(row["peak_va"])

    reg_skew = float(info["RegSkew"])
    return WymtSite(
        site_no=key,
        name=info["station_nm"],
        peaks=peaks,
        begyear=int(info["BegYear"]),
        endyear=int(info["EndYear"]),
        # -999 is the "station skew, no regional value" sentinel.
        regional_skew=None if reg_skew <= -998.0 else reg_skew,
        regional_skew_mse=float(info["RegMSEG"]),
        skew_option=info["SkewOption"],
        n_historical=int(info["HistPeaks"]),
        expected_74={
            "skew": float(info["Skew"]),
            "mean": float(info["Mean"]),
            "std_dev": float(info["StandDev"]),
            "at_site_skew": float(info["AtSiteSkew"]),
            "at_site_mseg": float(info["AtSiteMSEG"]),
            "pilf_threshold": float(info["PILF_Thresh"]),
            "n_pilf": float(info["PILFs"]),
        },
    )
