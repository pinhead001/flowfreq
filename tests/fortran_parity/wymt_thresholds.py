"""Wyoming/Montana stations whose perception thresholds are not all before the record.

Each station's ``PCPT_Thresh`` lines from ``vendor/peakfqr/inst/testdata/
wymt_ffa_2022A.psf``, in file order (``siteQT`` applies them in that order,
the later one winning). Every station here has a threshold period that
overlaps or lies inside its systematic record, or several separate periods --
exactly what ``ExpectedMomentsAlgorithm`` used to drop or collapse.

Stations with zero-flow years are left out: ``siteQT`` gives a zero an exact
row at ``Qmin`` that the native engine does not build, a separate gap.
"""

from __future__ import annotations

from typing import Dict, Tuple

THRESHOLDS: Dict[str, Dict[Tuple[int, int], float]] = {
    # A gap inside the record, 1952-1958, known to be below 104,000 cfs.
    "06185500.10": {(1942, 2022): 0.0, (1952, 1958): 104000.0},
    # A historical period running into the record, plus an in-record gap.
    "06324500.00": {(1878, 2022): 0.0, (1878, 1928): 100000.0, (1973, 1974): 33000.0},
    # The historic 1923 flood lies outside the only nonzero period.
    "06324500.01": {(1878, 2022): 0.0, (1878, 1918): 100000.0},
    # The whole second half of the record, 1993-2022, is a threshold period.
    "06324710.00": {(1976, 2022): 0.0, (1993, 2022): 30000.0},
    # Two separate in-record periods, listed out of order.
    "06325500.00": {(1947, 1978): 0.0, (1973, 1978): 3160.0, (1954, 1955): 3160.0},
    # Three in-record periods.
    "06327700.00": {
        (1955, 2004): 0.0,
        (1961, 1961): 14600.0,
        (1964, 1964): 14600.0,
        (1968, 2004): 14600.0,
    },
}

#: The subset whose weighted fit does not match ``emafitpr`` yet, for a reason
#: unrelated to perception thresholds: ``emafitb`` (``emafit.f`` lines 706-710)
#: computes the at-site skew MSE with the B17B formula alone -- no ADJE
#: censoring adjustment -- whenever MGBT is run and finds low outliers, and the
#: native engine always uses ADJE. Both stations have 17 MGBT low outliers.
B17B_MSE_SWITCH = {"06324500.00", "06324500.01"}


def bulletin17c_inputs(site_no: str) -> dict:
    """``ExpectedMomentsAlgorithm``/``build_emafit_arrays`` keyword arguments.

    Parameters
    ----------
    site_no : str
        A key of :data:`THRESHOLDS`.

    Returns
    -------
    dict
        ``peak_flows``, ``water_years``, ``historical_peaks``,
        ``perception_thresholds``, ``regional_skew``, ``regional_skew_mse``.
    """
    from tests.fixtures.wymt_peaks import load_site

    site = load_site(site_no, allow_historic=True)
    years = sorted(site.peaks)
    weighted = site.regional_skew is not None
    return {
        "peak_flows": [site.peaks[y] for y in years],
        "water_years": years,
        "historical_peaks": sorted(site.historical.items()) or None,
        "perception_thresholds": dict(THRESHOLDS[site_no]),
        "regional_skew": site.regional_skew if weighted else None,
        "regional_skew_mse": site.regional_skew_mse if weighted else None,
    }
