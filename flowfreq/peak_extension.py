"""Extension of a short annual-peak record with a nearby long one (Bulletin 17C Appendix 8).

A short-record gage ``y`` and a long-record index gage ``x`` have ``n1`` concurrent years
and the index gage ``n2`` more. Appendix 8 does two things, both in base-10 logarithms:

1. **Record augmentation.** The Matalas-Jacobs estimators (equations 8-7 and 8-8) give the
   mean and variance of the short gage as improved by the index gage's longer record, and
   equation 8-19 gives the effective record length ``n1 + ne`` of that variance.
2. **Record extension.** ``ne`` values are generated from ``ne`` of the index gage's
   non-concurrent years by the line of equation 8-20, with the intercept and slope of
   equations 8-23 and 8-24, so that the short record plus the ``ne`` values has exactly the
   augmented mean and variance. Those values are added to the record as if observed, and the
   frequency analysis (EMA, with its low-outlier test and regional skew) is run on the
   whole.

The ``ne`` values stand for the *information* the index gage adds, which is why there are
fewer of them than there are index years: a correlation short of 1 is worth less than a
year of record per year.

Reference
---------
England, J.F., Jr., Cohn, T.A., Faber, B.A., Stedinger, J.R., Thomas, W.O., Jr., Veilleux,
A.G., Kiang, J.E., and Mason, R.R., Jr., 2019, Guidelines for determining flood flow
frequency, Bulletin 17C (ver. 1.1): U.S. Geological Survey Techniques and Methods, book 4,
chap. B5, appendix 8 (after Matalas and Jacobs, 1964, and Vogel and Stedinger, 1985).

What this module cannot do
--------------------------
- **It cannot say the two records belong together.** The equations assume the logarithms
  of the concurrent peaks are jointly normal with one fixed relation. A relation that
  changed during the concurrent period (a dam, a gage move, a shift the short record's
  stationarity screen finds) is averaged over, not seen.
- **Years the short gage has and the index gage lacks take no part in the fit.** They stay
  in the extended record, so that record's mean and variance then differ slightly from the
  augmented ones; :attr:`PeakExtension.n_short_only` counts them.
- **Zero and negative peaks have no logarithm** and are left out of the fit. A year the
  short gage observed at zero flow stays observed: no peak is generated for it.
- **Censored and historical peaks are not treated.** Every value is taken as an exact
  annual peak.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

#: Fewest concurrent years Appendix 8 allows.
MIN_CONCURRENT_YEARS = 10
#: The correlation above which Appendix 8 says extension "may be suitable".
RECOMMENDED_CORRELATION = 0.80


@dataclass(frozen=True)
class PeakExtension:
    """A short peak record extended from an index gage, with every intermediate value.

    All statistics are of base-10 logarithms of flow.

    Attributes
    ----------
    n_concurrent, n_additional : int
        ``n1`` and ``n2``: concurrent years, and index-gage years outside them.
    n_short_only : int
        Years the short gage has that the index gage lacks; not in the fit.
    mean_short, mean_index_concurrent, mean_index_additional : float
        Equations 8-1 to 8-3.
    sd_short, sd_index_concurrent, sd_index_additional : float
        Square roots of equations 8-4 to 8-6.
    correlation, beta : float
        Equations 8-9 and 8-10.
    mean_augmented, variance_augmented : float
        The Matalas-Jacobs estimators, equations 8-7 and 8-8.
    effective_length_mean, effective_length_variance : float
        ``n1 + ne`` by equation 8-18 and by equation 8-19.
    n_extended : int
        ``ne``: the values added.
    mean_index_extended, intercept, slope : float
        Equation 8-20's ``x_e`` mean, ``a`` and ``b``.
    extended : pandas.DataFrame
        One row per value added: ``water_year``, ``index_flow_cfs``, ``flow_cfs``.
    notes : tuple of str
        Cautions the numbers call for (a low correlation, a small gain).
    """

    n_concurrent: int
    n_additional: int
    n_short_only: int
    mean_short: float
    mean_index_concurrent: float
    mean_index_additional: float
    sd_short: float
    sd_index_concurrent: float
    sd_index_additional: float
    correlation: float
    beta: float
    mean_augmented: float
    variance_augmented: float
    effective_length_mean: float
    effective_length_variance: float
    n_extended: int
    mean_index_extended: float
    intercept: float
    slope: float
    extended: pd.DataFrame
    notes: tuple

    def predict(self, index_flows: Sequence[float]) -> np.ndarray:
        """Equation 8-20 in flow units, for index-gage peaks in cfs."""
        x = np.log10(np.asarray(index_flows, dtype=float))
        return np.asarray(10.0 ** (self.intercept + self.slope * (x - self.mean_index_extended)))


def _positive(series: pd.Series, name: str) -> pd.Series:
    if not isinstance(series, pd.Series):
        raise TypeError(f"{name} must be a pandas Series of peak flow indexed by water year")
    out = pd.to_numeric(series, errors="coerce").dropna()
    if out.index.has_duplicates:
        raise ValueError(f"{name} has a water year more than once")
    return out[out > 0].astype(float).sort_index()


def variance_terms(n1: int, n2: int) -> tuple:
    """``A``, ``B`` and ``C`` of equations 8-14 to 8-16."""
    a = (n2 + 2) * (n1 - 6) * (n1 - 8) / (n1 - 5) + (n1 - 4) * (
        n1 * n2 * (n1 - 4) / ((n1 - 3) * (n1 - 2)) - 2 * n2 * (n1 - 4) / (n1 - 3) - 4
    )
    b = (
        6 * (n2 + 2) * (n1 - 6) / (n1 - 5)
        + 2 * (n1**2 - n1 - 14)
        + (n1 - 4)
        * (
            2 * n2 * (n1 - 5) / (n1 - 3)
            - 2 * (n1 + 3)
            - 2 * n1 * n2 * (n1 - 4) / ((n1 - 3) * (n1 - 2))
        )
    )
    c = (
        2 * (n1 + 1)
        + 3 * (n2 + 2) / (n1 - 5)
        - (n1 + 1) * (2 * n1 + n2 - 2) * (n1 - 3) / (n1 - 1)
        + (n1 - 4) * (2 * n2 / (n1 - 3) + 2 * (n1 + 1) + n1 * n2 * (n1 - 4) / ((n1 - 3) * (n1 - 2)))
    )
    return float(a), float(b), float(c)


def extend_peaks(
    short: pd.Series,
    index: pd.Series,
    *,
    n_extend: Optional[int] = None,
    years: Optional[Sequence[int]] = None,
) -> PeakExtension:
    """Extend a short annual-peak record from an index gage, as Bulletin 17C Appendix 8 does.

    Parameters
    ----------
    short, index : pandas.Series
        Annual peak flow (cfs) indexed by water year, at the short-record gage and at the
        long-record index gage. Zero and negative values are left out.
    n_extend : int, optional
        How many values to add. Defaults to the effective record gained in the variance,
        equation 8-19, rounded to the nearest whole year.
    years : sequence of int, optional
        The index-gage years to generate from, all outside the concurrent period. Defaults
        to the most recent ``n_extend`` of them, as the Bulletin's procedure does; give
        others where the most recent hold a flood that would distort the extension.

    Returns
    -------
    PeakExtension

    Raises
    ------
    TypeError
        If either record is not a Series.
    ValueError
        Fewer than :data:`MIN_CONCURRENT_YEARS` concurrent years; no index years outside
        them; a record with no variance; a correlation that is not positive; no gain in
        effective record (``ne`` below 1); ``years`` that are concurrent, missing or the
        wrong number; or a negative ``b**2`` in equation 8-24, where the Bulletin says a
        useful extension "is not possible".

    Examples
    --------
    >>> ext = extend_peaks(suwanee, etowah)             # doctest: +SKIP
    >>> ext.n_extended, round(ext.correlation, 4)       # doctest: +SKIP
    (13, 0.8519)
    >>> ext.extended                                    # doctest: +SKIP
    """
    y_all, x_all = _positive(short, "short"), _positive(index, "index")
    concurrent = y_all.index.intersection(x_all.index)
    n1 = len(concurrent)
    if n1 < MIN_CONCURRENT_YEARS:
        raise ValueError(
            f"only {n1} concurrent year(s) of positive peaks; Bulletin 17C Appendix 8 needs "
            f"at least {MIN_CONCURRENT_YEARS}"
        )
    # A year the short gage observed, even at zero flow, is not a year to generate.
    observed_years = pd.to_numeric(short, errors="coerce").dropna().index
    additional = x_all.index.difference(observed_years)
    n2 = len(additional)
    if n2 < 2:
        raise ValueError("the index gage has fewer than 2 years outside the concurrent period")
    y1 = np.log10(y_all.loc[concurrent].to_numpy())
    x1 = np.log10(x_all.loc[concurrent].to_numpy())
    x2 = np.log10(x_all.loc[additional].to_numpy())
    y1_bar, x1_bar, x2_bar = float(y1.mean()), float(x1.mean()), float(x2.mean())
    s2_y1, s2_x1, s2_x2 = (float(np.var(v, ddof=1)) for v in (y1, x1, x2))
    if s2_y1 <= 0 or s2_x1 <= 0:
        raise ValueError("a concurrent record has no variance, so there is no relation to fit")
    beta = float(np.sum((x1 - x1_bar) * (y1 - y1_bar)) / np.sum((x1 - x1_bar) ** 2))
    rho = beta * math.sqrt(s2_x1 / s2_y1)
    if rho <= 0:
        raise ValueError(
            f"the concurrent peaks are not positively correlated (correlation {rho:.3f})"
        )

    # Equations 8-7, 8-8 and 8-11.
    mu = y1_bar + n2 / (n1 + n2) * beta * (x2_bar - x1_bar)
    alpha2 = n2 * (n1 - 4) * (n1 - 1) / ((n2 - 1) * (n1 - 3) * (n1 - 2))
    variance = (
        (n1 - 1) * s2_y1
        + (n2 - 1) * beta**2 * s2_x2
        + (n2 - 1) * alpha2 * (1 - rho**2) * s2_y1
        + n1 * n2 / (n1 + n2) * beta**2 * (x2_bar - x1_bar) ** 2
    ) / (n1 + n2 - 1)

    # Equations 8-18 and 8-19.
    length_mean = n1 / (1 - n2 / (n1 + n2) * (rho**2 - (1 - rho**2) / (n1 - 3)))
    a_term, b_term, c_term = variance_terms(n1, n2)
    length_variance = (
        2.0
        / (
            2.0 / (n1 - 1)
            + n2 / ((n1 + n2 - 1) ** 2 * (n1 - 3)) * (a_term * rho**4 + b_term * rho**2 + c_term)
        )
        + 1.0
    )

    if years is not None:
        chosen = [int(year) for year in years]
        if len(set(chosen)) != len(chosen):
            raise ValueError("years repeats a water year")
        outside = [year for year in chosen if year not in set(additional)]
        if outside:
            raise ValueError(
                f"years {outside} are not index-gage years outside the concurrent period"
            )
        if n_extend is not None and n_extend != len(chosen):
            raise ValueError(f"n_extend is {n_extend} but {len(chosen)} years were given")
        ne = len(chosen)
    else:
        ne = int(round(length_variance - n1)) if n_extend is None else int(n_extend)
        if ne > n2:
            raise ValueError(f"n_extend {ne} is more than the {n2} index years available")
        chosen = [int(year) for year in additional[-ne:]] if ne > 0 else []
    if ne < 2:
        raise ValueError(
            f"the index gage adds an effective {length_variance - n1:.1f} year(s) to the "
            "variance of the short record, too little to extend it (correlation "
            f"{rho:.3f}, {n1} concurrent years)"
        )

    # Equations 8-21 to 8-24.
    xe = np.log10(x_all.loc[chosen].to_numpy())
    xe_bar, s2_xe = float(xe.mean()), float(np.var(xe, ddof=1))
    if s2_xe <= 0:
        raise ValueError("the index years chosen have no variance")
    intercept = ((n1 + ne) * mu - n1 * y1_bar) / ne
    b2 = (
        (n1 + ne - 1) * variance
        - (n1 - 1) * s2_y1
        - n1 * (y1_bar - mu) ** 2
        - ne * (intercept - mu) ** 2
    ) / ((ne - 1) * s2_xe)
    if b2 < 0:
        raise ValueError(
            "equation 8-24 gives a negative squared slope, so a useful extension with this "
            "method is not possible (Bulletin 17C Appendix 8)"
        )
    slope = math.sqrt(b2)
    flows = 10.0 ** (intercept + slope * (xe - xe_bar))

    notes = []
    if rho <= RECOMMENDED_CORRELATION:
        notes.append(
            f"The correlation of the concurrent peaks, {rho:.3f}, is not above the "
            f"{RECOMMENDED_CORRELATION:.2f} at which Bulletin 17C says record extension may be "
            "suitable."
        )
    if length_variance - n1 < 4:
        notes.append(
            f"The extension adds an effective {length_variance - n1:.1f} years, fewer than "
            "the 4 or 5 the Bulletin calls a substantial improvement."
        )
    for note in notes:
        logger.warning(note)
    return PeakExtension(
        n_concurrent=n1,
        n_additional=n2,
        n_short_only=len(y_all.index.difference(x_all.index)),
        mean_short=y1_bar,
        mean_index_concurrent=x1_bar,
        mean_index_additional=x2_bar,
        sd_short=math.sqrt(s2_y1),
        sd_index_concurrent=math.sqrt(s2_x1),
        sd_index_additional=math.sqrt(s2_x2),
        correlation=rho,
        beta=beta,
        mean_augmented=mu,
        variance_augmented=variance,
        effective_length_mean=length_mean,
        effective_length_variance=length_variance,
        n_extended=ne,
        mean_index_extended=xe_bar,
        intercept=intercept,
        slope=slope,
        extended=pd.DataFrame(
            {
                "water_year": chosen,
                "index_flow_cfs": x_all.loc[chosen].to_numpy(),
                "flow_cfs": flows,
            }
        ),
        notes=tuple(notes),
    )


def concurrent_relation(
    short: pd.Series, index: pd.Series, *, split_year: Optional[int] = None
) -> pd.DataFrame:
    """How the two gages' concurrent peaks relate in each part of the concurrent period.

    Record extension assumes one fixed relation between the gages. This sets the first
    part of the concurrent period beside the second: if the ratio of the short gage's
    peaks to the index gage's, or their correlation, differs between the parts, the
    relation moved and the extension averages over two conditions.

    Parameters
    ----------
    short, index : pandas.Series
        Annual peak flow (cfs) indexed by water year. Zero and negative values are left out.
    split_year : int, optional
        The first water year of the second part. Defaults to the middle of the concurrent
        years. Give the year of a change point a stationarity screen found.

    Returns
    -------
    pandas.DataFrame
        One row for each part and one for the whole: ``period`` (``"1930-1945"``),
        ``n``, ``median_ratio`` (short over index), ``ratio_log_sd`` (the standard
        deviation of log10 of the ratio) and ``correlation`` (of the log10 peaks; NaN with
        fewer than 3 years).

    Raises
    ------
    ValueError
        No concurrent years, or a ``split_year`` that leaves a part empty.
    """
    y_all, x_all = _positive(short, "short"), _positive(index, "index")
    years = y_all.index.intersection(x_all.index)
    if len(years) < 2:
        raise ValueError("the two records have fewer than 2 concurrent years")
    split = int(years[len(years) // 2]) if split_year is None else int(split_year)
    parts = [years[years < split], years[years >= split], years]
    if len(parts[0]) == 0 or len(parts[1]) == 0:
        raise ValueError(
            f"split_year {split} leaves one part of the concurrent period "
            f"({int(years[0])} to {int(years[-1])}) empty"
        )
    rows = []
    for part in parts:
        y, x = np.log10(y_all.loc[part].to_numpy()), np.log10(x_all.loc[part].to_numpy())
        ratio = y - x
        correlation = float(np.corrcoef(x, y)[0, 1]) if len(part) >= 3 else float("nan")
        rows.append(
            {
                "period": f"{int(part[0])}-{int(part[-1])}",
                "n": len(part),
                "median_ratio": float(10.0 ** np.median(ratio)),
                "ratio_log_sd": float(np.std(ratio, ddof=1)) if len(part) >= 2 else float("nan"),
                "correlation": correlation,
            }
        )
    return pd.DataFrame(rows)


__all__ = [
    "MIN_CONCURRENT_YEARS",
    "RECOMMENDED_CORRELATION",
    "PeakExtension",
    "concurrent_relation",
    "extend_peaks",
    "variance_terms",
]
