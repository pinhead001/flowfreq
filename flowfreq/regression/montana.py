"""Montana: weighting basin-characteristics and channel-width regression estimates.

Chase and others (2021, USGS SIR 2020-5142) publish channel-width equations for
each of Montana's eight hydrologic regions, one per method:

====== ============================================ ===================
Method Explanatory variable                         Region code in MT.json
====== ============================================ ===================
``BC`` basin characteristics (SIR 2015-5019-F)      ``W``, ``NW``, ...
``AC`` active-channel width, ``WACTCH`` (ft)        ``W-AC``, ...
``BF`` bankfull width, ``WBANKFULL`` (ft)           ``W-BF``, ...
``RS`` width from aerial photographs, ``CHANWD_RS`` ``W-RS``, ...
====== ============================================ ===================

Estimates from different methods at one site are combined in log space with the
minimum-variance weights of Parrett and Johnson (2004), which account for each
estimate's standard error of prediction (SEP) and the cross-correlation between the
methods' residuals (the report's Table 6, stored under ``cross_correlations`` in
``MT.json``). For two methods (report eqs. 11 and 12)::

    S12 = r12 SEP1 SEP2
    a1  = (SEP2^2 - S12) / (SEP1^2 + SEP2^2 - 2 S12),   a2 = 1 - a1
    Z   = a1 x1 + a2 x2
    SEPz = sqrt((SEP1^2 SEP2^2 - S12^2) / (SEP1^2 + SEP2^2 - 2 S12))

and for three (eqs. 5 to 8 and 10) the analogous closed form. ``x_i`` are log10 of
the estimates, and each ``SEP_i`` is the site's standard error of prediction in log
units, ``sqrt(model error variance + x' (X^T Lambda^-1 X)^-1 x)``. For a channel-width
equation the model error variance includes the measurement error variance (MEV;
report pp. 41-44), which ``MT.json`` already folds into ``model_error_variance``.

:func:`flowfreq.regression.evaluate_weighted` cannot express this: its weights are
fixed fractions, while these depend on the site's SEPs and the methods'
correlation.

The report's version 1.1 adds a caution: weighting can give unreasonable results when
two methods are highly correlated (r above about 0.5) and their SEPs differ (a ratio
above about 1.4 or below about 0.7). If the weighted estimate falls outside the
range of the individual estimates, use a different combination or a single equation.
:attr:`WeightedEstimate.outside_range` flags that case.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np

from flowfreq.regression.equations import RegressionEquation, RegressionEstimate, evaluate
from flowfreq.regression.library import DATA_DIR, StateLibrary, load_state

#: Method tag -> suffix appended to the hydrologic-region code in ``MT.json``.
METHOD_SUFFIX: Dict[str, str] = {"BC": "", "AC": "-AC", "BF": "-BF", "RS": "-RS"}

#: Montana's eight hydrologic regions (SIR 2015-5019-F and SIR 2020-5142).
REGIONS: Tuple[str, ...] = ("W", "NW", "NWF", "NEP", "ECP", "SEP", "UYCM", "SW")


def region_code(region: str, method: str) -> str:
    """Return the ``MT.json`` region code for a hydrologic region and method.

    Raises
    ------
    ValueError
        If the region or method is unknown.
    """
    if region not in REGIONS:
        raise ValueError(f"unknown Montana hydrologic region {region!r}; choose from {REGIONS}")
    if method not in METHOD_SUFFIX:
        raise ValueError(f"unknown method {method!r}; choose from {sorted(METHOD_SUFFIX)}")
    return region + METHOD_SUFFIX[method]


def _load_cross_correlations() -> Dict[str, Dict[str, Dict[str, float]]]:
    data = json.loads((DATA_DIR / "MT.json").read_text(encoding="utf-8"))
    regions: Dict[str, Dict[str, Dict[str, float]]] = data["cross_correlations"]["regions"]
    return regions


def cross_correlation(region: str, aep: float, method1: str, method2: str) -> float:
    """Return Table 6's cross-correlation between two methods' residuals.

    Parameters
    ----------
    region : str
        Hydrologic region code (``"W"``, ``"NW"``, ...).
    aep : float
        One of the report's AEPs (0.667, 0.5, 0.429, 0.2, 0.1, 0.04, 0.02, 0.01,
        0.005, 0.002).
    method1, method2 : {"BC", "AC", "BF", "RS"}
        Two different methods, in either order.

    Returns
    -------
    float

    Raises
    ------
    KeyError
        If the region, AEP or pair is not tabulated.
    """
    if method1 == method2:
        raise KeyError(f"cross-correlation needs two different methods, got {method1!r} twice")
    table = _load_cross_correlations()
    if region not in table:
        raise KeyError(f"no cross-correlations for region {region!r}")
    by_aep = table[region]
    key = next((k for k in by_aep if abs(float(k) - aep) < 1e-9), None)
    if key is None:
        raise KeyError(f"no cross-correlations for region {region!r} at AEP {aep}")
    row = by_aep[key]
    for pair in (f"{method1},{method2}", f"{method2},{method1}"):
        if pair in row:
            return float(row[pair])
    raise KeyError(f"no cross-correlation for methods {method1!r} and {method2!r}")


def prediction_sep(eq: RegressionEquation, characteristics: Mapping[str, float]) -> float:
    """Site standard error of prediction, log10 units: ``sqrt(MEV_total + x'Ux)``.

    Raises
    ------
    ValueError
        If the equation has no covariance matrix or model error variance.
    KeyError
        If a characteristic is missing.
    """
    if eq.covariance is None or eq.model_error_variance is None:
        raise ValueError(f"{eq.region_code} AEP {eq.aep}: no covariance to compute a site SEP")
    x = np.array([1.0] + [v.apply(float(characteristics[v.code])) for v in eq.variables])
    return math.sqrt(eq.model_error_variance + float(x @ np.asarray(eq.covariance) @ x))


def combine_two(
    x1: float, sep1: float, x2: float, sep2: float, r12: float
) -> Tuple[float, float, Tuple[float, float]]:
    """Report eqs. 11 and 12: weighted log estimate of two correlated estimates.

    Parameters
    ----------
    x1, x2 : float
        log10 of the two estimates.
    sep1, sep2 : float
        Their standard errors of prediction, log10 units.
    r12 : float
        Cross-correlation of the two methods' residuals.

    Returns
    -------
    (Z, SEPz, (a1, a2))

    Raises
    ------
    ValueError
        If an SEP is not positive or the weights are undefined.
    """
    if sep1 <= 0 or sep2 <= 0:
        raise ValueError("SEPs must be positive")
    s12 = r12 * sep1 * sep2
    denom = sep1**2 + sep2**2 - 2 * s12
    if denom <= 0:
        raise ValueError("weights undefined: the two estimates are perfectly correlated")
    a1 = (sep2**2 - s12) / denom
    a2 = 1.0 - a1
    z = a1 * x1 + a2 * x2
    sepz = math.sqrt((sep1**2 * sep2**2 - s12**2) / denom)
    return z, sepz, (a1, a2)


def combine_three(
    x: Sequence[float], sep: Sequence[float], r12: float, r13: float, r23: float
) -> Tuple[float, float, Tuple[float, float, float]]:
    """Report eqs. 5 to 8 and 10: weighted log estimate of three correlated estimates.

    Parameters
    ----------
    x : sequence of 3 floats
        log10 of the estimates from methods 1, 2 and 3.
    sep : sequence of 3 floats
        Their standard errors of prediction, log10 units.
    r12, r13, r23 : float
        Cross-correlations of the methods' residuals.

    Returns
    -------
    (Z, SEPz, (a1, a2, a3))

    Raises
    ------
    ValueError
        If an SEP is not positive or the weights are undefined.
    """
    if len(x) != 3 or len(sep) != 3:
        raise ValueError("combine_three needs exactly three estimates")
    s1, s2, s3 = sep
    if min(sep) <= 0:
        raise ValueError("SEPs must be positive")
    s12, s13, s23 = r12 * s1 * s2, r13 * s1 * s3, r23 * s2 * s3
    a_ = s1**2 + s3**2 - 2 * s13
    b_ = s3**2 + s12 - s13 - s23
    c_ = s2**2 + s3**2 - 2 * s23
    det = a_ * c_ - b_**2
    if det <= 0:
        raise ValueError("weights undefined: the estimates' covariance is singular")
    a1 = (c_ * (s3**2 - s13) - b_ * (s3**2 - s23)) / det
    a2 = (a_ * (s3**2 - s23) - b_ * (s3**2 - s13)) / det
    a3 = 1.0 - a1 - a2
    z = a1 * x[0] + a2 * x[1] + a3 * x[2]
    var = (
        (a1 * s1) ** 2
        + (a2 * s2) ** 2
        + (a3 * s3) ** 2
        + 2 * a1 * a2 * s12
        + 2 * a1 * a3 * s13
        + 2 * a2 * a3 * s23
    )
    return z, math.sqrt(var), (a1, a2, a3)


@dataclass(frozen=True)
class WeightedEstimate:
    """A Montana estimate weighted across methods (SIR 2020-5142 eqs. 5-12).

    Attributes
    ----------
    region, aep :
        Hydrologic region and AEP.
    flow_cfs : float
        ``10**log_flow``.
    log_flow : float
        The weighted log10 estimate ``Z``.
    sep_log : float
        ``SEPz``, log10 units.
    weights : dict of str to float
        Weight on each method's log estimate; they sum to 1.
    estimates : dict of str to RegressionEstimate
        Each method's estimate.
    seps : dict of str to float
        Each method's site SEP, log10 units.
    """

    region: str
    aep: float
    flow_cfs: float
    log_flow: float
    sep_log: float
    weights: Dict[str, float]
    estimates: Dict[str, RegressionEstimate]
    seps: Dict[str, float]

    @property
    def outside_range(self) -> bool:
        """Whether the weighted flow lies outside the individual estimates' range.

        The report (ver. 1.1) advises another combination, or one equation, when it does.
        """
        flows = [e.flow_cfs for e in self.estimates.values()]
        return not (min(flows) * (1 - 1e-12) <= self.flow_cfs <= max(flows) * (1 + 1e-12))


def estimate_weighted(
    region: str,
    aep: float,
    characteristics: Mapping[str, float],
    methods: Sequence[str],
    *,
    lib: Optional[StateLibrary] = None,
    allow_extrapolation: bool = False,
) -> WeightedEstimate:
    """Weight two or three Montana methods at one site (SIR 2020-5142 eqs. 5-12).

    Parameters
    ----------
    region : str
        Hydrologic region (``"W"``, ``"NW"``, ``"NWF"``, ``"NEP"``, ``"ECP"``,
        ``"SEP"``, ``"UYCM"``, ``"SW"``).
    aep : float
        AEP, one of the report's ten.
    characteristics : mapping of str to float
        Everything the chosen methods need, keyed by StreamStats code: the
        basin characteristics for ``"BC"`` and ``WACTCH``, ``WBANKFULL`` or
        ``CHANWD_RS`` (feet) for the channel-width methods.
    methods : sequence of 2 or 3 distinct method tags
        From ``"BC"``, ``"AC"``, ``"BF"`` and ``"RS"``, in the order the weights
        are reported.
    lib : StateLibrary, optional
        The Montana library; ``load_state("MT")`` when omitted.
    allow_extrapolation : bool, default False
        Passed to :func:`flowfreq.regression.evaluate`.

    Returns
    -------
    WeightedEstimate

    Raises
    ------
    ValueError
        If fewer than two or more than three methods are given, or a method repeats.
    KeyError, OutOfRangeError, EquationsUnavailable
        As :func:`flowfreq.regression.evaluate` and the library lookup raise them.
    """
    methods = list(methods)
    if len(methods) not in (2, 3) or len(set(methods)) != len(methods):
        raise ValueError(f"give two or three distinct methods, got {methods}")
    lib = lib if lib is not None else load_state("MT")
    if lib.state != "MT":
        raise ValueError(f"estimate_weighted needs the Montana library, got {lib.state}")
    estimates: Dict[str, RegressionEstimate] = {}
    seps: Dict[str, float] = {}
    for m in methods:
        eq = lib.equation(region_code(region, m), aep)
        estimates[m] = evaluate(eq, characteristics, allow_extrapolation=allow_extrapolation)
        seps[m] = prediction_sep(eq, characteristics)
    r = {(m1, m2): cross_correlation(region, aep, m1, m2) for m1, m2 in combinations(methods, 2)}
    xs = [estimates[m].log_flow for m in methods]
    ss = [seps[m] for m in methods]
    w: Tuple[float, ...]
    if len(methods) == 2:
        z, sepz, w = combine_two(xs[0], ss[0], xs[1], ss[1], r[(methods[0], methods[1])])
    else:
        m1, m2, m3 = methods
        z, sepz, w = combine_three(xs, ss, r[(m1, m2)], r[(m1, m3)], r[(m2, m3)])
    return WeightedEstimate(
        region=region,
        aep=aep,
        flow_cfs=10.0**z,
        log_flow=z,
        sep_log=sepz,
        weights=dict(zip(methods, w)),
        estimates=estimates,
        seps=seps,
    )
