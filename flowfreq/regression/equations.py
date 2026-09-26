"""Regional regression equations: schema, evaluation, prediction intervals.

An equation has the log-linear form USGS state reports use::

    log10(Q_p) = b0 + sum_i b_i * T_i(X_i)

where each ``T_i`` is the transform the report applies to basin
characteristic ``X_i`` (``log10``, ``log10(X+1)``, or none), optionally after a
linear rescaling ``scale * X_i + offset`` (see :class:`Variable`). Power-form
equations, ``Q = 10**b0 * A**b1 * ...``, are the same thing written out, so a
published term such as ``(FOREST/100 + 1)**b`` is ``log10`` with ``scale=0.01``
and ``offset=1``.

The design follows :class:`flowfreq.transpose.RegressionExponents`. An
equation cannot be built without a citation. Inputs outside the calibrated
range raise by default, and ``allow_extrapolation=True`` computes the value
anyway but records every violation. That second rule matters because NSS will
return a plausible number for a wildly out-of-range input with no warning at
all (``docs/STREAMSTATS_NSS_ADDENDUM.md``).

Roadmap: ``docs/MASTER_ROADMAP.md`` §3.1, issue #35.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy import stats

logger = logging.getLogger(__name__)

TRANSFORMS: Dict[str, Callable[[float], float]] = {
    "identity": lambda x: x,
    "log10": math.log10,
    "log10_plus1": lambda x: math.log10(x + 1.0),
}

# Exclusive lower bound on the argument each transform accepts. ``log10_plus1``
# is ``log10`` with its argument shifted by one, so its bound is -1.
_DOMAIN_LOWER: Dict[str, Optional[float]] = {
    "identity": None,
    "log10": 0.0,
    "log10_plus1": -1.0,
}


class OutOfRangeError(ValueError):
    """A basin characteristic lies outside the equation's calibrated range."""


@dataclass(frozen=True)
class Citation:
    """Where an equation came from. ``publication`` and ``table`` are required."""

    publication: str
    table: str
    year: Optional[int] = None
    url: str = ""

    def __post_init__(self) -> None:
        if not self.publication.strip() or not self.table.strip():
            raise ValueError("Citation requires both a publication and a table number")


@dataclass(frozen=True)
class Variable:
    """An explanatory variable: a StreamStats characteristic code with its limits.

    The value entering the regression is ``T(scale * x + offset)``, where ``x``
    is the characteristic in its published units and ``T`` is ``transform``.
    The rescaling lets an equation be stored exactly as the report prints it:

    =========================  ==========  =====  ======
    Published term             transform   scale  offset
    =========================  ==========  =====  ======
    ``A**b``                   log10       1      0
    ``(FOREST + 1)**b``        log10       1      1
    ``(FOREST/100 + 1)**b``    log10       0.01   1
    ``(ELEV/1000)**b``         log10       0.001  0
    ``(GUTTER + 0.1)**b``      log10       1      0.1
    ``b * (X - 20)``           identity    1      -20
    =========================  ==========  =====  ======

    ``transform="log10_plus1"`` is kept as a backward-compatible alias for
    ``log10`` with ``offset`` increased by 1: it evaluates
    ``log10(scale * x + offset + 1)``.

    Parameters
    ----------
    code : str
        StreamStats characteristic code.
    transform : {"identity", "log10", "log10_plus1"}, default "log10"
    units : str
        Published units of ``x``.
    minimum, maximum : float, optional
        Calibrated range of the **raw** characteristic ``x`` in its published
        units, before ``scale`` and ``offset`` are applied. That is how reports
        tabulate ranges.
    scale : float, default 1.0
        Positive, finite multiplier applied to ``x`` before the transform.
    offset : float, default 0.0
        Finite constant added after ``scale`` and before the transform.

    Notes
    -----
    A published covariance matrix ``(X^T Lambda^-1 X)^-1`` is in the basis of the
    *transformed* predictors, ``T(scale * x + offset)``, because that is the
    design matrix the regression was fitted on. An equation stored as published,
    with its rescaling in ``scale`` and ``offset``, therefore takes the covariance
    matrix exactly as printed. Folding a ``log10`` scale into the intercept
    instead (``b * log10(x / k) = b * log10(x) - b * log10(k)``) moves the
    predictor basis, and the intercept's row and column of the covariance matrix
    would have to be transformed to match.
    """

    code: str
    transform: str = "log10"
    units: str = ""
    minimum: Optional[float] = None
    maximum: Optional[float] = None
    scale: float = 1.0
    offset: float = 0.0

    def __post_init__(self) -> None:
        if self.transform not in TRANSFORMS:
            raise ValueError(
                f"Unknown transform {self.transform!r}; choose from {sorted(TRANSFORMS)}"
            )
        if self.minimum is not None and self.maximum is not None and self.minimum > self.maximum:
            raise ValueError(f"{self.code}: minimum {self.minimum} > maximum {self.maximum}")
        if not (math.isfinite(self.scale) and self.scale > 0):
            raise ValueError(f"{self.code}: scale must be positive and finite, got {self.scale}")
        if not math.isfinite(self.offset):
            raise ValueError(f"{self.code}: offset must be finite, got {self.offset}")

    def apply(self, x: float) -> float:
        """Return ``T(scale * x + offset)``, the value entering the regression.

        Parameters
        ----------
        x : float
            The characteristic in its published units.

        Returns
        -------
        float

        Raises
        ------
        ValueError
            If the transform's argument is outside its domain: ``<= 0`` for
            ``log10``, ``<= -1`` for ``log10_plus1``. ``allow_extrapolation``
            cannot override this, because there is no value to extrapolate to.
        """
        arg = self.scale * x + self.offset
        lower = _DOMAIN_LOWER[self.transform]
        if lower is not None and not arg > lower:
            raise ValueError(
                f"{self.code}={x}: {self.transform} needs {self._argument_text()} > "
                f"{lower:g}, got {arg:g}"
            )
        return TRANSFORMS[self.transform](arg)

    def _argument_text(self) -> str:
        text = self.code if self.scale == 1.0 else f"{self.scale:g}*{self.code}"
        if self.offset:
            text += f" {'+' if self.offset > 0 else '-'} {abs(self.offset):g}"
        return text

    def to_dict(self) -> Dict[str, Any]:
        """Return the JSON form, omitting ``scale`` and ``offset`` at their defaults.

        Omitting the defaults keeps a variable that does not rescale identical
        to the form written before ``scale`` and ``offset`` existed.
        """
        d: Dict[str, Any] = {"code": self.code, "transform": self.transform}
        if self.units:
            d["units"] = self.units
        if self.minimum is not None:
            d["minimum"] = self.minimum
        if self.maximum is not None:
            d["maximum"] = self.maximum
        if self.scale != 1.0:
            d["scale"] = self.scale
        if self.offset != 0.0:
            d["offset"] = self.offset
        return d


@dataclass(frozen=True)
class RegressionEquation:
    """One published equation for one AEP in one regression region.

    Parameters
    ----------
    region_code : str
        The report's or NSS's region identifier.
    aep : float
        Annual exceedance probability, strictly between 0 and 1.
    intercept : float
        ``b0`` in log10 space.
    variables : tuple of Variable
        Explanatory variables, in the order the covariance matrix uses.
    coefficients : tuple of float
        ``b_i``, one per variable.
    citation : Citation
        Required.
    sep_log : float, optional
        Standard error of prediction, log10 units.
    avp : float, optional
        Average variance of prediction, log10 units squared.
    model_error_variance : float, optional
        ``gamma**2`` of a GLS fit, log10 units squared.
    covariance : tuple of tuple of float, optional
        ``(X^T Lambda^-1 X)^-1``, of size ``(1 + p) x (1 + p)`` with the intercept first.
    n_sites : int, optional
        Number of gages in the calibration, which sets the t-distribution's degrees of freedom.
    bias_correction : float, default 1.0
        Multiplier applied to the retransformed estimate, when the report prescribes one.
    """

    region_code: str
    aep: float
    intercept: float
    variables: Tuple[Variable, ...]
    coefficients: Tuple[float, ...]
    citation: Citation
    region_name: str = ""
    sep_log: Optional[float] = None
    avp: Optional[float] = None
    model_error_variance: Optional[float] = None
    covariance: Optional[Tuple[Tuple[float, ...], ...]] = None
    n_sites: Optional[int] = None
    bias_correction: float = 1.0

    def __post_init__(self) -> None:
        if not 0.0 < self.aep < 1.0:
            raise ValueError(f"aep must be in (0, 1), got {self.aep}")
        if len(self.variables) != len(self.coefficients):
            raise ValueError("variables and coefficients differ in length")
        codes = [v.code for v in self.variables]
        if len(set(codes)) != len(codes):
            raise ValueError(f"duplicate variable codes: {codes}")
        if self.covariance is not None:
            k = 1 + len(self.variables)
            if len(self.covariance) != k or any(len(r) != k for r in self.covariance):
                raise ValueError(f"covariance must be {k}x{k}")
        if self.bias_correction <= 0:
            raise ValueError("bias_correction must be positive")

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "RegressionEquation":
        """Build from the JSON form stored under ``flowfreq/data/regression/``."""
        cov = d.get("covariance")
        return cls(
            region_code=str(d["region_code"]),
            region_name=str(d.get("region_name", "")),
            aep=float(d["aep"]),
            intercept=float(d["intercept"]),
            variables=tuple(Variable(**v) for v in d["variables"]),
            coefficients=tuple(float(c) for c in d["coefficients"]),
            citation=Citation(**d["citation"]),
            sep_log=d.get("sep_log"),
            avp=d.get("avp"),
            model_error_variance=d.get("model_error_variance"),
            covariance=tuple(tuple(float(x) for x in row) for row in cov) if cov else None,
            n_sites=d.get("n_sites"),
            bias_correction=float(d.get("bias_correction", 1.0)),
        )


@dataclass(frozen=True)
class RegressionEstimate:
    """An evaluated equation, with its provenance."""

    equation: RegressionEquation
    flow_cfs: float
    log_flow: float
    out_of_range: Tuple[str, ...] = ()
    interval: Optional[Tuple[float, float]] = None
    interval_level: Optional[float] = None
    interval_method: Optional[str] = None
    inputs: Dict[str, float] = field(default_factory=dict)

    @property
    def extrapolated(self) -> bool:
        """Whether any input was outside the calibrated range."""
        return bool(self.out_of_range)


def _check_range(eq: RegressionEquation, values: Mapping[str, float]) -> List[str]:
    problems = []
    for v in eq.variables:
        x = values[v.code]
        if v.minimum is not None and x < v.minimum:
            problems.append(f"{v.code}={x} < min {v.minimum}")
        if v.maximum is not None and x > v.maximum:
            problems.append(f"{v.code}={x} > max {v.maximum}")
    return problems


def _prediction_sd(
    eq: RegressionEquation, x_row: np.ndarray
) -> Tuple[Optional[float], Optional[str]]:
    if eq.covariance is not None and eq.model_error_variance is not None:
        u = np.asarray(eq.covariance, dtype=float)
        var = eq.model_error_variance + float(x_row @ u @ x_row)
        return math.sqrt(var), "site-specific (model error + sampling covariance)"
    if eq.avp is not None:
        return math.sqrt(eq.avp), "average variance of prediction"
    if eq.sep_log is not None:
        return eq.sep_log, "average standard error of prediction"
    return None, None


def evaluate(
    eq: RegressionEquation,
    characteristics: Mapping[str, float],
    *,
    allow_extrapolation: bool = False,
    interval_level: float = 0.90,
) -> RegressionEstimate:
    """Evaluate an equation at one site.

    Parameters
    ----------
    eq : RegressionEquation
    characteristics : mapping of str to float
        Basin characteristics keyed by StreamStats code.
    allow_extrapolation : bool, default False
        Compute the value even outside the calibrated range, and record it.
    interval_level : float, default 0.90
        Two-sided prediction-interval level.

    Returns
    -------
    RegressionEstimate

    Raises
    ------
    KeyError
        If a required characteristic is missing.
    OutOfRangeError
        If an input is outside its limits and ``allow_extrapolation`` is False.
    ValueError
        If a transform's argument, ``scale * x + offset``, is outside its domain
        (``<= 0`` for ``log10``), whatever ``allow_extrapolation`` says.
    """
    missing = [v.code for v in eq.variables if v.code not in characteristics]
    if missing:
        raise KeyError(f"Missing basin characteristic(s) for {eq.region_code}: {missing}")
    values = {v.code: float(characteristics[v.code]) for v in eq.variables}
    problems = _check_range(eq, values)
    if problems and not allow_extrapolation:
        raise OutOfRangeError(f"{eq.region_code} AEP {eq.aep}: " + "; ".join(problems))
    if problems:
        logger.warning("Extrapolating %s AEP %s: %s", eq.region_code, eq.aep, "; ".join(problems))

    try:
        transformed = np.array([v.apply(values[v.code]) for v in eq.variables])
    except ValueError as exc:
        raise ValueError(f"{eq.region_code} AEP {eq.aep}: {exc}") from None
    log_q = eq.intercept + float(np.dot(eq.coefficients, transformed))
    flow = eq.bias_correction * 10.0**log_q

    x_row = np.concatenate([[1.0], transformed])
    sd, method = _prediction_sd(eq, x_row)
    interval = None
    if sd is not None:
        alpha = 1.0 - interval_level
        if eq.n_sites is not None and eq.n_sites > len(eq.variables) + 1:
            crit = float(stats.t.ppf(1 - alpha / 2, eq.n_sites - len(eq.variables) - 1))
        else:
            crit = float(stats.norm.ppf(1 - alpha / 2))
        t_factor = 10.0 ** (crit * sd)
        interval = (flow / t_factor, flow * t_factor)
    return RegressionEstimate(
        equation=eq,
        flow_cfs=flow,
        log_flow=log_q,
        out_of_range=tuple(problems),
        interval=interval,
        interval_level=interval_level if interval else None,
        interval_method=method,
        inputs=values,
    )


def evaluate_weighted(
    estimates_by_region: Sequence[Tuple[float, RegressionEstimate]],
    *,
    space: str = "log",
) -> float:
    """Area-weight estimates for a basin that spans several regions.

    Reports differ on whether to weight the flows or their logarithms, so the
    caller must say which: ``space="log"`` or ``space="linear"``.

    Parameters
    ----------
    estimates_by_region : sequence of (area_fraction, RegressionEstimate)
        The area fractions must sum to 1.
    space : {"log", "linear"}

    Returns
    -------
    float
        The weighted flow in cfs.

    Raises
    ------
    ValueError
        If the fractions do not sum to 1, or the space is unknown.
    """
    fracs = np.array([f for f, _ in estimates_by_region], dtype=float)
    if not np.isclose(fracs.sum(), 1.0, atol=1e-6) or (fracs < 0).any():
        raise ValueError(f"area fractions must be non-negative and sum to 1, got {fracs.sum()}")
    if space == "log":
        return float(10.0 ** sum(f * math.log10(e.flow_cfs) for f, e in estimates_by_region))
    if space == "linear":
        return float(sum(f * e.flow_cfs for f, e in estimates_by_region))
    raise ValueError(f"space must be 'log' or 'linear', got {space!r}")
