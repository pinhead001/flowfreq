"""Regional regression equations: schema, evaluation, prediction intervals.

An equation has the log-linear form USGS state reports use::

    log10(Q_p) = b0 + sum_i b_i * T_i(X_i)

where each ``T_i`` is the transform the report applies to basin
characteristic ``X_i`` (``log10``, ``log10(X+1)``, a power ``X**k``, or none),
optionally after a linear rescaling ``scale * X_i + offset`` (see :class:`Variable`). Power-form
equations, ``Q = 10**b0 * A**b1 * ...``, are the same thing written out, so a
published term such as ``(FOREST/100 + 1)**b`` is ``log10`` with ``scale=0.01``
and ``offset=1``.

Two published things are not log-linear, and the schema stores them explicitly:

- **A printed constant.** Nevada's Region 6 prints its 2-year equation as ``Q=0``
  (WSP 2433, Table 10). It is stored with ``form="constant"`` and ``constant_cfs``,
  without an ``intercept``, so a reader that predates the form fails with
  ``KeyError: 'intercept'`` instead of misreading it. Its ``variables`` set the
  calibrated range only, and every coefficient is zero.
- **An error statistic that is not a prediction error.** The hybrid-method regions'
  "estimated average standard error of regression" is, the report says, "not
  comparable" to a regression standard error. It goes in ``published_errors``
  (:class:`PublishedError`), which the prediction-interval code never reads.

See ``docs/NEVADA_REGION6_SCHEMA.md`` for that design.

An equation can also carry **applicability limits** (``applicability``): ranges
on basin characteristics it does not use, outside which it does not apply. A
lake-storage limit is one, where a basin with more lake area than the
calibration's is outside the equations even though lake area is not a term.
They are checked with the variables' own ranges and need the characteristic as
an input.

The schema follows :class:`flowfreq.transpose.RegressionExponents`. An
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

#: Transforms that take a parameter, so cannot live in :data:`TRANSFORMS`:
#: ``power`` is ``x**exponent`` and needs :attr:`Variable.exponent`.
PARAMETRIC_TRANSFORMS = frozenset({"power"})

# Exclusive lower bound on the argument each transform accepts. ``log10_plus1``
# is ``log10`` with its argument shifted by one, so its bound is -1. ``power``
# needs a positive base, because a real power of a non-positive number is
# undefined or not real for a general (negative, fractional) exponent.
_DOMAIN_LOWER: Dict[str, Optional[float]] = {
    "identity": None,
    "log10": 0.0,
    "log10_plus1": -1.0,
    "power": 0.0,
}


#: Equation forms. ``loglinear`` is ``log10 Q = b0 + sum b_i T_i(X_i)``; ``constant``
#: is a printed constant flow (``Q=0``) with no regression behind it.
FORMS = frozenset({"loglinear", "constant"})

#: Kinds of :class:`PublishedError`. Each is a statistic a report prints that is *not*
#: a standard error of prediction, so none of them may set a prediction interval.
PUBLISHED_ERROR_KINDS = frozenset({"hybrid_standard_error_of_regression"})


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
    ``10**(b * A**-c)``        power       1      0       (``exponent=-c``)
    =========================  ==========  =====  ======

    ``transform="power"`` evaluates ``(scale * x + offset) ** exponent``. It is
    for a raw power of a characteristic *inside* the log-linear sum, as in
    Arizona's region 5, ``Q = 10**(a - b * A**-c)`` (SIR 2014-5211, Table 9):
    intercept ``a``, coefficient ``-b``, exponent ``-c``. An ordinary power-form
    term ``A**b`` is still ``log10`` with coefficient ``b``.

    ``transform="log10_plus1"`` is kept as a backward-compatible alias for
    ``log10`` with ``offset`` increased by 1: it evaluates
    ``log10(scale * x + offset + 1)``.

    Parameters
    ----------
    code : str
        StreamStats characteristic code.
    transform : {"identity", "log10", "log10_plus1", "power"}, default "log10"
    units : str
        Published units of ``x``.
    minimum, maximum : float, optional
        Calibrated range of the **raw** characteristic ``x`` in its published
        units, before ``scale`` and ``offset`` are applied. That is how reports
        tabulate ranges.
    scale : float, default 1.0
        Finite, non-zero multiplier applied to ``x`` before the transform. A
        negative scale is NSS's negated base, ``(-1*(x - 13))^b``.
    offset : float, default 0.0
        Finite constant added after ``scale`` and before the transform.
    exponent : float, optional
        Required with ``transform="power"`` (finite and non-zero) and not
        allowed with any other transform.

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
    exponent: Optional[float] = None

    def __post_init__(self) -> None:
        if self.transform not in TRANSFORMS and self.transform not in PARAMETRIC_TRANSFORMS:
            raise ValueError(
                f"Unknown transform {self.transform!r}; choose from "
                f"{sorted(set(TRANSFORMS) | PARAMETRIC_TRANSFORMS)}"
            )
        if self.transform == "power":
            if self.exponent is None or not math.isfinite(self.exponent) or self.exponent == 0:
                raise ValueError(
                    f"{self.code}: transform 'power' needs a finite, non-zero exponent, "
                    f"got {self.exponent}"
                )
        elif self.exponent is not None:
            raise ValueError(
                f"{self.code}: exponent is only used with transform 'power', "
                f"not {self.transform!r}"
            )
        if self.minimum is not None and self.maximum is not None and self.minimum > self.maximum:
            raise ValueError(f"{self.code}: minimum {self.minimum} > maximum {self.maximum}")
        if not (math.isfinite(self.scale) and self.scale != 0):
            raise ValueError(f"{self.code}: scale must be finite and non-zero, got {self.scale}")
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
            ``log10`` and ``power``, ``<= -1`` for ``log10_plus1``. ``allow_extrapolation``
            cannot override this, because there is no value to extrapolate to.
        """
        arg = self.scale * x + self.offset
        lower = _DOMAIN_LOWER[self.transform]
        if lower is not None and not arg > lower:
            bound = "positive" if lower == 0 else f"above {lower:g}"
            raise ValueError(
                f"{self.code}: {self.transform} argument {arg:g} is not {bound}; "
                f"{self.transform} needs {self._argument_text()} > {lower:g} at {self.code}={x}"
            )
        if self.transform == "power":
            assert self.exponent is not None  # enforced in __post_init__
            return float(arg**self.exponent)
        return TRANSFORMS[self.transform](arg)

    def _argument_text(self) -> str:
        text = self.code if self.scale == 1.0 else f"{self.scale:g}*{self.code}"
        if self.offset:
            text += f" {'+' if self.offset > 0 else '-'} {abs(self.offset):g}"
        return text

    def to_dict(self) -> Dict[str, Any]:
        """Return the JSON form, omitting ``scale``, ``offset`` and ``exponent`` when unset.

        Omitting the defaults keeps a variable that does not rescale identical
        to the form written before ``scale``, ``offset`` and ``exponent`` existed.
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
        if self.exponent is not None:
            d["exponent"] = self.exponent
        return d


@dataclass(frozen=True)
class PublishedError:
    """An error statistic a report prints that is not a standard error of prediction.

    It is stored so that everything the report publishes is machine-readable, and it
    is **never** used for a prediction interval: :func:`evaluate` reads only
    ``sep_log``, ``avp`` and the covariance. A statistic that *is* comparable to a
    standard error of prediction belongs in :attr:`RegressionEquation.sep_log`.

    Parameters
    ----------
    kind : str
        One of :data:`PUBLISHED_ERROR_KINDS`.
    value : float
        The printed value, unconverted. It must be finite and non-negative.
    units : str
        The printed units, verbatim (for example ``"log units"``). Required.
    source : str
        Where it is printed, such as ``"Table 10, p. 46"``.
    """

    kind: str
    value: float
    units: str
    source: str = ""

    def __post_init__(self) -> None:
        if self.kind not in PUBLISHED_ERROR_KINDS:
            raise ValueError(
                f"Unknown published-error kind {self.kind!r}; choose from "
                f"{sorted(PUBLISHED_ERROR_KINDS)}"
            )
        if not (math.isfinite(self.value) and self.value >= 0):
            raise ValueError(f"{self.kind}: value must be finite and >= 0, got {self.value}")
        if not self.units.strip():
            raise ValueError(f"{self.kind}: units are required (store them as printed)")

    def to_dict(self) -> Dict[str, Any]:
        """Return the JSON form, omitting ``source`` when empty."""
        d: Dict[str, Any] = {"kind": self.kind, "value": self.value, "units": self.units}
        if self.source:
            d["source"] = self.source
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
    intercept : float or None
        ``b0`` in log10 space. Required for ``form="loglinear"``, and must be ``None``
        for ``form="constant"``.
    variables : tuple of Variable
        Explanatory variables, in the order the covariance matrix uses. For a constant
        equation they set the calibrated range only.
    coefficients : tuple of float
        ``b_i``, one per variable. All zero for a constant equation.
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
    form : {"loglinear", "constant"}, default "loglinear"
        ``constant`` is a printed constant flow, such as Nevada Region 6's ``Q=0``. It
        takes ``constant_cfs``, no intercept, no error statistics and no bias correction.
    constant_cfs : float, optional
        The constant flow, finite and non-negative. Required for ``form="constant"`` and
        not allowed otherwise.
    published_errors : tuple of PublishedError
        Printed error statistics that are not standard errors of prediction. They never
        set an interval.
    equivalent_years : float, optional
        The report's equivalent years of record for the estimate (Hardison, 1971), finite
        and positive. Informational: no evaluator reads it. A report that weights a gaged
        site's station estimate with the regression estimate by years of record (Sauer,
        1974) uses it there, as :func:`flowfreq.regression.nevada.weight_gaged_estimate`
        does. Not allowed on a constant equation.
    applicability : tuple of Variable
        Applicability limits: characteristics the equation does not use, each with a
        ``minimum`` and/or ``maximum`` outside which the equation does not apply. They
        carry no coefficient and no covariance entry; :func:`evaluate` range-checks
        them like the variables, so each must be supplied. A code may not repeat a
        variable's (set that variable's range instead).
    """

    region_code: str
    aep: float
    intercept: Optional[float]
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
    form: str = "loglinear"
    constant_cfs: Optional[float] = None
    published_errors: Tuple[PublishedError, ...] = ()
    equivalent_years: Optional[float] = None
    applicability: Tuple[Variable, ...] = ()

    def __post_init__(self) -> None:
        if self.form not in FORMS:
            raise ValueError(f"Unknown form {self.form!r}; choose from {sorted(FORMS)}")
        if self.form == "constant":
            self._check_constant()
        else:
            if self.intercept is None or not math.isfinite(self.intercept):
                raise ValueError(
                    f"a loglinear equation needs a finite intercept, got {self.intercept}"
                )
            if self.constant_cfs is not None:
                raise ValueError("constant_cfs is only used with form 'constant'")
        if not 0.0 < self.aep < 1.0:
            raise ValueError(f"aep must be in (0, 1), got {self.aep}")
        if len(self.variables) != len(self.coefficients):
            raise ValueError("variables and coefficients differ in length")
        codes = [v.code for v in self.variables]
        if len(set(codes)) != len(codes):
            raise ValueError(f"duplicate variable codes: {codes}")
        limit_codes = [v.code for v in self.applicability]
        if len(set(limit_codes)) != len(limit_codes) or set(limit_codes) & set(codes):
            raise ValueError(
                f"applicability limits {limit_codes} repeat a code or a variable's ({codes})"
            )
        for v in self.applicability:
            if v.minimum is None and v.maximum is None:
                raise ValueError(f"applicability limit {v.code} has no minimum or maximum")
        if self.covariance is not None:
            k = 1 + len(self.variables)
            if len(self.covariance) != k or any(len(r) != k for r in self.covariance):
                raise ValueError(f"covariance must be {k}x{k}")
        if self.bias_correction <= 0:
            raise ValueError("bias_correction must be positive")
        ey = self.equivalent_years
        if ey is not None and not (math.isfinite(ey) and ey > 0):
            raise ValueError(f"equivalent_years must be finite and positive, got {ey}")

    def _check_constant(self) -> None:
        where = f"{self.region_code} AEP {self.aep}"
        c = self.constant_cfs
        if c is None or not (math.isfinite(c) and c >= 0):
            raise ValueError(f"{where}: form 'constant' needs a finite constant_cfs >= 0, got {c}")
        if self.intercept is not None:
            raise ValueError(f"{where}: a constant equation has no intercept")
        if any(b != 0 for b in self.coefficients):
            raise ValueError(f"{where}: a constant equation's coefficients must all be zero")
        stats_set = [
            name
            for name in (
                "sep_log",
                "avp",
                "model_error_variance",
                "covariance",
                "equivalent_years",
            )
            if getattr(self, name) is not None
        ]
        if stats_set:
            raise ValueError(f"{where}: a constant equation has no error statistics: {stats_set}")
        if self.bias_correction != 1.0:
            raise ValueError(f"{where}: a constant equation takes no bias correction")

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "RegressionEquation":
        """Build from the JSON form stored under ``flowfreq/data/regression/``.

        A ``form="constant"`` entry has no ``intercept`` key (so a reader that predates
        the form raises ``KeyError: 'intercept'`` rather than misreading it) and may omit
        ``coefficients``, which are then zeros. A ``loglinear`` entry without an
        ``intercept`` raises ``KeyError``.
        """
        cov = d.get("covariance")
        form = str(d.get("form", "loglinear"))
        variables = tuple(Variable(**v) for v in d["variables"])
        intercept: Optional[float]
        if form == "constant":
            if "intercept" in d:
                raise ValueError(
                    f"{d['region_code']} AEP {d['aep']}: a constant equation has no intercept"
                )
            intercept = None
            coefficients = tuple(float(c) for c in d.get("coefficients", [0.0] * len(variables)))
        else:
            intercept = float(d["intercept"])
            coefficients = tuple(float(c) for c in d["coefficients"])
        constant = d.get("constant_cfs")
        ey = d.get("equivalent_years")
        return cls(
            region_code=str(d["region_code"]),
            region_name=str(d.get("region_name", "")),
            aep=float(d["aep"]),
            intercept=intercept,
            variables=variables,
            coefficients=coefficients,
            citation=Citation(**d["citation"]),
            sep_log=d.get("sep_log"),
            avp=d.get("avp"),
            model_error_variance=d.get("model_error_variance"),
            covariance=tuple(tuple(float(x) for x in row) for row in cov) if cov else None,
            n_sites=d.get("n_sites"),
            bias_correction=float(d.get("bias_correction", 1.0)),
            form=form,
            constant_cfs=None if constant is None else float(constant),
            published_errors=tuple(PublishedError(**p) for p in d.get("published_errors", [])),
            equivalent_years=None if ey is None else float(ey),
            applicability=tuple(Variable(**v) for v in d.get("applicability", [])),
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
    #: Variance of prediction, log10 units squared: ``interval_method``'s quantity,
    #: squared. ``None`` when the equation publishes no standard error of prediction.
    prediction_variance: Optional[float] = None

    @property
    def extrapolated(self) -> bool:
        """Whether any input was outside the calibrated range."""
        return bool(self.out_of_range)


def _check_range(eq: RegressionEquation, values: Mapping[str, float]) -> List[str]:
    problems = []
    for v in (*eq.variables, *eq.applicability):
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
        Basin characteristics keyed by StreamStats code, including those of the
        equation's ``applicability`` limits.
    allow_extrapolation : bool, default False
        Compute the value even outside the calibrated range, and record it.
    interval_level : float, default 0.90
        Two-sided prediction-interval level.

    Returns
    -------
    RegressionEstimate

    Notes
    -----
    A ``form="constant"`` equation is still range-checked against its variables, and
    then returns ``constant_cfs``, with ``log_flow = log10(constant_cfs)`` (``-inf``
    for a printed ``Q=0``) and no interval. An equation whose only error statistics
    are :class:`PublishedError` entries also gets no interval, and
    ``interval_method`` says why.

    Raises
    ------
    KeyError
        If a required characteristic is missing.
    OutOfRangeError
        If an input is outside its limits and ``allow_extrapolation`` is False.
    ValueError
        If a transform's argument, ``scale * x + offset``, is outside its domain
        (``<= 0`` for ``log10`` and ``power``), whatever ``allow_extrapolation`` says.
    """
    needed = (*eq.variables, *eq.applicability)
    missing = [v.code for v in needed if v.code not in characteristics]
    if missing:
        raise KeyError(f"Missing basin characteristic(s) for {eq.region_code}: {missing}")
    values = {v.code: float(characteristics[v.code]) for v in eq.variables}
    problems = _check_range(eq, {v.code: float(characteristics[v.code]) for v in needed})
    if problems and not allow_extrapolation:
        raise OutOfRangeError(f"{eq.region_code} AEP {eq.aep}: " + "; ".join(problems))
    if problems:
        logger.warning("Extrapolating %s AEP %s: %s", eq.region_code, eq.aep, "; ".join(problems))

    try:
        transformed = np.array([v.apply(values[v.code]) for v in eq.variables])
    except ValueError as exc:
        raise ValueError(f"{eq.region_code} AEP {eq.aep}: {exc}") from None
    if eq.form == "constant":
        assert eq.constant_cfs is not None  # enforced in __post_init__
        return RegressionEstimate(
            equation=eq,
            flow_cfs=eq.constant_cfs,
            log_flow=math.log10(eq.constant_cfs) if eq.constant_cfs > 0 else -math.inf,
            out_of_range=tuple(problems),
            interval_method=f"none: printed constant (Q={eq.constant_cfs:g})",
            inputs=values,
        )
    assert eq.intercept is not None  # enforced in __post_init__
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
    elif eq.published_errors:
        kinds = ", ".join(sorted({p.kind for p in eq.published_errors}))
        method = f"none: the published error ({kinds}) is not a standard error of prediction"
    return RegressionEstimate(
        equation=eq,
        flow_cfs=flow,
        log_flow=log_q,
        out_of_range=tuple(problems),
        interval=interval,
        interval_level=interval_level if interval else None,
        interval_method=method,
        inputs=values,
        prediction_variance=None if sd is None else sd**2,
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
        If the fractions do not sum to 1, or the space is unknown, or ``space="log"``
        and a component's flow is zero (a printed ``Q=0``): a logarithm of zero has no
        weighted mean, so weight such estimates in linear space.
    """
    fracs = np.array([f for f, _ in estimates_by_region], dtype=float)
    if not np.isclose(fracs.sum(), 1.0, atol=1e-6) or (fracs < 0).any():
        raise ValueError(f"area fractions must be non-negative and sum to 1, got {fracs.sum()}")
    if space == "log":
        zero = [e for _, e in estimates_by_region if not e.flow_cfs > 0]
        if zero:
            names = ", ".join(f"{e.equation.region_code} AEP {e.equation.aep}" for e in zero)
            raise ValueError(
                f"cannot weight in log space: {names} has flow {zero[0].flow_cfs:g} cfs; "
                "weight in linear space"
            )
        return float(10.0 ** sum(f * math.log10(e.flow_cfs) for f, e in estimates_by_region))
    if space == "linear":
        return float(sum(f * e.flow_cfs for f, e in estimates_by_region))
    raise ValueError(f"space must be 'log' or 'linear', got {space!r}")
