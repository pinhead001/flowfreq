"""
Variance-based weighting of gage transfers against a regional regression.

A drainage-area transfer from a donor gage is usually allowed inside a fixed
area-ratio band (0.5-1.5) and the regression used outside it. The band has no
stated basis (Ries and Friesz, 2000). This module replaces it with a measured
error model. A transfer's variance is

    V_T = V_i + L^2 Var(b) + s2_d,        L = log10(A_t / A_i)

and ``s2_d``, the dissimilarity variance, is fitted from flow-connected gage
pairs. Transfers and the regression are then combined by generalized least
squares on their covariance. The area ratio at which ``V_T`` equals the
regression's variance becomes a derived band, with a bootstrap interval.

All quantities are log10 discharges at one AEP. See
``docs/TRANSPOSITION_DESIGN.md`` section 11 for the specification.

The pair fit
------------
For a connected pair with at-site estimates ``y_i`` and ``y_j``::

    e_ij   = y_i + b L_ij - y_j
    E[e^2] = V_i + V_j - 2 C_ij + L^2 Var(b) + s2_d

The receiving gage's own sampling variance ``V_j`` is subtracted because the
"truth" there is also an estimate. Nested gages share most of their flood
years, so the sampling covariance ``C_ij`` is subtracted too
(:func:`flowfreq.regression.gls.sampling_covariance`, labelled by
:func:`sampling_covariance_frame`). ``s2_d = s0 + sum_k s_k x_k^2`` is fitted
to ``e^2 - known`` by iterated weighted least squares under ``s_k >= 0``.

What this module cannot do
--------------------------
- **It does not find the nesting.** Connected pairs are an input.
- **Correlated errors are approximated.** The sampling covariance follows
  WREG's distance model unless the caller supplies a better matrix. Two
  donors' dissimilarity errors share one pooled correlation. The transfer and
  the regression are correlated near a gage that calibrated the regression,
  and that correlation is a scalar the caller sets.
- **``Var(b)`` is treated as independent per pair.** It is common to every
  pair, but it is small (``<= 0.0025`` against ``s2_d`` of order 0.1).
- **One AEP per call.**
- **Log space throughout.** ``10 ** log_q`` is a median.

References
----------
Ries, K.G., III, and Friesz, P.J., 2000, Methods for estimating low-flow
statistics for Massachusetts streams: U.S. Geological Survey Water-Resources
Investigations Report 00-4135.

Eng, K., Chen, Y.-Y., and Kiang, J.E., 2009, User's guide to the weighted-
multiple-linear-regression program (WREG version 1.0): U.S. Geological Survey
Techniques and Methods 4-A8.
"""

from __future__ import annotations

import logging
import math
from collections import Counter
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import norm

from flowfreq.regression.gls import WregCorrelation, sampling_covariance
from flowfreq.transpose import DEFAULT_B_BOUNDS, RegressionExponents, _require_kind

logger = logging.getLogger(__name__)

#: Name of the intercept ``s0`` in :attr:`DissimilarityModel.coefficients`.
CONST = "const"

#: Label of the regression row and column in :func:`combine`'s output.
REGRESSION = "regression"

_MAX_ITER = 50
_REL_TOL = 1e-8
#: Floor on a pair's expected squared error when weighting, so a pair whose
#: known variance happens to be near zero cannot take all the weight.
_MIN_EXPECTED = 1e-8


# --------------------------------------------------------------------------- inputs


def sampling_covariance_frame(
    gages: pd.DataFrame,
    aep: float,
    skew_mse: float,
    corr: WregCorrelation,
    concurrent: pd.DataFrame,
) -> pd.DataFrame:
    """:func:`~flowfreq.regression.gls.sampling_covariance`, labelled by site.

    Parameters
    ----------
    gages : pd.DataFrame
        Columns ``site_no``, ``sd``, ``skew``, ``n``, ``lat`` and ``lon``.
    aep, skew_mse : float
        As :func:`~flowfreq.regression.gls.lp3_quantile_sampling_variance`.
    corr : WregCorrelation
        Correlation-distance model.
    concurrent : pd.DataFrame
        Concurrent years of record, indexed and columned by ``site_no``. It is
        aligned to ``gages`` by label.

    Returns
    -------
    pd.DataFrame
        Sampling covariance, indexed and columned by ``site_no`` in ``gages`` order.

    Raises
    ------
    KeyError
        A site in ``gages`` is missing from ``concurrent``.
    """
    sites = gages["site_no"].astype(str).tolist()
    c = concurrent.rename(index=str, columns=str)
    missing = sorted(set(sites) - set(c.index) | set(sites) - set(c.columns))
    if missing:
        raise KeyError(f"concurrent has no row or column for {missing}")
    sigma = sampling_covariance(
        gages["sd"],
        gages["skew"],
        gages["n"],
        aep,
        skew_mse,
        gages["lat"],
        gages["lon"],
        c.loc[sites, sites].to_numpy(dtype=float),
        corr,
    )
    return pd.DataFrame(sigma, index=sites, columns=sites)


def _exponent(exponents: RegressionExponents, aep: float, caller: str) -> float:
    _require_kind(exponents, "aep", caller)
    frame = exponents.exponent_at([aep], b_bounds=DEFAULT_B_BOUNDS)
    return float(frame["exponent"].iloc[0])


def _column(gages: pd.DataFrame, sites: List[str], name: str) -> Dict[str, float]:
    return dict(zip(sites, gages[name].to_numpy(dtype=float).tolist()))


class _Network:
    """The gage table and its sampling covariance, as float lookups by site."""

    def __init__(
        self,
        gages: pd.DataFrame,
        sampling_cov: Optional[pd.DataFrame],
        attributes: Sequence[str],
    ) -> None:
        self.sites: List[str] = gages["site_no"].astype(str).tolist()
        self.attributes = tuple(attributes)
        self.area = _column(gages, self.sites, "drainage_area_sqmi")
        self.log_q = _column(gages, self.sites, "log_q")
        self.attr = {a: _column(gages, self.sites, a) for a in self.attributes}
        if sampling_cov is not None:
            cov = sampling_cov.rename(index=str, columns=str)
            missing = sorted(set(self.sites) - set(cov.index) | set(self.sites) - set(cov.columns))
            if missing:
                raise KeyError(f"sampling_cov has no row or column for {missing}")
            matrix = cov.loc[self.sites, self.sites].to_numpy(dtype=float)
        elif "variance" in gages.columns:
            matrix = np.diag(gages["variance"].to_numpy(dtype=float))
        else:
            raise ValueError(
                "the gage table needs a 'variance' column (each at-site estimate's sampling "
                "variance) or a sampling_cov frame"
            )
        self._pos = {s: n for n, s in enumerate(self.sites)}
        self._cov = matrix

    def cov(self, a: str, b: str) -> float:
        return float(self._cov[self._pos[a], self._pos[b]])

    def cov_frame(self, sites: Sequence[str]) -> pd.DataFrame:
        idx = [self._pos[s] for s in sites]
        return pd.DataFrame(self._cov[np.ix_(idx, idx)], index=list(sites), columns=list(sites))

    def features(self, donor: str, receiver: str) -> Dict[str, float]:
        """``log_ar`` and each ``d_<attribute>`` for a transfer donor -> receiver."""
        x = {"log_ar": math.log10(self.area[receiver] / self.area[donor])}
        for a in self.attributes:
            x[f"d_{a}"] = self.attr[a][receiver] - self.attr[a][donor]
        return x


def _check_connections(
    connections: Sequence[Tuple[str, str]], sites: Sequence[str]
) -> List[Tuple[str, str]]:
    known = set(sites)
    seen = set()
    out = []
    for a, b in connections:
        a, b = str(a), str(b)
        for s in (a, b):
            if s not in known:
                raise KeyError(f"connection names {s}, which is not in the gage table")
        if a == b:
            raise ValueError(f"gage {a} is connected to itself")
        key = frozenset((a, b))
        if key in seen:
            raise ValueError(
                f"the pair ({a}, {b}) is listed more than once. Each connection enters the "
                "fit once: the reversed pair has the same squared error."
            )
        seen.add(key)
        out.append((a, b))
    return out


def build_pairs(
    gages: pd.DataFrame,
    connections: Sequence[Tuple[str, str]],
    exponents: RegressionExponents,
    aep: float,
    *,
    exponent_variance: float = 0.0,
    sampling_cov: Optional[pd.DataFrame] = None,
    attributes: Sequence[str] = (),
) -> pd.DataFrame:
    """The transfer error of every connected gage pair, at one AEP.

    Parameters
    ----------
    gages : pd.DataFrame
        Columns ``site_no``, ``drainage_area_sqmi`` and ``log_q`` (the at-site
        log10 quantile at ``aep``), plus ``variance`` when ``sampling_cov`` is
        not given, plus every column named in ``attributes``.
    connections : sequence of (str, str)
        Flow-connected pairs as ``(donor, receiver)``. The direction only sets
        the sign of ``error`` and ``log_ar``; each pair is listed once.
    exponents : RegressionExponents
        Drainage-area exponents indexed by AEP.
    aep : float
        The annual exceedance probability ``log_q`` is at.
    exponent_variance : float, default 0.0
        ``Var(b)`` at ``aep``.
    sampling_cov : pd.DataFrame, optional
        Sampling covariance of the at-site estimates, indexed and columned by
        ``site_no`` (e.g. :func:`sampling_covariance_frame`). Without it the
        ``variance`` column is used and the covariances are taken as zero.
    attributes : sequence of str
        Gage columns whose receiver-minus-donor difference is carried as
        ``d_<attribute>`` (e.g. mean annual precipitation).

    Returns
    -------
    pd.DataFrame
        One row per connection: ``donor``, ``receiver``, ``area_ratio``,
        ``log_ar``, ``exponent``, ``transfer``, ``observed``, ``error``,
        ``known_variance`` (``V_i + V_j - 2 C_ij + L^2 Var(b)``), and the
        ``d_`` columns.

    Raises
    ------
    KeyError
        A connection or ``sampling_cov`` names an unknown site.
    ValueError
        A gage is connected to itself, a pair is listed twice, the exponent
        set is not indexed by AEP, ``exponent_variance`` is negative, or no
        variance is available.
    """
    if exponent_variance < 0:
        raise ValueError(f"exponent_variance must be >= 0, got {exponent_variance}")
    b = _exponent(exponents, aep, "build_pairs")
    net = _Network(gages, sampling_cov, attributes)
    return _pairs(net, _check_connections(connections, net.sites), b, exponent_variance)


def _pairs(
    net: _Network, conns: Sequence[Tuple[str, str]], b: float, exponent_variance: float
) -> pd.DataFrame:
    rows = []
    for donor, receiver in conns:
        x = net.features(donor, receiver)
        log_ar = x["log_ar"]
        transfer = net.log_q[donor] + b * log_ar
        observed = net.log_q[receiver]
        known = (
            net.cov(donor, donor)
            + net.cov(receiver, receiver)
            - 2.0 * net.cov(donor, receiver)
            + log_ar**2 * exponent_variance
        )
        row: Dict[str, object] = {
            "donor": donor,
            "receiver": receiver,
            "area_ratio": 10.0**log_ar,
            "log_ar": log_ar,
            "exponent": b,
            "transfer": transfer,
            "observed": observed,
            "error": transfer - observed,
            "known_variance": known,
        }
        row.update({k: v for k, v in x.items() if k != "log_ar"})
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- fit


@dataclass(frozen=True)
class DissimilarityModel:
    """A fitted dissimilarity variance ``s2_d = s0 + sum_k s_k x_k^2``.

    Attributes
    ----------
    terms : tuple of str
        The pair-table columns ``x_k`` (``"log_ar"``, ``"d_precip"``, ...).
    coefficients : dict of str to float
        ``s0`` under :data:`CONST` and each ``s_k`` under its term. All ``>= 0``.
    bound_terms : tuple of str
        Coefficients held at zero by the non-negativity constraint.
    n_pairs, n_gages : int
        Pairs fitted, and distinct gages among them.
    iterations : int
        Reweighting iterations.
    converged : bool
        Whether the coefficients settled within the iteration limit.
    residuals : pd.DataFrame
        Per pair: ``donor``, ``receiver``, ``excess`` (``e^2 - known``),
        ``fitted`` (``s2_d``), ``expected`` (``known + s2_d``) and
        ``standardized`` (``e / sqrt(expected)``, mean square about 1 when the
        model is right).
    by_ratio : pd.DataFrame
        A model-free check: pairs binned by ``|log_ar|``, with
        ``abs_log_ar_mid``, ``n_pairs``, ``mean_excess`` and ``mean_fitted``.
    """

    terms: Tuple[str, ...]
    coefficients: Dict[str, float]
    bound_terms: Tuple[str, ...]
    n_pairs: int
    n_gages: int
    iterations: int
    converged: bool
    residuals: pd.DataFrame
    by_ratio: pd.DataFrame

    def variance(self, **features: float) -> float:
        """``s2_d`` for one transfer.

        Parameters
        ----------
        **features : float
            A value for every term, signed or not (each enters squared).

        Raises
        ------
        KeyError
            A term has no value.
        """
        missing = [t for t in self.terms if t not in features]
        if missing:
            raise KeyError(f"no value for dissimilarity term(s) {missing}")
        return float(
            self.coefficients[CONST]
            + sum(self.coefficients[t] * float(features[t]) ** 2 for t in self.terms)
        )


def _design(pairs: pd.DataFrame, terms: Sequence[str]) -> np.ndarray:
    cols = [np.ones(len(pairs))] + [pairs[t].to_numpy(dtype=float) ** 2 for t in terms]
    return np.column_stack(cols)


def _nnls_irls(
    F: np.ndarray, z: np.ndarray, known: np.ndarray, max_iter: int
) -> Tuple[np.ndarray, int, bool]:
    """Weighted NNLS of ``z`` on ``F``, reweighted by ``1 / (known + F s)^2``."""
    start = max(float(np.mean(z)), 0.0)
    expected = np.maximum(known + start, _MIN_EXPECTED)
    coef = np.zeros(F.shape[1])
    for it in range(1, max_iter + 1):
        w = 1.0 / expected
        new, _ = nnls(F * w[:, None], z * w)
        settled = np.allclose(
            new, coef, rtol=_REL_TOL, atol=_REL_TOL * max(1.0, float(np.max(new)))
        )
        coef = new
        expected = np.maximum(known + F @ coef, _MIN_EXPECTED)
        if settled:
            return coef, it, True
    return coef, max_iter, False


def _by_ratio(
    pairs: pd.DataFrame, excess: np.ndarray, fitted: np.ndarray, n_bins: int
) -> pd.DataFrame:
    a = np.abs(pairs["log_ar"].to_numpy(dtype=float))
    edges = np.unique(np.quantile(a, np.linspace(0.0, 1.0, n_bins + 1)))
    if edges.size < 2:
        edges = np.array([a.min(), a.max() + 1e-12])
    idx = np.clip(np.searchsorted(edges, a, side="right") - 1, 0, edges.size - 2)
    rows = []
    for k in range(edges.size - 1):
        m = idx == k
        if not m.any():
            continue
        rows.append(
            {
                "abs_log_ar_mid": float((edges[k] + edges[k + 1]) / 2.0),
                "n_pairs": int(m.sum()),
                "mean_excess": float(excess[m].mean()),
                "mean_fitted": float(fitted[m].mean()),
            }
        )
    return pd.DataFrame(rows)


def fit_dissimilarity(
    pairs: pd.DataFrame,
    terms: Sequence[str] = ("log_ar",),
    *,
    max_iter: int = _MAX_ITER,
    n_bins: int = 5,
) -> DissimilarityModel:
    """Fit ``s2_d = s0 + sum_k s_k x_k^2`` to the pairs' squared errors.

    Weighted least squares on ``e^2 - known_variance``, weighted by
    ``1 / E[e^2]^2`` (for normal errors ``Var(e^2) = 2 E[e^2]^2``) and
    reweighted until the coefficients settle, with every coefficient
    constrained ``>= 0``. A coefficient held at zero is reported in
    ``bound_terms`` and logged: it means the known variance already explains
    more than that term's share of the scatter.

    Parameters
    ----------
    pairs : pd.DataFrame
        From :func:`build_pairs`.
    terms : sequence of str, default ("log_ar",)
        Pair columns entering squared, after the intercept.
    max_iter : int, default 50
        Reweighting iterations.
    n_bins : int, default 5
        Bins of ``|log_ar|`` for the model-free check.

    Returns
    -------
    DissimilarityModel

    Raises
    ------
    KeyError
        A term, or ``error`` / ``known_variance``, is not a column.
    ValueError
        Fewer pairs than coefficients plus one.
    """
    return _fit(pairs, tuple(terms), max_iter=max_iter, n_bins=n_bins, warn=True)


def _fit(
    pairs: pd.DataFrame,
    terms: Tuple[str, ...],
    *,
    max_iter: int = _MAX_ITER,
    n_bins: int = 5,
    warn: bool,
) -> DissimilarityModel:
    """:func:`fit_dissimilarity`; ``warn=False`` for refits inside a resampling loop."""
    needed = ["error", "known_variance", "log_ar", *terms]
    missing = [c for c in needed if c not in pairs.columns]
    if missing:
        raise KeyError(f"pair table has no column(s) {missing}")
    p = len(terms) + 1
    if len(pairs) <= p:
        raise ValueError(f"need more pairs than coefficients: {len(pairs)} pairs, {p} coefficients")

    e = pairs["error"].to_numpy(dtype=float)
    known = pairs["known_variance"].to_numpy(dtype=float)
    z = e**2 - known
    F = _design(pairs, terms)
    coef, iterations, converged = _nnls_irls(F, z, known, max_iter)
    if warn and not converged:
        logger.warning("fit_dissimilarity did not settle in %d iterations", max_iter)

    names = (CONST, *terms)
    bound = tuple(n for n, c in zip(names, coef) if c == 0.0)
    if warn and bound:
        logger.warning(
            "dissimilarity coefficient(s) %s held at zero: the gages' known sampling "
            "variance explains at least that share of the pairs' squared errors",
            ", ".join(bound),
        )
    fitted = F @ coef
    expected = np.maximum(known + fitted, _MIN_EXPECTED)
    residuals = pd.DataFrame(
        {
            "donor": pairs["donor"].to_numpy() if "donor" in pairs else None,
            "receiver": pairs["receiver"].to_numpy() if "receiver" in pairs else None,
            "excess": z,
            "fitted": fitted,
            "expected": expected,
            "standardized": e / np.sqrt(expected),
        }
    )
    gage_ids: set[str] = set()
    for col in ("donor", "receiver"):
        if col in pairs:
            gage_ids.update(pairs[col].astype(str))
    return DissimilarityModel(
        terms=terms,
        coefficients={n: float(c) for n, c in zip(names, coef)},
        bound_terms=bound,
        n_pairs=len(pairs),
        n_gages=len(gage_ids),
        iterations=iterations,
        converged=converged,
        residuals=residuals,
        by_ratio=_by_ratio(pairs, z, fitted, n_bins),
    )


def estimate_donor_correlation(
    gages: pd.DataFrame,
    connections: Sequence[Tuple[str, str]],
    exponents: RegressionExponents,
    aep: float,
    model: DissimilarityModel,
    *,
    exponent_variance: float = 0.0,
    sampling_cov: Optional[pd.DataFrame] = None,
    attributes: Sequence[str] = (),
) -> Tuple[Optional[float], int]:
    """Pooled correlation of two donors' dissimilarity errors at one receiver.

    For every gage ``j`` connected to two others ``i`` and ``k``::

        E[e_ij e_kj] = C_ik - C_ij - C_kj + V_j + L_ij L_kj Var(b) + rho s_ij s_kj

    and ``rho`` is the pooled moment estimate
    ``sum(e_ij e_kj - known) / sum(s_ij s_kj)``, clipped to ``[-1, 1]``.

    The two errors share the receiver's own deviation, so under the model
    ``rho`` is not negative, and a negative estimate is sampling noise. It is
    returned as estimated. :func:`combine` refuses ``rho <= -1/(k - 1)`` with
    ``k`` donors, and :func:`leave_one_out` floors its fold estimates at 0.

    Parameters
    ----------
    gages, connections, exponents, aep, exponent_variance, sampling_cov, attributes
        As :func:`build_pairs`. ``attributes`` must supply the model's ``d_`` terms.
    model : DissimilarityModel
        The fitted ``s2_d``.

    Returns
    -------
    (float or None, int)
        ``rho``, and how many donor pairs it pools. ``None`` when no gage has
        two connections, or the model gives no dissimilarity variance.
    """
    b = _exponent(exponents, aep, "estimate_donor_correlation")
    net = _Network(gages, sampling_cov, attributes)
    tri = _triples(net, _check_connections(connections, net.sites), b, exponent_variance)
    return _rho_from_triples(tri, model, np.ones(len(tri["product"]), dtype=bool))


def _neighbours(conns: Sequence[Tuple[str, str]]) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for a, c in conns:
        out.setdefault(a, []).append(c)
        out.setdefault(c, []).append(a)
    return out


def _triples(
    net: _Network, conns: Sequence[Tuple[str, str]], b: float, exponent_variance: float
) -> Dict[str, np.ndarray]:
    """Every (donor i, donor k, receiver j) with both donors connected to j.

    Returns arrays: ``product`` (``e_ij e_kj - known``), the three site
    labels, and for each feature ``x`` its value for both transfers as
    ``<x>_i`` and ``<x>_k``. Computed once, so :func:`leave_one_out` can mask
    folds without rebuilding it.
    """
    feats = ["log_ar", *(f"d_{a}" for a in net.attributes)]
    num: Dict[str, List[float]] = {"product": []}
    labels: Dict[str, List[str]] = {"i": [], "k": [], "j": []}
    for f in feats:
        num[f"{f}_i"] = []
        num[f"{f}_k"] = []
    for j, nbrs in _neighbours(conns).items():
        stats = {}
        for i in nbrs:
            x = net.features(i, j)
            stats[i] = (net.log_q[i] + b * x["log_ar"] - net.log_q[j], x)
        for a in range(len(nbrs)):
            for c in range(a + 1, len(nbrs)):
                i, k = nbrs[a], nbrs[c]
                (ei, xi), (ek, xk) = stats[i], stats[k]
                known = (
                    net.cov(i, k)
                    - net.cov(i, j)
                    - net.cov(k, j)
                    + net.cov(j, j)
                    + xi["log_ar"] * xk["log_ar"] * exponent_variance
                )
                num["product"].append(ei * ek - known)
                labels["i"].append(i)
                labels["k"].append(k)
                labels["j"].append(j)
                for f in feats:
                    num[f"{f}_i"].append(xi[f])
                    num[f"{f}_k"].append(xk[f])
    out = {key: np.asarray(v, dtype=float) for key, v in num.items()}
    out.update({key: np.asarray(v, dtype=object) for key, v in labels.items()})
    return out


def _rho_from_triples(
    tri: Dict[str, np.ndarray], model: DissimilarityModel, mask: np.ndarray
) -> Tuple[Optional[float], int]:
    n = int(mask.sum())
    if n == 0:
        return None, 0

    def s2(side: str) -> np.ndarray:
        missing = [t for t in model.terms if f"{t}_{side}" not in tri]
        if missing:
            raise KeyError(f"no value for dissimilarity term(s) {missing}")
        total = np.full(n, model.coefficients[CONST])
        for t in model.terms:
            total = total + model.coefficients[t] * tri[f"{t}_{side}"][mask].astype(float) ** 2
        return total

    den = float(np.sum(np.sqrt(s2("i") * s2("k"))))
    if den <= 0.0:
        return None, n
    num = float(np.sum(tri["product"][mask].astype(float)))
    return float(np.clip(num / den, -1.0, 1.0)), n


# --------------------------------------------------------------------------- blend


@dataclass(frozen=True)
class BlendResult:
    """A GLS combination of transfers and the regression at one target.

    Attributes
    ----------
    log_q : float
        Combined log10 estimate.
    variance : float
        Its variance, ``1 / (1' C^-1 1)``.
    level : float
        Two-sided interval level.
    interval : (float, float)
        ``10 ** (log_q -/+ z sqrt(variance))``, cfs.
    weights : pd.Series
        Per component (donor labels, then :data:`REGRESSION`). They sum to 1
        and can be negative when components are strongly correlated.
    covariance : pd.DataFrame
        The components' error covariance ``C``.
    components : pd.DataFrame
        ``estimate`` and ``variance`` per component.
    """

    log_q: float
    variance: float
    level: float
    interval: Tuple[float, float]
    weights: pd.Series
    covariance: pd.DataFrame
    components: pd.DataFrame

    @property
    def flow_cfs(self) -> float:
        """``10 ** log_q`` (a median in flow space)."""
        return float(10.0**self.log_q)


def _check_rho(name: str, rho: float) -> None:
    if not -1.0 <= rho <= 1.0:
        raise ValueError(f"{name} correlation must be in [-1, 1], got {rho}")


def _check_donor_bound(rho: float, k: int) -> None:
    # rho * s s' + (1 - rho) diag(s^2) is S [(1 - rho) I + rho 1 1'] S, positive
    # definite exactly when rho > -1/(k - 1). The other blocks are covariances
    # and cannot repair it.
    if k >= 2 and rho <= -1.0 / (k - 1):
        raise ValueError(
            f"donor correlation {rho:.4g} is not above -1/(k - 1) = {-1.0 / (k - 1):.4f} "
            f"for {k} donors with dissimilarity variance, so their errors have no "
            "positive definite covariance"
        )


def combine(
    model: DissimilarityModel,
    transfers: pd.DataFrame,
    regression: Optional[Tuple[float, float]],
    *,
    exponent_variance: float = 0.0,
    donor_cov: Optional[pd.DataFrame] = None,
    donor_correlation: float = 0.0,
    regression_correlation: float = 0.0,
    level: float = 0.90,
) -> BlendResult:
    """Combine one target's transfers, and its regression, by GLS.

    Covariance of the components' errors::

        C_ii = V_i + L_i^2 Var(b) + s2_d,i
        C_ik = Cov_ik + L_i L_k Var(b) + donor_correlation * sqrt(s2_d,i s2_d,k)
        C_iR = regression_correlation * sqrt(C_ii V_R)

    Parameters
    ----------
    model : DissimilarityModel
    transfers : pd.DataFrame
        One row per donor: ``donor``, ``transfer`` (log10), ``donor_variance``,
        ``log_ar`` (``log10(A_target / A_donor)``) and every other model term.
    regression : (float, float) or None
        The regression's log10 estimate at the target and its prediction
        variance ``V_R``; ``None`` for transfers alone.
    exponent_variance : float, default 0.0
        ``Var(b)``.
    donor_cov : pd.DataFrame, optional
        Covariance of the donors' estimates, labelled by donor. Its diagonal
        replaces ``donor_variance``. Without it, donors' estimates are
        uncorrelated.
    donor_correlation : float, default 0.0
        Correlation of two donors' dissimilarity errors
        (:func:`estimate_donor_correlation`). With ``k`` donors whose
        ``s2_d`` is positive it must exceed ``-1/(k - 1)``.
    regression_correlation : float, default 0.0
        Correlation of a transfer's error with the regression's.
    level : float, default 0.90
        Two-sided interval level.

    Returns
    -------
    BlendResult

    Raises
    ------
    ValueError
        Nothing to combine, a correlation outside ``[-1, 1]``, a donor
        correlation at or below ``-1/(k - 1)``, or a covariance matrix that is
        not positive definite.
    """
    _check_rho("donor", donor_correlation)
    _check_rho("regression", regression_correlation)
    donors = transfers["donor"].astype(str).tolist()
    if not donors and regression is None:
        raise ValueError("nothing to combine: no transfers and no regression")

    L = transfers["log_ar"].to_numpy(dtype=float) if donors else np.zeros(0)
    if donor_cov is not None:
        dc = donor_cov.rename(index=str, columns=str).loc[donors, donors].to_numpy(dtype=float)
    else:
        dc = (
            np.diag(transfers["donor_variance"].to_numpy(dtype=float))
            if donors
            else np.zeros((0, 0))
        )
    s = np.array(
        [
            model.variance(**{t: float(row[t]) for t in model.terms})
            for _, row in transfers.iterrows()
        ]
    )
    k = len(donors)
    _check_donor_bound(donor_correlation, int(np.count_nonzero(s > 0.0)))
    m = k + (regression is not None)
    C = np.zeros((m, m))
    C[:k, :k] = (
        dc + np.outer(L, L) * exponent_variance + donor_correlation * np.sqrt(np.outer(s, s))
    )
    C[np.arange(k), np.arange(k)] = np.diag(dc) + L**2 * exponent_variance + s
    est = list(transfers["transfer"].to_numpy(dtype=float)) if donors else []
    labels = list(donors)
    if regression is not None:
        r, v_r = regression
        C[k, k] = v_r
        C[:k, k] = C[k, :k] = regression_correlation * np.sqrt(np.diag(C)[:k] * v_r)
        est.append(float(r))
        labels.append(REGRESSION)

    try:
        np.linalg.cholesky(C)
    except np.linalg.LinAlgError:
        raise ValueError("the components' covariance is not positive definite") from None
    ci = np.linalg.inv(C)
    one = np.ones(m)
    denom = float(one @ ci @ one)
    w = ci @ one / denom
    log_q = float(w @ np.asarray(est))
    var = 1.0 / denom
    half = float(norm.ppf(0.5 + level / 2.0)) * math.sqrt(var)
    return BlendResult(
        log_q=log_q,
        variance=var,
        level=level,
        interval=(float(10.0 ** (log_q - half)), float(10.0 ** (log_q + half))),
        weights=pd.Series(w, index=labels),
        covariance=pd.DataFrame(C, index=labels, columns=labels),
        components=pd.DataFrame({"estimate": est, "variance": np.diag(C)}, index=labels),
    )


# --------------------------------------------------------------------------- band


@dataclass(frozen=True)
class Band:
    """The area-ratio band inside which a transfer beats the regression.

    Symmetric in log space, because ``s2_d`` is even in ``log_ar``.

    Attributes
    ----------
    half_width_log : float
        ``log10(AR*)``: 0 when empty, ``inf`` when unbounded.
    lower, upper : float
        ``1/AR*`` and ``AR*``.
    empty : bool
        The transfer's variance at ratio 1 already reaches ``V_R``.
    unbounded : bool
        Nothing in ``V_T`` grows with the ratio.
    """

    half_width_log: float
    lower: float
    upper: float
    empty: bool
    unbounded: bool


def _band_from(L: float) -> Band:
    if L == 0.0:
        return Band(0.0, 1.0, 1.0, empty=True, unbounded=False)
    if math.isinf(L):
        return Band(math.inf, 0.0, math.inf, empty=False, unbounded=True)
    return Band(L, 10.0**-L, 10.0**L, empty=False, unbounded=False)


def solve_band(
    model: DissimilarityModel,
    *,
    gage_variance: float,
    regression_variance: float,
    exponent_variance: float = 0.0,
    **features: float,
) -> Band:
    """The ratio ``AR*`` at which a transfer's variance equals the regression's.

    ``V_T(L) = V_g + L^2 (Var(b) + s_log_ar) + s0 + sum_other s_k x_k^2``, so
    ``L* = sqrt((V_R - V_g - s0 - other) / (Var(b) + s_log_ar))``.

    Parameters
    ----------
    model : DissimilarityModel
    gage_variance : float
        ``V_g``, the donor estimate's variance.
    regression_variance : float
        ``V_R`` at the target.
    exponent_variance : float, default 0.0
        ``Var(b)``.
    **features : float
        Values for every model term other than ``log_ar``.

    Returns
    -------
    Band

    Raises
    ------
    KeyError
        A term other than ``log_ar`` has no value.
    """
    others = {t: features[t] if t in features else None for t in model.terms if t != "log_ar"}
    missing = [t for t, v in others.items() if v is None]
    if missing:
        raise KeyError(f"no value for dissimilarity term(s) {missing}")
    fixed = model.variance(log_ar=0.0, **{t: float(v) for t, v in others.items() if v is not None})
    slope = exponent_variance + model.coefficients.get("log_ar", 0.0)
    gap = regression_variance - gage_variance - fixed
    if gap <= 0.0:
        return _band_from(0.0)
    if slope <= 0.0:
        return _band_from(math.inf)
    return _band_from(math.sqrt(gap / slope))


def _resample_pairs_by_gage(pairs: pd.DataFrame, counts: Mapping[str, int]) -> pd.DataFrame:
    """Keep each pair ``m_donor * m_receiver`` times."""
    reps = np.array(
        [
            counts.get(str(d), 0) * counts.get(str(r), 0)
            for d, r in zip(pairs["donor"], pairs["receiver"])
        ],
        dtype=int,
    )
    return pairs.loc[pairs.index.repeat(reps)].reset_index(drop=True)


@dataclass(frozen=True)
class BandBootstrap:
    """:func:`solve_band` on the full fit, with a gage-resampled interval.

    Attributes
    ----------
    estimate : Band
        From the full set of pairs.
    interval : (float, float)
        Percentile interval of ``AR*`` (1 for an empty band, ``inf`` for an
        unbounded one).
    level : float
    replicates : numpy.ndarray
        ``log10(AR*)`` per successful replicate.
    n_boot, n_failed : int
        Replicates drawn, and how many had too few pairs to fit.
    fraction_empty, fraction_unbounded : float
        Shares of successful replicates with an empty or unbounded band.
    """

    estimate: Band
    interval: Tuple[float, float]
    level: float
    replicates: np.ndarray
    n_boot: int
    n_failed: int
    fraction_empty: float
    fraction_unbounded: float


def bootstrap_band(
    pairs: pd.DataFrame,
    terms: Sequence[str],
    *,
    gage_variance: float,
    regression_variance: float,
    exponent_variance: float = 0.0,
    n_boot: int = 1000,
    level: float = 0.90,
    seed: Optional[int] = None,
    **features: float,
) -> BandBootstrap:
    """:func:`solve_band` with an interval from resampling gages.

    Each replicate draws as many gages as there are, with replacement, and
    keeps a pair ``m_i * m_j`` times, where ``m`` counts how often each end was
    drawn. Pairs share gages, so resampling pairs would understate the spread.

    Parameters
    ----------
    pairs : pd.DataFrame
        From :func:`build_pairs`.
    terms : sequence of str
        Passed to :func:`fit_dissimilarity`.
    gage_variance, regression_variance, exponent_variance, **features
        Passed to :func:`solve_band`.
    n_boot : int, default 1000
    level : float, default 0.90
    seed : int, optional

    Returns
    -------
    BandBootstrap
    """
    estimate = solve_band(
        fit_dissimilarity(pairs, terms),
        gage_variance=gage_variance,
        regression_variance=regression_variance,
        exponent_variance=exponent_variance,
        **features,
    )
    rng = np.random.default_rng(seed)
    gages = sorted(set(pairs["donor"].astype(str)) | set(pairs["receiver"].astype(str)))
    reps: List[float] = []
    failed = 0
    for _ in range(n_boot):
        drawn = Counter(rng.choice(gages, size=len(gages), replace=True).tolist())
        sample = _resample_pairs_by_gage(pairs, drawn)
        try:
            # A replicate's bound coefficient is not news; its band records it.
            model = _fit(sample, tuple(terms), warn=False)
        except ValueError:
            failed += 1
            continue
        band = solve_band(
            model,
            gage_variance=gage_variance,
            regression_variance=regression_variance,
            exponent_variance=exponent_variance,
            **features,
        )
        reps.append(band.half_width_log)
    arr = np.asarray(reps, dtype=float)
    if arr.size == 0:
        raise ValueError("every bootstrap replicate had too few pairs to fit")
    alpha = (1.0 - level) / 2.0
    lo, hi = np.quantile(arr, [alpha, 1.0 - alpha], method="inverted_cdf")
    return BandBootstrap(
        estimate=estimate,
        interval=(float(10.0**lo), float(10.0**hi)),
        level=level,
        replicates=arr,
        n_boot=n_boot,
        n_failed=failed,
        fraction_empty=float(np.mean(arr == 0.0)),
        fraction_unbounded=float(np.mean(np.isinf(arr))),
    )


# --------------------------------------------------------------------------- leave one out


def leave_one_out(
    gages: pd.DataFrame,
    connections: Sequence[Tuple[str, str]],
    exponents: RegressionExponents,
    aep: float,
    terms: Sequence[str] = ("log_ar",),
    *,
    exponent_variance: float = 0.0,
    sampling_cov: Optional[pd.DataFrame] = None,
    attributes: Sequence[str] = (),
    donor_correlation: Optional[float] = None,
    min_donor_correlation: Optional[float] = 0.0,
    regression_correlation: float = 0.0,
    level: float = 0.90,
) -> pd.DataFrame:
    """Score the weighting by holding out each connected gage in turn.

    For each gage ``j`` with at least one connection: ``s2_d`` (and, unless
    given, the donor correlation) is refitted from the pairs that do not
    involve ``j``; its neighbours' transfers to ``j`` are combined with its
    regression value; and the result is compared with its at-site estimate.

    The at-site value is itself an estimate, so coverage is checked against
    the variance of ``blend - y_j``: ``V_blend + V_j - 2 sum_k w_k C_kj``. The
    naive check against ``V_blend`` alone is reported as
    ``blend_covered_naive``.

    Parameters
    ----------
    gages : pd.DataFrame
        As :func:`build_pairs`, plus optional ``regression_log_q`` and
        ``regression_variance``. Without them the blend is transfers only.
    connections, exponents, aep, exponent_variance, sampling_cov, attributes
        As :func:`build_pairs`.
    terms : sequence of str
        Passed to :func:`fit_dissimilarity`.
    donor_correlation : float, optional
        Fixed donor correlation; ``None`` estimates it in each fold
        (:func:`estimate_donor_correlation`), and uses 0 where it cannot.
    min_donor_correlation : float or None, default 0.0
        Floor on a fold's estimated correlation. A negative estimate is noise
        (see :func:`estimate_donor_correlation`), and with several donors it
        can leave no valid covariance (#177). ``None`` uses the estimate as it
        is. A fixed ``donor_correlation`` is never floored.
    regression_correlation : float, default 0.0
        Passed to :func:`combine`.
    level : float, default 0.90

    Returns
    -------
    pd.DataFrame
        One row per connected gage: ``site_no``, ``n_donors``,
        ``donor_correlation`` (the value the fold used),
        ``donor_correlation_estimate`` (the fold's estimate before the floor,
        NaN when fixed or not estimable), ``donor_correlation_clipped``,
        ``observed``, ``blend``, ``blend_variance``, ``blend_error``,
        ``blend_covered``, ``blend_covered_naive``, ``transfer_only``,
        ``transfer_error``, ``transfer_covered``, ``regression``,
        ``regression_error``, ``regression_covered`` and ``note``. Errors are
        estimate minus observed, log10. A gage whose fold cannot be refitted or
        combined is not scored: its estimates are NaN, ``note`` says why
        (it is empty otherwise), and the run carries on.
    """
    b = _exponent(exponents, aep, "leave_one_out")
    net = _Network(gages, sampling_cov, attributes)
    conns = _check_connections(connections, net.sites)
    pairs = _pairs(net, conns, b, exponent_variance)
    tri = _triples(net, conns, b, exponent_variance)
    neighbours = _neighbours(conns)
    has_reg = {"regression_log_q", "regression_variance"} <= set(gages.columns)
    reg_q = _column(gages, net.sites, "regression_log_q") if has_reg else {}
    reg_v = _column(gages, net.sites, "regression_variance") if has_reg else {}
    zcrit = float(norm.ppf(0.5 + level / 2.0))
    donor = pairs["donor"].to_numpy()
    receiver = pairs["receiver"].to_numpy()

    unscored = dict.fromkeys(
        [
            "blend",
            "blend_variance",
            "blend_error",
            "blend_covered",
            "blend_covered_naive",
            "transfer_only",
            "transfer_error",
            "transfer_covered",
            "regression",
            "regression_error",
            "regression_covered",
        ],
        np.nan,
    )
    rows: List[Dict[str, object]] = []
    for j in net.sites:
        nbrs = neighbours.get(j, [])
        if not nbrs:
            continue
        row: Dict[str, object] = {
            "site_no": j,
            "n_donors": len(nbrs),
            "donor_correlation": np.nan,
            "donor_correlation_estimate": np.nan,
            "donor_correlation_clipped": False,
            "observed": net.log_q[j],
            **unscored,
            "note": "",
        }
        try:
            model = _fit(pairs.loc[(donor != j) & (receiver != j)], tuple(terms), warn=False)
        except ValueError as exc:
            row["note"] = f"s2_d refit failed: {exc}"
            rows.append(row)
            continue
        rho = donor_correlation
        if rho is None:
            outside = (tri["i"] != j) & (tri["k"] != j) & (tri["j"] != j)
            rho = _rho_from_triples(tri, model, outside)[0]
            if rho is not None:
                row["donor_correlation_estimate"] = rho
                if min_donor_correlation is not None and rho < min_donor_correlation:
                    rho = min_donor_correlation
                    row["donor_correlation_clipped"] = True
        rho_d = 0.0 if rho is None else rho
        row["donor_correlation"] = rho_d
        trows = []
        for i in nbrs:
            x = net.features(i, j)
            trows.append(
                {
                    "donor": i,
                    "transfer": net.log_q[i] + b * x["log_ar"],
                    "donor_variance": net.cov(i, i),
                    **x,
                }
            )
        transfers = pd.DataFrame(trows)
        dcov = net.cov_frame(nbrs)
        reg = (reg_q[j], reg_v[j]) if has_reg else None
        try:
            t_only = combine(
                model,
                transfers,
                None,
                exponent_variance=exponent_variance,
                donor_cov=dcov,
                donor_correlation=rho_d,
                level=level,
            )
            blend = (
                t_only
                if reg is None
                else combine(
                    model,
                    transfers,
                    reg,
                    exponent_variance=exponent_variance,
                    donor_cov=dcov,
                    donor_correlation=rho_d,
                    regression_correlation=regression_correlation,
                    level=level,
                )
            )
        except ValueError as exc:
            row["note"] = f"combine failed: {exc}"
            rows.append(row)
            continue
        y = net.log_q[j]
        v_j = net.cov(j, j)
        c_kj = np.array([net.cov(i, j) for i in nbrs])

        def scored(res: BlendResult) -> Tuple[float, bool, bool]:
            w = res.weights.reindex(nbrs).to_numpy(dtype=float)
            v_err = res.variance + v_j - 2.0 * float(w @ c_kj)
            err = res.log_q - y
            return (
                err,
                bool(abs(err) <= zcrit * math.sqrt(max(v_err, 0.0))),
                bool(abs(err) <= zcrit * math.sqrt(res.variance)),
            )

        b_err, b_cov, b_naive = scored(blend)
        t_err, t_cov, _ = scored(t_only)
        row.update(
            {
                "blend": blend.log_q,
                "blend_variance": blend.variance,
                "blend_error": b_err,
                "blend_covered": b_cov,
                "blend_covered_naive": b_naive,
                "transfer_only": t_only.log_q,
                "transfer_error": t_err,
                "transfer_covered": t_cov,
            }
        )
        if reg is not None:
            r_err = reg[0] - y
            row["regression"] = reg[0]
            row["regression_error"] = r_err
            row["regression_covered"] = bool(abs(r_err) <= zcrit * math.sqrt(reg[1] + v_j))
        rows.append(row)
    return pd.DataFrame(rows)
