"""
Transposing computed flows from a gaged donor basin to an ungaged target.

The drainage-area ratio method:

    Q_target(p) = Q_donor(p) * (A_target / A_donor) ** b(p)

where ``b(p)`` is the exponent on drainage area in the applicable published
USGS regional regression for exceedance probability ``p``. ``b = 1`` is the
naive area ratio and is almost never what the published regression says.

See ``docs/TRANSPOSITION_DESIGN.md`` for the full specification, including
why the exponent's provenance is mandatory rather than optional.

**The failure mode is silent.** A wrong exponent does not raise. It returns a
plausible discharge -- right order of magnitude, monotone across return
periods, passing every internal consistency check -- that is simply wrong by
10-30%. Two consequences run through this module:

- :class:`RegressionExponents` cannot be constructed without a citation. A
  transposition whose exponent has no stated source is not a result worth
  having.
- Every transposed quantile records whether its exponent was *published*,
  *interpolated* between published points, or *extrapolated* beyond them, so
  a report can say which is which instead of presenting all of them alike.

Published flood regressions typically give ``b`` at 50%, 20%, 10%, 4%, 2%,
1%, 0.5% and 0.2% AEP (Q2 through Q500). ``Bulletin17C.STANDARD_AEP`` runs
from 0.995 to 0.002, so the frequent end of a standard analysis sits well
below any flood regression's published range and is necessarily
extrapolated. Extrapolation defaults to ``"clamp"`` (hold the endpoint
exponent) for exactly that reason: linearly extending a trend 2.6 normal
deviates past the last observation is a curve fit, not an estimate. A caller
who has duration-regression exponents for the frequent end should supply
them in the same :class:`RegressionExponents`, which removes the
extrapolation entirely.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.special import ndtri

from flowfreq.core import FrequencyResults

logger = logging.getLogger(__name__)

#: Commonly cited USGS guidance for the drainage-area ratio method: the
#: target basin between half and one-and-a-half times the donor's area, on
#: the same stream, in the same hydrologic region. Outside this band the
#: published regression should be used directly instead of transposing.
DEFAULT_AREA_RATIO_RANGE: Tuple[float, float] = (0.5, 1.5)

#: Hard bounds applied to every exponent after interpolation or
#: extrapolation, whatever those produced. A drainage-area exponent outside
#: this range is not a number any published flood regression reports; getting
#: one means the interpolation walked somewhere it should not have.
DEFAULT_B_BOUNDS: Tuple[float, float] = (0.5, 1.0)

#: Tolerance for calling a query AEP "the same as" a published one, so that a
#: query at a published point is reported as published rather than
#: interpolated.
_AEP_MATCH_TOL: float = 1e-9

_SOURCE_PUBLISHED = "published"
_SOURCE_INTERPOLATED = "interpolated"
_SOURCE_EXTRAPOLATED = "extrapolated"
_SOURCE_CONSTANT = "constant"


def _normal_deviate(aep: np.ndarray) -> np.ndarray:
    """Standard normal deviate for an exceedance probability.

    Exponents are interpolated against this rather than against AEP
    directly: it is the abscissa the LP3 frequency curve itself is nearly
    linear in, published regression points are roughly evenly spaced along
    it, and it treats both tails symmetrically.

    Parameters
    ----------
    aep : np.ndarray
        Annual exceedance probabilities, strictly between 0 and 1.

    Returns
    -------
    np.ndarray
        ``z = ndtri(1 - aep)``; increasing as AEP decreases.
    """
    return ndtri(1.0 - np.asarray(aep, dtype=float))


@dataclass(frozen=True, eq=False)
class RegressionExponents:
    """Drainage-area exponents by AEP, from a published regression.

    Parameters
    ----------
    aeps : sequence of float
        Annual exceedance probabilities the exponents were published at,
        strictly between 0 and 1.
    exponents : sequence of float
        The exponent on drainage area at each of those AEPs.
    citation : str
        Where these came from, specifically enough that a reviewer can find
        it -- report number, table, region. **Required.** A transposition
        whose exponent has no stated source cannot be defended, so this class
        will not construct without one.
    region : str, optional
        The hydrologic region the exponents apply to, carried into
        provenance. This module cannot check that the donor and target are
        actually in it.
    valid_area_range_sqmi : tuple of float, optional
        The drainage-area range the regression was developed over. Used only
        to warn when a basin falls outside it.

    Raises
    ------
    ValueError
        If the citation is empty, the arrays disagree in length, any AEP is
        outside (0, 1), any exponent is non-finite or non-positive, or an
        AEP is repeated.
    """

    aeps: np.ndarray
    exponents: np.ndarray
    citation: str
    region: str = ""
    valid_area_range_sqmi: Optional[Tuple[float, float]] = None

    def __post_init__(self) -> None:
        if not self.citation or not self.citation.strip():
            raise ValueError(
                "RegressionExponents requires a citation naming where the exponents came "
                "from. A transposition whose exponent has no stated source is not a result "
                "worth having; see docs/TRANSPOSITION_DESIGN.md."
            )

        aeps = np.asarray(self.aeps, dtype=float).ravel()
        exponents = np.asarray(self.exponents, dtype=float).ravel()

        if aeps.size != exponents.size:
            raise ValueError(
                f"aeps and exponents must be the same length; got {aeps.size} and "
                f"{exponents.size}"
            )
        if aeps.size == 0:
            raise ValueError("aeps is empty; at least one published exponent is required")
        if np.any((aeps <= 0.0) | (aeps >= 1.0)):
            raise ValueError("every AEP must lie strictly between 0 and 1")
        if not np.all(np.isfinite(exponents)) or np.any(exponents <= 0.0):
            raise ValueError("every exponent must be finite and positive")

        # Sort by normal deviate ascending, i.e. AEP descending -- frequent
        # events first, so index 0 is always the frequent end and index -1 the
        # rare end wherever this is indexed below.
        order = np.argsort(_normal_deviate(aeps))
        aeps = aeps[order]
        exponents = exponents[order]

        if np.any(np.diff(aeps) == 0.0):
            raise ValueError("aeps contains a repeated value")

        object.__setattr__(self, "aeps", aeps)
        object.__setattr__(self, "exponents", exponents)

    @classmethod
    def constant(
        cls,
        exponent: float,
        citation: str,
        region: str = "",
        valid_area_range_sqmi: Optional[Tuple[float, float]] = None,
    ) -> "RegressionExponents":
        """A single exponent applied at every AEP.

        For the naive area ratio (``exponent=1.0``) or a sensitivity check.
        Still requires a citation -- "assumed 1.0 for sensitivity, not from a
        regression" is a perfectly good one, and saying so is the point.

        Parameters
        ----------
        exponent : float
            The exponent to use at every AEP.
        citation : str
            Where it came from, or why it was assumed.
        region : str, optional
            Carried into provenance.
        valid_area_range_sqmi : tuple of float, optional
            Carried into provenance.

        Returns
        -------
        RegressionExponents
        """
        return cls(
            aeps=np.array([0.5]),
            exponents=np.array([float(exponent)]),
            citation=citation,
            region=region,
            valid_area_range_sqmi=valid_area_range_sqmi,
        )

    @property
    def is_constant(self) -> bool:
        """Whether this is a single exponent applied at every AEP."""
        return bool(self.aeps.size == 1)

    def exponent_at(
        self,
        aeps: Sequence[float],
        *,
        interpolation: str = "linear",
        extrapolation: str = "clamp",
        b_bounds: Tuple[float, float] = DEFAULT_B_BOUNDS,
    ) -> pd.DataFrame:
        """Exponents at arbitrary AEPs, with each one's provenance.

        Parameters
        ----------
        aeps : sequence of float
            The AEPs to produce exponents for.
        interpolation : {"linear", "pchip"}, optional
            How to interpolate between published points, in normal-deviate
            space. ``"linear"`` is the default because a reviewer can
            reproduce it with a calculator; ``"pchip"`` is smoother and
            shape-preserving. Over a typical eight-point published set the
            two differ immaterially (see the tests, which measure it).
        extrapolation : {"clamp", "linear", "error"}, optional
            What to do beyond the published range. ``"clamp"`` (default)
            holds the endpoint exponent. ``"linear"`` extends the slope
            through the two nearest published points. ``"error"`` raises.
        b_bounds : tuple of float, optional
            Hard bounds applied after interpolation or extrapolation.

        Returns
        -------
        pd.DataFrame
            Columns ``aep``, ``exponent``, ``source`` (published /
            interpolated / extrapolated / constant) and ``clamped`` (whether
            ``b_bounds`` altered the value).

        Raises
        ------
        ValueError
            If an AEP is outside (0, 1), a mode is unrecognised, or
            ``extrapolation="error"`` and an AEP falls outside the published
            range.
        """
        query = np.asarray(aeps, dtype=float).ravel()
        if query.size == 0:
            raise ValueError("aeps is empty; nothing to compute an exponent for")
        if np.any((query <= 0.0) | (query >= 1.0)):
            raise ValueError("every AEP must lie strictly between 0 and 1")
        if interpolation not in ("linear", "pchip"):
            raise ValueError(f"unknown interpolation {interpolation!r}; use 'linear' or 'pchip'")
        if extrapolation not in ("clamp", "linear", "error"):
            raise ValueError(
                f"unknown extrapolation {extrapolation!r}; use 'clamp', 'linear' or 'error'"
            )

        if self.is_constant:
            values = np.full(query.shape, float(self.exponents[0]))
            sources = np.full(query.shape, _SOURCE_CONSTANT, dtype=object)
            return self._finish(query, values, sources, b_bounds)

        z_query = _normal_deviate(query)
        z_published = _normal_deviate(self.aeps)

        below = z_query < z_published[0]
        above = z_query > z_published[-1]
        outside = below | above

        if extrapolation == "error" and np.any(outside):
            raise ValueError(
                "extrapolation='error' and these AEPs fall outside the published range "
                f"[{self.aeps[0]:g}, {self.aeps[-1]:g}]: "
                f"{np.array2string(query[outside], precision=4)}"
            )

        values = self._interpolate(z_query, z_published, interpolation)

        if np.any(outside):
            values = self._extrapolate(values, z_query, z_published, below, above, extrapolation)

        sources = np.where(outside, _SOURCE_EXTRAPOLATED, _SOURCE_INTERPOLATED).astype(object)
        # A query landing on a published point is published, not interpolated.
        for i, aep_value in enumerate(query):
            if np.any(np.isclose(self.aeps, aep_value, rtol=0.0, atol=_AEP_MATCH_TOL)):
                sources[i] = _SOURCE_PUBLISHED

        return self._finish(query, values, sources, b_bounds)

    def _interpolate(
        self, z_query: np.ndarray, z_published: np.ndarray, interpolation: str
    ) -> np.ndarray:
        """Interpolated exponents; values outside the published range are
        placeholders that :meth:`_extrapolate` replaces."""
        if interpolation == "linear":
            # np.interp clamps outside the range on its own, which is the
            # right placeholder: it is already the "clamp" answer, and
            # _extrapolate overwrites it for the other modes.
            return np.interp(z_query, z_published, self.exponents)

        from scipy.interpolate import PchipInterpolator

        interpolator = PchipInterpolator(z_published, self.exponents, extrapolate=False)
        interpolated = np.asarray(interpolator(z_query), dtype=float)
        # PCHIP returns NaN outside; fill with the endpoint so the array is
        # finite before _extrapolate decides what those points deserve.
        interpolated = np.where(z_query < z_published[0], self.exponents[0], interpolated)
        interpolated = np.where(z_query > z_published[-1], self.exponents[-1], interpolated)
        return interpolated

    def _extrapolate(
        self,
        values: np.ndarray,
        z_query: np.ndarray,
        z_published: np.ndarray,
        below: np.ndarray,
        above: np.ndarray,
        extrapolation: str,
    ) -> np.ndarray:
        """Replace out-of-range values according to the extrapolation mode."""
        extrapolated = values.copy()

        if extrapolation == "clamp":
            extrapolated[below] = self.exponents[0]
            extrapolated[above] = self.exponents[-1]
            return extrapolated

        # extrapolation == "linear": extend the slope through the two nearest
        # published points at whichever end is being left.
        if np.any(below):
            slope_low = (self.exponents[1] - self.exponents[0]) / (z_published[1] - z_published[0])
            extrapolated[below] = self.exponents[0] + slope_low * (z_query[below] - z_published[0])
        if np.any(above):
            slope_high = (self.exponents[-1] - self.exponents[-2]) / (
                z_published[-1] - z_published[-2]
            )
            extrapolated[above] = self.exponents[-1] + slope_high * (
                z_query[above] - z_published[-1]
            )
        return extrapolated

    @staticmethod
    def _finish(
        query: np.ndarray,
        values: np.ndarray,
        sources: np.ndarray,
        b_bounds: Tuple[float, float],
    ) -> pd.DataFrame:
        """Apply the hard bounds and assemble the result frame."""
        low, high = float(b_bounds[0]), float(b_bounds[1])
        if low > high:
            raise ValueError(f"b_bounds is inverted: {b_bounds!r}")
        clamped = (values < low) | (values > high)
        bounded = np.clip(values, low, high)
        return pd.DataFrame(
            {
                "aep": query,
                "exponent": bounded,
                "source": sources,
                "clamped": clamped,
            }
        )


@dataclass
class TranspositionProvenance:
    """Everything needed to reproduce and defend a transposition.

    Attributes
    ----------
    donor_area_sqmi, target_area_sqmi : float
        The two drainage areas.
    area_ratio : float
        ``target_area_sqmi / donor_area_sqmi``.
    citation, region : str
        Carried through from the :class:`RegressionExponents` used.
    exponents : pd.DataFrame
        Per-AEP exponent and its source, as returned by
        :meth:`RegressionExponents.exponent_at`.
    interpolation, extrapolation : str
        The modes used.
    b_bounds, area_ratio_range : tuple of float
        The bounds applied.
    area_ratio_in_range : bool
        Whether ``area_ratio`` fell inside ``area_ratio_range``. False only
        when the caller passed ``allow_out_of_range=True``.
    donor_site_no, target_name : str, optional
        Labels for the report.
    """

    donor_area_sqmi: float
    target_area_sqmi: float
    area_ratio: float
    citation: str
    region: str
    exponents: pd.DataFrame
    interpolation: str
    extrapolation: str
    b_bounds: Tuple[float, float]
    area_ratio_range: Tuple[float, float]
    area_ratio_in_range: bool
    donor_site_no: Optional[str] = None
    target_name: Optional[str] = None

    @property
    def extrapolated_aeps(self) -> np.ndarray:
        """AEPs whose exponent was extrapolated beyond the published range."""
        mask = self.exponents["source"] == _SOURCE_EXTRAPOLATED
        return self.exponents.loc[mask, "aep"].to_numpy()

    @property
    def clamped_aeps(self) -> np.ndarray:
        """AEPs whose exponent was altered by the ``b_bounds`` clamp."""
        return self.exponents.loc[self.exponents["clamped"], "aep"].to_numpy()


@dataclass
class TransposedResults:
    """Flows transposed to a target basin, with their provenance.

    Deliberately **not** a :class:`~flowfreq.core.FrequencyResults`. That
    class describes a distribution fitted to a record; no distribution was
    fitted at the target, so there are no moments to report and reporting
    the donor's would be inventing a fit that was never performed. What
    exists here is a set of scaled quantiles and the arithmetic that
    produced them.

    Attributes
    ----------
    quantiles : pd.DataFrame
        Columns ``aep``, ``return_period``, ``donor_flow_cfs``, ``exponent``,
        ``area_ratio_factor`` (the multiplier actually applied) and
        ``flow_cfs`` (the transposed discharge).
    confidence_limits : pd.DataFrame
        The donor's limits scaled by the same factor. **First order only**:
        this treats the transposition itself as exact and therefore
        understates the true uncertainty, which also carries the regression's
        own error. Empty if the donor had no limits.
    provenance : TranspositionProvenance
    """

    quantiles: pd.DataFrame
    confidence_limits: pd.DataFrame
    provenance: TranspositionProvenance

    def to_markdown(self) -> str:
        """A report-ready account of the transposition, arithmetic included.

        Returns
        -------
        str
            Markdown: the inputs, the per-quantile arithmetic, and every
            caveat that applies to this particular transposition.
        """
        p = self.provenance
        target = p.target_name or "ungaged target"
        donor = p.donor_site_no or "donor gage"

        lines = [
            f"# Transposed flows: {target}",
            "",
            "Drainage-area ratio transposition, "
            "`Q_target(p) = Q_donor(p) * (A_target/A_donor) ** b(p)`.",
            "",
            f"- Donor: {donor}, {p.donor_area_sqmi:,.1f} sq mi",
            f"- Target: {target}, {p.target_area_sqmi:,.1f} sq mi",
            f"- Area ratio: {p.area_ratio:.4f}",
            f"- Exponent source: {p.citation}",
        ]
        if p.region:
            lines.append(f"- Region: {p.region}")
        lines.extend(
            [
                f"- Interpolation: {p.interpolation} (in normal-deviate space); "
                f"extrapolation: {p.extrapolation}",
                "",
                "## Quantiles",
                "",
                "| AEP | Return period (yr) | Donor Q (cfs) | b | Ratio^b | Target Q (cfs) | Exponent source |",
                "|---:|---:|---:|---:|---:|---:|---|",
            ]
        )

        for _, row in self.quantiles.iterrows():
            lines.append(
                f"| {row['aep']:.4g} | {row['return_period']:.4g} | "
                f"{row['donor_flow_cfs']:,.0f} | {row['exponent']:.4f} | "
                f"{row['area_ratio_factor']:.4f} | {row['flow_cfs']:,.0f} | "
                f"{row['source']} |"
            )

        lines.extend(["", "## Caveats", ""])

        extrapolated = p.extrapolated_aeps
        if extrapolated.size:
            lines.append(
                f"- **{extrapolated.size} of {len(self.quantiles)} quantiles use an "
                "extrapolated exponent**, at AEP "
                + ", ".join(f"{a:g}" for a in extrapolated)
                + ". These fall outside the range the cited regression publishes, so "
                "their exponent is an assumption of this analysis rather than a "
                "published value."
            )
        else:
            lines.append(
                "- Every exponent came from within the cited regression's published range."
            )

        clamped = p.clamped_aeps
        if clamped.size:
            lines.append(
                f"- Exponents at AEP "
                + ", ".join(f"{a:g}" for a in clamped)
                + f" were clamped to the bounds {p.b_bounds}."
            )

        if not p.area_ratio_in_range:
            lines.append(
                f"- **The area ratio {p.area_ratio:.4f} is outside the applicable range "
                f"{p.area_ratio_range}.** This transposition was run with an explicit "
                "override; the published regression should be preferred at the target."
            )

        if not self.confidence_limits.empty:
            lines.append(
                "- Confidence limits are the donor's, scaled by the same factor. That is "
                "first order and understates the true uncertainty: it treats the "
                "transposition as exact and carries none of the regression's own error."
            )

        return "\n".join(lines) + "\n"


def transpose_frequency(
    results: FrequencyResults,
    donor_area_sqmi: float,
    target_area_sqmi: float,
    exponents: RegressionExponents,
    *,
    interpolation: str = "linear",
    extrapolation: str = "clamp",
    b_bounds: Tuple[float, float] = DEFAULT_B_BOUNDS,
    area_ratio_range: Tuple[float, float] = DEFAULT_AREA_RATIO_RANGE,
    allow_out_of_range: bool = False,
    max_aep: Optional[float] = None,
    donor_site_no: Optional[str] = None,
    target_name: Optional[str] = None,
) -> TransposedResults:
    """Transpose a donor gage's fitted flood quantiles to an ungaged target.

    Every quantile present in ``results.quantiles`` is transposed, including
    the frequent end (``Bulletin17C.STANDARD_AEP`` runs to 0.995). Those
    quantiles sit below the range any flood regression publishes, so their
    exponents are extrapolated and flagged as such in the result; pass
    ``max_aep`` to leave them out instead, or supply exponents covering the
    frequent end to remove the extrapolation.

    Parameters
    ----------
    results : FrequencyResults
        The donor's analysis. Its ``quantiles`` frame supplies the AEPs.
    donor_area_sqmi, target_area_sqmi : float
        Drainage areas, same units. Both must be positive.
    exponents : RegressionExponents
        Drainage-area exponents and their citation.
    interpolation : {"linear", "pchip"}, optional
        Passed to :meth:`RegressionExponents.exponent_at`.
    extrapolation : {"clamp", "linear", "error"}, optional
        Passed to :meth:`RegressionExponents.exponent_at`. Defaults to
        ``"clamp"``.
    b_bounds : tuple of float, optional
        Hard bounds on the exponent.
    area_ratio_range : tuple of float, optional
        The applicable range for the drainage-area ratio method.
    allow_out_of_range : bool, optional
        Transpose anyway when the area ratio falls outside
        ``area_ratio_range``, recording the violation in provenance. False by
        default: outside that band the published regression should be used
        directly, and a plausible wrong number is worse than a refusal.
    max_aep : float, optional
        Drop quantiles more frequent than this before transposing. ``None``
        (default) transposes everything present.
    donor_site_no, target_name : str, optional
        Labels for the report.

    Returns
    -------
    TransposedResults

    Raises
    ------
    ValueError
        If an area is not positive, the donor has no quantiles, ``max_aep``
        leaves nothing to transpose, or the area ratio is outside
        ``area_ratio_range`` without ``allow_out_of_range``.
    """
    if not np.isfinite(donor_area_sqmi) or donor_area_sqmi <= 0:
        raise ValueError(f"donor_area_sqmi must be positive; got {donor_area_sqmi!r}")
    if not np.isfinite(target_area_sqmi) or target_area_sqmi <= 0:
        raise ValueError(f"target_area_sqmi must be positive; got {target_area_sqmi!r}")

    donor_quantiles = results.quantiles
    if donor_quantiles is None or donor_quantiles.empty:
        raise ValueError(
            "the donor FrequencyResults carries no quantiles; run the analysis with "
            "compute_quantiles before transposing"
        )

    selected = donor_quantiles
    if max_aep is not None:
        selected = selected.loc[selected["aep"] <= max_aep + _AEP_MATCH_TOL]
        if selected.empty:
            raise ValueError(
                f"max_aep={max_aep!r} excluded every quantile the donor carries "
                f"(AEPs {donor_quantiles['aep'].min():g} to "
                f"{donor_quantiles['aep'].max():g})"
            )

    area_ratio = float(target_area_sqmi) / float(donor_area_sqmi)
    ratio_low, ratio_high = area_ratio_range
    in_range = bool(ratio_low <= area_ratio <= ratio_high)
    if not in_range:
        if not allow_out_of_range:
            raise ValueError(
                f"area ratio {area_ratio:.4f} is outside the applicable range "
                f"({ratio_low}, {ratio_high}) for the drainage-area ratio method. "
                "Use the published regression at the target directly, or pass "
                "allow_out_of_range=True to transpose anyway and record the violation."
            )
        logger.warning(
            "Transposing at an area ratio of %.4f, outside the applicable range %s. "
            "The published regression should be preferred at the target.",
            area_ratio,
            area_ratio_range,
        )

    valid_area = exponents.valid_area_range_sqmi
    if valid_area is not None:
        for label, area in (("donor", donor_area_sqmi), ("target", target_area_sqmi)):
            if not valid_area[0] <= area <= valid_area[1]:
                logger.warning(
                    "%s drainage area %.1f sq mi is outside the range %s the cited "
                    "regression was developed over.",
                    label,
                    area,
                    valid_area,
                )

    aeps = selected["aep"].to_numpy(dtype=float)
    exponent_frame = exponents.exponent_at(
        aeps,
        interpolation=interpolation,
        extrapolation=extrapolation,
        b_bounds=b_bounds,
    )
    b_values = exponent_frame["exponent"].to_numpy(dtype=float)
    factors = area_ratio**b_values

    extrapolated = exponent_frame.loc[
        exponent_frame["source"] == _SOURCE_EXTRAPOLATED, "aep"
    ].to_numpy()
    if extrapolated.size:
        logger.warning(
            "%d of %d quantiles use an exponent extrapolated beyond the cited "
            "regression's published range (AEP %s); their exponent is an assumption "
            "of this analysis rather than a published value.",
            extrapolated.size,
            len(aeps),
            np.array2string(extrapolated, precision=4),
        )

    donor_flows = selected["flow_cfs"].to_numpy(dtype=float)
    transposed_quantiles = pd.DataFrame(
        {
            "aep": aeps,
            "return_period": selected["return_period"].to_numpy(dtype=float),
            "donor_flow_cfs": donor_flows,
            "exponent": b_values,
            "area_ratio_factor": factors,
            "flow_cfs": donor_flows * factors,
            "source": exponent_frame["source"].to_numpy(),
        }
    )

    donor_limits = results.confidence_limits
    if donor_limits is None or donor_limits.empty:
        transposed_limits = pd.DataFrame()
    else:
        limits = donor_limits.loc[donor_limits["aep"].isin(aeps)].copy()
        # Match each retained limit row to the factor computed for its AEP,
        # rather than assuming the two frames are already aligned.
        factor_by_aep = dict(zip(aeps, factors))
        limit_factors = limits["aep"].map(factor_by_aep).to_numpy(dtype=float)
        transposed_limits = pd.DataFrame(
            {
                "aep": limits["aep"].to_numpy(dtype=float),
                "return_period": limits["return_period"].to_numpy(dtype=float),
                "flow_cfs": limits["flow_cfs"].to_numpy(dtype=float) * limit_factors,
                "lower_5pct": limits["lower_5pct"].to_numpy(dtype=float) * limit_factors,
                "upper_5pct": limits["upper_5pct"].to_numpy(dtype=float) * limit_factors,
            }
        )

    provenance = TranspositionProvenance(
        donor_area_sqmi=float(donor_area_sqmi),
        target_area_sqmi=float(target_area_sqmi),
        area_ratio=area_ratio,
        citation=exponents.citation,
        region=exponents.region,
        exponents=exponent_frame,
        interpolation=interpolation,
        extrapolation=extrapolation,
        b_bounds=b_bounds,
        area_ratio_range=area_ratio_range,
        area_ratio_in_range=in_range,
        donor_site_no=donor_site_no,
        target_name=target_name,
    )

    return TransposedResults(
        quantiles=transposed_quantiles,
        confidence_limits=transposed_limits,
        provenance=provenance,
    )
