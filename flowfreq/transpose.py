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
from typing import Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from scipy.special import ndtri

from flowfreq.core import FrequencyResults, LowFlowResults

logger = logging.getLogger(__name__)

#: Commonly cited USGS guidance for the drainage-area ratio method: the
#: target basin between half and one-and-a-half times the donor's area, on
#: the same stream, in the same hydrologic region. Outside this band the
#: published regression should be used directly instead of transposing.
DEFAULT_AREA_RATIO_RANGE: Tuple[float, float] = (0.5, 1.5)

#: Tighter band for low flows. A low-flow statistic is controlled by
#: baseflow storage rather than area, so the range over which a pure
#: area-ratio transfer is defensible is narrower than it is for floods.
LOW_FLOW_AREA_RATIO_RANGE: Tuple[float, float] = (0.7, 1.3)

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

#: What kind of probability an exponent set is indexed by. These are not
#: interchangeable, and mixing them up is the category error
#: ``docs/TRANSPOSITION_DESIGN.md`` §7 warns about -- a flood regression's
#: exponent describes flood response and says nothing about 7Q10, which is
#: controlled by baseflow storage and geology. Each transpose function
#: requires its own kind, so the mistake raises instead of returning a
#: plausible wrong number.
#:
#: - ``"aep"`` -- annual exceedance probability of an annual maximum (floods)
#: - ``"exceedance"`` -- fraction of time exceeded (flow-duration curves)
#: - ``"non_exceedance"`` -- annual non-exceedance probability of an annual
#:   minimum (low flows: 7Q10 is non-exceedance 0.1)
PROBABILITY_KINDS: Tuple[str, ...] = ("aep", "exceedance", "non_exceedance")


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
    probability_kind : {"aep", "exceedance", "non_exceedance"}, optional
        What ``aeps`` actually holds. Defaults to ``"aep"`` (floods). See
        :data:`PROBABILITY_KINDS`: each transpose function requires the kind
        that matches the statistic it transposes, so handing a flood
        exponent set to :func:`transpose_low_flow` raises rather than
        returning a plausible wrong number.

        The field is named ``aeps`` for all three because the rest of this
        library indexes probabilities that way and a second name would be
        worse; for duration and low-flow sets it holds exceedance fractions
        and non-exceedance probabilities respectively, both of which occupy
        the same (0, 1) domain.

    Raises
    ------
    ValueError
        If the citation is empty, the arrays disagree in length, any
        probability is outside (0, 1), any exponent is non-finite or
        non-positive, an entry is repeated, or ``probability_kind`` is not
        one of :data:`PROBABILITY_KINDS`.
    """

    aeps: np.ndarray
    exponents: np.ndarray
    citation: str
    region: str = ""
    valid_area_range_sqmi: Optional[Tuple[float, float]] = None
    probability_kind: str = "aep"

    def __post_init__(self) -> None:
        if not self.citation or not self.citation.strip():
            raise ValueError(
                "RegressionExponents requires a citation naming where the exponents came "
                "from. A transposition whose exponent has no stated source is not a result "
                "worth having; see docs/TRANSPOSITION_DESIGN.md."
            )

        if self.probability_kind not in PROBABILITY_KINDS:
            raise ValueError(
                f"unknown probability_kind {self.probability_kind!r}; use one of "
                f"{PROBABILITY_KINDS}"
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
        probability_kind: str = "aep",
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
            probability_kind=probability_kind,
        )

    @property
    def is_constant(self) -> bool:
        """Whether this is a single exponent applied at every AEP."""
        return bool(self.aeps.size == 1)

    def exponent_at(
        self,
        aeps: Union[Sequence[float], np.ndarray],
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
    statistic: str = "flood frequency"
    degenerate_count: int = 0
    monotonic: bool = True
    hydrogeologic_setting: Optional[str] = None
    bfi_difference: Optional[float] = None

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
        if p.hydrogeologic_setting:
            lines.append(f"- Shared setting asserted by the analyst: {p.hydrogeologic_setting}")
        lines.extend(
            [
                f"- Statistic: {p.statistic}",
                f"- Interpolation: {p.interpolation} (in normal-deviate space); "
                f"extrapolation: {p.extrapolation}",
                "",
                "## Quantiles",
                "",
            ]
        )

        # Which probability the rows are indexed by depends on the statistic:
        # AEP for floods, exceedance fraction for duration, non-exceedance for
        # low flow. Pick whichever the frame actually carries rather than
        # assuming the flood one.
        prob_column, prob_label = _probability_column(self.quantiles)
        has_return_period = "return_period" in self.quantiles.columns

        header = f"| {prob_label} |"
        rule = "|---:|"
        if has_return_period:
            header += " Return period (yr) |"
            rule += "---:|"
        header += " Donor Q (cfs) | b | Ratio^b | Target Q (cfs) | Exponent source |"
        rule += "---:|---:|---:|---:|---|"
        lines.extend([header, rule])

        for _, row in self.quantiles.iterrows():
            cells = f"| {row[prob_column]:.4g} |"
            if has_return_period:
                cells += f" {row['return_period']:.4g} |"
            target_flow = row["flow_cfs"]
            target_text = "not transposable" if pd.isna(target_flow) else f"{target_flow:,.0f}"
            cells += (
                f" {row['donor_flow_cfs']:,.0f} | {row['exponent']:.4f} | "
                f"{row['area_ratio_factor']:.4f} | {target_text} | {row['source']} |"
            )
            lines.append(cells)

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

        if p.degenerate_count:
            lines.append(
                f"- **{p.degenerate_count} donor statistic(s) were zero** and could not be "
                "transposed: scaling a zero by an area ratio returns zero, which would "
                "assert the target is dry rather than estimate it. Those rows are blank, "
                "not zero."
            )

        if p.hydrogeologic_setting:
            lines.append(
                "- Low-flow transposition rests on the analyst's assertion that donor and "
                f"target share a setting ({p.hydrogeologic_setting}). Low flows are "
                "controlled by baseflow storage and geology rather than drainage area; "
                "where the target's own basin characteristics are available, a published "
                "low-flow regression is the better tool than transposition."
            )

        if p.bfi_difference is not None:
            lines.append(f"- Donor and target baseflow indices differ by {p.bfi_difference:.3f}.")

        if not self.confidence_limits.empty:
            lines.append(
                "- Confidence limits are the donor's, scaled by the same factor. That is "
                "first order and understates the true uncertainty: it treats the "
                "transposition as exact and carries none of the regression's own error."
            )

        return "\n".join(lines) + "\n"


def _probability_column(frame: pd.DataFrame) -> Tuple[str, str]:
    """The probability column a transposed frame is indexed by, and its label."""
    for column, label in (
        ("aep", "AEP"),
        ("exceedance_prob", "Exceedance"),
        ("non_exceedance_prob", "Non-exceedance"),
    ):
        if column in frame.columns:
            return column, label
    raise KeyError(f"no probability column in {list(frame.columns)}")


def _require_kind(exponents: RegressionExponents, expected: str, caller: str) -> None:
    """Refuse an exponent set indexed by the wrong kind of probability.

    Structural, not advisory: handing a flood regression's exponents to
    :func:`transpose_low_flow` is a category error rather than an
    approximation (``docs/TRANSPOSITION_DESIGN.md`` §7), and it is invisible
    in the output if it is allowed through.
    """
    if exponents.probability_kind != expected:
        raise ValueError(
            f"{caller} needs exponents indexed by {expected!r}, but this set is indexed "
            f"by {exponents.probability_kind!r}. These are not interchangeable: a flood "
            "regression's drainage-area exponent describes flood response and says "
            "nothing about a low-flow or duration statistic, which are controlled by "
            "baseflow storage and geology. Use the regression published for the "
            "statistic being transposed."
        )


def _check_areas(
    donor_area_sqmi: float,
    target_area_sqmi: float,
    exponents: RegressionExponents,
    area_ratio_range: Tuple[float, float],
    allow_out_of_range: bool,
) -> Tuple[float, bool]:
    """Validate the two areas and the ratio between them.

    Shared by every transpose function so the guardrail is one
    implementation rather than three that drift apart.

    Returns
    -------
    tuple of (float, bool)
        The area ratio, and whether it fell inside ``area_ratio_range``.
    """
    if not np.isfinite(donor_area_sqmi) or donor_area_sqmi <= 0:
        raise ValueError(f"donor_area_sqmi must be positive; got {donor_area_sqmi!r}")
    if not np.isfinite(target_area_sqmi) or target_area_sqmi <= 0:
        raise ValueError(f"target_area_sqmi must be positive; got {target_area_sqmi!r}")

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

    return area_ratio, in_range


def _check_monotonic(
    probs: np.ndarray,
    flows: np.ndarray,
    flow_rises_with_prob: bool,
    on_non_monotonic: str,
    context: str,
) -> bool:
    """Verify the transposed curve is still physically ordered.

    Transposing with an exponent that varies across the curve can invert it.
    Writing ``Q_t(p) = Q_d(p) * r**b(p)`` in logs,

        d/dp [ln Q_t] = d/dp [ln Q_d] + ln(r) * db/dp

    so when ``r < 1`` (``ln r`` negative) and ``b`` falls toward the dry end
    (``db/dp`` negative), the second term is *positive* and fights the first.
    If the donor's own curve is flat enough at that end, the second term
    wins and the transposed curve turns back on itself -- a flow exceeded
    99% of the time coming out larger than the flow exceeded 95% of the
    time, which cannot happen.

    That combination is not exotic. A flat dry end is exactly the signature
    of a groundwater-dominated basin, and a ``b`` that falls toward the dry
    end is exactly what the duration-regression literature reports, because
    low flows scale with storage rather than area. So this is checked rather
    than assumed, and it is checked on the transposed result rather than
    argued about in advance.

    Returns
    -------
    bool
        Whether the curve came out monotone.

    Raises
    ------
    ValueError
        If it did not and ``on_non_monotonic="raise"``.
    """
    if on_non_monotonic not in ("raise", "warn"):
        raise ValueError(f"unknown on_non_monotonic {on_non_monotonic!r}; use 'raise' or 'warn'")

    finite = np.isfinite(flows)
    if finite.sum() < 2:
        return True

    order = np.argsort(probs[finite])
    ordered = flows[finite][order]
    steps = np.diff(ordered)
    monotonic = bool(np.all(steps > 0) if flow_rises_with_prob else np.all(steps < 0))

    if not monotonic:
        direction = "rise" if flow_rises_with_prob else "fall"
        message = (
            f"the transposed {context} is not monotone: flows do not {direction} "
            "consistently with probability, which is physically impossible. The "
            "exponent set varies too steeply across the curve relative to how flat "
            "the donor's own curve is at that end -- see _check_monotonic. Use a "
            "smoother exponent set (or one exponent across the flat tail) rather "
            "than accepting the result; clipping the flows would invent numbers "
            "with no basis."
        )
        if on_non_monotonic == "raise":
            raise ValueError(message)
        logger.warning(message)

    return monotonic


def _warn_if_extrapolated(exponent_frame: pd.DataFrame, n_total: int) -> np.ndarray:
    """Log once naming the probabilities whose exponent was extrapolated."""
    extrapolated = exponent_frame.loc[
        exponent_frame["source"] == _SOURCE_EXTRAPOLATED, "aep"
    ].to_numpy()
    if extrapolated.size:
        logger.warning(
            "%d of %d quantiles use an exponent extrapolated beyond the cited "
            "regression's published range (probability %s); their exponent is an "
            "assumption of this analysis rather than a published value.",
            extrapolated.size,
            n_total,
            np.array2string(extrapolated, precision=4),
        )
    return extrapolated


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
    on_non_monotonic: str = "raise",
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
    _require_kind(exponents, "aep", "transpose_frequency")
    area_ratio, in_range = _check_areas(
        donor_area_sqmi,
        target_area_sqmi,
        exponents,
        area_ratio_range,
        allow_out_of_range,
    )

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

    aeps = selected["aep"].to_numpy(dtype=float)
    exponent_frame = exponents.exponent_at(
        aeps,
        interpolation=interpolation,
        extrapolation=extrapolation,
        b_bounds=b_bounds,
    )
    b_values = exponent_frame["exponent"].to_numpy(dtype=float)
    factors = area_ratio**b_values

    _warn_if_extrapolated(exponent_frame, len(aeps))

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

    monotonic = _check_monotonic(
        aeps,
        transposed_quantiles["flow_cfs"].to_numpy(dtype=float),
        flow_rises_with_prob=False,
        on_non_monotonic=on_non_monotonic,
        context="frequency curve",
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
        monotonic=monotonic,
    )

    return TransposedResults(
        quantiles=transposed_quantiles,
        confidence_limits=transposed_limits,
        provenance=provenance,
    )


def transpose_duration(
    curve: pd.DataFrame,
    donor_area_sqmi: float,
    target_area_sqmi: float,
    exponents: RegressionExponents,
    *,
    interpolation: str = "linear",
    extrapolation: str = "clamp",
    b_bounds: Tuple[float, float] = DEFAULT_B_BOUNDS,
    area_ratio_range: Tuple[float, float] = DEFAULT_AREA_RATIO_RANGE,
    allow_out_of_range: bool = False,
    on_non_monotonic: str = "raise",
    donor_site_no: Optional[str] = None,
    target_name: Optional[str] = None,
) -> TransposedResults:
    """Transpose flow-duration statistics to an ungaged target.

    Same drainage-area ratio machinery as :func:`transpose_frequency`, with
    exponents indexed by *exceedance fraction* rather than AEP, because ``b``
    is not constant across a duration curve: it is near-linear at the wet
    end, where contributing area largely sets the flow, and falls away at the
    dry end, where geology and storage dominate. A single exponent across the
    whole curve is wrong at one end or the other, which is why this takes a
    whole exponent set rather than one number.

    Parameters
    ----------
    curve : pd.DataFrame
        The donor's duration statistics as
        :func:`flowfreq.regime.flow_duration_curve` returns them: a
        ``flow_cfs`` column plus ``exceedance_prob`` or ``exceedance_pct``.
    donor_area_sqmi, target_area_sqmi : float
        Drainage areas, same units.
    exponents : RegressionExponents
        Must have ``probability_kind="exceedance"``.
    interpolation, extrapolation, b_bounds : optional
        As :func:`transpose_frequency`.
    area_ratio_range, allow_out_of_range : optional
        As :func:`transpose_frequency`.
    donor_site_no, target_name : str, optional
        Labels for the report.

    Returns
    -------
    TransposedResults
        ``confidence_limits`` is always empty: a duration curve is an
        empirical description of a record, not a fitted distribution, so
        there are no limits to carry.

    Raises
    ------
    ValueError
        If the exponent set is indexed by the wrong kind of probability, the
        curve is empty, or the areas fail their checks.
    KeyError
        If the curve lacks a flow or probability column.

    Notes
    -----
    A donor duration statistic of zero cannot be transposed: ``0 * anything``
    is ``0``, which asserts something about the target that the donor's
    record does not support. Those rows come back as ``NaN`` and are counted
    in the report rather than quietly reported as zero flow.
    """
    _require_kind(exponents, "exceedance", "transpose_duration")
    area_ratio, in_range = _check_areas(
        donor_area_sqmi,
        target_area_sqmi,
        exponents,
        area_ratio_range,
        allow_out_of_range,
    )

    if curve is None or curve.empty:
        raise ValueError("curve is empty; nothing to transpose")
    if "flow_cfs" not in curve.columns:
        raise KeyError("curve has no 'flow_cfs' column")

    if "exceedance_prob" in curve.columns:
        probs = curve["exceedance_prob"].to_numpy(dtype=float)
    elif "exceedance_pct" in curve.columns:
        probs = curve["exceedance_pct"].to_numpy(dtype=float) / 100.0
    else:
        raise KeyError("curve has neither 'exceedance_prob' nor 'exceedance_pct'")

    exponent_frame = exponents.exponent_at(
        probs,
        interpolation=interpolation,
        extrapolation=extrapolation,
        b_bounds=b_bounds,
    )
    _warn_if_extrapolated(exponent_frame, len(probs))

    b_values = exponent_frame["exponent"].to_numpy(dtype=float)
    factors = area_ratio**b_values
    donor_flows = curve["flow_cfs"].to_numpy(dtype=float)

    # A zero donor statistic is degenerate under a ratio method: scaling it
    # returns zero, which is an assertion about the target rather than an
    # estimate of it. NaN says "this could not be transposed", which is true.
    degenerate = ~(donor_flows > 0.0)
    transposed = np.where(degenerate, np.nan, donor_flows * factors)
    n_degenerate = int(np.count_nonzero(degenerate))
    if n_degenerate:
        logger.warning(
            "%d duration statistic(s) are zero or non-positive at the donor and cannot "
            "be transposed by a ratio method; reported as NaN, not as zero flow.",
            n_degenerate,
        )

    quantiles = pd.DataFrame(
        {
            "exceedance_prob": probs,
            "exceedance_pct": probs * 100.0,
            "donor_flow_cfs": donor_flows,
            "exponent": b_values,
            "area_ratio_factor": factors,
            "flow_cfs": transposed,
            "source": exponent_frame["source"].to_numpy(),
        }
    )

    monotonic = _check_monotonic(
        probs,
        transposed,
        flow_rises_with_prob=False,
        on_non_monotonic=on_non_monotonic,
        context="duration curve",
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
        statistic="flow duration",
        degenerate_count=n_degenerate,
        monotonic=monotonic,
    )

    return TransposedResults(
        quantiles=quantiles,
        confidence_limits=pd.DataFrame(),
        provenance=provenance,
    )


def transpose_low_flow(
    results: LowFlowResults,
    donor_area_sqmi: float,
    target_area_sqmi: float,
    exponents: RegressionExponents,
    *,
    hydrogeologic_setting: str,
    interpolation: str = "linear",
    extrapolation: str = "clamp",
    b_bounds: Tuple[float, float] = DEFAULT_B_BOUNDS,
    area_ratio_range: Tuple[float, float] = LOW_FLOW_AREA_RATIO_RANGE,
    allow_out_of_range: bool = False,
    on_non_monotonic: str = "raise",
    donor_bfi: Optional[float] = None,
    target_bfi: Optional[float] = None,
    max_bfi_difference: float = 0.15,
    max_p_zero: float = 0.0,
    donor_site_no: Optional[str] = None,
    target_name: Optional[str] = None,
) -> TransposedResults:
    """Transpose low-flow statistics to an ungaged target, under protest.

    Deliberately a separate function from :func:`transpose_frequency` rather
    than a flag on it. Low flows are controlled by baseflow storage --
    surficial geology, aquifer transmissivity, soil permeability -- not by
    drainage area. Two adjacent basins of identical area can differ
    severalfold in 7Q10, and one of them can be zero. That is why published
    low-flow regressions carry a geology or baseflow term and flood
    regressions do not, and why the guardrails here are stricter:

    - The exponent set must be indexed by ``"non_exceedance"``. A flood
      exponent is not an approximation here, it is the wrong quantity.
    - The area-ratio band defaults to :data:`LOW_FLOW_AREA_RATIO_RANGE`
      (0.7-1.3), tighter than the flood band.
    - ``hydrogeologic_setting`` is required. This module cannot verify that
      donor and target share an aquifer or physiographic province, so it
      makes the caller assert it in writing and records the assertion.
    - A donor statistic of zero, or a donor that goes dry at all
      (``p_zero > max_p_zero``), refuses rather than scaling a zero.

    **Where the target's own basin characteristics are available, prefer the
    published low-flow regression directly.** Transposition is the fallback
    here, not the better tool.

    Parameters
    ----------
    results : LowFlowResults
        The donor's low-flow analysis.
    donor_area_sqmi, target_area_sqmi : float
        Drainage areas, same units.
    exponents : RegressionExponents
        Must have ``probability_kind="non_exceedance"``.
    hydrogeologic_setting : str
        The setting the caller asserts donor and target share, e.g.
        "Valley and Ridge carbonate, same aquifer". Required and recorded;
        an empty string raises.
    donor_bfi, target_bfi : float, optional
        Baseflow indices from :func:`flowfreq.regime.baseflow_index`. When
        both are given and differ by more than ``max_bfi_difference``, a
        warning is logged and the difference recorded. Note ``regime``'s own
        caveat: two BFIs are comparable only when computed by the same
        method with the same parameters, so compute both sides yourself.
    max_bfi_difference : float, optional
        How far the two may differ before it is worth saying so.
    max_p_zero : float, optional
        The largest fraction of zero-flow years the donor may have. Defaults
        to 0.0: a donor that ever goes dry cannot tell you what the target
        does.
    interpolation, extrapolation, b_bounds : optional
        As :func:`transpose_frequency`.
    area_ratio_range, allow_out_of_range : optional
        As :func:`transpose_frequency`, with a tighter default band.
    donor_site_no, target_name : str, optional
        Labels for the report.

    Returns
    -------
    TransposedResults
        ``quantiles`` is indexed by ``non_exceedance_prob``.

    Raises
    ------
    ValueError
        If the exponent kind is wrong, ``hydrogeologic_setting`` is empty,
        the donor carries no quantiles, any donor statistic is zero, or the
        donor's ``p_zero`` exceeds ``max_p_zero``.
    """
    if not hydrogeologic_setting or not hydrogeologic_setting.strip():
        raise ValueError(
            "transpose_low_flow requires hydrogeologic_setting: a written assertion that "
            "donor and target share an aquifer or physiographic setting. Low flows are "
            "controlled by geology rather than area, this module cannot check the "
            "assertion itself, and leaving it unstated is what makes a transposed 7Q10 "
            "indefensible."
        )

    _require_kind(exponents, "non_exceedance", "transpose_low_flow")
    area_ratio, in_range = _check_areas(
        donor_area_sqmi,
        target_area_sqmi,
        exponents,
        area_ratio_range,
        allow_out_of_range,
    )

    donor_quantiles = results.quantiles
    if donor_quantiles is None or donor_quantiles.empty:
        raise ValueError("the donor LowFlowResults carries no quantiles")

    p_zero = float(results.p_zero or 0.0)
    if p_zero > max_p_zero:
        raise ValueError(
            f"the donor has zero flow in {p_zero:.1%} of years (p_zero={p_zero:g}, "
            f"max_p_zero={max_p_zero:g}). A ratio method cannot transpose a record that "
            "goes dry: the statistic is then controlled by whether the channel holds "
            "water at all, which is a property of the donor's geology, not its area."
        )

    donor_flows = donor_quantiles["flow_cfs"].to_numpy(dtype=float)
    if np.any(~(donor_flows > 0.0)):
        raise ValueError(
            "at least one donor low-flow statistic is zero or non-positive, and "
            "0 * (area ratio) is 0 -- which would assert that the target is dry rather "
            "than estimate it. Refusing; the target needs its own analysis or a "
            "regional low-flow regression."
        )

    probs = donor_quantiles["non_exceedance_prob"].to_numpy(dtype=float)
    exponent_frame = exponents.exponent_at(
        probs,
        interpolation=interpolation,
        extrapolation=extrapolation,
        b_bounds=b_bounds,
    )
    _warn_if_extrapolated(exponent_frame, len(probs))

    bfi_difference: Optional[float] = None
    if donor_bfi is not None and target_bfi is not None:
        bfi_difference = abs(float(donor_bfi) - float(target_bfi))
        if bfi_difference > max_bfi_difference:
            logger.warning(
                "Donor and target baseflow indices differ by %.3f (%.3f vs %.3f), more "
                "than max_bfi_difference=%.3f. Low-flow transposition assumes similar "
                "baseflow behaviour; this pair may not have it.",
                bfi_difference,
                donor_bfi,
                target_bfi,
                max_bfi_difference,
            )

    b_values = exponent_frame["exponent"].to_numpy(dtype=float)
    factors = area_ratio**b_values

    monotonic = _check_monotonic(
        probs,
        donor_flows * factors,
        flow_rises_with_prob=True,
        on_non_monotonic=on_non_monotonic,
        context="low-flow curve",
    )

    quantiles = pd.DataFrame(
        {
            "non_exceedance_prob": probs,
            "return_period": donor_quantiles["return_period"].to_numpy(dtype=float),
            "donor_flow_cfs": donor_flows,
            "exponent": b_values,
            "area_ratio_factor": factors,
            "flow_cfs": donor_flows * factors,
            "source": exponent_frame["source"].to_numpy(),
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
        statistic=f"{results.n_day}-day low flow",
        hydrogeologic_setting=hydrogeologic_setting,
        bfi_difference=bfi_difference,
        monotonic=monotonic,
    )

    return TransposedResults(
        quantiles=quantiles,
        confidence_limits=pd.DataFrame(),
        provenance=provenance,
    )
