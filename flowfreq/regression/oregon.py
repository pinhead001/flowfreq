"""Oregon-specific regression procedures that the equation schema does not express.

Western Oregon Regions 2A and 2B (Cooper, 2005, USGS SIR 2005-5116) are split on a
watershed's mean elevation at 3,000 ft: Region 2A at or above it, Region 2B below.
Near the boundary the report blends the two regions' estimates by linear
interpolation on mean elevation across a transition zone of width ``W`` centred on
3,000 ft (report p. 31, equation 7)::

    QT = Q2B * ((3,000 + W/2) - E) / W  +  Q2A * (E - (3,000 - W/2)) / W

where ``E`` is the mean watershed elevation in feet and ``W`` the width of the
zone in feet of elevation. The report evaluated widths of 100 to 1,000 ft (Table
13) and selected ``W = 250`` ft, for which equation 7 becomes equation 8 (p. 44)::

    QT = Q2B * (3,125 - E) / 250  +  Q2A * (E - 2,875) / 250

Below 2,875 ft Region 2B applies alone and above 3,125 ft Region 2A alone.
The blend is of flows, not logarithms, which is exactly
:func:`flowfreq.regression.evaluate_weighted` with ``space="linear"`` and the
equation 7 weights, so :func:`estimate_region2` computes the weights and delegates
to it. The elevation split is not a regression variable, so nothing in ``OR.json``
can express it. That is why this lives in a helper rather than the schema.

The report gives no prediction interval for a blended estimate. Each component
estimate keeps its own interval.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Mapping, Optional

from flowfreq.regression.equations import RegressionEstimate, evaluate, evaluate_weighted
from flowfreq.regression.library import StateLibrary, load_state

#: Mean watershed elevation, in feet, dividing Region 2A (above) from 2B (below).
REGION2_BOUNDARY_FT = 3000.0

#: Transition-zone width, in feet of elevation, that the report selected (equation 8).
TRANSITION_WIDTH_FT = 250.0

#: StreamStats/NSS code for mean basin elevation, in feet.
ELEVATION_CODE = "ELEV"


def region2_weights(mean_elevation: float, width: float = TRANSITION_WIDTH_FT) -> Dict[str, float]:
    """Return the equation 7 weights on Regions 2A and 2B for a mean elevation.

    Parameters
    ----------
    mean_elevation : float
        Mean watershed elevation ``E``, in feet.
    width : float, default 250
        Transition-zone width ``W``, in feet of elevation, centred on 3,000 ft.
        The report's selected value is 250 (equation 8).

    Returns
    -------
    dict
        ``{"2A": w, "2B": 1 - w}``, with ``w = (E - (3,000 - W/2)) / W`` clipped to
        ``[0, 1]``: 0 at or below the zone's lower edge, 1 at or above its upper edge.

    Raises
    ------
    ValueError
        If ``mean_elevation`` is not finite, or ``width`` is not positive and finite.
    """
    if not math.isfinite(mean_elevation):
        raise ValueError(f"mean_elevation must be finite, got {mean_elevation}")
    if not (math.isfinite(width) and width > 0):
        raise ValueError(f"width must be positive and finite, got {width}")
    lower = REGION2_BOUNDARY_FT - width / 2.0
    w_2a = min(1.0, max(0.0, (mean_elevation - lower) / width))
    return {"2A": w_2a, "2B": 1.0 - w_2a}


@dataclass(frozen=True)
class Region2Estimate:
    """A Region 2A/2B estimate, blended when the site lies in the transition zone.

    Attributes
    ----------
    aep : float
    flow_cfs : float
        ``Q2B * weights["2B"] + Q2A * weights["2A"]``, or the single region's flow.
    mean_elevation : float
        Feet.
    width : float
        Transition-zone width used, in feet.
    weights : dict of str to float
        Equation 7 weights on ``"2A"`` and ``"2B"``; they sum to 1.
    estimates : dict of str to RegressionEstimate
        The component estimates, only for regions with a positive weight.
    """

    aep: float
    flow_cfs: float
    mean_elevation: float
    width: float
    weights: Dict[str, float]
    estimates: Dict[str, RegressionEstimate]

    @property
    def blended(self) -> bool:
        """Whether both regions contributed (the site is strictly inside the zone)."""
        return len(self.estimates) == 2

    @property
    def extrapolated(self) -> bool:
        """Whether any contributing estimate used an input outside its calibrated range."""
        return any(e.extrapolated for e in self.estimates.values())


def estimate_region2(
    aep: float,
    characteristics: Mapping[str, float],
    *,
    mean_elevation: Optional[float] = None,
    width: float = TRANSITION_WIDTH_FT,
    lib: Optional[StateLibrary] = None,
    allow_extrapolation: bool = False,
    interval_level: float = 0.90,
) -> Region2Estimate:
    """Estimate a western Oregon Region 2 peak flow, applying the transition-zone blend.

    Parameters
    ----------
    aep : float
        Annual exceedance probability (the report publishes 0.5, 0.2, 0.1, 0.04,
        0.02, 0.01 and 0.002).
    characteristics : mapping of str to float
        Basin characteristics keyed by StreamStats code. Only the characteristics
        of the regions with a positive weight are required.
    mean_elevation : float, optional
        Mean watershed elevation, in feet. Taken from ``characteristics["ELEV"]``
        when omitted.
    width : float, default 250
        Transition-zone width in feet of elevation (equation 7's ``W``). Leave it
        at the report's selected 250 ft (equation 8) unless reproducing Table 13.
    lib : StateLibrary, optional
        The Oregon library; ``load_state("OR")`` when omitted.
    allow_extrapolation : bool, default False
        Passed to :func:`flowfreq.regression.evaluate` for each component.
    interval_level : float, default 0.90
        Passed to :func:`flowfreq.regression.evaluate` for each component.

    Returns
    -------
    Region2Estimate

    Raises
    ------
    KeyError
        If no mean elevation is given, or a contributing region's characteristic is missing.
    ValueError
        If ``lib`` is not Oregon's, or the elevation or width is invalid.
    OutOfRangeError
        If a contributing region's input is out of range and extrapolation is not allowed.
    EquationsUnavailable
        If the AEP is not published.
    """
    lib = lib if lib is not None else load_state("OR")
    if lib.state != "OR":
        raise ValueError(f"estimate_region2 needs the Oregon library, got {lib.state}")
    if mean_elevation is None:
        if ELEVATION_CODE not in characteristics:
            raise KeyError(
                f"mean elevation needed to choose Region 2A or 2B: pass mean_elevation "
                f"or characteristics[{ELEVATION_CODE!r}]"
            )
        mean_elevation = float(characteristics[ELEVATION_CODE])
    weights = region2_weights(mean_elevation, width)
    estimates = {
        region: evaluate(
            lib.equation(region, aep),
            characteristics,
            allow_extrapolation=allow_extrapolation,
            interval_level=interval_level,
        )
        for region in ("2B", "2A")
        if weights[region] > 0.0
    }
    flow = evaluate_weighted([(weights[r], e) for r, e in estimates.items()], space="linear")
    return Region2Estimate(
        aep=aep,
        flow_cfs=flow,
        mean_elevation=mean_elevation,
        width=width,
        weights=weights,
        estimates=estimates,
    )
