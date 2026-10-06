"""Nevada-specific regression procedures that the equation schema does not express.

Nevada's peak-flow equations are those of Thomas, Hjalmarson and Waltemeyer (1997,
USGS Water-Supply Paper 2433). The report defines two transition zones in which an
estimate is a weighted average of two regions' estimates ("Transition Zones", p. 19;
worked examples pp. 39-44). Both average flows, not logarithms, so both delegate to
:func:`flowfreq.regression.evaluate_weighted` with ``space="linear"``.

**Drainage area in two low- to middle-elevation regions** (equation 6)::

    QT(w) = (QT(a) * AREA(a) + QT(b) * AREA(b)) / AREA

Each region's equation is evaluated with the characteristics of the *entire* basin,
and the results are weighted by the drainage area in each region.

**Site near the High-Elevation Region 1 boundary** (equation 7)::

    QT(w) = QT(l) * (B - E) / 700  +  QT(h) * (1 - (B - E) / 700)

where ``E`` is the elevation of the study *site* (not the mean basin elevation), ``B``
the elevation of the lower boundary of High-Elevation Region 1, ``QT(l)`` the estimate
for the low- to middle-elevation region the site is in and ``QT(h)`` the Region 1
estimate, both with the characteristics of the entire basin. The zone extends 700 ft
below ``B``. A site at or above ``B`` is in Region 1 (p. 16: the site elevation decides
whether a site is in the high-elevation region).

``B`` comes from figure 5 (p. 12). The text gives it numerically only south of 41°
latitude: "South of 41° latitude, all study sites with an elevation between 6,800 and
7,500 ft are in the transition zone" (p. 19), so ``B`` = 7,500 ft there. North of 41°
the boundary falls with latitude, but the report gives it only as a line on figure 5,
so :func:`high_elevation_boundary_ft` refuses a latitude north of 41° and the caller
must read ``B`` off figure 5 and pass it explicitly.

The report's third transition zone (a site near a boundary whose basin lies in one
region, for which "a straight average ... may be appropriate") is a judgement call,
not a formula, and is not implemented. The report gives no prediction interval for a
weighted estimate; each component keeps its own.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Mapping, Optional

from flowfreq.regression.equations import RegressionEstimate, evaluate, evaluate_weighted
from flowfreq.regression.library import StateLibrary, load_state

#: Region code of High-Elevation Region 1.
HIGH_ELEVATION_REGION = "1"

#: Width of the elevation transition zone below the Region 1 boundary, in feet (eq. 7).
TRANSITION_ZONE_FT = 700.0

#: Lower boundary of High-Elevation Region 1 south of 41° latitude, in feet (p. 19).
BOUNDARY_SOUTH_OF_41_FT = 7500.0

#: Latitude, in decimal degrees, north of which the boundary is given only graphically.
BOUNDARY_KNEE_LATITUDE = 41.0


def high_elevation_boundary_ft(latitude: float) -> float:
    """Lower boundary ``B`` of High-Elevation Region 1 at a latitude, in feet.

    Parameters
    ----------
    latitude : float
        Latitude of the study site, decimal degrees.

    Returns
    -------
    float
        7,500 ft at or south of 41° latitude (WSP 2433, p. 19).

    Raises
    ------
    ValueError
        North of 41°, where the report gives ``B`` only as a line on figure 5: read it
        from the figure and pass ``boundary_ft`` to the estimator instead. Also if the
        latitude is not finite.
    """
    if not math.isfinite(latitude):
        raise ValueError(f"latitude must be finite, got {latitude}")
    if latitude > BOUNDARY_KNEE_LATITUDE:
        raise ValueError(
            f"latitude {latitude} is north of 41°: WSP 2433 gives the High-Elevation "
            "Region 1 boundary there only graphically (figure 5); pass boundary_ft"
        )
    return BOUNDARY_SOUTH_OF_41_FT


def elevation_transition_weights(site_elevation: float, boundary_ft: float) -> Dict[str, float]:
    """Equation 7 weights on the low- to middle-elevation region and on Region 1.

    Parameters
    ----------
    site_elevation : float
        Elevation ``E`` of the study site, in feet.
    boundary_ft : float
        Lower boundary ``B`` of High-Elevation Region 1, in feet.

    Returns
    -------
    dict
        ``{"low": w, "high": 1 - w}`` with ``w = (B - E) / 700`` clipped to ``[0, 1]``:
        1 at or below ``B - 700`` (low region alone), 0 at or above ``B`` (Region 1 alone).

    Raises
    ------
    ValueError
        If either elevation is not finite.
    """
    if not (math.isfinite(site_elevation) and math.isfinite(boundary_ft)):
        raise ValueError(f"elevations must be finite, got E={site_elevation}, B={boundary_ft}")
    w_low = min(1.0, max(0.0, (boundary_ft - site_elevation) / TRANSITION_ZONE_FT))
    return {"low": w_low, "high": 1.0 - w_low}


@dataclass(frozen=True)
class WeightedEstimate:
    """A Nevada estimate weighted across two regions (equation 6 or 7).

    Attributes
    ----------
    aep : float
    flow_cfs : float
        The weighted flow, or the single contributing region's flow.
    weights : dict of str to float
        Weight on each region code; they sum to 1.
    estimates : dict of str to RegressionEstimate
        The component estimates, only for regions with a positive weight.
    method : str
        ``"area (eq. 6)"`` or ``"site elevation (eq. 7)"``.
    """

    aep: float
    flow_cfs: float
    weights: Dict[str, float]
    estimates: Dict[str, RegressionEstimate]
    method: str

    @property
    def blended(self) -> bool:
        """Whether more than one region contributed."""
        return len(self.estimates) > 1

    @property
    def extrapolated(self) -> bool:
        """Whether any contributing estimate used an input outside its calibrated range."""
        return any(e.extrapolated for e in self.estimates.values())


def _library(lib: Optional[StateLibrary]) -> StateLibrary:
    lib = lib if lib is not None else load_state("NV")
    if lib.state != "NV":
        raise ValueError(f"needs the Nevada library, got {lib.state}")
    return lib


def _weighted(
    lib: StateLibrary,
    aep: float,
    weights: Mapping[str, float],
    characteristics: Mapping[str, float],
    method: str,
    allow_extrapolation: bool,
    interval_level: float,
) -> WeightedEstimate:
    estimates = {
        region: evaluate(
            lib.equation(region, aep),
            characteristics,
            allow_extrapolation=allow_extrapolation,
            interval_level=interval_level,
        )
        for region, w in weights.items()
        if w > 0.0
    }
    flow = evaluate_weighted([(weights[r], e) for r, e in estimates.items()], space="linear")
    return WeightedEstimate(
        aep=aep, flow_cfs=flow, weights=dict(weights), estimates=estimates, method=method
    )


def estimate_area_weighted(
    aep: float,
    region_areas: Mapping[str, float],
    characteristics: Mapping[str, float],
    *,
    lib: Optional[StateLibrary] = None,
    allow_extrapolation: bool = False,
    interval_level: float = 0.90,
) -> WeightedEstimate:
    """Equation 6: a basin whose drainage area lies in two low- to middle-elevation regions.

    Parameters
    ----------
    aep : float
        Annual exceedance probability (0.5 to 0.01; the report has no 200- or 500-year
        equations, and Region 6's 2-year equation is printed ``Q=0``, so is not stored).
    region_areas : mapping of str to float
        Drainage area in each of the two regions, square miles, keyed by region code.
    characteristics : mapping of str to float
        Characteristics of the **entire** basin, used for both regions.
    lib : StateLibrary, optional
        The Nevada library; ``load_state("NV")`` when omitted.
    allow_extrapolation : bool, default False
    interval_level : float, default 0.90
        Passed to :func:`flowfreq.regression.evaluate` for each component.

    Returns
    -------
    WeightedEstimate

    Raises
    ------
    ValueError
        If there are not exactly two regions, one is High-Elevation Region 1 (use
        :func:`estimate_elevation_transition`), an area is negative or the total is not
        positive, or ``lib`` is not Nevada's.
    EquationsUnavailable
        If a region has no equation for the AEP.
    """
    lib = _library(lib)
    if len(region_areas) != 2:
        raise ValueError(f"equation 6 weights exactly two regions, got {sorted(region_areas)}")
    if HIGH_ELEVATION_REGION in region_areas:
        raise ValueError(
            "equation 6 is for two low- to middle-elevation regions; near High-Elevation "
            "Region 1 use estimate_elevation_transition (equation 7)"
        )
    areas = {r: float(a) for r, a in region_areas.items()}
    if any(not (math.isfinite(a) and a >= 0.0) for a in areas.values()):
        raise ValueError(f"region areas must be finite and non-negative, got {areas}")
    total = sum(areas.values())
    if not total > 0.0:
        raise ValueError("total drainage area must be positive")
    weights = {r: a / total for r, a in areas.items()}
    return _weighted(
        lib, aep, weights, characteristics, "area (eq. 6)", allow_extrapolation, interval_level
    )


def estimate_elevation_transition(
    aep: float,
    low_region: str,
    characteristics: Mapping[str, float],
    *,
    site_elevation: float,
    latitude: Optional[float] = None,
    boundary_ft: Optional[float] = None,
    lib: Optional[StateLibrary] = None,
    allow_extrapolation: bool = False,
    interval_level: float = 0.90,
) -> WeightedEstimate:
    """Equation 7: a site in a low- to middle-elevation region near the Region 1 boundary.

    Parameters
    ----------
    aep : float
        Annual exceedance probability.
    low_region : str
        Code of the low- to middle-elevation region the site is in.
    characteristics : mapping of str to float
        Characteristics of the entire basin, used for both regions. Only the
        characteristics of the regions with a positive weight are required.
    site_elevation : float
        Elevation ``E`` of the study site, in feet (not the mean basin elevation).
    latitude : float, optional
        Latitude of the site, decimal degrees; gives ``B`` = 7,500 ft at or south of 41°.
    boundary_ft : float, optional
        ``B`` read from figure 5, in feet. Required north of 41°; overrides ``latitude``.
    lib : StateLibrary, optional
        The Nevada library; ``load_state("NV")`` when omitted.
    allow_extrapolation : bool, default False
    interval_level : float, default 0.90
        Passed to :func:`flowfreq.regression.evaluate` for each component.

    Returns
    -------
    WeightedEstimate
        Weights keyed by ``low_region`` and ``"1"``. At or below ``B - 700`` ft only the
        low region contributes; at or above ``B`` only Region 1 does.

    Raises
    ------
    ValueError
        If neither ``latitude`` nor ``boundary_ft`` is given, the latitude is north of
        41° without ``boundary_ft``, ``low_region`` is Region 1, or ``lib`` is not Nevada's.
    EquationsUnavailable
        If a contributing region has no equation for the AEP.
    """
    lib = _library(lib)
    if low_region == HIGH_ELEVATION_REGION:
        raise ValueError("low_region must be a low- to middle-elevation region, not Region 1")
    if boundary_ft is None:
        if latitude is None:
            raise ValueError("pass latitude (at or south of 41°) or boundary_ft (figure 5)")
        boundary_ft = high_elevation_boundary_ft(latitude)
    w = elevation_transition_weights(site_elevation, boundary_ft)
    weights = {low_region: w["low"], HIGH_ELEVATION_REGION: w["high"]}
    return _weighted(
        lib,
        aep,
        weights,
        characteristics,
        "site elevation (eq. 7)",
        allow_extrapolation,
        interval_level,
    )
