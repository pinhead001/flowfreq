"""Regional (generalized) skew lookup for Bulletin 17C.

Bulletin 17C weights the at-site skew with a regional skew and its MSE. Until
now flowfreq took that pair as two bare floats from the caller. This module
gives them a versioned table, with a citation for every value, and a lookup
that raises rather than guessing.

Two rules the table enforces:

- A row is usable only when ``status == "verified"``, meaning it has been
  transcribed from its cited report and second-checked. ``pending`` rows exist
  so that the gaps are visible, and they are never returned as values.
- There is no silent fallback to the Bulletin 17B Plate I map. A caller who
  wants a fallback must pass one explicitly, and it is labelled as such.

Roadmap: ``docs/MASTER_ROADMAP.md`` §1.3, issue #34.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union

import pandas as pd
import requests

logger = logging.getLogger(__name__)

DEFAULT_TABLE = Path(__file__).parent / "data" / "regional_skew.csv"
COLUMNS = (
    "state",
    "skew_region",
    "skew",
    "skew_mse",
    "effective_record_years",
    "citation",
    "table",
    "status",
    "notes",
)
STATUSES = frozenset({"pending", "verified"})


class RegionalSkewUnavailable(LookupError):
    """No verified regional skew for the requested state or region."""


@dataclass(frozen=True)
class RegionalSkew:
    """A regional skew with its uncertainty and source."""

    state: str
    skew_region: str
    skew: float
    skew_mse: float
    citation: str
    table: str = ""
    effective_record_years: Optional[float] = None
    #: Set when the value comes from a study whose region covers the point but
    #: which the state itself has not adopted (e.g. the Pacific Northwest study
    #: in western Montana). Empty otherwise.
    advisory: str = ""

    def __post_init__(self) -> None:
        if not self.citation.strip():
            raise ValueError("RegionalSkew requires a citation")
        if not math.isfinite(self.skew):
            raise ValueError(f"Non-finite skew {self.skew}")
        if not (self.skew_mse > 0 and math.isfinite(self.skew_mse)):
            raise ValueError(f"skew_mse must be positive and finite, got {self.skew_mse}")


def load_table(path: Union[str, Path, None] = None) -> pd.DataFrame:
    """Load and validate the regional skew table.

    Parameters
    ----------
    path : str or Path, optional
        Defaults to the packaged ``flowfreq/data/regional_skew.csv``.

    Returns
    -------
    pandas.DataFrame

    Raises
    ------
    ValueError
        On missing columns, an unknown status, a duplicate (state, region)
        pair, or a verified row missing a value or citation.
    """
    df = pd.read_csv(path or DEFAULT_TABLE, dtype={"state": str, "skew_region": str, "table": str})
    missing = [c for c in COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"regional skew table missing columns: {missing}")
    bad = set(df["status"]) - STATUSES
    if bad:
        raise ValueError(f"unknown status value(s): {sorted(bad)}")
    dup = df.duplicated(subset=["state", "skew_region"])
    if dup.any():
        raise ValueError(
            f"duplicate (state, skew_region): {df.loc[dup, ['state', 'skew_region']].values.tolist()}"
        )
    verified = df[df["status"] == "verified"]
    for col in ("skew", "skew_mse", "citation"):
        if verified[col].isna().any():
            raise ValueError(f"verified row(s) missing {col}")
    return df


def available(state: str, path: Union[str, Path, None] = None) -> List[str]:
    """List the verified skew regions for a state.

    Parameters
    ----------
    state : str
        Two-letter postal code.
    path : str or Path, optional
        Table override.

    Returns
    -------
    list of str
    """
    df = load_table(path)
    rows = df[(df["state"] == state.upper()) & (df["status"] == "verified")]
    return sorted(rows["skew_region"].tolist())


def regional_skew_for(
    state: str,
    skew_region: Optional[str] = None,
    *,
    path: Union[str, Path, None] = None,
) -> RegionalSkew:
    """Return the verified regional skew for a state and skew region.

    Parameters
    ----------
    state : str
        Two-letter postal code.
    skew_region : str, optional
        Required whenever the state has more than one skew-region row, of
        *any* status. A pending row marks an area the verified value does not
        cover (e.g. Idaho's Snake River Plain), so returning the verified row
        for an unnamed region could apply it where its study says it is invalid.
    path : str or Path, optional
        Table override.

    Returns
    -------
    RegionalSkew

    Raises
    ------
    RegionalSkewUnavailable
        If there is no verified row, or the region is ambiguous.
    """
    df = load_table(path)
    st = state.upper()
    rows = df[df["state"] == st]
    if skew_region is None and len(rows) > 1:
        listing = ", ".join(f"{r.skew_region!r} ({r.status})" for r in rows.itertuples())
        raise RegionalSkewUnavailable(
            f"{st} has {len(rows)} skew regions: {listing}. Name one with skew_region -- "
            "a verified value does not apply inside a region still pending."
        )
    if skew_region is not None:
        rows = rows[rows["skew_region"] == skew_region]
    verified = rows[rows["status"] == "verified"]
    if verified.empty:
        pending = rows[rows["status"] == "pending"]
        why = "pending transcription" if not pending.empty else "no entry"
        raise RegionalSkewUnavailable(
            f"No verified regional skew for {st}"
            + (f" region {skew_region!r}" if skew_region else "")
            + f" ({why}). Supply regional_skew/regional_skew_mse explicitly."
        )
    if len(verified) > 1:
        raise RegionalSkewUnavailable(
            f"{st} has {len(verified)} skew regions {sorted(verified['skew_region'])}; name one"
        )
    r = verified.iloc[0]
    erl = r["effective_record_years"]
    return RegionalSkew(
        state=st,
        skew_region=str(r["skew_region"]),
        skew=float(r["skew"]),
        skew_mse=float(r["skew_mse"]),
        citation=str(r["citation"]),
        table="" if pd.isna(r["table"]) else str(r["table"]),
        effective_record_years=None if pd.isna(erl) else float(erl),
    )


# ----------------------------------------------------------------------------
# Location lookup: hydrologic-unit membership
# ----------------------------------------------------------------------------

#: Which hydrologic units each skew study covers, excludes, or leaves open.
#: Regenerate with ``tools/build_regional_skew_hucs.py``; never edit by hand.
HUC_TABLE = Path(__file__).parent / "data" / "regional_skew_hucs.csv"

#: States each HUC8 of the tabulated studies touches (Watershed Boundary Dataset).
HUC8_STATES_TABLE = Path(__file__).parent / "data" / "regional_skew_huc8_states.csv"

HUC_COLUMNS = ("huc", "skew_region", "rule", "valid_states", "source", "notes")
HUC_RULES = frozenset({"member", "excluded", "boundary", "unresolved"})

#: StreamStats layer holding Idaho's SIR 2016-5083 peak-flow region polygons.
#: The Snake River Plain ("region 0") is the part of Idaho outside all of
#: them, so a point query answering no ``Peak_Flow_Region_*_2016_5083``
#: feature is on the Plain. Used only to decide ``boundary`` HUs.
IDAHO_REGIONS_QUERY_URL = (
    "https://gis.streamstats.usgs.gov/arcgis/rest/services/nss/regions/MapServer/11/query"
)

#: USGS Watershed Boundary Dataset map service, 12-digit HU layer. Answers a
#: point query with the HU12 containing it. Live-verified 2026-09-27: the
#: point of 12340500 (Clark Fork above Missoula) returns 170102040104, the
#: HU12 the Water Data API gives as that site's ``hydrologic_unit_code``.
WBD_HU12_QUERY_URL = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/6/query"

#: FIPS state code to postal code, for Water Data API monitoring locations.
STATE_FIPS: Dict[str, str] = {
    "01": "AL", "02": "AK", "04": "AZ", "05": "AR", "06": "CA", "08": "CO", "09": "CT",
    "10": "DE", "11": "DC", "12": "FL", "13": "GA", "15": "HI", "16": "ID", "17": "IL",
    "18": "IN", "19": "IA", "20": "KS", "21": "KY", "22": "LA", "23": "ME", "24": "MD",
    "25": "MA", "26": "MI", "27": "MN", "28": "MS", "29": "MO", "30": "MT", "31": "NE",
    "32": "NV", "33": "NH", "34": "NJ", "35": "NM", "36": "NY", "37": "NC", "38": "ND",
    "39": "OH", "40": "OK", "41": "OR", "42": "PA", "44": "RI", "45": "SC", "46": "SD",
    "47": "TN", "48": "TX", "49": "UT", "50": "VT", "51": "VA", "53": "WA", "54": "WV",
    "55": "WI", "56": "WY", "60": "AS", "66": "GU", "69": "MP", "72": "PR", "78": "VI",
}  # fmt: skip


def _is_huc(huc: str) -> bool:
    return huc.isdigit() and len(huc) in (2, 4, 6, 8, 10, 12)


def load_huc_table(path: Union[str, Path, None] = None) -> pd.DataFrame:
    """Load and validate the hydrologic-unit definitions of the skew studies.

    Parameters
    ----------
    path : str or Path, optional
        Defaults to the packaged ``flowfreq/data/regional_skew_hucs.csv``.

    Returns
    -------
    pandas.DataFrame

    Raises
    ------
    ValueError
        On missing columns, an unknown rule, a malformed or duplicated HUC, a
        member row without valid states, or a row without a source.
    """
    df = pd.read_csv(path or HUC_TABLE, dtype=str, keep_default_na=False)
    missing = [c for c in HUC_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"regional skew HUC table missing columns: {missing}")
    bad = set(df["rule"]) - HUC_RULES
    if bad:
        raise ValueError(f"unknown rule value(s): {sorted(bad)}")
    malformed = [h for h in df["huc"] if not _is_huc(h)]
    if malformed:
        raise ValueError(f"malformed HUC(s): {malformed[:5]}")
    if df["huc"].duplicated().any():
        raise ValueError(f"duplicate HUC(s): {sorted(df.loc[df['huc'].duplicated(), 'huc'])}")
    if (df.loc[df["rule"] == "member", "valid_states"].str.strip() == "").any():
        raise ValueError("member rows need valid_states")
    if (df["source"].str.strip() == "").any():
        raise ValueError("every HUC row needs a source")
    return df


def _states_of(huc: str, path: Union[str, Path, None]) -> Set[str]:
    """States the HUC touches, from the packaged HUC8 table; empty if unknown."""
    df = pd.read_csv(path or HUC8_STATES_TABLE, dtype=str, keep_default_na=False)
    if len(huc) < 8:
        rows = df[df["huc8"].str.startswith(huc)]
    else:
        rows = df[df["huc8"] == huc[:8]]
    return {s for states in rows["states"] for s in states.split(",") if s}


def regional_skew_at_huc(
    huc: str,
    state: Optional[str] = None,
    *,
    path: Union[str, Path, None] = None,
    huc_path: Union[str, Path, None] = None,
    states_path: Union[str, Path, None] = None,
    point: Optional[Tuple[float, float]] = None,
    timeout: int = 60,
) -> RegionalSkew:
    """Return the verified regional skew for a hydrologic unit.

    Skew studies define their regions by drainage, so membership is decided
    by hydrologic unit: the most specific (longest) prefix of ``huc`` listed
    in ``regional_skew_hucs.csv`` decides. Nothing is defaulted -- a unit no
    verified study covers raises.

    Parameters
    ----------
    huc : str
        Hydrologic unit code, 2 to 12 digits. A 12-digit HU is the reliable
        input: a coarser unit that is only partly inside a study (or partly
        inside one of its exclusions) raises rather than guess.
    state : str, optional
        Two-letter postal code of the site. Needed only where the HU crosses
        into a state the study does not cover (e.g. Snake headwaters units
        shared by Idaho and Wyoming); otherwise inferred from the HUC8.
    path : str or Path, optional
        Override for the regional skew value table.
    huc_path, states_path : str or Path, optional
        Overrides for the HUC membership and HUC8-states tables.
    point : (float, float), optional
        The site's ``(lat, lon)``. Only consulted for a ``boundary`` HU -- one
        the Snake River Plain exclusion covers in part -- which is then decided
        by a live query of the SIR 2016-5083 region polygons
        (:data:`IDAHO_REGIONS_QUERY_URL`). Without it such an HU raises.
    timeout : int
        Request timeout for that query, seconds.

    Returns
    -------
    RegionalSkew
        With ``advisory`` set when the study covers the site but the state has
        no verified row adopting it (western Montana: USGS Montana applies the
        Bulletin 17B map statewide instead; see the MT row's notes).

    Raises
    ------
    ValueError
        ``huc`` is not 2-12 digits of even length.
    RegionalSkewUnavailable
        No verified study covers the HU; it is excluded by its study (the
        Snake River Plain), unresolved by the report, only partly covered at
        the resolution given, or in a state the study does not cover.
    """
    huc = str(huc).strip()
    if not _is_huc(huc):
        raise ValueError(f"huc must be 2-12 digits of even length, got {huc!r}")
    st = state.upper() if state else None
    hucs = load_huc_table(huc_path)

    matches = hucs[[huc.startswith(h) for h in hucs["huc"]]]
    # A listed unit *inside* the query with a different outcome makes a coarse
    # query ambiguous (e.g. a HUC8 part of which is Snake River Plain).
    finer = hucs[hucs["huc"].str.startswith(huc) & (hucs["huc"].str.len() > len(huc))]
    if matches.empty and finer.empty:
        raise RegionalSkewUnavailable(
            f"No verified regional skew study covers HUC {huc}. Supply "
            "regional_skew/regional_skew_mse explicitly (Bulletin 17C p. 31: consult the USGS)."
        )
    best = matches.loc[matches["huc"].str.len().idxmax()] if not matches.empty else None
    outcomes = {(str(r.skew_region), str(r.rule)) for r in finer.itertuples()}
    if best is not None:
        outcomes.discard((str(best["skew_region"]), str(best["rule"])))
    if outcomes or best is None:
        listing = ", ".join(sorted(f"{r} ({rule})" for r, rule in outcomes))
        raise RegionalSkewUnavailable(
            f"HUC {huc} is only partly covered: units inside it are {listing}. "
            "Give the site's 12-digit HUC."
        )

    region, rule = str(best["skew_region"]), str(best["rule"])
    table = load_table(path)
    if rule == "boundary":
        if point is None:
            raise RegionalSkewUnavailable(
                f"HUC {huc} straddles the edge of the {region} ({best['notes']}); its code "
                "alone cannot decide. Use regional_skew_at(lat, lon) with the site location."
            )
        if _on_snake_river_plain(*point, timeout=timeout):
            rule = "excluded"
        else:
            outer = matches[matches["rule"] != "boundary"]
            if outer.empty:
                raise RegionalSkewUnavailable(f"No verified regional skew study covers HUC {huc}")
            best = outer.loc[outer["huc"].str.len().idxmax()]
            region, rule = str(best["skew_region"]), str(best["rule"])
    if rule == "excluded":
        rows = table[table["skew_region"] == region]
        why = str(rows["notes"].iloc[0]) if not rows.empty else str(best["notes"])
        raise RegionalSkewUnavailable(
            f"HUC {huc} is in the {region}, where no regional skew applies: {why}"
        )
    if rule == "unresolved":
        raise RegionalSkewUnavailable(
            f"HUC {huc} (unit {best['huc']}) is at the edge of the {region} study and "
            f"unresolved: {best['notes']}"
        )

    valid = set(str(best["valid_states"]).split(";"))
    touched = _states_of(huc, states_path)
    if st is None:
        if touched and not touched & valid:
            raise RegionalSkewUnavailable(
                f"HUC {huc} lies in {', '.join(sorted(touched))}, outside the {region} study "
                f"(which covers {', '.join(sorted(valid))})."
            )
        if not touched or not touched <= valid:
            raise RegionalSkewUnavailable(
                f"HUC {huc} touches {', '.join(sorted(touched)) or 'unknown states'}; the "
                f"{region} study covers only {', '.join(sorted(valid))}. Pass state=."
            )
        if len(touched) == 1:
            st = next(iter(touched))
    elif st not in valid:
        raise RegionalSkewUnavailable(
            f"{st} is outside the {region} study (which covers {', '.join(sorted(valid))})."
        )

    verified = table[(table["skew_region"] == region) & (table["status"] == "verified")]
    if verified.empty:
        raise RegionalSkewUnavailable(f"No verified row for skew region {region!r}")
    values = verified[["skew", "skew_mse", "effective_record_years"]].drop_duplicates()
    if len(values) != 1:
        raise ValueError(f"verified rows for {region!r} disagree: {values.values.tolist()}")
    own = verified[verified["state"] == st] if st else verified.iloc[0:0]
    advisory = ""
    if own.empty:
        r = verified.sort_values("state").iloc[0]
        unadopted = sorted({st} if st else touched - set(verified["state"]))
        if unadopted:
            advisory = (
                f"{', '.join(unadopted)} has no verified row adopting the {region} study; the "
                f"value is the study's own. Read the {', '.join(unadopted)} notes in "
                "regional_skew.csv before using it."
            )
            logger.warning("HUC %s: %s", huc, advisory)
    else:
        r = own.iloc[0]
    erl = r["effective_record_years"]
    return RegionalSkew(
        state=st or ",".join(sorted(touched)),
        skew_region=region,
        skew=float(r["skew"]),
        skew_mse=float(r["skew_mse"]),
        citation=str(r["citation"]),
        table="" if pd.isna(r["table"]) else str(r["table"]),
        effective_record_years=None if pd.isna(erl) else float(erl),
        advisory=advisory,
    )


def regional_skew_for_site(site_no: str, *, timeout: int = 60) -> RegionalSkew:
    """Regional skew for a USGS monitoring location, from its own HUC.

    Uses the Water Data API monitoring-locations record
    (:func:`flowfreq.waterdata.fetch_monitoring_location`): its 12-digit
    ``hydrologic_unit_code`` and its ``state_code``.

    Parameters
    ----------
    site_no : str
        USGS site number.
    timeout : int
        Request timeout, seconds.

    Returns
    -------
    RegionalSkew

    Raises
    ------
    RegionalSkewUnavailable
        The site has no HUC on record, or :func:`regional_skew_at_huc` raises.
    requests.RequestException
        The monitoring-location request failed.
    """
    from flowfreq.waterdata import fetch_monitoring_location

    props = fetch_monitoring_location(str(site_no), timeout=timeout)
    huc = str(props.get("hydrologic_unit_code") or "").strip()
    if not huc:
        raise RegionalSkewUnavailable(f"Site {site_no} has no hydrologic_unit_code on record")
    state = STATE_FIPS.get(str(props.get("state_code") or "").zfill(2))
    return regional_skew_at_huc(huc, state)


def huc12_at(lat: float, lon: float, *, timeout: int = 60) -> Tuple[str, List[str]]:
    """The 12-digit hydrologic unit containing a point, and the states it touches.

    Queries the USGS Watershed Boundary Dataset map service
    (:data:`WBD_HU12_QUERY_URL`).

    Parameters
    ----------
    lat, lon : float
        WGS84 decimal degrees.
    timeout : int
        Request timeout, seconds.

    Returns
    -------
    (str, list of str)
        HU12 code and its WBD ``states`` (postal codes; ``CN`` is Canada).

    Raises
    ------
    ValueError
        Coordinates out of range.
    RegionalSkewUnavailable
        No HU12 contains the point.
    requests.RequestException
        The request failed or the service returned an error.
    """
    if not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0):
        raise ValueError(f"({lat}, {lon}) is not a WGS84 latitude/longitude")
    params = {
        "geometry": f"{lon},{lat}",
        "geometryType": "esriGeometryPoint",
        "inSR": "4326",
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "huc12,states",
        "returnGeometry": "false",
        "f": "json",
    }
    try:
        response = requests.get(WBD_HU12_QUERY_URL, params=params, timeout=timeout)
        response.raise_for_status()
        payload = response.json()
    except (requests.RequestException, ValueError) as exc:
        raise requests.RequestException(f"WBD HU12 query failed at ({lat}, {lon}): {exc}") from exc
    if not isinstance(payload, dict) or "error" in payload:
        raise requests.RequestException(f"WBD HU12 query error at ({lat}, {lon}): {payload}")
    feats = payload.get("features") or []
    if not feats:
        raise RegionalSkewUnavailable(f"No hydrologic unit contains ({lat}, {lon})")
    attrs = feats[0].get("attributes") or {}
    states = [s for s in str(attrs.get("states") or "").split(",") if s]
    return str(attrs.get("huc12")), states


def regional_skew_at(
    lat: float, lon: float, state: Optional[str] = None, *, timeout: int = 60
) -> RegionalSkew:
    """Regional skew for a point, via the Watershed Boundary Dataset HU12 there.

    The point's HU12 decides membership (:func:`regional_skew_at_huc`); an
    HU12 only partly on the Snake River Plain is decided by the point itself.

    Parameters
    ----------
    lat, lon : float
        WGS84 decimal degrees of the site (the gage or basin outlet).
    state : str, optional
        Postal code of the site. Taken from the HU12 when it lies in one
        state; needed only for an HU12 on a state line the study stops at.
    timeout : int
        Request timeout, seconds.

    Returns
    -------
    RegionalSkew

    Raises
    ------
    RegionalSkewUnavailable
        No verified study covers the point (see :func:`regional_skew_at_huc`),
        or its HU12 crosses into a state the study does not cover.
    requests.RequestException
        The WBD request failed.
    """
    huc, states = huc12_at(lat, lon, timeout=timeout)
    if state is None and len(states) == 1:
        state = states[0]
    return regional_skew_at_huc(huc, state, point=(lat, lon), timeout=timeout)


def _on_snake_river_plain(lat: float, lon: float, *, timeout: int = 60) -> bool:
    """True if the point is in none of Idaho's SIR 2016-5083 peak-flow regions."""
    params = {
        "geometry": f"{lon},{lat}",
        "geometryType": "esriGeometryPoint",
        "inSR": "4326",
        "spatialRel": "esriSpatialRelIntersects",
        "where": "NAME LIKE 'Peak_Flow_Region%2016_5083'",
        "outFields": "NAME",
        "returnGeometry": "false",
        "f": "json",
    }
    try:
        response = requests.get(IDAHO_REGIONS_QUERY_URL, params=params, timeout=timeout)
        response.raise_for_status()
        payload = response.json()
    except (requests.RequestException, ValueError) as exc:
        raise requests.RequestException(
            f"StreamStats Idaho region query failed at ({lat}, {lon}): {exc}"
        ) from exc
    if not isinstance(payload, dict) or "error" in payload or "features" not in payload:
        raise requests.RequestException(
            f"StreamStats Idaho region query error at ({lat}, {lon}): {payload}"
        )
    return not payload["features"]
