"""Captured-shape StreamStats API payloads for delineation/characteristics tests.

Values are plausible for the Methow River basin, WA (the design doc's own
verification points, ``docs/STREAMSTATS_MODULE_DESIGN.md`` appendix). The snap shape is
confirmed live (2026-09-11); basin characteristics per the design doc (2026-09-09/10).
The ``ss-delineate`` ``sshydro`` chaining response's overall shape
(``{"stateAbbreviation", "bcrequest"}``) is confirmed live; the exact nesting of a
``WarningMsg`` within it is not, hence ``_find_warning_msg``'s recursive scan rather
than a fixed path.
"""

# A successful pourpoint snap: the point lies close to the stream network. `output` is
# a GeoJSON Point -- coordinates are [lon, lat], confirmed live 2026-09-11.
SNAP_GOOD = {
    "region": "WA",
    "input": {"type": "Point", "coordinates": [-120.37893, 48.57426]},
    "output": {"type": "Point", "coordinates": [-120.3789336506848, 48.57425501781417]},
    "couldSnap": True,
}

# The point is off the flowline -- the service will not snap it. The whole point of
# this fixture is that FR-1 forbids proceeding to delineation from here.
SNAP_UNSNAPPABLE = {
    "couldSnap": False,
    "output": {},
}

# The ss-delineate 'sshydro' chaining response: the bcrequest body to POST to ss-hydro.
# This is the one ss-delineate call the shipped pipeline actually makes -- an earlier
# version also called `delineate/features/{region}` for a watershed polygon, but that
# was found live to return an unrelated point feature, not a polygon, and was removed;
# see flowfreq/streamstats.py's module docstring.
DELINEATE_SSHYDRO_GOOD = {
    "stateAbbreviation": "WA",
    "bcrequest": {
        "cid": "WA20260910120000123",
        "rcode": "WA",
        "workspaceID": "WA20260910120000123",
    },
}

# The design doc's own captured failure mode (S4): HTTP 200 carrying the exact
# WarningMsg string it recorded live for an unsnappable point. The precise nesting of
# WarningMsg within this chaining variant's response was not independently
# re-confirmed (module docstring), so it's placed inside bcrequest here defensively --
# _find_warning_msg scans the whole structure, not one fixed path. This fixture -- and
# the test asserting it raises -- is the whole point of this module.
DELINEATE_SSHYDRO_WITH_WARNING = {
    "stateAbbreviation": "WA",
    "bcrequest": {
        "cid": "WA20260910120000124",
        "rcode": "WA",
        "workspaceID": "WA20260910120000124",
        "WarningMsg": (
            ", Point not snappable using ss-pourpoint API service; " "results may be inaccurate."
        ),
    },
}

# No bcrequest payload -- ss-hydro has nothing usable to compute from.
DELINEATE_SSHYDRO_MALFORMED = {"stateAbbreviation": "WA"}

# ss-hydro basin-characteristics response for tp_goat_creek (design doc appendix):
# DRNAREA 412.0 sq mi, PRECPRIS10 45.62 in, CANOPY_PCT 45.242 %.
HYDRO_CHARACTERISTICS_GOOD = [
    {
        "code": "DRNAREA",
        "name": "Drainage Area",
        "description": "Area that drains to a point on a stream",
        "value": 412.0,
        "unit": "square miles",
        "msg": "",
    },
    {
        "code": "PRECPRIS10",
        "name": "Mean Annual Precipitation",
        "description": "Mean annual precipitation, PRISM 1971-2000",
        "value": 45.62,
        "unit": "inches",
        "msg": "",
    },
    {
        "code": "CANOPY_PCT",
        "name": "Percent Canopy",
        "description": "Percentage of canopy cover",
        "value": 45.242,
        "unit": "percent",
        "msg": "Value computed from 2016 NLCD canopy layer",
    },
]

# The same characteristics wrapped in a {"parameters": [...]} envelope, the shape the
# ff-idea02 transcript's script expected -- supported defensively alongside the design
# doc's bare-list shape.
HYDRO_CHARACTERISTICS_WRAPPED = {"parameters": HYDRO_CHARACTERISTICS_GOOD}

# Neither a list nor a {"parameters": [...]} dict -- an unexpected response shape.
HYDRO_CHARACTERISTICS_MALFORMED = {"error": "internal error", "code": 500}

# A characteristic entry missing its value -- present but unusable.
HYDRO_CHARACTERISTICS_MISSING_VALUE = [
    {"code": "DRNAREA", "name": "Drainage Area", "description": "", "unit": "square miles"}
]

# Body of the 422 ss-hydro (or snap) returns when a required parameter -- here,
# region -- is missing from the request.
RESPONSE_422_MISSING_REGION = '{"errors": ["region is required"]}'
