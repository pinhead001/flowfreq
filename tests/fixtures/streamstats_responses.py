"""Captured-shape StreamStats API payloads for delineation/characteristics tests.

Values are plausible for the Methow River basin, WA (the design doc's own
verification points, ``docs/STREAMSTATS_MODULE_DESIGN.md`` appendix). The snap and
basin-characteristics shapes below are confirmed against the live service (2026-09-11
for the snap shape, 2026-09-09/10 for basin characteristics per the design doc); the
``ss-delineate``/``ss-hydro`` chaining shapes are this module's best-effort
interpretation of the documented protocol and not independently re-verified here.
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

# A successful ss-delineate 'features' response for a good point: no WarningMsg, a
# real (if simplified) watershed polygon on the globalwatershed feature.
DELINEATE_FEATURES_GOOD = {
    "type": "FeatureCollection",
    "features": [
        {
            "id": "globalwatershed",
            "properties": {"Name": "globalwatershed", "WarningMsg": ""},
            "geometry": {
                "type": "Polygon",
                "coordinates": [
                    [
                        [-120.40, 48.55],
                        [-120.35, 48.55],
                        [-120.35, 48.60],
                        [-120.40, 48.60],
                        [-120.40, 48.55],
                    ]
                ],
            },
        }
    ],
}

# The design doc's own captured failure mode (S4): HTTP 200, a structurally valid
# hillslope-sliver polygon, and the exact WarningMsg string it recorded live. This
# fixture -- and the test asserting it raises -- is the whole point of this module.
DELINEATE_FEATURES_UNSNAPPABLE_WARNING = {
    "type": "FeatureCollection",
    "features": [
        {
            "id": "globalwatershed",
            "properties": {
                "Name": "globalwatershed",
                "WarningMsg": (
                    ", Point not snappable using ss-pourpoint API service; "
                    "results may be inaccurate."
                ),
            },
            "geometry": {
                "type": "Polygon",
                "coordinates": [
                    [
                        [-120.3701, 48.5840],
                        [-120.3700, 48.5840],
                        [-120.3700, 48.5841],
                        [-120.3701, 48.5840],
                    ]
                ],
            },
        }
    ],
}

# No 'features' key at all -- an unexpected response shape, not a warned delineation.
DELINEATE_FEATURES_MALFORMED = {"type": "FeatureCollection"}

# A degenerate polygon with no WarningMsg to explain it -- must still be refused
# (FR-3's validation does not depend on the service always labeling its own failures).
DELINEATE_FEATURES_DEGENERATE_NO_WARNING = {
    "type": "FeatureCollection",
    "features": [
        {
            "id": "globalwatershed",
            "properties": {"Name": "globalwatershed", "WarningMsg": ""},
            "geometry": {"type": "Polygon", "coordinates": [[[-120.37, 48.58]] * 4]},
        }
    ],
}

# The ss-delineate 'sshydro' chaining response: the bcrequest body to POST to ss-hydro.
DELINEATE_SSHYDRO_GOOD = {
    "stateAbbreviation": "WA",
    "bcrequest": {
        "cid": "WA20260910120000123",
        "rcode": "WA",
        "workspaceID": "WA20260910120000123",
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
