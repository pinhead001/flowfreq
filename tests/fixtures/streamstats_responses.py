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

# --- Phase 2 (NSS) fixtures. Trimmed from real live responses captured 2026-09-11
# --- against WA / Goat Creek (docs/STREAMSTATS_NSS_ADDENDUM.md). Shapes are real;
# --- values are the real ones, region lists trimmed to what each test needs.

# GET /nssservices/statisticgroups (global, no region) -- trimmed to 3 of 52 entries.
NSS_STATISTIC_GROUPS_ALL = [
    {"id": 1, "name": "Physical Characteristics", "code": "PC", "defType": "BC"},
    {"id": 2, "name": "Peak-Flow Statistics", "code": "PFS", "defType": "FS"},
    {"id": 4, "name": "Low-Flow Statistics", "code": "LFS", "defType": "FS"},
]

# GET /nssservices/regions/WA/statisticgroups -- WA supports only these two.
NSS_STATISTIC_GROUPS_WA = [
    {"id": 2, "name": "Peak-Flow Statistics", "code": "PFS", "defType": "FS"},
    {"id": 4, "name": "Low-Flow Statistics", "code": "LFS", "defType": "FS"},
]

# GET /nssservices/regions/WA/Scenarios?statisticgroups=PFS&unitsystem=2 -- the
# "scenario template": one scenario per statistic group, each with one or more
# independently-calibrated regressionRegions carrying placeholder (-999.99) parameter
# values and the region's own valid [min, max] range. Trimmed to one region (of WA
# Peak-Flow's real four) for the mocked tests; the live test exercises all four.
NSS_SCENARIO_TEMPLATE_PFS = [
    {
        "statisticGroupID": 2,
        "statisticGroupName": "Peak-Flow Statistics",
        "regressionRegions": [
            {
                "id": 717,
                "name": "Peak_Region_1_2016_5118",
                "code": "GC1750",
                "description": "",
                "statusID": 4,
                "methodID": 1,
                "citationID": 150,
                "parameters": [
                    {
                        "id": 19532,
                        "name": "Drainage Area",
                        "description": "Area that drains to a point on a stream",
                        "code": "DRNAREA",
                        "unitType": {"id": 35, "unit": "square miles", "abbr": "mi^2"},
                        "value": -999.99,
                        "limits": {"max": 3310.0, "min": 0.25},
                    },
                    {
                        "id": 19533,
                        "name": "Mean Annual Precip PRISM 1981 2010",
                        "description": "Basin average mean annual precipitation for 1981 to 2010",
                        "code": "PRECPRIS10",
                        "unitType": {"id": 36, "unit": "inches", "abbr": "in"},
                        "value": -999.99,
                        "limits": {"max": 52.5, "min": 9.82},
                    },
                    {
                        "id": 19534,
                        "name": "Percent Area Under Canopy",
                        "description": "Percentage of drainage area covered by canopy",
                        "code": "CANOPY_PCT",
                        "unitType": {"id": 0, "unit": "percent", "abbr": "%"},
                        "value": -999.99,
                        "limits": {"max": 77.4, "min": 0.0},
                    },
                ],
            }
        ],
        "links": [
            {
                "rel": "Citations",
                "href": "https://streamstats.usgs.gov/nssservices/citations?regressionregions=717",
                "method": "GET",
            }
        ],
    }
]

# POST /nssservices/Scenarios/Estimate response -- the same shape as the template,
# with each regressionRegion's "results" now populated. Real values for Goat Creek
# (DRNAREA=412.0, PRECPRIS10=45.62, CANOPY_PCT=45.242), trimmed to 2 of the real 8
# statistics.
NSS_ESTIMATE_RESPONSE_PFS = [
    {
        "statisticGroupID": 2,
        "statisticGroupName": "Peak-Flow Statistics",
        "regressionRegions": [
            {
                "id": 717,
                "name": "Peak_Region_1_2016_5118",
                "code": "GC1750",
                "citationID": 150,
                "parameters": [
                    {"code": "DRNAREA", "value": 412.0},
                    {"code": "PRECPRIS10", "value": 45.62},
                    {"code": "CANOPY_PCT", "value": 45.242},
                ],
                "results": [
                    {
                        "equivalentYears": 0.0,
                        "id": 0,
                        "name": "50-percent AEP flood",
                        "code": "PK50AEP",
                        "description": "Maximum instantaneous flow with a 50% AEP",
                        "value": 4370.0,
                        "sep": 0.350768271099932,
                        "errors": [
                            {
                                "id": 0,
                                "name": "Average standard error of prediction",
                                "code": "ASEp",
                                "value": 95.0,
                            }
                        ],
                        "unit": {"id": 44, "unit": "cubic feet per second", "abbr": "ft^3/s"},
                        "equation": "3.846*DRNAREA^0.745*10^(0.032*PRECPRIS10)/10^(0.0078*CANOPY_PCT)",
                        "intervalBounds": {"lower": 1140.0, "upper": 16700.0},
                    },
                    {
                        "equivalentYears": 0.0,
                        "id": 0,
                        "name": "20-percent AEP flood",
                        "code": "PK20AEP",
                        "description": "Maximum instantaneous flow with a 20% AEP",
                        "value": 6050.0,
                        "sep": 0.28093883640234,
                        "errors": [
                            {
                                "id": 0,
                                "name": "Average standard error of prediction",
                                "code": "ASEp",
                                "value": 71.9,
                            }
                        ],
                        "unit": {"id": 44, "unit": "cubic feet per second", "abbr": "ft^3/s"},
                        "equation": "12.106*DRNAREA^0.713*10^(0.028*PRECPRIS10)/10^(0.0098*CANOPY_PCT)",
                        "intervalBounds": {"lower": 2060.0, "upper": 17700.0},
                    },
                ],
            }
        ],
    }
]

# GET /nssservices/citations?regressionregions=717 -- real citation for WA Peak-Flow.
NSS_CITATIONS_RESPONSE = [
    {
        "id": 150,
        "title": (
            "2016, Magnitude, frequency, and trends of floods at gaged and ungaged "
            "sites in Washington: U.S. Geological Survey Scientific Investigations "
            "Report 2016-5118, 70 p."
        ),
        "author": "Mastin, M.C., Konrad, C.P., Veilleux, A.G., and Tecca, A.E.,",
        "citationURL": "http://dx.doi.org/10.3133/sir20165118",
        "lastYearOfData": 2014,
    }
]

# The real, confirmed-live failure mode: POST /nssservices/Scenarios/Estimate with the
# scenario array wrapped in {"scenarioList": [...]} (matching the api-config metadata's
# literal parameter name) rather than posted bare -- a 500 with no other diagnostic.
NSS_ESTIMATE_500_BODY = (
    '{"code":500,"message":"An error occured while processing your request. '
    'See messages for more information.","content":"Internal Server Error Occured"}'
)
