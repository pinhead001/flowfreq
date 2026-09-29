"""StreamStats API payloads for delineation/characteristics tests.

Values are for the Methow River basin, WA (the design doc's own verification points,
``docs/STREAMSTATS_MODULE_DESIGN.md`` appendix). The snap shape is confirmed live
(2026-09-11); basin characteristics per the design doc (2026-09-09/10).

The ``delineate/sshydro`` responses are live captures (2026-09-27) stored as JSON under
``tests/fixtures/streamstats/`` (see the comment on each constant for how it was
trimmed). They carry the watershed polygon (design doc S10) and, for the off-network
point, the real nesting of ``WarningMsg``.
"""

import json
from pathlib import Path
from typing import Any

_CAPTURES = Path(__file__).parent / "streamstats"


def load_capture(name: str) -> Any:
    """A fresh copy of a committed live capture, so a test may mutate it freely."""
    return json.loads((_CAPTURES / name).read_text(encoding="utf-8"))


# A successful pourpoint snap: the point lies close to the stream network. `output` is
# a GeoJSON Point -- coordinates are [lon, lat], confirmed live 2026-09-11.
SNAP_GOOD = {
    "region": "WA",
    "input": {"type": "Point", "coordinates": [-120.37893, 48.57426]},
    "output": {"type": "Point", "coordinates": [-120.3789336506848, 48.57425501781417]},
    "couldSnap": True,
}

# Methow River near Pateros (USGS 12449950's site), live 2026-09-27.
SNAP_METHOW_PATEROS = {
    "region": "WA",
    "input": {"type": "Point", "coordinates": [-119.9837, 48.0776]},
    "output": {"type": "Point", "coordinates": [-119.98382843100944, 48.0772686409578]},
    "couldSnap": True,
}

# The point is off the flowline -- the service will not snap it. The whole point of
# this fixture is that FR-1 forbids proceeding to delineation from here.
SNAP_UNSNAPPABLE = {
    "couldSnap": False,
    "output": {},
}

# The ss-delineate 'sshydro' response for Goat Creek: the bcrequest body to POST to
# ss-hydro, which also carries the globalwatershed polygon. Live capture 2026-09-27,
# polygon trimmed from 4590 to 332 vertices (every vertex within 2 km of the pour point
# kept, every 20th elsewhere): area 411.72 mi^2 against DRNAREA 412.0.
DELINEATE_SSHYDRO_GOOD = load_capture("sshydro_WA_goat_creek_trimmed.json")

# Methow near Pateros: 1793 mi^2, and its polygon has a hole (an interior ring). Live
# capture 2026-09-27, exterior ring trimmed the same way from 10309 to 358 vertices,
# hole kept whole: area 1792.80 mi^2 against DRNAREA 1793.0.
DELINEATE_SSHYDRO_METHOW_PATEROS = load_capture("sshydro_WA_methow_pateros_trimmed.json")

# ss-hydro's answer for the Methow at Pateros, trimmed to DRNAREA (live 2026-09-27).
HYDRO_CHARACTERISTICS_METHOW_PATEROS = [
    {
        "code": "DRNAREA",
        "name": "Drainage Area",
        "description": "Area that drains to a point on a stream",
        "value": 1793.0,
        "unit": "square miles",
        "msg": "Local AreaOp successful",
    }
]

# The design doc's own failure mode (S4), captured live 2026-09-27 and verbatim: HTTP
# 200, a 13-vertex hillslope sliver, and the WarningMsg in both features' properties.
# This fixture -- and the test asserting it raises -- is the whole point of this module.
DELINEATE_SSHYDRO_WITH_WARNING = load_capture("sshydro_WA_offnetwork.json")

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

# --- NSS region selection by watershed polygon (addendum S5). Live captures,
# --- 2026-09-27, under tests/fixtures/streamstats/.

# GET /nssservices/regions/WA/Scenarios?statisticgroups=PFS,LFS&unitsystem=2 -- the
# full WA template: four Peak-Flow regions (GC1750-GC1753), four Low-Flow regions.
NSS_SCENARIOS_WA_PFS_LFS = load_capture("nss_scenarios_WA_PFS_LFS.json")

# POST /nssservices/regions/WA/regressionregions/bylocation with Goat Creek's polygon:
# one region, Peak Region 2 (GC1751), 100% -- not GC1750, which the older fixtures use.
NSS_BYLOCATION_WA_GOAT_CREEK = load_capture("nss_bylocation_WA_goat_creek.json")

# POST /nssservices/Scenarios/Estimate for Goat Creek's real characteristics: GC1751
# only (50% AEP 3,290 cfs), and all four PFS regions (4,370 / 3,290 / 5,280 / 7,970).
NSS_ESTIMATE_WA_GOAT_CREEK_GC1751 = load_capture("nss_estimate_WA_goat_creek_GC1751.json")
NSS_ESTIMATE_WA_GOAT_CREEK_ALL_PFS = load_capture("nss_estimate_WA_goat_creek_all_PFS.json")

# Ogeechee River near Louisville, GA (USGS 02200500, DRNAREA 805): the one live basin
# found that spans several regions of a state's own study. The sshydro capture is a
# MultiPolygon (three one-cell islands), exterior ring trimmed from 22119 to 722
# vertices; region-scoped bylocation puts it 44/32/24% in SIR 2014-5030 regions 1/3/4
# (ss-hydro's own PCTREG1/3/4 are 44.05/31.87/24.08), and 100% in GC1934.
DELINEATE_SSHYDRO_GA_OGEECHEE = load_capture("sshydro_GA_ogeechee_louisville_trimmed.json")
NSS_BYLOCATION_GA_OGEECHEE = load_capture("nss_bylocation_GA_ogeechee_louisville.json")

# GET /nssservices/regions/GA/Scenarios?statisticgroups=PFS&unitsystem=2&
# regressionregions=GC1572,GC1573 -- two GA rural (<1 mi^2) regions, used to exercise
# NSS's area averaging.
NSS_SCENARIOS_GA_RURAL_UNDER_1 = load_capture("nss_scenarios_GA_PFS_rural_under_1sqmi.json")

# POST /nssservices/Scenarios/Estimate with those two regions carrying percentWeight
# 60 and 40 (chosen for the test, not from a basin) and DRNAREA=0.5, LC06IMP=5,
# LC06DEV=20: NSS appends a third region, code "areaave", name "Area-Averaged", whose
# values are the weighted mean (0.6*116 + 0.4*33.3 = 82.92 for PK50AEP).
NSS_ESTIMATE_GA_AREA_AVERAGED = load_capture("nss_estimate_GA_area_averaged.json")

# --- Low-Flow Statistics (LFS), live 2026-09-28 (addendum S6).

# ss-hydro's answer when asked (BCs=DRNAREA;ELEV1000) for a characteristic the region
# does not compute: HTTP 200 and a -999 sentinel. code/value/unit/msg are as returned
# live for Skookumchuck River near Vail, WA; the ELEV1000 entry's name/description were
# not recorded and are filled in here. WA computes 11 characteristics; ELEV1000 is not
# one of them.
HYDRO_CHARACTERISTICS_WITH_UNAVAILABLE = [
    {
        "name": "Drainage Area",
        "description": "Area that drains to a point on a stream",
        "code": "DRNAREA",
        "unit": "square miles",
        "value": 39.9,
        "msg": "Local AreaOp successful",
    },
    {
        "name": "ELEV1000",
        "description": "",
        "code": "ELEV1000",
        "unit": "",
        "value": -999.0,
        "msg": "Basin Characteristic not found in database",
    },
]

# GET /nssservices/regions/WA/Scenarios?statisticgroups=LFS&unitsystem=2 -- WA's four
# Low-Flow regions: GC1434 (Nooksack, DRNAREA+ELEV1000), GC1556 (DRNAREA),
# GC1557 (DRNAREA+PRECIP), GC1558 (DRNAREA+PRECIP+TAU_ANN_G).
NSS_SCENARIOS_WA_LFS = load_capture("nss_scenarios_WA_LFS.json")

# POST /nssservices/Scenarios/Estimate for GC1556 and GC1557 with Skookumchuck River
# near Vail's live characteristics (DRNAREA 39.9, PRECIP 71.07): M7D10Y 16.2 and 12.0
# cfs, standard error reported under code "SE" (133 and 114), not "ASEp".
NSS_ESTIMATE_WA_LFS_SKOOKUMCHUCK = load_capture("nss_estimate_WA_LFS_skookumchuck_vail.json")
