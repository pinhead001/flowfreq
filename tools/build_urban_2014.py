#!/usr/bin/env python
"""Store SIR 2014-5030's urban equations (Georgia, North and South Carolina) at their full range.

The state files already hold the report's equations (Feaster, Gotvald and Weaver, 2014, ver.
1.1) for small rural basins, to 1 mi2, as regions HR1, HR3, HR4 and (Georgia) HR5. The report's
Table 7 takes the same equations to 3 mi2 in hydrologic region 1, 459 in region 3, 53.5 in
region 4 and 10 in region 5, and has another column for region 1 from 3 to 436 mi2. This tool
stores all of them as sets used only when named (``named_regions``): the owner's decision of
2026-11-24. The first four are copies of the stored equations with the report's drainage-area
range; the fifth is transcribed here from Tables 7, 9, 10 and 11.

Running it again changes nothing. Run it before or after ``build_urban_2026.py``.

    python tools/build_urban_2014.py
"""

from __future__ import annotations

import copy
import json
import math
from pathlib import Path
from typing import Any, Dict, List

DATA = Path(__file__).resolve().parent.parent / "flowfreq" / "data" / "regression"
STATES = ("GA", "NC", "SC")
PREFIX = "Urban 2014 "
LARGE = PREFIX + "HR1 large"
AEPS = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)

#: stored region -> (new region, its name, upper drainage area in mi2 of Table 7).
#: The lower limits stay the stored ones but for HR3, whose Table 7 starts at 0.22 mi2 (0.14
#: is Table 11's smallest gage of the region, rural).
LOWER = {"HR3": 0.22}
WIDENED = {
    "HR1": (
        PREFIX + "HR1 small",
        "Hydrologic region 1, Piedmont and Ridge and Valley, urban, 0.10 to 3 mi2 (SIR 2014-5030)",
        3.0,
    ),
    "HR3": (PREFIX + "HR3", "Hydrologic region 3, Sand Hills, urban (SIR 2014-5030)", 459.0),
    "HR4": (PREFIX + "HR4", "Hydrologic region 4, Coastal Plain, urban (SIR 2014-5030)", 53.5),
    "HR5": (
        PREFIX + "HR5",
        "Hydrologic region 5, Lower Tifton Upland, urban (SIR 2014-5030)",
        10.0,
    ),
}
LARGE_NAME = (
    "Hydrologic region 1, Piedmont and Ridge and Valley, urban, 3 to 436 mi2 (SIR 2014-5030)"
)

# Hydrologic region 1, 3 to 436 mi2: Table 7's constant, exponent of DRNAREA and factor of
# IMPNLCD06; Table 9's average standard error of prediction (percent) and AVP; Table 10's
# model error variance and the covariance's upper triangle by rows.
LARGE_TABLE = {
    0.5: (198.0, 0.5735, 0.0101, 34.4, 0.021, 0.020)
    + (1.51e-03, -6.19e-04, -1.79e-05, 3.58e-04, 3.53e-06, 1.31e-06),
    0.2: (359.0, 0.5605, 0.0074, 31.4, 0.018, 0.017)
    + (1.49e-03, -5.90e-04, -1.69e-05, 3.33e-04, 3.39e-06, 1.27e-06),
    0.1: (484.0, 0.5539, 0.0060, 30.7, 0.017, 0.016)
    + (1.59e-03, -6.19e-04, -1.78e-05, 3.45e-04, 3.54e-06, 1.35e-06),
    0.04: (657.0, 0.5470, 0.0046, 31.4, 0.018, 0.017)
    + (1.81e-03, -7.03e-04, -2.02e-05, 3.89e-04, 4.02e-06, 1.56e-06),
    0.02: (794.0, 0.5428, 0.0037, 32.4, 0.019, 0.018)
    + (2.02e-03, -7.85e-04, -2.26e-05, 4.33e-04, 4.49e-06, 1.75e-06),
    0.01: (941.0, 0.5386, 0.0028, 34.2, 0.021, 0.020)
    + (2.28e-03, -8.90e-04, -2.56e-05, 4.91e-04, 5.10e-06, 1.99e-06),
    0.005: (1096.0, 0.5351, 0.0021, 35.8, 0.023, 0.022)
    + (2.53e-03, -9.94e-04, -2.86e-05, 5.48e-04, 5.72e-06, 2.23e-06),
    0.002: (1319.0, 0.5305, 0.0011, 38.7, 0.026, 0.025)
    + (2.93e-03, -1.16e-03, -3.33e-05, 6.40e-04, 6.71e-06, 2.60e-06),
}
#: 175 rural and 60 urban streamgages (Table 6), and Table 11's range.
LARGE_SITES = 235
LARGE_AREA = (3.0, 436.0)
IMPERVIOUS = (0.0, 47.9)

NOTE = (
    " Urban 2014: SIR 2014-5030 again (Table 7), at its whole drainage area range, by the "
    "owner's decision of 2026-11-24, as sets used only when named (named_regions; in no "
    "region_preference tier and not keyed to the urban regions NSS serves): 'Urban 2014 HR1 "
    "small' (0.10 to 3 mi2), 'Urban 2014 HR1 large' (3 to 436 mi2), 'Urban 2014 HR3' (0.22 to 459), "
    "'Urban 2014 HR4' (to 53.5){hr5}. All but 'HR1 large' are the equations stored as HR1, HR3"
    ", HR4{hr5_short}, whose numbers they copy, by tools/build_urban_2014.py. 'HR1 large' is "
    "Tables 7, 9, 10 and 11 (235 streamgages, 60 of them urban): Table 5's regression flows "
    "reproduce within 0.8 percent at the 61 gages it prints over 3 mi2, the coefficients and "
    "matrices are the report's spreadsheet tool's to rounding, and MEV + mean(xUx) over Table "
    "8's 235 gages returns Table 9's AVP. Reedy River example (30.0 mi2, 18 percent, p. 46): "
    "6,600 ft3/s reproduces. The xUx the report prints for its HR1 example (0.001447) is what "
    "Table 10 gives, to rounding, with the impervious area entered as 0.30 where the example "
    "says 30.0, as its spreadsheet tool enters it (percent over 100); the coefficients and "
    "matrices are for percent (the tests of significance on p. 25 follow only so), and "
    "0.00196 is used here. A named set is used for the whole basin: the report's average by area over "
    "hydrologic regions (p. 47) is not made."
)
#: What the note on the small-basin sets said of the rest of the report, and says now.
SCOPE = (
    "NSS serves those as separate urban regions, which are out of scope with HR1's equations "
    "over 3 mi2.",
    "NSS serves those, and HR1's equations over 3 mi2, as separate urban regions; here they are "
    "the Urban 2014 sets.",
)
NOTE_HR5 = (", 'Urban 2014 HR5' (0.2 to 10)", " and HR5")

#: What every one of the sets carries; ``{range}`` is the set's own sentence.
CAUTION = (
    "SIR 2014-5030 (Feaster, Gotvald and Weaver, 2014, Table 7) is for urban and small rural "
    "streams, with peaks through water year 2011 at the latest. At 1 mi2 or more and under 10 "
    "percent impervious the report says to use the rural equations (pp. 8, 34; it names the "
    "2009 rural reports, and the rural set stored here is HR1-5, SIR 2023-5006). {range} "
    "Named, this set is used for the whole basin: the report averages a basin in two hydrologic "
    "regions by area (p. 47), and has no equations for urban streams in the Blue Ridge (region "
    "2). Nothing here checks which region the basin is in. Land cover is from NLCD 2006 (a "
    "later NLCD by the same methods may give the impervious area), not NLCD 1992. Not for streams substantially affected by "
    "regulation, channelization or levees, nor by storm surge or tides."
)
RANGES = {
    PREFIX
    + "HR1 small": "Its 62 gages reach 3 mi2, 21 of them urban, the least 13.2 percent impervious.",
    LARGE: "It starts at 3 mi2. The largest urban basin among its 235 gages (60 urban) is 378 mi2.",
    PREFIX
    + "HR3": (
        "Only 8 of its 34 gages are urban, at 67 to 98.5 percent developed, and the largest of "
        "those drains 5.67 mi2; its standard error of prediction is 42 to 73 percent."
    ),
    PREFIX
    + "HR4": (
        "Of its 157 gages 27 are urban, 16 of them in New Jersey, and the largest urban basin "
        "is 43.9 mi2."
    ),
    PREFIX
    + "HR5": (
        "It is Gotvald and Knaak (2011), not refitted, and has no urban term; the basins it was "
        "fitted to were 0 to 36 percent impervious."
    ),
}


def _log_error(percent: float) -> float:
    """A standard error in percent as one in log10 units."""
    return math.sqrt(math.log(1.0 + (percent / 100.0) ** 2)) / math.log(10.0)


def large_equations(template: Dict[str, Any]) -> List[Dict[str, Any]]:
    """The 3 to 436 mi2 equations of hydrologic region 1, cited as ``template`` is."""
    equations = []
    for aep, row in LARGE_TABLE.items():
        constant, exponent, factor, percent, avp, mev, aa, ab, ac, bb, bc, cc = row
        variables = [
            ("DRNAREA", "log10", "mi2", *LARGE_AREA),
            ("LC06IMP", "identity", "percent", *IMPERVIOUS),
        ]
        keys = ("code", "transform", "units", "minimum", "maximum")
        equations.append(
            {
                "region_code": LARGE,
                "region_name": LARGE_NAME,
                "aep": aep,
                "intercept": math.log10(constant),
                "variables": [dict(zip(keys, v)) for v in variables],
                "coefficients": [exponent, factor],
                "citation": copy.deepcopy(template["citation"]),
                "sep_log": round(_log_error(percent), 6),
                "avp": avp,
                "model_error_variance": mev,
                "covariance": [[aa, ab, ac], [ab, bb, bc], [ac, bc, cc]],
                "n_sites": LARGE_SITES,
            }
        )
    return equations


def widened(equation: Dict[str, Any]) -> Dict[str, Any]:
    """A stored small-basin equation with the report's drainage-area range."""
    region, name, maximum = WIDENED[equation["region_code"]]
    copied = copy.deepcopy(equation)
    copied["region_code"], copied["region_name"] = region, name
    area = copied["variables"][0]
    if area["code"] != "DRNAREA" or area["maximum"] != 1.0:
        raise ValueError(f"unexpected drainage area in {equation['region_code']}: {area}")
    area["maximum"] = maximum
    area["minimum"] = LOWER.get(equation["region_code"], area["minimum"])
    return copied


def store(state: str) -> None:
    """Rewrite the state's file with the sets."""
    path = DATA / f"{state}.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    stored = [e for e in data["equations"] if not e["region_code"].startswith(PREFIX)]
    small = [e for e in stored if e["region_code"] in WIDENED]
    added = [widened(e) for e in small]
    template = next(e for e in small if e["region_code"] == "HR1")
    added += large_equations(template)
    added.sort(key=lambda e: e["region_code"])  # stable: AEPs stay in their order
    data["equations"] = stored + added
    regions = sorted({e["region_code"] for e in added})

    kept = {k: v for k, v in data["region_cautions"].items() if not k.startswith(PREFIX)}
    data["region_cautions"] = {
        **kept,
        **{r: CAUTION.replace("{range}", RANGES[r]) for r in regions},
    }
    named = [r for r in data["named_regions"] if not r.startswith(PREFIX)]
    data["named_regions"] = named + regions

    # Before the 2026 set's note, which build_urban_2026.py replaces to the end.
    head, marker, tail = data["notes"].partition(" Urban 2026: ")
    head = head.partition(" Urban 2014: ")[0].replace(*SCOPE)
    hr5, short = NOTE_HR5 if "HR5" in {e["region_code"] for e in small} else ("", "")
    data["notes"] = head + NOTE.format(hr5=hr5, hr5_short=short) + marker + tail

    # The state files are indent-1 ASCII JSON with CRLF line ends.
    text = json.dumps(data, indent=1) + "\n"
    path.write_bytes(text.replace("\n", "\r\n").encode("ascii"))


def main() -> int:
    for state in STATES:
        store(state)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
