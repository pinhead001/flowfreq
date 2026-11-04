"""Nebraska's verified equations (WRI 99-4032, Tables 2-8, appendix B tables B1 and B2).

The report has no worked example and no covariance. The cross-check is table B2's regression flow
at every gage, from table B1's characteristics, which are printed with few digits. So each flow is
compared with the band that the printed rounding of the characteristics, coefficients, constant and
flow allows: nothing here is a chosen tolerance.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from flowfreq.regression import evaluate
from flowfreq.regression.library import load_state
from flowfreq.ungaged import region_for_nss_code

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "flowfreq" / "data" / "regression"
NSS = {
    "HPS": "GC141",
    "HPC": "GC142",
    "N&W": "GC143",
    "NE": "GC144",
    "C&SC": "GC145",
    "East": "GC146",
    "UR": "GC147",
    "BB": "GC148",
}
AEPS = (0.5, 0.2, 0.1, 0.04, 0.02, 0.01, 0.005, 0.002)
GAGES = json.loads((ROOT / "fixtures" / "regression" / "nebraska_gages.json").read_text("utf-8"))

#: ASE and SEE columns of Tables 2-8, which have no schema field.
ASE = {
    "HPS": (0.003, 0.004, 0.005, 0.007, 0.007, 0.008, 0.009, 0.011),
    "HPC": (0.006, 0.008, 0.011, 0.014, 0.016, 0.019, 0.022, 0.026),
    "N&W": (0.032, 0.014, 0.014, 0.016, 0.018, 0.018, 0.020, 0.023),
    "NE": (0.007, 0.006, 0.006, 0.007, 0.008, 0.010, 0.009, 0.011),
    "C&SC": (0.016, 0.011, 0.012, 0.014, 0.016, 0.019, 0.021, 0.025),
    "East": (0.006, 0.004, 0.004, 0.005, 0.005, 0.006, 0.007, 0.008),
    "UR": (0.008, 0.008, 0.008, 0.010, 0.012, 0.014, 0.016, 0.018),
    "BB": (0.007, 0.004, 0.004, 0.004, 0.005, 0.006, 0.006, 0.008),
}
SEE = {
    "HPS": (0.174, 0.172, 0.176, 0.181, 0.189, 0.195, 0.203, 0.217),
    "HPC": (0.149, 0.175, 0.182, 0.195, 0.201, 0.210, 0.223, 0.244),
    "N&W": (0.424, 0.247, 0.222, 0.224, 0.236, 0.254, 0.259, 0.274),
    "NE": (0.191, 0.153, 0.147, 0.151, 0.157, 0.166, 0.176, 0.185),
    "C&SC": (0.269, 0.196, 0.187, 0.198, 0.212, 0.228, 0.245, 0.268),
    "East": (0.191, 0.126, 0.107, 0.104, 0.109, 0.116, 0.124, 0.136),
    "UR": (0.211, 0.192, 0.196, 0.211, 0.224, 0.239, 0.255, 0.276),
    "BB": (0.164, 0.079, 0.044, 0.041, 0.045, 0.052, 0.061, 0.074),
}
#: Table 5 at 1 percent: SEP .263 is below sqrt(ASE + AME) even at the edges of their rounding.
SEP_DISAGREES_WITH_ASE_PLUS_AME = {("C&SC", 0.01)}

#: Flows outside the band. Table B1 prints AWC 1.64 and PLP 41.61 at 06464500 (see its own test).
B1_MISPRINT = {("06464500", "N&W"): {0.5, 0.2, 0.1, 0.04, 0.02}}
#: Consistent across a gage's low-AEP rows with B2 having used other TTP and RR than B1 prints.
C_AND_SC_INPUTS = {
    ("06768050", "C&SC"): {0.5, 0.2, 0.1, 0.04, 0.02},
    ("06769000", "C&SC"): {0.5, 0.2, 0.1},
    ("06769500", "C&SC"): {0.5, 0.2, 0.1, 0.04, 0.02},
    ("06784700", "C&SC"): {0.2},
    ("06849600", "C&SC"): {0.5, 0.2, 0.1},
}
#: Isolated cells, each at a gage whose other flows reproduce. Printed as read from the page image.
UNEXPLAINED_CELLS = {
    ("06396490", "N&W"): {0.2},
    ("06447500", "HPC"): {0.2},
    ("06449500", "N&W"): {0.005, 0.002},
    ("06453500", "N&W"): {0.005, 0.002},
    ("06456400", "N&W"): {0.5, 0.002},
    ("06810060", "East"): {0.005},
    ("06880710", "BB"): {0.5},
    ("06880720", "BB"): {0.1},
    ("06883000", "BB"): {0.5},
}
PINNED = {**B1_MISPRINT, **C_AND_SC_INPUTS, **UNEXPLAINED_CELLS}
#: Printed characteristics just past a printed range, by the range's rounding, and not marked in B2.
ROUNDS_PAST_RANGE = {
    ("06464500", "N&W"),  # RR 4.17 for 4.2, and the PLP misprint
    ("06791500", "NE"),  # PLP 5.57 for 5.56
    ("06799080", "NE"),  # DF 0.007 for 0.01
    ("06800500", "HPS"),  # MCS 4.4 for 4.41
    ("06883600", "BB"),  # SD 0.138 for 0.14
}


def _ulp(token: str) -> float:
    """The unit of a printed number's last digit; these tables print integers to 3 figures."""
    text = token.replace(",", "")
    if "." in text:
        return 10.0 ** -(len(text) - text.index(".") - 1)
    return 10.0 ** min(len(text) - len(text.rstrip("0")), max(len(text) - 3, 0))


def _half_width(eq, printed: dict[str, str]) -> float:
    """Half-width, in log10 units, of the estimate from the printed rounding of its inputs."""
    constant = 10.0**eq.intercept
    width = 0.5 * 10.0 ** (math.floor(math.log10(constant)) - 2) / constant / math.log(10.0)
    for v, b in zip(eq.variables, eq.coefficients):
        x = float(printed[v.code].replace(",", "")) + v.offset
        width += abs(b) * (0.5 * _ulp(printed[v.code]) / x) / math.log(10.0)
        width += (0.005 if abs(b) >= 1.0 else 0.0005) * abs(math.log10(x))
    return width


def _inside(eq, printed: dict[str, str], flow: int) -> bool:
    inputs = {k: float(t.replace(",", "")) for k, t in printed.items()}
    est = evaluate(eq, inputs, allow_extrapolation=True)
    allowed = _half_width(eq, printed) + 0.5 * _ulp(str(flow)) / flow / math.log(10.0)
    return abs(est.log_flow - math.log10(flow)) <= allowed


@pytest.fixture(scope="module")
def lib():
    return load_state("NE")


def test_nebraska_is_verified(lib):
    assert lib.status == "verified"
    assert set(lib.regions) == set(NSS)
    assert len(lib.equations) == 64
    for region, code in NSS.items():
        assert region_for_nss_code(lib, code) == region
        assert set(lib.aeps(region)) == set(AEPS)


def test_every_regression_flow_reproduces_within_its_printed_rounding_but_for_the_pins(lib):
    assert len(GAGES["rows"]) == 323
    cells = outside = 0
    found: dict[tuple[str, str], set[float]] = {}
    for row in GAGES["rows"]:
        for aep, flow in zip(AEPS, row["flows"]):
            cells += 1
            if not _inside(lib.equation(row["region"], aep), row["characteristics"], flow):
                outside += 1
                found.setdefault((row["station"], row["region"]), set()).add(aep)
    assert cells == 2584
    assert outside == 34
    assert found == PINNED


def test_table_b1_misprints_permeability_at_keya_paha_river_at_wewela(lib):
    """AWC 1.64 and PLP 41.61 are impossible, and five PLP equations imply about 2.13."""
    row = next(r for r in GAGES["rows"] if r["station"] == "06464500" and r["region"] == "N&W")
    assert row["characteristics"]["LSTPERM"] == "41.61"
    fixed = {**row["characteristics"], "LSTPERM": "2.13"}
    for aep, flow in zip(AEPS, row["flows"]):
        assert _inside(lib.equation("N&W", aep), fixed, flow)


def test_evaluate_gives_a_published_flow(lib):
    """Big Blue River at Barneston (06882000), to B2's three figures at 1 percent AEP."""
    row = next(r for r in GAGES["rows"] if r["station"] == "06882000")
    inputs = {k: float(t.replace(",", "")) for k, t in row["characteristics"].items()}
    est = evaluate(lib.equation("BB", 0.01), inputs)
    assert est.flow_cfs == pytest.approx(58100, rel=0.01)
    assert row["flows"][5] == 58100


def test_evaluate_refuses_a_flow_outside_the_published_range(lib):
    inputs = {"DRNAREA": 7000.0, "AVMXSS": 6.5, "BSHAPE": 2.0, "STRDEN": 0.6}
    with pytest.raises(ValueError):
        evaluate(lib.equation("BB", 0.01), inputs)


def test_gages_b2_marks_out_of_range_are_outside_the_tables(lib):
    assert len(GAGES["outside_range"]) == 12
    for row in GAGES["outside_range"]:
        eqs = [lib.equation(row["region"], aep) for aep in AEPS]
        limits = {v.code: (v.minimum, v.maximum) for e in eqs for v in e.variables}
        values = {k: float(t.replace(",", "")) for k, t in row["characteristics"].items()}
        assert any(not lo <= values[k] <= hi for k, (lo, hi) in limits.items())


def test_other_gages_are_inside_the_tables_but_for_rounding(lib):
    past = set()
    for row in GAGES["rows"]:
        eqs = [lib.equation(row["region"], aep) for aep in AEPS]
        limits = {v.code: (v.minimum, v.maximum) for e in eqs for v in e.variables}
        values = {k: float(t.replace(",", "")) for k, t in row["characteristics"].items()}
        if any(not lo <= values[k] <= hi for k, (lo, hi) in limits.items()):
            past.add((row["station"], row["region"]))
    assert past == ROUNDS_PAST_RANGE


def test_station_counts_reproduce_from_table_b2(lib):
    """Regression gages: B2's rows with flows, less those it marks as new network stations."""
    for region in NSS:
        used = {
            r["station"] for r in GAGES["rows"] if r["region"] == region and not r["new_station"]
        }
        assert len(used) == lib.equation(region, 0.01).n_sites, region


def test_error_statistics_are_consistent_with_the_reports_definitions(lib):
    """SEP = sqrt(ASE + AME) and SEE = sqrt(AME), to the rounding of three decimals."""
    bad = set()
    for region in NSS:
        for i, aep in enumerate(AEPS):
            eq = lib.equation(region, aep)
            sep, ame = eq.sep_log, eq.model_error_variance
            if abs(sep**2 - (ASE[region][i] + ame)) > 0.001 + 2 * sep * 0.0005 + 1e-12:
                bad.add((region, aep))
            see = SEE[region][i]
            assert abs(see**2 - ame) <= 0.0005 + 2 * see * 0.0005 + 1e-12, (region, aep)
    assert bad == SEP_DISAGREES_WITH_ASE_PLUS_AME


def test_every_equation_equals_the_nss_snapshot(lib):
    snap = json.loads((DATA / "nss_snapshot" / "NE.json").read_text(encoding="utf-8"))
    by_key = {(s["region_code"], s["aep"]): s for s in snap["equations"]}
    assert {(NSS[e.region_code], e.aep) for e in lib.equations} == set(by_key)
    for e in lib.equations:
        s = by_key[(NSS[e.region_code], e.aep)]
        assert e.intercept == pytest.approx(s["intercept"], abs=1e-12)
        mine = {
            v.code: (b, v.minimum, v.maximum, v.offset) for v, b in zip(e.variables, e.coefficients)
        }
        theirs = {
            v["code"]: (b, v["minimum"], v["maximum"], v.get("offset", 0.0))
            for v, b in zip(s["variables"], s["coefficients"])
        }
        assert mine == theirs, (e.region_code, e.aep)


def test_high_permeability_keeps_both_analyses_with_no_preference(lib):
    assert lib.equation("HPS", 0.01).n_sites == 49
    assert lib.equation("HPC", 0.01).n_sites == 23
