"""Colorado SIR 2025-5047 peak-flow equations, stored as an opt-in second set.

Kohn, Mast and Gross, 2026, USGS SIR 2025-5047, Table 1 (pp. 20-22): GLS peak-flow
equations for four elevation regions of central and western Colorado. NSS still serves
the 2009/2016/2006 equations in ``CO.json``, so this USGS set lives in
``provisional/CO.json`` and is returned only with ``allow_provisional=True``.

Every literal number below is typed from the report's Table 1 or copied from a dated
NSS staging response, never from the JSON itself: the tests check the file against the
source, not against itself. The report prints no worked example, so the evaluation
check uses NSS's staging draft of this report (test.streamstats.usgs.gov, citation 200,
regression regions 995-998, status "Review", fetched 2026-10-13).
"""

from __future__ import annotations

import math

import pytest

from flowfreq.regression import EquationsUnavailable, OutOfRangeError, evaluate, load_state

AEP_OF_RP = {2: 0.5, 5: 0.2, 10: 0.1, 25: 0.04, 50: 0.02, 100: 0.01, 200: 0.005, 500: 0.002}

# Table 1: (return period, b0, AREA exponent, ((variable, coefficient), ...), N,
# pseudo-R2 %, SEP %, SME %). Variables use the StreamStats codes printed in the
# table's "approximate range" row.
TABLE1 = {
    "Alpine": [
        (2, 0.244, 0.807, (("BSLDEM10M", 0.00898), ("APR1_SWE", 0.0537)), 105, 89.6, 39, 38),
        (5, 0.526, 0.786, (("BSLDEM10M", 0.00841), ("APR1_SWE", 0.0472)), 105, 90.8, 34, 33),
        (10, 0.664, 0.775, (("BSLDEM10M", 0.00819), ("APR1_SWE", 0.0440)), 105, 90.7, 33, 32),
        (25, 0.803, 0.763, (("BSLDEM10M", 0.00801), ("APR1_SWE", 0.0410)), 105, 89.7, 34, 33),
        (50, 0.888, 0.756, (("BSLDEM10M", 0.00794), ("APR1_SWE", 0.0392)), 105, 88.8, 35, 34),
        (100, 0.961, 0.749, (("BSLDEM10M", 0.00789), ("APR1_SWE", 0.0377)), 105, 87.5, 37, 36),
        (200, 1.03, 0.742, (("BSLDEM10M", 0.00787), ("APR1_SWE", 0.0365)), 105, 85.9, 39, 38),
        (500, 1.10, 0.735, (("BSLDEM10M", 0.00787), ("APR1_SWE", 0.0350)), 105, 83.8, 42, 40),
    ],
    "SubAlpine": [
        (
            2,
            -6.29,
            0.895,
            (("LONG_CENT", -0.0604), ("I24H100Y", 0.115), ("APR1_SWE", 0.0363)),
            104,
            87.1,
            54,
            52,
        ),
        (
            5,
            -6.05,
            0.874,
            (("LONG_CENT", -0.0606), ("I24H100Y", 0.119), ("APR1_SWE", 0.0303)),
            104,
            88.4,
            49,
            47,
        ),
        (
            10,
            -5.94,
            0.864,
            (("LONG_CENT", -0.0607), ("I24H100Y", 0.124), ("APR1_SWE", 0.0273)),
            104,
            88.9,
            47,
            46,
        ),
        (
            25,
            -5.79,
            0.854,
            (("LONG_CENT", -0.0604), ("I24H100Y", 0.131), ("APR1_SWE", 0.0241)),
            104,
            88.5,
            48,
            46,
        ),
        (
            50,
            -5.67,
            0.847,
            (("LONG_CENT", -0.0600), ("I24H100Y", 0.137), ("APR1_SWE", 0.0222)),
            104,
            88.0,
            49,
            47,
        ),
        (
            100,
            -5.56,
            0.841,
            (("LONG_CENT", -0.0595), ("I24H100Y", 0.142), ("APR1_SWE", 0.0205)),
            104,
            87.4,
            50,
            48,
        ),
        (
            200,
            -5.43,
            0.835,
            (("LONG_CENT", -0.0588), ("I24H100Y", 0.148), ("APR1_SWE", 0.0189)),
            104,
            86.6,
            52,
            50,
        ),
        (
            500,
            -5.27,
            0.829,
            (("LONG_CENT", -0.0578), ("I24H100Y", 0.154), ("APR1_SWE", 0.0171)),
            104,
            85.5,
            54,
            52,
        ),
    ],
    "MidElevation": [
        (2, -6.73, 0.919, (("LONG_CENT", -0.0594), ("PRECIP", 0.0418)), 104, 87.3, 65, 63),
        (5, -7.21, 0.893, (("LONG_CENT", -0.0687), ("PRECIP", 0.0332)), 104, 88.2, 59, 57),
        (10, -7.33, 0.878, (("LONG_CENT", -0.0723), ("PRECIP", 0.0286)), 104, 88.1, 58, 56),
        (25, -7.36, 0.860, (("LONG_CENT", -0.0753), ("PRECIP", 0.0239)), 104, 87.1, 60, 58),
        (50, -7.36, 0.848, (("LONG_CENT", -0.0769), ("PRECIP", 0.0208)), 104, 85.7, 63, 61),
        (100, -7.33, 0.837, (("LONG_CENT", -0.0782), ("PRECIP", 0.0182)), 104, 84.4, 66, 64),
        (200, -7.29, 0.826, (("LONG_CENT", -0.0791), ("PRECIP", 0.0158)), 104, 83.0, 69, 67),
        (500, -7.23, 0.814, (("LONG_CENT", -0.0802), ("PRECIP", 0.0129)), 104, 80.9, 74, 72),
    ],
    "Plateau": [
        (2, 2.37, 0.646, (("ELEV", -0.000364), ("I24H100Y", 0.478)), 105, 73.9, 112, 108),
        (5, 2.77, 0.573, (("ELEV", -0.000368), ("I24H100Y", 0.514)), 105, 75.8, 89, 86),
        (10, 2.98, 0.539, (("ELEV", -0.000378), ("I24H100Y", 0.548)), 105, 73.4, 90, 86),
        (25, 3.21, 0.506, (("ELEV", -0.000397), ("I24H100Y", 0.598)), 105, 67.8, 99, 94),
        (50, 3.36, 0.486, (("ELEV", -0.000412), ("I24H100Y", 0.635)), 105, 62.8, 109, 104),
        (100, 3.50, 0.469, (("ELEV", -0.000427), ("I24H100Y", 0.669)), 105, 58.2, 120, 114),
        (200, 3.63, 0.453, (("ELEV", -0.000441), ("I24H100Y", 0.702)), 105, 53.9, 133, 126),
        (500, 3.78, 0.435, (("ELEV", -0.000459), ("I24H100Y", 0.742)), 105, 48.4, 153, 144),
    ],
}

# Table 1 headnote: one pooled range for all four regions, in published units.
RANGES = {
    "DRNAREA": (0.22, 22629.0, "log10"),
    "LONG_CENT": (-110.219, -105.093, "identity"),
    "BSLDEM10M": (3.37, 69.4, "identity"),
    "ELEV": (4808.0, 11954.0, "identity"),
    "APR1_SWE": (0.0, 33.8, "identity"),
    "I24H100Y": (2.19, 5.21, "identity"),
    "PRECIP": (7.71, 61.6, "identity"),
}

# NSS staging codes for the report's variables (same quantity, different name).
NSS_CODE = {
    "BSLDEM10M": "ALTBSLDEM",
    "ELEV": "ALTELEV",
    "I24H100Y": "ALI24H100Y",
    "PRECIP": "PRECPRIS10",
}

# NSS staging Estimate, 2026-10-13: inputs (NSS codes) and returned values, 2 to 500 years.
NSS_STAGING = {
    "Alpine": [
        (
            {"DRNAREA": 70.55763601482124, "ALTBSLDEM": 15.291486154141461, "APR1_SWE": 16.9},
            [603, 804, 924, 1070, 1170, 1270, 1380, 1480],
        ),
        (
            {"DRNAREA": 5.0, "ALTBSLDEM": 40.0, "APR1_SWE": 25.0},
            [323, 391, 430, 481, 518, 553, 597, 636],
        ),
    ],
    "SubAlpine": [
        (
            {
                "DRNAREA": 70.55763601482124,
                "LONG_CENT": -107.6559929093485,
                "ALI24H100Y": 3.3812462281535303,
                "APR1_SWE": 16.9,
            },
            [740, 1010, 1180, 1380, 1560, 1690, 1830, 1970],
        ),
        (
            {"DRNAREA": 350.0, "LONG_CENT": -106.2, "ALI24H100Y": 2.6, "APR1_SWE": 8.0},
            [979, 1450, 1760, 2140, 2460, 2700, 2980, 3270],
        ),
    ],
    "MidElevation": [
        (
            {
                "DRNAREA": 70.55763601482124,
                "LONG_CENT": -107.6559929093485,
                "PRECPRIS10": 21.793026407546062,
            },
            [188, 363, 501, 720, 870, 1080, 1250, 1550],
        ),
        (
            {"DRNAREA": 1200.0, "LONG_CENT": -108.9, "PRECPRIS10": 30.0},
            [6640, 10400, 12700, 16000, 17800, 20400, 21900, 25000],
        ),
    ],
    "Plateau": [
        (
            {"DRNAREA": 70.55763601482124, "ALTELEV": 7000.0, "ALI24H100Y": 3.3812462281535303},
            [429, 980, 1530, 2450, 3320, 4370, 5680, 7590],
        ),
        (
            {"DRNAREA": 15.0, "ALTELEV": 6200.0, "ALI24H100Y": 4.4},
            [947, 2650, 4800, 9440, 14800, 22300, 32900, 51400],
        ),
    ],
}


def _pct_to_log(p: float) -> float:
    return math.sqrt(math.log(1.0 + (p / 100.0) ** 2)) / math.log(10.0)


def _from_nss(inputs: dict) -> dict:
    back = {v: k for k, v in NSS_CODE.items()}
    return {back.get(k, k): v for k, v in inputs.items()}


@pytest.fixture(scope="module")
def lib():
    return load_state("CO", allow_provisional=True)


def test_default_lookup_still_returns_the_nss_served_set():
    published = load_state("CO")
    assert published.status == "verified" and not published.provisional
    assert "Mountain" in published.regions and "Alpine" not in published.regions
    assert len(published.equations) == 62


def test_opt_in_returns_sir_2025_5047(lib):
    assert lib.state == "CO" and lib.provisional
    assert lib.regions == ["Alpine", "MidElevation", "Plateau", "SubAlpine"]
    assert len(lib.equations) == 32
    assert lib.source_reports[0]["doi"] == "10.3133/sir20255047"
    assert lib.source_reports[0]["nss"]["adopted"] is False
    assert "NOT a flowfreq-developed set" in lib.notes
    for region in TABLE1:
        assert lib.aeps(region) == sorted(AEP_OF_RP.values(), reverse=True)


@pytest.mark.parametrize("region", sorted(TABLE1))
def test_every_value_matches_table1(lib, region):
    for rp, b0, b1, terms, n, _pr2, sep, sme in TABLE1[region]:
        eq = lib.equation(region, AEP_OF_RP[rp])
        assert eq.intercept == b0
        assert [v.code for v in eq.variables] == ["DRNAREA"] + [c for c, _ in terms]
        assert list(eq.coefficients) == [b1] + [k for _, k in terms]
        assert eq.n_sites == n
        assert eq.sep_log == pytest.approx(_pct_to_log(sep), abs=5e-7)
        assert eq.model_error_variance == pytest.approx(_pct_to_log(sme) ** 2, abs=5e-8)
        assert eq.covariance is None and eq.avp is None
        assert eq.citation.table.startswith("Table 1") and eq.citation.year == 2026
        assert "2025-5047" in eq.citation.publication


def test_variables_carry_the_pooled_published_ranges(lib):
    for eq in lib.equations:
        for v in eq.variables:
            lo, hi, transform = RANGES[v.code]
            assert (v.minimum, v.maximum, v.transform) == (lo, hi, transform)
            assert v.scale == 1.0 and v.offset == 0.0


@pytest.mark.parametrize("region", sorted(NSS_STAGING))
def test_reproduces_nss_staging_draft(lib, region):
    for inputs, values in NSS_STAGING[region]:
        basin = _from_nss(inputs)
        for rp, expected in zip(AEP_OF_RP, values):
            got = evaluate(lib.equation(region, AEP_OF_RP[rp]), basin).flow_cfs
            # NSS rounds to three significant figures.
            half_unit = 0.5 * 10.0 ** (math.floor(math.log10(expected)) - 2)
            assert abs(got - expected) <= half_unit * 1.0001, (region, rp, got, expected)


@pytest.mark.parametrize("region", sorted(NSS_STAGING))
def test_quantiles_rise_with_rarer_aeps(lib, region):
    for inputs, _ in NSS_STAGING[region]:
        basin = _from_nss(inputs)
        flows = [evaluate(lib.equation(region, a), basin).flow_cfs for a in AEP_OF_RP.values()]
        assert flows == sorted(flows)


def test_intervals_fall_back_to_sep(lib):
    est = evaluate(lib.equation("Alpine", 0.01), _from_nss(NSS_STAGING["Alpine"][0][0]))
    assert est.interval_method == "average standard error of prediction"
    lo, hi = est.interval
    assert lo < est.flow_cfs < hi


def test_out_of_range_input_raises(lib):
    basin = {"DRNAREA": 100.0, "ELEV": 7000.0, "I24H100Y": 6.0}
    with pytest.raises(OutOfRangeError, match="I24H100Y"):
        evaluate(lib.equation("Plateau", 0.01), basin)
    est = evaluate(lib.equation("Plateau", 0.01), basin, allow_extrapolation=True)
    assert est.extrapolated


def test_unpublished_aep_raises(lib):
    with pytest.raises(EquationsUnavailable):
        lib.equation("Alpine", 0.3)
