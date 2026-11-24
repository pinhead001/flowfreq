"""Texas's OmegaEM equations: SIR 2009-5087 Table 3, region ``Statewide OmegaEM``.

Two sources in the report are compared: Table 3, typed from the page image and stored, and
Appendix 3's R output, parsed from the PDF text into
``fixtures/regression/texas_omegaem_appendix3.json``. The stations of
``texas_omegaem_gages.csv`` are the report's own (its Appendix 1 data files): the 502 of 638
that the gage catalog locates in a quadrangle of figure 2 with one OmegaEM value.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd
import pytest
from scipy.stats import t

from flowfreq.regression import evaluate
from flowfreq.regression.library import EquationsUnavailable, load_state
from flowfreq.ungaged import estimate_from_characteristics

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "regression"
REGION = "Statewide OmegaEM"
YEARS = (2, 5, 10, 25, 50, 100, 200, 250, 500)
#: The report's worked example: 08190000, Nueces River at Laguna (p. 17).
LAGUNA = {"CONTDA": 737.0, "PRECIPCENT": 24.5, "MCSLOPE": 0.00326, "OmegaEM": 0.33}


@pytest.fixture(scope="module")
def lib() -> Any:
    return load_state("TX")


@pytest.fixture(scope="module")
def appendix() -> Dict[str, Any]:
    path = FIXTURES / "texas_omegaem_appendix3.json"
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def gages() -> pd.DataFrame:
    return pd.read_csv(FIXTURES / "texas_omegaem_gages.csv", dtype={"STATION": str})


def _equation(lib: Any, years: int) -> Any:
    return lib.equation(REGION, 1.0 / years)


def test_nine_equations_and_no_nss_region(lib: Any) -> None:
    mine = [e for e in lib.equations if e.region_code == REGION]
    assert [e.aep for e in mine] == [1.0 / t for t in YEARS]
    assert "NSS" not in mine[0].region_name
    assert all(e.n_sites == 638 for e in mine)
    assert all("2009-5087" in e.citation.publication for e in mine)


def test_the_report_example(lib: Any) -> None:
    """P. 17: about 145,000 cubic feet per second."""
    flow = evaluate(_equation(lib, 100), LAGUNA).flow_cfs
    assert round(flow, -3) == 145000


@pytest.mark.parametrize("years", YEARS)
def test_table_3_is_appendix_3_rounded(lib: Any, appendix: Dict[str, Any], years: int) -> None:
    eq, printed = _equation(lib, years), appendix[str(years)]
    b0, area, precip, slope, omega = printed["coef"]
    assert eq.intercept == pytest.approx(b0, abs=0.00501)
    assert eq.coefficients[0] == float(f"{area:.4g}")
    assert eq.coefficients[1:] == (round(precip, 3), round(slope, 3), round(omega, 3))
    assert eq.variables[0].exponent == printed["power"]


@pytest.mark.parametrize("years", YEARS)
def test_matrix_gives_the_printed_standard_errors(
    lib: Any, appendix: Dict[str, Any], years: int
) -> None:
    eq, printed = _equation(lib, years), appendix[str(years)]
    assert eq.model_error_variance == pytest.approx(printed["rse"] ** 2)
    covariance = np.array(eq.covariance)
    errors = np.sqrt(np.diag(covariance))
    np.testing.assert_allclose(errors, printed["se"], rtol=4e-4)
    np.testing.assert_allclose(covariance, covariance.T)
    # Every term, not only the diagonal.
    matrix = np.array(printed["matrix"])
    np.testing.assert_allclose(covariance, printed["rse"] ** 2 * matrix, rtol=2e-6)


@pytest.mark.parametrize("years", YEARS)
def test_limits_are_table_1_and_figure_2(lib: Any, years: int) -> None:
    """Table 1: 0.1 to 9,329 mi2, 8 to 57 in, slope 0.00023 to 0.0703; figure 2's least and most."""
    limits = [(v.code, v.units, v.minimum, v.maximum) for v in _equation(lib, years).variables]
    assert limits == [
        ("CONTDA", "mi2", 0.1, 9329.0),
        ("PRECIPCENT", "in", 8.0, 57.0),
        ("MCSLOPE", "ft/ft", 0.00023, 0.0703),
        ("OmegaEM", "dimensionless", -0.503, 0.386),
    ]


def test_the_interval_is_students_t_on_633_degrees_of_freedom(lib: Any) -> None:
    """Equation 7 at the report's example, 90 percent: 638 stations less five coefficients."""
    est = evaluate(_equation(lib, 100), LAGUNA)
    assert est.interval is not None and est.interval_level == 0.9
    half = t.ppf(0.95, 633) * math.sqrt(est.prediction_variance)
    assert est.interval[0] == pytest.approx(est.flow_cfs / 10.0**half, rel=1e-9)
    assert est.interval[1] == pytest.approx(est.flow_cfs * 10.0**half, rel=1e-9)
    assert [round(x, -1) for x in est.interval] == [45780.0, 457040.0]


@pytest.mark.parametrize("years", YEARS)
def test_stations_scatter_about_the_equation_as_the_report_says(
    lib: Any, appendix: Dict[str, Any], gages: pd.DataFrame, years: int
) -> None:
    """Residuals from the report's own station flows: no bias, and the printed standard error."""
    eq = _equation(lib, years)
    residuals = np.array(
        [
            math.log10(row[f"trimQ{years}"])
            - math.log10(
                evaluate(
                    eq,
                    {
                        "CONTDA": row["CDA"],
                        "PRECIPCENT": row["MAP"],
                        "MCSLOPE": row["MCS"],
                        "OmegaEM": row["OMEGA"],
                    },
                ).flow_cfs
            )
            for _, row in gages.iterrows()
        ]
    )
    assert len(residuals) == 502
    assert abs(residuals.mean()) < 0.025
    assert residuals.std() == pytest.approx(appendix[str(years)]["rse"], rel=0.06)


@pytest.mark.parametrize("years", YEARS)
def test_rounded_coefficients_move_a_station_flow_by_two_percent_at_most(
    lib: Any, appendix: Dict[str, Any], gages: pd.DataFrame, years: int
) -> None:
    eq, printed = _equation(lib, years), appendix[str(years)]
    b = printed["coef"]
    area, precip = gages["CDA"] ** printed["power"], np.log10(gages["MAP"])
    slope, omega = np.log10(gages["MCS"]), gages["OMEGA"]
    unrounded = b[0] + b[1] * area + b[2] * precip + b[3] * slope + b[4] * omega
    c = eq.coefficients
    stored = eq.intercept + c[0] * area + c[1] * precip + c[2] * slope + c[3] * omega
    assert 10 ** abs(unrounded - stored).max() - 1 < 0.02


def test_interval_is_the_report_form(lib: Any, appendix: Dict[str, Any]) -> None:
    """Equations 7 and 8: s times the root of 1 + h0, with h0 from the printed matrix."""
    eq, printed = _equation(lib, 100), appendix["100"]
    x = np.array([1.0, 737.0 ** printed["power"], math.log10(24.5), math.log10(0.00326), 0.33])
    h0 = float(x @ np.array(printed["matrix"]) @ x)
    assert 0.0 < h0 < 0.06  # the printed greatest leverage of a station is 0.05
    variance = eq.model_error_variance + float(x @ np.array(eq.covariance) @ x)
    assert variance == pytest.approx(printed["rse"] ** 2 * (1.0 + h0), rel=1e-4)


def test_it_is_used_only_when_asked_for(lib: Any) -> None:
    tiers = [code for tier in lib.region_preference.tiers for code in tier]
    assert REGION not in tiers
    assert REGION not in lib.region_limits


def test_the_total_drainage_area_is_not_taken_for_the_contributing_one(lib: Any) -> None:
    """P. 3: area is contributing. StreamStats' DRNAREA is total, many times as large in places."""
    total = {**{k: v for k, v in LAGUNA.items() if k != "CONTDA"}, "DRNAREA": 737.0}
    with pytest.raises(EquationsUnavailable, match="CONTDA"):
        estimate_from_characteristics("TX", total, region=REGION, lat=29.43, lon=-100.0)


def test_an_estimate_carries_what_the_set_assumes() -> None:
    estimate = estimate_from_characteristics("TX", LAGUNA, region=REGION, lat=29.43, lon=-100.0)
    for words in (
        "undeveloped",
        "contributing drainage area (CONTDA), not the total",
        "feet per foot",
        "not a 10-85 slope",
        "too narrow",
    ):
        assert words in estimate.advisory, words


def test_the_statewide_set_says_the_authors_prefer_this_one() -> None:
    estimate = estimate_from_characteristics(
        "TX", {"DRNAREA": 737.0, "PRECIPCENT": 24.5}, region="Statewide", lat=29.43, lon=-100.0
    )
    assert "SIR 2009-5087" in estimate.advisory and "'Statewide OmegaEM'" in estimate.advisory


def test_other_texas_sets_carry_no_caution(lib: Any) -> None:
    assert set(lib.region_cautions) == {"Statewide", REGION}


def test_estimate_from_characteristics(lib: Any) -> None:
    estimate = estimate_from_characteristics("TX", LAGUNA, region=REGION, lat=29.43, lon=-100.0)
    table = estimate.quantiles.set_index("aep")
    flow = evaluate(_equation(lib, 100), LAGUNA).flow_cfs
    assert table.loc[0.01, "flow_cfs"] == pytest.approx(flow)


@pytest.mark.parametrize(
    "code,value", [("CONTDA", 9400.0), ("MCSLOPE", 0.1), ("OmegaEM", 0.5), ("PRECIPCENT", 60.0)]
)
def test_out_of_range_is_refused(lib: Any, code: str, value: float) -> None:
    with pytest.raises(ValueError, match=code):
        evaluate(_equation(lib, 100), {**LAGUNA, code: value})


def test_the_notes_say_the_slope_must_be_supplied(lib: Any) -> None:
    assert "STATEWIDE OMEGAEM: it needs MCSLOPE" in lib.notes
    assert "StreamStats returns no channel slope in Texas" in lib.notes
