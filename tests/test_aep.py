"""AEPs that two sources write differently are the same flood: 0.67 and 0.667, 0.43 and 0.429.

A Bulletin 17C fit tabulates the former and a regression report the latter. Each place that
pairs the two used to compare exactly, and dropped, mislabelled or left unweighted the
1.5-year flood without saying so.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq import estimate_from_characteristics
from flowfreq.aep import AEP_MATCH, match_aep, same_aep
from flowfreq.future_flow import apply_change_factors, factor_set_from_dict
from flowfreq.gage_weighting import weight_with_regression
from flowfreq.reporting import Paragraph, Table, build_section
from flowfreq.transpose import RegressionExponents
from flowfreq.workflow import run_ffa
from tests.test_reporting_sites import FLOWS, GAGE, SKEW, YEARS

IDAHO = {
    "state": "ID",
    "lat": 48.8,
    "lon": -116.2,
    "region": "1_2",
    "characteristics": {"DRNAREA": 205.0, "PRECPRIS10": 34.0},
}


@pytest.mark.parametrize(
    "a,b,same",
    [
        (0.67, 0.667, True),
        (0.43, 0.429, True),
        (0.5, 0.5, True),
        (0.005, 0.002, False),
        (0.01, 0.005, False),
        (0.2, 0.1, False),
        (0.99, 0.995, True),  # within 0.75 percent; no source tabulates both
    ],
)
def test_same_aep(a, b, same):
    assert same_aep(a, b) is same and same_aep(b, a) is same


def test_match_aep_takes_the_nearest_and_refuses_the_far():
    assert AEP_MATCH == 0.0075
    assert match_aep([0.8, 0.667, 0.5], 0.67) == 0.667
    assert match_aep([0.667, 0.67], 0.67) == 0.67
    assert match_aep([0.005, 0.01], 0.002) is None
    assert match_aep([], 0.5) is None


@pytest.fixture(scope="module")
def fit():
    return run_ffa(np.array(FLOWS, dtype=float), np.array(YEARS), stationarity=False, **SKEW)


def test_the_1_5_year_flood_is_weighted_with_the_regression(fit):
    estimate = estimate_from_characteristics(
        "ID", IDAHO["characteristics"], region="1_2", lat=48.8, lon=-116.2
    )
    regression = [per_region[0][1] for per_region in estimate.estimates.values()]
    frame = weight_with_regression(fit["b17c"].results, regression)
    row = frame[np.isclose(frame["aep"], 0.67)].iloc[0]
    assert np.isfinite(row["weighted_flow_cfs"])
    assert row["note"] == "regression at AEP 0.667"
    exact = frame[np.isclose(frame["aep"], 0.5)].iloc[0]
    assert exact["note"] == "" and np.isfinite(exact["weighted_flow_cfs"])


def test_the_weighted_report_section_keeps_and_explains_the_row():
    built = build_section("weighted", {"gage": dict(GAGE), "regression": IDAHO, "detail": False})
    (table,) = [
        b
        for b in built.blocks
        if isinstance(b, Table) and b.caption.startswith("At-site and regression flood quantiles")
    ]
    rows = {row[0]: row for row in table.rows}
    assert "67" in rows and rows["67"][-1] == "regression at AEP 0.667"
    assert table.header[-1] == "Note"
    assert any(isinstance(b, Paragraph) for b in built.blocks)


def test_a_published_endpoint_written_differently_is_not_an_extrapolation():
    exponents = RegressionExponents(
        aeps=[0.667, 0.5, 0.1, 0.01], exponents=[0.95, 0.93, 0.88, 0.84], citation="made up"
    )
    result = exponents.exponent_at([0.67, 0.5, 0.3, 0.002], extrapolation="clamp")
    assert result.loc[0, "exponent"] == pytest.approx(0.95)
    assert list(result["source"]) == ["published", "published", "interpolated", "extrapolated"]
    assert list(result["aep"]) == [0.67, 0.5, 0.3, 0.002]  # reported as asked
    # 'error' no longer refuses the 1.5-year flood, and still refuses a real one.
    exponents.exponent_at([0.67], extrapolation="error")
    with pytest.raises(ValueError, match="outside the published range"):
        exponents.exponent_at([0.002], extrapolation="error")


FACTORS = {
    "name": "made-up",
    "aeps": [0.667, 0.5, 0.1, 0.01],
    "factors": [1.05, 1.06, 1.10, 1.20],
    "scenario": "S",
    "horizon": "2050",
    "citation": "Example, A.B., 2024, made up for the tests.",
}


def test_a_change_factor_at_0_667_serves_a_quantile_at_0_67():
    factors = factor_set_from_dict(FACTORS)
    values, extrapolated = factors.factor_at([0.67, 0.5, 0.002])
    assert values[0] == pytest.approx(1.05) and not extrapolated[0]
    assert extrapolated[2]
    table = apply_change_factors(
        pd.DataFrame({"aep": [0.67, 0.01], "flow_cfs": [100.0, 400.0]}), factors
    ).table
    assert list(table["extrapolated"]) == [False, False]
    assert table.loc[0, "future_flow_cfs"] == pytest.approx(105.0)


@pytest.mark.parametrize("bad", [[50.0, 10.0], [1.0], [0.0], [float("nan")]])
def test_an_aep_given_in_percent_is_refused(bad):
    with pytest.raises(ValueError, match="fraction in \\(0, 1\\)"):
        factor_set_from_dict(FACTORS).factor_at(bad)


def test_the_future_flow_section_keeps_the_1_5_year_row_as_published():
    built = build_section(
        "future_flow", {"gage": dict(GAGE), "factor_set": dict(FACTORS), "detail": False}
    )
    (table,) = [
        b
        for b in built.blocks
        if isinstance(b, Table) and b.caption == "Flood quantiles under the future condition."
    ]
    rows = {row[0]: row for row in table.rows}
    assert rows["67"][2] == "1.050" and rows["67"][4] == "published"
    asked = build_section(
        "future_flow",
        {"gage": dict(GAGE), "factor_set": dict(FACTORS), "aeps": [0.67], "detail": False},
    )
    assert not any(
        isinstance(b, Paragraph) and "nearest published factor" in b.text for b in asked.blocks
    )
