"""B-WLS/B-GLS regional skew tooling (``flowfreq.skew_study``).

The validation reproduces the Pacific Northwest CONSTANT model from its own
published inputs: SIR 2016-5083 appendix B, Table B1 (station skew, pseudo
record length and basin centroid of the 290 gages used), eq. B11 (the study's
correlation-distance model), and Tables B2/B3 (the result). The only input the
report does not print -- each gage's historical period, needed for the pseudo
concurrent record length (eq. B10) -- comes from the NWIS peak record through
WY2012 (``tools/build_pnw_skew_validation_data.py``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from flowfreq.skew_study import (
    PNW_CORRELATION_MODEL,
    CorrelationDistanceModel,
    concurrent_pseudo_record_length,
    distance_miles,
    effective_record_length,
    fit_regional_skew,
    gs_skew_variance,
    pseudo_record_length,
    skew_correlation_matrix,
    unbias_skew,
    unbiased_skew_variance,
)

B1 = Path(__file__).parent / "fixtures" / "skew_study" / "pnw_table_b1.csv"


# ----------------------------------------------------------------------------
# Building blocks
# ----------------------------------------------------------------------------


def test_gs_variance_limits():
    # Large n: Var[G] -> 6/n at G = 0.
    assert gs_skew_variance(10000, 0.0) == pytest.approx(6 / 10000, rel=1e-2)
    assert gs_skew_variance(50, 0.5) > gs_skew_variance(50, 0.0)
    with pytest.raises(ValueError):
        gs_skew_variance(0, 0.0)


def test_unbiasing():
    assert unbias_skew(-0.5, 60) == pytest.approx(-0.55)
    assert unbiased_skew_variance(60, 0.1) == pytest.approx(1.1**2 * gs_skew_variance(60, 0.1))


def test_effective_record_length_reproduces_pnw():
    # SIR 2016-5083 p. 52: AVPnew 0.18 "corresponds to an effective record length of 41 years".
    assert effective_record_length(0.18, -0.07) == pytest.approx(41.0, abs=0.05)
    with pytest.raises(ValueError):
        effective_record_length(0.0, 0.0)


def test_pseudo_record_length_clipping():
    assert pseudo_record_length(40, 0.2, 0.1) == 80
    assert pseudo_record_length(40, 0.2, 0.1, historical_period=60) == 60
    assert pseudo_record_length(40, 0.1, 0.2) == 40  # never below Ps
    with pytest.raises(ValueError):
        pseudo_record_length(0, 0.1, 0.1)


def test_concurrent_record_length_b10():
    cy = concurrent_pseudo_record_length([1950, 1970], [2012, 2012], [50, 30], [63, 43])
    assert cy[0, 1] == pytest.approx(43 * (50 / 63) * (30 / 43))
    assert concurrent_pseudo_record_length([1900, 1960], [1950, 2000], [40, 40], [51, 41])[
        0, 1
    ] == pytest.approx(0.0)


def test_distance():
    d = distance_miles([46.87, 46.87], [-113.93, -112.93])
    assert d[0, 1] == pytest.approx(47.3, abs=0.3) and d[0, 0] == 0


def test_correlation_model_fit_recovers_parameters():
    d = np.linspace(0, 600, 200)
    r = PNW_CORRELATION_MODEL.rho(d)
    m = CorrelationDistanceModel.fit(d, r)
    assert (m.a, m.b, m.c) == pytest.approx((0.21, -0.17, -0.0058), abs=1e-3)
    with pytest.raises(ValueError):
        CorrelationDistanceModel.fit([1.0], [0.5])


def test_skew_correlation_matrix_properties():
    c = skew_correlation_matrix(
        [45, 45.1, 48], [-110, -110.1, -115], [50, 50, 50], [1960] * 3, [2009] * 3, [50] * 3,
        PNW_CORRELATION_MODEL, 2.8,
    )  # fmt: skip
    assert np.allclose(np.diag(c), 1) and np.allclose(c, c.T)
    assert c[0, 1] > c[0, 2] > 0


# ----------------------------------------------------------------------------
# fit_regional_skew
# ----------------------------------------------------------------------------


def test_synthetic_recovers_constant():
    rng = np.random.default_rng(1)
    n = 400
    prl = rng.integers(40, 100, n).astype(float)
    true = -0.2 + rng.normal(0, np.sqrt(0.1), n)
    g = true + rng.normal(0, np.sqrt(unbiased_skew_variance(prl, -0.2)))
    g = g / (1 + 6 / prl)  # back to the biased scale the function expects
    m = fit_regional_skew(g, prl)
    assert m.beta[0] == pytest.approx(-0.2, abs=0.05)
    assert m.sigma2 == pytest.approx(0.1, abs=0.04)
    assert m.mbv == pytest.approx(1.0, abs=0.01)  # independent: no misrepresentation


def test_synthetic_explanatory_variable():
    rng = np.random.default_rng(2)
    n = 300
    x = rng.uniform(0, 3, n)
    prl = np.full(n, 60.0)
    g = (0.3 - 0.2 * x + rng.normal(0, 0.2, n) + rng.normal(0, 0.3, n)) / 1.1
    m = fit_regional_skew(g, prl, x[:, None], ["x"])
    assert m.names == ["constant", "x"]
    assert m.beta[1] == pytest.approx(-0.2, abs=0.06)
    assert 0 < m.pseudo_r2 <= 1
    assert m.regional_skew([1.0]) == pytest.approx(m.beta[0] + m.beta[1])
    with pytest.raises(ValueError):
        m.regional_skew([])


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"skew": [0.1, 0.2], "prl": [30.0]}, "aligned"),
        ({"skew": [0.1, 0.2, 0.3], "prl": [30.0, -1.0, 40.0]}, "positive"),
        ({"skew": [0.1, 0.2], "prl": [30.0, 30.0]}, "cannot fit"),
    ],
)
def test_fit_errors(kwargs, match):
    with pytest.raises(ValueError, match=match):
        fit_regional_skew(**kwargs)


def test_non_psd_correlation_raises():
    c = np.array([[1, 0.99, -0.99], [0.99, 1, 0.99], [-0.99, 0.99, 1]])
    with pytest.raises(ValueError, match="semidefinite"):
        fit_regional_skew([0.1, 0.2, 0.3], [50.0, 50.0, 50.0], correlation=c)


# ----------------------------------------------------------------------------
# Validation: Pacific Northwest CONSTANT model (SIR 2016-5083 Tables B2, B3)
# ----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def pnw():
    d = pd.read_csv(B1, dtype={"site_no": str})
    u = d[d["used"] == "Yes"].reset_index(drop=True)
    c = skew_correlation_matrix(
        u["lat"], u["lon"], u["prl"], u["yb"], u["ye"], u["ph"], PNW_CORRELATION_MODEL, 2.8
    )
    return u, fit_regional_skew(u["skew"], u["prl"], correlation=c, kappa=2.8)


def test_pnw_table_b1_fixture():
    d = pd.read_csv(B1, dtype={"site_no": str})
    assert len(d) == 461
    assert d["used"].value_counts().to_dict() == {"Yes": 290, "no-P": 135, "no-R": 36}
    assert (d.loc[d["used"] == "Yes", "prl"] >= 35).all()


def test_pnw_constant_model_table_b2(pnw):
    _, m = pnw
    # Table B2: -0.07 (0.10); sigma2 0.17 (0.022); ASEV 0.010; AVPnew 0.18; pseudo-R2 0.
    assert round(m.beta[0], 2) == -0.07
    assert round(m.beta_sd[0], 2) == 0.10
    assert round(m.sigma2, 2) == 0.17
    assert round(m.sigma2_sd, 3) == 0.022
    assert m.asev == pytest.approx(0.010, abs=0.001)
    assert round(m.avp_new, 2) == 0.18
    assert m.pseudo_r2 == 0.0
    # p. 52: effective record length 41 years.
    assert round(m.effective_record_length) == 41


def test_pnw_pseudo_anova_table_b3(pnw):
    _, m = pnw
    assert round(float(np.sum(m.sampling_variance))) == 35  # sampling error SS
    assert round(m.n_sites * m.sigma2) in (49, 50)  # model error SS: 50
    assert round(m.evr, 1) == 0.7
    # MBV* published as 10; the historical periods here are approximated from
    # NWIS (see module docstring), which moves the cross-correlations a little.
    assert 8.5 <= m.mbv <= 11


def test_pnw_high_influence_sites(pnw):
    u, m = pnw
    published = [
        "12354000", "12311000", "12358500", "12325500", "14353500", "14353000", "13336500",
        "14158790", "12043000", "12344000", "14091500", "13112000", "12342500", "12354500",
        "12039500", "12121600", "14032000", "14166500", "13120000",
    ]  # fmt: skip
    ours = list(u["site_no"][np.argsort(-m.influence)[:19]])
    assert ours[:2] == published[:2]
    assert len(set(ours) & set(published)) >= 17


def test_pnw_wls_alone_understates_the_constant_variance(pnw):
    u, m = pnw
    wls = fit_regional_skew(u["skew"], u["prl"])
    assert wls.mbv == pytest.approx(1.0, abs=0.01)
    assert m.beta_sd[0] ** 2 / wls.beta_sd[0] ** 2 > 5


# ----------------------------------------------------------------------------
# station_skew (flowfreq's own EMA)
# ----------------------------------------------------------------------------


def _record(n: int = 60, seed: int = 3) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {
            "water_year": np.arange(1950, 1950 + n),
            "peak_flow_cfs": 10 ** rng.normal(3.0, 0.25, n),
            "qualification_code": [""] * n,
        }
    )


def test_station_skew_systematic_record():
    from flowfreq.bulletin17c import Bulletin17C
    from flowfreq.skew_study import station_skew

    rec = _record()
    s = station_skew("00000001", rec)
    res = Bulletin17C(rec["peak_flow_cfs"].to_numpy(), rec["water_year"].to_numpy()).run_analysis()
    assert s.skew == pytest.approx(res.skew_station)
    assert (s.yb, s.ye, s.ph) == (1950, 2009, 60)
    assert s.n_systematic <= s.prl <= s.ph
    assert s.skew_mse > 0 and s.at_site_option in ("ADJE", "B17B")


def test_station_skew_historic_peak_lengthens_prl():
    from flowfreq.skew_study import station_skew

    rec = _record()
    hist = pd.DataFrame(
        {"water_year": [1900], "peak_flow_cfs": [20000.0], "qualification_code": ["7"]}
    )
    s = station_skew("00000002", pd.concat([hist, rec], ignore_index=True))
    assert s.yb == 1900 and s.ph == 110
    assert s.prl > s.n_systematic


MT_DIR = Path(__file__).resolve().parents[1] / "data" / "skew_study"


def test_montana_provisional_study_reproduces_from_its_inputs():
    """The committed Montana results follow from the committed inputs."""
    import json

    res = json.loads((MT_DIR / "montana_provisional_v1_results.json").read_text(encoding="utf-8"))
    assert "NOT a USGS study" in res["status"]
    d = pd.read_csv(MT_DIR / "montana_provisional_v1_inputs.csv", dtype={"site_no": str})
    u = d[d["excluded"].isna()].reset_index(drop=True)
    assert len(u) == res["n_sites"]
    cm = res["correlation_model"]
    model = CorrelationDistanceModel(cm["a"], cm["b"], cm["c"])
    c = skew_correlation_matrix(
        u["lat"], u["lon"], u["prl"], u["yb"], u["ye"], u["ph"], model, res["kappa"]
    )
    m = fit_regional_skew(u["skew"], u["prl"], correlation=c)
    pub = res["models"]["CONSTANT"]
    assert m.beta[0] == pytest.approx(pub["beta_constant"], abs=1e-4)
    assert m.sigma2 == pytest.approx(pub["sigma2"], abs=1e-4)
    assert m.avp_new == pytest.approx(pub["avp_new"], abs=1e-4)


def test_montana_provisional_study_is_not_a_table_value():
    from flowfreq.regional_skew import RegionalSkewUnavailable, load_table, regional_skew_for

    df = load_table()
    assert not ((df["state"] == "MT") & (df["status"] == "verified")).any()
    with pytest.raises(RegionalSkewUnavailable):
        regional_skew_for("MT")


def test_station_skew_empty_record_raises():
    from flowfreq.skew_study import station_skew

    with pytest.raises(ValueError):
        station_skew("00000003", _record().iloc[0:0])
