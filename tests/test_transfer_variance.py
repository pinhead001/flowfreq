"""Tests for flowfreq.transfer_variance on synthetic gage networks.

Every network here is built so the dissimilarity variance is known exactly.
Gages sit on chains (one stream each). Gage ``g`` on chain ``c`` has true
log quantile

    mu_g = a + b * log10(A_g) + beta_c * log10(A_g) + gamma_c * P_g + eta_g

with ``beta_c ~ N(0, s1)``, ``gamma_c ~ N(0, s2)`` and ``eta_g ~ N(0, s0/2)``.
For two gages on one chain, the transfer error ``mu_i + b L - mu_j`` then has
variance ``s0 + s1 L^2 + s2 dP^2``: the model form ``fit_dissimilarity``
fits, with known coefficients. At-site estimates add sampling error with
known (co)variance.
"""

from __future__ import annotations

import itertools
from types import SimpleNamespace
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import pytest

import flowfreq.transfer_variance as tv
from flowfreq.regression.gls import WregCorrelation
from flowfreq.transfer_variance import (
    BlendResult,
    DissimilarityModel,
    PowerProfile,
    bootstrap_band,
    build_pairs,
    combine,
    estimate_donor_correlation,
    fit_dissimilarity,
    leave_one_out,
    profile_power,
    sampling_covariance_frame,
    solve_band,
)
from flowfreq.transpose import RegressionExponents, transpose_frequency

B = 0.8
EXPONENTS = RegressionExponents.constant(B, citation="synthetic: b = 0.8 by construction")


def make_network(
    n_chains: int,
    gages_per_chain: int,
    *,
    s0: float,
    s1: float,
    s2: float = 0.0,
    sampling_var: float = 0.004,
    sampling_rho: float = 0.6,
    regression_var: float = 0.05,
    seed: int = 0,
) -> Tuple[pd.DataFrame, List[Tuple[str, str]], pd.DataFrame]:
    """Gages, their connections (every pair on a chain), and the sampling covariance."""
    rng = np.random.default_rng(seed)
    rows = []
    connections: List[Tuple[str, str]] = []
    for c in range(n_chains):
        beta = rng.normal(0.0, np.sqrt(s1))
        gamma = rng.normal(0.0, np.sqrt(s2)) if s2 > 0 else 0.0
        log_a = np.sort(rng.uniform(0.5, 3.0, gages_per_chain))
        precip = rng.normal(30.0, 3.0, gages_per_chain)
        ids = [f"c{c:03d}g{k}" for k in range(gages_per_chain)]
        for k, site in enumerate(ids):
            mu = (
                2.0
                + B * log_a[k]
                + beta * log_a[k]
                + gamma * precip[k]
                + rng.normal(0.0, np.sqrt(s0 / 2.0))
            )
            rows.append(
                {
                    "site_no": site,
                    "chain": c,
                    "drainage_area_sqmi": 10.0 ** log_a[k],
                    "precip": precip[k],
                    "mu": mu,
                    "variance": sampling_var,
                }
            )
        connections.extend(itertools.combinations(ids, 2))
    gages = pd.DataFrame(rows)
    sites = gages["site_no"].tolist()
    chain = gages["chain"].to_numpy()
    same = chain[:, None] == chain[None, :]
    cov = np.where(same, sampling_rho * sampling_var, 0.0)
    np.fill_diagonal(cov, sampling_var)
    eps = rng.multivariate_normal(np.zeros(len(sites)), cov)
    gages["log_q"] = gages["mu"] + eps
    gages["regression_log_q"] = gages["mu"] + rng.normal(0.0, np.sqrt(regression_var), len(sites))
    gages["regression_variance"] = regression_var
    return gages, connections, pd.DataFrame(cov, index=sites, columns=sites)


# --------------------------------------------------------------------------- pairs


class TestBuildPairs:
    def test_hand_computed_pair(self):
        gages = pd.DataFrame(
            {
                "site_no": ["a", "b"],
                "drainage_area_sqmi": [100.0, 400.0],
                "log_q": [3.0, 3.5],
                "variance": [0.01, 0.02],
                "precip": [30.0, 25.0],
            }
        )
        pairs = build_pairs(
            gages, [("a", "b")], EXPONENTS, 0.01, exponent_variance=0.002, attributes=["precip"]
        )
        assert len(pairs) == 1
        row = pairs.iloc[0]
        L = np.log10(4.0)
        assert row["donor"] == "a" and row["receiver"] == "b"
        assert row["log_ar"] == pytest.approx(L)
        assert row["transfer"] == pytest.approx(3.0 + B * L)
        assert row["error"] == pytest.approx(3.0 + B * L - 3.5)
        assert row["known_variance"] == pytest.approx(0.01 + 0.02 + L**2 * 0.002)
        assert row["d_precip"] == pytest.approx(-5.0)

    def test_sampling_covariance_is_subtracted_twice(self):
        gages = pd.DataFrame(
            {"site_no": ["a", "b"], "drainage_area_sqmi": [100.0, 200.0], "log_q": [3.0, 3.2]}
        )
        cov = pd.DataFrame([[0.01, 0.006], [0.006, 0.02]], index=["a", "b"], columns=["a", "b"])
        row = build_pairs(gages, [("a", "b")], EXPONENTS, 0.01, sampling_cov=cov).iloc[0]
        assert row["known_variance"] == pytest.approx(0.01 + 0.02 - 2 * 0.006)

    def test_covariance_order_follows_labels_not_rows(self):
        gages = pd.DataFrame(
            {"site_no": ["a", "b"], "drainage_area_sqmi": [100.0, 200.0], "log_q": [3.0, 3.2]}
        )
        cov = pd.DataFrame([[0.02, 0.006], [0.006, 0.01]], index=["b", "a"], columns=["b", "a"])
        row = build_pairs(gages, [("a", "b")], EXPONENTS, 0.01, sampling_cov=cov).iloc[0]
        assert row["known_variance"] == pytest.approx(0.01 + 0.02 - 2 * 0.006)

    def test_matches_transpose_frequency(self):
        """The pair's transfer is transpose_frequency's arithmetic, in logs."""
        donor = SimpleNamespace(
            quantiles=pd.DataFrame({"aep": [0.01], "return_period": [100.0], "flow_cfs": [5000.0]}),
            confidence_limits=None,
        )
        out = transpose_frequency(donor, 100.0, 140.0, EXPONENTS)  # type: ignore[arg-type]
        gages = pd.DataFrame(
            {
                "site_no": ["a", "b"],
                "drainage_area_sqmi": [100.0, 140.0],
                "log_q": [np.log10(5000.0), 0.0],
                "variance": [0.01, 0.01],
            }
        )
        row = build_pairs(gages, [("a", "b")], EXPONENTS, 0.01).iloc[0]
        assert 10.0 ** row["transfer"] == pytest.approx(out.quantiles["flow_cfs"].iloc[0])

    def test_reversed_pair_has_the_same_squared_error(self):
        gages = pd.DataFrame(
            {
                "site_no": ["a", "b"],
                "drainage_area_sqmi": [100.0, 300.0],
                "log_q": [3.0, 3.3],
                "variance": [0.01, 0.01],
            }
        )
        ab = build_pairs(gages, [("a", "b")], EXPONENTS, 0.01).iloc[0]
        ba = build_pairs(gages, [("b", "a")], EXPONENTS, 0.01).iloc[0]
        assert ba["error"] == pytest.approx(-ab["error"])
        assert ba["log_ar"] == pytest.approx(-ab["log_ar"])

    def test_duplicate_connection_raises(self):
        gages = pd.DataFrame(
            {
                "site_no": ["a", "b"],
                "drainage_area_sqmi": [1.0, 2.0],
                "log_q": [1.0, 1.0],
                "variance": [0.01, 0.01],
            }
        )
        with pytest.raises(ValueError, match="more than once"):
            build_pairs(gages, [("a", "b"), ("b", "a")], EXPONENTS, 0.01)

    def test_unknown_and_self_connections_raise(self):
        gages = pd.DataFrame(
            {"site_no": ["a"], "drainage_area_sqmi": [1.0], "log_q": [1.0], "variance": [0.01]}
        )
        with pytest.raises(KeyError, match="zz"):
            build_pairs(gages, [("a", "zz")], EXPONENTS, 0.01)
        with pytest.raises(ValueError, match="itself"):
            build_pairs(gages, [("a", "a")], EXPONENTS, 0.01)

    def test_needs_variance_or_covariance(self):
        gages = pd.DataFrame(
            {"site_no": ["a", "b"], "drainage_area_sqmi": [1.0, 2.0], "log_q": [1.0, 1.1]}
        )
        with pytest.raises(ValueError, match="variance"):
            build_pairs(gages, [("a", "b")], EXPONENTS, 0.01)

    def test_wrong_probability_kind_raises(self):
        gages = pd.DataFrame(
            {
                "site_no": ["a", "b"],
                "drainage_area_sqmi": [1.0, 2.0],
                "log_q": [1.0, 1.1],
                "variance": [0.01, 0.01],
            }
        )
        duration = RegressionExponents.constant(0.9, citation="x", probability_kind="exceedance")
        with pytest.raises(ValueError, match="indexed by"):
            build_pairs(gages, [("a", "b")], duration, 0.01)


class TestSamplingCovarianceFrame:
    def test_labels_and_values_match_gls(self):
        from flowfreq.regression.gls import sampling_covariance

        gages = pd.DataFrame(
            {
                "site_no": ["a", "b", "c"],
                "sd": [0.25, 0.3, 0.28],
                "skew": [0.0, 0.1, -0.1],
                "n": [40, 60, 50],
                "lat": [48.0, 48.1, 48.5],
                "lon": [-120.0, -120.1, -120.3],
            }
        )
        concurrent = pd.DataFrame(
            [[40, 35, 20], [35, 60, 45], [20, 45, 50]],
            index=["a", "b", "c"],
            columns=["a", "b", "c"],
        )
        corr = WregCorrelation(alpha=0.002, theta=0.98)
        frame = sampling_covariance_frame(gages, 0.01, 0.12, corr, concurrent)
        expected = sampling_covariance(
            gages["sd"],
            gages["skew"],
            gages["n"],
            0.01,
            0.12,
            gages["lat"],
            gages["lon"],
            concurrent.to_numpy(),
            corr,
        )
        assert list(frame.index) == ["a", "b", "c"] and list(frame.columns) == ["a", "b", "c"]
        np.testing.assert_allclose(frame.to_numpy(), expected)

    def test_concurrent_is_aligned_by_label(self):
        gages = pd.DataFrame(
            {
                "site_no": ["a", "b"],
                "sd": [0.25, 0.3],
                "skew": [0.0, 0.0],
                "n": [40, 60],
                "lat": [48.0, 48.1],
                "lon": [-120.0, -120.1],
            }
        )
        corr = WregCorrelation(alpha=0.002, theta=0.98)
        forward = pd.DataFrame([[40, 30], [30, 60]], index=["a", "b"], columns=["a", "b"])
        backward = forward.loc[["b", "a"], ["b", "a"]]
        a = sampling_covariance_frame(gages, 0.01, 0.0, corr, forward)
        b = sampling_covariance_frame(gages, 0.01, 0.0, corr, backward)
        pd.testing.assert_frame_equal(a, b)


# --------------------------------------------------------------------------- fit


class TestFitDissimilarity:
    def test_recovers_known_coefficients(self):
        gages, conns, cov = make_network(600, 4, s0=0.02, s1=0.03, seed=1)
        pairs = build_pairs(gages, conns, EXPONENTS, 0.01, sampling_cov=cov)
        model = fit_dissimilarity(pairs, terms=["log_ar"])
        assert isinstance(model, DissimilarityModel)
        assert model.coefficients["const"] == pytest.approx(0.02, abs=0.006)
        assert model.coefficients["log_ar"] == pytest.approx(0.03, abs=0.01)
        assert model.n_pairs == 600 * 6
        assert model.n_gages == 600 * 4
        assert model.converged
        assert not model.bound_terms

    def test_recovers_a_precipitation_term(self):
        gages, conns, cov = make_network(600, 4, s0=0.02, s1=0.02, s2=0.0004, seed=2)
        pairs = build_pairs(gages, conns, EXPONENTS, 0.01, sampling_cov=cov, attributes=["precip"])
        model = fit_dissimilarity(pairs, terms=["log_ar", "d_precip"])
        assert model.coefficients["d_precip"] == pytest.approx(0.0004, rel=0.35)
        assert model.coefficients["log_ar"] == pytest.approx(0.02, abs=0.012)

    def test_ignoring_the_sampling_covariance_biases_s0_low(self):
        """Subtracting V_i + V_j without -2 C_ij over-subtracts the truth noise."""
        gages, conns, cov = make_network(600, 4, s0=0.02, s1=0.03, sampling_rho=0.8, seed=3)
        right = fit_dissimilarity(
            build_pairs(gages, conns, EXPONENTS, 0.01, sampling_cov=cov), terms=["log_ar"]
        )
        wrong = fit_dissimilarity(build_pairs(gages, conns, EXPONENTS, 0.01), terms=["log_ar"])
        assert wrong.coefficients["const"] < right.coefficients["const"] - 0.003

    def test_never_negative_and_reports_the_bind(self, caplog):
        gages, conns, cov = make_network(15, 3, s0=0.0, s1=0.0, sampling_var=0.02, seed=4)
        pairs = build_pairs(gages, conns, EXPONENTS, 0.01, sampling_cov=cov)
        # Force negatives: report a larger known variance than the truth.
        pairs["known_variance"] *= 3.0
        with caplog.at_level("WARNING", logger="flowfreq.transfer_variance"):
            model = fit_dissimilarity(pairs, terms=["log_ar"])
        assert all(v >= 0.0 for v in model.coefficients.values())
        assert "const" in model.bound_terms
        assert "held at zero" in caplog.text

    def test_variance_evaluates_the_fitted_form(self):
        model = DissimilarityModel(
            terms=("log_ar", "d_precip"),
            coefficients={"const": 0.01, "log_ar": 0.04, "d_precip": 0.001},
            bound_terms=(),
            n_pairs=10,
            n_gages=8,
            iterations=3,
            converged=True,
            residuals=pd.DataFrame(),
            by_ratio=pd.DataFrame(),
        )
        assert model.variance(log_ar=-0.5, d_precip=2.0) == pytest.approx(
            0.01 + 0.04 * 0.25 + 0.001 * 4.0
        )
        with pytest.raises(KeyError, match="d_precip"):
            model.variance(log_ar=0.5)

    def test_diagnostics(self):
        gages, conns, cov = make_network(300, 4, s0=0.02, s1=0.03, seed=5)
        pairs = build_pairs(gages, conns, EXPONENTS, 0.01, sampling_cov=cov)
        model = fit_dissimilarity(pairs, terms=["log_ar"])
        z = model.residuals["standardized"]
        assert len(z) == len(pairs)
        # Standardized by the fitted E[e^2], so the mean square is about 1.
        assert float(np.mean(z**2)) == pytest.approx(1.0, abs=0.1)
        assert {"abs_log_ar_mid", "mean_excess", "mean_fitted", "n_pairs"} <= set(
            model.by_ratio.columns
        )
        assert model.by_ratio["n_pairs"].sum() == len(pairs)

    def test_too_few_pairs_raises(self):
        pairs = pd.DataFrame(
            {
                "donor": ["a"],
                "receiver": ["b"],
                "log_ar": [0.1],
                "error": [0.1],
                "known_variance": [0.01],
            }
        )
        with pytest.raises(ValueError, match="pairs"):
            fit_dissimilarity(pairs, terms=["log_ar"])

    def test_unknown_term_raises(self):
        gages, conns, cov = make_network(20, 3, s0=0.02, s1=0.03, seed=6)
        pairs = build_pairs(gages, conns, EXPONENTS, 0.01, sampling_cov=cov)
        with pytest.raises(KeyError, match="d_canopy"):
            fit_dissimilarity(pairs, terms=["log_ar", "d_canopy"])


# --------------------------------------------------------------------------- blend


def _model(s0: float, s1: float) -> DissimilarityModel:
    return DissimilarityModel(
        terms=("log_ar",),
        coefficients={"const": s0, "log_ar": s1},
        bound_terms=(),
        n_pairs=100,
        n_gages=80,
        iterations=1,
        converged=True,
        residuals=pd.DataFrame(),
        by_ratio=pd.DataFrame(),
    )


class TestCombine:
    def test_one_donor_and_regression_is_inverse_variance(self):
        model = _model(0.02, 0.04)
        transfers = pd.DataFrame(
            {"donor": ["a"], "transfer": [3.0], "donor_variance": [0.01], "log_ar": [0.5]}
        )
        out = combine(model, transfers, regression=(3.3, 0.06), exponent_variance=0.002)
        assert isinstance(out, BlendResult)
        v_t = 0.01 + 0.25 * 0.002 + 0.02 + 0.04 * 0.25
        expected = (3.0 / v_t + 3.3 / 0.06) / (1 / v_t + 1 / 0.06)
        assert out.log_q == pytest.approx(expected)
        assert out.variance == pytest.approx(1 / (1 / v_t + 1 / 0.06))
        assert out.weights.sum() == pytest.approx(1.0)
        assert list(out.weights.index) == ["a", "regression"]
        lo, hi = out.interval
        assert lo < out.flow_cfs < hi
        assert np.log10(hi) - out.log_q == pytest.approx(1.6448536 * np.sqrt(out.variance))

    def test_two_donors_use_their_covariance(self):
        model = _model(0.02, 0.04)
        transfers = pd.DataFrame(
            {
                "donor": ["up", "down"],
                "transfer": [3.0, 3.2],
                "donor_variance": [0.01, 0.012],
                "log_ar": [0.3, -0.2],
            }
        )
        cov = pd.DataFrame(
            [[0.01, 0.007], [0.007, 0.012]], index=["up", "down"], columns=["up", "down"]
        )
        out = combine(
            model,
            transfers,
            regression=None,
            exponent_variance=0.002,
            donor_cov=cov,
            donor_correlation=0.3,
        )
        s_up = 0.02 + 0.04 * 0.09
        s_dn = 0.02 + 0.04 * 0.04
        off = 0.007 + 0.3 * (-0.2) * 0.002 + 0.3 * np.sqrt(s_up * s_dn)
        c = np.array(
            [
                [0.01 + 0.09 * 0.002 + s_up, off],
                [off, 0.012 + 0.04 * 0.002 + s_dn],
            ]
        )
        ci = np.linalg.inv(c)
        one = np.ones(2)
        w = ci @ one / (one @ ci @ one)
        np.testing.assert_allclose(out.covariance.to_numpy(), c)
        np.testing.assert_allclose(out.weights.to_numpy(), w)
        assert out.log_q == pytest.approx(w @ np.array([3.0, 3.2]))
        assert out.variance == pytest.approx(1.0 / (one @ ci @ one))

    def test_regression_correlation_enters_the_covariance(self):
        model = _model(0.02, 0.0)
        transfers = pd.DataFrame(
            {"donor": ["a"], "transfer": [3.0], "donor_variance": [0.01], "log_ar": [0.0]}
        )
        out = combine(model, transfers, regression=(3.1, 0.05), regression_correlation=0.4)
        assert out.covariance.loc["a", "regression"] == pytest.approx(0.4 * np.sqrt(0.03 * 0.05))

    def test_coverage_by_simulation(self):
        """A 90 % interval from the blend covers the truth about 90 % of the time."""
        rng = np.random.default_rng(11)
        model = _model(0.02, 0.04)
        hits = 0
        n = 2000
        for _ in range(n):
            truth = 3.0
            L = np.array([0.4, -0.3])
            s = 0.02 + 0.04 * L**2
            t = truth + rng.normal(0.0, np.sqrt(0.01 + s))
            r = truth + rng.normal(0.0, np.sqrt(0.05))
            transfers = pd.DataFrame(
                {"donor": ["u", "d"], "transfer": t, "donor_variance": [0.01, 0.01], "log_ar": L}
            )
            out = combine(model, transfers, regression=(r, 0.05))
            lo, hi = out.interval
            hits += lo <= 10.0**truth <= hi
        assert hits / n == pytest.approx(0.90, abs=0.02)

    def test_needs_an_estimate(self):
        with pytest.raises(ValueError, match="nothing to combine"):
            combine(
                _model(0.02, 0.04),
                pd.DataFrame(columns=["donor", "transfer", "donor_variance", "log_ar"]),
                regression=None,
            )

    def test_bad_correlation_raises(self):
        transfers = pd.DataFrame(
            {"donor": ["a"], "transfer": [3.0], "donor_variance": [0.01], "log_ar": [0.0]}
        )
        with pytest.raises(ValueError, match="correlation"):
            combine(
                _model(0.02, 0.0), transfers, regression=(3.0, 0.05), regression_correlation=1.5
            )


class TestEstimateDonorCorrelation:
    def test_shared_receiver_noise_gives_one_half(self):
        # eta is independent per gage, so with s1 = 0 two donors' dissimilarity
        # errors at one receiver share only the receiver's own eta: correlation 1/2.
        gages, conns, cov = make_network(800, 3, s0=0.04, s1=0.0, seed=12)
        pairs = build_pairs(gages, conns, EXPONENTS, 0.01, sampling_cov=cov)
        model = fit_dissimilarity(pairs, terms=["log_ar"])
        rho, n = estimate_donor_correlation(gages, conns, EXPONENTS, 0.01, model, sampling_cov=cov)
        assert n == 800 * 3
        assert rho == pytest.approx(0.5, abs=0.1)

    def test_none_without_triples(self):
        gages = pd.DataFrame(
            {
                "site_no": ["a", "b"],
                "drainage_area_sqmi": [1.0, 2.0],
                "log_q": [1.0, 1.1],
                "variance": [0.01, 0.01],
            }
        )
        rho, n = estimate_donor_correlation(gages, [("a", "b")], EXPONENTS, 0.01, _model(0.02, 0.0))
        assert rho is None and n == 0


# --------------------------------------------------------------------------- band


class TestSolveBand:
    def test_closed_form(self):
        model = _model(0.02, 0.04)
        band = solve_band(
            model, gage_variance=0.01, regression_variance=0.09, exponent_variance=0.002
        )
        L = np.sqrt((0.09 - 0.01 - 0.02) / (0.04 + 0.002))
        assert band.upper == pytest.approx(10.0**L)
        assert band.lower == pytest.approx(10.0**-L)
        assert band.half_width_log == pytest.approx(L)
        # At the edge the two variances are equal.
        v_t = 0.01 + L**2 * 0.002 + model.variance(log_ar=L)
        assert v_t == pytest.approx(0.09)

    def test_empty_when_the_transfer_never_wins(self):
        band = solve_band(_model(0.05, 0.04), gage_variance=0.05, regression_variance=0.08)
        assert band.empty and band.half_width_log == 0.0 and band.upper == 1.0

    def test_unbounded_when_nothing_grows_with_ratio(self):
        band = solve_band(_model(0.02, 0.0), gage_variance=0.01, regression_variance=0.08)
        assert band.unbounded and band.upper == np.inf and band.lower == 0.0

    def test_other_features_shift_the_band(self):
        model = DissimilarityModel(
            terms=("log_ar", "d_precip"),
            coefficients={"const": 0.02, "log_ar": 0.04, "d_precip": 0.001},
            bound_terms=(),
            n_pairs=10,
            n_gages=8,
            iterations=1,
            converged=True,
            residuals=pd.DataFrame(),
            by_ratio=pd.DataFrame(),
        )
        wet = solve_band(model, gage_variance=0.01, regression_variance=0.09, d_precip=5.0)
        dry = solve_band(model, gage_variance=0.01, regression_variance=0.09, d_precip=0.0)
        assert wet.upper < dry.upper
        with pytest.raises(KeyError, match="d_precip"):
            solve_band(model, gage_variance=0.01, regression_variance=0.09)


class TestBootstrapBand:
    def test_interval_covers_the_true_band(self):
        gages, conns, cov = make_network(300, 4, s0=0.02, s1=0.04, seed=21)
        pairs = build_pairs(gages, conns, EXPONENTS, 0.01, sampling_cov=cov)
        result = bootstrap_band(
            pairs,
            terms=["log_ar"],
            gage_variance=0.01,
            regression_variance=0.09,
            n_boot=200,
            seed=3,
        )
        true_upper = 10.0 ** np.sqrt((0.09 - 0.01 - 0.02) / 0.04)
        lo, hi = result.interval
        assert lo <= true_upper <= hi
        assert lo <= result.estimate.upper <= hi
        assert result.n_boot == 200
        assert result.n_failed == 0

    def test_resamples_gages_not_pairs(self):
        """A pair is kept m_i * m_j times; one drawn gage alone keeps none of its pairs."""
        from flowfreq.transfer_variance import _resample_pairs_by_gage

        pairs = pd.DataFrame(
            {"donor": ["a", "a", "b"], "receiver": ["b", "c", "c"], "x": [1, 2, 3]}
        )
        out = _resample_pairs_by_gage(pairs, {"a": 2, "b": 3})
        assert sorted(out["x"].tolist()) == [1] * 6

    def test_reproducible_with_a_seed(self):
        gages, conns, cov = make_network(60, 3, s0=0.02, s1=0.04, seed=22)
        pairs = build_pairs(gages, conns, EXPONENTS, 0.01, sampling_cov=cov)
        kw = dict(terms=["log_ar"], gage_variance=0.01, regression_variance=0.09, n_boot=30)
        a = bootstrap_band(pairs, seed=5, **kw)
        b = bootstrap_band(pairs, seed=5, **kw)
        assert a.interval == b.interval


# --------------------------------------------------------------------------- leave one out


@pytest.fixture(scope="module")
def scored() -> pd.DataFrame:
    gages, conns, cov = make_network(
        150, 3, s0=0.01, s1=0.02, sampling_var=0.004, regression_var=0.06, seed=31
    )
    return leave_one_out(gages, conns, EXPONENTS, 0.01, terms=["log_ar"], sampling_cov=cov)


class TestLeaveOneOut:

    def test_one_row_per_connected_gage(self, scored):
        assert len(scored) == 450
        assert (scored["n_donors"] == 2).all()
        assert {
            "site_no",
            "observed",
            "blend",
            "blend_variance",
            "blend_error",
            "blend_covered",
            "blend_covered_naive",
            "transfer_error",
            "regression_error",
            "regression_covered",
        } <= set(scored.columns)

    def test_coverage_is_near_nominal(self, scored):
        assert scored["blend_covered"].mean() == pytest.approx(0.90, abs=0.05)
        assert scored["regression_covered"].mean() == pytest.approx(0.90, abs=0.05)

    def test_blend_beats_regression_when_the_regression_is_noisy(self, scored):
        assert scored["blend_error"].abs().mean() < scored["regression_error"].abs().mean()

    def test_without_regression_columns_scores_transfers_only(self):
        gages, conns, cov = make_network(40, 3, s0=0.01, s1=0.02, seed=32)
        gages = gages.drop(columns=["regression_log_q", "regression_variance"])
        out = leave_one_out(gages, conns, EXPONENTS, 0.01, terms=["log_ar"], sampling_cov=cov)
        assert out["regression_error"].isna().all()
        np.testing.assert_allclose(out["blend"], out["transfer_only"])


def test_leave_one_out_uses_a_fixed_donor_correlation():
    gages, conns, cov = make_network(30, 3, s0=0.01, s1=0.02, seed=33)
    out = leave_one_out(
        gages, conns, EXPONENTS, 0.01, terms=["log_ar"], sampling_cov=cov, donor_correlation=0.25
    )
    assert (out["donor_correlation"] == 0.25).all()


def _two_networks() -> Tuple[pd.DataFrame, List[Tuple[str, str]], pd.DataFrame]:
    """Chains of three gages (two donors each) beside chains of five (four donors each)."""
    small, c_small, v_small = make_network(20, 3, s0=0.01, s1=0.02, seed=37)
    big, c_big, v_big = make_network(20, 5, s0=0.01, s1=0.02, seed=38)
    big["site_no"] = "x" + big["site_no"]
    c_big = [("x" + a, "x" + b) for a, b in c_big]
    v_big = v_big.rename(index=lambda s: "x" + s, columns=lambda s: "x" + s)
    cov = pd.concat([v_small, v_big]).fillna(0.0)
    return pd.concat([small, big], ignore_index=True), c_small + c_big, cov[cov.index]


class TestNegativeDonorCorrelation:
    """#177: a negative donor correlation with several donors is not a covariance."""

    def _transfers(self, k: int) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "donor": [f"d{i}" for i in range(k)],
                "transfer": 3.0 + 0.01 * np.arange(k),
                "donor_variance": 0.01,
                "log_ar": 0.1 * (np.arange(k) + 1),
            }
        )

    def test_combine_names_the_bound(self):
        # Four donors allow rho > -1/3, whatever the sampling covariance adds.
        with pytest.raises(ValueError, match=r"-1/\(k - 1\) = -0.3333 for 4 donors"):
            combine(_model(0.02, 0.04), self._transfers(4), None, donor_correlation=-0.4)

    def test_combine_accepts_a_negative_correlation_inside_the_bound(self):
        out = combine(_model(0.02, 0.04), self._transfers(4), None, donor_correlation=-0.3)
        assert np.all(np.linalg.eigvalsh(out.covariance.to_numpy()) > 0)

    def test_two_donors_take_any_correlation_above_minus_one(self):
        out = combine(_model(0.02, 0.04), self._transfers(2), None, donor_correlation=-0.9)
        assert out.weights.sum() == pytest.approx(1.0)

    def test_the_bound_ignores_donors_without_dissimilarity_variance(self):
        # With s2_d = 0 the correlation multiplies nothing.
        out = combine(_model(0.0, 0.0), self._transfers(4), None, donor_correlation=-0.9)
        assert out.weights.sum() == pytest.approx(1.0)

    def test_leave_one_out_floors_a_negative_fold_estimate(self, monkeypatch):
        gages, conns, cov = make_network(30, 5, s0=0.01, s1=0.02, seed=34)
        monkeypatch.setattr(tv, "_rho_from_triples", lambda tri, model, mask: (-0.6, 10))
        out = leave_one_out(gages, conns, EXPONENTS, 0.01, terms=["log_ar"], sampling_cov=cov)
        assert len(out) == 150 and out["note"].eq("").all()
        assert (out["donor_correlation"] == 0.0).all()
        assert (out["donor_correlation_estimate"] == -0.6).all()
        assert out["donor_correlation_clipped"].all()

    def test_min_donor_correlation_none_keeps_the_estimate(self, monkeypatch):
        gages, conns, cov = make_network(30, 3, s0=0.01, s1=0.02, seed=35)
        monkeypatch.setattr(tv, "_rho_from_triples", lambda tri, model, mask: (-0.2, 10))
        out = leave_one_out(
            gages,
            conns,
            EXPONENTS,
            0.01,
            terms=["log_ar"],
            sampling_cov=cov,
            min_donor_correlation=None,
        )
        assert (out["donor_correlation"] == -0.2).all()
        assert not out["donor_correlation_clipped"].any()

    def test_a_failing_gage_is_recorded_and_the_run_continues(self):
        # A fixed rho of -0.5 is a covariance for two donors but not for four.
        gages, conns, cov = _two_networks()
        out = leave_one_out(
            gages,
            conns,
            EXPONENTS,
            0.01,
            terms=["log_ar"],
            sampling_cov=cov,
            donor_correlation=-0.5,
        )
        assert len(out) == 160
        four = out["n_donors"] == 4
        assert out.loc[four, "note"].str.contains(r"-1/\(k - 1\)").all()
        assert out.loc[four, ["blend", "transfer_only", "blend_covered"]].isna().all().all()
        assert out.loc[~four, "note"].eq("").all() and out.loc[~four, "blend"].notna().all()

    def test_a_fixed_correlation_is_not_floored(self):
        gages, conns, cov = make_network(30, 3, s0=0.01, s1=0.02, seed=36)
        out = leave_one_out(
            gages,
            conns,
            EXPONENTS,
            0.01,
            terms=["log_ar"],
            sampling_cov=cov,
            donor_correlation=-0.2,
        )
        assert (out["donor_correlation"] == -0.2).all()
        assert out["donor_correlation_estimate"].isna().all()
        assert not out["donor_correlation_clipped"].any()


def test_combine_refuses_a_covariance_that_is_not_positive_definite():
    transfers = pd.DataFrame(
        {
            "donor": ["a", "b"],
            "transfer": [3.0, 3.1],
            "donor_variance": [0.01, 0.01],
            "log_ar": [0.1, -0.1],
        }
    )
    bad = pd.DataFrame([[0.01, 0.5], [0.5, 0.01]], index=["a", "b"], columns=["a", "b"])
    with pytest.raises(ValueError, match="positive definite"):
        combine(_model(0.0, 0.0), transfers, regression=None, donor_cov=bad)


# --------------------------------------------------------------------------- powers (#181)


def make_pairs(
    n: int, *, s0: float, s1: float, power: float, known: float = 0.01, seed: int = 0
) -> pd.DataFrame:
    """A pair table drawn directly from ``e ~ N(0, known + s0 + s1 |L|^power)``."""
    rng = np.random.default_rng(seed)
    L = rng.uniform(-3.5, 3.5, n)
    e = rng.normal(0.0, np.sqrt(known + s0 + s1 * np.abs(L) ** power))
    return pd.DataFrame(
        {
            "donor": [f"d{k}" for k in range(n)],
            "receiver": [f"r{k}" for k in range(n)],
            "log_ar": L,
            "error": e,
            "known_variance": known,
        }
    )


def _model_p(s0: float, s1: float, power: float) -> DissimilarityModel:
    return DissimilarityModel(
        terms=("log_ar",),
        coefficients={"const": s0, "log_ar": s1},
        bound_terms=(),
        n_pairs=100,
        n_gages=80,
        iterations=1,
        converged=True,
        residuals=pd.DataFrame(),
        by_ratio=pd.DataFrame(),
        powers={"log_ar": power},
    )


class TestPowers:
    def test_default_power_is_two(self):
        assert _model(0.01, 0.04).power("log_ar") == 2.0
        model = fit_dissimilarity(make_pairs(200, s0=0.02, s1=0.03, power=2.0))
        assert model.powers == {"log_ar": 2.0}

    def test_variance_applies_the_power_to_the_absolute_value(self):
        assert _model_p(0.01, 0.04, 1.0).variance(log_ar=-0.5) == pytest.approx(0.01 + 0.04 * 0.5)
        assert _model_p(0.01, 0.04, 1.5).variance(log_ar=-0.25) == pytest.approx(
            0.01 + 0.04 * 0.25**1.5
        )

    def test_recovers_a_linear_term(self):
        pairs = make_pairs(6000, s0=0.02, s1=0.05, power=1.0, seed=41)
        model = fit_dissimilarity(pairs, powers={"log_ar": 1.0})
        assert model.powers == {"log_ar": 1.0}
        assert model.coefficients["const"] == pytest.approx(0.02, abs=0.006)
        assert model.coefficients["log_ar"] == pytest.approx(0.05, abs=0.008)

    def test_power_two_is_the_default_fit(self):
        pairs = make_pairs(300, s0=0.02, s1=0.03, power=2.0, seed=42)
        a = fit_dissimilarity(pairs)
        b = fit_dissimilarity(pairs, powers={"log_ar": 2.0})
        assert a.coefficients == b.coefficients

    def test_bad_powers_raise(self):
        pairs = make_pairs(50, s0=0.02, s1=0.03, power=2.0)
        with pytest.raises(ValueError, match="> 0"):
            fit_dissimilarity(pairs, powers={"log_ar": 0.0})
        with pytest.raises(KeyError, match="d_precip"):
            fit_dissimilarity(pairs, powers={"d_precip": 1.0})


class TestProfilePower:
    @pytest.mark.parametrize("true_power", [1.0, 2.0])
    def test_picks_the_generating_power(self, true_power):
        pairs = make_pairs(4000, s0=0.02, s1=0.05, power=true_power, seed=43)
        prof = profile_power(pairs, grid=(0.5, 1.0, 1.5, 2.0, 3.0))
        assert isinstance(prof, PowerProfile)
        assert prof.best_power == true_power
        assert prof.best.powers == {"log_ar": true_power}
        assert list(prof.table["power"]) == [0.5, 1.0, 1.5, 2.0, 3.0]
        assert {"log_likelihood", "converged", "const", "log_ar"} <= set(prof.table.columns)
        assert set(prof.models) == {0.5, 1.0, 1.5, 2.0, 3.0}
        assert prof.models[1.5].by_ratio["n_pairs"].sum() == len(pairs)

    def test_log_likelihood_is_gaussian(self):
        pairs = make_pairs(100, s0=0.02, s1=0.05, power=1.0, seed=44)
        prof = profile_power(pairs, grid=(1.0,))
        expected = prof.models[1.0].residuals["expected"].to_numpy()
        e = pairs["error"].to_numpy()
        ll = -0.5 * float(np.sum(np.log(2.0 * np.pi * expected) + e**2 / expected))
        assert prof.table["log_likelihood"].iloc[0] == pytest.approx(ll)

    def test_term_and_grid_are_checked(self):
        pairs = make_pairs(50, s0=0.02, s1=0.03, power=2.0)
        with pytest.raises(ValueError, match="d_precip"):
            profile_power(pairs, term="d_precip")
        with pytest.raises(ValueError, match="grid"):
            profile_power(pairs, grid=())


class TestSolveBandPowers:
    @staticmethod
    def _check_edge(model: DissimilarityModel, L: float, vb: float) -> None:
        v_t = 0.01 + L**2 * vb + model.variance(log_ar=L)
        assert v_t == pytest.approx(0.09, rel=1e-9)

    def test_linear_term_has_a_closed_form(self):
        model = _model_p(0.02, 0.04, 1.0)
        band = solve_band(
            model, gage_variance=0.01, regression_variance=0.09, exponent_variance=0.002
        )
        gap = 0.09 - 0.01 - 0.02
        L = (-0.04 + np.sqrt(0.04**2 + 4 * 0.002 * gap)) / (2 * 0.002)
        assert band.half_width_log == pytest.approx(L)
        self._check_edge(model, band.half_width_log, 0.002)

    def test_linear_term_without_exponent_variance(self):
        band = solve_band(_model_p(0.02, 0.04, 1.0), gage_variance=0.01, regression_variance=0.09)
        assert band.half_width_log == pytest.approx(0.06 / 0.04)

    @pytest.mark.parametrize("power", [0.5, 1.5, 3.0])
    def test_other_powers_by_root(self, power):
        model = _model_p(0.02, 0.04, power)
        band = solve_band(
            model, gage_variance=0.01, regression_variance=0.09, exponent_variance=0.002
        )
        assert 0.0 < band.half_width_log < np.inf
        self._check_edge(model, band.half_width_log, 0.002)

    def test_bootstrap_uses_the_powers(self):
        pairs = make_pairs(400, s0=0.02, s1=0.04, power=1.0, seed=45)
        result = bootstrap_band(
            pairs,
            ["log_ar"],
            gage_variance=0.01,
            regression_variance=0.09,
            n_boot=20,
            seed=1,
            powers={"log_ar": 1.0},
        )
        direct = solve_band(
            fit_dissimilarity(pairs, powers={"log_ar": 1.0}),
            gage_variance=0.01,
            regression_variance=0.09,
        )
        assert result.estimate == direct


def test_donor_correlation_uses_the_powers():
    # c receives from a (log ratio +1) and b (log ratio -2): one triple.
    gages = pd.DataFrame(
        {
            "site_no": ["a", "b", "c"],
            "drainage_area_sqmi": [10.0, 10000.0, 100.0],
            "log_q": [2.0, 3.4, 2.7],
            "variance": [0.01, 0.01, 0.01],
        }
    )
    conns = [("a", "c"), ("b", "c")]
    product = (2.0 + B * 1.0 - 2.7) * (3.4 + B * -2.0 - 2.7) - 0.01
    for power in (1.0, 2.0):
        rho, n = estimate_donor_correlation(
            gages, conns, EXPONENTS, 0.01, _model_p(0.1, 0.4, power)
        )
        s_a, s_b = 0.1 + 0.4, 0.1 + 0.4 * 2.0**power
        assert n == 1
        assert rho == pytest.approx(product / np.sqrt(s_a * s_b))


def test_leave_one_out_refits_with_the_powers():
    gages, conns, cov = make_network(40, 3, s0=0.01, s1=0.02, seed=34)
    default = leave_one_out(gages, conns, EXPONENTS, 0.01, sampling_cov=cov)
    squared = leave_one_out(gages, conns, EXPONENTS, 0.01, sampling_cov=cov, powers={"log_ar": 2.0})
    linear = leave_one_out(gages, conns, EXPONENTS, 0.01, sampling_cov=cov, powers={"log_ar": 1.0})
    pd.testing.assert_frame_equal(default, squared)
    assert not np.allclose(default["blend"], linear["blend"])
