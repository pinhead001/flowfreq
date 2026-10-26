"""``bootstrap_band(profile=...)``: re-choose the dissimilarity power in each replicate.

``bootstrap_band`` held a ``powers=`` choice fixed, so its interval left out the
uncertainty in choosing the power (#181 follow-up; TODO "Priorities for methow-flow"
item 12). With ``profile=<term>`` each replicate runs :func:`profile_power` on its own
resampled pairs and solves the band at the winning power; the result reports how often
each grid power won.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from flowfreq.transfer_variance import BandBootstrap, bootstrap_band, profile_power, solve_band

KW = dict(gage_variance=0.01, regression_variance=0.09)


def make_pairs(n: int, *, power: float, seed: int, s0: float = 0.02, s1: float = 0.05):
    rng = np.random.default_rng(seed)
    L = rng.uniform(-3.5, 3.5, n)
    e = rng.normal(0.0, np.sqrt(0.01 + s0 + s1 * np.abs(L) ** power))
    return pd.DataFrame(
        {
            "donor": [f"d{k}" for k in range(n)],
            "receiver": [f"r{k}" for k in range(n)],
            "log_ar": L,
            "error": e,
            "known_variance": 0.01,
        }
    )


def test_default_holds_the_power_and_reports_no_counts() -> None:
    res = bootstrap_band(make_pairs(300, power=2.0, seed=1), ["log_ar"], n_boot=5, seed=2, **KW)
    assert isinstance(res, BandBootstrap)
    assert res.power_counts is None and res.estimate_power is None


def test_profile_rechooses_the_power_per_replicate() -> None:
    pairs = make_pairs(2000, power=1.0, seed=3)
    res = bootstrap_band(
        pairs, ["log_ar"], n_boot=20, seed=4, profile="log_ar", grid=(1.0, 2.0), **KW
    )
    assert res.power_counts is not None
    assert set(res.power_counts) == {1.0, 2.0}
    assert sum(res.power_counts.values()) == res.n_boot - res.n_failed == len(res.replicates)
    assert res.power_counts[1.0] >= 15  # the generating power wins most replicates
    full = profile_power(pairs, grid=(1.0, 2.0))
    assert res.estimate_power == full.best_power == 1.0
    assert res.estimate == solve_band(full.best, **KW)


def test_profile_is_reproducible_with_a_seed() -> None:
    pairs = make_pairs(400, power=1.5, seed=5)
    kw = dict(n_boot=10, seed=6, profile="log_ar", grid=(1.0, 1.5, 2.0), **KW)
    a = bootstrap_band(pairs, ["log_ar"], **kw)
    b = bootstrap_band(pairs, ["log_ar"], **kw)
    assert a.power_counts == b.power_counts
    np.testing.assert_array_equal(a.replicates, b.replicates)


def test_profile_argument_checks() -> None:
    pairs = make_pairs(100, power=2.0, seed=7)
    with pytest.raises(ValueError, match="d_precip"):
        bootstrap_band(pairs, ["log_ar"], n_boot=2, profile="d_precip", **KW)
    with pytest.raises(ValueError, match="log_ar"):  # fixed and profiled at once
        bootstrap_band(pairs, ["log_ar"], n_boot=2, profile="log_ar", powers={"log_ar": 1.0}, **KW)
    with pytest.raises(ValueError, match="grid"):
        bootstrap_band(pairs, ["log_ar"], n_boot=2, profile="log_ar", grid=(), **KW)
