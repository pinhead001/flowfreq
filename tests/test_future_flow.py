"""Tests for flowfreq.future_flow. All factor values are synthetic."""

from __future__ import annotations

import pandas as pd
import pytest

from flowfreq.future_flow import ChangeFactorSet, apply_change_factors, select_factor_set


def _set(**kw):
    base = dict(
        name="synthetic-national",
        aeps=(0.5, 0.1, 0.01),
        factors=(1.05, 1.10, 1.20),
        scenario="test",
        horizon="2050",
        citation="Synthetic test source",
    )
    base.update(kw)
    return ChangeFactorSet(**base)


QUANTILES = pd.DataFrame({"aep": [0.5, 0.1, 0.01, 0.002], "flow_cfs": [100.0, 200.0, 400.0, 600.0]})


def test_apply_keeps_regulatory_and_flags_extrapolation():
    fq = apply_change_factors(QUANTILES, _set())
    t = fq.table
    assert list(t["regulatory_flow_cfs"]) == [100.0, 200.0, 400.0, 600.0]
    assert list(t["change_factor"][:3]) == pytest.approx([1.05, 1.10, 1.20])
    assert t["future_flow_cfs"].iloc[2] == pytest.approx(480.0)
    assert list(t["extrapolated"]) == [False, False, False, True]
    assert t["change_factor"].iloc[3] == pytest.approx(1.20)  # endpoint held
    assert fq.provenance["n_extrapolated"] == 1 and fq.notes


def test_interpolates_between_published_points():
    f, extrap = _set().factor_at([0.05])
    assert 1.10 < f[0] < 1.20 and not extrap[0]


def test_missing_column():
    with pytest.raises(KeyError):
        apply_change_factors(QUANTILES.rename(columns={"flow_cfs": "q"}), _set())


@pytest.mark.parametrize(
    "kw",
    [
        {"citation": ""},
        {"scenario": " "},
        {"factors": (1.0, 1.0)},
        {"factors": (1.0, -1.0, 1.0)},
        {"aeps": (0.5, 0.5, 0.01)},
        {"legal_status": "mandatory"},
        {"level": "state"},  # no geography
    ],
)
def test_invalid_sets(kw):
    with pytest.raises(ValueError):
        _set(**kw)


def test_state_overrides_national_only_where_it_applies():
    nat = _set()
    wa = _set(name="synthetic-wa", level="state", geography=("WA",), legal_status="recommended")
    chosen, overrode = select_factor_set(nat, [wa], state="wa")
    assert chosen is wa and overrode == "synthetic-national"
    chosen, overrode = select_factor_set(nat, [wa], state="OR")
    assert chosen is nat and overrode is None
    fq = apply_change_factors(QUANTILES, chosen, overrode=overrode)
    assert fq.provenance["overrode"] is None


def test_select_rejects_bad_default_and_ambiguity():
    wa1 = _set(name="a", level="state", geography=("WA",))
    wa2 = _set(name="b", level="state", geography=("WA",))
    with pytest.raises(ValueError, match="national"):
        select_factor_set(wa1)
    with pytest.raises(ValueError, match="2 state sets"):
        select_factor_set(_set(), [wa1, wa2], state="WA")
