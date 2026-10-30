"""``region=`` names a library region when NSS locates no peak-flow region (#207).

Checked against NSS on 2026-10-30: ``regressionregions/bylocation`` returns none of Wyoming's
six peak-flow regions anywhere in the state (only the SIR 2026-5120 flow-duration regions)
and nothing at all in Nevada, so a caller who reads the report's region map had no way to
use the verified WY or NV equations through ``estimate_at_site``. With ``region=`` naming a
library region (its ``region_code``, or an NSS code a ``region_name`` carries) and no located
region matching the library, that region is used alone and the assembly says the caller chose
it.
"""

from __future__ import annotations

import pytest

from flowfreq.regression.library import EquationsUnavailable
from flowfreq.ungaged import assemble

from .test_ungaged import CHARS, KW, LIB, _flow, _region

CALLER = "caller; NSS located no peak-flow region for this basin"


def _check_single(a, code: str, lib_code: str) -> None:
    assert [r.code for r in a.used] == [code]
    assert a.alternatives == []
    ((fraction, est),) = a.estimates[0.01]
    assert fraction == 1.0
    row = a.quantiles.set_index("aep").loc[0.01]
    assert row["flow_cfs"] == pytest.approx(_flow(lib_code, 0.01))
    assert est.interval is not None
    assert (row["lower"], row["upper"]) == pytest.approx(est.interval)
    assert a.region_choice_source == CALLER


def test_no_located_region_with_an_nss_code_choice_uses_that_region() -> None:
    a = assemble(LIB, [], CHARS, region="GC2", **KW)
    _check_single(a, "GC2", "R2")


def test_no_located_region_with_a_library_code_choice_uses_that_region() -> None:
    a = assemble(LIB, [], CHARS, region="r2", **KW)
    _check_single(a, "R2", "R2")


def test_unmatched_located_regions_with_a_choice_use_the_chosen_region() -> None:
    located = [_region("GC9", 100.0)]
    a = assemble(LIB, located, CHARS, region="GC2", **KW)
    _check_single(a, "GC2", "R2")
    assert a.unmatched == located


def test_unknown_choice_still_raises_naming_the_located_codes() -> None:
    with pytest.raises(EquationsUnavailable, match="GC9") as info:
        assemble(LIB, [_region("GC9", 100.0)], CHARS, region="NOPE", **KW)
    assert "'NOPE' is not a library region" in str(info.value)


def test_unknown_choice_with_nothing_located_still_raises() -> None:
    with pytest.raises(EquationsUnavailable, match="NSS located no regression region") as info:
        assemble(LIB, [], CHARS, region="NOPE", **KW)
    assert "'NOPE' is not a library region" in str(info.value)


def test_without_a_choice_both_refusals_are_unchanged() -> None:
    with pytest.raises(EquationsUnavailable, match="NSS located no regression region"):
        assemble(LIB, [], CHARS, **KW)
    with pytest.raises(EquationsUnavailable, match="GC9"):
        assemble(LIB, [_region("GC9", 100.0)], CHARS, **KW)


def test_region_choice_source_is_none_on_the_normal_path() -> None:
    a = assemble(LIB, [_region("GC1", 100.0)], CHARS, **KW)
    assert a.region_choice_source is None
    b = assemble(LIB, [_region("GC1", 100.0), _region("GC2", 100.0)], CHARS, region="GC2", **KW)
    assert b.region_choice_source is None
