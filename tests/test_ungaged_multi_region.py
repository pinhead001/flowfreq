"""``region=`` as a list of codes: one equation set per side of a boundary basin (#203).

When a basin straddles a region boundary and each side has its own overlapping sets
(superseded and current), one code covers only its side and is refused. A list of codes,
one per side, whose ``percent_weight`` together covers the basin (100 % within NSS's
rounding tolerance), is used as a partition: area-weighted across the chosen regions,
with every other region reported in ``alternatives``.
"""

from __future__ import annotations

import math

import pytest
from click.testing import CliRunner

from flowfreq.cli import cli
from flowfreq.regression.library import EquationsUnavailable
from flowfreq.ungaged import OverlappingRegionsError, assemble

from .test_ungaged import KW, LIB, _flow, _region

CHARS = {"DRNAREA": 100.0, "PRECIP": 30.0}


def _boundary():
    # GC1 60 % / GC2 40 % partition the basin; GC4 is a set that covers all of it.
    return [_region("GC1", 60.0), _region("GC2", 40.0), _region("GC4", 100.0)]


def test_list_covering_the_basin_is_used_as_a_partition() -> None:
    regions = _boundary()
    a = assemble(LIB, regions, CHARS, region=["GC1", "GC2"], **KW)
    assert [r.code for r in a.used] == ["GC1", "GC2"]
    assert a.alternatives == [regions[2]]
    assert [f for f, _ in a.estimates[0.1]] == pytest.approx([0.6, 0.4])
    row = a.quantiles.set_index("aep").loc[0.1]
    expected = 10 ** (0.6 * math.log10(_flow("GC1", 0.1)) + 0.4 * math.log10(_flow("R2", 0.1)))
    assert row["flow_cfs"] == pytest.approx(expected)
    assert row["interval_method"] == "none: area-weighted across 2 regions"
    assert row["regions"] == "GC1, GC2"


def test_list_accepts_library_codes_case_insensitively() -> None:
    a = assemble(LIB, _boundary(), CHARS, region=["gc1", "r2"], **KW)
    assert [r.code for r in a.used] == ["GC1", "GC2"]


def test_single_element_list_behaves_like_a_string() -> None:
    a = assemble(LIB, _boundary(), CHARS, region=["GC4"], **KW)
    assert [r.code for r in a.used] == ["GC4"]
    assert [r.code for r in a.alternatives] == ["GC1", "GC2"]
    with pytest.raises(OverlappingRegionsError, match="only 60 % of the basin"):
        assemble(LIB, _boundary(), CHARS, region=["GC1"], **KW)


def test_list_covering_part_of_the_basin_is_refused() -> None:
    regions = [_region("GC1", 60.0), _region("GC2", 40.0), _region("GC4", 100.0)]
    with pytest.raises(OverlappingRegionsError, match="40 %") as info:
        assemble(LIB, regions, CHARS, region=["GC2"], **KW)
    assert [c[0] for c in info.value.candidates] == ["GC1", "GC2", "GC4"]


def test_list_whose_regions_overlap_each_other_is_refused() -> None:
    with pytest.raises(OverlappingRegionsError, match="200 %"):
        assemble(LIB, _boundary(), CHARS, region=["GC1", "GC2", "GC4"], **KW)


def test_list_with_a_code_not_in_use_is_refused() -> None:
    with pytest.raises(OverlappingRegionsError, match="GC9"):
        assemble(LIB, _boundary(), CHARS, region=["GC1", "GC9"], **KW)


@pytest.mark.parametrize("bad", [[], ["GC1", "gc1"], ["GC1", "GC1"]])
def test_empty_or_duplicate_list_is_a_value_error(bad: list) -> None:
    with pytest.raises(ValueError):
        assemble(LIB, _boundary(), CHARS, region=bad, **KW)


def test_list_on_a_partition_is_a_no_op_when_every_code_is_in_use() -> None:
    regions = [_region("GC1", 60.0), _region("GC2", 40.0)]
    a = assemble(LIB, regions, CHARS, region=["GC2", "GC1"], **KW)
    assert [r.code for r in a.used] == ["GC1", "GC2"] and a.alternatives == []
    with pytest.raises(EquationsUnavailable, match="GC7"):
        assemble(LIB, regions, CHARS, region=["GC1", "GC7"], **KW)


def test_cli_region_accepts_a_comma_separated_list(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: dict = {}

    class _Est:
        def to_markdown(self) -> str:
            return "SYNTHETIC-MARKDOWN"

    def fake(lat: float, lon: float, **kw: object) -> _Est:
        calls.update(kw)
        return _Est()

    monkeypatch.setattr("flowfreq.ungaged.estimate_at_site", fake)
    base = ["ungaged", "--lat", "36.3", "--lon", "-87.0", "--state", "TN"]
    result = CliRunner().invoke(cli, base + ["--region", "GC1, GC2"])
    assert result.exit_code == 0, result.output
    assert calls["region"] == ["GC1", "GC2"]
    result = CliRunner().invoke(cli, base + ["--region", "GC2016"])
    assert calls["region"] == "GC2016"
