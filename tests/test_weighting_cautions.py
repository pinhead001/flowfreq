"""``StateLibrary.weighting_cautions``: a report's limit on area-weighting across regions.

Montana's SIR 2015-5019-F (p. 13) limits weighting a Northwest basin that runs into the
Northwest Foothills to AEPs greater than 4 percent.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from flowfreq import estimate_from_characteristics
from flowfreq.regression.library import WeightingCaution, load_state, parse_state

CHARS = {"CONTDA": 100.0, "PRECIP": 15.0}
AT = {"lat": 48.0, "lon": -112.5}


def test_montana_cautions_northwest_into_the_foothills_at_4_percent_and_rarer():
    est = estimate_from_characteristics(
        "MT", CHARS, region={"NW": 40.0, "NWF": 60.0}, aeps=[0.5, 0.1, 0.04, 0.01], **AT
    )
    assert "appropriate only for AEPs greater than 4 percent" in est.advisory
    assert est.advisory.endswith("AEPs this reaches here, in percent: 4, 1.")
    assert est.provenance["weighting_cautions"] == [
        {
            "regions": ["NW", "NWF"],
            "aeps": [0.04, 0.01],
            "source": "SIR 2015-5019-F ver. 1.1, p. 13",
        }
    ]
    # A caution, not a refusal: every AEP is still weighted.
    assert est.quantiles["flow_cfs"].notna().all()
    assert "greater than 4 percent" in est.to_markdown()


def test_no_caution_above_4_percent():
    est = estimate_from_characteristics(
        "MT", CHARS, region={"NW": 40.0, "NWF": 60.0}, aeps=[0.5, 0.1], **AT
    )
    assert "greater than 4 percent" not in est.advisory
    assert "weighting_cautions" not in est.provenance


@pytest.mark.parametrize(
    "region",
    ["NW", "NWF", {"NWF": 50.0, "NEP": 50.0}],
    ids=["northwest", "foothills", "foothills-and-plains"],
)
def test_no_caution_for_one_region_or_another_pair(region):
    chars = {**CHARS, "EL5000": 10.0, "FOREST": 10.0, "SLOP30_30M": 5.0, "ET0306MOD": 2.0}
    est = estimate_from_characteristics(
        "MT", chars, region=region, aeps=[0.01], allow_extrapolation=True, **AT
    )
    assert "greater than 4 percent" not in est.advisory
    assert "weighting_cautions" not in est.provenance


def test_affected_needs_every_region_and_keeps_the_limit():
    wc = WeightingCaution(regions=("A", "B"), text="T.", aep_at_or_below=0.04)
    assert wc.affected({"A", "B", "C"}, [0.01, 0.5, 0.04]) == [0.04, 0.01]
    assert wc.affected({"A"}, [0.01]) == []
    every = WeightingCaution(regions=("A", "B"), text="T.")
    assert every.affected({"A", "B"}, [0.01, 0.5]) == [0.5, 0.01]
    assert wc.caution([0.04, 0.002]) == "T. AEPs this reaches here, in percent: 4, 0.2."


def test_unknown_or_single_regions_are_refused():
    data = json.loads(
        (Path(__file__).parent.parent / "flowfreq/data/regression/MT.json").read_text("utf-8")
    )
    bad = dict(data, weighting_cautions=[{"regions": ["NW", "NOPE"], "text": "x"}])
    with pytest.raises(ValueError, match="weighting_cautions"):
        parse_state(bad)
    bad = dict(data, weighting_cautions=[{"regions": ["NW"], "text": "x"}])
    with pytest.raises(ValueError, match="two or more regions"):
        parse_state(bad)


def test_only_montana_has_one():
    files = sorted((Path(__file__).parent.parent / "flowfreq/data/regression").glob("*.json"))
    have = [f.stem for f in files if json.loads(f.read_text("utf-8")).get("weighting_cautions")]
    assert len(files) > 40 and have == ["MT"]


@pytest.mark.parametrize("width", ["BF", "AC", "RS"])
def test_montanas_channel_width_sets_carry_the_same_limit(width):
    """SIR 2020-5142 (pp. 38 and 47) repeats it for the sets it adds."""
    lib = load_state("MT")
    pair = (f"NW-{width}", f"NWF-{width}")
    (found,) = [c for c in lib.weighting_cautions if set(c.regions) == set(pair)]
    assert found.aep_at_or_below == 0.04 and "SIR 2020-5142" in found.text
    assert found.affected(pair, [0.5, 0.1, 0.04, 0.01]) == [0.04, 0.01]
    assert found.affected(pair[:1], [0.01]) == []
