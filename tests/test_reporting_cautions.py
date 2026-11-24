"""A regression's cautions reach every report table its flows are in, and a chat's reply.

Nothing here makes a request: StreamStats is faked as in ``test_ungaged``.
"""

from __future__ import annotations

from typing import Any, Dict, List

import pytest

from flowfreq.regression.library import load_state
from flowfreq.reporting import Paragraph, ReportError, build_section
from flowfreq.reporting.blocks import Bullets, Table
from flowfreq.reporting.session import ReportSession
from tests.test_reporting_daily import FACTORS
from tests.test_reporting_sites import GAGE, TENNESSEE
from tests.test_reporting_summary import LOCATIONS
from tests.test_ungaged import Fakes, _wa_site
from tests.test_ungaged_named_region import SUPPLIED, _other

REGION = "Urban 2026"
URBAN: Dict[str, Any] = {
    "state": "TN",
    "characteristics": {
        "DRNAREA": 20.0,
        "LC11DEV": 60.0,
        "LC11IMP": 30.0,
        "L3_PIEDMNT": 0.0,
        "L3_RDGVLY": 0.0,
    },
    "region": REGION,
    "lat": 36.04,
    "lon": -88.23,
}
CAUTION = load_state("TN").region_cautions[REGION]
ESTIMATES: List[Dict[str, Any]] = [
    {"label": "Gage", "kind": "gaged_ffa", "options": dict(GAGE)},
    {"label": "Urban regression", "kind": "ungaged", "options": dict(URBAN)},
]


def _options(kind: str, detail: bool, regression: Dict[str, Any]) -> Dict[str, Any]:
    if kind == "weighted":
        return {"gage": dict(GAGE), "regression": dict(regression), "detail": detail}
    if kind == "future_flow":
        return {"regression": dict(regression), "factor_set": dict(FACTORS), "detail": detail}
    estimates = [ESTIMATES[0], {**ESTIMATES[1], "options": dict(regression)}]
    if kind == "site_summary":
        return {"locations": LOCATIONS, "estimates": estimates}
    return {"estimates": estimates, "detail": detail}


def _all_text(blocks: List[Any]) -> str:
    out = []
    for b in blocks:
        out.extend(b.items if isinstance(b, Bullets) else [getattr(b, "text", "")])
    return "\n".join(out)


KINDS = [
    ("weighted", True),
    ("weighted", False),
    ("comparison", True),
    ("comparison", False),
    ("future_flow", True),
    ("future_flow", False),
    ("site_summary", True),
]


@pytest.mark.parametrize("kind,detail", KINDS)
def test_every_kind_prints_the_caution_once(kind: str, detail: bool) -> None:
    built = build_section(kind, _options(kind, detail, URBAN))
    assert _all_text(built.blocks).count(CAUTION) == 1
    assert len(built.cautions) == 1 and CAUTION in built.cautions[0]


@pytest.mark.parametrize("kind,detail", KINDS)
def test_no_caution_no_words(kind: str, detail: bool) -> None:
    built = build_section(kind, _options(kind, detail, TENNESSEE))
    assert built.cautions == [] and "Cautions" not in _all_text(built.blocks)


@pytest.mark.parametrize("kind", ["comparison", "site_summary"])
def test_the_caution_names_its_estimate_and_follows_the_flows(kind: str) -> None:
    blocks = build_section(kind, _options(kind, False, URBAN)).blocks
    (at,) = [i for i, b in enumerate(blocks) if isinstance(b, Paragraph) and CAUTION in b.text]
    assert blocks[at].text.startswith("Cautions on Urban regression: ")
    assert any(isinstance(b, Table) for b in blocks[:at])


def test_the_ungaged_section_itself() -> None:
    built = build_section("ungaged", dict(URBAN))
    assert _all_text(built.blocks).count(CAUTION) == 1 and built.cautions == [CAUTION]


class TestSession:
    def test_replies_carry_the_cautions(self, tmp_path: Any) -> None:
        session = ReportSession("Urban", output_dir=tmp_path)
        assert "cautions" not in session.add_section("gaged_ffa", dict(GAGE))
        assert session.add_section("ungaged", dict(URBAN))["cautions"] == [CAUTION]
        for reply in (
            session.weight("s1", "s2", keep_parts=True),
            session.compare(["s1", "s2"], labels=["Gage", "Urban"]),
            session.summarize(["s1", "s2"], LOCATIONS, labels=["Gage", "Urban"]),
        ):
            (said,) = reply["cautions"]
            assert CAUTION in said


# -- an ungaged section through StreamStats with values the analyst supplies -----------------

TEXAS: Dict[str, Any] = {
    "state": "TX",
    "lat": 29.43,
    "lon": -100.0,
    "region": "Statewide OmegaEM",
    "supplied_characteristics": dict(SUPPLIED),
}


@pytest.fixture
def texas(monkeypatch: pytest.MonkeyPatch) -> Fakes:
    site = {
        **_wa_site(),
        "streamstats_region": "TX",
        "characteristics": {"DRNAREA": 760.0, "PRECIPCENT": 24.5, "OmegaEM": 0.33},
    }
    f = Fakes(site)
    f.regions = [_other("TX")]
    monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", f.delineate)
    monkeypatch.setattr("flowfreq.ungaged.locate_regression_regions", f.locate)
    return f


class TestSupplied:
    def test_stored_with_streamstats_values_and_named_as_the_analysts(self, texas: Fakes) -> None:
        built = build_section("ungaged", dict(TEXAS))
        options = built.options
        assert options["region"] == {"Statewide OmegaEM": 100.0}
        assert options["characteristics"]["MCSLOPE"] == 0.00326
        assert "supplied_characteristics" not in options
        assert options["streamstats"]["supplied"] == ["CONTDA", "MCSLOPE"]
        opening = built.blocks[0].text
        assert "from USGS StreamStats" in opening
        assert "CONTDA and MCSLOPE were supplied by the analyst" in opening
        captions = [b.caption for b in built.blocks if isinstance(b, Table)]
        assert (
            "Basin characteristics the equations use, from StreamStats except CONTDA and MCSLOPE, "
            "supplied by the analyst." in captions
        )
        notes = _all_text(built.blocks)
        assert "NSS does not serve it" in notes and texas.regions[0].code in notes

    @pytest.mark.parametrize(
        "codes,said",
        [
            (["MCSLOPE"], ". MCSLOPE was supplied by the analyst"),
            (["A", "B", "C"], ". A, B and C were supplied by the analyst"),
        ],
    )
    def test_the_words_fit_the_number_of_codes(
        self, texas: Fakes, codes: List[str], said: str
    ) -> None:
        options = build_section("ungaged", dict(TEXAS)).options
        options["streamstats"]["supplied"] = codes
        assert said in build_section("ungaged", options).blocks[0].text

    def test_rebuilds_without_streamstats(self, texas: Fakes) -> None:
        first = build_section("ungaged", dict(TEXAS))
        again = build_section("ungaged", first.options)
        assert len(texas.delineate_calls) == 1
        assert again.options == first.options and again.blocks == first.blocks

    def test_only_without_characteristics(self) -> None:
        with pytest.raises(ReportError, match="supplied_characteristics is for a section without"):
            build_section("ungaged", {**URBAN, "supplied_characteristics": {"LC11DEV": 60.0}})

    @pytest.mark.parametrize("bad", [[1], "MCSLOPE=0.003", {}])
    def test_must_be_an_object(self, texas: Fakes, bad: Any) -> None:
        with pytest.raises(ReportError, match="supplied_characteristics must be an object"):
            build_section("ungaged", {**TEXAS, "supplied_characteristics": bad})
        assert texas.delineate_calls == []

    @pytest.mark.parametrize("bad", ["MCSLOPE", [1]])
    def test_the_stored_record_is_checked(self, texas: Fakes, bad: Any) -> None:
        options = build_section("ungaged", dict(TEXAS)).options
        options["streamstats"]["supplied"] = bad
        with pytest.raises(ReportError, match="streamstats supplied must be a list of codes"):
            build_section("ungaged", options)
