"""Report sections that take a site: `ungaged` through StreamStats, and `weighted`.

StreamStats is faked at the two names ``estimate_at_site`` calls, as in ``test_ungaged``;
nothing here makes a request.
"""

from __future__ import annotations

import datetime as dt
from typing import Any, Dict, List

import numpy as np
import pytest

from flowfreq import estimate_from_characteristics
from flowfreq.gage_weighting import weight_with_regression
from flowfreq.reporting import (
    SECTION_KINDS,
    Bullets,
    Heading,
    Paragraph,
    ReportError,
    ReportSpec,
    SectionSpec,
    Table,
    build_report,
    build_section,
)
from flowfreq.reporting.sections import BULLETIN_17C, cfs
from flowfreq.workflow import run_ffa
from tests.fixtures import big_sandy as bs
from tests.test_ungaged import WA_LAT, WA_LON, Fakes, _wa_site

YEARS = sorted(bs.SYSTEMATIC_PEAKS)
FLOWS = [bs.SYSTEMATIC_PEAKS[y] for y in YEARS]
SKEW = {"regional_skew": bs.REGIONAL_SKEW, "regional_skew_se": bs.REGIONAL_SKEW_SD}
GAGE: Dict[str, Any] = {
    "site_no": "03606500",
    "site_name": "Big Sandy River at Bruceton, TN",
    "drainage_area_sqmi": 205.0,
    "peaks": {"water_year": YEARS, "peak_flow_cfs": FLOWS},
    "stationarity": False,
    "peaks_table": False,
    **SKEW,
}
#: Tennessee's drainage-area-only equations for hydrologic area 4, at the Big Sandy gage.
TENNESSEE: Dict[str, Any] = {
    "state": "TN",
    "characteristics": {"DRNAREA": 205.0},
    "region": "HA4 DA only",
    "lat": 36.04,
    "lon": -88.23,
}
WASHINGTON: Dict[str, Any] = {"state": "WA", "lat": WA_LAT, "lon": WA_LON}
#: A characteristic StreamStats returns for Washington that no equation there reads.
UNUSED = "ELEVMAX"
WEIGHTED_CAPTION = "At-site and regression flood quantiles and their weighted estimate."


def _table(blocks: List[Any], caption: str) -> Table:
    (table,) = [b for b in blocks if isinstance(b, Table) and b.caption == caption]
    return table


def _text(blocks: List[Any]) -> str:
    return "\n".join(b.text for b in blocks if isinstance(b, Paragraph))


def _caveats(blocks: List[Any]) -> Bullets:
    last = blocks[-1]
    assert isinstance(last, Bullets)
    return last


@pytest.fixture
def fakes(monkeypatch) -> Fakes:
    site = _wa_site()
    site["characteristics"][UNUSED] = 5000.0
    f = Fakes(site)
    monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", f.delineate)
    monkeypatch.setattr("flowfreq.ungaged.locate_regression_regions", f.locate)
    return f


# -- ungaged, from StreamStats ------------------------------------------------------------


class TestUngagedFromStreamStats:
    def test_streamstats_values_are_written_into_the_options(self, fakes):
        content = build_section("ungaged", dict(WASHINGTON))
        options = content.options
        assert len(fakes.delineate_calls) == 1 and len(fakes.locate_calls) == 1
        assert options["characteristics"] == fakes.site["characteristics"]
        assert options["region"] == pytest.approx({"1": 0.92, "4": 99.08})
        record = options["streamstats"]
        assert record["retrieved"] == dt.date.today().isoformat()
        assert record["snapped_lat"] == pytest.approx(WA_LAT + 1e-4)
        assert record["snapped_lon"] == pytest.approx(WA_LON - 1e-4)
        assert record["snap_distance_m"] == 13.2
        assert record["polygon_area_sq_mi"] == 82.2
        assert record["regions_located"] == [
            {"code": "GC1750", "name": "Regression Region 1", "percent": 0.92},
            {"code": "GC1753", "name": "Regression Region 4", "percent": 99.08},
        ]
        assert record["notes"] == [] and "region_choice" not in record

    def test_the_section_is_the_stored_inputs_estimate(self, fakes):
        content = build_section("ungaged", dict(WASHINGTON))
        direct = estimate_from_characteristics(
            "WA",
            content.options["characteristics"],
            region=content.options["region"],
            lat=WA_LAT,
            lon=WA_LON,
        )
        table = _table(content.blocks, "Regression flood quantiles with prediction limits.")
        assert [row[1] for row in table.rows] == [cfs(q) for q in direct.quantiles["flow_cfs"]]
        assert len(table.rows) == len(direct.quantiles) > 0

    def test_the_stored_options_rebuild_without_streamstats(self, fakes):
        first = build_section("ungaged", dict(WASHINGTON))
        again = build_section("ungaged", first.options)
        assert len(fakes.delineate_calls) == 1 and len(fakes.locate_calls) == 1
        assert again.options == first.options
        assert again.blocks == first.blocks and again.references == first.references

    def test_the_text_says_where_the_values_came_from(self, fakes):
        content = build_section("ungaged", dict(WASHINGTON))
        opening = _text(content.blocks)
        assert f"from USGS StreamStats, retrieved {dt.date.today().isoformat()}" in opening
        assert "which it moved to 45.84422, -120.79649 on its stream network (13 m away)" in opening
        assert "the delineated basin is 82.2 square miles" in opening
        assert "does not change if StreamStats later does" in opening
        assert "supplied by the analyst" not in opening

    def test_only_the_characteristics_the_equations_use_are_listed(self, fakes):
        content = build_section("ungaged", dict(WASHINGTON))
        assert content.options["characteristics"][UNUSED] == 5000.0
        table = _table(content.blocks, "Basin characteristics the equations use, from StreamStats.")
        assert [row[0] for row in table.rows] == ["CANOPY_PCT", "DRNAREA", "PRECPRIS10"]

    def test_all_characteristics_lists_every_one(self, fakes):
        content = build_section("ungaged", {**WASHINGTON, "all_characteristics": True})
        (table,) = [
            b for b in content.blocks if isinstance(b, Table) and b.header == ("Code", "Value")
        ]
        assert UNUSED in [row[0] for row in table.rows]
        assert table.caption != "Basin characteristics the equations use, from StreamStats."

    def test_a_streamstats_advisory_is_stored_and_printed(self, fakes, monkeypatch):
        delineate = fakes.delineate

        def advised(*args: Any, **kwargs: Any) -> Any:
            ws = delineate(*args, **kwargs)
            ws.advisories = ["Delineated basin may contain karst"]
            return ws

        monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", advised)
        first = build_section("ungaged", dict(WASHINGTON))
        (note,) = first.options["streamstats"]["notes"]
        assert note.startswith("StreamStats warned about this basin: Delineated basin may contain")
        assert note in _caveats(first.blocks).items
        # ... and the stored note still prints when StreamStats is not asked again
        again = build_section("ungaged", first.options)
        assert note in _caveats(again.blocks).items

    def test_a_region_choice_is_passed_on_and_recorded(self, fakes, monkeypatch):
        import flowfreq.ungaged as ungaged

        seen: Dict[str, Any] = {}
        real = ungaged.estimate_at_site

        def recording(*args: Any, **kwargs: Any) -> Any:
            seen.update(kwargs)
            return real(*args, **{**kwargs, "region": None})

        monkeypatch.setattr(ungaged, "estimate_at_site", recording)
        content = build_section("ungaged", {**WASHINGTON, "region": ["GC1750", "GC1753"]})
        assert seen["region"] == ["GC1750", "GC1753"]
        assert content.options["streamstats"]["region_choice"] == ["GC1750", "GC1753"]
        assert content.options["region"] == pytest.approx({"1": 0.92, "4": 99.08})

    def test_a_region_choice_must_be_a_code_or_a_list(self, fakes):
        with pytest.raises(ReportError, match="region is the choice among overlapping NSS"):
            build_section("ungaged", {**WASHINGTON, "region": {"1": 100.0}})
        assert fakes.delineate_calls == []

    def test_an_estimate_the_stored_inputs_do_not_reproduce_is_refused(self, fakes, monkeypatch):
        import flowfreq.ungaged as ungaged

        real = ungaged.estimate_from_characteristics

        def drifted(*args: Any, **kwargs: Any) -> Any:
            estimate = real(*args, **kwargs)
            estimate.quantiles["flow_cfs"] *= 1.001
            return estimate

        monkeypatch.setattr(ungaged, "estimate_from_characteristics", drifted)
        with pytest.raises(ReportError, match="does not reproduce the StreamStats estimate"):
            build_section("ungaged", dict(WASHINGTON))

    def test_streamstats_errors_pass_through(self, fakes, monkeypatch):
        from flowfreq.streamstats import UnsnappablePointError

        def unsnappable(region: str, lat: float, lon: float, **kw: Any) -> Any:
            raise UnsnappablePointError("no stream within tolerance")

        monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", unsnappable)
        with pytest.raises(UnsnappablePointError, match="no stream within tolerance"):
            build_section("ungaged", dict(WASHINGTON))

    @pytest.mark.parametrize(
        "record, message",
        [
            ("yesterday", "streamstats must be an object"),
            ({"retrieved": "2026-11-15", "server": "prodweba"}, "unknown streamstats key.*server"),
            ({"notes": "karst"}, "notes must be a list of sentences"),
            ({"notes": [1]}, "notes must be a list of sentences"),
        ],
    )
    def test_a_bad_streamstats_record(self, record, message):
        with pytest.raises(ReportError, match=message):
            build_section("ungaged", {**TENNESSEE, "streamstats": record})

    def test_a_streamstats_record_needs_its_characteristics(self, fakes):
        with pytest.raises(ReportError, match="a streamstats record needs the characteristics"):
            build_section("ungaged", {**WASHINGTON, "streamstats": {"retrieved": "2026-11-15"}})
        assert fakes.delineate_calls == []

    def test_a_report_keeps_the_frozen_section(self, fakes, tmp_path):
        spec = ReportSpec("Site", [SectionSpec("ungaged", options=dict(WASHINGTON))])
        report = build_report(spec)
        assert spec.sections[0].options == WASHINGTON  # the caller's spec is not changed
        path = tmp_path / "site.json"
        report.spec.save(path)
        rebuilt = build_report(ReportSpec.load(path))
        assert len(fakes.delineate_calls) == 1
        assert rebuilt.to_markdown() == report.to_markdown()


# -- weighted: the at-site fit with the regression ----------------------------------------


@pytest.fixture(scope="module")
def weighted() -> Any:
    return build_section("weighted", {"gage": dict(GAGE), "regression": dict(TENNESSEE)})


class TestWeighted:
    def test_kind_is_registered(self):
        kind = SECTION_KINDS["weighted"]
        assert kind.required == ("gage", "regression") and kind.optional == ("detail",)

    def test_table_is_weight_with_regression(self, weighted):
        fit = run_ffa(
            np.array(FLOWS), np.array(YEARS), site_no="03606500", stationarity=False, **SKEW
        )
        estimate = estimate_from_characteristics(
            "TN", {"DRNAREA": 205.0}, region="HA4 DA only", lat=36.04, lon=-88.23
        )
        frame = weight_with_regression(
            fit["b17c"].results, [per[0][1] for per in estimate.estimates.values()]
        )
        frame = frame[frame["regression_flow_cfs"].notna()]
        table = _table(weighted.blocks, WEIGHTED_CAPTION)
        assert len(table.header) == 7 and len(table.rows) == len(estimate.estimates) > 0
        for row, (_, expected) in zip(table.rows, frame.iterrows()):
            assert row == (
                f"{expected['aep'] * 100:.3g}",
                cfs(expected["flow_cfs"]),
                cfs(expected["regression_flow_cfs"]),
                f"{expected['weight_at_site']:.2f}",
                cfs(expected["weighted_flow_cfs"]),
                cfs(expected["weighted_lower_cfs"]),
                cfs(expected["weighted_upper_cfs"]),
            )

    def test_the_one_percent_flow(self, weighted):
        by_aep = {row[0]: row for row in _table(weighted.blocks, WEIGHTED_CAPTION).rows}
        assert by_aep["1"][1:] == ("18,069", "24,552", "0.83", "19,019", "14,333", "25,237")

    def test_title_method_and_references(self, weighted):
        assert weighted.title == (
            "Weighted flood frequency at USGS 03606500 Big Sandy River at Bruceton, TN"
        )
        assert "Bulletin 17C Appendix 9 (England and others, 2019)" in _text(weighted.blocks)
        assert weighted.references[0] == BULLETIN_17C
        assert any("sir20245130" in r for r in weighted.references[1:])
        assert len(set(weighted.references)) == len(weighted.references)

    def test_detail_carries_both_analyses_one_level_down(self, weighted):
        gage = build_section("gaged_ffa", dict(GAGE))
        regression = build_section("ungaged", dict(TENNESSEE))
        headings = [(b.level, b.text) for b in weighted.blocks if isinstance(b, Heading)]

        def below(blocks: List[Any]) -> List[Any]:
            return [(b.level + 1, b.text) for b in blocks if isinstance(b, Heading)]

        assert (3, "Quantiles") in below(regression.blocks)
        assert headings == [
            (2, "At-site analysis"),
            *below(gage.blocks),
            (2, "Regression estimate"),
            *below(regression.blocks),
        ]
        assert len(weighted.blocks) == 4 + len(gage.blocks) + len(regression.blocks)
        captions = [b.caption for b in weighted.blocks if isinstance(b, Table) and b.caption]
        assert "Flood quantiles with 90-percent confidence limits." in captions
        assert "Regression flood quantiles with prediction limits." in captions

    def test_without_detail_the_section_is_the_method_and_the_table(self):
        content = build_section(
            "weighted", {"gage": dict(GAGE), "regression": dict(TENNESSEE), "detail": False}
        )
        assert [type(b) for b in content.blocks] == [Paragraph, Table]
        assert content.references[0] == BULLETIN_17C and len(content.references) == 2

    def test_the_options_as_built_hold_both_resolved_sets(self, weighted):
        assert weighted.options["gage"] == build_section("gaged_ffa", dict(GAGE)).options
        assert weighted.options["regression"] == TENNESSEE
        again = build_section("weighted", weighted.options)
        assert again.options == weighted.options
        assert _table(again.blocks, WEIGHTED_CAPTION) == _table(weighted.blocks, WEIGHTED_CAPTION)

    def test_a_regression_through_streamstats_is_frozen_in_the_section(self, monkeypatch):
        site = _wa_site()
        site["region_weights"] = {"4": 1.0}
        f = Fakes(site)
        monkeypatch.setattr("flowfreq.ungaged.delineate_and_get_characteristics", f.delineate)
        monkeypatch.setattr("flowfreq.ungaged.locate_regression_regions", f.locate)
        options = {"gage": dict(GAGE), "regression": dict(WASHINGTON), "detail": False}
        content = build_section("weighted", options)
        assert options["regression"] == WASHINGTON  # the caller's options are not changed
        stored = content.options["regression"]
        assert stored["region"] == {"4": 100.0}
        assert stored["streamstats"]["polygon_area_sq_mi"] == 82.2
        build_section("weighted", content.options)
        assert len(f.delineate_calls) == 1

    def test_an_extrapolated_regression_is_noted_per_row(self):
        regression = {
            **TENNESSEE,
            "characteristics": {"DRNAREA": 90000.0},
            "allow_extrapolation": True,
        }
        content = build_section(
            "weighted", {"gage": dict(GAGE), "regression": regression, "detail": False}
        )
        table = _table(content.blocks, WEIGHTED_CAPTION)
        assert table.header[-1] == "Note" and len(table.align) == len(table.header) == 8
        assert {row[-1] for row in table.rows} == {"extrapolated regression"}

    def test_a_basin_in_two_regions_is_refused(self):
        regression = {
            "state": "MT",
            "characteristics": {
                "CONTDA": 17.27,
                "EL6000": 0.0,
                "SLOP30_30M": 0.05,
                "ET0306MOD": 1.15,
            },
            "region": {"UYCM": 36.0, "ECP": 64.0},
            "lat": 46.0,
            "lon": -109.0,
        }
        with pytest.raises(ReportError, match="spans 2 regression regions"):
            build_section("weighted", {"gage": dict(GAGE), "regression": regression})

    def test_a_regression_without_a_variance_of_prediction_is_refused(self):
        regression = {
            "state": "NV",
            "characteristics": {"DRNAREA": 57, "ELEV": 6500},
            "region": "10",
            "lat": 39.5,
            "lon": -117.0,
        }
        with pytest.raises(ReportError, match="no AEP could be weighted.*no variance of predict"):
            build_section("weighted", {"gage": dict(GAGE), "regression": regression})

    @pytest.mark.parametrize(
        "change, message",
        [
            ({"gage": "03606500"}, "weighted: gage must be an object of gaged_ffa options"),
            ({"regression": ["TN"]}, "weighted: regression must be an object of ungaged options"),
            ({"gage": {**GAGE, "regional_skw": 0.1}}, "gaged_ffa: unknown option.*regional_skw"),
            ({"regression": {"state": "TN"}}, "ungaged: missing option.*lat"),
        ],
    )
    def test_bad_nested_options(self, change, message):
        options = {"gage": dict(GAGE), "regression": dict(TENNESSEE), **change}
        with pytest.raises(ReportError, match=message):
            build_section("weighted", options)

    def test_missing_part(self):
        with pytest.raises(ReportError, match="weighted: missing option.*regression"):
            build_section("weighted", {"gage": dict(GAGE)})

    def test_in_a_report(self, weighted, tmp_path):
        spec = ReportSpec(
            "Big Sandy River",
            [SectionSpec("weighted", options={"gage": dict(GAGE), "regression": dict(TENNESSEE)})],
            date="2026-11-15",
        )
        report = build_report(spec)
        text = report.to_markdown()
        assert "## 1. Weighted flood frequency at USGS 03606500" in text
        assert "### At-site analysis" in text and "### Regression estimate" in text
        assert "#### Quantiles" in text
        assert f"Table 1. {WEIGHTED_CAPTION}" in text
        path = tmp_path / "weighted.json"
        report.spec.save(path)
        assert build_report(ReportSpec.load(path)).to_markdown() == text

    def test_docx(self, tmp_path):
        docx = pytest.importorskip("docx")
        spec = ReportSpec(
            "Big Sandy River",
            [SectionSpec("weighted", options={"gage": dict(GAGE), "regression": dict(TENNESSEE)})],
        )
        path = build_report(spec).save(tmp_path / "weighted.docx")
        document = docx.Document(str(path))
        styles = {p.style.name for p in document.paragraphs}
        assert {"Heading 1", "Heading 2", "Heading 3"} <= styles
        cells = [c.text for t in document.tables for c in t.rows[0].cells]
        assert "Weight on the at-site flow" in cells
