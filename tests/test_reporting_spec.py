"""Tests for flowfreq.reporting.spec: the report specification and its JSON form."""

from __future__ import annotations

import json

import pytest

from flowfreq.reporting import SCHEMA, ReportError, ReportSpec, SectionSpec


def _spec() -> ReportSpec:
    return ReportSpec(
        title="Example",
        subtitle="Draft",
        author="A. Analyst",
        project="Bridge 12",
        date="2026-11-14",
        sections=[
            SectionSpec("text", title="Purpose", narrative="Why.", narrative_status="reviewed"),
            SectionSpec("gaged_ffa", options={"site_no": "03606500", "station_skew_only": True}),
        ],
    )


class TestRoundTrip:
    def test_dict_form(self):
        data = _spec().to_dict()
        assert data["schema"] == SCHEMA
        assert data["sections"][0] == {
            "kind": "text",
            "title": "Purpose",
            "narrative": "Why.",
            "narrative_status": "reviewed",
        }
        # What is unset is left out, so a hand-written spec stays short.
        assert data["sections"][1] == {
            "kind": "gaged_ffa",
            "options": {"site_no": "03606500", "station_skew_only": True},
        }
        assert ReportSpec.from_dict(data) == _spec()

    def test_json_and_file(self, tmp_path):
        spec = _spec()
        assert ReportSpec.from_json(spec.to_json()) == spec
        path = spec.save(tmp_path / "spec.json")
        assert ReportSpec.load(path) == spec

    def test_non_ascii_survives(self, tmp_path):
        spec = ReportSpec(title="Río Grande", sections=[])
        path = spec.save(tmp_path / "spec.json")
        assert "Río Grande" in path.read_text(encoding="utf-8")
        assert ReportSpec.load(path).title == "Río Grande"

    def test_minimal_spec(self):
        spec = ReportSpec.from_dict({"title": "T"})
        assert spec == ReportSpec(title="T")
        assert spec.to_dict() == {"schema": SCHEMA, "title": "T", "sections": []}


class TestRefusals:
    def test_unknown_report_key(self):
        with pytest.raises(ReportError, match="unknown key.*titel"):
            ReportSpec.from_dict({"title": "T", "titel": "x"})

    def test_unknown_section_key(self):
        with pytest.raises(ReportError, match="section 'text'.*unknown key.*body"):
            ReportSpec.from_dict({"title": "T", "sections": [{"kind": "text", "body": "x"}]})

    def test_other_schema(self):
        with pytest.raises(ReportError, match="schema 2"):
            ReportSpec.from_dict({"schema": 2, "title": "T"})

    def test_missing_title(self):
        with pytest.raises(ReportError, match="needs a title"):
            ReportSpec.from_dict({"sections": []})

    def test_not_an_object(self):
        with pytest.raises(ReportError, match="JSON object"):
            ReportSpec.from_dict([])  # type: ignore[arg-type]

    def test_sections_not_a_list(self):
        with pytest.raises(ReportError, match="sections must be a list"):
            ReportSpec.from_dict({"title": "T", "sections": {"kind": "text"}})

    def test_section_without_kind(self):
        with pytest.raises(ReportError, match="needs a 'kind'"):
            ReportSpec.from_dict({"title": "T", "sections": [{"title": "x"}]})

    def test_options_not_an_object(self):
        with pytest.raises(ReportError, match="options must be an object"):
            SectionSpec.from_dict({"kind": "ungaged", "options": [1, 2]})

    def test_bad_narrative_status(self):
        with pytest.raises(ReportError, match="draft, reviewed"):
            SectionSpec("text", narrative="x", narrative_status="final")
        with pytest.raises(ReportError, match="draft, reviewed"):
            SectionSpec.from_dict({"kind": "text", "narrative": "x", "narrative_status": "final"})

    def test_invalid_json(self):
        with pytest.raises(ReportError, match="not valid JSON"):
            ReportSpec.from_json("{title: T}")

    def test_nan_is_not_written(self):
        spec = ReportSpec("T", [SectionSpec("gaged_ffa", options={"regional_skew": float("nan")})])
        with pytest.raises(ValueError):
            spec.to_json()


def test_report_error_is_a_value_error():
    assert issubclass(ReportError, ValueError)


def test_json_is_plain_and_ends_with_a_newline():
    text = _spec().to_json()
    assert text.endswith("}\n")
    assert json.loads(text)["title"] == "Example"
