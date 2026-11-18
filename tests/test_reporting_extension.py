"""The `extension` option of a `gaged_ffa` section: a record extended from an index gage.

The records are Bulletin 17C Appendix 8's own example (Suwanee Creek extended from the
Etowah River), carried inline, so nothing here makes a request.
"""

from __future__ import annotations

from typing import Any, Dict, List

import pytest

from flowfreq.reporting import (
    Bullets,
    Heading,
    Paragraph,
    ReportError,
    Table,
    build_section,
    sections,
)
from tests.test_peak_extension import CONCURRENT, ETOWAH_BEFORE, EXTENDED

ETOWAH: Dict[str, Any] = {
    "water_year": list(range(1892, 1985)) + [y for y, _, _ in CONCURRENT],
    "peak_flow_cfs": [float(q) for q in ETOWAH_BEFORE] + [float(e) for _, e, _ in CONCURRENT],
}
SUWANEE: Dict[str, Any] = {
    "site_no": "02334885",
    "site_name": "Suwanee Creek at Suwanee, GA",
    "drainage_area_sqmi": 47.0,
    "peaks": {
        "water_year": [y for y, _, _ in CONCURRENT],
        "peak_flow_cfs": [float(s) for _, _, s in CONCURRENT],
    },
    "station_skew_only": True,
    "apply_peak_codes": False,
}
EXTENSION: Dict[str, Any] = {
    "index": {
        "site_no": "02392000",
        "site_name": "Etowah River at Canton, GA",
        "peaks": ETOWAH,
    }
}


def _table(blocks: List[Any], caption: str) -> Table:
    (table,) = [b for b in blocks if isinstance(b, Table) and b.caption == caption]
    return table


@pytest.fixture(scope="module")
def built() -> Any:
    return build_section("gaged_ffa", {**SUWANEE, "extension": EXTENSION})


def test_the_extension_is_the_bulletins(built):
    facts = dict(_table(built.blocks, "Record extension from the index gage.").rows)
    assert facts["Concurrent years"] == "20" and facts["Peaks added"] == "13"
    assert facts["Correlation of the concurrent log10 peaks"] == "0.8519"
    assert facts["Mean of log10 flow, improved"] == "3.3024"
    assert facts["Extension line"] == "log10 Q = 3.4364 + 1.4416 (log10 Q_index - 4.1414)"
    text = "\n".join(b.text for b in built.blocks if isinstance(b, Paragraph))
    assert "index gage USGS 02392000 Etowah River at Canton, GA" in text
    assert "from the index gage's peaks of 1972 to 1984" in text
    assert "not estimates of what happened in those years" in text
    assert built.title.endswith("record extended")


def test_the_fit_and_the_table_of_peaks_use_the_extended_record(built):
    record = dict(_table(built.blocks, "Station and peak-flow record.").rows)
    assert record["Annual peaks"] == "20" and record["Peaks added by record extension"] == "13"
    assert record["Water years"] == "1985 to 2004"
    peaks = _table(built.blocks, "Annual peak flows used in the analysis.")
    assert peaks.header[-1] == "Source" and len(peaks.rows) == 33
    extended = [row for row in peaks.rows if row[-1] == "extended"]
    assert [row[0] for row in extended] == [str(y) for y in range(1972, 1985)]
    for row, printed in zip(extended, EXTENDED):
        assert float(row[1].replace(",", "")) == pytest.approx(printed, rel=6e-3)
    alone = build_section("gaged_ffa", dict(SUWANEE))
    caption = "Fitted log-Pearson Type III parameters and low-outlier screening."
    own = dict(_table(alone.blocks, caption).rows)
    fitted = dict(_table(built.blocks, caption).rows)
    # The Bulletin's point: the 20 dry years alone sit below the extended record.
    assert float(fitted["Mean of log10 flow"]) > float(own["Mean of log10 flow"])


def test_no_stationarity_screen_on_an_extended_record(built):
    assert not any(isinstance(b, Heading) and b.text == "Stationarity screen" for b in built.blocks)
    alone = build_section("gaged_ffa", dict(SUWANEE))
    assert any(isinstance(b, Heading) and b.text == "Stationarity screen" for b in alone.blocks)


def test_the_index_record_is_downloaded_once_and_stored(monkeypatch):
    asked: List[str] = []

    def fake(site_no: str) -> Dict[str, Any]:
        asked.append(site_no)
        return {"peaks": dict(ETOWAH), "site_name": "ETOWAH RIVER AT CANTON, GA"}

    monkeypatch.setattr(sections, "download_peaks", fake)
    built = build_section(
        "gaged_ffa", {**SUWANEE, "extension": {"index": {"site_no": "02392000"}, "n_extend": 8}}
    )
    assert asked == ["02392000"]
    stored = built.options["extension"]
    assert stored["index"]["peaks"]["water_year"][0] == 1892 and stored["n_extend"] == 8
    assert stored["index"]["site_name"] == "ETOWAH RIVER AT CANTON, GA"
    again = build_section("gaged_ffa", built.options)
    assert asked == ["02392000"]
    assert (
        dict(_table(again.blocks, "Record extension from the index gage.").rows)["Peaks added"]
        == "8"
    )


def test_a_weak_relation_carries_its_caution():
    import numpy as np

    rng = np.random.default_rng(11)
    noisy = dict(SUWANEE["peaks"])
    noisy["peak_flow_cfs"] = [
        float(q * 10 ** rng.normal(0.0, 0.18)) for q in noisy["peak_flow_cfs"]
    ]
    try:
        built = build_section("gaged_ffa", {**SUWANEE, "peaks": noisy, "extension": EXTENSION})
    except ReportError as exc:
        assert "could not be extended" in str(exc)
    else:
        facts = dict(_table(built.blocks, "Record extension from the index gage.").rows)
        if float(facts["Correlation of the concurrent log10 peaks"]) <= 0.80:
            assert any(isinstance(b, Bullets) and "0.80" in b.items[0] for b in built.blocks)


@pytest.mark.parametrize(
    "extension,said",
    [
        ({"index": "02392000"}, "must be an object with an index gage"),
        ({"index": {}}, "needs site_no or peaks"),
        ({"index": {"peaks": ETOWAH}, "method": "MOVE.1"}, "unknown extension key"),
        ({"index": {"peaks": ETOWAH}, "years": [1990]}, "could not be extended"),
        (
            {"index": {"peaks": {"water_year": [1985, 1986], "peak_flow_cfs": [1.0, 2.0]}}},
            "could not be extended: only 2 concurrent",
        ),
    ],
)
def test_refusals(extension, said):
    with pytest.raises(ReportError, match=said):
        build_section("gaged_ffa", {**SUWANEE, "extension": extension})
