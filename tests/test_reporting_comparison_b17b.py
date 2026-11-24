"""Bulletin 17B (superseded) as an estimate of a ``comparison`` section.

Everything here is offline: the peaks are inline and the regression carries its
characteristics.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterator, List

import matplotlib.pyplot as plt
import pytest

import flowfreq
from flowfreq import Bulletin17B, Bulletin17BError
from flowfreq.reporting import (
    SECTION_KINDS,
    ReportError,
    ReportSpec,
    SectionSpec,
    _base,
    build_report,
    build_section,
    comparison,
    render_markdown,
)
from flowfreq.reporting._base import _version, section_title, with_label
from flowfreq.reporting.b17b import B17B_CAPTION, BULLETIN_17B, _b17b_ffa
from flowfreq.reporting.blocks import Bullets, Figure, Heading, Paragraph, Table
from flowfreq.reporting.sections import COMPARISON_CAPTION, cfs
from flowfreq.reporting.session import ReportSession
from tests.test_reporting_b17b import LABEL, REFUSED, REGULATED_SITE, peaks_spec
from tests.test_reporting_cautions import CAUTION, URBAN
from tests.test_reporting_comparison import ESTIMATES as ORDINARY
from tests.test_reporting_comparison import SITE, TRANSPOSED
from tests.test_reporting_sites import TENNESSEE
from tests.test_reporting_summary import LOCATIONS, SUMMARY

SKEW: Dict[str, Any] = {"regional_skew": 0.6, "regional_skew_se": 0.55}
LIMITS_CAPTION = "Lower and upper limits of each estimate."
CURVES = f"Bulletin 17C and {LABEL} frequency curves."
STORED = Path(__file__).parent / "fixtures" / "comparison_without_b17b.json"


def gage_opts(name: str) -> Dict[str, Any]:
    """Options that a ``gaged_ffa`` and a ``b17b_ffa`` section both take."""
    return {"peaks": peaks_spec(name), **SKEW}


def b17c(label: str = "Bulletin 17C", **extra: Any) -> Dict[str, Any]:
    options = {**gage_opts("example1"), "stationarity": False, "peaks_table": False, **extra}
    return {"label": label, "kind": "gaged_ffa", "options": options}


def b17b(label: str = "Procedure of 1982", **extra: Any) -> Dict[str, Any]:
    return {"label": label, "kind": "b17b_ffa", "options": {**gage_opts("example1"), **extra}}


def refused(label: str = "Old procedure") -> Dict[str, Any]:
    return b17b(label, low_outlier_threshold=2190.0)


def is_figure(block: Any) -> bool:
    return isinstance(block, Figure)


def compare(*estimates: Dict[str, Any], **options: Any) -> Any:
    return build_section("comparison", {"estimates": list(estimates), **options})


def table(content: Any, caption: str = COMPARISON_CAPTION) -> Table:
    (found,) = [b for b in content.blocks if isinstance(b, Table) and b.caption == caption]
    return found


def text_of(blocks: List[Any]) -> str:
    out = []
    for b in blocks:
        out.extend(b.items if isinstance(b, Bullets) else [getattr(b, "text", "")])
    return "\n".join(out)


def headings(content: Any) -> List[str]:
    return [b.text for b in content.blocks if isinstance(b, Heading) and b.level == 2]


@pytest.fixture
def drawn(monkeypatch: pytest.MonkeyPatch) -> List[Any]:
    """What each curve figure was asked to draw, in order."""
    calls: List[Any] = []

    def record(curves: Any, aeps: Any, where: str, **style: Any) -> bytes:
        calls.append((list(curves), list(aeps), where))
        return _base._curve_figure(curves, aeps, where, **style)

    monkeypatch.setattr(comparison, "_curve_figure", record)
    return calls


# -- the table ----------------------------------------------------------------------------


def test_comparison_puts_b17b_beside_b17c():
    spec = {
        "estimates": [
            {"label": "Bulletin 17C", "kind": "gaged_ffa", "options": gage_opts("example1")},
            {"label": "Bulletin 17B", "kind": "b17b_ffa", "options": gage_opts("example1")},
        ],
        "detail": False,
    }
    content = build_section("comparison", spec)
    assert sum(1 for b in content.blocks if is_figure(b)) == 1
    assert table(content).header == (
        "Annual exceedance probability (percent)",
        "Bulletin 17C (cfs)",
        f"{LABEL} (cfs)",
        "Largest over smallest",
    )
    _, res = _b17b_ffa(gage_opts("example1"))
    own = res.quantiles.set_index("aep")
    rows = {r[0]: r for r in table(content).rows}
    assert rows["1"][2] == cfs(own.loc[0.01, "flow_cfs"])
    limits = {r[0]: r for r in table(content, LIMITS_CAPTION).rows}
    assert limits["1"][2] == f"{cfs(own.loc[0.01, 'lower'])} to {cfs(own.loc[0.01, 'upper'])}"
    flows = [float(c.replace(",", "")) for c in rows["1"][1:3]]
    assert float(rows["1"][3]) == pytest.approx(max(flows) / min(flows), abs=0.006)
    assert BULLETIN_17B in content.references
    said = text_of(content.blocks)
    assert f"; {LABEL}, 90-percent confidence limits (Bulletin 17B)" in said
    assert "superseded by Bulletin 17C" in said and "estimated by 2 methods" in said
    assert "gives no answer" not in said and content.cautions == []


def test_comparison_matches_aeps_with_b17c():
    content = compare(b17c(), b17b(), detail=False)
    shown = [r[0] for r in table(content).rows]
    assert shown == ["99", "90", "50", "20", "10", "4", "2", "1", "0.5", "0.2"]
    assert all("-" not in r for r in table(content).rows)
    # 99.5 percent is Bulletin 17C's alone, and is not taken for 99.
    own = comparison._compared(0, b17c())[2].set_index("aep")["flow_cfs"]
    assert 0.995 in own.index and table(content).rows[0][1] == cfs(own[0.99])
    with pytest.raises(ReportError, match="not every estimate has a quantile at AEP 0.995"):
        compare(b17c(), b17b(), aeps=[0.995])


def test_three_columns_split_procedure_and_skew():
    low_outlier = {"peaks": peaks_spec("example3")}
    content = compare(
        b17c("B17C, today", **low_outlier),
        b17b("B17B, the same skew", **low_outlier),
        b17b("B17B, the old skew", **low_outlier, generalized_skew=-0.3),
        detail=False,
    )
    assert table(content).header[1:4] == (
        "B17C, today (cfs)",
        f"B17B, the same skew, {LABEL} (cfs)",
        f"B17B, the old skew, {LABEL} (cfs)",
    )
    row = {r[0]: r for r in table(content).rows}["1"]
    flows = [float(c.replace(",", "")) for c in row[1:4]]
    assert len(set(flows)) == 3 and flows[2] < flows[1]
    assert "estimated by 3 methods" in text_of(content.blocks)
    assert content.title == f"Comparison of flood estimates, with {LABEL}"


def test_the_procedures_agree_on_a_complete_record():
    """No outlier, no history, the same skew: the moments and weights are the same."""
    content = compare(b17c(), b17b(), detail=False)
    assert {r[3] for r in table(content).rows} == {"1.00"}
    assert all(r[1] == r[2] for r in table(content).rows)


def test_options_as_built_build_the_same_section_again():
    content = compare(b17c(), b17b(), refused(), detail=False)
    kept = content.options["estimates"]
    assert [e["label"] for e in kept] == ["Bulletin 17C", "Procedure of 1982", "Old procedure"]
    again = build_section("comparison", json.loads(json.dumps(content.options)))
    assert again.blocks == content.blocks and again.options == content.options
    assert again.cautions == content.cautions


def test_a_bad_option_still_raises():
    with pytest.raises(ReportError, match="low_outlier_threshold must not be negative"):
        compare(b17c(), b17b(low_outlier_threshold=-1.0))
    with pytest.raises(ReportError, match=r"unknown option\(s\) stationarity"):
        compare(b17c(), b17b(stationarity=False))
    with pytest.raises(ReportError, match="kind must be one of gaged_ffa.*b17b_ffa; got 'mean'"):
        compare(b17c(), {"label": "x", "kind": "mean", "options": {}})


# -- labels -------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "label, shown",
    [
        ("Gage", f"Gage, {LABEL}"),
        (f"Old study, {LABEL}", f"Old study, {LABEL}"),
        ("bulletin 17b (SUPERSEDED), old skew", "bulletin 17b (SUPERSEDED), old skew"),
        ("Bulletin 17B", LABEL),
        ("bulletin 17b, old skew", "bulletin 17b, old skew (superseded)"),
        ("Superseded bulletin 17B fit", "Superseded bulletin 17B fit"),
        ("Superseded", f"Superseded, {LABEL}"),
        ("Bulletin 17C", f"Bulletin 17C, {LABEL}"),
    ],
)
def test_column_headers_and_headings_say_superseded_once(label, shown):
    content = compare(b17c("Gage fit"), b17b(label))
    assert table(content).header[2] == f"{shown} (cfs)"
    assert table(content, LIMITS_CAPTION).header[2] == f"{shown} (cfs)"
    assert headings(content) == ["Gage fit", shown]


def test_the_subsection_keeps_its_own_labels():
    content = compare(b17c(), b17b())
    (own,) = [b for b in content.blocks if isinstance(b, Table) and b.caption == B17B_CAPTION]
    assert own.rows
    assert any(b.level == 3 and LABEL in b.text for b in content.blocks if isinstance(b, Heading))


def test_the_same_shown_label_is_refused():
    with pytest.raises(ReportError, match="two estimates have the same label"):
        compare(b17b("Gage"), b17b(f"Gage, {LABEL}"))
    with pytest.raises(ReportError, match="two estimates have the same label"):
        compare(b17b("Bulletin 17B"), b17b(LABEL))
    with pytest.raises(ReportError, match="two estimates have the same label"):
        compare(b17c("Gage"), b17c("Gage"))


def test_title_says_superseded():
    assert compare(b17c(), b17b(), site_name=SITE, detail=False).title == (
        f"Flood quantiles at {SITE} by 2 methods, with {LABEL}"
    )
    given = compare(b17c(), b17b(), site_name="Bulletin 17B (Superseded) Creek", detail=False)
    assert given.title == "Flood quantiles at Bulletin 17B (Superseded) Creek by 2 methods"


def test_a_given_title_keeps_the_label_in_a_report():
    both = {"estimates": [b17c(), b17b()], "detail": False}
    spec = ReportSpec(
        title="Report",
        sections=[
            SectionSpec("comparison", title="Old and new", options=both),
            SectionSpec("comparison", title="old and new, bulletin 17b (superseded)", options=both),
            SectionSpec("comparison", options=both),
            SectionSpec(
                "comparison",
                title="Two skews",
                options={"estimates": [b17c("A"), b17c("B", regional_skew=0.1)], "detail": False},
            ),
        ],
    )
    built = build_report(spec)
    shown = [b.text for b in built.blocks if isinstance(b, Heading) and b.level == 1][:4]
    assert shown == [
        f"1. Old and new, with {LABEL}",
        "2. old and new, bulletin 17b (superseded)",
        f"3. Comparison of flood estimates, with {LABEL}",
        "4. Two skews",
    ]
    assert built.spec.sections[0].title == "Old and new"


def test_section_title_label_follows_the_options():
    with_b17b = {"estimates": [b17c(), b17b()]}
    without = {"estimates": [b17c()]}
    assert section_title("comparison", "Given", "Own", with_b17b) == f"Given, with {LABEL}"
    assert section_title("comparison", None, "Own", with_b17b) == f"Own, with {LABEL}"
    assert section_title("comparison", f"X, {LABEL}", "Own", with_b17b) == f"X, {LABEL}"
    assert section_title("comparison", "Given", "Own", without) == "Given"
    assert section_title("b17b_ffa", "Given", "Own", {}) == f"Given, {LABEL}"
    assert with_label("", LABEL) == "" and with_label("A", "") == "A"


# -- refusals -----------------------------------------------------------------------------


def test_refused_b17b_estimate_does_not_fail_the_comparison():
    content = compare(b17c(), refused(), detail=False)
    assert table(content).header[2] == f"Old procedure, {LABEL} (cfs)"
    rows = table(content).rows
    assert [r[0] for r in rows][:2] == ["99.5", "99"]  # the AEPs of the estimate with numbers
    assert all(r[2] == "-" and r[3] == "-" and r[1] != "-" for r in rows)
    limits = table(content, LIMITS_CAPTION)
    assert all(r[2] == "-" and " to " in r[1] for r in limits.rows)
    said = text_of(content.blocks)
    assert "1 of them gives no answer" in said
    assert "90-percent confidence limits (Bulletin 17B)" not in said
    assert content.title == f"Comparison of flood estimates, with {LABEL}"


def test_a_refused_estimate_is_left_out_of_the_ratio():
    alone = compare(b17c(), b17b("Old"), detail=False)
    content = compare(b17c(), refused(), b17b("Old"), detail=False)
    assert [(r[0], r[1], r[3], r[4]) for r in table(content).rows] == [
        tuple(r) for r in table(alone).rows
    ]
    assert all(r[2] == "-" for r in table(content).rows)


@pytest.mark.parametrize("detail", [False, True])
def test_the_refusal_is_printed_once(detail):
    content = compare(b17c(), refused(), b17b("Old"), detail=detail)
    said = text_of(content.blocks)
    assert said.count(REFUSED) == 1
    assert "Cautions on" not in said
    (caution,) = content.cautions
    assert caution.startswith(f"Old procedure, {LABEL}: {REFUSED}")
    if detail:
        assert headings(content)[1] == f"Old procedure, {LABEL}"
    else:
        (shown,) = [b for b in content.blocks if REFUSED in getattr(b, "text", "")]
        assert isinstance(shown, Paragraph) and shown.text == caution


def test_a_regression_caution_still_prints_beside_a_refusal():
    urban = {"label": "Urban regression", "kind": "ungaged", "options": dict(URBAN)}
    content = compare(b17c(), refused(), urban, detail=False)
    said = text_of(content.blocks)
    assert said.count(f"Cautions on Urban regression: {CAUTION}") == 1
    assert said.count(REFUSED) == 1 and len(content.cautions) == 2


def test_a_regulated_gage_is_a_refusal_too():
    content = compare(b17c(), b17b("Old", site_no=REGULATED_SITE), detail=False)
    assert all(r[2] == "-" for r in table(content).rows)
    assert f"USGS {REGULATED_SITE} is classed regulated" in text_of(content.blocks)


@pytest.mark.parametrize("detail", [False, True])
def test_every_estimate_refused(detail, drawn):
    content = compare(refused("A"), refused("B"), site_name=SITE, detail=detail, aeps=[0.01])
    tables = [b.caption for b in content.blocks if isinstance(b, Table)]
    assert COMPARISON_CAPTION not in tables and LIMITS_CAPTION not in tables
    assert not any(is_figure(b) for b in content.blocks) and drawn == []
    said = text_of(content.blocks)
    assert "no table" in said and said.count(REFUSED) == 2
    assert content.title == f"Flood quantiles at {SITE} by 2 methods, with {LABEL}"
    assert len(content.cautions) == 2
    compare(refused("A"), refused("B"), detail=detail)
    with pytest.raises(ReportError, match="aeps must be a list of numbers"):
        compare(refused("A"), refused("B"), aeps=["x"])


# -- the figure ---------------------------------------------------------------------------


def caption_of(content: Any) -> str:
    (figure,) = [b for b in content.blocks if is_figure(b)]
    assert figure.png.startswith(b"\x89PNG")
    return figure.caption


def test_the_figure_draws_every_estimate_with_its_limits(drawn):
    estimates = [
        b17c(),
        b17b(),
        {"label": "Regression", "kind": "ungaged", "options": dict(TENNESSEE)},
    ]
    content = compare(*estimates, site_name=SITE, detail=False)
    ((curves, aeps, where),) = drawn
    assert [label for label, _ in curves] == [
        "Bulletin 17C",
        f"Procedure of 1982, {LABEL}",
        "Regression",
    ]
    assert where == SITE and [f"{a * 100:.3g}" for a in aeps] == [r[0] for r in table(content).rows]
    for _, flows in curves:
        assert flows[["lower", "upper"]].notna().all(axis=None)
    caption = caption_of(content)
    assert caption.startswith(CURVES)
    assert f"Superseded, dashed: Procedure of 1982, {LABEL}." in caption
    assert "Regression, 90-percent prediction limits" in caption
    assert "90-percent confidence limits (Bulletin 17B)" in caption


def test_no_curve_goes_beyond_the_table(drawn):
    compare(b17c(), b17b(), detail=False)
    ((curves, _, _),) = drawn
    assert all(0.0019 < f["aep"].min() and f["aep"].max() < 0.9901 for _, f in curves)


def test_a_refused_estimate_is_not_drawn(drawn):
    content = compare(b17c(), refused(), b17b("Old"), detail=False)
    ((curves, _, _),) = drawn
    assert [label for label, _ in curves] == ["Bulletin 17C", f"Old, {LABEL}"]
    caption = caption_of(content)
    assert caption.startswith(CURVES) and f"Superseded, dashed: Old, {LABEL}." in caption
    assert f"Not drawn: Old procedure, {LABEL}." in caption


def test_the_caption_names_only_what_is_drawn():
    only_17c = caption_of(compare(b17c(), refused(), detail=False))
    assert only_17c.startswith("Bulletin 17C frequency curves.")
    assert "Superseded, dashed:" not in only_17c
    two = caption_of(compare(b17b("A"), b17b("B", generalized_skew=-0.3), detail=False))
    assert two.startswith(f"{LABEL} frequency curves.")
    assert f"Superseded, dashed: A, {LABEL}; B, {LABEL}." in two
    moved = {"label": "Moved", "kind": "transposed", "options": dict(TRANSPOSED)}
    assert caption_of(compare(moved, b17b(), detail=False)).startswith(CURVES)
    regression = {"label": "R", "kind": "ungaged", "options": dict(TENNESSEE)}
    assert caption_of(compare(regression, refused(), detail=False)).startswith("Frequency curves.")


@pytest.fixture
def figures(monkeypatch: pytest.MonkeyPatch) -> Iterator[List[Any]]:
    """The figures drawn, left open."""
    shown: List[Any] = []

    def keep(fig: Any) -> bytes:
        fig.canvas.draw()
        shown.append(fig)
        return b"png"

    monkeypatch.setattr(_base, "_figure_png", keep)
    yield shown
    for fig in shown:
        plt.close(fig)


def legend_of(fig: Any) -> List[str]:
    return [t.get_text() for t in fig.axes[0].get_legend().get_texts()]


def curves_of(fig: Any) -> List[Any]:
    """The curves, without the limits."""
    return [line for line in fig.axes[0].get_lines() if line.get_marker() == "o"]


def test_the_figure_is_on_the_summary_axis(figures):
    compare(b17c(), b17b(), detail=False)
    (fig,) = figures
    (ax,) = fig.axes
    assert ax.get_yscale() == "log"
    assert [t.get_text() for t in ax.get_xticklabels()][-2:] == ["0.5", "0.2"]
    assert legend_of(fig) == ["Bulletin 17C", f"Procedure of 1982, {LABEL}"]
    assert len(ax.collections) == 2  # a band for each


@pytest.mark.parametrize("label", ["_hidden", "$x^{$ in dollars", "$_{$", "50% $"])
def test_any_label_is_drawn_as_written(figures, label):
    compare(b17c(label), b17b(), site_name="Cash $^$ Creek", detail=False)
    (fig,) = figures
    assert legend_of(fig) == [label, f"Procedure of 1982, {LABEL}"]
    assert fig.axes[0].get_title() == "Flood frequency at Cash $^$ Creek"


@pytest.mark.parametrize("label", ["_hidden", "$x^{$ in dollars"])
def test_any_label_is_drawn_as_written_in_a_summary(figures, label):
    estimates = [{**ORDINARY[0], "label": label}]
    build_section("site_summary", {**SUMMARY, "estimates": estimates})
    assert legend_of(figures[-1]) == [label]


def test_the_legend_covers_nothing(figures):
    long = "Regional regression equations of the study, drainage area only, " * 3
    regression = {"label": long, "kind": "ungaged", "options": dict(TENNESSEE)}
    compare(b17c(), b17b(), b17b("Old skew", generalized_skew=-0.3), regression, detail=False)
    (fig,) = figures
    (ax,) = fig.axes
    box, legend = ax.get_window_extent(), ax.get_legend().get_window_extent()
    assert legend.y1 < min(t.get_window_extent().y0 for t in ax.get_xticklabels())
    assert legend.y1 < ax.xaxis.label.get_window_extent().y0
    assert legend.width < 1.3 * box.width  # long labels are wrapped
    assert " ".join(legend_of(fig)[-1].split()) == long.strip()


def test_a_summary_keeps_its_legend_inside(figures):
    build_section("site_summary", dict(SUMMARY))
    ax = figures[-1].axes[0]
    assert ax.get_window_extent().contains(*ax.get_legend().get_window_extent().p0)
    assert {line.get_linestyle() for line in curves_of(figures[-1])} == {"-"}


def test_coinciding_curves_both_show(figures):
    compare(b17b("First"), b17c(), b17b("Second", generalized_skew=-0.3), detail=False)
    old, new, other = curves_of(figures[0])
    at = dict(zip(new.get_xdata().round(6), new.get_ydata()))
    both = [(at[x], y) for x, y in zip(old.get_xdata().round(6), old.get_ydata()) if x in at]
    assert len(both) == 10 and all(a == pytest.approx(b) for a, b in both)
    assert new.get_linestyle() == "-"
    assert old.get_linestyle() == "--" and other.get_linestyle() == "--"
    assert min(old.get_zorder(), other.get_zorder()) > new.get_zorder()


def test_the_figure_leaves_no_pyplot_figure_open():
    before = plt.get_fignums()
    compare(b17c(), b17b(), detail=False)
    assert plt.get_fignums() == before


# -- nothing else changes -----------------------------------------------------------------


def test_comparison_without_b17b_is_unchanged():
    spec = {"estimates": [b17c("A"), b17c("B", regional_skew=0.1)], "detail": False}
    content = build_section("comparison", spec)
    assert not any(is_figure(b) for b in content.blocks)
    assert "17B" not in text_of(content.blocks) and "17B" not in content.title


def _plain(block: Any) -> List[Any]:
    fields = {k: v for k, v in vars(block).items() if not isinstance(v, bytes)}
    return [type(block).__name__, fields]


@pytest.mark.parametrize("name", ["detail", "brief", "cautioned"])
def test_comparison_without_b17b_is_as_it_was_stored(name):
    """The whole section, as commit 0518b01 built it (figures without their pixels)."""
    urban = {"label": "Urban", "kind": "ungaged", "options": dict(URBAN)}
    options: Dict[str, Dict[str, Any]] = {
        "detail": {"site_name": SITE, "estimates": ORDINARY},
        "brief": {"site_name": SITE, "estimates": ORDINARY, "detail": False},
        "cautioned": {"estimates": [ORDINARY[0], urban], "detail": False, "aeps": [0.1, 0.01]},
    }
    content = build_section("comparison", options[name])
    stored = json.loads(STORED.read_text(encoding="utf-8").replace("{version}", _version()))[name]
    built = {
        "title": content.title,
        "blocks": [_plain(b) for b in content.blocks],
        "options": content.options,
        "references": content.references,
        "cautions": content.cautions,
        "markdown": render_markdown([b for b in content.blocks if not is_figure(b)]),
    }
    assert json.loads(json.dumps(built)) == stored


def test_nothing_computes_b17b_without_the_kind(monkeypatch):
    monkeypatch.setattr(
        "flowfreq.reporting.comparison._b17b_ffa", lambda *a: pytest.fail("computed")
    )
    compare(b17c("A"), b17c("B", regional_skew=0.1), detail=False)
    build_section("site_summary", dict(SUMMARY))


def test_site_summary_refuses_b17b():
    with pytest.raises(ReportError, match="site_summary: estimate 2: .*Bulletin 17B") as said:
        build_section("site_summary", {**SUMMARY, "estimates": [b17c(), b17b()]})
    message = str(said.value)
    assert "kind must be one of gaged_ffa, transposed, ungaged; got 'b17b_ffa'" in message
    assert "a b17b_ffa section or a comparison, not in a summary" in message
    with pytest.raises(ReportError) as said:
        estimates = [{**ORDINARY[0], "kind": "weighted"}]
        build_section("site_summary", {**SUMMARY, "estimates": estimates})
    assert "17B" not in str(said.value)


def test_only_the_labelled_kinds_have_a_label():
    labelled = {name for name, kind in SECTION_KINDS.items() if kind.title_label}
    assert labelled == {"b17b_ffa", "comparison"}


# -- a session ----------------------------------------------------------------------------


def test_session_compares_and_will_not_summarize(tmp_path):
    session = ReportSession("Report", output_dir=tmp_path)
    session.add_section("gaged_ffa", b17c()["options"])
    session.add_section("b17b_ffa", b17b()["options"])
    added = session.compare(["s1", "s2"], labels=["17C", "Old"], title="Old and new")
    assert added["title"] == f"Old and new, with {LABEL}"
    assert COMPARISON_CAPTION in added["result"] or "17C (cfs)" in added["result"]
    assert f"Old, {LABEL} (cfs)" in added["result"] and "cautions" not in added
    default = session.compare(["s1", "s2"])
    assert f"Flood frequency, {LABEL} (cfs)" in default["result"]
    session.render("report.md")
    assert f"Old and new, with {LABEL}" in (tmp_path / "report.md").read_text(encoding="utf-8")
    before = session.outline()
    with pytest.raises(ReportError, match="section s2 is a b17b_ffa section") as said:
        session.summarize(["s1", "s2"], LOCATIONS)
    assert "not in a summary" in str(said.value)
    session.add_text("Words", "Words.")
    with pytest.raises(ReportError, match="gaged_ffa, transposed, ungaged or b17b_ffa section"):
        session.compare(["s1", "s5"])
    assert session.outline()["sections"][:-1] == before["sections"]


def test_session_reply_carries_the_refusal(tmp_path):
    session = ReportSession("Report", output_dir=tmp_path)
    session.add_section("gaged_ffa", b17c()["options"])
    session.add_section("b17b_ffa", refused()["options"])
    session.add_section("b17b_ffa", refused()["options"])
    reply = session.compare(["s1", "s2"], labels=["17C", "Old"])
    (said,) = reply["cautions"]
    assert said.startswith(f"Old, {LABEL}: {REFUSED}")
    none = session.compare(["s2", "s3"], labels=["Old", "Older"])
    assert none["result"] is None and len(none["cautions"]) == 2


# -- the package --------------------------------------------------------------------------


def test_the_package_exports_the_class():
    assert flowfreq.Bulletin17B is Bulletin17B and issubclass(Bulletin17BError, ValueError)
    for name in ("Bulletin17B", "Bulletin17BError", "Bulletin17BResults", "bulletin17b_kwargs"):
        assert name in flowfreq.__all__
        assert getattr(flowfreq, name) is getattr(flowfreq.bulletin17b, name)
