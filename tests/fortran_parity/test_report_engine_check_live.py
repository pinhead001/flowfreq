"""The `engine_check` report section with the reference it makes itself.

``tests/test_reporting_engine_check.py`` covers the section from stored numbers, without
the extension. Here the section calls ``emafitpr`` on the fixture's record, as it does when
no reference is given, and the report it writes has to say the two engines agree.
"""

from __future__ import annotations

import datetime as dt

import pytest

pytest.importorskip(
    "flowfreq.peakfqr",
    reason="Fortran extension not built; run python build_fortran/build.py "
    "(see docs/FORTRAN_UPLOAD.md)",
    # See tests/fortran_parity/test_fortran_oracles.py's exc_type comment: pytest
    # 9.1 narrowed importorskip's default exc_type to ModuleNotFoundError, but this
    # module deliberately re-raises as ImportError.
    exc_type=ImportError,
)

pytestmark = pytest.mark.requires_fortran


@pytest.fixture(scope="module")
def content():
    from flowfreq.reporting import build_section
    from tests.test_reporting_sites import GAGE

    return build_section("engine_check", {"gage": dict(GAGE), "detail": False})


def test_the_engines_agree_on_big_sandy(content):
    from flowfreq.reporting import Paragraph
    from tests.test_reporting_engine_check import AGREE

    verdict = content.blocks[0]
    assert isinstance(verdict, Paragraph)
    assert AGREE in verdict.text and "Outside tolerance" not in verdict.text
    assert "Neither censors a value as a low outlier." in verdict.text


def test_the_reference_is_stored_as_the_tests_without_the_extension_have_it(content):
    """The numbers ``tests/test_reporting_engine_check.py`` carries are this call's.

    Held to 1e-4 relative rather than exactly: the stored flows are rounded to a
    thousandth of a cfs, and a different gfortran moves the fit's last digits (CLAUDE.md,
    Validation Status).
    """
    from tests.test_reporting_engine_check import REFERENCE

    stored = content.options["reference"]
    assert stored["source"].startswith("USGS peakfq 8.1.0, its Fortran routine emafitpr")
    assert stored["computed"] == dt.date.today().isoformat()
    assert stored["record"] == REFERENCE["record"]
    assert stored["low_outlier_count"] == REFERENCE["low_outlier_count"]
    assert stored["parameters"] == pytest.approx(REFERENCE["parameters"], rel=1e-4)
    assert [row["aep"] for row in stored["quantiles"]] == [
        row["aep"] for row in REFERENCE["quantiles"]
    ]
    for row, expected in zip(stored["quantiles"], REFERENCE["quantiles"]):
        assert row == pytest.approx(expected, rel=1e-4), row["aep"]


def test_it_rebuilds_from_what_it_stored(content, monkeypatch):
    from flowfreq.bulletin17c import Bulletin17C
    from flowfreq.reporting import build_section

    def absent(self, *args, **kwargs):
        raise AssertionError("the stored reference should have been used")

    monkeypatch.setattr(Bulletin17C, "fortran_reference", absent)
    again = build_section("engine_check", content.options)
    assert again.blocks == content.blocks


def test_through_the_session(tmp_path):
    from flowfreq.reporting.session import ReportSession
    from tests.test_reporting_engine_check import AGREE
    from tests.test_reporting_sites import GAGE

    session = ReportSession("Big Sandy River", output_dir=tmp_path)
    gage = session.add_section("gaged_ffa", GAGE)
    added = session.engine_check(gage["section_id"])
    assert added["kind"] == "engine_check" and AGREE in added["result"]
