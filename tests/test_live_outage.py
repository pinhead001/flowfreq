"""What the live run counts as a service outage (tests/live_outage.py)."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest
import requests

from flowfreq.streamstats import StreamStatsTransportError
from tests.live_outage import OUTAGE_PREFIX, outage_reason

REPO_ROOT = Path(__file__).resolve().parent.parent


def _http_error(status: int) -> requests.HTTPError:
    response = requests.Response()
    response.status_code = status
    return requests.HTTPError(f"{status} error", response=response)


def _wrapped(inner: BaseException) -> requests.RequestException:
    """The WBD lookup's shape: a bare RequestException raised from the real error."""
    try:
        try:
            raise inner
        except BaseException as exc:
            raise requests.RequestException("WBD HU12 query failed") from exc
    except requests.RequestException as outer:
        return outer


class TestOutageReason:
    @pytest.mark.parametrize(
        "exc",
        [
            _http_error(500),
            _http_error(504),
            requests.ConnectionError("refused"),
            requests.Timeout("slow"),
            StreamStatsTransportError("failed after 3 attempts"),
            _wrapped(_http_error(504)),
        ],
    )
    def test_outages(self, exc: BaseException) -> None:
        reason = outage_reason(exc)
        assert reason is not None and reason.startswith(OUTAGE_PREFIX)

    @pytest.mark.parametrize(
        "exc",
        [
            None,
            _http_error(403),
            _http_error(404),
            requests.HTTPError("no response attached"),
            ValueError("unexpected JSON shape"),
            AssertionError("assert 503 == 200"),
            _wrapped(_http_error(404)),
        ],
    )
    def test_not_outages(self, exc: BaseException) -> None:
        assert outage_reason(exc) is None

    def test_assertion_raised_while_handling_an_outage_still_fails(self) -> None:
        try:
            try:
                raise _http_error(503)
            except requests.HTTPError:
                raise AssertionError("wrong answer")
        except AssertionError as exc:
            assert outage_reason(exc) is None

    def test_a_cycle_in_the_chain_terminates(self) -> None:
        a, b = ValueError("a"), ValueError("b")
        a.__cause__, b.__cause__ = b, a
        assert outage_reason(a) is None

    def test_long_reason_is_cut(self) -> None:
        reason = outage_reason(requests.Timeout("x" * 5000))
        assert reason is not None and len(reason) <= 300


SAMPLE = """
import pytest, requests

def _503():
    response = requests.Response()
    response.status_code = 503
    return requests.HTTPError("503 error", response=response)

@pytest.fixture
def needs_service():
    raise _503()

@pytest.mark.requires_network
def test_outage():
    raise _503()

@pytest.mark.requires_network
def test_outage_in_fixture(needs_service):
    pass

@pytest.mark.requires_network
def test_wrong_answer():
    assert 1 == 2

@pytest.mark.requires_network
def test_client_error():
    response = requests.Response()
    response.status_code = 404
    raise requests.HTTPError("404 error", response=response)

def test_offline_test_keeps_its_failure():
    raise _503()
"""


def test_hook_end_to_end(tmp_path: Path) -> None:
    """Outages skip with the reason listed and summarized; everything else still fails."""
    (tmp_path / "test_sample.py").write_text(SAMPLE, encoding="utf-8")
    (tmp_path / "conftest.py").write_text(
        "from tests.live_outage import pytest_runtest_makereport  # noqa: F401\n",
        encoding="utf-8",
    )
    (tmp_path / "pytest.ini").write_text(
        "[pytest]\nmarkers =\n    requires_network: live\n", encoding="utf-8"
    )
    summary = tmp_path / "summary.md"
    env = {
        **os.environ,
        "PYTHONPATH": str(REPO_ROOT),
        "GITHUB_STEP_SUMMARY": str(summary),
    }
    done = subprocess.run(
        [sys.executable, "-m", "pytest", "-rs", "-p", "no:cacheprovider", str(tmp_path)],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert "3 failed, 2 skipped" in done.stdout, done.stdout
    assert done.stdout.count(OUTAGE_PREFIX) == 2, done.stdout
    assert done.returncode == 1

    written = summary.read_text(encoding="utf-8")
    assert "test_outage`" in written and "test_outage_in_fixture`" in written
    assert "test_wrong_answer" not in written
    assert "test_offline_test_keeps_its_failure" not in written
