"""Hooks shared by the whole suite."""

import pytest

from tests.live_outage import pytest_runtest_makereport  # noqa: F401


@pytest.fixture(autouse=True)
def _no_bylocation_pause(request, monkeypatch):
    """Stubbed NSS ``bylocation`` calls skip the pause between requests.

    A live test keeps it: the pause is for the service (issue #365).
    """
    # Each test reads NSS's statistic-group list from its own stub, not a neighbour's.
    monkeypatch.setattr("flowfreq.streamstats._GROUP_CODES_BY_ID", {})
    if request.node.get_closest_marker("requires_network") is None:
        monkeypatch.setattr("flowfreq.streamstats.BYLOCATION_PAUSE", 0.0)


@pytest.fixture(autouse=True)
def _cold_fit_cache():
    """Each test fits its gages afresh: a fit memoised by a neighbour is not reused."""
    from flowfreq.reporting import gaged

    gaged.clear_fit_cache()
    yield
    gaged.clear_fit_cache()
