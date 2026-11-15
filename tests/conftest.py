"""Hooks shared by the whole suite."""

import pytest

from tests.live_outage import pytest_runtest_makereport  # noqa: F401


@pytest.fixture(autouse=True)
def _no_bylocation_pause(request, monkeypatch):
    """Stubbed NSS ``bylocation`` calls skip the pause between requests.

    A live test keeps it: the pause is for the service (issue #365).
    """
    if request.node.get_closest_marker("requires_network") is None:
        monkeypatch.setattr("flowfreq.streamstats.BYLOCATION_PAUSE", 0.0)
