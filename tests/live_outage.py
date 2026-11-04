"""Report a live test that hit a service outage as an outage, not as a failure.

All six ``requires_network`` runs from 2026-10-05 to 2026-11-02 were red, not on one test but
on whichever service was down that week: WBD 504s, an NSS 500, a missing ss-hydro raster. A run
that is always red notices nothing, and skipping test by test only ever covered the test that
had just failed. So the rule for transport failures lives here, once.

Only a transport failure or a 5xx is an outage. A 4xx, a wrong answer, or any
``AssertionError`` (``assert response.status_code == 200`` on a 503 included) stays a
failure: those are how a service changing shape under this code gets noticed.

Outages are not silent: each is a skip whose reason starts with :data:`OUTAGE_PREFIX`, listed
by ``-rs`` and, under GitHub Actions, written to the job summary.
"""

from __future__ import annotations

import os
from typing import Any, Generator, Optional

import pytest
import requests

from flowfreq.streamstats import StreamStatsTransportError

OUTAGE_PREFIX = "service outage: "

_REASON_CHARS = 300


def outage_reason(exc: Optional[BaseException]) -> Optional[str]:
    """Return a skip reason if ``exc`` is a service outage, else ``None``.

    The cause chain is followed because library code wraps the transport error: the WBD
    lookup raises a bare ``RequestException`` from the ``HTTPError`` that carries the 504.
    An ``AssertionError`` is never an outage, whatever it was raised from.

    Parameters
    ----------
    exc : BaseException or None
        The exception a test raised.

    Returns
    -------
    str or None
        ``OUTAGE_PREFIX`` and the message, or ``None`` when it is not an outage.
    """
    if isinstance(exc, AssertionError):
        return None
    seen = set()
    current = exc
    while current is not None and id(current) not in seen:
        seen.add(id(current))
        if isinstance(
            current, (StreamStatsTransportError, requests.ConnectionError, requests.Timeout)
        ) or (
            isinstance(current, requests.HTTPError)
            and current.response is not None
            and current.response.status_code >= 500
        ):
            return f"{OUTAGE_PREFIX}{type(current).__name__}: {current}"[:_REASON_CHARS]
        current = current.__cause__ or current.__context__
    return None


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item: pytest.Item, call: Any) -> Generator[None, Any, None]:
    """Turn a ``requires_network`` test's outage into a skip, in its setup or its body."""
    outcome = yield
    report = outcome.get_result()
    if (
        call.excinfo is None
        or report.when not in ("setup", "call")
        or item.get_closest_marker("requires_network") is None
    ):
        return
    reason = outage_reason(call.excinfo.value)
    if reason is None:
        return
    report.outcome = "skipped"
    report.longrepr = (str(item.path), item.location[1] or 0, reason)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as fh:
            fh.write(f"- `{item.nodeid}`: {reason}\n")
