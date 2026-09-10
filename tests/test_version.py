"""The package version must agree with pyproject.toml.

``flowfreq.__version__`` was a hand-maintained literal for most of this
project's life and drifted every time: it read "0.3.0" through the 0.4.0
release and "0.4.0" through 0.5.0, 0.6.0 and 0.6.1. It now derives from the
installed distribution metadata, which cannot disagree with what pip
installed -- but that only helps if what pip installed matches the tree, so
this checks the round trip.

The version is read out of pyproject.toml with a regex rather than a TOML
parser. Historical note: this was forced when the runtime floor was 3.9
(``tomllib`` is 3.11+ and ``tomli`` was never a declared dev dependency); the
floor is 3.11 now, so ``tomllib`` would parse fine here too, but there is no
strong reason to switch a working regex for a stdlib import that would only
matter if this file's ``version = "..."`` line ever became ambiguous to
match, which it is not -- ``[tool.black]`` has ``target-version``, which does
not match the anchor.
"""

from __future__ import annotations

import pathlib
import re

import flowfreq


def _pyproject_version() -> str:
    root = pathlib.Path(__file__).resolve().parents[1]
    text = (root / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version = "([^"]+)"', text, re.MULTILINE)
    assert match is not None, 'no `version = "..."` line found in pyproject.toml'
    return match.group(1)


def test_version_is_reported():
    """ "unknown" means the fallback fired, which should not happen in a test
    run against an installed package."""
    assert flowfreq.__version__
    assert flowfreq.__version__ != "unknown"


def test_version_matches_pyproject():
    """If this fails locally after a version bump, the installed package is
    stale rather than the code being wrong: reinstall with
    ``pip install -e ".[dev]"``. CI installs fresh every run, so a failure
    there is a genuine disagreement."""
    expected = _pyproject_version()
    assert flowfreq.__version__ == expected, (
        f"flowfreq.__version__ is {flowfreq.__version__!r} but pyproject.toml says "
        f'{expected!r}; if you just bumped the version, reinstall with pip install -e ".[dev]"'
    )


def test_version_is_not_a_hardcoded_literal_again():
    """Guards the actual regression. A literal reintroduced here would sit at
    whatever was current when it was typed and then quietly stop being true,
    which is exactly what happened three times."""
    source = pathlib.Path(flowfreq.__file__).read_text(encoding="utf-8")
    assert '__version__ = "0.' not in source, (
        "__version__ looks like a hardcoded literal again; derive it from "
        "importlib.metadata so it cannot drift from pyproject.toml"
    )
