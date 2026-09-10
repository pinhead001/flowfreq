"""The formatter config must be one the pinned formatter understands.

``[tool.black]`` gained ``py314`` in 82b1dab when the test matrix rose to
3.14, while the dev extra still pinned ``black>=24.0,<25``. black 24 has no
``py314`` target, so it exited with a usage error before formatting anything:
``make lint`` was not reporting a style problem, it was not running, and every
file went unchecked. That broke CI and the v0.7.0 release run until 3dba003
dropped the target again.

Dropping it was the right call -- black has no py314-specific formatting rules,
so the target changes no output and only decides whether black agrees to start.
But the pressure to re-add it comes back every time someone reads the 3.14 test
matrix and notices the formatter does not mention 3.14, and the failure it
causes does not look like its cause.

So this is the tie between the two settings, which nothing else provides:
re-adding ``py314`` is fine the moment the pin allows a black that knows it,
and fails loudly before then. Same kind of check as ``test_version.py`` --
two places in one file that have to agree, asserted rather than remembered.

The target list is read with a regex rather than a TOML parser to match
``test_version.py``'s reasoning about dev-dependency weight; ``target-version``
is unambiguous in this file.
"""

from __future__ import annotations

import pathlib
import re
from typing import Set

import pytest


def _declared_target_versions() -> Set[str]:
    root = pathlib.Path(__file__).resolve().parents[1]
    text = (root / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r"^target-version = \[([^\]]*)\]", text, re.MULTILINE)
    assert match is not None, "no `target-version = [...]` line found in pyproject.toml"
    return set(re.findall(r"['\"]([^'\"]+)['\"]", match.group(1)))


def test_pyproject_declares_target_versions() -> None:
    assert _declared_target_versions()


def test_every_declared_target_is_supported_by_the_pinned_black() -> None:
    """The assertion that would have caught the broken lint job.

    A target black does not know is not a style disagreement -- it stops the
    formatter from starting, so every file goes unchecked.
    """
    black_mode = pytest.importorskip("black.mode")
    supported = {version.name.lower() for version in black_mode.TargetVersion}

    unsupported = _declared_target_versions() - supported
    assert not unsupported, (
        f"[tool.black] target-version declares {sorted(unsupported)}, which the "
        f"installed black does not support. black exits with a usage error rather "
        f"than skipping the unknown target, so `make lint` checks nothing. Raise the "
        f"black pin in the dev extra, or drop the target."
    )
