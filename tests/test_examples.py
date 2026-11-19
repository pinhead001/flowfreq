"""The shipped examples and the clean script.

Two examples imported names that no longer existed (``USGSGage``, ``B17CEngine`` from the
package root) and failed on their first line. Running them needs the network, so this
checks what can be checked without it: each parses, and every name it imports from
flowfreq exists.
"""

from __future__ import annotations

import ast
import importlib
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = sorted((ROOT / "examples").glob("*.py"))


def test_there_are_examples_to_check():
    assert len(EXAMPLES) >= 3


@pytest.mark.parametrize("path", EXAMPLES, ids=lambda p: p.name)
def test_an_examples_flowfreq_imports_resolve(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    checked = 0
    for node in ast.walk(tree):
        if not isinstance(node, ast.ImportFrom) or not (node.module or "").startswith("flowfreq"):
            continue
        module = importlib.import_module(node.module)
        for alias in node.names:
            assert hasattr(
                module, alias.name
            ), f"{path.name} line {node.lineno}: {node.module} has no {alias.name}"
            checked += 1
    assert checked, f"{path.name} imports nothing from flowfreq"


def _clean():
    spec = importlib.util.spec_from_file_location(
        "flowfreq_tools_clean", ROOT / "tools" / "clean.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_clean_finds_artifacts_and_leaves_the_vendored_tree(tmp_path):
    clean = _clean()
    for made in (
        "build_fortran/mbuild/x.o",
        "flowfreq/peakfqr/_emafort.cp312-win_amd64.pyd",
        ".pytest_cache/v",
        "build/lib/a.py",
        "flowfreq/__pycache__/a.pyc",
        "vendor/peakfqr/__pycache__/a.pyc",
        ".venv/lib/__pycache__/a.pyc",
    ):
        target = tmp_path / made
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("x", encoding="utf-8")
    plain = {p.relative_to(tmp_path).as_posix() for p in clean.targets(tmp_path, False)}
    assert plain == {
        "build_fortran/mbuild",
        "flowfreq/peakfqr/_emafort.cp312-win_amd64.pyd",
        ".pytest_cache",
    }
    everything = {p.relative_to(tmp_path).as_posix() for p in clean.targets(tmp_path, True)}
    assert everything == plain | {"build", "flowfreq/__pycache__"}
