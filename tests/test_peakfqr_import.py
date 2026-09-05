"""The guidance flowfreq.peakfqr gives when its extension is not built.

That is the normal state for anyone who pip-installed flowfreq: the .so is
gitignored and built on demand. The bare ImportError named a private submodule
and said nothing about how to get one, so this pins the replacement.
"""

from __future__ import annotations

import builtins
import importlib
import sys

import pytest


@pytest.fixture
def without_extension(monkeypatch):
    """Import flowfreq.peakfqr with _emafort unavailable, however it is built."""
    for name in [n for n in sys.modules if n.startswith("flowfreq.peakfqr")]:
        monkeypatch.delitem(sys.modules, name, raising=False)

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "flowfreq.peakfqr._emafort" or name.endswith("._emafort"):
            raise ImportError("No module named 'flowfreq.peakfqr._emafort'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    return None


class TestMissingExtensionMessage:
    def test_raises_import_error(self, without_extension):
        with pytest.raises(ImportError):
            importlib.import_module("flowfreq.peakfqr")

    def test_names_the_build_command(self, without_extension):
        """A user hitting this needs the command, not the missing symbol."""
        with pytest.raises(ImportError) as exc:
            importlib.import_module("flowfreq.peakfqr")
        assert "build_fortran/build.py" in str(exc.value)

    def test_names_the_toolchain(self, without_extension):
        with pytest.raises(ImportError) as exc:
            importlib.import_module("flowfreq.peakfqr")
        assert "gfortran" in str(exc.value) and "meson" in str(exc.value)

    def test_says_the_rest_of_the_library_is_unaffected(self, without_extension):
        """Otherwise this reads like flowfreq itself is broken."""
        with pytest.raises(ImportError) as exc:
            importlib.import_module("flowfreq.peakfqr")
        assert "flowfreq.bulletin17c" in str(exc.value)


def test_docstring_points_at_a_path_that_exists():
    """It used to cite _shared/peakfqr/, which is not in this repository."""
    import pathlib

    src = pathlib.Path("flowfreq/peakfqr/__init__.py").read_text()
    assert "_shared/peakfqr" not in src
    assert "vendor/peakfqr/src/emafit.f" in src
