"""Import smoke test for ``app/streamlit_app.py``.

The app is a top-level script, not a module with a ``main()``: importing it
executes the entire body -- every widget call, every ``st.session_state``
default, and every hydrolib call reachable before the first button press.
Streamlit calls this *bare mode*, widgets return their declared defaults, and
``download_data`` is therefore ``False``, so nothing here contacts NWIS.

That makes a plain import a real check, and it is coverage the app did not
have: CI linted ``hydrolib/`` and ``tests/`` only and ran no app tests at all,
which is why Streamlit changes could not be reviewed with any confidence. This
catches a stale import, a ``NameError`` on a branch no test takes, and -- the
one that matters -- an app call that no longer matches a hydrolib signature.

Skipped when Streamlit is absent, which includes the 3.9 matrix job: Streamlit
requires Python >= 3.10. The ``app`` job in ``.github/workflows/tests.yml``
installs ``app/requirements.txt`` so that this actually runs somewhere.
"""

from __future__ import annotations

import importlib
import logging

import pytest

pytest.importorskip(
    "streamlit",
    reason="Streamlit not installed; pip install -r app/requirements.txt",
)

# One 'missing ScriptRunContext!' warning per widget call, and the app builds
# dozens. Expected in bare mode, and it would bury a real failure in the log.
logging.getLogger("streamlit.runtime.scriptrunner_utils.script_run_context").setLevel(logging.ERROR)

# Names the app pulls out of hydrolib and app.ffa_*. If one of these stops
# resolving, the app is broken for every user even though hydrolib's own tests
# are green -- which is exactly the failure this module exists to catch.
WIRED_CALLABLES = (
    "run_ffa",
    "format_parameters_df",
    "format_quantile_df",
    "build_skew_curves_dict",
    "compute_skew_tables",
    "export_comparison_csv",
    "export_ffa_to_zip",
    "plot_frequency_curve_streamlit",
    "USGSgage",
    "Hydrograph",
)


@pytest.fixture(scope="module")
def app_module():
    """The imported app. Importing it *is* the smoke test; assertions follow."""
    return importlib.import_module("app.streamlit_app")


def test_app_imports(app_module):
    assert app_module.__name__ == "app.streamlit_app"


@pytest.mark.parametrize("name", WIRED_CALLABLES)
def test_wired_entry_point_resolves(app_module, name):
    attr = getattr(app_module, name, None)
    assert attr is not None, f"app/streamlit_app.py no longer exposes {name}"
    assert callable(attr), f"{name} is not callable"


def test_app_reports_the_installed_hydrolib_version(app_module):
    from hydrolib import __version__

    assert app_module.__version__ == __version__


def test_skew_options_are_offered(app_module):
    """SKEW_OPTIONS drives a sidebar selectbox; an empty one would render a dead control."""
    assert app_module.SKEW_OPTIONS
