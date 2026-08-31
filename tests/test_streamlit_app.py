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
import pathlib

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


class TestPilfOverrideWiring:
    """The PILF override was reachable only from Python; now it is in the UI.

    run_ffa has always accepted low_outlier_threshold_override, but
    streamlit_app.py never passed it. These check the wiring end to end: the
    sidebar value reaches run_ffa, and the applied cut is drawn on the record
    rather than only changing a number in the parameter table.
    """

    @staticmethod
    def _peaks():
        import pandas as pd

        from tests.fixtures.big_sandy import SYSTEMATIC_PEAKS

        return pd.DataFrame(
            {
                "water_year": sorted(SYSTEMATIC_PEAKS),
                "peak_flow_cfs": [SYSTEMATIC_PEAKS[y] for y in sorted(SYSTEMATIC_PEAKS)],
            }
        )

    def test_app_passes_the_override_to_run_ffa(self, app_module):
        """Both call sites, not just the download path."""
        source = pathlib.Path(app_module.__file__).read_text()
        assert source.count("low_outlier_threshold_override=pilf_override or None") == 2

    def test_override_participates_in_the_refit_trigger(self, app_module):
        """Changing it without refitting would leave the curve silently stale."""
        source = pathlib.Path(app_module.__file__).read_text()
        start = source.index("current_ffa_inputs = (")
        assert "pilf_override" in source[start : start + 300]

    def test_plot_marks_censored_peaks_when_a_threshold_applies(self, app_module):
        fig = app_module.plot_peak_timeseries(
            self._peaks(), "Big Sandy", "03606500", pilf_threshold=2000.0, pilf_source="override"
        )
        labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
        assert any("censored" in label for label in labels)
        assert any("2,000 cfs (override)" in label for label in labels)

    def test_plot_is_unchanged_when_no_threshold_applies(self, app_module):
        fig = app_module.plot_peak_timeseries(self._peaks(), "Big Sandy", "03606500")
        assert fig.axes[0].get_legend() is None

    def test_a_threshold_below_every_peak_censors_nothing(self, app_module):
        """Guard against drawing an empty 'censored' series and a stray legend."""
        fig = app_module.plot_peak_timeseries(
            self._peaks(), "Big Sandy", "03606500", pilf_threshold=1.0
        )
        assert fig.axes[0].get_legend() is None
