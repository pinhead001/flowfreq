"""Regional regression equations for flood quantiles at ungaged sites.

See ``docs/MASTER_ROADMAP.md`` §3. This package is the offline evaluator.
:mod:`flowfreq.streamstats` is the online (NSS) one.
"""

from flowfreq.regression.equations import (
    Citation,
    OutOfRangeError,
    PublishedError,
    RegressionEquation,
    RegressionEstimate,
    Variable,
    evaluate,
    evaluate_weighted,
)
from flowfreq.regression.jurisdictions import (
    JURISDICTIONS,
    WAVE_NAMES,
    Jurisdiction,
    jurisdiction,
    wave,
)
from flowfreq.regression.library import (
    EquationsUnavailable,
    StateLibrary,
    available_states,
    awaiting_nss_states,
    load_state,
    nss_snapshot_states,
    provisional_states,
)
from flowfreq.regression.nss_library import asep_to_sep_log, snapshot_to_library

__all__ = [
    "Citation",
    "EquationsUnavailable",
    "JURISDICTIONS",
    "Jurisdiction",
    "OutOfRangeError",
    "PublishedError",
    "RegressionEquation",
    "RegressionEstimate",
    "StateLibrary",
    "Variable",
    "WAVE_NAMES",
    "asep_to_sep_log",
    "available_states",
    "awaiting_nss_states",
    "evaluate",
    "evaluate_weighted",
    "jurisdiction",
    "load_state",
    "nss_snapshot_states",
    "provisional_states",
    "snapshot_to_library",
    "wave",
]
