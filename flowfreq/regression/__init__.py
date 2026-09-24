"""Regional regression equations for flood quantiles at ungaged sites.

See ``docs/MASTER_ROADMAP.md`` §3. This package is the offline evaluator.
:mod:`flowfreq.streamstats` is the online (NSS) one.
"""

from flowfreq.regression.equations import (
    Citation,
    OutOfRangeError,
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
    load_state,
)

__all__ = [
    "Citation",
    "EquationsUnavailable",
    "JURISDICTIONS",
    "Jurisdiction",
    "OutOfRangeError",
    "RegressionEquation",
    "RegressionEstimate",
    "StateLibrary",
    "Variable",
    "WAVE_NAMES",
    "available_states",
    "evaluate",
    "evaluate_weighted",
    "jurisdiction",
    "load_state",
    "wave",
]
