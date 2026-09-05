# Design: the vendored Fortran as a selectable analysis engine

Letting a caller run the same record through the vendored USGS `emafitpr`
instead of, or alongside, the native Python EMA — and compare the two.

**Status: design only.** Nothing here is implemented. This document exists so
the work can start from a specification rather than from a conversation, and so
the one part with a silent failure mode is written down before anyone builds it.

---

## 1. Why

The motivating case is a CLOMR/LOMR submittal. A reviewer can challenge the
analysis, and "this reproduces USGS PeakFQ 8.1.0 to within 0.04% at Q100 on
this gage" is evidence in a way that "the library has a test suite" is not.

Be clear about what it does and does not add. Parity is already established: the
committed goldens under `tests/fortran_parity/golden/` show agreement to 1e-8 on
the moments and better than 0.11% on quantiles for four reference gages. Running
both engines on an ordinary record mostly re-proves that.

The value is in the records where the two might genuinely diverge — heavy
censoring, historic peaks, extreme skew, short records. The standing
`xfail(strict=True)` on Cains Coulee's `skew_weighted` is a real, unexplained
0.058-skew-unit gap that exists today. It is not known to be a port defect —
the investigation in `TestCainsCouleeAsGMseDiscrepancy` traces it to
`emafitpr`'s own internally-reported `as_G_mse` (0.2212) disagreeing with a
standalone call to the same `mseg_all` routine on identical inputs (0.0749) —
which is precisely the point: a disagreement worth seeing is not the same as a
verdict about which side is wrong. A per-run comparison surfaces that class of
case on a user's own data, which nothing currently does.

It also catches defects in the native path that no unit test would. The biased
station-skew estimator in `B17CEngine.fit` — see the strict xfail in
`tests/test_engine.py` — sat undetected in a public API for exactly as long as
nothing compared it against a reference.

---

## 2. What already exists

Most of the bridge is built and works.

| Piece | Where | State |
|---|---|---|
| f2py extension | `flowfreq/peakfqr/`, built by `build_fortran/build.py` | Works. 13 s to build, 812 KB `.so` |
| Signature file | `build_fortran/_emafort.pyf` | Works. **Must** be passed or f2py wraps QUADPACK's `dqag` and the build fails |
| Fortran call + output unpacking | `ReferenceResult.from_emafit()` | Works. Returns parameters, quantiles, CIs, variance |
| Interval construction | `tests/fortran_parity/cases.py::build_emafit_inputs` | **Test code.** Handles the parity cases, not the full input space |
| Golden capture | `tools/gen_fortran_golden.py` | Works |

The gap is one item: there is no library-side, supported translation from a
`Bulletin17C` input set to `emafitpr` arguments. That is the whole of the risk.

---

## 3. The interval builder, and why it is the hard part

`emafitpr` does not take peaks and water years. It takes five parallel arrays
describing, for every row, a flow interval and the perception thresholds it was
observed against:

```
ql, qu    flow interval bounds, log10   (ql == qu for an exact observation)
tl, tu    perception threshold bounds, log10
dtype     1 only for a peak carrying the USGS historic flag; 0 otherwise
```

Construction follows `siteQT` in `vendor/peakfqr/R/readInputs.R`. `Bulletin17C`
by contrast accepts `peak_flows`, `water_years`, `historical_peaks`,
`perception_thresholds`, `user_low_outlier_threshold` and `ema_params`, and
derives the rest internally. Something has to translate, and the translation is
where meaning is assigned.

### The failure mode

A wrong translation does not raise. It produces a **different, valid analysis**,
and the comparison then reports disagreement between two implementations that
are both correct. Inside a LOMR package that is worse than not offering the
feature: it manufactures doubt about numbers that were fine.

### The worked example

USGS 12363000 has 98 peaks spanning 102 water years; 1924–1927 are unmeasured.
Two defensible-looking constructions:

| Construction | rows | at-site skew | Q100 |
|---|---|---|---|
| 98 observations, gap years omitted | 98 | **+0.435** | **120,064** |
| 102 rows, gap years censored at the lower threshold | 102 | +0.250 | 119,473 |

Neither errors. The first is right *here*, because no perception threshold
asserts what would have been recorded in those years — they are simply
unmeasured, and `Bulletin17C` given peaks and water years fits the 98. The
second would be right if a threshold had been declared, because then the absence
of a peak is itself an observation.

`ParityCase.fill_missing_years` exists for exactly this distinction, and
`tests/fixtures/site_12363000.py` documents it. **The library-side builder needs
the same switch, driven by whether the caller supplied perception thresholds
covering the gap.**

### Cases the builder must handle

`build_emafit_inputs` today handles systematic peaks, historic peaks and
threshold-fill. A supported builder additionally has to decide, and be tested
on, each of:

1. **Gap years with no threshold declared** — omit (the 12363000 case above)
2. **Gap years inside a declared threshold period** — censor at the lower bound
3. **Historic peaks** — `dtype = 1`, and only for the USGS historic flag (peak
   code 7), not for every peak inside the historical period
4. **Zero flows** — `Bulletin17C` tracks `n_zeros`; `siteQT` and `QMIN = 1e-20`
   need checking against what the native path does with them
5. **A user-supplied PILF threshold** — `gbthrsh0 > -6` uses the value, `<= -6`
   runs MGBT. `run_ffa`'s `low_outlier_threshold_override` must map onto this
6. **Multiple perception threshold periods** — the dict is keyed
   `(start, end) -> lower`; overlapping or adjacent periods need a defined rule
7. **Records where `water_years` is omitted** — `Bulletin17C` allows it

Items 4, 6 and 7 have no precedent in the parity cases. They need answering
against `readInputs.R`, not invented.

### Encodings to get right

From `CLAUDE.md`, restated because they are easy to invert:

- **Skew MSE:** `0` = generalized, no error; `< 0` = generalized with
  `MSE = -value`; `> 0` = weighted; `> 1e10` = station-only
- **MGBT:** `gbthrsh0 <= -6` computes MGBT; `> -6` is used as the threshold
- **Weight options:** 1 = HWN, 2 = ERL, 3 = INV
- **Everything is log10** at the `emafitpr` boundary
- `QMIN = 1e-20`, `QMAX = 1e20` from `peakfqr/R/main.R`

---

## 4. Result mapping

`ReferenceResult` and `FrequencyResults` overlap in five field *names*, and even
those differ in shape — `quantiles` is `Dict[float, float]` on one and a
`DataFrame` on the other. The adapter is real work, not a rename.

Available from `emafitpr`, whose output groups are
`n`, `mgbt`, `skew`, `cmoms`, `quantiles`:

| `FrequencyResults` field | Source |
|---|---|
| `mean_log`, `std_log`, `skew_station`, `skew_weighted` | `parameters` dict via `cmoms_to_parameters` (`skew_at_site` is the station skew) |
| `n_peaks`, `n_systematic`, `n_historical` | direct |
| `low_outlier_threshold` | `10 ** mgbt.gbval`, guarded by the `> -6` sentinel |
| `n_low_outliers` | `mgbt.gbnlow` |
| `n_zeros` | `mgbt.gbnzero` |
| `quantiles`, `confidence_limits` | `quantiles.yp / ci_low / ci_high`, log10 → cfs, into DataFrames matching the native column names |
| `skew_used_mse` | `skew.as_G_mse_o` |
| `method` | constant `AnalysisMethod.EMA` |

Not reported by `emafitpr`, so the adapter must decide rather than guess:

| Field | Options |
|---|---|
| `mgb_critical_value` | Not in the outputs. Either recompute natively, or leave `None` and let the comparison skip it |
| `ema_iterations`, `ema_converged` | Not reported. Leave `None`; a `False` would be a fabrication |
| `pilf_flows` | Derivable from the input peaks and the returned threshold, not returned directly |
| `n_censored` | Count from the constructed intervals, not from the Fortran |
| `skew_regional`, `skew_used` | Echo the caller's input; `emafitpr` does not return them |

**Rule: never synthesise a field the Fortran did not report.** `None` that the
comparison skips is honest; a plausible-looking default is not, and this output
may end up in a submittal.

---

## 5. Proposed API

Three layers, each usable without the ones above it.

```python
# 1. library
b = Bulletin17C(peak_flows=..., water_years=...)
b.run_analysis(method="ema", engine="native")    # default, unchanged
b.run_analysis(method="ema", engine="fortran")   # requires the built extension

# 2. comparison — the thing a submittal actually wants
from flowfreq.workflow import compare_engines
report = compare_engines(peak_flows=..., water_years=..., aeps=[...])
report.max_quantile_deviation_pct   # 0.106 for 12363000 (at the 500-yr)
report.to_markdown()                # a table that can go in an appendix

# 3. CLI
flowfreq compare --site 12363000 --peaks peaks.csv
```

`engine="fortran"` raises the existing actionable `ImportError` from
`flowfreq/peakfqr/__init__.py` when the extension is absent. `engine` defaults
to `"native"` forever — the Fortran path is opt-in, never a silent substitution.

`compare_engines` is the feature. `engine=` alone makes a user run two analyses
and diff them by hand, which is the thing they wanted the library to do.

---

## 6. Distribution, and why it is deferred

The extension builds in 13 s and is 812 KB, so the cost is not the build. The
cost is *where* it can happen.

| Consumer | Feasible today? |
|---|---|
| Source checkout with gfortran + meson | Yes. `make fortran` |
| `pip install` from git | No — no compiler step, no wheel |
| Streamlit Community Cloud | Plausibly, via `packages.txt` for gfortran plus meson/ninja in requirements — untested, and a build failure there is opaque |
| `pip install flowfreq` from PyPI | Needs binary wheels |

Binary wheels mean `cibuildwheel`, `auditwheel`, and a matrix over Python
version × OS. It also converts flowfreq from a pure-Python wheel into a
platform-specific one, which is a distribution change with its own maintenance
cost.

**Recommendation: build layers 1–3 for source checkouts first.** That already
serves the LOMR case, where the analysis is run by an engineer on a workstation,
not by an anonymous visitor to a hosted app. Decide on wheels only once the
feature has been used in anger and the demand is known.

The deployed Streamlit app is explicitly **out of scope** for the first pass.

---

## 7. Correctness plan

The feature is a verification tool, so its own verification has to be stronger
than usual. Each risk with the test that retires it:

| Risk | Test |
|---|---|
| Interval builder disagrees with `siteQT` | Build inputs for all four existing parity cases via the **library** builder and assert byte-equality with `build_emafit_inputs`. The goldens then re-validate them |
| Gap years handled wrongly | 12363000 with and without a declared threshold, asserting the +0.435 / +0.250 split above |
| Historic peaks mis-flagged | Big Sandy, which has them, must reproduce its golden through the library builder |
| Zero flows | Needs a fixture; none of the current parity sites has one |
| Adapter fabricates a field | Assert `ema_converged is None` and `mgb_critical_value is None` on a Fortran result — the absence is the contract |
| Comparison reports false disagreement | On all four parity sites, `max_quantile_deviation_pct` must be under the tolerances already measured in the goldens |

**Do not add the CLI or app surface until the builder passes the first row.**
Everything else is presentation over a translation that either is or is not
faithful.

---

## 8. Estimate

| Piece | Estimate |
|---|---|
| Library interval builder + tests | 1 day — the risk lives here |
| `ReferenceResult` → `FrequencyResults` adapter | 2 h |
| `engine=` on `run_analysis` | 2 h |
| `compare_engines` + markdown output | 0.5 day |
| CLI `compare` subcommand | 2 h |
| **Total, source checkouts** | **~2 days** |
| Binary wheels, if pursued later | ~1 week |

---

## 9. Open questions

1. **Zero flows.** What does `siteQT` do with them, and does it match
   `Bulletin17C`'s `n_zeros` handling? No parity site exercises this.
2. **Overlapping perception threshold periods.** Is the dict key rule
   documented in `readInputs.R`, or is it undefined behaviour to reject?
3. **`mgb_critical_value`.** Recompute natively for the comparison, or report
   the field as unavailable on the Fortran path?
4. **Should `compare_engines` be able to run without the extension**, falling
   back to the committed goldens for the four reference sites? That would give
   the deployed app *something* honest to show without any build step.

Question 4 is worth deciding early: it is the cheap 80% of the feature, and it
changes what layer 2 should look like.
