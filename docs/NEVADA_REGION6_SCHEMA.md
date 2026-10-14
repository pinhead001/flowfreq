# Nevada region 6 and the hybrid regions: a schema proposal

Issue #104. Design only; nothing here is implemented yet.

`NV.json` is `partial` for two things WSP 2433 publishes that the equation schema
(`flowfreq/regression/equations.py`) cannot hold:

1. **Region 6's 2-year equation is printed `Q=0`** (Table 10). A log-linear equation
   `log10 Q = b0 + Σ b_i T_i(X_i)` cannot give zero.
2. **The hybrid regions (6 and 10) print error statistics that are not standard errors of
   prediction.** The report says they "include much of the within-station residual variance"
   and are "not comparable" to a regression standard error. `sep_log` is the only error field
   that could hold them, and putting them there would be wrong.

The `verified` definition in `flowfreq/regression/library.py` is: "everything the cited
reports publish for the stored regions is stored, double-entered against the source, and
cross-checked". NV can only become `verified` once both items are stored.

**Recommendation (details in §5).** Add an explicit `form` to `RegressionEquation`, with
`"constant"` next to the default `"loglinear"`, and store region 6's 2-year equation as
`{"form": "constant", "constant_cfs": 0}`. Add a separate `published_errors` list, which the
evaluator never reads, for statistics the report itself says are not prediction errors. Both
changes are additive, so every other state file stays byte-identical.

---

## 1. What WSP 2433 publishes

Source: Thomas, Hjalmarson and Waltemeyer, 1997, USGS Water-Supply Paper 2433,
<https://pubs.usgs.gov/wsp/2433/report.pdf>. The page numbers below are the printed page
numbers; the PDF page is 9 higher. Each was checked on the rendered page.

| Item | Where | Content |
|---|---|---|
| Region 6 equations | Table 10, p. 46 | 2-yr **`Q=0`**; 5-100 yr `Q = a·AREA^b·(ELEV/1,000)^c`. "Data were based on 80 stations. Dashes indicate no data." |
| Region 6 errors | Table 10, p. 46, column "Estimated average standard error of regression, in log units" | 2-yr `--`; 5, 10, 25, 50, 100 yr: **1.47, 1.12, .796, 1.10, 1.84** |
| Region 6 equivalent years | Table 10, p. 46 | `--`, 0.233, .748, 2.52, 1.75, .794 |
| Region 10 equations | Table 14, p. 54 | `Q = a·AREA^b` at 2-100 yr, 104 stations |
| Region 10 errors | Table 14, p. 54, same column heading | 2, 5, 10, 25, 50, 100 yr: **1.14, .602, .675, .949, .928, 1.23** |
| Region 10 equivalent years | Table 14, p. 54 | 0.618, 3.13, 3.45, 2.49, 3.22, 2.22 |
| Table headnote (both) | pp. 46, 54 | "Estimated average standard error of regression for the hybrid method includes much of the within-station residual variance and therefore is not comparable to standard error of estimate from an ordinary least-squares regression." |
| Hybrid method | "Hybrid Method", pp. 95-97 | Station records are composited by basin-characteristic level and plotted with Weibull plotting positions; the method "addresses the problem of no-flow years" (p. 97). |
| How the error is computed | p. 99 | A "direct assessment procedure": random subsets of station records, the sample variance of the T-year estimates, 30 subsets per step. Details are in Hjalmarson and Thomas (1992). |
| How to read the error | pp. 99-100 | The relations "must be considered as only rough estimates. The magnitude of the errors is too large to have much meaning, but it does show the relative differences between relations for the different recurrence intervals." |
| Transition zones | p. 19 (eqs. 6-7), p. 100; examples pp. 39-44 | Weighted averages of flows. Already in `flowfreq/regression/nevada.py`. |

The report never says why region 6's 2-year flood is zero. The likely reason, inferred and not
stated, is the hybrid method's handling of no-flow years (p. 97) combined with the high share of
undefined station relations in regions 6 and 10 ("about 55 percent", p. 99). If more than half
of the composited annual peaks are zero, the median peak is zero. Read this way, `Q=0` is a
genuine published estimate, not a missing value. NSS treats it the same way: its region GC164
`PK50AEP` equation is `0*(0*(DRNAREA))^(0)*(0*(ELEV))^(0)` and returns 0
(`data/nss_snapshots/NV_2026-10-01.json`).

**The unit "log units" is not pinned down.** On p. 99 the average station standard deviation is
"about 0.9 log units for Regions 6 and 10 ... and 0.4 log units for the 12 GLS regions", which
suggests log10. The same page's ratios ("more than 10 for Region 6, about 6 for Region 10")
cannot be reproduced from Table 29's GLS errors either as ratios of standard errors or as ratios
of variances. The proposal therefore stores the values as printed, with the printed unit
string, and converts nothing.

**NSS carries some of the hybrid errors, under the wrong code.** The snapshot lists them under
the same `SE` code as the GLS regions' *percent* errors (region 1: `SE 59.0`; region 6:
`SE 1.47`). It also has gaps:

- region 6 has 4 of 5 (it is missing 0.796 at 25 years);
- region 10 has 2 of 6 (1.14 and 1.23 only).

So NSS can cross-check 6 of the 11 values; the other 5 must be double-entered from the
rendered page. The snapshot's `nss_sep_note` already refuses to map NSS's `sep` to `sep_log`
for the same reason.

## 2. Current state and who reads it

- `NV.json` has 35 equations. Region 6 has AEPs 0.2-0.01 and no 0.5. All 11 hybrid equations
  have `sep_log = None`, so `evaluate()` gives them `interval = None` and
  `interval_method = None`.
- `nevada.estimate_area_weighted(0.5, {"6": ..., "10": ...})` raises `EquationsUnavailable`
  (asserted in `tests/test_regression_nevada.py`). By contrast, the report's eq. 6 and NSS's
  area averaging give `Q2 = Q2_region10 · A10 / A`.
- `tests/test_regression_nevada.py::test_matches_live_nss` skips region 6 `PK50AEP` with
  `assert est.value == 0`.

These code paths assume every equation is log-linear with a finite, positive flow:

| Consumer | Assumption | What breaks on a zero flow |
|---|---|---|
| `equations.evaluate` | `flow = bias · 10**log_q` | nothing, if handled explicitly |
| `equations.evaluate_weighted(space="log")` | `math.log10(e.flow_cfs)` | `ValueError: math domain error` with no context |
| `nevada._weighted` | `space="linear"` | nothing: `0·w` is correct for eq. 6 and eq. 7 |
| `validation/wave_qa.py` (`estimate_state`, `discrepancy`, l. 563-579, 767, 779) | reads `eq.sep_log` directly; compares `log_flow` | `log_flow = -inf` gives an infinite difference and a spurious flag |
| `future_flow.py` l. 736-743 | `change_factor = q_proj / q_hist` | `ZeroDivisionError` |
| `tests/test_regression_equations.py::test_packaged_files_evaluate_as_before` | recomputes `intercept + Σ b·T(x)` for every packaged equation | needs to skip non-log-linear forms |
| `tools/gen_regression_coverage.py` | counts equations | NV goes from 35 to 36; regenerate `docs/REGRESSION_COVERAGE.md` |

`wave_qa` reading `eq.sep_log` directly is the strongest argument against storing the hybrid
errors in `sep_log`, even with a flag. One consumer that forgets to check the flag turns 1.47
log units into a "standard error" a factor of about 30 wide, and that wrongly loosens the
cross-border QA threshold. As things stand, the hybrid regions' missing SEP makes that
threshold stricter (`docs/WAVE_QA_W1_W2.md`).

## 3. Problem A: a constant equation

### A1. An explicit `form` with a `"constant"` variant (recommended)

JSON (the shared `region_name` and `citation` are omitted here):

```json
{
  "region_code": "6", "aep": 0.5,
  "form": "constant", "constant_cfs": 0.0,
  "variables": [
    {"code": "DRNAREA", "transform": "log10", "minimum": ..., "maximum": ...},
    {"code": "ELEV", "transform": "log10", "scale": 0.001, "minimum": ..., "maximum": ...}
  ],
  "n_sites": 80
}
```

- `RegressionEquation` gains `form: str = "loglinear"` and `constant_cfs: Optional[float] = None`.
  `intercept` becomes `Optional[float]`: required for `loglinear`, and **absent** for
  `constant`.
- In the constant form, `variables` set the range check only (the region's calibrated ranges,
  the same ones its 5-year equation uses). `coefficients` is omitted and stored internally as
  zeros, so `len(variables) == len(coefficients)` still holds. `__post_init__` rejects these
  for a constant equation: a negative or non-finite constant; any non-zero coefficient;
  `sep_log`, `avp`, `model_error_variance` or `covariance`; a `bias_correction` other than 1.
- `evaluate()` short-circuits after the missing-characteristic check and the range check:
  `flow_cfs = constant_cfs`; `log_flow = log10(constant_cfs)`, which is `-inf` for 0;
  `interval = None`; `interval_method = "none: printed constant (Q=0)"`.

Trade-offs:

- **Backward compatibility.** Additive. Every existing file has no `form` key and loads
  unchanged. An older flowfreq given the new `NV.json` through `data_dir` fails loudly with
  `KeyError: 'intercept'` instead of silently misreading it, because the constant form omits
  `intercept`. That is the reason for omitting it rather than storing a dummy. `SCHEMA_VERSION`
  can stay 1 (but see decision D3).
- **Evaluator.** One branch in `evaluate`. The intercept-is-required check moves into
  `__post_init__`.
- **Callers.** `lib.aeps("6")` returns all six AEPs, and `lib.equation("6", 0.5)` and
  `evaluate()` work as for any other equation. Eq. 6 and eq. 7 at AEP 0.5 now return the
  report's weighted value instead of raising. `log_flow = -inf` is honest, but it is a value
  some consumers will not expect. `evaluate_weighted(space="log")` must raise a clear
  `ValueError` naming the region when a component is zero, because log-space weighting of a
  zero has no meaning. `future_flow` must refuse a zero historical flow explicitly.
- **Typing.** Making `intercept` `Optional` touches every reader of `eq.intercept`. In library
  code there are none outside `equations.py` (`nss.py` reads its own parsed-equation
  intercept). A few tests read it (AZ, ID, provisional). mypy enforces `flowfreq/regression/`
  (it is not on the `[tool.mypy]` override list), so each reader needs a narrowing check.

### A2. A sentinel inside the log-linear form

Store `"intercept": null`, read as `-inf`, with no variables. The evaluator already gives
`10**-inf = 0`, so there are no evaluator changes.

- **For:** the smallest diff.
- **Against:** the intent is opaque. `null` reads as "missing", not "zero". Strict JSON cannot
  carry `-Infinity`, so a mapping convention is needed. There is no range check, because there
  are no variables. Every consumer has to detect `isinf(intercept)` for itself. It also widens
  the meaning of `intercept` for every state just to cover one cell.

### A3. Keep it out of the equation list

Add a library-level `printed_values: [{"region_code": "6", "aep": 0.5, "value_cfs": 0,
"citation": ...}]`. Raise a subclass of `EquationsUnavailable` that carries the printed value
from `lib.equation("6", 0.5)`. `nevada.py`'s transition functions catch it and use the value.

- **For:** "every equation is log-linear" stays true. Generic consumers (`wave_qa`,
  `future_flow`) need no changes, because they keep getting a refusal.
- **Against:** a new top-level schema field and a second lookup path. `lib.aeps("6")` either
  still omits 0.5 or has to merge two sources. A generic caller of
  `evaluate(lib.equation(...))` sees an exception for a number the report prints, which is
  arguably still "not stored" in the sense the `verified` definition means.

| | A1 `form: constant` | A2 sentinel | A3 printed values |
|---|---|---|---|
| Existing JSON unchanged | yes | yes | yes |
| Old reader on new NV.json | loud `KeyError` | silent `-inf` intercept, or `TypeError` | silently ignores |
| Range check on the Q=0 cell | yes | no | no |
| Generic `evaluate()` returns 0 | yes | yes | no (raises) |
| Consumers needing guards | 3 (`evaluate_weighted` log, `wave_qa`, `future_flow`) | same 3, plus `isinf` checks | 0 generic, 1 Nevada |
| Extends to other odd forms later | yes (the `form` discriminator) | no | no |

## 4. Problem B: error statistics that are not prediction errors

### B1. A `published_errors` list the evaluator never reads (recommended)

```json
"published_errors": [
  {"kind": "hybrid_standard_error_of_regression", "value": 1.47, "units": "log units",
   "comparable_to_sep": false, "source": "Table 10, p. 46"}
]
```

- `RegressionEquation.published_errors: Tuple[PublishedError, ...] = ()`. `PublishedError` is a
  frozen dataclass with `kind: str`, `value: float`, `units: str` (required, non-empty),
  `comparable_to_sep: bool = False` and `source: str = ""`. `kind` comes from a closed set,
  starting with `{"hybrid_standard_error_of_regression"}`, so a typo fails at load time.
- `_prediction_sd` never reads it, so hybrid equations keep `interval = None`. `evaluate` sets
  `interval_method = "none: published error is not a prediction error (WSP 2433 p. 99)"`
  whenever an equation has no comparable statistic but does have `published_errors`. Callers
  then see *why* there is no interval instead of a bare `None`.
- Region 6's 2-year constant has no entry, because the report prints `--`.

Trade-offs:

- **Backward compatibility.** Additive, and omitted when empty. An older reader ignores the key
  (`from_dict` reads only the keys it knows). That is harmless, because the values change no
  number.
- **Evaluator and prediction intervals.** No behaviour change for any state.
- **Callers.** The values become machine-readable (for reports and the app's display) without
  any risk of being used as a SEP. `wave_qa` keeps treating the hybrid regions as having no SEP.
- **Reuse.** The same list fits the other published-but-unschematized statistics the
  verified files mention only in their notes. Examples are WY's 95-percent prediction-interval
  factors and the equivalent years of record that WY, OR and NV all print. Those kinds would be
  added only on request (decision D2).

### B2. Store them in `sep_log` with a kind flag

Add `sep_kind: "prediction" | "hybrid_regression"` and make `_prediction_sd` skip non-prediction
kinds.

- **For:** one field, and the obvious place to look.
- **Against:** `wave_qa` already reads `eq.sep_log` without going through `_prediction_sd`, and
  any future consumer can do the same. A forgotten check silently produces 30x-wide intervals
  or loosens QA thresholds. A field whose meaning depends on a second field is easy to misread.

### B3. Notes only, with the `verified` definition amended

Leave the values in `NV.json` `notes` (where they already are, verbatim). Amend `library.py`'s
definition so that a statistic the report itself disclaims as not a prediction error may be
recorded in notes rather than in fields.

- **For:** no code changes.
- **Against:** it lowers the bar for every state to fix one, and it leaves the values readable
  only by humans. Every other `verified` file stores what it publishes in fields, or names a
  missing schema field in its notes (WY). A notes-only route would make "verified" mean less
  for NV than for WY.

| | B1 `published_errors` | B2 `sep_kind` | B3 notes only |
|---|---|---|---|
| Schema change | additive list | additive flag on an existing field | none |
| Risk of misuse as a SEP | none (separate field) | high (`wave_qa` reads `sep_log` today) | none |
| Machine-readable | yes | yes | no |
| Changes the `verified` bar | no | no | yes |

## 5. Recommendation: A1 + B1

- **Schema** (`equations.py`). Add `form` and `constant_cfs` (A1), and `PublishedError` and
  `published_errors` (B1). `intercept` becomes `Optional[float]`, required for `loglinear`.
  `from_dict` reads the new keys, and its defaults reproduce today's objects exactly.
- **Evaluator.**
  - The constant short-circuit in `evaluate`.
  - The `interval_method` explanation.
  - `evaluate_weighted(space="log")` raises `ValueError("... region 6 AEP 0.5 is a printed
    constant 0; weight in linear space")` on a zero component.
- **Nevada.**
  - `NV.json` gains region 6 AEP 0.5 (36 equations) with region 6's DRNAREA and ELEV ranges.
  - It gains 11 `published_errors` entries, from the Table 10 and Table 14 values in §1.
  - `notes` is updated, and `status` changes to `verified` once the cross-checks in §7 pass.
  - In `nevada.py`, the `estimate_area_weighted` docstring drops "is not stored". At AEP 0.5,
    eq. 6 with region 6 returns `Q2_10 · A10/A`, and eq. 7 with `low_region="6"` returns
    `w_high · Q2_1`.
- **Consumers.**
  - `wave_qa` skips equations with `form == "constant"` in its discrepancy screen and reports
    them as not evaluable, the way it already treats a missing SEP.
  - `future_flow` raises a `ValueError` on a constant (or zero) historical flow; see D4.
- **What a caller sees.** `evaluate(load_state("NV").equation("6", 0.5), basin)` returns
  `flow_cfs == 0.0` with `interval is None` and an `interval_method` string saying why. For
  hybrid equations, `eq.published_errors` lists the printed values while `sep_log` stays
  `None`. No other state's numbers change.

## 6. Files affected

| File | Change |
|---|---|
| `flowfreq/regression/equations.py` | `form`, `constant_cfs`, `PublishedError`, `published_errors`, Optional `intercept`, `evaluate` branch, `evaluate_weighted` zero guard, module docstring |
| `flowfreq/regression/__init__.py` | export `PublishedError` |
| `flowfreq/regression/library.py` | docstring only, noting that `aeps()` includes constant-form cells; the status text is unchanged under A1 + B1 |
| `flowfreq/regression/nevada.py` | docstrings (module, `estimate_area_weighted`) |
| `flowfreq/data/regression/NV.json` | +1 equation, +11 `published_errors`, notes, status |
| `flowfreq/validation/wave_qa.py` | skip constant-form equations |
| `flowfreq/future_flow.py` | refuse a zero historical flow |
| `tools/snapshot_nss.py` (optional) | recognise NSS's `0*(...)^(0)` as the constant form instead of a `parse_error` |
| `docs/REGRESSION_COVERAGE.md` | regenerate with `tools/gen_regression_coverage.py` |
| `TODO.md`, `CHANGELOG.md`, `docs/MASTER_ROADMAP.md` | at implementation time, not in this PR |

## 7. Tests

`tests/test_regression_equations.py`:

- The constant form round-trips through `from_dict`. A file without `form` gives
  `form == "loglinear"` and an object equal to today's.
- `evaluate` on a constant equation:
  - returns `flow_cfs == 0.0`, `log_flow == -inf`, `interval is None` and the
    `interval_method` text;
  - raises `OutOfRangeError` outside the variables' range, and with `allow_extrapolation=True`
    computes the value and records the violation;
  - raises `KeyError` on a missing characteristic.
- Rejected at construction: a constant with a negative or NaN value, a non-zero coefficient,
  `sep_log`, `covariance`, or `bias_correction != 1`; and a log-linear equation without an
  `intercept`.
- `published_errors`: parses; `units` is required; an unknown `kind` raises; an entry never
  sets an interval (an equation with only `published_errors` has `interval is None`); an
  equation with both `sep_log` and `published_errors` still uses `sep_log`.
- `evaluate_weighted`: a zero component works in linear space; in log space it raises a
  `ValueError` that names the region.
- `test_packaged_files_evaluate_as_before` skips `form == "constant"` explicitly. A companion
  test asserts that NV region 6 AEP 0.5 is the *only* packaged constant-form equation.

`tests/test_regression_nevada.py`:

- `status == "verified"`; `lib.aeps("6") == AEPS`; 36 equations. Replace
  `test_region_6_two_year_is_not_stored` with a test that evaluates it to 0 at `BASIN`.
- The hybrid published errors equal the Table 10 and Table 14 literals copied from the report
  (as the existing tests copy them). They are `comparable_to_sep is False`, `sep_log` stays
  `None`, and `evaluate(...).interval is None` for all 11 equations.
- `estimate_area_weighted(0.5, {"6": 21.0, "10": 36.0}, ...)` equals the region 10 Q2 at the
  whole-basin characteristics times 36/57, replacing the `EquationsUnavailable` assertion in
  `test_errors`. Eq. 7 with `low_region="6"` at AEP 0.5 gives `w_high · Q2_1`.
- `test_matches_live_nss` (`requires_network`) evaluates region 6 `PK50AEP` through the
  library instead of skipping it (36 compared). Where NSS has an `SE` for a hybrid cell,
  it equals the stored `published_errors` value (6 of 11 cells).

Elsewhere:

- `tests/test_wave_qa.py`: a constant-form equation is excluded from the discrepancy screen and
  not flagged.
- `tests/test_future_flow.py`: a change factor on a constant-zero equation raises `ValueError`.

## 8. Decisions for the maintainer

- **D1. Range-check the `Q=0` cell?** The proposal says yes, using region 6's DRNAREA and ELEV
  ranges, so that a site outside region 6's calibration is refused at every AEP and not only
  at 5-100 years. NSS also range-checks the cell (its GC164 parameters carry limits). The
  alternative is a constant with no variables, which is accepted anywhere.
- **D2. Equivalent years of record.** All six NV tables print them, and so do the WY and OR
  reports, and no state stores them. WY is `verified` without them, so the precedent is that
  they are outside the `verified` bar. Should `published_errors` (or a sibling list) take them
  now, or should the precedent be written into `library.py`'s status definition?
- **D3. `SCHEMA_VERSION`.** It can stay 1, because the change is additive and an old reader
  fails loudly on the constant form. Bumping to 2 would mean `parse_state` accepting {1, 2}.
- **D4. `future_flow` on a zero flow:** raise (as proposed) or return `change_factor = None`?
- **D5. Is `log_flow = -inf` acceptable** on `RegressionEstimate`, or should a constant-zero
  estimate carry `log_flow = None`? `None` would make the field `Optional` for every caller.
- **D6. Release line.** This is a schema addition, so it would land on `main` (0.12.0). It
  would not be backported to `release/0.11.x`.
