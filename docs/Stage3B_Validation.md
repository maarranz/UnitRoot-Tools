# Stage 3B: ADF and DF-GLS integration validation

27 October 2026. Stage 3B is complete and stops for instructor review. No PP/KPSS adapter, commit or push was made.

## Implemented behavior

`adf()` and `dfgls()` are exported with the approved keyword-only signatures. Fixed order, including zero, is passed explicitly to public released ARCH 8.0.0. Automatic mode uses the unchanged independent selector with default MAIC or explicit MBIC/AIC/BIC, then passes the selected integer order to that same fixed-lag backend path.

The backend receives a validated copy of the **original series**. DF-GLS selection uses complete-series OLS detrending; its final GLS detrending happens exactly once inside ARCH. A constructor-spy test verifies the raw values, explicit lag argument and single public constructor call. ADF's final regression uses original levels and its requested deterministics, not selection residuals.

ADF supports `n`, `c`, `ct`, with `ctt` restricted to fixed lags. DF-GLS supports `c`/`ct`. Conflicting fixed/automatic options, invalid types/values and insufficient degrees of freedom fail explicitly. Automatic K uses the approved test-specific bound, warning and result note; the immutable selector result retains the original default/request metadata and the effective common selection sample. No candidate equation or penalty is changed.

`UnitRootResult` supplies frozen inference/sample/selection metadata, copied read-only critical-value/version maps, calibration provenance and public regression access. The regression is a documented mutable Statsmodels backend object, not deeply frozen. The critical values, p-value and sample size are read from public ARCH properties without substitutions. `reject()` uses the strict left-tail critical-value rule at 1%, 5% and 10%; exact equality does not reject. P-values never override that rule. `summary()` displays H0, terms, lag choice, distinct selection/final observation counts, approximate p-value, critical values and a nominal 5% decision, explicitly distinguishing failure to reject from establishing the null.

Only left-tail inference is implemented. Shared bandwidth/type fields remain None; they implement no PP/KPSS behavior.

## Test outcome and coverage

**797 passed, 0 failures, 0 errors, 0 skipped**, including all 337 Stage 3A cases. The original Stage 3A API/import-boundary guard was necessarily expanded to permit `adf`, `dfgls`, `UnitRootResult` and public ARCH imports **in the adapter module only**. No numerical selector test or fixture was removed or relaxed.

| Measure | Result |
|---|---:|
| Statements | 341 / 342 = **99.71%** |
| Branches | 107 / 108 = **99.07%** |
| Combined coverage | **99.56%** |
| New adapter module | **100%** statement and branch coverage |
| Shared result module | **100%** statement and branch coverage |

The single uncovered statement/branch is the unchanged defensive non-finite-criterion guard in the Stage 3A selector. No coverage exclusion was added.

The first integration run exposed two real error-handling gaps: Statsmodels could raise its rank-deficiency warning before the backend produced a public result, for a zero regressor and an extreme-scale design. The adapter now converts that **specific** `SingularMatrixWarning`, plus numerical RuntimeWarnings, into contextual, chained `EstimationError`. It does not suppress unrelated warnings or swallow unexpected programming errors. All final tests pass without changing numerical tolerances.

Tests exercise:

- Fixed p=0/4 against saved Stage 2 ARCH evidence and fresh public ARCH calls for every supported deterministic specification, including fixed ADF `ctt`.
- Every automatic method, supported trend and K=4/12 on the six reviewed vectors, with exact selected-order checks, final statistics, p-values, critical values and N−1−p sample counts.
- Explicit K=0, uncapped defaults, capped defaults, test-specific bounds, shortest accepted fixed/automatic samples and one below, and invalid explicit orders.
- Original-level DF-GLS handoff, immutable result/maps, public regression identity/labels, calibration/version provenance, approximate endpoint p-values and Series names.
- Strict comparisons immediately below, exactly at and immediately above each critical value, including intentionally conflicting p-values; unsupported alpha values are rejected rather than interpolated.
- Conflicting arguments, missing/complex/boolean/object/2-D inputs, constant/deterministic/singular designs, poor conditioning, unsupported backend versions, numeric-range failures, eager backend failures, malformed backend results and error chaining. Public backend symbols are temporarily replaced in fault-injection tests; no upstream source files are modified.

## Numerical comparisons

[stage3b_comparisons.csv](../validation/stage3b_comparisons.csv) contains **312** recorded wrapper/reference cases: 72 fixed and 240 automatic. [stage3b_numerics.json](../validation/stage3b_numerics.json) records the measured maxima, versions, coverage and source/preservation hashes.

| Quantity | Cases | Maximum absolute difference |
|---|---:|---:|
| Fixed statistic, p-value, each critical value versus public ARCH | 72 | **0** |
| Fixed statistic, p-value, each critical value versus saved Stage 2 ARCH | 72 | **0** |
| Automatic statistic, p-value, each critical value versus public ARCH at independently selected p | 240 | **0** |
| Automatic selected order versus independent selector | 240 | **0 lags** |
| DF-GLS statistic versus saved urca | 24 | **4.71e-14** |

“Zero” describes this run's identical numerical outputs, not a promise of bitwise portability across platforms. Saved-reference assertions retain the reviewed absolute **1e-8** tolerance; fresh-public-call assertions use **1e-12**. Integer choices and observation counts are exact. The 24 urca checks are additional to the 312 recorded wrapper cases. Their differing critical values are explicitly tested and **not** substituted into ARCH inference.

All modes return the final backend sample N−1−p. Automatic results separately retain selection sample N−1−K. In the executable example, N=240 and K=12 give 227 selection observations, while selected p=0 gives 239 final observations. Conventional AIC/BIC remain the approved two-step OLS/common-sample scores, not aliases for the modified criteria or a backend one-step autolag request.

## Preservation of Stage 3A and prior evidence

Byte hashes confirm `selection.py`, `_validation.py`, `errors.py`, existing numerical-fixture/input-validation tests, and every one of the eleven Stage 2 result files are unchanged. The `CandidateDiagnostics` and `LagSelectionResult` class definitions are AST-identical to the Stage 3A wheel. `results.py` adds the new result type; `__init__.py` adds approved exports and a shared version constant module. Historical Stage 3A/Stage 2 reports and result artifacts are retained, rather than rewritten to match the new API.

All Stage 3A cases continue to pass; the sole existing test adaptation is the explicitly staged API/import guard described above. The independent selector still imports no ARCH API, and package imports load the backend lazily only when an adapter is called. Production imports use `arch` and `arch.unitroot.ADF/DFGLS`; no private ARCH modules, selectors, calibration arrays, development snapshot or R runtime is required.

## Examples, packaging and reproduction

The updated [README](../README.md) has short complete examples. [examples/adf_dfgls.py](../examples/adf_dfgls.py) was executed successfully, displaying fixed ADF, automatic MAIC DF-GLS and automatic BIC ADF summaries, all three nominal rejection decisions, and public regression coefficients.

A regular wheel built successfully and its packaged exports were smoke-tested from the wheel itself. Build intermediates remain in the scratch workspace/ignored directories. No dependency update, package publication or upstream source modification was performed.

```sh
python -m pip install -e '.[test]'
python -m pytest -q -p no:cacheprovider \
  --cov=unitroot_tools --cov-branch --cov-report=term-missing \
  --cov-report=json:validation/stage3b_coverage.json \
  --junitxml=validation/stage3b_pytest.xml
python examples/adf_dfgls.py
```

The numerical CSV repeats the fixed/automatic grids covered by `tests/test_adapters.py`; it is a review artifact, not an additional inference implementation. Machine-readable test results are in [stage3b_pytest.xml](../validation/stage3b_pytest.xml) and coverage in [stage3b_coverage.json](../validation/stage3b_coverage.json).

## Tested environment and unresolved limits

Executed on Python **3.14.10**, macOS arm64, NumPy **2.4.4**, SciPy **1.17.1**, Statsmodels **0.15.0**, released ARCH **8.0.0**. Python 3.11+ remains a target; Python 3.11/3.12/3.13, lower dependency bounds, Linux and other BLAS environments have not been tested here and are not claimed compatible.

No unexplained numerical discrepancy remains for this scope. Backend estimation retains the original units: unsuitable extreme scales, rank reduction and effectively zero variance produce errors rather than silent rescaling, dropping regressors or misleading finite inference. Compatibility bounds and cross-platform behavior remain release gates.

The approved lag-selection contract and ARCH's approximations are not exact finite-sample or post-selection-adjusted inference. The prior full-text Perron–Qu provenance limitation remains documented; no paper-wide replication or size/power claim is added. PP/KPSS and their bandwidth/right-tail behavior are deferred to a separately authorized stage.
