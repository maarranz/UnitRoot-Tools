# Scaled-SVD ADF integration and validation

31 October 2026 · Mac M1 Pro · Python 3.13.18

**ADF integration is complete and passes both source and installed-wheel validation.** No unexplained rejection-decision changes occurred. Work stops at ADF; no Enders restricted/auxiliary solver, AIC/BIC formula, RNG or Colab integration was undertaken.

## Implementation and files

Repository: `/Users/maarranz/GitHub/UnitRoot-Tools`, branch `main`, HEAD `37fa6efab21a5f332c062c263007411cdb46623e`.

Three existing files changed:

- `src/unitroot_tools/_adapters.py`: adds `_adf_regression` and `_adf_calibration`; the ADF branch now constructs the exact released difference regression and uses the existing `scaled_ols` helper. It never constructs an unused raw ARCH ADF. DFGLS fitting remains in its unchanged branch.
- `tests/test_adapters.py`: applies the explicitly authorized four-case accuracy-oracle policy and retargets affected ADF fit/calibration mocks. Ordinary ARCH parity and historical release-CSV assertions remain.
- `tests/test_unitroot_results.py`: retargets the endpoint mock to calibration, retaining exact 0/1 values, boundary labels and decision checks.

Five files added:

- `tests/test_adf_scaled_integration.py`: 28 integration checks covering exact construction/sample/df, public regression behavior, numerical accuracy, independent calibration, no raw ADF construction and numerical-range errors.
- `validation/scaled_svd/adf_ctt_references.json`: unchanged reference values and raw-ARCH diagnostic differences copied from preserved 120-digit evidence.
- `validation/scaled_svd/adf_inputs.json`: portable hexadecimal binary64 inputs from the saved nine numerical cases.
- `validation/scaled_svd/adf_calibration_references.json`: frozen ARCH-8 calibration evaluations, with source-evidence hash; no response-surface coefficients copied.
- `validation/scaled_svd/verify_adf_references.py`: independent 80/120-digit Decimal reproduction of the four references from the unchanged, hash-checked CSV.

All other source files, including `_regression.py`, `_enders.py`, `_deterministic.py`, selection, PP and KPSS, are unchanged. The complete reviewable task patch is `ADF_Integration.patch` beside this report.

## Coherent regression and calibration

The genuine Statsmodels `OLSResults`/`RegressionResultsWrapper` comes from one column-scaled full-rank SVD fit. Its original-unit coefficients, residuals, fitted values, RSS, covariance, standard errors and t statistics are mutually consistent. The ADF statistic is exactly its first coefficient t statistic. No reduced-rank fallback or raw-backend statistic is substituted.

Tests establish exact equality to ARCH's design/response arrays and named column order for n/c/ct/ctt at lags 0/4. Existing lag selection, requested metadata, effective samples and degrees of freedom remain intact. Genuine-result checks include prediction, covariance, conventional coefficient p-values, confidence intervals, t/F restrictions and likelihood. Rank/degeneracy, poor-conditioning warnings, zero variance, underflow/overflow and invalid-inference guards pass.

The private calibration seam imports only public `statsmodels.tsa.adfvalues.mackinnonp` and `mackinnoncrit`, with N=1, the requested trend and actual regression nobs. Conventional OLS coefficient p-values are not used as ADF p-values. **All 685 frozen p-value evaluations and 40 critical-value triplets match exactly**, including surface boundaries, adjacent floating-point values, tails, finite samples and infinity. Calibration equivalence is tested separately from regression accuracy. Critical values also remain exactly unchanged in all 13 numerical fixture results.

Existing `backend='arch'`, backend version and calibration label are retained as calibration-compatibility identifiers. The ADF docstring, regression description and p-value/critical-value provenance explicitly identify UnitRoot-Tools scaled-SVD fitting and public Statsmodels evaluation; they do not claim raw ARCH estimation. No result field or public signature changed. No unconditional fitting-method note was appended. Default notes-equality and architecture tests remain intact and pass.

Existing estimation-error categories, contextual messages and exception chaining are retained. Expected numerical failures are wrapped; unexpected programming errors still propagate. Only ADF mocks affected by the new boundary were retargeted; DFGLS constructor checks remain applicable.

## Precisely bounded reference changes

Only the statistic assertion for these four fixed-lag ADF cases now uses the preserved independent reference at **absolute tolerance 1e-12, rtol=0**:

| Series / trend / lag | Stable statistic | Stable minus independent reference | Stable minus raw ARCH |
|---|---:|---:|---:|
| rw_negma / ctt / 0 | −11.752256901055627 | −7.11e-15 | +1.995e-12 |
| hetero_rw / ctt / 0 | −1.2008196950031496 | +2.22e-16 | +1.052e-12 |
| rw_negma / ctt / 4 | −4.282434791175040 | −1.78e-15 | +3.079e-12 |
| trend_ar08 / ctt / 4 | −4.127040912512630 | −1.78e-15 | +2.741e-12 |

The reference values reproduced **exactly at both 80 and 120 digits**, using separate high-precision elimination rather than the production SVD helper. The reference procedure converts solved Decimal coefficients/covariance to binary64 before forming its reported t statistic, matching the preserved evidence. Its diagnostic normal equations are high-precision arithmetic, not a proposed binary64 production method.

The original fixture CSV, exact-binary64 helper reference file and helper tests are unchanged. All other raw-backend statistic assertions retain 1e-12; p-value parity, exact critical values and release-CSV checks retain their original requirements. There is no adaptive exception, relaxed tolerance or selectively inaccurate fallback.

## Validation results

| Gate | Result |
|---|---:|
| Before edits: helper + adapter + result contracts | 481 passed |
| Focused initial integration + architecture/selection | 538 passed |
| Full staged package suite, including 28 additions | 1,416 passed |
| Existing research suite | 26 passed |
| Final actual-repository package + research suites | **1,442 passed** |
| Clean wheel, installed package + research suites | **1,442 passed** |
| Existing independent Enders validator | Passed, 18 regression cases |
| Existing independent deterministic validator | Passed, 12 selections / 18 distinct fits |
| Independent ctt reference reproduction | All four exact |

No failures, skips, warning-policy changes or meaningful test bypasses were needed. The original 21 helper tests pass within the full suites. Every installed runtime module matches final repository source bytes; installed import/distribution paths were verified to use the separate wheel target, not repository or historical installed code. Build and install used no dependency resolution or build isolation. Numerical threads were pinned to one, bytecode and pytest cache writes disabled, and warnings retained under the existing error policy.

Nine original exact-binary64 designs through the **public integrated ADF** achieved:

| Accuracy measure against saved high-precision reference | Maximum |
|---|---:|
| Relative RSS error | 2.22e-16 |
| Coefficient error / reference standard error | 1.03e-14 |
| Covariance error / sqrt(reference Cii Cjj) | 1.14e-15 |
| Relative standard-error error | 6.66e-16 |
| Absolute t-statistic error | 7.11e-15 |

Across 13 saved numerical fixtures, statistics/p-values changed as expected from stable fitting. Maximum absolute statistic change was **0.00181364** on the extreme-scale ct cases; maximum absolute p-value change was **1.79e-13**. All 1%, 5%, 10% decisions, samples, lags, notes and warnings were unchanged. Exact changes, full regression arrays and raw/stable comparisons are in `baseline-bounded.json`, `staged-bounded.json` and `bounded-comparison.json`.

Nine saved Linux portable vectors were evaluated without generating new simulations. All selection/Enders paths, rejection decisions, input hashes and warnings matched; maximum floating difference versus the previous Mac path was 1.14e-13 and versus Linux evidence 1.99e-13, within the unchanged 1e-8 portable requirement. The 27 DFGLS/PP/KPSS comparison records were exactly unchanged. Installed-wheel fixture records are exactly equal to source records. Enders unrestricted quantities inherit the intended tiny ADF arithmetic changes; its rules, restricted/auxiliary implementations and checked behavior are unchanged.

## Preservation, limits and next step

The initial inventory covers **826 existing repository files, including 39 documentation files**. Only the three authorized existing files changed; all other 823 files retain bytes and modification times. All pre-existing research changes, notebooks, ignored evidence, documentation, original fixture CSV and exact references are preserved. The qualified environment's **8,106 files retain bytes and modification times, with no added files**. Its dependencies and historical installed package were not replaced. Tests used explicit source paths or the new external installed-wheel target. Python 3.12 was not used or modified.

`git diff --check` passes; branch and HEAD are unchanged. Inventories, logs/XML, commands, copied baseline, installed-wheel checks and checksum manifest remain in this external work directory. No commit, push, tag, original pilot, full simulation or Colab qualification occurred.

Qualification is for this Mac and its established Python 3.13.18 numerical stack. Frozen calibration/portable evidence supports the tested domain; fresh Linux/Colab integration runs and the known innovation-portability issue are not resolved by this assignment. ARCH 7.2 compatibility was not revisited. Near-singular identification and nonrepresentable squared units continue to fail explicitly.

Recommended next step: review this ADF patch and evidence in the planning thread. Separate authorization is required before integrating stable restricted/auxiliary Enders fits, addressing AIC/BIC integration or starting platform qualification/simulations.
