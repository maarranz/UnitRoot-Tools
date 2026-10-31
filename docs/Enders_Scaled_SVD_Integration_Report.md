# Scaled-SVD Enders integration

31 October 2026 · Mac M1 Pro · qualified Python 3.13.18

**Original and Conservative Enders now use the validated scaled-SVD layer for their restricted and auxiliary regressions, while reusing the integrated ADF for unrestricted fits. All 1,597 source and installed-wheel tests pass.** The existing deterministic-selection procedure and public API are preserved.

## Changes and fitting architecture

Repository: `/Users/maarranz/GitHub/UnitRoot-Tools`, branch `main`, HEAD `37fa6efab21a5f332c062c263007411cdb46623e`.

Only one existing file changes: `src/unitroot_tools/_enders.py`.

- `_enders_regression` constructs a genuine named Statsmodels result through the existing `scaled_ols`, applies existing full-rank/variance/conditioning checks, and translates expected numerical failures to contextual `EstimationError` with exception chaining. Nonfinite arrays are rejected before Statsmodels model construction. Runtime and singular-matrix warnings are not suppressed.
- `_joint` fits the original retained columns through that wrapper. The empty restricted model (constant specification, lag zero) retains its exact residual equal to the original response and uses the helper's checked residual sum of squares. Restricted/unrestricted RSS use original response units. The formula `(RSS_R-RSS_U)/2 / (RSS_U/df_U)`, the two restrictions, negative-roundoff convention and existing invalid-statistic guard are retained. There is no replacement with a Wald test.
- `_auxiliary` uses the same stable wrapper on the original design with only `Level.L1` removed. Its coefficients, covariance, residual variance and deterministic-term t statistic come from that one fit.

Unrestricted estimation and final ADF tau/conditional statistics continue to use the existing integrated ADF. `_adapters.py`, `_regression.py`, calibration, lag selection, `_deterministic.py`, DFGLS, PP and KPSS are unchanged. There is no reduced-rank or raw-OLS fallback.

New files:

- `tests/test_enders_scaled_integration.py`: 155 tests covering independent numerical accuracy, exact restricted/auxiliary construction, all Conservative branches with both no-constant settings, exact boundaries at every supported alpha for both variants, genuine regression objects, errors, conditioning and no raw fallback.
- `validation/enders/scaled_references.json`: approximately 75 KiB of independent references for 18 saved-input cases: ordinary, positive 1e-12 and negative 1e-12 scales, c/ct, lags 0/1/4.
- `validation/enders/verify_scaled_references.py`: independently reproduces the references from saved hexadecimal binary64 inputs using 80/120-digit Decimal elimination, without importing production fitting code.
- This permanent report, `docs/Enders_Scaled_SVD_Integration_Report.md`; no earlier document was overwritten.

## Unchanged decision procedure

The entire `enders()` decision-procedure body is byte-for-byte unchanged; its hash is recorded in the external evidence. Result dataclasses, public signatures and finite phi tables are unchanged.

Both variants still automatically choose deterministic terms from n/c/ct. The user supplies the series, fixed lag, alpha and variant, without supplying a deterministic specification. Original retains its full sequence; Conservative retains its default constant-only policy stop and optional permission to continue to n. Sequential order, hypotheses, strict rejection boundaries, Student-t conventions, phi-table lookup, notes, sample alignment and chronological histories remain intact. No quadratic-trend selection or new result schema was introduced. The pre-existing optional lag-selection behavior was neither implemented anew nor changed; new numerical probes use fixed lags only.

## Independent accuracy and regression validation

The 18 cases include unrestricted, restricted and auxiliary fits, including the empty restricted model. Exact design/response arrays, names, samples, restrictions and residual degrees of freedom are checked separately from numerical values. Independent references agree at binary64 output precision between 80 and 120 digits; phi uses unrounded 120-digit RSS values. High-precision normal equations are reference arithmetic only, not the production solver.

Maximum errors over these cases:

| Quantity | Previous restricted/auxiliary fitting | Scaled-SVD integration |
|---|---:|---:|
| Absolute phi error | 2.22e-5 | 2.27e-13 |
| Relative restricted RSS error | 4.82e-8 | 3.33e-16 |
| Relative auxiliary RSS error | 3.30e-4 | 3.33e-16 |
| Absolute auxiliary t error | 0.1141 | 1.24e-14 |
| Auxiliary coefficient error / reference SE | 0.1142 | 2.08e-14 |
| Auxiliary covariance error / sqrt(reference Cii Cjj) | 7.07e-3 | 1.65e-15 |

Focused requirements remain 1e-12 absolute for relevant t/phi statistics and 1e-12 relative/standardized for RSS, coefficients and covariance. No existing test, fixture, expected value or tolerance was changed. Existing independent ordinary-scale QR/Wald-versus-RSS checks and the published GDP example pass unchanged.

Singular/zero-column designs, invalid residual df, exact deterministic fits, zero variance, nonfinite arrays and squared-unit overflow/underflow fail explicitly. A tested near-singular full-rank design keeps every column and emits the established conditioning warning; a design below the existing normalized rank threshold is rejected. Error-chain and no-raw-fit tests pass.

## Validation results

| Gate | Result |
|---|---:|
| Baseline Enders, ADF/helper and deterministic-selection tests | 689 passed |
| Focused Enders tests before additional all-alpha boundaries | 254 passed |
| Final complete package suite | 1,571 passed |
| Existing research suite | 26 passed |
| Full source total | **1,597 passed** |
| Separate installed-wheel total | **1,597 passed** |
| Existing independent Enders validator | Passed: 18 regression cases and GDP reproduction |
| Existing independent deterministic validator | Passed: 12 selections / 18 distinct fits |
| New independent reference reproduction | All 18 cases reproduced exactly |

All final tests pass without skips or warning suppression. Early new-test development exposed a pytest iterator-deprecation error, corrected by materializing its parameter list, and the nonfinite Statsmodels constructor escape, corrected by the input guard above. Original failed logs remain preserved. Neither issue required changing a scientific tolerance or decision rule.

A clean wheel was built externally with no dependency resolution/build isolation and installed into a separate target. Import and distribution paths were verified. Every runtime module matches staged source bytes. The complete tests, bounded fixtures and independent validators run against that installed target; installed bounded records equal source records exactly. Numerical threads are pinned to one; bytecode/cache writing is disabled. The qualified environment is not replaced by the wheel installation.

## Saved-series pre/post comparison

The bounded comparison evaluates 19 already saved vectors (including GDP and the nine portable Linux vectors), at lags 0/1/4, alpha 1%/5%/10%, Original, Conservative default and Conservative with no-constant permission: **513 complete Enders paths**.

Every selected specification, rejection decision, intermediate branch, critical value, sample/lag field, note and warning is unchanged. All n/c/ct outcomes occur. The largest changed statistic in these actually executed histories is **5.51e-14**. Direct restricted/auxiliary probes additionally expose the much larger previous extreme-scale errors in the accuracy table, even when a public procedure terminates before those branches. This distinction avoids claiming that all intermediate raw fits were already accurate merely because decisions agreed.

The comparison also records **171 ADF results and 114 AIC/BIC selections/scores: all are exactly unchanged**. ADF MacKinnon calibration and result contracts remain untouched. Existing controlled Original-path tests plus the new Conservative-path tests exercise all documented branch types; separate equality checks confirm strict nonrejection at tau, phi, auxiliary and conditional boundaries for both variants at all three alpha levels.

## Preservation, evidence and limits

The initial inventory covers 832 pre-existing repository files, including tracked, untracked and ignored work. Only `_enders.py` is authorized to change among them; existing ADF integration, research changes, documentation, notebooks, CSVs and saved numerical evidence are preserved. The qualified environment's 8,106 original files and dependencies are preserved. Final verification records hashes and timestamps, branch/HEAD, task-only changes and `git diff --check` outside the repository.

Large logs, complete regression arrays, baseline/staged snapshots, wheels and installed target remain at `/Users/maarranz/Documents/ChatGPT/UnitRoot-Tools/work/enders-scaled-svd-integration-2026-10-31-01/`. See `bounded-comparison.json`, `decision-procedure-preservation.json`, test XML/logs, wheel checks and preservation inventories there.

This qualifies numerical integration on the current Mac stack. It is not a size/power study, a fresh Linux/Colab qualification, a solution to the known cross-platform innovation discrepancy, or proof that intrinsic near-collinearity is harmless. The published GDP intercept-t discrepancy already documented by the independent validator remains explicitly recorded; no reference was rewritten to hide it.

No original pilot, full Monte Carlo run, RNG modification, automatic-lag implementation, AIC/BIC integration, commit, push or tag occurred. Stop here for planning-thread review before any further integration or platform qualification.
