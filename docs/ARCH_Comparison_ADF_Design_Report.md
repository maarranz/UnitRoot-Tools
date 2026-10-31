# ARCH comparison and ADF integration design

30 October 2026 · Mac M1 Pro · Investigation only

Recommendation: use the public `statsmodels.tsa.adfvalues.mackinnonp` and `mackinnoncrit` functions to calibrate the stable ADF statistic. Retain independent high-precision accuracy requirements and ordinary-scale ARCH parity. No production code or tests were changed.

## A. Versions, environment and comparison limits

The repository remains on `main`, HEAD `37fa6efab21a5f332c062c263007411cdb46623e`. The initial audit included existing tracked modifications, untracked/ignored files, afternoon reports, the scaled-SVD helper/tests, exact reference evidence and the preserved unsuccessful ADF integration. All 825 repository files, including 38 documentation files, retained their bytes and modification times; Git status and tracked diff are identical. See `preservation-before.json` and `preservation-after-containment.json`.

There is no PyPI release ARCH 7.9 or 7.9.0: both release endpoints returned 404. The latest actual 7.x release is **7.2.0**, selected for comparison with **8.0.0**. Original release metadata, source archives and verified SHA256 values are saved here. [PyPI release history](https://pypi.org/project/arch/), [ARCH 7.2.0](https://pypi.org/project/arch/7.2.0/).

Two separate environments were installed with identical dependencies except ARCH, using Python 3.13.18, NumPy 2.4.4, SciPy 1.17.1, Statsmodels 0.15.0 and pandas 3.1.0 on macOS arm64/Accelerate. Both pass dependency metadata checks. Full freezes, build configurations, commands and installation logs are retained. UnitRoot-Tools was not installed into these comparison environments, avoiding its ARCH-8 constraint being misrepresented as satisfied by ARCH 7.

**ARCH 7.2.0 fails at import**: its use of pandas' `deprecate_kwarg("type", "plot_type")` raises `TypeError` because pandas 3.1 requires an additional warning-class argument. ARCH 8 uses its vendored decorator. See `probe-7.2.0.log` and `failed-7.2.0-environment.json`. No downgrade or compatibility patch was applied. Therefore no genuine runtime ARCH-7-versus-8 numerical comparison is claimed.

Source inspection provides narrower evidence:

| Component | 7.2.0 versus 8.0.0 |
|---|---|
| `_estimate_df_regression` | Identical source: same differencing, lag construction, trends and Statsmodels OLS fit |
| `ADF._compute_statistic` | Typing/string changes; no relevant numerical change |
| `mackinnonp` | Typing/error-text changes; same numerical logic |
| `mackinnoncrit` | Same response surfaces; 8 adds `.astype(float)` to finite-sample output |
| Dickey–Fuller coefficient module | Same expressions excluding the typing import |

Thus inspected source does not identify an ADF algorithm/calibration change explaining the observed discrepancy. This is not proof of cross-version runtime identity.

The isolated ARCH-8 environment and qualified ARCH-8 environment produced **exactly identical complete numerical records for all 13 saved fixtures**: designs, responses, samples/lags, coefficients, residuals, RSS, covariance, standard errors, statistics, p-values and critical values. All recorded warnings were empty. Numerical maximum differences were zero. Full arrays and hashes are in `results-8.0.0/adf-fixtures.json`; the comparison is in `isolated8-qualified8-numerical-comparison.json`.

## B. Calibration architecture

**A — public ARCH:** Neither inspected release exports a public standalone MacKinnon evaluator from `arch.unitroot`. The public ADF constructor does not accept an external statistic; its read-only inference properties are not a supported calibration injection point. Mutating `_stat` or depending on private computation hooks would not solve the architecture requirement. [Public ADF interface](https://arch.readthedocs.io/en/latest/unitroot/generated/arch.unitroot.ADF.html).

**B — existing public implementation: recommended.** Statsmodels explicitly exports `mackinnonp` and `mackinnoncrit` in `statsmodels.tsa.adfvalues.__all__`. They implement the matching ADF conventions: MacKinnon 1994 asymptotic p-values and 2010 critical values, with the historical no-constant critical-value convention. [Public source and docstrings](https://github.com/statsmodels/statsmodels/blob/main/statsmodels/tsa/adfvalues.py).

Use a small private UnitRoot-Tools calibration seam, conceptually:

- `mackinnonp(stable_t, regression=trend, N=1)`;
- `mackinnoncrit(N=1, regression=trend, nobs=effective_nobs)`.

Pass the actual effective ADF regression sample, normally input length minus one minus selected lags. Preserve `n/c/ct/ctt` exactly. Do not special-case `n` by substituting infinite sample size. Do not use conventional OLS coefficient p-values as ADF p-values.

Verification found **685 p-value checks and 40 critical-value triplets exactly equal** to ARCH 8, including all deterministic cases, response-surface transitions and adjacent floating-point values, tails, finite samples and infinity. Matching N=1 coefficient surfaces/cutoffs were inspected without implementing or copying coefficients. Applying both evaluators to the improved statistics for all 13 fixtures also yielded exact equality; all 1%, 5% and 10% rejection decisions remained unchanged. See `results-8.0.0/calibration-grid.json`, `calibration-threshold-equivalence.json` and `public-calibration-on-improved-statistics.json`. These results qualify the installed versions and tested domain, not every possible future dependency version.

**C — internal ARCH adapter:** unnecessary while B is supported and equivalent. If B ceases to qualify, an explicitly reviewed, centralized, pinned two-function dependency would require an architecture exception and version tests; no bypass is proposed now.

**D — independent calibration:** unnecessary coefficient provenance, licensing, convention and numerical maintenance burden. No coefficients were implemented or vendored.

All reported ADF regression quantities must come from the same stable fit; its statistic feeds the public calibration seam. Critical values use that fit's sample/specification. Keep provenance explicit about UnitRoot-Tools fitting and Statsmodels evaluation of ARCH-equivalent surfaces. Review existing backend-label semantics rather than silently claiming the regression was fitted by ARCH. Use regression description/provenance fields for method information, preserving default notes.

## C. Reference-oracle policy

The preserved 120-digit investigation establishes that the four quadratic-trend conflicts are raw-backend accuracy limits, not grounds for relaxing the accuracy threshold:

| Fixture | Current scaled minus raw statistic | Previously saved scaled minus 120-digit reference |
|---|---:|---:|
| `rw_negma`, ctt, lag 0 | 1.997e-12 | −7.1e-15 |
| `hetero_rw`, ctt, lag 0 | 1.050e-12 | +2.2e-16 |
| `rw_negma`, ctt, lag 4 | 3.078e-12 | −1.8e-15 |
| `trend_ar08`, ctt, lag 4 | 2.741e-12 | −1.8e-15 |

Small last-bit differences from the earlier report are retained in the separate current evidence rather than rewriting it. Extreme-scale positive/negative ct fixtures additionally show a raw-versus-stable statistic difference of about 0.001814; calibration equality does not repair a raw regression error.

Proposed assertion policy:

| Category | Preserve or revise |
|---|---|
| Econometric definition | Preserve exact design, deterministic columns, lag/sample, degrees-of-freedom and rejection-tail checks |
| Independent accuracy | Preserve existing scaled-SVD exact-binary64/120-digit evidence and tight statistic requirements; add the four saved ctt references as explicit fixtures |
| Backend parity | Keep ordinary-scale checks. For precisely the four demonstrated ctt cases, replace raw ARCH as accuracy oracle with independent references; retain raw differences as diagnostics |
| Calibration | Test both evaluators at the same supplied statistic, separately from regression accuracy; preserve exact critical-value checks |
| API/metadata | Preserve notes, lag requests, warnings, errors, endpoint labels and result contracts |

Locations requiring review are `/Users/maarranz/GitHub/UnitRoot-Tools/tests/test_adapters.py:19` (`check_inference`, particularly raw-statistic assertion at line 20) and `:31` (`test_fixed_reference_evidence`). Limit oracle substitution to ADF, ctt, the four series/lag combinations above. No automatic exception based on observed error, blanket tolerance increase or changed historical CSV is proposed. Keep the existing release-reference CSV and its 1e-8 checks. Expected reference values and current tests remain unchanged in this assignment.

## D. Five nonnumerical integration failures

| Failures | Cause and smallest coherent correction |
|---|---|
| Three notes-equality cases | The attempted adapter appended an unconditional fitting-method note, breaking `test_default_caps_preserve_request_metadata` (`tests/test_adapters.py:85`, equality at :95). Preserve equality and genuine warning notes; place method information in description/provenance. Implementation compatibility mistake. |
| One import-boundary failure | Importing `arch.unitroot.unitroot` violated the intended AST boundary at `tests/test_selection.py:141`. Public Statsmodels calibration avoids this dependency; retain the architecture test unchanged. |
| One endpoint-mock failure | A concrete `isinstance(..., ReleasedADF)` gate rejected the protocol mock in `tests/test_unitroot_results.py:81`. The stable path also needs endpoint calibration of its own statistic. Retarget this mock to the calibration seam, retaining exact 0/1 endpoint behavior and labels; do not special-case the mock or weaken the endpoint contract. |

These five are implementation/API-contract issues, not econometric conflicts. Related constructor/error/warning/regression-guard tests at `tests/test_adapters.py:176–246` will need ADF mocks directed at the stable fit/calibration seams while preserving validation and exception chains. DFGLS backend mocks remain applicable. Do not instantiate an unused raw ADF merely to satisfy old mocks.

## E. Bounded next implementation specification

After planning review, modify only the ADF path in `/Users/maarranz/GitHub/UnitRoot-Tools/src/unitroot_tools/_adapters.py` (`_backend`:23, `_check_regression`:65, `_estimate`:99, `adf`:176); preserve DFGLS at :190. Assemble the same design/sample, use the delivered scaled-SVD helper, derive coefficients/covariance/t from one fit, and calibrate that t through the public Statsmodels seam. Preserve selection, notes, warnings, metadata and failure handling.

Approve the explicit four-case reference-oracle policy and provenance semantics alongside implementation. Retarget only mocks affected by the new dependency seams. Validate exact reference fixtures, calibration boundaries, ordinary ARCH parity, all public contracts, full source/wheel suites and bounded saved fixtures. Investigate any new decision change before proceeding. Complete ADF before Enders and AIC/BIC; no RNG, pilot, Monte Carlo or Colab work is included.

### Evidence and preservation

All new material is in this directory. `installation-commands.json`, full freezes, source comparisons and probe scripts make the investigation reproducible; numerical probes pin BLAS/OpenMP threads to one. The qualified environment's 8,106 original files are unchanged and its dependency set was not altered. One inspection inadvertently generated 83 new bytecode cache files there; only those newly generated files were archived under `generated-bytecode-from-probe/`, restoring the original file inventory. The preliminary and corrected preservation records are both retained (`bytecode-containment.json`). No existing evidence was overwritten. No production/test edits, commits, pushes, tags or simulations were performed. The investigation stops here.
