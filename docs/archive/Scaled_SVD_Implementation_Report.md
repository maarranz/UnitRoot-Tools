# Column-scaled SVD implementation: staged result

30 October 2026 · qualified Python 3.13.18 · no commit/push/tag

**Stage A is implemented and validated. ADF integration stopped at the authorized boundaries; Enders/scoring integration did not proceed.** The attempted ADF change is preserved externally and its task-owned edits were restored to the exact starting bytes. Existing estimators therefore retain their starting behavior. A passing final suite must not be interpreted as completed estimator integration.

## Delivered files and fitting method

Three new repository files remain for review:

- `src/unitroot_tools/_regression.py`: `scaled_fit`, `residual_ssr`, `scaled_ols`, and private `ScaledFit` data.
- `tests/test_scaled_svd.py`: 21 focused numerical and genuine-result-object tests.
- `validation/scaled_svd/exact_references.json`: portable hexadecimal binary64 arrays from the original nine saved designs, original NPZ hashes, and independently computed 120-digit coefficient/covariance references. Original evidence was not changed or regenerated.

The helper validates shape, real/finite data, positive residual df, nonzero columns and response. Each column is divided by its maximum and then its normalized norm to avoid premature norm overflow. It performs a full SVD on the resulting unit-length columns, with **the existing EPS × max(shape) × largest-singular-value rank cutoff**, retaining every column or failing explicitly. No reduced-rank fallback, centering, detrending or sample modification occurs.

For X=ZD and y=cv, it solves g=V S^-1 U'v and maps beta=c D^-1 g. The normalized covariance is D^-1 V S^-2 V' D^-1; coefficient covariance multiplies this by **RSS/(n−k)**. No inverse of X'X is formed. Fitted values and residuals use the original design/response. Before accepting the ordinary residual dot product, `residual_ssr` cross-checks a scaled-norm calculation and rejects overflow, underflow, zero variance or inconsistent squared-unit arithmetic. Nonfinite coefficients, pseudoinverse, covariance, standard errors and t statistics fail explicitly.

`scaled_ols` constructs a **genuine Statsmodels OLSResults/RegressionResultsWrapper** from one coherent fit on the original named model. It does not edit a raw fitted result or return a look-alike object. Focused checks cover named parameters, unchanged arrays, fitted values/residuals, covariance, standard errors, t/p values, confidence intervals, prediction and ordinary restriction tests. This demonstrates a viable regression-result construction, not complete ADF calibration compatibility.

## Stage A accuracy and tests

All nine exact positive/negative/ordinary designs pass independent high-precision checks. Maximum observed errors:

| Quantity | Maximum error |
|---|---:|
| Relative RSS | 2.22e-16 |
| Coefficient error divided by reference standard error | 7.36e-15 |
| Covariance error divided by sqrt(reference Cii × Cjj) | 9.68e-16 |
| Relative standard error | 6.67e-16 |
| Absolute t statistic | 1.07e-14 |

The corresponding maximum raw-backend relative RSS error is 8.64784e-6. Tests also cover exact sample/df, singular and zero-column designs, zero response/residual variance, nonfinite inputs, numerical-range failure, and the existing normalized rank threshold on near-collinear designs. Ill-conditioned full-rank designs remain intrinsically sensitive; the helper does not change their statistical interpretation.

Final Stage A: **21 passed**. An initial new-test failure was a pandas DataFrame Boolean-reduction mistake in the test; converting its finite-value check to an array corrected that test bug without changing any numerical threshold. Initial and final logs are retained.

## Why Stage B stopped

The attempted ADF adapter constructed the released difference design, used the genuine stable OLS result, and calculated its matching p-value/critical values using ARCH's existing calibration routines. No raw ARCH statistic was combined with corrected RSS. The complete suite then reported **1,379 passed, 9 failed** (1,388 package tests including the 21 additions).

Four established `test_fixed_reference_evidence` cases require raw ARCH t-statistic parity within **1e-12**. Independently solving the same exact designs at 120 decimal digits shows a genuine reference-oracle conflict:

| Quadratic-trend case | p | Scaled − raw t | Scaled − high-precision t |
|---|---:|---:|---:|
| rw_negma | 0 | +1.995e-12 | −7.11e-15 |
| hetero_rw | 0 | +1.052e-12 | +2.22e-16 |
| rw_negma | 4 | +3.079e-12 | −1.78e-15 |
| trend_ar08 | 4 | +2.741e-12 | −1.78e-15 |

The stable outputs are accurate, but cannot satisfy the established raw-backend oracle in these cases. Designs/responses are exactly equal and all checked 1/5/10% rejection decisions remain unchanged. Neither expectations nor tolerances were changed, and no selectively inaccurate fallback was introduced to mimic raw ARCH.

The remaining five failures expose integration issues: three existing note-equality contracts, the explicit **public-ARCH-import boundary**, and the endpoint-calibration mock contract. In particular, the current architecture test permits only public `arch.unitroot` test classes; the attempted calls to `arch.unitroot.unitroot.mackinnonp/mackinnoncrit` violate that rule. These failures are preserved rather than hidden with dynamic imports or special test bypasses. Some metadata issues are straightforward implementation corrections, but the numerical oracle conflict and calibration boundary require review before continuing.

Per the stop instructions, only this task's `_adapters.py` edits were restored from its verified starting copy. `attempted-_adapters.py`, `attempted-adf-integration.patch`, the complete failing log/XML and `ctt-parity-conflict.json` preserve the attempt and both numerical results. No pre-existing research change was reverted.

## Final verification and incomplete stages

With ADF restored and the helper isolated:

- **1,414 passed**: original 1,367 package tests, 21 new helper tests and 26 research tests.
- A normal wheel was built from a new external source copy and installed into a separate target directory. Installed helper bytes match source; all **21 installed helper tests pass**. The qualified environment itself was not replaced by this test installation.
- `git diff --check` passes. Public exports/signatures and existing tests remain unchanged.

**Stage B is not complete. Stages C and D were not started.** No stable Enders restricted/auxiliary integration, AIC/BIC integration, new post-integration research replay, original-pilot rerun, full simulation or Colab qualification is claimed. The RNG and known Linux–Mac innovation hashes remain untouched. Source-level DFGLS, PP, KPSS and lag-selection behavior is unchanged.

## Preservation and next decision

The initial audit recorded branch `main`, HEAD `37fa6efab21a5f332c062c263007411cdb46623e`, all tracked/untracked state and existing documentation. All pre-existing repository file **bytes**, including every document, notebook, existing test, research file and validation artifact, match the starting inventory. Documentation/research/evidence modification times are unchanged; `_adapters.py` has a new modification time from applying/restoring this task's attempt, with original bytes verified. The tracked diff equals the original uncommitted research diff. Python 3.12's 9,670 inventoried files and original handoff/diagnostic hashes remain unchanged.

Recommendation: **NO-GO for further estimator integration or platform qualification until scientific/API review resolves these boundaries.** Review whether the accurate independent reference should replace raw-ARCH parity as the oracle for demonstrated biased cases, without loosening the numerical accuracy requirement. Separately authorize a coherent calibration bridge and its architecture tests. Do not bypass the public-backend test or silently relabel fitting/calibration provenance. Then resume Stage B and require a fully green source/installation gate before Enders/scoring and bounded deterministic research comparisons.

All new evidence is under `/Users/maarranz/Documents/ChatGPT/UnitRoot-Tools/work/scaled-svd-implementation-2026-10-30-01/`. `preservation-before.json`, `preservation-final.json`, test logs/XML, helper accuracy, exact conflict records, attempted patch and wheel-installation checks make the result reviewable. No existing report was updated, no dependency was added/upgraded, and no commit, push or tag was created.
