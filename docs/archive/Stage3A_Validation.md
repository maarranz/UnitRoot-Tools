# Stage 3A implementation and validation

27 October 2026. Implemented the package foundation and independent lag selector only, in the existing repository. No commit or push was made.

## Scope completed

Before implementation, `Implementation_Specification.md` was amended with explicit conventional AIC/BIC formulas and common-sample normalization. The instructor's decisions are incorporated: namespace `unitroot_tools`, MAIC default, explicit MAIC/MBIC/AIC/BIC methods, and `ctt` reserved for a future fixed-lag ADF adapter.

The package includes one implemented function, `select_lags`, public `EstimationError` and `SpecificationWarning`, and frozen `LagSelectionResult` / `CandidateDiagnostics` records. Candidate tables are tuples of frozen scalar records; they contain all four scores and common-sample diagnostics. Full-series OLS detrending, scaled SVD estimation, sample/cap checks, warnings, and numerical-range failures are independent NumPy/SciPy calculations. No ARCH code or private/development API is imported or used. The production dependency remains pinned to **released ARCH 8.0.0** for future adapters.

Minimal `pyproject.toml`, source layout, test configuration, README and ignore rules are present. ADF, DF-GLS, PP, KPSS adapters, their result/decision objects, and bandwidth estimation are **not implemented or exported**, including no placeholder adapters.

## Test results

**337 passed; 0 failed; 0 skipped**, final run on Python **3.14.10**, macOS arm64. Warnings are errors except where a test explicitly checks the documented warning.

| Coverage measure | Result |
|---|---:|
| Executable statements | 156 / 157 = **99.36%** |
| Branches | 49 / 50 = **98.00%** |
| Combined statement/branch coverage | **99.03%** |

All source files except `selection.py` have 100% measured coverage. The uncovered statement/branch is the defensive check rejecting non-finite coefficients or criterion scores after a successful finite, full-rank SVD fit. It was not marked as excluded or artificially triggered to inflate coverage.

The first run had one incorrect test expectation: a detrended linear sequence became perfectly predictable before the later rank-deficient candidate. The test was corrected to use geometric lagged regressors with a final-response perturbation, isolating rank deficiency while keeping the earlier candidate estimable. The final suite has no remaining failures. A zero matched-level sum now produces an explicit rank-deficiency error rather than misleading rescaling advice.

Tests cover:

- All 36 reviewed PQ fixture groups (six series × three deterministic specifications × two K values), and every method's selected order and candidate diagnostics.
- Saved development DFGLS MAIC selected orders/minima, plus frozen bootUR two-step MAIC/MBIC choices and the known indexing discrepancy. These tests read reviewed CSVs; they require neither R nor development ARCH.
- A rational-arithmetic hand example distinguishing modified/conventional scores and selections; public Statsmodels OLS likelihood criteria on matched preprocessed regressions, with the declared common constants removed.
- The complete matched lagged-level sum; SSR/n versus residual degrees of freedom; common selection sample and zero-order behavior; default MAIC, label normalization, explicit maximums and automatic cap warnings.
- Exact/near ties, level/trend shifts where appropriate, sign and unit changes, frozen results/input copies, supported minimum samples, zero columns, singular/deterministic/nearly deterministic inputs, poor conditioning, missing/non-finite/type errors, and numeric diagnostic underflow/overflow.
- Numerical error chaining and production API/import boundaries. No deferred adapter appears in the package API.

## Numerical differences

The companion [stage3a_numerics.json](../../validation/stage3a_numerics.json) records versions, coverage, source/test hashes and measured maxima. Historical Stage 2 CSVs, scripts and manifest remain unchanged; their recorded hashes were rechecked.

Comparison against all **324** saved `PQ_full` candidate rows:

| Quantity | Maximum absolute difference |
|---|---:|
| Level coefficient | 3.45e-15 |
| SSR | 1.03e-12 |
| Variance | 6.22e-15 |
| Lagged-level sum of squares | 2.33e-10 |
| Tau correction | 4.84e-13 |
| MAIC | 1.65e-14 |
| MBIC | 3.38e-14 |
| Conventional AIC | 1.34e-15 |
| Conventional BIC | 1.34e-15 |

Every quantity passes the reviewed absolute **1e-8** tolerance, and selected orders agree exactly with the respective independently computed convention. Conventional reference scores are derived from the saved candidate variance and common n, and are additionally checked against public OLS likelihood scores; they are not claimed to have been stored as native AIC/BIC outputs in the earlier study.

`mic_boundary`, intercept, K=4 selects **MAIC 3, MBIC 2, AIC 0, BIC 0**. The bootUR-style truncated-sum MAIC choice remains 2, intentionally not copied into production. Across the 36 fixture groups, AIC differs from MAIC in 21, and BIC differs from MBIC in 21. These counts illustrate distinct penalties on these fixtures; they are not size/power or disagreement-frequency estimates.

## Installation and build

An editable installation in the external validation environment succeeded. A regular pure-Python wheel also built successfully; its contents, MIT license and exact ARCH 8.0.0 dependency were inspected, and importing directly from the wheel produced the expected selector result in a hand-example smoke test. No package was published. Wheel/build intermediates remain in the chat's scratch workspace and ignored build directories.

The first editable-install attempt was not executed because automatic approval review timed out; retrying after installing build tools with build isolation disabled succeeded. The build's deprecated TOML license-table format was replaced with an SPDX MIT license expression and `setuptools>=77`, then the wheel was rebuilt successfully.

Reproduce from the repository root in an appropriately isolated environment:

```sh
python -m pip install -e '.[test]'
python -m pytest --cov=unitroot_tools --cov-branch \
  --cov-report=term-missing \
  --cov-report=json:validation/stage3a_coverage.json \
  --junitxml=validation/stage3a_pytest.xml -q -p no:cacheprovider
python -m build --wheel --outdir /tmp/unitroot-stage3a-wheel
```

The last command additionally requires the `build` developer tool. Editable installation/build tools were installed only in the external scratch environment, not a repository virtual environment. Saved [coverage JSON](../../validation/stage3a_coverage.json) and [pytest XML](../../validation/stage3a_pytest.xml) provide machine-readable evidence.

## Support and methodological limits

Tested numerical versions: NumPy 2.4.4, SciPy 1.17.1, Statsmodels 0.15.0; installed ARCH 8.0.0. Test tools: pytest 9.1.1, pytest-cov 7.1.0. Python **3.11+ is a target**, not a tested compatibility claim: Python 3.11/3.12/3.13, lower dependency bounds, Linux and alternative BLAS environments were not exercised in Stage 3A.

Conventional AIC/BIC are explicitly the unmodified penalties on the approved **two-step OLS-preprocessed common-sample** family. They are not aliases for MAIC/MBIC and need not equal one-step deterministic ADF autolag implementations. No new methodological blocker was found for this defined contract. The earlier full-text Perron–Qu provenance limitation persists; no exact paper-replication, full Ng–Perron family or finite-sample size/power claim is made.

All four scores are reported in the declared original units. Inputs whose squared-unit diagnostics cannot be represented as positive finite float64 values fail explicitly with rescaling advice; this is a numerical range limit, not silent missing-value handling. Numerical singularity/deterministic thresholds and cap rules are documented in the governing specification.

Stage 3A stops here for instructor review. The next authorized stage would need separate instruction before adding ADF/DF-GLS or PP/KPSS adapters.
