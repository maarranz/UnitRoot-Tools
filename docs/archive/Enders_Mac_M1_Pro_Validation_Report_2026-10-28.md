# Enders Mac M1 Pro Validation Report

28 October 2026 · UnitRoot-Tools 0.1.0rc1 · validation only

## Assessment

**PASS for the tested Mac configuration.** The current source checkout, a newly built wheel, and a fresh installation from the corresponding source distribution each passed all **1,287 package tests**, **123 separately rerun Enders tests**, and the existing independent numerical-validation script. There were no test failures, errors, skips or reported test warnings. All three numerical reports are identical, including the GDP path and numerical discrepancies.

The implementation is ready for the next separately scoped development stage on this configuration. No econometric correction is indicated by these checks. The planned experimental deterministic AIC/BIC feature can now be specified and discussed separately; it was not implemented. Existing repository distribution files are stale and exclude Enders, so they must not be used as distributions of this checkpoint. Fresh validation artifacts worked correctly but were not copied into the repository or published.

This is numerical/software validation, not evidence that either complete sequential procedure controls nominal statistical size. Conservative Enders remains a proposed modification. Monte Carlo and bootstrap work were not performed.

## Repository and preservation

- Repository: `/Users/maarranz/GitHub/UnitRoot-Tools`.
- Branch: `main`.
- HEAD, independently checked before and after: `8c4f1a17cf12207dba3c051be6f1de7943df895c`.
- History: Enders implementation commit `329468c`; successor `8c4f1a1` corrects the GDP CSV hash metadata.
- Initial and final status: `?? docs/Theoretical_Audit_Sequential_Strategies.md`.
- SHA-256 comparison of **780 existing files outside `.git`** found no additions, removals or content changes in the repository. This includes source, tests, docs, untracked audit, historical evidence, original distributions, and existing cache/coverage files.
- No pull, dependency upgrade, environment recreation, metadata/version edit in the repository, commit, tag or push occurred. No existing report was overwritten.

Testing instructions inspected: `pyproject.toml`, `docs/Enders.md`, `validation/README.md`, `validation/stage4/Reproduce.md`, `validation/validate_enders.py`, `validation/inspect_distributions.py`, and the installed-validation runner. The test configuration treats warnings as errors and requires strict markers.

## Machine and Python environment

| Item | Verified value |
|---|---|
| Hardware | MacBook Pro, MacBookPro18,3, Apple M1 Pro, 8 cores, 16 GB memory |
| OS | macOS 26.7.1, build 25G241 |
| Process architecture | native arm64 / aarch64 |
| Python used | CPython 3.12.16; Clang 23.1.3 |
| UnitRoot-Tools version | 0.1.0rc1 |
| ARCH | 8.0.0 |
| NumPy / SciPy | 2.4.4 / 1.17.1 |
| Statsmodels / Pandas | 0.15.0 / 3.1.0 |
| pytest / pytest-cov | 9.1.1 / 7.1.0 |
| setuptools / pip | 84.2.1 / 26.2.2 |
| NumPy BLAS / LAPACK | Apple Accelerate |

No environment was activated in the shell (`CONDA_PREFIX` and `VIRTUAL_ENV` absent). The shell-default Python 3.15.2 had none of the project/test dependencies. The existing conda `unitroot_release` environment also lacked runtime/test dependencies. Neither was changed. The previously used, dependency-complete project validation interpreter was selected explicitly:

`/Users/maarranz/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/stage4b/env312-current/bin/python`

Its previously installed UnitRoot-Tools copy predates Enders. Source validation therefore explicitly imported the current checkout via `PYTHONPATH=.../src` and the repository's source-test configuration. For artifact tests, newly installed target directories took precedence; import paths and byte identity of every runtime Python module were verified. No result is attributed to the old installed copy.

## Test execution

| Code under test | Complete suite | Focused Enders suite | Failures / errors / skips, each run | Independent script |
|---|---:|---:|---|---|
| source | 1287 pass (6.460s) | 123 pass (1.262s) | 0 / 0 / 0 | pass; 18 cases |
| wheel | 1287 pass (6.307s) | 123 pass (1.341s) | 0 / 0 / 0 | pass; 18 cases |
| sdist | 1287 pass (6.306s) | 123 pass (1.339s) | 0 / 0 / 0 | pass; 18 cases |

These are repeated executions of the same tests, not thousands of independent statistical experiments. The separate focused run deliberately repeats the 123 Enders cases already included in the full suite.

The focused suite includes 50 controlled Original-path combinations with immediate-stopping assertions; the published-GDP Conservative distinction and permitted-variant equivalence; Original's irrelevant-flag handling; all four lag-selection methods and common-sample checks; 30 sample-grid/significance-level combinations testing both phi tables; invalid options/series; tail/result metadata; default Conservative policy termination; too-few-effective-observations rejection; and equality-at-critical-value nonrejection for tau, phi, auxiliary and conditional tests.

An additional external, fixture-based consistency check performed **84 comparisons** across four existing validation/test fixtures, three alpha levels, fixed lags 0/1/4 and automatic MAIC/MBIC/AIC/BIC with maximum lag 8. Original and permitted Conservative matched in terminal results, scalar histories, coefficients, standard errors and residual arrays. Default Conservative rejection always implied Original rejection; regression observation counts matched. These are deterministic regression/consistency checks using existing fixed-seed fixtures, not a size/power simulation. Repository tests were not edited.

Coverage from each full run was identical: package statements **650/657 (98.93%)**, package branches **174/180 (96.67%)**; Enders statements **144/149 (96.64%)**, Enders branches **29/34 (85.29%)**. Unexecuted Enders lines are defensive checks for unexpected ARCH column names, invalid joint statistic, sample mismatch, and the unreachable final guard (lines 102, 110, 119, 181, 215). Passing branch tests does not mean 100% defensive-path coverage. No tests or tolerances were changed to improve coverage.

## Independent numerical reconstruction

The existing `validation/validate_enders.py` was inspected before execution. Required input is the checked-in `validation/enders/rgdp.csv` and `source.json`, beside the script; the script verifies the fixture's SHA-256. Its required `--report` destination was outside the repository. Python assertions were enabled.

It reconstructs 18 cases: three controlled vectors × ct/c × lags 0/1/4. Independent QR solutions supply coefficients, residuals, covariance matrices and t statistics. Phi is checked both as an RSS restriction statistic and through an independently computed two-restriction Wald quadratic form. Restricted-model RSS and auxiliary regressions are independently checked. No ordinary F quantiles replace the nonstandard phi tables.

All maximum absolute discrepancies below met the unchanged **1e-8** numerical threshold:

| Quantity | Maximum absolute discrepancy |
|---|---:|
| coefficient | 4.69624339416e-14 |
| standard_error | 5.56499291093e-15 |
| ssr | 5.55111512313e-17 |
| t_statistic | 7.66053886991e-13 |
| phi | 3.11217718263e-12 |
| aux_coefficient | 2.40529818285e-13 |
| aux_standard_error | 3.74977826567e-14 |
| aux_ssr | 2.84217094304e-14 |
| aux_t | 1.19548815292e-12 |

The largest was **3.112177182629239e-12** for phi, approximately 3,213 times below the threshold. Tests additionally use their existing tighter quantity-specific relative/absolute comparisons, including rtol 1e-10 and atol 1e-12 for regression quantities. All passed unchanged.

Tau calibration remained ARCH 8.0.0. Focused tests verified phi1/phi3 table entries, row transitions at 25/50/100/250/500 observations, the largest-grid-row convention, the cap beyond 500, refusal below 25, alpha 1%/5%/10%, and strict rejection inequalities. Auxiliary two-sided and conditional left-sided finite-df t critical values were checked by the existing tests. This confirms implementation of the documented conventions, not exact finite-sample or post-selection validity.

## Published GDP reproduction

Input: 264 quarterly RGDP levels, transformed to natural logs; augmentation p=1; **262 regression observations**; alpha=0.05. The complete Original path reproduced:

`tau3 → phi3 → tau2 → phi1 → auxiliary_intercept → conditional_t`.

| Test | Mac statistic | Published rounded statistic | Absolute difference | Mac critical value | Step decision |
|---|---:|---:|---:|---:|---|
| tau3 | -1.489116883600 | -1.49 | 0.000883116400 | -3.427381410512 | nonreject |
| phi3 | 2.965689732029 | 2.97 | 0.004310267971 | 6.340000000000 | nonreject |
| tau2 | -2.048529930013 | -2.05 | 0.001470069987 | -2.872635586277 | nonreject |
| phi1 | 25.485041627542 | 25.49 | 0.004958372458 | 4.630000000000 | reject |
| auxiliary_intercept | 6.797466523465 | 6.80 | 0.002533476535 | 1.969130003360 | reject |
| conditional_t | -2.048529930013 | -2.05 | 0.001470069987 | -1.650758133653 | reject |

All six path statistics are within **0.005** of the manual's two-decimal values. Both joint entries use the T=250 row. Phi/auxiliary rejection determines the route; it is not itself a unit-root rejection.

- **Original:** selects c and rejects the unit-root null through `conditional_t`.
- **Default Conservative:** selects c and stops at the `tau2` nonrejection; no phi1 or conditional-intercept branch is run.
- **Conservative with `allow_no_constant=True`:** reproduces Original's path and rejection.

All nine printed coefficient comparisons passed the existing **0.00005** tolerance, including the documented intercept transformation to the publication's original time index. Eight of nine additional printed t-statistic comparisons were within 0.005. The previously documented exception remains: A2.3 intercept t = **2.646058967065688**, printed **2.64**, difference **0.006058967065688**. The script explicitly records this exception and verifies its identity; it does not relax the general rounding threshold. Thus `all_passed=true` coexists correctly with `published_t_statistics_all_within_rounding=false`. This is not a newly discovered Mac failure and does not alter the decision path.

The source's inconsistent mention of 244 versus 262 usable observations remains documented. The reproduced sample follows the 264-level fixture and p=1, yielding 262 observations. ARCH tau response-surface critical values need not equal rounded textbook critical values.

## Distribution validation

The repository's existing wheel and sdist both lack `_enders.py`; they are pre-Enders artifacts despite sharing the version 0.1.0rc1. They were inspected without modification and were not used as current-checkpoint validation evidence.

A fresh committed-source copy was obtained with `git archive HEAD`. It intentionally excludes the untracked theoretical-audit document. Using existing setuptools 84.2.1, an sdist was built in that copy and a wheel was built from the extracted sdist. All generated metadata/build output stayed in this chat's validation directory. The existing distribution-inspection script passed against the clean committed-source copy: runtime-module bytes, documentation/examples/license, metadata and required sdist tests/fixtures matched.

- `unitroot_tools-0.1.0rc1-py3-none-any.whl`: 107,138 bytes; SHA-256 `8b13620bf525721709e43e6f445608127b5825d25b89e1a7bdbf67f4e14c4a9f`.
- `unitroot_tools-0.1.0rc1.tar.gz`: 257,972 bytes; SHA-256 `115ead7870b122ac8301e8ece05e6efd6615b1456f7282f0a4bba9b057bc45bc`.

Each artifact was installed with pip `--target` into a different temporary directory, using `--no-deps --no-build-isolation --no-index --no-cache-dir`. Tests ran from copied test/fixture directories outside the checkout, with an external pytest configuration containing **no source pythonpath**. Both package import locations and all runtime module bytes were verified before testing. Full suites, focused suites, independent scripts and pip dependency checks all passed. Source/wheel/sdist numerical JSON reports are exactly identical in this environment.

Isolation limitation: the two installed package targets share the unchanged existing interpreter and numerical dependencies; they are not freshly provisioned virtual environments or additional Python-version/platform tests. This validates artifact installation and import separation while honoring the restriction against environment recreation and dependency changes. No twine upload, release tag or publication was attempted.

## Linux comparison and limits

| Evidence dimension | Conclusion |
|---|---|
| Software test execution | New Mac result: 1,287 pass, matching the earlier **reported** Linux count of 1,287. Focused Mac result: 123 pass. |
| Independent numerical agreement | Mac's QR/RSS/Wald comparisons pass at 1e-8; maximum 3.1122e-12. Earlier conversation reported Linux discrepancies below 2.3e-12. These are internal validation discrepancies, **not direct Mac-minus-Linux differences**. |
| Published GDP example | Mac reproduces the same reported Original/Conservative distinction and the previously documented printed-intercept exception. |
| Direct cross-platform numeric equality | Not quantified: original Linux JSON/log/archive was not available in the inspected repository or searched local document/download locations. The historical attachment was referenced in earlier conversation but its contents were unavailable. |
| Whole-procedure nominal size | Not assessed or established by these checks; belongs to theoretical and later simulation research. |

The Mac result is independently executed evidence. It is not inferred from Linux history. The difference between the reported Linux error bound and the Mac internal maximum is not a failed tolerance and cannot identify a platform difference without matched raw results. Retain that direct comparison as an evidence limitation rather than manufacturing a cross-platform discrepancy figure.

## Warnings and failures

- **Tests:** zero failures, errors, skips and reported test warnings in all full/focused runs. Warnings-as-errors remained enabled.
- **Independent scripts / consistency checks:** all passed; known publication discrepancy preserved as described above.
- **Build:** harmless manifest warning that no `validation/*.xml` files exist in the committed-source copy; byte-compilation warnings reflect intentionally setting `PYTHONDONTWRITEBYTECODE=1`. Artifact inspection and installations passed.
- **Environment diagnostics:** NumPy `show_config()` recommends optional `pyyaml` for better display; its fallback still reports Accelerate configuration. No dependency was added. An initial pip-check command disabled an unwritable default cache; the recorded no-cache rerun passed without that warning. A CPU sysctl query was permission-restricted and system_profiler warned about `hw.cpufamily`, but its hardware overview identified Apple M1 Pro. No permissions or machine settings were changed for these diagnostics.
- **Saved artifact staleness:** the old repository distributions omit Enders. Recommendation: replace them only in a separately authorized packaging/release step. No source correction is needed for this issue; fresh artifacts from current committed sources passed.

## Reproduction commands and evidence

All paths below are literal paths used in this validation. `PYTHONPATH` intentionally selects source for source checks and the artifact target for installed checks. No shell activation is required.

```sh
TASK_ROOT=/Users/maarranz/Documents/Codex/2026-10-28/continue-the-unfinished-theoretical-audit-of
REPO=/Users/maarranz/GitHub/UnitRoot-Tools
PYTHON=/Users/maarranz/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/stage4b/env312-current/bin/python
OUT="$TASK_ROOT/work/enders-mac-validation"

cd "$REPO"
git branch --show-current
git rev-parse HEAD
git status --porcelain
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export PYTHONPATH="$REPO/src"
export COVERAGE_FILE="$OUT/source.coverage"
"$PYTHON" -m pytest -p no:cacheprovider --junitxml="$OUT/source-full.xml" \
  --cov=unitroot_tools --cov-branch --cov-report=json:"$OUT/source-coverage.json"
"$PYTHON" -m pytest tests/test_enders.py -p no:cacheprovider \
  --junitxml="$OUT/source-focused.xml"
"$PYTHON" validation/validate_enders.py --report "$OUT/source-numerical.json"
"$PYTHON" "$OUT/equivalence_checks.py"
"$PYTHON" -m pip --no-cache-dir check
```

Actual stdout/stderr were redirected to the corresponding `.log` files. The external consistency checker and installed-validation runner are preserved with the evidence; they do not modify repository tests. To reproduce without overwriting evidence, use a new output directory and update the checker/runner location accordingly.

For packaging, a fresh directory was populated from `git archive HEAD`, without untracked files. In `build-source/`, the backend invocation was `setuptools.build_meta.build_sdist(str(out))`, with `out=../fresh-dist`. The resulting sdist was extracted into `sdist-source/`; inside its package directory, `setuptools.build_meta.build_wheel(str(out))` used `out=../../fresh-dist`. Both invocations ran under the same Python with byte-compilation disabled. Equivalent exact backend expressions are:

```sh
cd "$OUT/build-source"
"$PYTHON" -c 'from pathlib import Path; from setuptools.build_meta import build_sdist; build_sdist(str(Path("../fresh-dist").resolve()))'
cd "$OUT/sdist-source/unitroot_tools-0.1.0rc1"
"$PYTHON" -c 'from pathlib import Path; from setuptools.build_meta import build_wheel; build_wheel(str(Path("../../fresh-dist").resolve()))'
cd "$OUT/build-source"
"$PYTHON" validation/inspect_distributions.py --dist ../fresh-dist \
  --report ../distribution-inspection.json
cd "$OUT"
"$PYTHON" -m pip install --target wheel/site-packages --no-deps \
  --no-build-isolation --no-index --no-cache-dir \
  fresh-dist/unitroot_tools-0.1.0rc1-py3-none-any.whl
"$PYTHON" -m pip install --target sdist/site-packages --no-deps \
  --no-build-isolation --no-index --no-cache-dir \
  fresh-dist/unitroot_tools-0.1.0rc1.tar.gz
"$PYTHON" validate_installed.py wheel
"$PYTHON" validate_installed.py sdist
```

The installed runner sets each target's import path, disables bytecode/cache writes, copies the existing committed tests/fixtures, writes a separate warnings-as-errors pytest configuration, and executes import provenance, full suite, focused suite, numerical validation and pip check. Exact argument vectors, working directories, return codes and environment overrides are in each target's `commands.json`. The script requires fresh suite directories for a rerun.

The accompanying evidence archive contains logs, XML results, JSON numerical/coverage/preservation reports, the external validation scripts, and the fresh distribution artifacts. The original source and installed environment were not replaced. **Ready for the next development stage on this tested Mac stack, with the direct Linux-output comparison and statistical-size assessment explicitly outstanding.**
