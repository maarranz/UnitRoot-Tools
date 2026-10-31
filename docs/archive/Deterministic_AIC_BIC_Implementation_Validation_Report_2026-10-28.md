# Experimental deterministic AIC/BIC: implementation and validation report

28 October 2026 · UnitRoot-Tools 0.1.0rc1 · MacBook Pro M1 Pro

## Assessment

The approved experimental selector is implemented and passes software and independent
numerical validation on the source checkout, a fresh wheel installation and a fresh
source-distribution installation. Each passes **1,363 tests**: the 1,287-test baseline
plus **76 new tests**. The existing export-list assertion was extended for the three
approved exports. No existing econometric implementation was changed.

Independent QR reconstruction agrees within the fixed 1e−9 absolute tolerances;
the largest likelihood/score/margin discrepancy is 1.4210854715202e-13.
All numerical selector reports agree exactly across the three paths. The complete
Enders numerical reports and 84-case variant-equivalence reports agree exactly with
the earlier saved Mac baseline for each path.

Ready for user review and a separately authorized commit as an **experimental
specification-selection feature**. Statistical size and power of the whole
selection-and-testing procedure remain unassessed. Nothing was committed, tagged,
pushed, published, or installed into the active environment.

## 1. Repository and environment

Repository: `/Users/maarranz/GitHub/UnitRoot-Tools`.
Initial and final branch: `main`; initial and final HEAD:
`8c4f1a17cf12207dba3c051be6f1de7943df895c` (following Enders commit `329468c`).
Initial tracked working tree was clean; the four untracked design/audit/validation
artifacts listed in §9 were already present and remain byte-identical.
The approved design is `docs/Deterministic_AIC_BIC_Design.md`; the user's implementation
assignment supplies the final fixed-specification terminology.

| Component | Observed version / configuration |
|---|---|
| Hardware | MacBook Pro MacBookPro18,3; Apple M1 Pro; 8 cores; 16 GB RAM |
| OS | macOS 26.7.1 (25G241), Darwin 25.6.0 |
| Architecture | native arm64 |
| Python | 3.12.16, Clang 23.1.3 |
| Package | 0.1.0rc1, unchanged |
| ARCH | 8.0.0 |
| NumPy / SciPy | 2.4.4 / 1.17.1 |
| Statsmodels / pandas | 0.15.0 / 3.1.0 |
| pytest / pytest-cov | 9.1.1 / 7.1.0 |
| setuptools / pip | 84.2.1 / 26.2.2 |
| NumPy BLAS/LAPACK | macOS Accelerate; library version not reported |

Interpreter used throughout:
`/Users/maarranz/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/stage4b/env312-current/bin/python`.
Its previously installed same-version package predates Enders; explicit import-path
checks therefore matter. Source checks used the repository's absolute `src` path;
installed checks used isolated `--target` directories and asserted import location
and runtime byte equality. No dependency upgrades or environment recreation occurred.

## 2. Files and additive architecture

| File | Change |
|---|---|
| `src/unitroot_tools/_deterministic.py` | New selector, frozen result types, exact comparator, scoring and alignment safeguards |
| `src/unitroot_tools/__init__.py` | Export selector and two result types |
| `tests/test_deterministic.py` | 76 new tests |
| `tests/test_selection.py` | Extend exact public-export assertion; retain all previous exports and backend-boundary checks |
| `validation/validate_deterministic.py` | Standalone independent QR/Gaussian-density validator with JSON report |
| `docs/Deterministic_Selection.md` | Methodology, API, reproducible example, interpretation and limits |
| `README.md` | Experimental feature overview and documentation link |
| `docs/User_Guide.md` | Usage section and inference qualification |

The new module composes the existing fixed-lag `adf()` for n/c/ct. It introduces no
second ADF engine. Existing `results.py`, `_adapters.py`, `_enders.py`, Enders tables,
`selection.py`, `_bandwidth.py`, validation helpers, old validation scripts and
fixtures are byte-identical. Version and dependency metadata are unchanged.

## 3. Public API and result contract

```python
select_deterministic(y, *, lags: int, criterion: str = "bic")
```

`lags` is mandatory and nonnegative; NumPy integral values are accepted, bool/None
and nonintegral values rejected. Criterion normalization follows existing
case-insensitive aic/bic validation. There are no automatic-lag, max_lags, method,
alpha or candidate-subset arguments.

Public `DeterministicCandidate` fields: specification, deterministic_dimension,
n_coefficients, n_parameters, df_resid, ssr, variance_mle, log_likelihood,
aic_score, bic_score and adf_result (the existing `UnitRootResult`).

Public `DeterministicSelectionResult` fields: criterion, selected_specification,
lags, input_nobs, nobs, sample_start, sample_stop, canonical n/c/ct candidates,
versions, notes, input_name, lag_policy, criterion_normalization, tie_break,
inference_status and experimental. Defaults identify the supplied-common-integer
lag policy, total reduced-score normalization including variance, exact-score tie
rule and fixed-specification calibration after selection.

Properties expose selected_candidate, criterion_value, runner_up_gap,
score_differences (read-only requested-score differences from the minimum), and
selected_fixed_spec_adf. The selected nested object is retained by identity;
it is never refitted. Records are frozen and maps copied/read-only; nested backend
regression objects remain mutable, consistently with the existing result contract.
There is **no selector-level pvalue, reject() or reject_unit_root**.

The summary shows the experimental status, criterion/fixed p, common sample,
all scores/counts/differences, winner/margin and fixed-specification ADF statistic,
approximate p-value, critical values and calibration/provenance. Beside those
numbers it states:

**Fixed-specification ADF p-value reported after selection.** Ordinary ADF calibration
ignores deterministic selection and any upstream lag selection. This is not a
selective p-value conditional on winning and establishes no nominal size control
for the whole procedure. Ordinary coefficient Student t inference is not substituted.

## 4. Methodological conventions and safeguards

All candidates fit unrestricted ADF regressions to identical raw-data response,
lagged-level and lagged-difference observations at fixed p. For input length N,
n=N−1−p, with zero-based response-level slice `[p+1:N]`. The ct trend is 1,…,n.
Checks compare actual values, named columns, rank, sample metadata and residual df.

With deterministic dimension d=0,1,2, q=p+d+1 and **k=q+1=p+d+2** including variance:

- ML variance = RSS/n, not RSS/df_resid.
- Full conditional Gaussian ℓ = −n/2 [log(2π)+1+log(RSS/n)].
- AIC score = n log(RSS/n)+2k.
- BIC score = n log(RSS/n)+log(n)k.

Both total scores omit n[1+log(2π)], common to all candidates. No lag-selection or
MAIC/MBIC score is reused. The global unrounded minimum wins; exact floating-point
ties prefer n, then c, then ct. One-ULP improvements win without a near-tie tolerance.

All candidates must succeed. Failures preserve the project's exception categories
and identify the offending candidate; no omission or lag reduction occurs.
Existing full-rank/positive-residual-variance and conditioning safeguards are
retained, with additional finite likelihood/score and alignment checks.
N≥2p+5 is only an algebraic feasibility boundary. No Enders phi-table n≥25 limit
is imposed. Inputs are copied; missing/nonreal/nonfinite/multidimensional data
follow existing validation. Calendar spacing remains the caller's responsibility.

## 5. Test execution

| Path | Suite | Tests | Failures | Errors | Skipped | Seconds |
|---|---|---:|---:|---:|---:|---:|
| source | full | 1363 | 0 | 0 | 0 | 5.990 |
| source | focused | 76 | 0 | 0 | 0 | 0.973 |
| source | enders | 123 | 0 | 0 | 0 | 1.078 |
| wheel | full | 1363 | 0 | 0 | 0 | 5.955 |
| wheel | focused | 76 | 0 | 0 | 0 | 1.010 |
| wheel | enders | 123 | 0 | 0 | 0 | 1.049 |
| sdist | full | 1363 | 0 | 0 | 0 | 5.981 |
| sdist | focused | 76 | 0 | 0 | 0 | 0.991 |
| sdist | enders | 123 | 0 | 0 | 0 | 1.045 |

Each final suite emitted zero uncaptured test warnings; warnings are configured as
errors. An explicit focused test confirms inherited SpecificationWarning propagation.
The 76 tests cover independent least-squares calculations; p=0/1/3; exact sample
values/names; parameter counts; real AIC/BIC disagreement and ct winners; analytic
scoring; exact/one-ULP ties; scale shifts; invalid arguments/data; degenerate fits;
positive-df sample boundaries; all-or-error candidate failures; corrupted response,
regressor and metadata guards; no input mutation; no winner refit; nested ADF
calibration consistency; provenance; summary wording and immutable scalar containers.

The first complete run had 1 failure and 1,362 passes because an existing test
hard-coded the old export set. Only that assertion was extended for the additive
API; no test was removed or weakened. The final complete run has 1,363 passes.
Initial failure logs/JUnit are retained separately in the evidence archive.

Source full-run coverage: new module 125/126 executable statements (99.21%) and
20/22 branches (90.91%); coverage is supporting software evidence, not statistical
validation. `git diff --check` passed.

## 6. Independent numerical evidence

`validation/validate_deterministic.py` independently constructs every design row by
explicit time indices, solves QR with NumPy, calculates residuals/RSS, and sums
SciPy Gaussian log densities at RSS/n. It converts the independently obtained
full likelihood to reduced IC scores. It does not use selector helpers to build
expected results. Two fixed analytic level vectors at p=0/1/3 give 18 distinct
regressions and 12 AIC/BIC selection calls per distribution. These are numerical
fixtures, not size/power simulations.

Absolute tolerances were fixed at **1e−9**, stricter than the design's proposed
1e−8 ceiling; no tolerance was loosened. Indices, counts, winner labels, alignment,
and nested fixed-specification ADF calibration are checked exactly. Maximum errors
below apply identically to source, wheel and sdist:

| Quantity | Maximum absolute discrepancy | Absolute tolerance |
|---|---:|---:|
| coefficients | 4.884981308350689e-15 | 1e−9 |
| ssr | 1.4210854715202e-14 | 1e−9 |
| variance_mle | 2.220446049250313e-16 | 1e−9 |
| log_likelihood | 5.684341886080801e-14 | 1e−9 |
| aic_score | 9.947598300641403e-14 | 1e−9 |
| bic_score | 9.237055564881302e-14 | 1e−9 |
| runner_up_gap | 1.13686837721616e-13 | 1e−9 |
| score_differences | 1.4210854715202e-13 | 1e−9 |

Actual fitted-score example at p=1 (all values rounded here only for presentation):

| Fixture | Criterion | n score | c score | ct score | Winner | Runner-up gap |
|---|---|---:|---:|---:|---|---:|
| disagreement | AIC | -118.499125490536 | -119.694237130157 | -117.728255367299 | c | 1.195111639621 |
| disagreement | BIC | -110.187071617139 | -108.611498632294 | -103.874832244971 | n | 1.575572984845 |
| trend | AIC | -35.786618470133 | -51.189167348475 | -117.728255367299 | ct | 66.539088018824 |
| trend | BIC | -27.474564596736 | -40.106428850612 | -103.874832244971 | ct | 63.768403394358 |

The reproducible disagreement fixture is yₜ=.1+sin(.8t)+.6cos(1.71t)+.15sin(.17t²),
t=0,…,119. Adding .03t gives the trend fixture. Independent results confirm ct
wins under both criteria at each tested lag order.

Section 9.1's separate arithmetic fixture uses n=100, p=1, k=(3,4,5) and
RSS=(100,100exp(−.03),100exp(−.035)). AIC=(6,5,6.5) selects c;
BIC=(3log100,4log100−3,5log100−3.5) selects n. The focused check uses 1e−12
absolute tolerance. These prescribed RSS values are not misrepresented as a
fitted dataset. Exact/near-tie tests use deliberately equal floating-point scores
and a one-ULP difference, independent of numerical-fit tolerances.

## 7. Enders and baseline preservation

Each path reran all 123 focused Enders tests and the existing standalone validator
(18 independent-regression cases, unchanged 1e−8 numerical tolerance). Every field
of its JSON report equals the saved earlier Mac report. Every field of the 84-case
Original/permitted-Conservative equivalence report also equals its saved baseline:
complete histories, selected lag orders, samples and fixture decisions are retained.
These extra fixture checks contain no Monte Carlo size/power experiment.

Published real-GDP example: n=262, p=1, selected c. Original rejects through its
conditional branch; default Conservative does not reject at tau2; Conservative
with allow_no_constant=True agrees with Original. Tau2 remains
−2.048529930012573 against the fixed-specification 5% value −2.872635586277424;
Original's conditional left-tail comparison uses −1.6507581336531618.

The pre-existing published rounding limitation remains explicit: A2.3 intercept
ordinary t=2.646058967065688 versus printed 2.64 differs by 0.00605896706568787,
outside the published-statistic rounding tolerance 0.005. The validator records
that exception; it is not a new implementation discrepancy and was not concealed
by changing tolerances. The six published branch statistics reproduce to the
existing 0.005 rounding convention.

The 1,287-test baseline is retained within the expanded suite; the historical Linux
baseline is not a new Linux run. Passing software tests and matching published
calculations do not establish nominal size of Enders or the new selector.

## 8. Fresh distribution validation

Artifacts were built from a clean external `git archive HEAD` snapshot overlaid
with the eight authorized changed/new feature paths in §2. The four pre-existing
untracked audit/design/evidence files were excluded from this validation-only
package snapshot. No repository dist artifact or active environment was overwritten.

First build the sdist using setuptools.build_meta in that snapshot; extract it;
then build the wheel from the extracted sdist. Install wheel and sdist separately
with --target, --no-deps, --no-build-isolation, --no-index and --no-cache-dir.
The sdist path therefore includes an independent local wheel build by pip.

For each installed path, tests and validation fixtures were copied to a separate
suite directory; its pytest configuration has **no pythonpath=src**. PYTHONPATH
names only that installation. Import assertions prove the active module and all
runtime module bytes match the tested snapshot. Both full suites, both focused
suites, both independent scripts, 84-case Enders equivalence checks and new
documentation examples passed. Both pip dependency checks passed.

Distribution inspection verified metadata/version/dependencies, runtime bytes,
license, documentation/examples and required sdist tests/fixtures. The new
methodology page is in wheel and sdist; the new validator and focused tests are in
the sdist and executed externally against the wheel. No release was published.

| Artifact (version unchanged) | Bytes | SHA-256 |
|---|---:|---|
| unitroot_tools-0.1.0rc1-py3-none-any.whl | 115168 | `7f96d670dc822184a23c79bc65379fddd7b357eaab89c481cf17f6cc863b49cb` |
| unitroot_tools-0.1.0rc1.tar.gz | 269052 | `df2f6e46081f490a8a4438afc0df312d8f0d8a78f82768f69a6cd0996a750bf0` |

Nonfatal diagnostic/build messages are retained: NumPy configuration display asks
for optional pyyaml; system_profiler reports an hw.cpufamily lookup error while
successfully reporting the M1 hardware; the existing manifest's `validation/*.xml`
glob matches no files; byte-compilation is disabled for clean external builds.
No dependency was added to remove these messages. They are separate from test
warnings/failures, which are zero in final runs.

## 9. Preservation and final working tree

Before/after SHA-256 comparison covered 719 existing non-Git/non-bytecode files.
Only four pre-existing files changed: README, User Guide, public exports and the
export-list assertion. No existing file was removed. The four pre-existing
untracked artifacts retain their exact bytes and hashes (`preservation.json`).
No existing econometric source, test fixture, historical distribution, audit or
validation evidence was overwritten. pyproject.toml and package version bytes
are unchanged. There was no pull, commit, tag, push or active-environment install.

Final working-tree status:

```text
## main...origin/main
 M README.md
 M docs/User_Guide.md
 M src/unitroot_tools/__init__.py
 M tests/test_selection.py
?? docs/Deterministic_AIC_BIC_Design.md
?? docs/Deterministic_Selection.md
?? docs/Enders_Mac_M1_Pro_Validation_Evidence_2026-10-28.zip
?? docs/Enders_Mac_M1_Pro_Validation_Report_2026-10-28.md
?? docs/Theoretical_Audit_Sequential_Strategies.md
?? src/unitroot_tools/_deterministic.py
?? tests/test_deterministic.py
?? validation/validate_deterministic.py
8c4f1a17cf12207dba3c051be6f1de7943df895c
```

## 10. Reproduction, evidence and next steps

Evidence archive: `Deterministic_AIC_BIC_Implementation_Validation_Evidence_2026-10-28.zip`.
It contains exact argument arrays, environment overrides and working directories
in commands.json and extra-commands.json, complete logs/JUnit/coverage JSON,
independent reports, preservation hashes, artifact inspection, tracked patch,
validation orchestration scripts and fresh wheel/sdist. Full feature source and
tests are recoverable from the sdist; uncommitted state is identified explicitly.
Hardware identifiers are redacted in the archive's hardware log.

Core source commands (P is the interpreter in §1; R the repository; W an external
writable evidence directory):

```sh
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export PYTHONPATH="$R/src" COVERAGE_FILE="$W/source-coverage.data"
cd "$R"
"$P" -m pytest tests -p no:cacheprovider --junitxml="$W/source-full.xml" --cov=unitroot_tools --cov-branch --cov-report=json:"$W/source-coverage.json"
"$P" -m pytest tests/test_deterministic.py -p no:cacheprovider --junitxml="$W/source-focused.xml"
"$P" -m pytest tests/test_enders.py -p no:cacheprovider --junitxml="$W/source-enders.xml"
"$P" validation/validate_deterministic.py --report "$W/source-numerical-deterministic.json"
"$P" validation/validate_enders.py --report "$W/source-numerical-enders.json"
```

For exact fresh-build and installed commands, use the archived command records;
run_validation.py and extra_checks.py document the complete procedure and paths.
They deliberately create new validation directories and are not intended to
silently overwrite an earlier run. Dependencies must already be present as in §1.

Review the additive API, numerical conventions, new documentation and working-tree
diff before authorizing a commit. A subsequent separately approved research
protocol should assess selection frequencies and size/power for the full procedure.
For now: mandatory supplied p, fixed n/c/ct, no selective calibration, no automatic
lag selection, no new sequential methods, no Monte Carlo run and no bootstrap.
The economic warning is explicit: including n may be inappropriate; selection
proves neither drift nor trend presence/absence, and strong theory favors a
prespecified ADF specification. These limitations do not prevent review of this
experimental software implementation but preclude treating it as a calibrated
new unit-root test.
