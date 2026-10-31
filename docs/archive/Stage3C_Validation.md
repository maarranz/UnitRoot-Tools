# Stage 3C validation: Phillips–Perron and KPSS

27 October 2026. **Stage 3C passes on the tested environment.** PP and KPSS
complete the approved initial test family. No methodological blocker was found;
this is not a release or cross-platform compatibility certification. No upstream
package was modified, and no commit or push was made.

## Implemented contract

```python
pp(y, *, trend="c", bandwidth=None, test_type="tau")
kpss(y, *, trend="c", bandwidth=None)
```

PP supports n/c/ct and tau/rho. KPSS supports c/ct. Both use only public released
ARCH 8.0.0 estimation and inference. The active version is checked at runtime.
Imports remain lazy, and the independent augmentation selector neither imports
nor calls ARCH. No development snapshot or private ARCH API is used.

Bandwidth selection has a separate implementation and argument from ADF/DF-GLS
augmentation-order selection. Neither adapter accepts `lags`, `method` or
`max_lags`; `lags` and `lag_selection` in their returned records are None.
For N input levels:

| Test | Final sample, zero-based half-open | Default bandwidth | Fixed feasibility |
|---|---|---|---|
| PP | [1,N), N−1 observations | min(ceil(12(N/100)^(1/4)), N−2) | integer 0≤b<N−1 |
| KPSS | [0,N), N observations | ARCH public Hobijn-style selection (`lags=None`) | integer 0≤b<N |

Minimum lengths are d+3 for PP and d+2 for KPSS, with d=0/1/2 deterministic
regressors for n/c/ct. Data-dependent rank and numerical feasibility also apply.
An explicit zero bandwidth is preserved. Booleans, floats, negative integers
and legacy −1 are rejected. A binding PP cap emits `SpecificationWarning` and
records the uncapped rule. A resolved invalid automatic KPSS bandwidth raises
`EstimationError`; it is not silently capped or replaced.

Result metadata includes effective bandwidth, method (`fixed`, `schwert_capped`,
`hobijn`), requested bandwidth, uncapped PP rule, sample, ARCH/backend versions,
calibration provenance and approximation notes. The Schwert method label names
the capped rule even when its cap does not bind. PP exposes ARCH's public HAC
regression; KPSS returns None without attempting private regression access.
The record and metadata maps are immutable; the backend regression object is
not deeply immutable, as in Stage 3B.

PP's null is a unit root, with strict left-tail rejection. KPSS's null is
stationarity around the stated deterministic terms, with strict right-tail
rejection. Decisions use ARCH critical values at exactly 1%, 5% or 10%, never a
p-value override. Exact critical-value equality does not reject. The summary
names the null and says that failure to reject does not establish it. KPSS
rejection alone does not establish a unit root.

## Test results and coverage

The full suite has **1,164 passing tests, zero failures**: all 797 prior tests
remain passing, plus 367 PP/KPSS cases. The only prior test edit is the necessary
export/public-import boundary update in `tests/test_selection.py`. No earlier
numerical assertion, fixture or selector convention was changed.

Coverage with branch measurement: **458/460 statements (99.57%)** and
**145/146 branches (99.32%)**; combined coverage **99.50%**. The two unexecuted
statements are the new adapter's direct re-raise of a pre-existing
`EstimationError` from regression checks and the pre-existing selector's
non-finite-coefficient guard. Coverage is supporting evidence, not a guarantee
against all numerical problems.

Machine-readable artifacts:

- [Test results](../../validation/stage3c_pytest.xml).
- [Coverage](../../validation/stage3c_coverage.json).
- [Numerical rows](../../validation/stage3c_comparisons.csv) and
  [numerical metadata/group maxima](../../validation/stage3c_numerics.json).
- [Preservation check](../../validation/stage3c_preservation.json).

Tests cover all supported specifications/types; fixed bandwidths 0, 4 and 12;
automatic rules; largest feasible fixed bandwidth and the next invalid value;
short samples; cap warnings; malformed/complex/boolean/nonfinite inputs;
deterministic degeneracy; mocked public numerical failures; invalid returned
bandwidth/sample/inference; wrong backend version; public regression access;
input copying; immutable result metadata; approximation endpoints; and summary
semantics. Independent OLS/Newey-West formulas additionally verify PP with no
deterministic regressors. Exact equality and adjacent floating-point values
exercise both tails, with deliberately conflicting p-values to ensure the
critical-value decision is authoritative. Unexpected programming errors remain
visible rather than being converted indiscriminately.

## Numerical comparisons

`validation/validate_stage3c.py` reproduces **1,944 scalar comparisons**, including
**1,800 expected-equality checks, all passing at absolute tolerance 1e−8**.
These counts are distinct from the historical Stage 2 totals. Repeated fields
on the same six vectors are regression checks, not independent experiments.

| Comparison | Model pairs | Findings |
|---|---:|---|
| Wrapper vs fresh public ARCH | 192 | All statistics, p-values, three critical values, bandwidths and nobs exactly equal; maximum difference 0 |
| Wrapper vs saved released ARCH | 72 | All six saved fields exactly equal; maximum difference 0 |
| PP vs saved R urca | 48 | All 24 tau and 24 rho statistics differ at 1e−8, as expected |
| KPSS vs saved Statsmodels | 24 | Statistics agree; p-values and each of the three critical values differ |

The 192 fresh pairs comprise 108 fixed PP, 36 automatic PP, 36 fixed KPSS and
12 automatic KPSS. The saved ARCH pairs comprise 48 PP and 24 KPSS. Saved
references have c/ct and bandwidths 4/12; new n, zero-bandwidth and automatic
cases are covered by fresh public calls, not invented historical fixtures.
Twelve fresh Statsmodels zero-bandwidth KPSS cases are also checked in the
unit suite; they are not included in the 1,944-row report.

PP/urca maximum absolute statistic differences remain **0.031518998709593404**
for tau and **1.297161193184138** for rho. These reflect the already audited
finite-sample formula conventions, including regression variance/degrees of
freedom, rather than mismatched input vectors or a tolerance failure. Stage 2's
independent package-specific reconstructions and formula explanations remain in
[the numerical study](../Numerical_Validation.md) and
[the implementation specification](Implementation_Specification.md).
No R formula is substituted for ARCH's statistic or calibration.

KPSS/Statsmodels maximum statistic difference is **1.887379141862766e−15**;
maximum p-value difference is **0.03093507074918886**. Their different tabulated
calibrations and interpolation/endpoints explain the inference differences.
The wrapper keeps ARCH's simulated table rather than Statsmodels' calibration.
PP retains statistic-specific ADF-t/tau and ADF-z/rho response surfaces.

## Warnings and numerical limits

The exact ARCH automatic-KPSS historical default-change `DeprecationWarning`
is suppressed because this API explicitly selects the current automatic
behavior. Only that full message/category is filtered; unrelated deprecations
and user warnings remain visible and are tested. Numerical runtime warnings
and Statsmodels' singular-matrix warning become chained `EstimationError`s.

ARCH's Hobijn calculation can fail on very short samples even when the stated
length bound passes; native numerical failures are preserved as causes and
reported rather than replaced with an undocumented rule. Tests reproduce such
native failures for short intercept/trend samples. Invalid automatic returned
bandwidths are additionally tested using controlled public-backend substitutes.
A feasible fixed bandwidth or more observations may resolve a particular case;
no fixed choice is guaranteed to repair degenerate data.

PP p-values at 0/1 and KPSS p-values at 0.0001/1 receive explicit approximation
endpoint notes. These are response-surface/table endpoints, not exact tail
probabilities. Missing values are rejected, never dropped. Equally spaced
observations remain an input assumption. These six simulated vectors do not
establish empirical test size, power or universal floating-point robustness.

## Reproduction and preservation

From the existing repository with the package and test dependencies installed:

```sh
python -m pytest --cov=unitroot_tools --cov-branch \
  --cov-report=term-missing --cov-report=json:validation/stage3c_coverage.json \
  --junitxml=validation/stage3c_pytest.xml
python validation/validate_stage3c.py
python examples/adf_dfgls.py
python examples/pp_kpss.py
python -m pip check
```

Tested: Python **3.14.10**, macOS **26.7.1 arm64**, ARCH **8.0.0**, NumPy
**2.4.4**, SciPy **1.17.1**, Statsmodels **0.15.0**, pytest **9.1.1** and
pytest-cov **7.1.0**. Dependency consistency passes. Both examples execute.
A pure-Python wheel builds without isolation or downloads; a separate import
from that wheel successfully calls ADF, DF-GLS, PP and KPSS. This checks the
built distribution on the existing environment, not a clean dependency install.

SHA-256 checks against the start-of-stage baseline verify all eleven original
`validation/results/` files are byte-identical. Selector source, validation,
errors, version and packaging configuration are unchanged. AST comparison
against the Stage 3B wheel verifies all existing ADF/DF-GLS adapter functions
except the public-class registry are unchanged; selector diagnostic record
classes are unchanged. The public registry now adds KPSS and PhillipsPerron.
The shared result record gains bandwidth metadata and right-tail support.
Historical reports, scripts, fixtures and the ADF/DF-GLS example are preserved.

No live R rerun was performed: urca comparisons use the reviewed saved Stage 2
results. Python 3.11, dependency lower bounds, other platforms/BLAS and a clean
installation matrix remain untested. There is no unresolved methodological
change requiring approval within Stage 3C. The next decision is instructor
review of this completed stage before teaching/release verification; no new
methodological family or inference method has been introduced.
