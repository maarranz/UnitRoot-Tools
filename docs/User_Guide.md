# UnitRoot-Tools user guide

UnitRoot-Tools provides `adf`, `dfgls`, `pp`, `kpss` and `select_lags` for
teaching applied time-series econometrics. ADF regression estimation uses UnitRoot-Tools' integrated **column-scaled SVD**
OLS implementation. Its p-values and critical values are evaluated by the public
Statsmodels MacKinnon calibration functions, using the ARCH 8.0.0-compatible ADF
conventions. DFGLS, PP and KPSS continue to use released **ARCH 8.0.0**; their
estimation and calibration were not changed by the scaled-SVD integration. The separate `enders`
function provides sequential deterministic-specification decisions; it does
not decide automatically whether economic differencing is appropriate.

## Installation and first run

Python 3.11 or later is required. Consult `Stage4B_Release_Candidate.md` for the
versions actually tested; a declared dependency range is not a tested matrix.
From a local checkout:

```sh
python -m venv .venv
# macOS/Linux:
source .venv/bin/activate
python -m pip install .
python examples/teaching_walkthrough.py
```

For development use `python -m pip install -e '.[test]'` and `python -m pytest`.
For a supplied wheel use `python -m pip install /path/to/unitroot_tools-0.1.0rc1-py3-none-any.whl`.
Then run `python -m pip check`. The version in that filename is the current
unpublished release candidate, not a publication announcement.

In Google Colab, upload the supplied wheel to `/content` using the Files pane,
then run this in a notebook cell:

```python
%pip install /content/unitroot_tools-0.1.0rc1-py3-none-any.whl
```

Restart the runtime if replacing already imported numerical dependencies.
The examples below use ordinary Python and can be pasted into later cells.
Colab/Linux itself was not available for validation in this stage. No PyPI
installation or GitHub release is assumed to exist.

## Start with data and a stated question

```python
import numpy as np
from unitroot_tools import adf, dfgls, pp, kpss, select_lags

y = np.random.default_rng(20261027).standard_normal(240).cumsum()
```

This simulated random walk is reproducible; it is not evidence about test size
or power. Supply a one-dimensional finite numeric series with equally spaced
observations. Missing values are rejected rather than dropped. Select the
sample and deterministic terms using the scientific context, plots and theory;
trying many specifications and reporting the most convenient rejection changes
the interpretation of nominal significance levels.

| Test | Null hypothesis | Supported `trend` | Rejection tail |
|---|---|---|---|
| ADF | stochastic component has a unit root | n, c, ct; ctt with fixed lags only | left |
| DF-GLS | stochastic component has a unit root | c, ct | left |
| Phillips–Perron | stochastic component has a unit root | n, c, ct | left |
| KPSS | stationarity around the stated deterministic terms | c, ct | right |

`n` means no deterministic terms, `c` an intercept, `ct` an intercept and linear
trend, and `ctt` adds a quadratic trend. Trend stationarity is not stationarity
around a constant. Tests under different deterministic specifications answer
different questions.

## ADF and DF-GLS: augmentation order

```python
fixed_adf = adf(y, trend="c", lags=4)
fixed_dfgls = dfgls(y, trend="ct", lags=0)
automatic_adf = adf(y, trend="c", max_lags=12)  # default MAIC
automatic_dfgls = dfgls(y, trend="ct", method="mbic", max_lags=12)
print(automatic_adf.summary())
```

`lags=p` fixes the number of lagged differences, including p=0. Do not combine
fixed `lags` with `method` or `max_lags`. Automatic selection compares candidate
p=0,…,K. An omitted K uses the approved Schwert-style rule and feasibility cap;
a reduced automatic maximum produces a warning. An explicit infeasible order
raises an error. ADF with `ctt` requires fixed `lags`.

The four explicit automatic methods are:

| Method | Candidate score on a common sample of size n |
|---|---|
| AIC | log(σ²ₚ) + 2p/n |
| BIC | log(σ²ₚ) + log(n)p/n |
| MAIC, default | log(σ²ₚ) + 2(τₚ+p)/n |
| MBIC | log(σ²ₚ) + log(n)(τₚ+p)/n |

All use the same full-series OLS detrending and response-level indices
K+1,…,N−1, so n=N−K−1 for every candidate. The fitted selection regression is
Δũₜ=γₚũₜ₋₁+Σⱼ₌₁ᵖφⱼΔũₜ₋ⱼ+eₜ. Here σ²ₚ=Σeₜ²/n and
τₚ=γₚ²Σũₜ₋₁²/σ²ₚ, with the **entire matched lagged-level sum**, including its
first observation. Constants common to candidates are omitted from scores.
Conventional AIC/BIC omit τₚ and are not aliases for MAIC/MBIC.

These approved Perron–Qu OLS selection conventions do not use volatility
rescaling or the bootUR first-level omission. They need not match another
package's one-step autolag. See `MIC_Methodology.md` and the implementation
specification for provenance and limits of paper replication.

The selected integer is passed to ARCH's public fixed-lag implementation on
original levels. DF-GLS performs GLS detrending once inside ARCH; the selector's
OLS transformation is not fed into that final regression. Final sample size
is N−1−p and can exceed the common selection sample size.

```python
selection = select_lags(y, trend="c", method="aic", max_lags=12)
print(selection.selected_lags, selection.selection_nobs)
for candidate in selection.candidates:
    print(candidate.lags, candidate.aic, candidate.bic,
          candidate.maic, candidate.mbic, candidate.tau)
print(automatic_adf.lags, automatic_adf.nobs)
print(automatic_adf.lag_selection.selection_nobs)
```

## PP and KPSS: covariance bandwidth

```python
pp_tau = pp(y, trend="c", bandwidth=4)
pp_rho = pp(y, trend="ct", test_type="rho")
kpss_fixed = kpss(y, trend="c", bandwidth=0)
kpss_auto = kpss(y, trend="ct")
print(pp_rho.summary())
print(kpss_auto.summary())
```

PP supports tau (default) and rho statistics. Their scales and calibrations
are different. Bandwidth controls the Bartlett/Newey–West long-run-variance
estimate; it is **not** an augmentation order. Neither function takes `lags`,
`method` or `max_lags`.

With N input levels, PP defaults to
min(ceil(12(N/100)^(1/4)), N−2); a binding cap warns. Its covariance sample has
N−1 observations. KPSS delegates automatic bandwidth to ARCH's public
Hobijn-style calculation on N observations. Fixed bandwidths are nonnegative
integers, including zero; b<N−1 for PP and b<N for KPSS. KPSS never silently
caps an infeasible automatic answer. Short or degenerate data can fail even
when length bounds pass; choose a justified fixed bandwidth or investigate
the data rather than interpreting an estimation failure as a test outcome.

```python
print(pp_rho.bandwidth, pp_rho.bandwidth_method, pp_rho.bandwidth_uncapped)
print(kpss_auto.bandwidth, kpss_auto.bandwidth_method)
assert pp_tau.lags is None and kpss_auto.lag_selection is None
```

## Inference and interpretation

```python
for result in (fixed_adf, fixed_dfgls, pp_tau, kpss_auto):
    print(result.test, result.null_hypothesis)
    print(result.statistic, result.pvalue, dict(result.critical_values))
    print("Reject stated H0 at 5%:", result.reject(.05))
    print(result.calibration, result.notes)
```

`reject(alpha)` accepts only .01, .05 or .10. It compares the statistic strictly
below the matching ARCH critical value for ADF/DF-GLS/PP, strictly above it for
KPSS. Exact equality does not reject. Approximate p-values remain visible but
do not override this critical-value decision; rounding and calibration
approximations can make the two indications differ near a boundary.

A rejection supplies evidence against the stated null under the maintained
assumptions. Failure to reject is insufficient evidence against it; it does
not prove a unit root or stationarity. KPSS rejection does not by itself prove
a unit root. Conflicting results invite investigation of sample, dynamics,
deterministic terms and assumptions, rather than a mechanical voting rule.
These tests are not an automatic model-selection procedure.

ARCH's own finite-sample critical values, response-surface p-values and KPSS
simulated-table interpolation are preserved. Different packages can produce
different PP statistics or KPSS calibrations even on matched samples.
Endpoint notes identify approximate saturated p-values, not exact probabilities.

## Regression access, diagnostics and errors

```python
print(fixed_adf.regression.params)
print(fixed_dfgls.regression.params)  # regression on GLS-detrended data
print(pp_tau.regression.cov_type)    # HAC regression from public ARCH
assert kpss_auto.regression is None
print(fixed_adf.sample_start, fixed_adf.sample_stop, fixed_adf.versions)
```

The exposed Statsmodels regressions are mutable backend objects. Their ordinary
coefficient p-values are not the calibrated unit-root test p-value. Result
scalars, notes, candidate records and metadata maps are immutable. Sample
boundaries are zero-based with exclusive stop; KPSS's slice denotes the levels
used in its residual/partial-sum calculation, not an exposed regression.

Malformed arguments raise TypeError/ValueError. Unusable numerical estimation
raises `EstimationError`, preserving a numerical cause where available.
`SpecificationWarning` and stored notes explain caps or poor conditioning.
Only ARCH's exact historical automatic-KPSS default-change deprecation notice
is suppressed. Inspect all other warnings and numerical issues before inference.

Run the complete notebook-style sequence locally with
`python examples/teaching_walkthrough.py`. No network or R is needed once the
Python dependencies are installed.

## Sequential Enders procedures

```python
from unitroot_tools import enders
result = enders(y, variant="conservative", max_lags=12)
print(result.summary())
print(result.trend, result.reject_unit_root, result.terminal_test)
```

See [Enders procedures](Enders.md) for the complete API, original versus
conservative rules, finite-table conventions and conditional-inference caveats.
Unrestricted Enders regressions reuse the integrated ADF estimator. Restricted
and auxiliary regressions use the same stable column-scaled SVD fitting layer;
the hypotheses, critical values and decision trees are unchanged.
The conservative variant is a proposed modification; its default preserves
the constant and the nonstandard tau2 decision.

## Experimental deterministic AIC/BIC comparison

Use `select_deterministic(y, lags=p, criterion="bic")` to compare n/c/ct at a
supplied integer p; BIC is the default and AIC is optional. The candidates reuse
the integrated scaled-SVD ADF fits without a separate raw-OLS refit. Every candidate uses
the same original-data response sample N−1−p. The total likelihood scores count
all regression coefficients and the innovation variance. Both scores, the
selected specification, sample boundaries and runner-up gap are retained.
This selects deterministic terms; it does not select augmentation lags or apply
Enders's sequential branch tests.

```python
from unitroot_tools import select_deterministic

deterministic = select_deterministic(y, lags=4, criterion="bic")
print(deterministic.summary())
print(deterministic.selected_specification, deterministic.score_differences)
```

**Fixed-specification ADF p-value reported after selection.** The nested result
uses ordinary ADF calibration, ignores deterministic and upstream lag selection,
and is not a selective p-value conditional on winning. No overall nominal size
guarantee follows. The selector exposes no pvalue or reject() method. Including
n may be inappropriate for economic series and selection does not prove drift
or trend presence/absence; prefer prespecified ADF where theory gives guidance.
See [the detailed methodology and result contract](Deterministic_Selection.md)
for formulas, errors, a reproducible disagreement example and research limits.
