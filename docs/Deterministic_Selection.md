# Experimental deterministic AIC/BIC selection

`select_deterministic` compares three unrestricted ADF regressions at a supplied,
common augmentation order. It is an experimental specification-selection tool,
not a newly calibrated unit-root test. The existing ADF adapter estimates every
candidate and retains its nonstandard ARCH ADF calibration.

## Model and scores

For original levels y, fit Δyₜ = a + γt + byₜ₋₁ + Σⱼ₌₁ᵖβⱼΔyₜ₋ⱼ + εₜ,
with deterministic dimension d=0,1,2:

| Candidate | Deterministic terms | Regression coefficients q | Likelihood parameters k |
|---|---|---|---|
| n | none | p+1 | p+2 |
| c | intercept | p+2 | p+3 |
| ct | intercept and linear trend | p+3 | p+4 |

The variance counts as an estimated likelihood parameter. All three regressions
use raw levels, the same response observations, lagged levels and differences,
and n=N−1−p effective observations. Response-level indices are `[p+1:N]`, with
zero-based, exclusive-stop boundaries. The ct trend is 1,…,n on this sample.
The winner is not refitted.

For candidate RSS, use the maximum-likelihood variance σ̂²=RSS/n:

- Full conditional Gaussian log likelihood: ℓ=−n/2 [log(2π)+1+log(RSS/n)].
- `aic_score` = n log(RSS/n)+2k.
- `bic_score` = n log(RSS/n)+log(n)k.

Scores are **total** scores omitting the common Gaussian constant n[1+log(2π)].
Thus full AIC/BIC equals the stored score plus that constant. RSS/df_resid is
not the likelihood variance. Both scores are retained regardless of the chosen
criterion. Minimize the unrounded requested score globally; exact floating-point
ties favor n, then c, then ct. There is no near-tie tolerance.

## Public API and example

```python
import numpy as np
from unitroot_tools import select_deterministic

t = np.arange(120, dtype=float)
y = .1 + np.sin(.8*t) + .6*np.cos(1.71*t) + .15*np.sin(.17*t*t)
selection = select_deterministic(y, lags=1)  # BIC default: n
alternative = select_deterministic(y, lags=1, criterion="aic")  # c
print(selection.summary())
print(selection.selected_specification, selection.runner_up_gap)
for row in selection.candidates:
    print(row.specification, row.n_parameters, row.aic_score, row.bic_score)
fixed_spec = selection.selected_fixed_spec_adf
# Fixed-specification ADF p-value reported after selection.
print(fixed_spec.statistic, fixed_spec.pvalue)
```

`lags` is mandatory, nonnegative and integral; zero is valid, booleans and None
are invalid. `criterion` accepts case-insensitive aic/bic. There is no automatic
lag selection, candidate subset, alpha, method or max_lags argument.

`DeterministicSelectionResult` records criterion, selected specification,
common lags and boundaries, input/effective sample sizes, the canonical tuple
of candidates, versions, notes, optional input name, and explicit experimental,
lag-policy, normalization, tie and inference-status metadata. Properties expose
`selected_candidate`, `criterion_value`, read-only `score_differences`,
`runner_up_gap`, and `selected_fixed_spec_adf` (the identical nested fitted result).
Differences use the requested criterion minus its winning score. The runner-up
gap is the second-smallest score minus the smallest, including zero for ties.

Each `DeterministicCandidate` stores specification, deterministic dimension,
coefficient/parameter counts, residual df, RSS, ML variance, full log likelihood,
both scores and the existing `UnitRootResult` in `adf_result`. Scalars are frozen;
version/difference maps are read-only copies. Exposed nested Statsmodels regression
objects remain mutable. Do not mutate them if retaining an auditable record.

## Interpretation and limitations

**Fixed-specification ADF p-value reported after selection.** The selected
model's approximate p-value and critical values use ordinary fixed-specification
ADF calibration. They ignore deterministic selection on the same data and any
upstream lag selection. This is not a selective p-value conditional on the event
that the model won. It does not establish nominal size control for the complete
selection-and-testing procedure. The summary places this qualification alongside
the ADF statistic and p-value, with critical values and calibration provenance.
There is no selector-level pvalue, reject() or unconditional unit-root decision.
Explicitly accessing the nested ADF's reject() remains a fixed-model calculation.
Ordinary regression Student t p-values do not replace ADF calibration.

Including n in every comparison may be inappropriate for economic series.
Selecting n/c/ct does not prove the absence or presence of drift or deterministic
trends. When economic theory gives strong guidance, use a prespecified `adf`
specification instead of mechanical information-criterion selection. Level-origin
changes can change comparisons involving n; translation invariance is not claimed.

Inputs follow the existing finite, real, one-dimensional validation and copy
conventions. Missing data are rejected. Equally spaced observations are required;
calendar spacing is the caller's responsibility. All three fits must be valid,
full rank and have positive finite residual variance. A failed candidate is
identified in the error; candidates are never silently omitted and p is never
reduced. Alignment checks inspect values and named columns as well as metadata.
N≥2p+5 is the algebraic positive-residual-df boundary for ct, not a recommendation
for reliable inference. Existing conditioning warnings remain visible.

## Relation to other procedures and validation

Enders uses sequential tau/phi and conditional deterministic-term tests, with
explicit branch rules; this selector uses a global unrestricted likelihood-score
comparison. It does not change Original or Conservative Enders. Existing
`select_lags` compares augmentation orders using full-series OLS detrending and
a maximum-lag common sample. Its AIC/BIC/MAIC/MBIC scores and modified tau penalty
are not reused here. This feature compares deterministic terms at fixed p.

Run `python -m pytest tests/test_deterministic.py` for focused tests and
`python validation/validate_deterministic.py --report /tmp/deterministic.json`
for independent QR regression, Gaussian-density likelihood and score checks.
Use `PYTHONPATH=src` for a checkout that has not been installed. The script uses
fixed numerical fixtures, including different AIC/BIC winners and ct winners;
it performs no Monte Carlo size/power study. The analytic scoring fixture is an
arithmetic check, not an asserted fitted dataset. Software/numerical agreement
does not establish statistical size. Monte Carlo research under a separately
approved protocol remains future work; bootstrap inference is outside this release.
