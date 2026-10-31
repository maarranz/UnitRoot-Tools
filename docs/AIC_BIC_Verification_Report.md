# AIC/BIC deterministic-specification verification

31 October 2026 · Mac M1 Pro · qualified Python 3.13.18

**Verification passes; no production correction is needed.** Fixed-lag AIC/BIC selection consistently uses the integrated scaled-SVD ADF fits. Scientific definitions, candidate sets, decisions and API remain unchanged. This is bounded numerical verification, not post-selection inference calibration.

## Scope and implementation review

Repository `/Users/maarranz/GitHub/UnitRoot-Tools`, branch `main`, HEAD `37fa6efab21a5f332c062c263007411cdb46623e`. Reviewed `_deterministic.py`, `_adapters.py`, `_regression.py`, deterministic/helper/ADF/Enders tests, saved fixtures, prior deterministic implementation/correction documentation and the ADF/Enders scaled-SVD integration reports.

The selector performs exactly three unrestricted fits, in order n/c/ct, at one explicitly supplied common nonnegative integer lag p. It calls the integrated public ADF, which uses column-scaled SVD and returns a genuine Statsmodels regression result. No raw-OLS fallback or re-estimation is introduced by selection. Alignment checks require the actual response, named columns, sample boundaries and residual degrees of freedom to agree. Quadratic trend is not a candidate. Automatic lag selection is not invoked.

For N input levels, all candidates have n=N−1−p response observations, level indices [p+1,N). Their response is the same first difference; regressors are lagged level, the same p lagged differences, then the candidate's constant and linear trend (1 through n). With deterministic dimension d=0/1/2, coefficient count q=1+p+d, total likelihood parameter count k=q+1 includes the estimated innovation variance, and residual df=n−q. All candidates must succeed; failures retain candidate context rather than silently dropping a model. The minimum algebraic sample N=2p+5 leaves ct with one residual degree of freedom.

## Independent definitions and conventions

For Gaussian conditional regression errors, maximization over variance gives sigma²_ML=RSS/n. Substitution in the Gaussian log density gives

- log L = −n/2 × [log(2π)+1+log(RSS/n)].
- AIC = n log(RSS/n)+2k.
- BIC = n log(RSS/n)+log(n)k.

The returned scores omit the common n[1+log(2π)] term from −2 log L. Their normalization is explicitly recorded; their ranking is unaffected because n and the response are identical across candidates. The full returned log likelihood includes this constant. Coefficient covariance and t statistics use RSS/(n−q), while IC likelihood uses RSS/n. This deliberate distinction is checked directly.

[R stats AIC documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/AIC.html) defines AIC/BIC penalties and requires comparable maximum-likelihood fits to the same observations. [R logLik documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/logLik.html) documents estimated scale and full constants for linear models. [Statsmodels OLS.loglike](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLS.loglike.html) documents concentrating out scale. [Statsmodels info_criteria](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.RegressionResults.info_criteria.html) documents that its default omits scale from the parameter count; `dk_params=1` includes it. Thus comparing the selector's scores directly with default backend AIC/BIC would mix conventions. New checks explicitly include scale and remove the common Gaussian constant. These independent sources confirm the existing contract; they do not motivate a formula or policy change.

Exact score ties favor smaller deterministic dimension, with no tolerance or rounding applied; a one-ULP advantage wins even for ct. Existing tests verify order independence, exact/near ties, invalid lags and inputs, alignment failures, sample boundaries, degeneracy, finite-positive variance and inherited warning propagation. Result fields correctly expose both scores, selected candidate/specification, requested-criterion value, unrounded differences and runner-up gap. The selected fixed-specification ADF object is reused. Frozen scalar records and read-only mappings retain their existing contracts; nested regression objects remain mutable as documented.

Notes and summary retain experimental status, fixed-specification calibration after selection, absence of a selector-level p-value or rejection decision, upstream lag-selection qualification, equally spaced data responsibility and the caution about including n for economic data. Conditioning warnings are inherited, not suppressed. No selective-inference or nominal-size claim is added.

## New numerical evidence

`validation/verify_deterministic_scaled.py` verifies the existing checksum-pinned hexadecimal binary64 saved inputs (`y`, `positive`, `negative`), at lags 0/1/4 and both criteria: **18 selection calls, 27 distinct candidate regressions**. Positive and negative fixtures have magnitude 1e-12. It independently constructs response and design by explicit time indices, uses the existing independent Decimal elimination reference routine at 80 and 120 digits, and computes variance, Gaussian likelihood and IC scores with Decimal logarithms. Both precisions produce identical binary64 reference outputs. High-precision normal equations are reference arithmetic only; production continues to use SVD.

Maximum errors across these cases:

| Quantity | Maximum error | Fixed test limit |
|---|---:|---:|
| RSS, relative | 3.031e-16 | 1e-12 |
| ML variance, relative | 2.723e-16 | 1e-12 |
| Full log likelihood, absolute | 4.547e-13 | 4e-12 |
| AIC score, absolute | 9.095e-13 | 4e-12 |
| BIC score, absolute | 9.095e-13 | 4e-12 |

Every selected specification matches independent ranking. Gaps and differences agree within 8e-12. Actual designs/responses and integer counts agree exactly. RSS equals the returned stable regression's RSS; variance conventions are checked separately. Backend concentrated likelihood agrees within 2e-12; backend IC with scale included agrees within 4e-12 after constant removal.

The existing small-scale tests retain their historical 1e-5 raw-fit RSS budget and explanatory comment unchanged. That budget predates ADF SVD integration; the 18 new strict tests supplement it with independent references rather than relaxing, replacing or rewriting any earlier fixture, expectation or tolerance.

## Baseline and bounded comparisons

Before additions, the complete package/research suite passes **1,597 tests** (1,571 package + 26 research). The unchanged existing independent QR deterministic validator passes 12 selections/18 fits, the independent Enders validator passes its 18 cases and GDP reproduction, and saved Enders/ADF high-precision references reproduce exactly (18 and 4 cases respectively).

The bounded probe uses 19 saved vectors: three exact binary64 numerical vectors, six reviewed CSV series, GDP and nine portable historical vectors. Each is evaluated at p=0/1/4 for AIC and BIC: **114 selections / 171 distinct fits**. Baseline and final staged records are exactly equal, including input hashes, all candidate scores/RSS/likelihoods/variance/counts, selected specifications, score differences/gaps, ADF statistics/p-values, sample and lag fields, notes and warnings. The 114 compact score/selection records also match the previous Enders integration's saved evidence exactly. Historical portable inputs are reused on this Mac; no Linux qualification or new simulation is performed.

## Final validation and delivery

Final full source and separate installed-wheel suites each pass **1,615 tests**: 1,589 package + 26 research, including the 18 new strict checks. Focused existing/new deterministic tests pass **98 tests**. Existing independent validators and new high-precision score verification pass against source and the separate wheel target. Installed bounded records equal source exactly. All installed runtime modules match source bytes. The new report is included in the wheel documentation.

The wheel is built from an external clean source copy with no dependency resolution or build isolation and installed into an external target, preserving the qualified environment and its historical installed package. Tests use explicit source or wheel paths, one numerical thread, no bytecode or pytest cache writes, and the existing warnings-as-errors policy. No warning suppression, skip or scientific tolerance change is used.

Exactly three repository additions belong to this assignment:

- `tests/test_deterministic_scaled_verification.py` — 18 focused independent checks.
- `validation/verify_deterministic_scaled.py` — reproducible independent score validator.
- `docs/AIC_BIC_Verification_Report.md` — this permanent report.

No pre-existing file is edited. All 836 initial repository files, including tracked/untracked/ignored evidence and permanent docs, retain their bytes and modification times; all 8,106 qualified-environment files are preserved with no additions. Branch and HEAD remain unchanged. Preservation inventories, baseline/staged copies, command records, JSON numerical evidence, test XML/logs and the wheel/installed target remain outside the repository at `/Users/maarranz/Documents/ChatGPT/UnitRoot-Tools/work/aic-bic-verification-2026-10-31-01/`. The preservation result and artifact checksums accompany that evidence. `git diff --check` passes.

No policy decision is required by the tested evidence. Next recommendation: review the combined ADF/Enders changes and verification evidence, then prepare a separately authorized scoped Git commit including permanent reports. No RNG changes, automatic-lag implementation, pilot/main Monte Carlo, Colab/Linux qualification, commit, push or tag occurred.
