# MAIC and MBIC methodology — Stage 2A

27 October 2026. This document specifies a recommended procedure; it does not implement a UnitRoot-Tools API. Source inspection precedes the numerical comparisons in `Numerical_Validation.md`.

## Recommendation

Adopt **Perron–Qu OLS-based lag selection, without volatility rescaling by default**, with MAIC as the default criterion and MBIC as an explicit alternative. Detrend the complete observed series by OLS, estimate all candidate difference regressions on one common sample, use the maximum-lag sample count in every residual-variance/criterion denominator, and sum the squared lagged levels over **all those same rows**. For DF-GLS, select the lag using this OLS series, then estimate the final statistic on GLS-detrended data at that fixed lag.

For a future implementation, use a small independently specified selector shared by MAIC and MBIC, followed by the released ARCH fixed-lag test backend. Use the pinned ARCH development DFGLS selector as a reference for MAIC, not as an unversioned production dependency. Its selector lacks MBIC. Do not adopt bootUR's implementation verbatim: the inspected release omits the first matched lagged-level value from its MIC penalty. Do not silently emulate that discrepancy. Keep original GLS-selection and volatility-rescaled criteria as separately named research options if later required.

This recommendation selects a reproducible convention; it does not assert universal finite-sample superiority. A small cross-validation set cannot establish size or power.

## 1. Sources and version boundary

- **Ng and Perron (2001)**, *Econometrica* 69(6), 1519–1554. Full published paper inspected, especially Section 3 (GLS detrending), Section 5 (MIC), equations (12)–(13), and discussion of candidate samples. Equation page 1529 was also rendered and visually checked. [Published-paper PDF](https://users.ssc.wisc.edu/~bhansen/718/NgPerron2001.pdf), [DOI](https://doi.org/10.1111/1468-0262.00256).
- **Perron and Qu (2007)**, *Economics Letters* 94(1), 12–19, DOI [10.1016/j.econlet.2006.06.009](https://doi.org/10.1016/j.econlet.2006.06.009). Publisher abstract explicitly identifies replacing GLS detrending with OLS detrending when constructing MIC. Author's [publication list](https://sites.bu.edu/qu/publications/) confirms 2007. The February 2006 working-paper abstract was located, but attempts to retrieve its full PDF failed (a nominal PDF download was actually HTML). This audit therefore does **not** claim to have checked every detail of the full Perron–Qu article. The OLS-selection convention is corroborated directly by both implementations below. [Publisher abstract](https://www.sciencedirect.com/science/article/abs/pii/S0165176506002138).
- **ARCH 8.0.0** is the released baseline; **development commit `09a6be5f1b490346d8ba5d3b5850991734aaaf3a`** is the MAIC baseline, not a release. Inspected `arch/unitroot/unitroot.py`, especially `_select_best_ic`, `_autolag_ols`, `_autolag_ols_low_memory`, `_df_select_lags`, `ADF` and `DFGLS`. [Pinned source](https://github.com/bashtage/arch/blob/09a6be5f1b490346d8ba5d3b5850991734aaaf3a/arch/unitroot/unitroot.py).
- **bootUR 1.0.5**, CRAN source archive: inspected `src/bootURfunctions.cpp`, `R/adftest.R`, `R/bootURtests.R`, `R/auxiliary.R`, `R/RcppExports.R`. Relevant C++ functions are `de_trend`, `adf_cpp`, `maic_cpp`, `mbic_cpp`, `adf_selectlags_cpp`, and the `adf_onestep_*` variants. [Versioned archive](https://cran.r-project.org/src/contrib/bootUR_1.0.5.tar.gz), [reference manual](https://cran.r-project.org/web/packages/bootUR/refman/bootUR.html). Archive digest is recorded below.

bootUR labels the simple modification “Perron and Qu (2008)” and misspells the journal in one source reference. The matching title, volume and pages belong to the **2007** article above. No distinct 2008 modification is assumed.

## 2. Original Ng–Perron definitions

In the paper's notation, for a fixed maximum candidate order $k_{max}$, fit the augmented regression on detrended data $\widetilde y$:

$$
\Delta\widetilde y_t
=\widehat\beta_0(k)\widetilde y_{t-1}
 +\sum_{j=1}^{k}\widehat\beta_j(k)\Delta\widetilde y_{t-j}
 +\widehat e_{t,k}.
$$

On the common candidate sample $t=k_{max}+1,\ldots,T$, define

$$
n_c=T-k_{max},\qquad
\widehat\sigma_k^2=\frac{1}{n_c}\sum_{t=k_{max}+1}^{T}\widehat e_{t,k}^2,
$$

$$
\tau_T(k)=\frac{\widehat\beta_0(k)^2}{\widehat\sigma_k^2}
             \sum_{t=k_{max}+1}^{T}\widetilde y_{t-1}^2.
$$

Then equations (12) and (13) give

$$
\operatorname{MIC}(k)=\log\widehat\sigma_k^2
 +\frac{C_T}{n_c}\{\tau_T(k)+k\},
\qquad
\widehat k=\arg\min_{0\leq k\leq k_{max}}\operatorname{MIC}(k).
$$

**MAIC uses $C_T=2$; MBIC uses $C_T=\log(T-k_{max})=\log n_c$.** Using $\log T$ is a different finite-sample convention. The variance divisor is the common sample count, not residual degrees of freedom and not a count changing with candidate k. The level coefficient is estimated freely; do not fit a restricted regression with that coefficient set to zero when computing the correction term.

The paper develops MIC for detrended series generally and combines it with GLS detrending in its recommended efficient tests. “Original Ng–Perron” should therefore be qualified: the **GLS-selection** version uses GLS-detrended data for both selection and the final GLS-based statistic. An MIC formula alone does not determine which detrending procedure was used.

## 3. Detrending and first observations

Let $d_t=(1)$ or $(1,t)'$. OLS detrending is

$$
\widehat\delta_{OLS}=\arg\min_\delta\sum_t(y_t-d_t'\delta)^2,
\qquad v_t^{OLS}=y_t-d_t'\widehat\delta_{OLS}.
$$

GLS detrending uses $\bar\alpha=1+\bar c/N$, with $\bar c=-7$ for an intercept and $-13.5$ for intercept/trend in the inspected implementations. Transform the first observed y and deterministic row **without differencing**, and for subsequent rows use $y_t-\bar\alpha y_{t-1}$ and $d_t-\bar\alpha d_{t-1}$. Regress transformed y on transformed d, then subtract the resulting deterministic fit from the original levels.

The paper introduces its transformation on a sequence indexed from 0 through T. ARCH, urca and bootUR operate on the supplied N observed levels. Their finite-sample coefficient uses N, and the first supplied value is the retained initial row. Do not confuse a paper's asymptotic/index notation with a claim that N supplied values provide an extra unobserved initial level.

**Perron–Qu operational convention:** replace the selection series $v^{GLS}$ with $v^{OLS}$; retain GLS detrending for the final DF-GLS/M statistic and its nuisance estimation as specified for that statistic. Changing the selection series does not turn the final test into an OLS ADF test. bootUR's public `adf(..., two_step=TRUE)` reports a two-step **OLS** ADF statistic; its selected order can still be compared with ARCH DFGLS. Those final test statistics are not the same comparison.

## 4. Exact array/sample convention for adoption

For N supplied observations $y_1,\ldots,y_N$, choose a declared maximum K. Detrend **all N levels**, then compare k = 0,…,K on

$$
 t=K+2,\ldots,N,\qquad n=N-1-K.
$$

At every k use the same response $\Delta v_t$, level $v_{t-1}$, and first k lagged differences. The matched level sum is

$$
 S=\sum_{j=K+1}^{N-1}v_j^2.
$$

In zero-based Python slicing this is `v[K:-1]`; response is `diff(v)[K:]`. These two vectors must have equal length n. Set $\widehat\sigma_k^2=SSR_k/n$, $\tau_k=\widehat\beta_{0,k}^2S/\widehat\sigma_k^2$, and use 2 or log(n) in MIC. Resolve exact ties toward the smallest k.

After selecting k, refit the final test using **N−1−k** observations. That final sample is generally larger than the selection sample. Do not replace the selection denominator by the final sample size. Fix K explicitly for cross-package validation: defaults differ (ARCH ceiling rule; bootUR rounded rule with a short-sample adjustment).

For future production input handling, reject non-finite data, rank-deficient candidate designs and nonpositive residual variance; specify a feasible K rather than silently changing it mid-comparison. This is a recommendation, not implemented library behavior.

## 5. What ARCH development implements

- `_select_best_ic`: `log(sigma2) + 2*(tau+k)/nobs`; **MAIC only**, no MBIC.
- `DFGLS`: OLS-detrends the complete series for selection, passes `trend="n"` to the selector, and runs the final fixed-lag regression on GLS-detrended levels. This matches the recommended PQ selection structure.
- Full-memory selector: the squared QR diagonal of the level column gives its sum of squares after preceding deterministic regressors are projected out. With pre-detrended input and no deterministic regressors, this is precisely the full matched sum S above.
- Low-memory selector: uses a Schur complement for the corresponding level sum; normalization is compensated in coefficients/variance, and the criterion's variance scale is restored. Both memory routes need numerical checks, not only symbolic inspection.
- `ADF`: estimates deterministic terms jointly inside each common-sample regression; its penalty projects the level on those trimmed-sample terms. This differs from full-series OLS detrending followed by a regression without deterministic terms. Consequently `ADF(method="maic", trend="ct")` should not be presented as an alias for PQ-selected DFGLS.

Released ARCH 8.0.0 has neither MAIC nor MBIC. Passing an externally selected integer `lags` to its public fixed-lag backend is available, but that alone is not a released MIC feature.

## 6. What bootUR 1.0.5 implements

With `criterion_scale=FALSE` and `two_step=TRUE`, bootUR OLS-detrends all N levels for selection. It estimates the candidate regression on zero-based rows K+1,…,N−1, equivalent to the n rows specified above. Residual variance uses SSR/n. Its `maic_cpp` and `mbic_cpp` implement the same criterion forms with multipliers 2 and log(n), respectively. Selection is OLS-based even when a GLS final regression is requested by the lower-level test routines.

**Penalty alignment discrepancy:** `adf_selectlags_cpp` sets

```text
ylag = de_trend(...).subvec(pmax + 1, n - 2)
```

Armadillo endpoints are inclusive. This gives zero-based levels K+1,…,N−2, whereas the regression's lagged-level vector uses K,…,N−2. Thus its penalty uses

$$
 S_{bootUR}=\sum_{j=K+2}^{N-1}v_j^2=S-v_{K+1}^2,
$$

one value fewer, while its SSR divisor remains n. It is not simply a different full-sample count. With rescaling off, the difference between the full-sum and bootUR criterion at a given k is

$$
 MIC_{full}(k)-MIC_{bootUR}(k)
 =\frac{C}{n}\frac{\widehat\beta_{0,k}^2 v_{K+1}^2}{\widehat\sigma_k^2}\geq0.
$$

It depends on k, so even selected orders can differ. Numerical validation checks this source reading against the installed binary; it does not patch bootUR.

**Additional output-label issue:** `adf_onestep_tests_parest_all_units_cpp` reports `par.n_elem - 1` as the selected lag count. Because this parameter vector includes the deterministic coefficients, the public one-step result reports actual k plus 1 for an intercept, or plus 2 for intercept/trend. The validation harness preserves this raw field as `reported_lags` and separately labels the inferred augmentation order; it does not modify the package.

With `two_step=FALSE`, deterministic terms are fitted jointly and the penalty instead uses raw levels, again starting one observation late. ARCH's one-step ADF penalty residualizes those levels on deterministic terms. Those two routines are not fully matched when deterministic terms are present.

`criterion_scale=TRUE` is bootUR's default: it rescales differences by a nonparametric volatility estimate and cumulatively sums them before selection. ARCH has no matching switch. Disable this for the main comparison and report rescaling separately. Choosing MAIC by name while leaving this default unchanged does not match ARCH.

## 7. Adoption table and remaining limits

| Candidate | MAIC | MBIC | Selection preprocessing | Adopt? |
|---|---|---|---|---|
| ARCH 8.0.0 | No | No | Standard AIC/BIC/t-stat | Reuse fixed-lag tests |
| Pinned ARCH development DFGLS | Yes | No | Full-series OLS, complete matched level sum | MAIC reference; no moving dependency |
| bootUR 1.0.5 two-step, unscaled | Yes | Yes | Full-series OLS, penalty omits first matched level | Comparator with documented discrepancy |
| Recommended later selector | Equation above | Equation above with log(n) | Full-series OLS, full common sample | Adopt this explicit contract after validation |
| Original-style GLS selection | Equation above | Equation above | GLS series in selection | Retain only as a clearly named research comparison |

This decision concerns lag selection, not implementation of the complete Ng–Perron M-test family, bootstrap inference, or bounded-series calibration. Stage 2's numerical study checks a small reproducible collection and diagnoses differences; it cannot certify original-paper Monte Carlo claims or establish performance for all data. Full Perron–Qu article review remains a bibliographic limitation; the central OLS substitution is verified from primary abstract and implementation source.

## Source digests

SHA-256 of downloaded research inputs (kept outside the repository):

| Input | SHA-256 |
|---|---|
| Ng–Perron published PDF | `442a6690d32c88b3b2ff28d0f708ba049a8c7bc958b9e7f2b34295dd0ff29787` |
| bootUR 1.0.5 CRAN source archive | `34ba4c50648b15ec98d68552efcdd520a5b458341bcf30854e6fb032f41e4a4a` |
| ARCH pinned GitHub source archive | `4b95b6e724abb045fcbfcd0cbfe7f7bfb297b0665d62b33ad81073b3d137e23f` |
| ARCH development `unitroot.py` | `b6be3547403ab880d425c780e1358c86eb9073a7b903ad8975650b45ad445fc9` |
