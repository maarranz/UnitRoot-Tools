# Deterministic AIC/BIC Selection — Methodology and API Design

28 October 2026 · UnitRoot-Tools (Python) · historical approved design preceding implementation

**Status: historical approved design. The experimental feature is now implemented.**

This document preserves the methodological contract approved before implementation.
Prospective instructions and validation plans below describe that design stage;
see [the current implementation, API and interpretation](Deterministic_Selection.md).

## 1. Recommended scope

Introduce an experimental comparison of ADF regressions with deterministic specifications `n`, `c` and `ct`. Estimate all three with a user-supplied common lag order and identical observations; choose the global minimum of AIC or BIC; return the full comparison and the selected model's existing fixed-specification ADF result.

Recommended initial interface:

```python
select_deterministic(y, *, lags: int, criterion: str = "bic")
```

Require `lags`, including an explicit `lags=0` when no augmentation is wanted. Defer internal automatic lag selection, `max_lags`, `method`, joint lag/specification searches, candidate-set options and model averaging. BIC is the proposed default for a parsimonious starting policy, not a claim of superior unit-root size or power. Compute and preserve both criteria for every candidate, even though only the requested criterion determines the winner.

The feature complements Enders by making a global penalized-likelihood comparison instead of following a sequence of unit-root and deterministic pretests. It is an experimental application of established criteria, not a claim that AIC/BIC or post-selection unit-root inference is a new established test.

## 2. Starting checkpoint and evidence

Read-only inspection confirmed branch `main`, HEAD `8c4f1a17cf12207dba3c051be6f1de7943df895c`, version `0.1.0rc1`. Existing untracked documents/evidence were preserved:

```text
?? docs/Enders_Mac_M1_Pro_Validation_Evidence_2026-10-28.zip
?? docs/Enders_Mac_M1_Pro_Validation_Report_2026-10-28.md
?? docs/Theoretical_Audit_Sequential_Strategies.md
```

The design follows §9.1 of the completed research assessment in the separately retained theoretical audit, with the initial lag policy narrowed to a mandatory supplied integer. The earlier separately retained Mac validation report records 1,287 full-suite and 123 focused Enders passes for source, fresh wheel and sdist installations on Python 3.12.16/macOS arm64. Those are prior validation results; no tests were rerun for this design assignment.

Repository sources inspected include `docs/Enders.md`, `docs/MIC_Methodology.md`, `docs/User_Guide.md`, `selection.py`, `_adapters.py`, `_validation.py`, `results.py`, and the public exports. The validated ADF and Enders paths are the preservation baseline.

## 3. Econometric specification and likelihood

For each deterministic choice fit

\[
\Delta y_t=a+\gamma t+b y_{t-1}
 +\sum_{j=1}^{p}\beta_j\Delta y_{t-j}+\varepsilon_t.
\]

| Specification | Restrictions on regression | Number of deterministic coefficients |
|---|---|---:|
| `n` | a=0 and gamma=0 | 0 |
| `c` | gamma=0 | 1 |
| `ct` | a and gamma unrestricted | 2 |

Estimate b and all included coefficients freely in every candidate. Do not impose b=0, constrain the fitted root to be stationary, demean/detrend the response separately, or choose a specification according to which ADF test rejects. This is a comparison of unrestricted estimating equations, not an IC comparison between separately fitted unit-root and stationary models.

Under independent homoskedastic Gaussian innovations, use the **conditional innovation likelihood**, conditioning on the first p+1 supplied levels needed to initialize these regressions. Ignore the initial-observation density consistently in all candidates. Do not claim to evaluate the exact unconditional likelihood of a stationary AR model. Lagged dependent regressors do not justify treating the entire design as fixed and strictly exogenous; the likelihood is the product of the sequential innovation densities given the past.

For common effective sample size n, q_d regression coefficients and residual sum of squares RSS_d, the likelihood is

\[
\ell_d(\theta_d,\sigma_d^2)
=-\frac n2\log(2\pi)-\frac n2\log\sigma_d^2
 -\frac{RSS_d(\theta_d)}{2\sigma_d^2}.
\]

OLS minimizes RSS and maximizes this conditional Gaussian likelihood. The variance MLE is

\[
\widehat\sigma_d^2=RSS_d/n,
\qquad
\widehat\ell_d=-\frac n2\{\log(2\pi)+1+\log(RSS_d/n)\}.
\]

Use `RSS/n` for likelihood, not the residual-df variance used for conventional OLS standard errors. The Gaussian working score remains computable outside the iid Gaussian model, but neither an appropriate likelihood nor good selection/inference properties follow automatically with residual serial correlation, heteroskedasticity, breaks or misspecified dynamics. These assumptions belong in the documentation.

### 3.1 Interpretation under null and alternative

The ADF null is b=0 against b<0, subject to the other dynamics satisfying the conditions needed for unit-root/stationary interpretation. A negative b alone is insufficient to ensure stability for an arbitrary augmented model.

- Under b=0 with stable difference dynamics, c allows drift in levels; its long-run increment mean is `a/(1−sum(beta_j))` when this denominator is nonzero. A genuine nonzero gamma in ct can generate quadratic mean growth in levels.
- Under a stable alternative, n corresponds to a zero deterministic equilibrium mean; c permits a nonzero constant mean; ct permits a linear deterministic mean. A difference-equation trend is not automatically a quadratic trend in stationary levels.
- Selecting c does not prove drift; selecting ct does not prove a genuine nonzero difference-equation trend. The selector does not decide integration order or the truth of deterministic coefficients.
- If prior economics excludes the no-constant model, this initial three-candidate selector may be inappropriate. Use a prespecified ADF specification; do not invent an undocumented candidate exclusion or choose n because it has a smaller p-value.

Standard regular-model interpretations of AIC/BIC do not establish uniform properties at a unit root, near-zero deterministic coefficients or after model/lag selection. No consistency, optimality or nominal-size guarantee is proposed here for the complete algorithm.

## 4. Exact information criteria and parameter count

Let d denote the deterministic dimension 0, 1 or 2. With a full-rank design,

\[
q_d=1+p+d,\qquad k_d=q_d+1=p+d+2.
\]

The first term in q_d is the lagged-level coefficient; the final +1 in k_d is the estimated innovation variance. Included intercept/trend coefficients are counted even when their estimates happen to be zero. Conditioned initial values and the supplied integer p are not additional estimated parameters. Each of the p augmentation coefficients is counted.

| Candidate | Regression coefficients q_d | Likelihood parameters k_d | Residual df |
|---|---:|---:|---:|
| n | p+1 | p+2 | n−p−1 |
| c | p+2 | p+3 | n−p−2 |
| ct | p+3 | p+4 | n−p−3 |

Define the stored comparison scores, using natural logarithms, as

\[
\mathrm{AICscore}_d=n\log(RSS_d/n)+2k_d,
\]
\[
\mathrm{BICscore}_d=n\log(RSS_d/n)+\log(n)k_d.
\]

These are **total**, not per-observation, scores. They omit the common Gaussian constant

\[
C_n=n\{1+\log(2\pi)\}.
\]

The full conventional values are `−2*log_likelihood + 2*k` and `−2*log_likelihood + log(n)*k`, respectively; add C_n to our scores to recover them. Store the full maximized log likelihood as well, so the convention is independently auditable. Scores may be negative; only differences on the matched comparison matter. Name the fields `aic_score` and `bic_score` and expose the normalization explicitly.

Proposed normalization label: `total_omit_common_gaussian_constant_include_variance`. No AICc, REML, per-observation score, MIC correction or estimated-scale omission is supported initially.

For nested candidates adding h coefficients, the larger model wins precisely when its gain `n*log(RSS_small/RSS_large)` exceeds `2*h` for AIC or `log(n)*h` for BIC; equality invokes the specified tie rule. This explains why a lower RSS alone is insufficient.

### 4.1 Independent convention checks

R's official AIC documentation specifies `−2*logLik + penalty*n_parameters` and comparison of maximized likelihoods on the same data. Its `logLik.lm(..., REML=FALSE)` implementation includes the full Gaussian constant and sets parameter df to regression rank +1. This agrees with the derivation above. [R AIC documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/AIC.html), [R logLik documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/logLik.html), [R implementation, logLik.lm](https://raw.githubusercontent.com/wch/r-source/trunk/src/library/stats/R/logLik.R).

Statsmodels ordinarily omits the scale parameter from OLS information-criterion counts. Consequently its default `.aic` and `.bic` are not our absolute scores. For a matched full-rank fit, our full AIC equals `.aic + 2`; our full BIC equals `.bic + log(n)`. Where supported, `info_criteria(..., dk_params=1)` includes scale. After subtracting C_n these are comparison-score checks; simply asserting equality to `.aic` would test the wrong convention. [Statsmodels OLS AIC](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.OLSResults.aic.html), [information-criterion parameter correction](https://www.statsmodels.org/stable/generated/statsmodels.regression.linear_model.RegressionResults.info_criteria.html).

These are references for definitions and future independent comparisons, not new runtime dependencies or claims of a completed numerical check.

## 5. Common sample, lag policy and numerical failures

Given N supplied levels in zero-based indexing and integer p:

- Response indices: `p+1, ..., N−1`, represented by `[p+1:N]`.
- Effective n: `N−1−p`.
- Common response: `diff(y)[p:]`.
- Common lagged level: `y[p:-1]`.
- Common lagged difference j, for j=1,...,p: `diff(y)[p−j : N−1−j]`.
- In ct, use the existing ARCH estimation-sample trend index `1,...,n`. Keep the time origin convention in the report; do not reset it independently between computations of scores and ADF results.

Call the existing fixed-lag `adf(y, trend=specification, lags=p)` for each candidate. Use its fitted raw-response RSS, regression design, t statistic and fixed-model calibration. Validate equality of response values, observation count, boundaries, lagged-level/difference columns and p across all candidates, rather than checking lengths only. The design matrices differ only by deterministic columns. Compute scores and expose inference from those same fitted objects; do not refit the winner on another sample.

The largest design is ct with p+3 coefficients. Positive residual df requires `N >= 2*p+5`, in addition to full rank, nonzero residual variance and existing adapter checks. This is an algebraic feasibility boundary, not a recommended minimum sample for reliable inference. Do not inherit Enders's n≥25 phi-table restriction: this procedure uses no phi table. Do not silently reduce p to make a comparison feasible.

Apply current real, finite, one-dimensional input validation and copy semantics. Reject invalid p, including booleans, negatives, None and nonintegral values. Follow existing case-insensitive label validation for AIC/BIC. Missing observations are not dropped or imputed. Equally spaced observations remain a researcher responsibility: the current array validator does not verify calendar spacing from a pandas index. Preserve the optional series name, and distinguish numerical row boundaries from dates.

All three candidates must be valid. If any candidate fails existing validation, has rank deficiency, nonpositive/nonfinite RSS, insufficient residual df or nonfinite likelihood/score, fail the comparison with the offending specification identified; do not silently choose among survivors. Use existing exception categories and preserve existing conditioning warnings. Do not suppress a failed `adf()` inference calibration merely to return an IC winner. No partial-success result is proposed initially.

### 5.1 Why mandatory lags

The provisional signature's `lags=None`, `max_lags` and `method="maic"` imply automatic lag selection, which adds another selection step, a detrending convention and potentially a different selection sample. Selecting p independently in each candidate also defeats the intended fixed-p experiment.

Defer those arguments completely. A caller may supply an integer obtained elsewhere, including from existing `select_lags()`, but the result must record `lag_policy="supplied_common_integer"`, not assert that p was selected independently of the data. No internal lag selection occurred; external lag-selection provenance is not inferred or certified. Documentation must state that upstream selection adds further uncertainty. A future automatic mode would need a separate approved convention and result provenance.

The existing `select_lags()` MAIC/MBIC/AIC/BIC scores concern lag comparison on separately detrended series, with omitted common quantities appropriate to that task. None of them, including its ordinary AIC/BIC fields, should be compared across deterministic specifications for this feature.

### 5.2 Tie-breaking and precision

Evaluate and display candidates in fixed order n, c, ct. Choose by the lexicographic pair `(unrounded_selected_criterion_score, deterministic_dimension)`. Exact score ties therefore favor n over c over ct. No near-tie tolerance, rounding before selection, or backend-default tie policy applies.

Report each candidate's difference from the minimum for the requested criterion and the gap to the runner-up. A zero gap records an exact tie. A small nonzero gap warns the reader of numerical sensitivity without changing the rule. Near-boundary winners can differ across platforms because of rounding; a deterministic tie rule does not establish universal cross-platform agreement for almost-ties. Do not interpret score gaps as model probabilities or unit-root p-values.

## 6. Public API and inference wording

### 6.1 Proposed call contract

```python
select_deterministic(y, *, lags: int, criterion: str = "bic")
```

Return a `DeterministicSelectionResult`. There is no alpha argument because significance levels do not determine the IC selection. Unsupported arguments, omitted lags and `lags=None` fail clearly; no compatibility shim for the provisional automatic signature is needed because it has not been released. No `reject` method, global p-value, or automatic 5% unit-root conclusion is added to the selection result.

### 6.2 The essential inference distinction

Use this reusable documentation wording:

> The displayed ADF p-value uses the calibration for the selected specification treated as fixed. It does not account for choosing that specification from these data, or for any upstream lag selection. It is neither a calibrated p-value for the complete selection-and-testing algorithm nor a selective p-value conditional on the event that this model won. The procedure has no established overall nominal size.

“Conditional on the selected specification” is acceptable shorthand only with this explanation. It must not suggest that the sampling distribution conditional on selection has been derived. Even the fixed-specification interpretation retains that ADF calibration's maintained assumptions and approximation limits.

Selecting an intercept or trend by IC never licenses replacing the nonstandard ADF calibration with ordinary Student t. Use the selected candidate's unchanged ARCH 8.0.0 tau calibration; no Enders phi or conditional conventional-t branches enter this selector.

### 6.3 Summary behavior

The first line of `summary()` should read:

`Experimental deterministic selection by BIC — supplied common lags: p`

Then show common sample boundaries, the three candidate score rows, parameter counts, selected specification and runner-up gap. End with a section labelled **Fixed-specification ADF calibration for the selected model**, showing statistic, approximate p-value, critical values and provenance. Place the inference limitation immediately beside these numbers, not only in a footnote. Do not print “BIC rejects the unit root” or an unconditional “decision at 5%”.

Do not reuse the existing `UnitRootResult.summary()` wholesale inside the parent summary: it prints a nominal 5% fixed-model decision and could obscure the distinction. Existing result behavior remains unchanged for users who explicitly access the nested object.

## 7. Result structure and illustrative usage

Use new, small frozen scalar-record types alongside existing types, not subclasses that redefine their meanings. This is a data contract, not an implementation.

### 7.1 DeterministicCandidate

| Field | Meaning |
|---|---|
| `specification` | n, c or ct |
| `deterministic_dimension` | 0, 1 or 2 |
| `n_coefficients` | q_d, including the lagged-level coefficient |
| `n_parameters` | k_d=q_d+1, including variance |
| `df_resid` | n−q_d |
| `ssr` | Raw-response residual sum of squares |
| `variance_mle` | SSR/n |
| `log_likelihood` | Maximized conditional Gaussian log likelihood, with its common constant |
| `aic_score`, `bic_score` | The reduced total scores defined in §4 |
| `adf_result` | Existing, unchanged `UnitRootResult` for this fit |

### 7.2 DeterministicSelectionResult

| Field | Meaning |
|---|---|
| `criterion` | Normalized `aic` or `bic` |
| `selected_specification` | The winning n/c/ct code |
| `lags` | The supplied common augmentation order |
| `lag_policy` | `supplied_common_integer` |
| `input_nobs`, `nobs` | N and the common n |
| `sample_start`, `sample_stop` | Zero-based, stop-exclusive response boundaries p+1 and N |
| `candidates` | Tuple of all three `DeterministicCandidate` objects, in n/c/ct order |
| `criterion_normalization` | Explicit label from §4 |
| `tie_break` | `exact_score_then_smaller_deterministic_dimension` |
| `inference_status` | `fixed_specification_calibration_after_selection` |
| `experimental` | True |
| `versions`, `notes`, `input_name` | Provenance using existing package conventions; no inferred external lag-selection provenance |

Convenience properties: `selected_candidate`, `criterion_value`, `runner_up_gap`, and `selected_fixed_spec_adf`. The last returns the **same** `UnitRootResult` as `selected_candidate.adf_result`; it does not construct another fit. Calculate score gaps from stored candidates. Preserve one source of truth for the winner and avoid duplicated flattened p-values/decisions at the top level. Frozen records do not make the nested Statsmodels regression objects deeply immutable, as already documented for current results.

Illustrative usage from the approved design; see the current implementation page for a self-contained, validated example:

```python
from unitroot_tools import select_deterministic

comparison = select_deterministic(y, lags=4, criterion="bic")
print(comparison.summary())
print(comparison.selected_specification)

for candidate in comparison.candidates:
    print(candidate.specification, candidate.n_parameters,
          candidate.aic_score, candidate.bic_score)

# This is the existing fixed-model calibration, without selection adjustment.
fixed_spec = comparison.selected_fixed_spec_adf
print("Selected-model fixed-specification ADF statistic:", fixed_spec.statistic)
print("Fixed-specification ADF p-value reported after selection:", fixed_spec.pvalue)
```

An explicitly requested `fixed_spec.reject(alpha)` continues to use the existing left-tail critical-value rule and supported levels, including its strict equality convention. That remains a fixed-model calculation, not a selector-level decision. No post-selection p-value is synthesized from the three component results or from the IC scores.

## 8. Architecture and compatibility

The existing public fixed-lag `adf()` adapter provides all three candidate regressions, calibration, fitted objects and sample metadata. A small new selector module can compose these results without changing `_enders.py`, the lag selector or current result classes. The approved implementation adds the new function and result types to public exports as an additive API change.

Do not reuse the lag-oriented `CandidateDiagnostics` or `LagSelectionResult`: their axes, normalization and detrending semantics differ. In particular, the existing lag selector's candidate df/RSS refer to its selection regression, not necessarily the final raw-level ADF fit. Name the new comparison explicitly rather than overloading those types.

Potential conflicts to resolve during implementation review are: Statsmodels' omitted-scale criterion convention; accidentally leaving `lags=None` and triggering backend selection; different response slices; silently skipped invalid fits; and a parent summary that inherits a misleading unconditional decision label. The contract above resolves each without a broad refactor. No general strategy framework, new numerical solver, new dependencies, regression-storage switch or bootstrap infrastructure is needed.

## 9. Historical independent software-validation plan

1. **Regression reconstruction.** On fixed, reproducible nondegenerate level vectors and the published GDP fixture, independently build the common response/design matrices for p=0 and positive p. Use NumPy QR or a separate SciPy least-squares calculation to verify coefficients, residual arrays, SSR, residual df and full-rank column counts. Verify every candidate against the existing fixed-lag ADF call. Do not use the selector itself to generate expected regressions.
2. **Likelihood and criterion arithmetic.** Independently sum Gaussian log densities at fitted means using variance SSR/n and verify the concentrated likelihood. Check reduced and full scores, all three k values, both penalties, and the common-constant conversion. Include tests that would fail if SSR/df_resid, omission of the level coefficient/variance, log(N) instead of log(n), or a per-observation score were used. Compare with independently specified OLS criteria using the documented variance-parameter correction; optional R reference calculations are not a package dependency.
3. **Sample and provenance invariants.** Check exact row values and boundaries, identical lag columns across candidates, trend origin and no winner refit. Verify p is supplied to each call and no lag-selector call occurs. Ensure the original input is not changed. Show that multiplying the entire series by nonzero s shifts every score by `2*n*log(abs(s))` without changing the ranking, within numerical tolerance. Do not assert translation invariance: including n makes the level origin substantively relevant.
4. **Selection and ties.** Verify the global minimum, canonical candidate ordering, both score columns, minimum-relative differences and runner-up gap. Use constructed equal **floating-point score values** to exercise the exact tie comparator; algebraic real-number ties alone may round differently. Check n/c/ct ties, candidate-order independence and near-ties without a tolerance. Add real regression fixtures where AIC and BIC disagree and where ct wins; lock their independently checked margins, not just labels.
5. **Invalid data/fits.** Test the supplied-lag contract, invalid labels, nonfinite/nonreal/multidimensional input, structural sample boundary, rank deficiency, degenerate RSS and an invalid single candidate. Verify all-or-error behavior, candidate identification and inherited warnings. No silent candidate removal, sample change or score normalization fallback is allowed.
6. **Inference/result consistency.** Check each nested ADF statistic, approximate p-value, critical values, calibration and versions against the unchanged fixed-model adapter. Verify identity of the selected nested object, immutable scalar containers and absence of selector-level `.pvalue`/`.reject`. Test the summary's selection warning and its lack of an unconditional 5% conclusion. Preserve the documented limits of nested regression mutability.
7. **Regression and distribution preservation.** Rerun the complete existing suite and independent Enders script after any approved implementation. Compare Original/Conservative histories, selected lags, samples and GDP decisions with the validated checkpoint. Build fresh artifacts in a committed-source copy outside the repository; inspect runtime bytes and package contents. Run full/focused checks on source, a wheel installation, and a fresh sdist installation from directories where checkout source cannot shadow installed code. Record environment/BLAS versions, failures/skips/warnings and numerical differences. Reuse the Mac-validation preservation safeguards; existing saved repository distributions predate Enders and must not be treated as current feature artifacts.

Proposed initial numerical tolerances: retain the existing fixture-specific regression checks where applicable; use a separately stated absolute 1e-8 ceiling for independent scalar likelihood/score comparisons on appropriately scaled fixtures. Exact tie tests compare deliberately equal scores, not values within that ceiling. Calibrated probability bounds, decision labels and sample indices must not be validated by loose numerical tolerances. Review tolerances against actual fixture scales before accepting a new validation report; do not loosen them to conceal a disagreement.

### 9.1 Analytic scoring example for AIC/BIC disagreement

This is a designed **arithmetic fixture**, not results from a fitted series or a simulation. Set n=100, p=1 and candidate SSRs to `100`, `100*exp(−0.03)` and `100*exp(−0.035)` for n/c/ct. Then:

| Candidate | k | n log(SSR/n) | AIC score | BIC score |
|---|---:|---:|---:|---:|
| n | 3 | 0 | 6 | 3 log(100) |
| c | 4 | −3 | 5 | 4 log(100)−3 |
| ct | 5 | −3.5 | 6.5 | 5 log(100)−3.5 |

AIC selects c while BIC selects n. This separates scoring logic from future regression fixtures and gives a non-boundary expected result. It does not establish that these SSRs arise from a particular supplied ADF series; real-series disagreement fixtures still require independent construction and verification.

## 10. Later statistical-validation design — no simulations now

After specification approval, implementation and independent numerical validation, design a separate Monte Carlo study of the **whole selection-and-test algorithm**. Begin with correctly specified iid Gaussian fixed-p models, then introduce controlled stable AR/MA short-run dynamics and fixed-order misspecification. Separate this from any future automatic-lag extension.

Use null DGPs with b=0 and zero, small/local and fixed nonzero drift; include genuine difference-equation trends only as an explicitly labelled expanded null permitting quadratic mean growth. For alternatives, use stationary and local-to-unity processes around zero/nonzero means and deterministic linear trends, with all remaining roots stable. Vary sample size, innovation scale relative to deterministic terms, initial conditions and short-run parameters. Define local sequences and initial-value rules before running; do not equate an estimating-regression coefficient with a level-trend coefficient across DGPs.

Compare AIC and BIC selectors against prespecified/oracle ADF specifications, always-retained c/ct benchmarks, and the existing distinct Original/default Conservative Enders procedures with matched lag/sample conventions. Do not duplicate Original by treating permitted Conservative as a separate research algorithm. Include a calibrated non-bootstrap union benchmark only after its precise source definition has been verified.

Report deterministic-selection frequencies, unconditional unit-root rejection rates under each declared null, raw and size-adjusted power, score-margin distributions, numerical-failure rates and Monte Carlo uncertainty. Conditional-on-selection rejection frequencies may be descriptive diagnostics, not evidence that the displayed component p-values are valid selective p-values. Fix seeds, repetitions, calibration/training versus evaluation separation, and interval conventions in a separately approved protocol. Report uncertainty especially where selection events are rare. This study would assess empirical performance for its specified DGPs, not prove general validity.

## 11. Decisions for approval before coding

| Decision | Recommended initial choice | Alternative requiring more scope |
|---|---|---|
| Lag handling | Mandatory `lags` integer; no automatic arguments | Internal common-ct lag selection with explicit additional provenance |
| Public return | Complete `DeterministicSelectionResult` with selected-specification convenience properties | Returning only a code would lose the audit trail |
| Default criterion | BIC, with explicit AIC option; both scores always stored | Mandatory criterion if no default is desired |
| Candidate set | Fixed n/c/ct; fail if any fit is invalid | Configurable candidate subsets or partial-success comparisons |
| Criterion convention | Total reduced scores, variance counted, full ML log likelihood stored | Full AIC/BIC values could be stored instead if field names/conversion are revised consistently |
| Ties | Exact unrounded score, then smaller deterministic dimension | Near-tie tolerance would define a different selection policy |
| Inference exposure | Nested unchanged `UnitRootResult` under `selected_fixed_spec_adf`; prominent fixed-calibration warning; no selector p-value/rejection method | A separate inference wrapper would add types without solving post-selection calibration |
| Release status | Experimental model-selection feature, subject to independent tests; no overall-size claim | Inclusion as an established calibrated unit-root test is unsupported |

The methodological contract and API above were approved before implementation. The current implementation follows these choices; see the linked implementation page. Existing ADF/Enders behavior, dependencies and version remain the baseline. Monte Carlo and bootstrap research require a separate assignment; repository commits and pushes are governed by subsequent explicit authorization.
