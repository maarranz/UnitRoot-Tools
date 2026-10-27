# UnitRoot-Tools: implementation specification and final validation audit

Stage 2C audit with Stage 3A/3B/3C instructor-approved amendments, 27 October 2026. Stages 3A and 3B completed the selector and ADF/DF-GLS integration. Stage 3C completes PP/KPSS within the approved scope. Historical stage-specific descriptions below record the scope at that stage. Work remains in the existing repository; no commit or push is authorized.

## 1. Audit verdict and evidence

**Stage 2 passes its stated numerical validation scope.** There is no unexplained numerical discrepancy blocking implementation of the architecture confirmed by the instructor. This is permission to recommend the next stage, not evidence that a production implementation already exists or that all future defaults and edge cases have been validated.

Read and inspected: `ARCH_Audit.md`, `MIC_Methodology.md`, `Numerical_Validation.md`, both `validation/cross_validate` scripts, and every file in `validation/results/` (six CSVs and five JSON files). The pre-existing `validation/results.zip` was also checked for CRC integrity and byte-for-byte agreement with its eleven result files; it was left unchanged.

The independent audit is reproducible with [audit_saved_results.py](../validation/audit_saved_results.py), and its output is [stage2c_audit.json](../validation/stage2c_audit.json). It does **not** import the Stage 2B harness or its MIC/PP reference functions. It reads the saved evidence, uses fresh QR-based least-squares calculations with explicit observation indices, and additionally tests public fixed-lag calls in released ARCH 8.0.0. Run with NumPy, SciPy and released ARCH installed:

```sh
python validation/audit_saved_results.py --report /tmp/unitroot-stage2c-audit.json
```

### Independently established results

| Check | Finding |
|---|---|
| Recorded comparisons | Exactly **1,188**, with differences and tolerance flags recalculated from numeric values |
| Expected-equality checks | Exactly **840**, all pass absolute tolerance **1e-8**; integer lag choices agree exactly |
| Summary consistency | Every group count, maximum difference and failure count agrees with `summary.json` |
| Original manifest | All four raw-data/result hashes and both Stage 2B script hashes match |
| Saved observations | 216 Python-release rows, 120 Python-development rows, 252 R rows; unique keys, finite reported numbers, valid probabilities and correct effective sample counts |
| Inputs | Five length-200 series and one length-64 diagnostic; all six regenerated exactly from the documented PCG64 seeds and recurrences |
| MIC candidate evidence | All **1,620 rows** independently reconstructed: response/lag alignment, deterministic treatment, variance denominator, level sum, tau, MAIC and MBIC |
| MIC output alignment | **264** unscaled saved selector outputs agree with independently selected orders, including both memory routes and corrected/raw bootUR lag fields |
| Public released-ARCH composition | **120** additional checks: choose MAIC/MBIC from aligned PQ curves, pass p through public `ADF(..., lags=p)` / `DFGLS(..., lags=p)`, verify final statistics and N−1−p sample counts; maximum statistic error **7.11e-14** |

The 840 are scalar comparisons, with shared data and some repeated diagnostics; they are not 840 independent econometric experiments. The additional 120 checks are reported separately, not added to the historical 1,188/840 totals.

Fresh QR reconstruction versus stored candidate fields had maximum errors: coefficient 1.88e-14, variance 4.45e-15, level sum 6.01e-11, tau 8.82e-11, MAIC 9.44e-13, MBIC 2.47e-12. Fresh fixed-test reconstruction had maximum statistic errors: ADF 3.09e-12, KPSS 1.73e-14, DFGLS 4.71e-14, PP tau 7.20e-14 and PP rho 2.01e-12. These all satisfy the declared absolute 1e-8 tolerance. Original between-package differences remain those in the saved CSVs, rather than these fresh-reconstruction rounding errors.

### Documented differences verified in the saved data

- ADF: 48 statistic/p-value/critical-value comparisons per field; all agree at tolerance. KPSS: 24 statistics agree, but p-values and calibration triples differ as documented. DF-GLS: 24 statistics agree while all three critical-value fields differ.
- PP: all 24 tau and all 24 rho package pairs differ; source-specific finite-sample corrections reconcile both packages. This is not a reason to replace released ARCH's public statistic with an R formula.
- Matched two-step, unscaled MAIC: ARCH development and bootUR agree in 20/20 base cases and 23/24 including the diagnostic. One-step ADF differs in 7/36 after accounting for bootUR's lag-label issue. Rescaling changes bootUR choices in 3/36. Aligned MBIC and bootUR choose the same order in 36/36, without equal criterion values.
- `mic_boundary`, seed 457, intercept, K=4: n=59; full level sum 42.200588309 versus bootUR-style 40.265549277. At p=2/3, aligned MAIC is 0.552905055/0.552592109 (select **3**); bootUR-style is 0.508145123/0.509366759 (select **2**). The difference at **each** candidate is exactly the omitted first matched squared-level contribution, multiplied by its coefficient/variance factor. Both MAIC and MBIC obey that identity.

### Corrections and limits

No numerical correction to the original CSVs or calculation scripts is required for the saved study. Historical scripts and results remain unchanged so their manifest remains meaningful. Two clarifications are added to `Numerical_Validation.md`:

1. Development metadata's `all_distributions["arch"]` says 8.0.0 because dictionary enumeration collapses two installed distributions. `arch_runtime_version`, the targeted `packages` entry and the development source SHA identify the imported code correctly. That dictionary is not an environment lock. Before a future regenerated evidence set, record distributions as a list with origins, or use a separate development environment.
2. Future input readers must reject interior gaps rather than copying the research harness's blanket blank/NA filtering. The present saved series have trailing padding only, verified independently, so this does not invalidate them. Future harness runs should hash derived CSVs/metadata too, validate candidate matrix rank, and write a fresh output directory. The 2C report now hashes all eleven existing files.

The complete Perron–Qu article was not obtained in Stage 2A. The explicit operational contract below is supported by the primary abstract, original Ng–Perron formulas and inspected implementations; it must not be advertised as a complete article replication. Automatic K and bandwidth choices, degenerate inputs, alternative Python environments and production error handling remain implementation acceptance tests. MBIC external evidence checks selected orders, not a native ARCH MBIC implementation, which does not exist. These are limits and release gates, not unexplained failures in the existing evidence.

## 2. Confirmed production architecture

Use **released `arch==8.0.0`**. Own the input checks, MIC selector, result metadata and teaching-oriented display. Delegate final estimation to public `arch.unitroot.ADF`, `DFGLS`, `PhillipsPerron` and `KPSS`. Always supply an explicit integer augmentation order to ADF/DFGLS after selection, bypassing all backend automatic lag-selection branches.

No production import from `arch.unitroot.unitroot`, `_df_select_lags`, `_ic_best`, internal calibration arrays, or development snapshots. No subclassing, monkey-patching, copied ARCH detrending routines or upstream modifications. The historical development harness remains a research artifact; normal package use and normal regression tests must not require that snapshot or R.

**Methodological distinction:** applying a PQ-selected p to the ordinary ADF regression is a valid explicitly named lag-selection convention. It does not reproduce ARCH development's one-step ADF MAIC selector, nor bootUR's final two-step OLS statistic. DF-GLS selects on OLS-detrended data and tests on GLS-detrended data; pass the original series and selected p to ARCH, rather than passing an already detrended series and detrending twice.

PP/KPSS use long-run-variance bandwidths, not MIC augmentation orders. Use ARCH calibration consistently with ARCH's final statistics. There is no complete Ng–Perron M-test family, structural-break test, cointegration, bootstrap, volatility-rescaling or bounded-series mode in this implementation stage.

## 3. Proposed small public API

Approved import namespace: **`unitroot_tools`**. The five functions are approved for the completed API; Stage 3A exported only `select_lags`, its result/diagnostic records, and the public error/warning types. Stage 3B adds `adf`, `dfgls` and `UnitRootResult`. Do not create adapter stubs for the deferred four functions. Five functions, two returned result types, and two documented error/warning types are sufficient. All options after y are keyword-only. `ArrayLike` below means the one-dimensional numeric input contract in Section 8.

```python
def select_lags(
    y: ArrayLike, *, trend: str = "c",
    max_lags: int | None = None, method: str = "maic"
) -> LagSelectionResult: ...

def adf(
    y: ArrayLike, *, trend: str = "c", lags: int | None = None,
    max_lags: int | None = None, method: str | None = None
) -> UnitRootResult: ...

def dfgls(
    y: ArrayLike, *, trend: str = "c", lags: int | None = None,
    max_lags: int | None = None, method: str | None = None
) -> UnitRootResult: ...

def pp(
    y: ArrayLike, *, trend: str = "c", bandwidth: int | None = None,
    test_type: str = "tau"
) -> UnitRootResult: ...

def kpss(
    y: ArrayLike, *, trend: str = "c", bandwidth: int | None = None
) -> UnitRootResult: ...
```

Resolved behavior:

- `adf(y)` / `dfgls(y)`: `lags=None, method=None` resolves to automatic **MAIC** with the declared K rule. `method="mbic"` changes the modified penalty multiplier; explicit `method="aic"` / `"bic"` use the distinct conventional penalties defined in Section 4.
- `lags=p` means fixed order. `method` and `max_lags` must then both be `None`; reject conflicting settings rather than ignoring them. `lags=0` is a real fixed-order request.
- `select_lags` supports `"maic"`, `"mbic"`, `"aic"` and `"bic"`, returning the full small candidate table. Normalize case once for method labels, trend labels and PP test type; reject unknown labels. AIC/BIC are explicit alternatives; MAIC remains the default. A t-stat selector is outside this first API.
- PP test types are `"tau"` and `"rho"`. There is no `lags` argument to PP/KPSS and no `bandwidth` argument to ADF/DFGLS.
- Materialize the statistic, calibration and regression (if available) before returning, so backend lazy errors occur inside the function call. No mutable specification setters or automatic differencing until significance is reached.

Teaching examples (proposed API only):

```python
adf(y, lags=4)                          # fixed augmentation order
adf(y, method="mbic", max_lags=12)      # PQ OLS selection, then ordinary ADF
dfgls(y, trend="ct", max_lags=12)       # PQ MAIC selection, then GLS test
pp(y, bandwidth=12, test_type="tau")    # fixed covariance bandwidth
kpss(y, trend="ct", bandwidth=12)       # stationarity null
```

### Deterministic specifications

| Function/mode | Supported terms | Interpretation |
|---|---|---|
| `select_lags`, automatic ADF | `n`, `c`, `ct` | None; intercept; intercept plus linear time trend in full-series OLS selection detrending |
| Fixed-order ADF | `n`, `c`, `ct`, `ctt` | `ctt` additionally includes quadratic time trend in the final ADF regression |
| DF-GLS, fixed/automatic | `c`, `ct` | Selection OLS trend matches final GLS deterministic specification |
| PP | `n`, `c`, `ct` | Deterministic terms in the level regression |
| KPSS | `c`, `ct` | Stationarity around an intercept or linear trend |

All automatic selection with `ctt` is deferred because Stage 2's MIC evidence does not cover it. `adf(..., trend="ctt", lags=p)` is supported; an automatic request with that trend raises an informative error. Do not infer that fixed-order ADF's larger deterministic menu applies to every selector/test.

## 4. Lag selection and sample contract

Let N be the supplied observation count, K the maximum candidate order, p the candidate order, and d the count of deterministic regressors (0/1/2 for `n`/`c`/`ct`). Use the **complete** series to estimate OLS deterministic coefficients. A normalized linear-time regressor t/N is permitted, because it spans exactly the same space as t. With `n`, leave y unchanged. Define resulting levels v and differences Δv.

For every p=0,…,K, estimate without additional deterministic terms:

$$
\Delta v_t=\beta_{0,p}v_{t-1}
          +\sum_{j=1}^{p}\beta_{j,p}\Delta v_{t-j}+e_{t,p},
\quad t=K+2,\ldots,N.
$$

Use the same n=N−1−K rows at every p, including p=0. The response is `diff(v)[K:]`, the complete matched level column is `v[K:-1]`, and the jth difference lag uses those same response indices shifted by j.

$$
S=\sum_{j=K+1}^{N-1}v_j^2,\quad
\widehat\sigma_p^2=\frac{\sum_{t=K+2}^{N}e_{t,p}^2}{n},\quad
\tau_p=\frac{\widehat\beta_{0,p}^{\,2}S}{\widehat\sigma_p^2}.
$$

$$
MAIC(p)=\log\widehat\sigma_p^2+\frac{2(\tau_p+p)}{n},
\qquad
MBIC(p)=\log\widehat\sigma_p^2+\frac{\log(n)(\tau_p+p)}{n}.
$$

### Conventional AIC/BIC amendment (approved before Stage 3A implementation)

All four methods use **the same complete-series OLS detrending**, common response rows t=K+2,…,N, unrestricted candidate regressions without further deterministic terms, variance SSR/n, and the same feasible candidate range 0,…,K. The conventional criteria are

$$
AIC(p)=\log\widehat\sigma_p^2+\frac{2p}{n},
\qquad
BIC(p)=\log\widehat\sigma_p^2+\frac{\log(n)p}{n}.
$$

They omit tau rather than estimating a different candidate regression. Therefore MAIC−AIC = 2·tau/n and MBIC−BIC = log(n)·tau/n **at each candidate**, with identical SSR, variance, level coefficient and level sum. They are not aliases and can select different orders. All four scores are always included in each diagnostic row, while `criterion_value`, selected p and runner-up gap refer only to the requested method.

Reported scores use normalization `per_observation_omit_common_constants`. For the no-deterministic conditional Gaussian candidate regression, its p+1 fitted coefficients include the lagged level. Dividing conventional full −2 log-likelihood AIC/BIC by n yields the above scores plus the common terms log(2π)+1 and, respectively, 2/n or log(n)/n. Counting an innovation-variance parameter would add another method-specific common term. These common terms do not affect argmin; they are deliberately not included in the reported scores. Deterministic coefficients estimated once in the preprocessing step likewise add only common counts if included.

“Conventional” refers to the unmodified AIC/BIC penalties **within this explicitly specified two-step OLS-preprocessed regression family**. This is not a promise to reproduce Statsmodels/ARCH automatic one-step ADF selection with deterministics fitted on the trimmed sample. It also differs from bootUR's volatility-rescaled defaults. Tests compare conventional scores to public Statsmodels OLS likelihood criteria on the matched preprocessed common-sample regressions, with the common terms removed; they do not use backend autolag or private selectors.

Use natural logarithms and unrestricted level coefficients. **Never** divide variance by residual degrees of freedom, replace n by N, use a p-dependent selection sample, omit the first level in S, use GLS selection detrending, or enable volatility rescaling in this contract. Choose the smallest actual floating-point criterion; exact ties choose the smallest p. A tolerance must not silently promote a near-tie to a tie. Record the runner-up criterion gap for teaching and diagnostics.

Final ADF/DFGLS estimation refits at chosen p on **N−1−p** observations. ADF uses original y and its requested deterministics; DFGLS receives original y and performs its own GLS detrending with ARCH's constants/initial-row convention. Backend covariance/standard-error normalization for the final statistic remains ARCH's, distinct from the MIC SSR/n divisor.

### Candidate cap and short-sample feasibility

Proposed automatic request is `K_rule = ceil(12*(N/100)**0.25)`. It is a reproducible convenience rule, not a data-optimal cap validated by Stage 2. For automatic defaults, cap to the smallest of the rule and the appropriate feasibility bounds; record uncapped and effective K and emit `SpecificationWarning` when a cap is reduced. An **explicit** infeasible K is an error and is never clipped or searched over a smaller range after a failure.

Necessary bounds (in addition to full rank and nonzero residual variance):

| Context | Bound |
|---|---|
| Standalone selector | `K <= floor((N-3)/2)` and `N >= max(3, d+2)` |
| Automatic ADF | Selector bound plus `K <= floor((N-3-d)/2)` and `N >= d+3` |
| Automatic DFGLS | Selector bound plus `K <= N-3-d` and `N >= d+3` |
| Fixed ADF | `N >= 3+d+p` and `N-1-p > 1+p+d` |
| Fixed DFGLS | `N >= 3+d+p` and `N-1-p > 1+p` |

The second-column constraints guarantee positive final residual degrees of freedom and honor released ARCH's minimum sample check. They do not guarantee rank for arbitrary data. Standalone selection has no knowledge of the eventual final test: ADF's automatic wrapper may therefore use a stricter K than a standalone selector. Report that distinction. K=0 remains valid when all conditions hold.

Use rank-revealing QR/SVD least squares and explicit residuals, not an inverse of X'X or subtraction of nearly equal quadratic forms for SSR. Column/time normalization must preserve estimates and criterion units after rescaling. Detect rank deficiency with a documented standard SVD threshold (machine epsilon × larger matrix dimension × largest singular value on the scaled design). Fail the entire requested search if any required candidate is infeasible; do not silently skip it. No low-memory alternative is needed initially.

## 5. Long-run-variance bandwidths

For PP, covariance residual count m=N−1; for KPSS, m=N. Fixed integer bandwidth b must satisfy **0 <= b < m**. Bartlett/Newey–West weights and statistic formulas remain ARCH's public implementation. Do not translate MIC-selected p into b.

Proposed defaults:

- PP `bandwidth=None`: resolve `b_rule=ceil(12*(N/100)**0.25)` and `b=min(b_rule, N−2)`, then pass integer `lags=b` to `PhillipsPerron`. Record and warn if capped. This bounded default is a transparent wrapper convention and may differ from unbounded released ARCH defaults in very short samples.
- KPSS `bandwidth=None`: request released ARCH's data-dependent Hobijn-style selector through public `KPSS(..., lags=None)`. Record the public resolved `.lags`. Validate `0 <= b < N` and finite results; if the backend produces an unacceptable bandwidth, fail with a request to specify b rather than silently recomputing a different statistic. Suppress only its exact historical lag-default DeprecationWarning, replacing it with the clear method label; preserve other warnings.

PP must have N >= d+3, and KPSS N >= d+2 under the wrapper's proposed residual-degree checks. Automatic bandwidths, b=0 and the upper valid bandwidth require implementation tests; the Stage 2 comparisons exercised only b=4 and 12. No generalized kernel selection or public legacy `-1` sentinel is proposed.

## 6. Results, regression access and teaching metadata

Use immutable result records with attributes, not backend-dependent tuple lengths. Read-only candidate arrays/mappings should not mutate subsequent summaries. A returned Statsmodels regression object remains a backend object; do not promise deep immutability or cross-version pickling of it.

`UnitRootResult` must provide:

- `test`, `statistic`, `pvalue`, `critical_values` keyed by nominal alpha **0.01/0.05/0.10**; `null_hypothesis`, `alternative_hypothesis`, `tail` (`left` for ADF/DFGLS/PP, `right` for KPSS).
- `trend`, `input_nobs`, final `nobs`, final response `sample_start` (zero-based) and `sample_stop` (exclusive). Starts are p+1 for ADF/DFGLS, 1 for PP, 0 for KPSS; stop is N.
- `lags` for ADF/DFGLS only, otherwise `None`; `bandwidth` for PP/KPSS only, otherwise `None`; `test_type` for PP only; `bandwidth_method` (`fixed`, `schwert_capped`, `hobijn`) and requested/resolved bandwidth details where applicable.
- `lag_selection: LagSelectionResult | None` (None for fixed p and PP/KPSS); `backend="arch"`, `backend_version="8.0.0"`, UnitRoot-Tools version, NumPy/SciPy/Statsmodels versions, calibration/method labels and explicit notes.
- `regression`: public `.regression` from ARCH ADF/DFGLS/PP, **None for KPSS**. Do not retrieve `_regression`/`_resids` or reconstruct a KPSS regression while calling it an ARCH result.
- `summary() -> str` and `reject(alpha=0.05) -> bool` as defined below. No generic serialized backend object or private-result escape hatch is needed.

`LagSelectionResult` must provide `method`, `selection_detrending="ols"` (or `"none"` for `n`), `trend`, `input_nobs`, requested/effective `max_lags`, uncapped rule value if used, `selected_lags`, `selection_nobs`, `sample_start=K+1`, `sample_stop=N`, `criterion_value`, runner-up gap (None with one candidate), and an immutable ordered candidate table. Each table row contains p, n, residual degrees of freedom, fitted level coefficient, SSR, variance, S, tau, the selected criterion value and all four AIC/BIC/MAIC/MBIC values for comparison. `CandidateDiagnostics` is a frozen scalar record, and `candidates` is a tuple of those records; `criterion_normalization` records the convention above. No backend-private values are required to produce it.

Regression labels must explain: DFGLS regression is on GLS-detrended data; PP regression is a **level** regression with HAC covariance. Its ordinary coefficient p-values do not test the PP unit-root hypothesis. ADF/DFGLS ordinary regression p-values similarly must not replace unit-root calibration. KPSS rejection is evidence against the stated stationarity null, not proof of a specific unit-root alternative.

## 7. Critical values, p-values and decisions

Read `stat`, `pvalue`, `critical_values` and `nobs` only from the instantiated public released-ARCH test object. Preserve its implemented sample-size convention and calibration; map its 1%/5%/10% labels to numeric alpha keys without changing numbers. Do not splice urca critical values into ARCH statistics or use Statsmodels KPSS p-values with ARCH KPSS results.

Record calibration provenance:

| Test | Result label |
|---|---|
| ADF / PP tau | ARCH 8.0.0 ADF-t MacKinnon p-value approximation and sample-size critical-value surfaces, with historical no-constant exception |
| PP rho | ARCH 8.0.0 ADF-z approximation/surfaces |
| DFGLS | ARCH 8.0.0 dedicated DFGLS simulation response surfaces, not ordinary intercept-ADF tables |
| KPSS | ARCH 8.0.0 simulated-quantile interpolation, not the Statsmodels four-point table |

`reject(alpha)` supports exactly 0.01, 0.05 and 0.10 initially. Compare the statistic directly with the returned critical value: **statistic < critical value** for left-tail tests, **statistic > critical value** for KPSS. Exact equality does not reject. Reject invalid alpha values; never interpolate extra critical values. P-values remain visible as approximate complementary information and are not used to override this decision rule near a calibration boundary. A future arbitrary-alpha mode would need a separate explicit contract.

Summaries display H0, deterministic terms, selection method, p or b with its meaning, final and selection sample sizes, statistic, approximate p-value and nominal critical-value decision. “Fail to reject” must not be printed as “accept H0.” No finite-sample exact-size or post-selection-adjusted inference claim is justified.

Expose `pvalue_status`/notes for approximation and boundary values: ARCH KPSS's supported table endpoints are 0.0001 and 1.0; flag equality as a table-boundary value rather than an exact tail probability. Likewise flag unit-root response-surface p-values at 0 or 1. These are documented version-specific presentation checks of public outputs, not imports of private calibration arrays. Unknown non-finite or out-of-range probabilities are errors.

## 8. Input validation and errors

Accept a nonempty one-dimensional real numeric list/tuple, NumPy array or pandas Series through numeric conversion; copy to a float64 working array. Reject scalars, multi-column/DataFrame and two-dimensional inputs (including N×1), booleans, complex/string/object data, NaN, infinity and overflow after conversion. A Series name may be retained for display; do not sort, reindex, drop gaps, infer frequency, aggregate or impute. Equally spaced observations are a documented modeling assumption, not something this vector API can establish.

Integer parameters accept Python/NumPy integer scalars, excluding booleans; reject float-valued integers such as 4.0, negative values and conflicting settings. Validate trend/type/method before estimation. Constant series, complete deterministic fits with effectively zero residual variance, singular designs and nonpositive long-run variance are infeasible tests, not valid p=0/p=1 results.

Use:

- `TypeError` for shape/type/category violations; `ValueError` for invalid labels, ranges, contradictory options and insufficient sample sizes.
- Public `EstimationError(RuntimeError)` for numerical/rank/zero-variance failures, with test/trend/p-or-K context and the original exception chained. Wrap known numerical `ValueError`, `RuntimeError` or linear-algebra errors only within the narrow backend estimation boundary; preserve unexpected programming errors. Do not depend on a private ARCH exception import.
- Public `SpecificationWarning(UserWarning)` when a default K/bandwidth is capped. The same note remains in the result. Warn about documented poor conditioning when estimation remains full-rank; do not silently change the requested model.

For copied/scaled internal working data, document zero-variance and rank thresholds relative to data scale, so changing measurement units does not change feasibility decisions. No arbitrary absolute tolerance should classify all very-small-unit economic series as zero. Reject unsupported ARCH runtime versions before estimation rather than assuming an 8.x or development release is compatible with the initial pinned backend.

## 9. Acceptance tests and dependencies

### Required tests before a first release

1. **MIC identities:** independently hand-checked tiny examples for every lag row and lagged-level sum; p=0, K=0 and selected K; common n versus deliberately wrong n/df divisors; MAIC multiplier 2 and MBIC log(n); conventional AIC/BIC score identities, independently fitted matched OLS likelihood comparisons, and examples with different conventional/modified choices. Include exact/near ties and the saved boundary fixture. MAIC must select 3 there; bootUR's 2 is a documented negative comparison, not a target to imitate.
2. **Detrending:** none/intercept/trend, full-series OLS before trimming, normalized-time invariance, level-shift invariance for `c`/`ct`, linear-trend invariance for `ct`, sign and unit-scale invariance with the expected 2·log|scale| shift for all four criteria. Constant/deterministic residual failures must be explicit. No unsupported invariance is required for `n`.
3. **Backend composition:** for both MIC methods and all supported trends, chosen p passed to public released ARCH on **original** y, fixed final n=N−1−p, and no accidental second detrending. Verify explicitly chosen p=0 and larger orders, including selected orders beyond 4. Include a guard that production imports never reference private/development ARCH selectors.
4. **Fixed-test regressions:** preserve the saved ADF/KPSS/DFGLS matches and the PP finite-sample differences; regression access/covariance labels, calibration values and KPSS None regression. Test PP `n` against released ARCH/independent formula because urca does not offer that specification.
5. **Decisions/results:** correct tails, reversed KPSS null, all three alpha values, exact critical-value equality, approximation/boundary labels, p-value/critical-value disagreements near thresholds, conflicting options and immutable metadata.
6. **Feasibility/errors:** shortest allowed samples and one below, all caps and explicit invalid K, b=0/b=m−1/b=m, fixed versus automatic orders, rank deficiency, nearly deterministic/scaled inputs, missing/complex/boolean/2-D data and backend numerical errors. Automatic Hobijn/Schwert behavior and reduced-cap warnings need independent checks; the old study did not establish them.
7. **Cross-language fixtures:** R urca fixed DFGLS/PP; bootUR 1.0.5 unscaled two-step MAIC/MBIC with its penalty and lag-label differences explicitly asserted. Use the saved vectors and original-style/PQ diagnostic distinction. Live R regeneration is an optional developer job with recorded versions; normal Python tests use reviewed frozen fixtures and released ARCH only. Extend MBIC fixtures to independently verify numerical criterion values, since old R results expose only selected orders.
8. **Platform/version matrix:** tested lower dependency bounds, Python minimum and current teaching version, macOS/Linux and at least one non-Apple BLAS. Start with absolute 1e-8 for saved statistic/criterion comparisons and exact integer choices; use documented scale-aware checks for variance/sum identities. Never loosen tolerances just to erase an unexplained package difference. Explain that the six fixture series are not a size/power experiment.

### Proposed runtime support

**Target minimum Python 3.11**, approved by the instructor; compatibility must not be claimed until tested. ARCH itself supports 3.10, but the exercised NumPy 2.4.4, SciPy 1.17.1 and Pandas 3.1.0 metadata require at least 3.11. No lower-Python validation was performed. The minimum is a proposed support promise requiring CI, not a claim that 3.11 has already been tested.

Proposed direct dependencies: `arch==8.0.0`, `numpy>=1.26.4,<3`, `scipy>=1.11.4,<2`, `statsmodels>=0.14.4,<0.16`. Statsmodels is named because the public regression-result contract exposes its objects. Pandas is brought in by ARCH/Statsmodels; no separate optional input adapter or additional plotting dependency is needed. Proposed lower bounds are conservative **untested bounds** until the acceptance matrix passes; the confirmed numerical stack remains the exact versions in Stage 2B metadata. Pin the teaching environment separately when packaging is eventually authorized; no such file is created here.

R, urca, bootUR, JSON tooling for R, pytest and the old development ARCH snapshot are validation/development tools, never runtime requirements. No external service, package compiler or internet connection is needed to run the eventual library once its dependencies are installed.

## 10. Proposed implementation sequence and approvals

Proceed only after this specification is reviewed; this task stops at the following plan:

1. **Contract and fixtures:** instructor confirms the choices below; freeze the reviewed Stage 2 fixtures and complete test cases for K/bandwidth defaults, errors and decisions. Establish the minimal package/test layout only when implementation is authorized.
2. **Independent MIC core:** input normalization, OLS deterministic removal, candidate construction, QR/SVD estimation, MAIC/MBIC, tie/cap rules and immutable selector result. Pass formula, invariance, boundary and independent fixture tests before integration.
3. **ADF/DFGLS public adapters:** fixed-order first, then PQ-selected order passed to released ARCH; regression access, sample/calibration metadata and decisions. Pass the public-composition regression tests.
4. **PP/KPSS public adapters:** distinct bandwidth handling, stationarity-null display, calibrated results and robust auto-bandwidth/error tests. Keep all source-documented inter-package differences visible.
5. **Teaching/release verification:** short examples/summary display, cross-language regeneration where available, lower-bound and cross-platform matrix, and final instructor review. No scope expansion is implied.

**Already decided by the instructor:** released ARCH 8.0.0; public final-test APIs; independent MAIC/MBIC; PQ OLS selection; no rescaling by default; no complete M family or excluded test classes. These decisions do not need to be requested again.

**Stage 3A instructor decisions:** namespace and five functions approved; default MAIC, explicit MBIC/AIC/BIC supported; fixed-only `ctt`; proposed maximum-lag/bandwidth rules, warnings and feasibility checks approved; critical-value decisions at 1%, 5%, 10% approved; Python 3.11+ is a testing target; released ARCH 8.0.0 exclusively. Dependency lower bounds still require testing. Only the foundation and selector are authorized now; stop for review before any test adapter.

**Recommendation:** move to implementation of this approved architecture after those API/default choices are accepted. No numerical blocker remains in the audited evidence. Production release still requires the acceptance tests above; the incomplete full-text Perron–Qu review remains a stated provenance limitation rather than a basis for claiming exact paper replication.

## Evidence and sources

The audit uses the local CSV/JSON files identified above, the two reviewed scripts, and sources already preserved during Stages 1–2. No upstream source was changed or new upstream release assumed.

- [ARCH audit](ARCH_Audit.md), [MIC methodology](MIC_Methodology.md), [numerical study](Numerical_Validation.md).
- [Original Stage 2 manifest](../validation/results/manifest.json), [comparison rows](../validation/results/comparisons.csv), [candidate equations/data](../validation/results/mic_candidates.csv), [independent audit report](../validation/stage2c_audit.json).
- Released ARCH 8.0.0 installed `METADATA`, public test signatures/source and `unitroot.py` digest `c9ada36525f293df31297085b4afbafdd6a75935af05ad098d4076cc568bb1e3`; [tagged source](https://github.com/bashtage/arch/blob/v8.0.0/arch/unitroot/unitroot.py).
- Saved bootUR 1.0.5 source `adf_selectlags_cpp`, `maic_cpp`, `mbic_cpp`, `adf_onestep_tests_parest_all_units_cpp`; archive digest and primary paper references are in the methodology document.
- Installed package metadata for the tested NumPy/SciPy/Pandas/Statsmodels versions supports the dependency discussion; it does not substitute for a minimum-version test run.

## Stage 3A numerical foundation details

The selector uses float64 working copies and independent scaled SVD least squares (SciPy `gelsd`). Scale levels before detrending, then normalize detrended levels again before constructing candidates; undo the scale in reported SSR, variance and level sums. Normalize each design column to Euclidean norm one before rank determination. The relative singular-value cutoff is machine epsilon times the larger design dimension. Explicit residual norm at or below the same relative cutoff times response norm is numerically zero. Deterministic-removal residual norm uses the analogous cutoff against full input norm. These checks are unit-scale relative, not absolute variance thresholds. A full-rank condition number above 1/sqrt(machine epsilon) emits `SpecificationWarning`, with the same note retained in the result.

Inputs with finite float64 levels but unrepresentable squared-unit diagnostics (variance/SSR/level sum underflow to zero or overflow) raise `EstimationError` with rescaling guidance; a finite log criterion alone must not hide an invalid zero/infinite diagnostic. The supported numeric range is thus explicitly bounded by representable scalar diagnostics. Exact geometric/polynomial sequences and deterministic fits can be infeasible even when the sample-length bound passes. No failed candidate is skipped. Zero-variance, rank and numeric-range errors must name the selection context and preserve a numerical cause when one exists.

Stage 3A adds only package/test configuration, the independent selector, immutable diagnostics, tests and a validation report. It implements no `adf`, `dfgls`, `pp`, `kpss`, `UnitRootResult` or rejection/long-run-bandwidth computation yet. ARCH is pinned for later adapters but is neither imported nor called by the selector.


## Stage 3B integration details and scope

The fixed-lag and automatic `adf` / `dfgls` adapters are now implemented. The Stage 3A selector calculations and its validation module are unchanged. Its existing API-availability test is updated to permit the newly authorized exports and public ARCH imports in the adapter module only; no selector/fixture test is removed. PP and KPSS functions remain absent.

Adapters validate and copy original levels, check deterministic degeneracy without performing GLS detrending, resolve fixed/automatic order and test-specific feasibility, then pass `data, lags=p, trend=trend` through the public released-ARCH constructor. Backend imports are lazy, so using only the independent selector does not import ARCH. The active `arch.__version__` must be exactly 8.0.0. Automatic defaults call the unchanged selector with the effective explicit K and restore requested/default cap metadata with a frozen-record replacement; candidate calculations are untouched.

The public final regression design is checked for finite data, expected dimensions/degrees of freedom, scaled rank and positive residual variance. A backend that drops columns due to numerical scale is rejected rather than silently accepted or re-estimated on rescaled input. Runtime numerical warnings and Statsmodels' specific `SingularMatrixWarning` are converted to chained `EstimationError` inside the estimation boundary; unrelated warnings and unexpected programming errors are preserved. No private ARCH exception/selector/calibration import is used.

`UnitRootResult` has frozen scalar attributes and copied read-only critical-value/version mappings; its public Statsmodels regression object is not deeply immutable. It reports all requested test/sample/selection/calibration metadata. Shared bandwidth/type fields remain None in Stage 3B and implement no PP/KPSS behavior. Only the left-tail rule is supported at present. Alpha is converted to a Python float then must equal 0.01, 0.05 or 0.10 exactly; nearby noncanonical values are rejected, not used for interpolation. Equality of the statistic to its critical value does not reject. Approximate p-values at 0 or 1 are flagged as boundary values without being used to override the decision.

`summary()` prints H0, the deterministic specification, fixed/automatic lag choice, distinct final/selection observation counts, statistic, approximate p-value, backend critical values and a nominal 5% decision. It explicitly says failure to reject does not establish the null. Public regression p-values are not substituted for unit-root inference.

See `Stage3B_Validation.md` for executable examples, test coverage, numerical comparisons and limits. Stage 3B stops for review before PP/KPSS or any scope expansion.


## Stage 3C integration details and scope

The approved public signatures are now implemented:

```python
pp(y, *, trend="c", bandwidth=None, test_type="tau") -> UnitRootResult
kpss(y, *, trend="c", bandwidth=None) -> UnitRootResult
```

PP supports n/c/ct and tau/rho; KPSS supports c/ct. Bandwidth is never an
augmentation order: both return `lags=None`, `lag_selection=None`, effective
`bandwidth`, method `fixed`, `schwert_capped` or `hobijn`, and distinct
`bandwidth_requested` / `bandwidth_uncapped` fields. Requested is the explicit
integer or None; uncapped is the PP automatic rule before the feasibility cap
or None. The `schwert_capped` method label denotes the capped rule even when
the cap does not bind; only a binding cap emits a warning and note.

For N levels, PP covariance sample m=N-1 and automatic rule
b=min(ceil(12*(N/100)**(1/4)), N-2). KPSS m=N and its default passes lags=None
to released ARCH's public Hobijn-style selector. Fixed b requires 0<=b<m;
KPSS's resolved automatic b must satisfy the same bound, otherwise raise
EstimationError with fixed-bandwidth guidance. ARCH automatic numerical
failures are chained, with no replacement selector or silent fallback.
Minimum lengths are d+3 for PP and d+2 for KPSS, where d is the number of
stated deterministic regressors. Additional data/rank/numeric checks apply.
Zero is valid, while legacy -1 is not. No selector methodology changed.

PP returns the public HAC/Newey-West regression, final sample [1,N), nobs=N-1,
and statistic-specific ARCH ADF-t/tau or ADF-z/rho calibration. KPSS returns
regression=None, sample [0,N), nobs=N, and ARCH's simulated quantile-table
interpolation. Its null is stationarity around the specified deterministic
terms; its right-tail decision uses statistic > critical_value. Other tests
remain left-tailed. Exact equality never rejects. Result summaries report
bandwidth rather than augmentation lags and retain the distinction between
failure to reject and acceptance. The regression object remains mutable.

Endpoint approximation notes are recorded for PP p=0/1 and KPSS p=0.0001/1.
They describe calibrated endpoints rather than exact tail probabilities.
Only the exact ARCH KPSS historical default-change DeprecationWarning is
suppressed for automatic calls; unrelated warnings remain visible. Public
numerical warnings/errors are handled consistently with Stage 3B.

See Stage3C_Validation.md for comparisons, coverage, reproducibility and
remaining compatibility limits. Stop after this stage for instructor review;
release verification and compatibility claims remain separate work.
