# ARCH unit-root audit — Stage 1

Audit date: **27 October 2026**. Scope: source and documentation review, with a targeted Statsmodels comparison and R follow-up agenda. This is a static audit, not numerical certification. No UnitRoot-Tools functions, packaging or environment files were created; no installation, commit or push was performed.

## Decision

Reuse ARCH's released public ADF, DFGLS, PhillipsPerron and KPSS classes as candidate backends, subject to later numerical validation. **ARCH 8.0.0 has no MAIC selection option; the inspected development snapshot does.** Do not describe either release DFGLS or development MAIC as a verified implementation of the complete Ng–Perron test family. Preserve deterministic specification, lag/bandwidth meaning, sample size and calibration provenance in any future results interface.

The main implementation gap for a release-only dependency is a separately specified and validated modified lag selector. A full Ng–Perron M-test implementation is a separate gap. Do not implement either until the intended original versus modified detrending procedure is settled.

## 1. Versions and evidence boundaries

| Surface | Observed version/revision | Interpretation |
|---|---|---|
| ARCH PyPI metadata | 8.0.0, first uploaded 2025-10-21 | Latest release reported by PyPI at inspection, not a statement about an installed local version. [A1] |
| ARCH GitHub latest release | v8.0.0, published 2025-10-21 | Release notes emphasize compatibility. Tagged source is the release baseline here. [A2] |
| ARCH `latest` documentation | Page banner says 7.2.0 | This moving URL is behind the release baseline; it cannot establish 8.0.0 or development functionality. [A3] |
| ARCH development | `09a6be5f1b490346d8ba5d3b5850991734aaaf3a`, commit dated 2026-10-26 | Source and tests fetched through `main` were re-fetched using this SHA and matched byte for byte. MAIC findings below apply to this snapshot only. [A4] |
| Statsmodels PyPI/GitHub release | 0.15.0 / v0.15.0, 2026-08-27 | Main comparison baseline; 0.14.6 source also read for API history. [S1–S3] |

ARCH's v8.0.0 change log raises minimum Python to 3.10, switches to Meson and records compatibility/typing/documentation work. It does not announce MAIC. Direct inspection of v7.2.0 and v8.0.0 selectors confirms three methods in both. A diff of these files shows typing, conversion and numerical-code maintenance, so unchanged method names must not be mistaken for proof of identical numerical behavior. No complete historical release bisect was performed. [A5, A6]

## 2. Public API and results

The v8.0.0 and inspected development `arch.unitroot.__init__` export `ADF`, `DFGLS`, `PhillipsPerron`, `KPSS`, `ZivotAndrews`, `VarianceRatio`, `auto_bandwidth`, `engle_granger` and `EngleGrangerTestResults`. There is no exported Ng–Perron M-test class. The unitroot implementation files inspected likewise contain no such implementation. This is a bounded source finding, not a search of every repository file. [A7]

Common test properties include `stat`, `pvalue`, `critical_values`, `lags`, `nobs`, `trend`, `valid_trends`, `y`, `null_hypothesis`, `alternative_hypothesis`, and `summary()`. Computation is lazy in released 8.0.0. Inspecting a result can trigger estimation and feasibility errors. ADF/DFGLS expose `max_lags`, but the selected information criterion is not a public result property in 8.0.0 (`_ic_best` is internal for ADF). [A6]

| Released class | Deterministic terms | Meaning of `lags` | Regression access |
|---|---|---|---|
| `ADF(y, lags=None, trend="c", max_lags=None, method="aic", low_memory=None)` | `n`, `c`, `ct`, `ctt` | Number of lagged differences | Public `.regression`: fitted Statsmodels OLS result |
| `DFGLS(...)`, same selection arguments, `trend="c"` | `c`, `ct` | Number of lagged differences after detrending | Public `.regression`: final no-deterministic ADF regression on GLS-detrended series |
| `PhillipsPerron(y, lags=None, trend="c", test_type="tau")` | `n`, `c`, `ct` | Newey–West bandwidth, not augmentation order | Public `.regression`: level AR(1) fit with HAC covariance |
| `KPSS(y, lags=None, trend="c")` | `c`, `ct` | Newey–West bandwidth | No public `.regression` in inspected class; deterministic OLS fit is local, residuals internal |

Here `n` means no deterministic component, `c` intercept, `ct` intercept plus linear trend, and `ctt` additionally quadratic trend. These sets are test-specific; a universal trend menu would be misleading. Constructor specifications should be treated as fixed: released source contains read-only properties despite examples/docstrings suggesting assignment such as `adf.trend = "ct"`. Reconstruct the object to change specification. Development has additional validation and method normalization; do not assume its input handling exists in 8.0.0. [A6, A8]

## 3. ADF and lag selection

The released regression is

`Δy_t = β y_(t−1) + deterministic terms + Σ(j=1..p) γ_j Δy_(t−j) + u_t`.

The test statistic is the first coefficient's OLS t-statistic, testing β = 0. Fixed `lags=p` bypasses automatic selection. With `lags=None`, candidates 0 through `max_lags` are considered using `aic`, `bic` or `t-stat`. The implementation selects on a common sample trimmed to the maximum candidate order, then refits the final regression using the chosen order; final `nobs` is normally `T−1−p`. [A6: `_df_select_lags`, `_estimate_df_regression`, `ADF`]

Default maximum order is `ceil(12*(T/100)**0.25)`, with a short-sample cap. ARCH's automatic cap is `max((T−1)//2−1, 0)` minus the number of deterministic terms, followed by nonnegative clipping. Statsmodels uses `T//2−ntrend−1`; these can differ by one for even T. Explicit maximum orders and infeasible samples need separate validation. [A6, S2]

ARCH minimizes Gaussian-likelihood criteria with penalties `2p` or `log(n)*p`; omitted common parameter penalties do not affect the minimizing p. Its t-stat selection uses the largest lag with absolute last-lag statistic at least `1.6448536269514722`, with order zero eligible. This is not an arbitrary user-selectable significance level. The low-memory selector uses a different numerical route; its existence is not a guarantee of equality in poorly conditioned cases. Released ADF automatically selects this route above 100,000 observations; released DFGLS uses a threshold of at least 100,000. [A6]

ADF critical values use `mackinnoncrit(..., num_unit_roots=1, nobs=regression.nobs)` and response surfaces in inverse sample size. P-values use the separate `mackinnonp` asymptotic approximation; they are not exact finite-sample probabilities. The source attributes critical-value surfaces to MacKinnon 2010, with the no-constant exception attributed to earlier work. Do not label every calibration “2010 p-values.” [A6, A9]

## 4. MAIC availability and Ng–Perron fidelity

### Availability is version-dependent

**v7.2.0 and v8.0.0:** the public signatures, `_select_best_ic`, and lag-selection routines only support AIC, BIC and t-stat. `method="maic"` is not a released supported selector; manually passing an externally selected integer lag is possible but does not make ARCH itself a MAIC implementation. Statsmodels 0.15.0's `adfuller` also validates only AIC/BIC/t-stat (or `None`). Its source comment citing Ng and Perron does not add MIC support. [A5, A6, S2]

**Development SHA above:** ADF and DFGLS accept `maic`, with uppercase normalized, and both full-memory and low-memory selectors calculate a modified penalty. ZivotAndrews also receives this selector, but break-test behavior was outside this audit. There is no MBIC branch in the inspected selector. [A8]

### What the development code actually computes

Writing `n = T−1−p_max` for the implemented common selection sample, it minimizes

`MAIC(p) = log(sigma2_p) + 2*(tau_p + p)/n`,

with `sigma2_p = SSR_p/n` and

`tau_p = beta_level,p**2 * SS_level / sigma2_p`.

In the full-memory route, `SS_level` is the squared QR diagonal for the lagged-level column after projection on preceding deterministic regressors. In the low-memory route it is the corresponding Schur complement after the level column is normalized; the coefficient and variance scaling compensate, and residual variance is restored to original units before criterion calculation. Thus this is a substantive modified penalty, not AIC renamed MAIC. The implemented `n` should be recorded explicitly rather than interpreting the source docstring's `T−k_max` without accounting for differencing. [A8: lines 148–410]

For ADF the penalty is formed inside the usual regression with the selected deterministic terms. For DFGLS, selection runs on a separately **OLS-detrended** full series with no deterministic regressors in the selection regression; the final test uses **GLS-detrended** data. These preprocessing choices are part of the procedure, not interchangeable details. [A8: `DFGLS._compute_statistic`]

The development test `test_dfgls_auto_maic` checks both memory paths against values attributed in its comment to R boundedur **1.0.3**, on demeaned macro inflation: maximum 14 selects 13 (criterion approximately −10.15394), and maximum 4 selects 2 (approximately −10.18967). These are upstream assertions inspected, not tests executed here. They cover an intercept case and two maxima, not original-paper replication, all trend cases, boundary conditions or size/power. The CRAN page retrieved during this audit showed boundedur 1.0.1, so the precise 1.0.3 comparator provenance remains unresolved. [A10, R3]

### Fidelity verdict

Ng and Perron (2001) distinguish modified information criteria and GLS-based M tests. Their paper's publisher abstract establishes that scope; it is not sufficient for equation-by-equation validation. ARCH development's explicit correction term supports calling it an implemented MAIC selector, but **full Ng–Perron fidelity is unverified**. [P1, A8]

Before claiming fidelity, specify and check: original GLS-based selection versus the Perron–Qu OLS-selection modification; initial-observation treatment; common selection sample and residual-variance divisor; lagged-level sum-of-squares projection; deterministic terms; finite-sample critical-value source; and the final statistic. The ARCH DFGLS documentation cites Perron and Qu (2007) for OLS detrending during lag selection. bootUR documentation refers to a Perron–Qu (2008) correction; resolve the precise publication/version rather than assuming these descriptions identify identical algorithms. [A8, R2]

The inspected ARCH API does not provide the full `MZα`, `MZt`, `MSB`, `MPT` family. A DFGLS statistic with MAIC is not those M statistics. Any later implementation requires its own formula audit and independent reference results; MAIC availability alone does not close that gap.

## 5. DF-GLS details and calibration discrepancy

Released DFGLS uses local-to-unity GLS detrending with `c_bar = −7` for an intercept and `−13.5` for intercept/trend, and transformation coefficient `1+c_bar/T`. First rows of the transformed y and deterministic matrix are retained untransformed. The final regression is an ADF without deterministic terms. Automatic lag selection instead uses OLS-detrended data. In 8.0.0 an initial `ADF(self._y, ...)` selection is made and then overwritten by that OLS-detrended selection; development removes the redundant initial selection. Fixed lags avoid that branch. [A6, A8]

**Documentation/source discrepancy:** the 7.2.0 documentation and released class docstring say intercept-case critical values and p-values are identical to ADF. Actual code routes both `c` and `ct` to distinct `dfgls_*` arrays. Their intercept asymptotic critical values are approximately **−2.56782, −1.94363, −1.61998** at 1%, 5%, 10%, versus ordinary intercept-ADF values approximately **−3.43035, −2.86154, −2.56677**. This is also not exact equality with the ADF no-constant arrays. Use implemented calibration provenance, and flag the prose claim rather than repeat it. The DFGLS table describes simulation-based approximation using MacKinnon-style methodology. Simulations were not rerun and simulation-generation code was not inspected. [A9, D1]

## 6. Phillips–Perron

PP fits `y_t = rho*y_(t−1) + deterministics + u_t` without lagged differences, then corrects its statistics using Bartlett/Newey–West long-run variance. `lags=None` defaults to `ceil(12*(T/100)**0.25)` using original series length. Both `test_type="tau"` (default) and `"rho"` exist: source computes both statistics and chooses the associated `adf-t` or `adf-z` calibration. Ordinary regression p-values are not PP unit-root p-values. [A6, A11]

`.regression` exposes the level regression with HAC standard errors using the same bandwidth. The PP correction itself also uses quantities from the conventional OLS fit, including its standard error and residual variance. Therefore neither the returned regression t-statistic nor a generic HAC coefficient test should be substituted for `.stat`. Explicit bandwidth, trend, statistic type and both effective/original sample sizes belong in later comparisons. [A6]

## 7. KPSS

KPSS reverses the testing direction: the null is stationarity around the selected intercept or linear trend; rejection points against that null. Its statistic is the sum of squared cumulative deterministic-regression residuals, divided by `T**2` and a Bartlett/Newey–West variance estimate. It uses all T observations after deterministic regression. [A6, A11]

In released ARCH, `lags=None` uses the data-dependent Hobijn-style bandwidth rule and emits a DeprecationWarning describing the historical change. `lags=-1` selects the legacy ceiling rule; an explicit nonnegative integer selects that bandwidth. These are covariance bandwidths, not ADF lag choices. The inspected release does not cap the legacy rule to T−1 in the same way as Statsmodels and permits bandwidth equal to residual count; small-sample behavior can therefore differ. [A6, S2]

ARCH interpolates an extensive simulated quantile table. The source reports 100 million replications with series length 2,000; this provenance was read, not independently reproduced. Returned 1%, 5%, 10% critical values for an intercept are **0.7428, 0.4614, 0.3475**, and for trend **0.2175, 0.1479, 0.1193**. `kpss_crit` uses interpolation with endpoint saturation, so extreme p-values are bounded by its table. A helper docstring incorrectly lists critical-value order as 10%, 5%, 1%; implementation uses 1%, 5%, 10% and the class labels match that implementation. [A6, A9]

Statsmodels 0.15.0 instead interpolates the original four quantiles at 10%, 5%, 2.5%, 1% and warns outside the supported 0.01–0.10 interval. Matching KPSS statistics does not imply matching p-values or critical values. [S2]

## 8. Statsmodels comparison and reuse choice

| Concern | ARCH 8.0.0 | Statsmodels 0.15.0 |
|---|---|---|
| Standard ADF | Class, lazy properties | `adfuller`, function; same deterministic menu and AIC/BIC/t-stat choices |
| Fixed order | `lags=p` | `maxlag=p, autolag=None` |
| DFGLS / PP | Public classes | Neither appears in inspected `tsa.stattools` exports; no claim about all third-party extensions |
| MAIC | Absent in release; present in inspected development | Absent from inspected ADF selector |
| Regression | ADF/DFGLS/PP public property | ADF `store=True` supplies `resstore.resols`; `regresults=True` additionally retains autolag fits |
| KPSS bandwidth | `None` automatic; `-1` legacy | `nlags="auto"` or `"legacy"`; explicit integer |
| Result interface | Consistent test objects, but test-specific regression availability | 0.15.0 offers `result_object=True` (`ADFullerResult` / `KPSSResult`); default still legacy tuples with warning |

Statsmodels 0.15.0's result objects unpack only statistic and p-value; other fields use attributes. The documented default migration is 0.16 or after July 2027, whichever is later. Do not assume that new result behavior already applies by default. Inspected 0.14.6 has the older tuple API and lacks this option. KPSS's stored result does not expose a fitted deterministic regression. [S2, S3]

ARCH already uses Statsmodels OLS internally. Use Statsmodels as an independent public ADF/KPSS comparison interface, while recognizing shared regression machinery and MacKinnon table ancestry. Later validation should start with explicitly matched lag/bandwidth, trend and sample before checking automatic selection. Differences in short-sample maximum orders and KPSS calibration are expected reasons for disagreement, not automatically bugs. [A6, S2, S4]

## 9. Reuse and implementation recommendations (no implementation in Stage 1)

1. Prefer composition over copying ARCH internals for standard released tests. Expose backend/version, deterministic terms, selection method, candidate cap, selected augmentation order or bandwidth, effective sample and calibration source. Keep regression access optional and test-specific.
2. For MAIC, choose explicitly between waiting for a release containing the inspected development work, pinning a reviewed development revision, or later implementing an independently validated selector and passing fixed lags to ARCH. No dependency pin is proposed or created by this audit.
3. Separate an original Ng–Perron mode from any OLS-selection/Perron–Qu mode. Define truncation, detrending, variance normalization and reference tables before coding. Never infer M-test support from a DFGLS name or MAIC option.
4. Later verification should compare fixed-order ADF/KPSS with Statsmodels, PP/DFGLS with urca under matched conventions, and modified selectors with a pinned R reference. Include intercept/trend, zero/max selected lags, short/constant/singular inputs, both memory paths, and negative-MA designs relevant to Ng–Perron size distortions. Simulation size/power claims require simulations, not agreement on one dataset.
5. If code is copied later, retain the actual upstream copyright/conditions/disclaimer. The inspected ARCH license permits reuse subject to notice and non-endorsement conditions; no upstream source was copied into this repository. R packages have different license terms, so do not assume their code can be copied under the repository's existing license. This is a source-notice observation, not a legal compatibility opinion. [A12, R1–R3]

## 10. R follow-ups, bounded to documentation review

- **urca 1.3-4:** manual entries `ur.df`, `ur.ers`, `ur.pp`, `ur.kpss` are relevant comparison targets. `ur.df` documents fixed/AIC/BIC selection; `ur.ers` distinguishes DF-GLS from its P-test, with different lag semantics and table provenance. Audit source/sample alignment and statistic variants before using it as an oracle. No R computation was performed. [R1]
- **bootUR 1.0.5:** `adf` documents AIC/BIC/MAIC/MBIC, default MAIC, a default rescaling option and two-step detrending. Its manual says modified criteria apply the Perron–Qu correction. Audit these switches and source before asserting equivalence to ARCH development. Bootstrap/union/panel procedures are later work, not supplied by this audit. [R2]
- **boundedur:** retrieved CRAN page describes bounded-series ADF/M-type tests and MAIC. Inspect the exact selector revision cited by ARCH's test, resolve the 1.0.1/1.0.3 provenance mismatch, and distinguish bounded Brownian-motion calibration from ordinary unbounded unit-root calibration. Its bounded tests are not drop-in calibration for ordinary ARCH results. [R3]

## 11. Inspection ledger and limitations

Local repository: checked `git status --short` (initially clean), remote URL (maarranz/UnitRoot-Tools), directory contents (only LICENSE plus Git metadata), and applicable ancestor/repository AGENTS.md locations (none found). The prior chat could be read, but its assistant messages resolved only to opaque content-reference markers; no missing plan text was reconstructed from them. This request therefore defines Stage 1 scope.

Upstream source files downloaded and inspected outside the target repository, in the chat workspace's `work/arch-audit/`:

- ARCH v7.2.0, v8.0.0 and development: `arch/unitroot/unitroot.py` and `arch/unitroot/__init__.py`. Read target classes, base properties, selectors, regression helper and calibration helpers; diffed 7.2.0 against 8.0.0 and searched for method/test availability. Other classes were scanned only for exports/method support.
- ARCH v8.0.0: `arch/unitroot/critical_values/{dfgls,dickey_fuller,kpss}.py`, `arch/utility/cov.py`, `arch/tests/unitroot/test_unitroot.py`, `doc/source/changes/8.0.rst`, `LICENSE.md`.
- ARCH development: `arch/tests/unitroot/test_unitroot.py`, `arch/unitroot/critical_values/{dfgls,kpss}.py`, `LICENSE.md`. DFGLS coefficient file matched the v8.0.0 file. Main unitroot implementation and test file were verified against immutable SHA downloads. Other main auxiliary files were not separately SHA-rechecked.
- Statsmodels v0.15.0: `statsmodels/tsa/stattools/{_stattools.py,__init__.py}`, `statsmodels/tsa/adfvalues.py`; v0.14.6: `statsmodels/tsa/stattools.py`. Targeted ADF/KPSS/result/selector sections and exports, not the entire time-series library.
- Registry/release metadata: PyPI JSON for both projects; GitHub latest-release responses and main-commit metadata. Statsmodels main was identified but not used as the implementation comparison baseline.
- Documentation: ARCH latest ADF, DFGLS, PP, KPSS and rendered unitroot module, change-log landing page; Statsmodels stable ADF/KPSS; urca manual's relevant test entries; bootUR manual's ADF/selection entries; boundedur CRAN description; Ng–Perron publisher abstract. Some CRAN index/HTML links and an attempted author-hosted paper PDF were inaccessible; successful manual/page alternatives are cited below.

**Not verified:** numerical outputs, package test execution, original-paper equations in full, simulation tables' generation/reproducibility, empirical size/power, R implementation equivalence, bootstrap methods, all release history, or all upstream modules. The default `python3` had no numpy/scipy/pandas/statsmodels/arch installations; no packages were installed and no environment was created. Upstream test expectations are evidence of intended coverage only. Development can change after the pinned snapshot; moving documentation/CRAN pages may be stale. Tagged Git URLs are reproducible locators but tags were not independently authenticated against release wheels/sdist hashes. No external-service writes were performed; the only repository change is this Markdown document.

## Sources

All URLs accessed 27 October 2026. Code symbol names above identify the exact sections supporting findings; development links are immutable.

- [A1: ARCH PyPI metadata](https://pypi.org/pypi/arch/json)
- [A2: ARCH v8.0.0 release](https://github.com/bashtage/arch/releases/tag/v8.0.0)
- [A3: ARCH latest ADF documentation, banner 7.2.0](https://arch.readthedocs.io/en/latest/unitroot/generated/arch.unitroot.ADF.html) and [rendered module](https://arch.readthedocs.io/en/latest/_modules/arch/unitroot/unitroot.html)
- [A4: audited development commit](https://github.com/bashtage/arch/commit/09a6be5f1b490346d8ba5d3b5850991734aaaf3a)
- [A5: v7.2.0 unitroot source](https://github.com/bashtage/arch/blob/v7.2.0/arch/unitroot/unitroot.py)
- [A6: v8.0.0 unitroot source](https://github.com/bashtage/arch/blob/v8.0.0/arch/unitroot/unitroot.py)
- [A7: released exports](https://github.com/bashtage/arch/blob/v8.0.0/arch/unitroot/__init__.py) and [development exports](https://github.com/bashtage/arch/blob/09a6be5f1b490346d8ba5d3b5850991734aaaf3a/arch/unitroot/__init__.py)
- [A8: pinned development unitroot source](https://github.com/bashtage/arch/blob/09a6be5f1b490346d8ba5d3b5850991734aaaf3a/arch/unitroot/unitroot.py)
- A9: released [DFGLS calibration](https://github.com/bashtage/arch/blob/v8.0.0/arch/unitroot/critical_values/dfgls.py), [Dickey–Fuller calibration](https://github.com/bashtage/arch/blob/v8.0.0/arch/unitroot/critical_values/dickey_fuller.py), [KPSS quantiles](https://github.com/bashtage/arch/blob/v8.0.0/arch/unitroot/critical_values/kpss.py)
- A10: [development tests](https://github.com/bashtage/arch/blob/09a6be5f1b490346d8ba5d3b5850991734aaaf3a/arch/tests/unitroot/test_unitroot.py) and [released tests](https://github.com/bashtage/arch/blob/v8.0.0/arch/tests/unitroot/test_unitroot.py)
- [A11: v8.0.0 Newey–West covariance](https://github.com/bashtage/arch/blob/v8.0.0/arch/utility/cov.py)
- [A12: ARCH license](https://github.com/bashtage/arch/blob/v8.0.0/LICENSE.md) and [8.0 change log](https://github.com/bashtage/arch/blob/v8.0.0/doc/source/changes/8.0.rst)
- D1: latest [DFGLS](https://arch.readthedocs.io/en/latest/unitroot/generated/arch.unitroot.DFGLS.html), [PP](https://arch.readthedocs.io/en/latest/unitroot/generated/arch.unitroot.PhillipsPerron.html), [KPSS](https://arch.readthedocs.io/en/latest/unitroot/generated/arch.unitroot.KPSS.html), [change logs](https://arch.readthedocs.io/en/latest/changes.html)
- S1: [Statsmodels PyPI](https://pypi.org/pypi/statsmodels/json), [v0.15.0 release](https://github.com/statsmodels/statsmodels/releases/tag/v0.15.0)
- S2: [v0.15.0 stattools implementation](https://github.com/statsmodels/statsmodels/blob/v0.15.0/statsmodels/tsa/stattools/_stattools.py) and [exports](https://github.com/statsmodels/statsmodels/blob/v0.15.0/statsmodels/tsa/stattools/__init__.py)
- [S3: v0.14.6 stattools](https://github.com/statsmodels/statsmodels/blob/v0.14.6/statsmodels/tsa/stattools.py)
- [S4: v0.15.0 ADF calibration](https://github.com/statsmodels/statsmodels/blob/v0.15.0/statsmodels/tsa/adfvalues.py)
- S5: stable [ADF documentation](https://www.statsmodels.org/stable/generated/statsmodels.tsa.stattools.adfuller.html) and [KPSS documentation](https://www.statsmodels.org/stable/generated/statsmodels.tsa.stattools.kpss.html)
- [P1: Ng and Perron (2001), Econometrica 69(6), 1519–1554, DOI 10.1111/1468-0262.00256](https://onlinelibrary.wiley.com/doi/pdf/10.1111/1468-0262.00256) — publisher abstract inspected; full equations not certified.
- [R1: urca 1.3-4 reference manual](https://cran.r-project.org/web/packages/urca/urca.pdf), relevant entries pp. 43–51
- [R2: bootUR 1.0.5 reference manual](https://cran.r-project.org/web/packages/bootUR/refman/bootUR.html), `adf` and selection options
- [R3: boundedur CRAN page, retrieved version 1.0.1](https://cran.r-project.org/web/packages/boundedur/index.html)
