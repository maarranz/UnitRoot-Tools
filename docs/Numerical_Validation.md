# Numerical cross-validation — Stage 2B

Executed 27 October 2026 on the existing local machine. This is a small reproducible implementation comparison, **not a size/power Monte Carlo study**. Read `MIC_Methodology.md` for the adopted definitions and source review.

## Outcome

Fixed-lag ADF and KPSS statistics agree between ARCH and Statsmodels to floating-point precision. Fixed-lag DF-GLS statistics also agree with R urca; critical values do not. Phillips–Perron statistics differ systematically in finite samples, and separate source-derived formulas reproduce both implementations. Development ARCH's PQ-style MAIC agrees with bootUR on the five base series under the matched specifications. A deliberately selected short diagnostic series demonstrates an actual selected-order difference from bootUR's penalty indexing.

The completed harness has **840 expected-equality checks, all passing at absolute tolerance 1e-8**, within 1,188 recorded comparisons. Selected-order checks are integer equality. “Expected equality” covers only matched quantities and explicit source-derived reconstructions, not all package outputs. Nonmatching quantities remain in the results rather than being forced to pass.

## 1. Code, data and exact versions

The reproducible entry point is [cross_validate.py](../validation/cross_validate.py); it invokes the companion [cross_validate.R](../validation/cross_validate.R). These are research/validation scripts, not production UnitRoot-Tools functions. They install nothing, download nothing, and do not alter either package. Release and development ARCH run in separate child processes so imports cannot accidentally mix.

| Component | Executed version |
|---|---|
| Python | 3.14.10 |
| NumPy / SciPy / Pandas | 2.4.4 / 1.17.1 / 3.1.0 |
| Statsmodels / patsy | 0.15.0 / 1.0.3 |
| Released ARCH | 8.0.0 |
| Development ARCH | Commit `09a6be5f1b490346d8ba5d3b5850991734aaaf3a` |
| R / urca / bootUR | 4.6.1 / 1.3-4 / 1.0.5 |
| Hardware | Apple arm64 |

The development archive was built locally with the **synthetic build label `8.0.1.dev0`** because a GitHub archive has no Git version metadata. This label is not a claimed upstream release; commit and source digest establish identity. The installed development `unitroot.py` is checked against SHA-256 `b6be3547403ab880d425c780e1358c86eb9073a7b903ad8975650b45ad445fc9` before execution. R package and Python release versions are checked too. All packages and temporary environments were kept in this chat's scratch workspace, outside the target repository. Existing R urca and compiler dependencies were used; bootUR was added to a temporary R library.

Source build of bootUR initially failed because the machine's R linker configuration referenced unavailable `/opt/gfortran` libraries. The official CRAN arm64 binary for R 4.6 installed and loaded successfully and was used for the results. No global R configuration was changed. The separately downloaded CRAN 1.0.5 source was inspected, and its source-based selector diagnostics agree with the installed binary.

## 2. Simulated series and settings

The script generates the data once in Python and writes **17 significant digits** to `series.csv`; R reads that same file. There is no attempt to match Python's and R's unrelated random generators. NumPy's generator is explicitly `Generator(PCG64(20261027))`, generating a 5×501 innovation array before constructing:

| Series | N | Construction |
|---|---:|---|
| `rw` | 200 | Last 200 values of a 500-increment Gaussian random walk |
| `rw_negma` | 200 | Last 200 integrated increments `e_t - 0.8 e_(t-1)` |
| `ar07` | 200 | AR(1) coefficient 0.7, initialized at zero; last 200 of 501 levels |
| `trend_ar08` | 200 | Analogous AR(1) coefficient 0.8 plus `2 + 0.03*t` over retained t=1,…,200 |
| `hetero_rw` | 200 | Cumulative innovations with standard deviation 1 for 100 increments and 3 for 100 |
| `mic_boundary` | 64 | Negative-MA integrated series, independent PCG64 seed 457; diagnostic selection described below |

The last series is **intentionally chosen to expose an implementation difference**: a deterministic search over seeds 0,…,999 retains the first series where the aligned and bootUR-style MAIC penalty choose different lags for K=4 and an intercept. It is not an unbiased Monte Carlo draw used to estimate a discrepancy frequency. The script records the search and seed in the manifest.

Settings:

- Fixed ADF: p=0 and 4; deterministic specifications `n`, `c`, `ct`, `ctt` (48 matched package pairs).
- Fixed KPSS: bandwidth 4 and 12; `c`, `ct` (24 pairs). These are long-run-variance bandwidths, not augmentation orders.
- Fixed DF-GLS: p=0 and 4; `c`, `ct`; urca `ur.ers(type="DF-GLS", lag.max=p)` (24 pairs).
- PP: bandwidth 4 and 12; `c`, `ct`; both tau and rho, mapped to urca Z-tau and Z-alpha (24 pairs of each statistic). Fixed positive bandwidth avoids urca's separate zero-bandwidth edge path.
- MIC: K=4 and 12, candidates 0,…,K, `n`, `c`, `ct`; bootUR MAIC/MBIC with `criterion_scale=FALSE`, one- and two-step variants; development ARCH MAIC in both memory modes. DFGLS supports only `c`, `ct`. bootUR's default rescaling is run separately as a sensitivity comparison.

The main script uses public test APIs. Its development MAIC minimum-value checks additionally inspect ADF's internal `_ic_best` and call the pinned `_df_select_lags` helper for DFGLS. Those private calls are diagnostics only, explicitly tied to the source hash; they are not proposed library dependencies. An independent least-squares calculation records all candidate criteria using aligned, bootUR-style, one-step, and original-style GLS-selection conventions.

## 3. Fixed-specification results

| Comparison | Pairs | Maximum absolute statistic difference | Interpretation |
|---|---:|---:|---|
| ARCH–Statsmodels ADF | 48 | 1.78e-15 | Same fixed-lag regression/statistic |
| ARCH–Statsmodels KPSS | 24 | 1.89e-15 | Same statistic at matched bandwidth |
| ARCH–urca DF-GLS | 24 | 4.71e-14 | Same GLS transformation and final fixed-lag regression |
| ARCH–urca PP tau | 24 | 0.0315190 | Different finite-sample correction |
| ARCH–urca PP rho | 24 | 1.2971612 | Different finite-sample correction |

Example: `rw`, intercept, augmentation p=4 or covariance bandwidth 4 as appropriate:

| Test | ARCH statistic | Comparator statistic | Comparator |
|---|---:|---:|---|
| ADF | -1.598060981325 | -1.598060981325 | Statsmodels |
| KPSS | 3.413879019162 | 3.413879019162 | Statsmodels |
| DFGLS | 0.053736157454 | 0.053736157454 | R |
| PP_tau | -1.540143168508 | -1.539734341087 | R |
| PP_rho | -3.950403495937 | -3.948314231579 | R |

### Calibration is a separate comparison

All 48 ADF critical-value triples match exactly; maximum ADF p-value difference is 2.22e-16.

All 24 KPSS statistic pairs agree, but all p-value pairs differ. The maximum p-value difference is **0.0309351**. For the `rw`, intercept, bandwidth-4 example, ARCH returns **0.0001**, Statsmodels **0.01**. Both are endpoint/table-based values; they do not imply that Statsmodels estimates a tail probability exactly equal to 0.01. The script suppresses only Statsmodels' `InterpolationWarning` and records the values. ARCH's 1%, 5%, 10% intercept critical values are **0.7428, 0.4614, 0.3475**; Statsmodels uses **0.739, 0.463, 0.347**.

DF-GLS critical values differ in all 24 pairs because ARCH uses its simulated response surfaces with final regression nobs, while urca uses older intercept formulas and sample-size buckets for trend, with values rounded to two decimals. Maximum absolute critical-value difference here is **0.318943** (10% level). For `rw`, intercept, p=4, ARCH's triple is **−2.668676, −2.048502, −1.731621**; urca's is **−2.58, −1.94, −1.62**. Agreement of statistics does not validate interchangeability of critical values. urca's tested `ur.ers`/`ur.pp` objects do not return corresponding p-values, so none were invented. urca does not supply critical values for the tested PP Z-alpha result; those cells remain blank.

### Why Phillips–Perron differs

Both packages fit the same level AR(1) with deterministic regressors and use the same positive-bandwidth Bartlett/Newey–West residual variance estimate. The difference here is the finite-sample level normalization in the correction.

Let n=N−1, rho be the fitted level coefficient, `gamma0=SSR/n`, lambda2 the long-run variance, `se(rho)` the ordinary OLS standard error, and D the deterministic matrix over these n rows. Put

`M_x = ||residualize(y_(t−1), D)||² / n²`,

`M_y = ||residualize(y_t, D)||² / n²`.

The source-derived reconstruction is

`Z_tau(M) = sqrt(gamma0/lambda2)*(rho−1)/se(rho) − (lambda2−gamma0)/(2*sqrt(lambda2*M))`,

`Z_rho(M) = n*(rho−1) − (lambda2−gamma0)/(2*M)`.

ARCH uses `M_x`. urca's intercept case uses `M_y`; its trend-case explicit M expression equals **`(1−1/n²)*M_y`**. Retaining that factor matters: simply substituting current-level residual sums left differences as large as 0.003114 in the diagnostic rho calculation. The final script includes it and reproduces ARCH to **7.39e-13** and urca to **2.14e-12** across all 48 statistic specifications for each backend. No discrepancy tolerance was relaxed. This is a reconciliation of the implementations, not a determination that one finite-sample formula is universally preferable. Their critical-value sources also differ.

## 4. MAIC/MBIC numerical findings

### Matched two-step OLS selection, rescaling off

For the five base series, ARCH development DFGLS MAIC and bootUR two-step MAIC choose the same order in **20/20** specifications (`c`/`ct`, K=4/12). Including the diagnostic series gives **23/24**. The differing case is `mic_boundary`, intercept, K=4: ARCH chooses **3**, bootUR **2**.

K=12 results (augmentation orders, not final test statistics):

| Series | Trend | ARCH DFGLS MAIC | bootUR two-step MAIC | bootUR two-step MBIC |
|---|---|---:|---:|---:|
| rw | c | 2 | 2 | 0 |
| rw | ct | 0 | 0 | 0 |
| rw_negma | c | 9 | 9 | 8 |
| rw_negma | ct | 9 | 9 | 9 |
| ar07 | c | 7 | 7 | 7 |
| ar07 | ct | 7 | 7 | 4 |
| trend_ar08 | c | 2 | 2 | 2 |
| trend_ar08 | ct | 2 | 2 | 2 |
| hetero_rw | c | 0 | 0 | 0 |
| hetero_rw | ct | 0 | 0 | 0 |
| mic_boundary | c | 9 | 9 | 9 |
| mic_boundary | ct | 9 | 9 | 9 |

Independent aligned MIC calculations match the pinned ARCH selected order in **all 60** ADF/DFGLS specifications. The minimum MAIC value agrees within **2.29e-14** across the full-memory and low-memory routes; both routes choose the same order in all 60 cases. This checks the modified penalty itself, not only a few coincidentally equal selected lags.

The independent calculation reproducing bootUR's penalty indexing matches its two-step MAIC and MBIC selected orders in **36/36 each**, including `n` cases. The recommended full-sum MBIC also selects the same orders in this set, but its criterion values are not identical; no native ARCH MBIC exists to compare. Equal selected orders do not erase the documented equation difference.

### An actual penalty-indexing counterexample

For `mic_boundary`, intercept and K=4, the common sample has n=59 rows. The aligned lagged-level sum is **42.200588309**; bootUR's truncated sum is **40.265549277**. The regressions, coefficients and residual-variance divisors are otherwise the same.

| Candidate p | Aligned MAIC | bootUR-style MAIC |
|---:|---:|---:|
| 0 | 0.584559279 | 0.535380835 |
| 1 | 0.615892442 | 0.566831550 |
| 2 | 0.552905055 | **0.508145123** |
| 3 | **0.552592109** | 0.509366759 |
| 4 | 0.840332034 | 0.783962503 |

Actual ARCH/bootUR outputs choose 3/2, respectively. This supports the source diagnosis that bootUR omits the first matched level from the penalty sum. It does not imply that 1/24 is a population disagreement rate, because the sixth series was selected to demonstrate the issue.

### Differences that must not be combined with that counterexample

- **Original GLS versus PQ OLS selection:** the independent GLS-selection and OLS-selection variants differ in **2/24** MAIC specifications (`rw_negma`, K=12, both `c` and `ct`: original-style GLS selects 8; PQ selects 9). These are script reference calculations, not a native bootUR option comparison.
- **One-step ADF:** after correcting only the reported-order label described next, ARCH ADF MAIC and bootUR one-step MAIC differ in **7/36** specifications. ARCH projects the penalty level on the deterministic regressors, whereas bootUR uses unprojected levels and drops the first matched value. This is not a matched PQ-DFGLS comparison.
- **bootUR one-step output label:** its public “selected lags” equals the actual augmentation order plus the number of deterministic regressors. The R results preserve that number in `reported_lags`; `lags` subtracts 0/1/2 for `n`/`c`/`ct`. Independent candidate curves reproduce the normalized choices in all 36 MAIC and all 36 MBIC cases. For example `rw`, intercept, K=4 selects actual p=4 but reports 5. No installed package code was patched.
- **Volatility rescaling:** turning bootUR's default `criterion_scale=TRUE` on changes selected MAIC orders in **3/36** cases: `rw/n` at K=4 and 12 goes 0→2; `rw_negma/c`, K=12 goes 9→8. This is a sensitivity result; there is no equivalent ARCH setting. Final bootUR OLS ADF statistics are not compared with ARCH's GLS final statistics.

## 5. Reproduction

Run from the repository root. Dependencies must be installed outside the repository. The following is a recipe for a separate temporary environment, not a packaging proposal:

```sh
VALIDATION_WORK=$(mktemp -d /tmp/unitroot-stage2.XXXXXX)
python3.14 -m venv "$VALIDATION_WORK/py"
"$VALIDATION_WORK/py/bin/python" -m pip install \
  numpy==2.4.4 scipy==1.17.1 pandas==3.1.0 \
  statsmodels==0.15.0 arch==8.0.0 patsy==1.0.3

curl -L 'https://codeload.github.com/bashtage/arch/tar.gz/09a6be5f1b490346d8ba5d3b5850991734aaaf3a' \
  -o "$VALIDATION_WORK/arch-dev.tar.gz"
shasum -a 256 "$VALIDATION_WORK/arch-dev.tar.gz"
# Expected: 4b95b6e724abb045fcbfcd0cbfe7f7bfb297b0665d62b33ad81073b3d137e23f

tar -xzf "$VALIDATION_WORK/arch-dev.tar.gz" -C "$VALIDATION_WORK"
SETUPTOOLS_SCM_PRETEND_VERSION=8.0.1.dev0 \
  "$VALIDATION_WORK/py/bin/python" -m pip install --no-deps \
  --target "$VALIDATION_WORK/arch-dev" \
  "$VALIDATION_WORK/arch-09a6be5f1b490346d8ba5d3b5850991734aaaf3a"
```

Install R **urca 1.3-4**, **bootUR 1.0.5** and **jsonlite** in a temporary R library, using compatible official binaries or a working source compiler. The supplied machine already had urca and bootUR's compiler dependencies; bootUR was installed with:

```sh
mkdir -p "$VALIDATION_WORK/R-lib"
Rscript --vanilla -e \
  'install.packages("bootUR", lib=commandArgs(TRUE)[1], repos="https://cloud.r-project.org", type="binary", dependencies=FALSE)' \
  "$VALIDATION_WORK/R-lib"
```

That command requires its dependencies to be installed already and is specific to a platform with suitable R binaries. On another system, install the versioned CRAN source and dependencies instead. Current moving CRAN repositories may later serve a different version; the harness refuses unreviewed urca/bootUR versions. Do not interpret a fresh install without checking versions as an exact reproduction.

Then execute:

```sh
"$VALIDATION_WORK/py/bin/python" validation/cross_validate.py \
  --arch-dev-path "$VALIDATION_WORK/arch-dev" \
  --r-library "$VALIDATION_WORK/R-lib" \
  --out "$VALIDATION_WORK/results"
```

Use `--rscript /path/to/Rscript` if needed. The final run recorded here used the same command with temporary paths in the chat workspace and the default repository `validation/results` output directory. Full Python dependency versions, including transitive libraries, appear in the worker metadata; exact build tools and wheel hashes are not locked. Numerical tolerance rather than bytewise portability across BLAS implementations is intended. The current saved raw results are authoritative for this run.

The process exits nonzero for failed expected-equality checks, unsupported package versions, missing dependencies, or a development-source hash mismatch. Known between-package differences are deliberately not equality assertions. No upstream package's entire test suite was run.

## 6. Saved evidence and limitations

The results directory contains:

- `series.csv`: six input vectors, with blank padding only after the short series ends; `manifest.json` records generator, seeds, lengths, script and raw-file SHA-256 hashes.
- `python_release.csv` (216 rows), `python_development.csv` (120 rows), `r_results.csv` (252 rows): raw numerical observations, selected orders, effective sample sizes, available p-values and critical values. Fields unavailable through a backend remain blank. R MIC `nobs` is inferred from the actual selected augmentation order; it is not a native bootUR output field.
- `mic_candidates.csv` (1,620 rows): every candidate order, coefficient, variance, level sum, tau, MAIC and MBIC for five source-derived diagnostic conventions.
- `comparisons.csv` (1,188 rows): each compared quantity, values, absolute difference and whether equality is expected; `summary.json` gives group counts and maxima. Groups such as DFGLS/KPSS include statistic and calibration comparisons, so their aggregate maximum is not a statistic-only error.
- Python worker metadata and `r_metadata.json`: exact runtime/package versions and development implementation digest. Metadata identifies imported development code separately from its synthetic build label.

Limits: six series, two augmentation orders/bandwidths, finite deterministic menus, and no exhaustive boundary/constant/missing-value coverage; no empirical size or power claims; no complete Ng–Perron M-test or bootstrap validation; no bounded-series test validation. Full Perron–Qu text could not be obtained, as recorded in the methodology document. Cross-package agreement also shares linear-regression/table ancestry and is not independent econometric proof. The candidate selector in this research harness is an explicit diagnostic calculation, not a tested production implementation.

The original Stage 1 audit was left unchanged. This stage adds only the two Markdown documents, the two research scripts and their numerical evidence. No package API, packaging/environment files, commit or push was created in the repository.

## Sources inspected for the numerical reconciliation

- [ARCH 8.0.0 unitroot source](https://github.com/bashtage/arch/blob/v8.0.0/arch/unitroot/unitroot.py), calibration files, and [Bartlett/Newey–West implementation](https://github.com/bashtage/arch/blob/v8.0.0/arch/utility/cov.py).
- [Pinned ARCH development source](https://github.com/bashtage/arch/blob/09a6be5f1b490346d8ba5d3b5850991734aaaf3a/arch/unitroot/unitroot.py); installed Python source digest verified before running.
- [Statsmodels 0.15.0 source](https://github.com/statsmodels/statsmodels/blob/v0.15.0/statsmodels/tsa/stattools/_stattools.py), ADF/KPSS calculation and result-object paths.
- Installed `urca::ur.ers` and `urca::ur.pp` function bodies, inspected directly in R; [urca reference manual](https://cran.r-project.org/web/packages/urca/urca.pdf), [versioned source](https://cran.r-project.org/src/contrib/urca_1.3-4.tar.gz).
- [bootUR 1.0.5 CRAN source](https://cran.r-project.org/src/contrib/bootUR_1.0.5.tar.gz), functions and source digest enumerated in `MIC_Methodology.md`; numerical runs used its official CRAN binary.


## Stage 2C audit clarification — 27 October 2026

An independent read-only audit of all eleven saved result files confirms 1,188 comparison rows and 840 passing expected-equality checks. Fresh QR-based calculations reproduce all 1,620 candidate rows, and 120 additional public fixed-lag composition checks pass using **released ARCH 8.0.0**. The 840 are scalar checks on shared fixtures, not independent simulation replications. See [Implementation_Specification.md](archive/Implementation_Specification.md), [audit script](../validation/audit_saved_results.py), and [audit report](../validation/stage2c_audit.json).

Two preservation/interpretation points apply to this historical run:

- `python_development.csv.metadata.json` has `all_distributions["arch"] = "8.0.0"` because its dictionary collapses duplicate installed distribution names. Its `packages["arch"]`, `arch_runtime_version` and verified development source digest identify the imported development code. The dictionary is not a dependency lock or proof that the development worker used the release. Future evidence generation should retain each distribution and origin or use separate environments.
- The research input readers filter blank/NA entries. Only trailing padding exists in the saved short series, independently verified; no interior gaps were dropped in this run. Future readers must reject interior gaps. Normal runtime input validation must never inherit this filtering behavior.

No numerical formula correction is required for this evidence. Original scripts, raw results and their manifest are preserved unchanged. The final production decision is **released ARCH 8.0.0 plus an independent PQ OLS MAIC/MBIC selector**; the development snapshot remains historical comparison evidence only. Automatic defaults, numerical/error boundaries, and minimum-version/platform support still require implementation acceptance tests.
