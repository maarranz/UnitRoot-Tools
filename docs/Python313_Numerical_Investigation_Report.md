# Python 3.13 qualification and numerical fitting proposal

30 October 2026 · Mac M1 Pro · investigation complete; no production replacement implemented

**Python 3.13 is qualified on this Mac with the preserved numerical stack.** All 1,367 package tests, 26 research tests and independent validators pass. The migration does not remove the extreme-scale fitting error or the Linux–Mac innovation discrepancy. Column-scaled QR and column-scaled SVD both achieve high accuracy; these experiments favor **scaled SVD as the initial default candidate**, with scaled QR as an independent reference/alternative, rather than attributing the remedy to QR alone. Implementation requires separate approval.

## A. Preservation and documentation audit

The mandatory read-only audit preceded environment creation, installation, tests and file writes. Branch `main`, HEAD `37fa6efab21a5f332c062c263007411cdb46623e`; six tracked research files have existing modifications, and eleven untracked paths contain prepared research files, reports/evidence and a notebook. Nothing was stashed, switched, staged or discarded.

The inventory covers **820 pre-existing non-Git repository files**, including ignored files, and **36 documentation files**. `initial-read-only-inventory-raw.json` records exact paths, sizes, timestamps, SHA256, complete status and tracked diff. `documentation-inventory.json` adds documentation creation dates; `files-created-or-modified-2026-10-30.json` identifies today's files. Final verification confirms all 820 file hashes **and modification times** unchanged, Git status/diff identical, and all **9,670 files** inventoried in Python 3.12 unchanged with no additions. See `preservation-verified.json`.

The relevant existing analysis was reused: [Stage 2 Colab report](archive/Stage2_Qualification_Report.md), [prepared Mac validation](archive/Prepared_Validation_Report.md), [innovation diagnosis](Innovation_Discrepancy_Report.md), [original numerical diagnosis](/Users/maarranz/Documents/ChatGPT/UnitRoot-Tools/handoff/original-linux-handoff/Diagnostic_Report.md), and the 28 October deterministic-selector correction/design documents. These already establish raw column imbalance and ordinary-scale agreement. New work targets interpreter compatibility, covariance mapping and solver cost; it does not rewrite their findings.

## B. Separate Python 3.13 environment

New executable:

`/Users/maarranz/Documents/ChatGPT/UnitRoot-Tools/work/python313-qualification-2026-10-30-01/env313/bin/python`

CPython **3.13.18**, macOS 26.7.1 arm64; NumPy **2.4.4**, SciPy **1.17.1**, statsmodels **0.15.0**, ARCH **8.0.0**, pandas **3.1.0**. NumPy/SciPy build metadata identifies Accelerate. Complete freezes, build configurations, executable paths and commands are recorded in `metadata.log`, `python313-freeze.txt`, `environment-creation.json` and `installation-commands.json`.

Following the existing Stage 4 convention, `uv venv` used the already available standalone 3.13 runtime, with **no system site packages**. All dependency versions match Python 3.12's full freeze; no unrelated ATSE environment was used. Runtime packages came from binary wheels. The unchanged current package was built from an **external source snapshot** and normally installed as a wheel with genuine metadata. `pip check` passes. Installed runtime bytes match current source exactly. No editable installation, repository build output or dependency-constraint change occurred.

This environment is the qualified development/testing target; no global interpreter setting or existing configuration file was changed. Its ordinary installation is a snapshot: future authorized source changes require a fresh external wheel installation or separately approved editable-install arrangements. The old 3.12 installed wheel predates Enders; historical tests used explicit source paths. Accordingly, the narrow 3.12 comparison used the preserved interpreter plus identical current source, while 3.13 used its genuine installed package. Python 3.12 itself was not upgraded or reinstalled.

## C. Compatibility results

Tests ran from the external snapshot against installed 3.13 code, with the existing strict marker/warning policy, bytecode/user-site disabled and numerical threads limited to one. No warnings were filtered to obtain success.

| Check | Python 3.13 result |
|---|---|
| Existing package suite, unchanged | **1,367 passed**, no failures/skips |
| Existing research suite, unchanged | **26 passed** |
| Independent arithmetic/tree, Enders, deterministic validators | All passed |
| Original nine portable Linux fixtures | Passed at unchanged 1e-8 tolerance and exact hash/discrete requirements |
| Nine portable records: fresh 3.12 versus 3.13 comparison | Entire records exactly equal, including innovation/input hashes, paths, decisions, floating fields and warnings |
| Previously discrepant validation key | 3.12/3.13 innovation bytes and complete record exactly equal; floating difference zero; warning lists both empty |
| Exact extreme-scale diagnostic: 3.12 versus 3.13 | Numerical fit records and selection results exactly equal |

The Mac innovation hash remains `44f1ffd79c4b73f67e60ea8fce1d338a1f881a971b083ec0ee0f53402c2379dc`, differing from original Linux `4381d297820c9e18827f3ac699f3007d8a56f99d72e602c488740963049ddfc5`. Both interpreters retain level hash `eea9b470f1a8688d1177d2bde39e1d2b2d0beb543f3ceb7b392d1a9adf4a3397`. The separately explained log1p discrepancy and strict gate remain unchanged. No 400-row rerun, pilot or main workload was needed. Supplied isolated main-phase portable keys were replayed only as fixtures.

No uncaptured test warning was reported. Pip's nonfatal cache-access warning is preserved in command/stderr logs; it is distinct from scientific warnings. The cold Matplotlib/pyparsing warning mechanism found on Colab was not reproduced in this matched Mac environment.

## D. Accuracy, covariance and cost

Actual path: adapter `_estimate` passes raw levels to ARCH; ARCH constructs the ADF difference design and calls Statsmodels OLS. Default `OLS.fit()` uses an SVD pseudoinverse of the **unscaled** design. It is not a normal-equation coefficient solve. Residuals/RSS derive from those coefficients; covariance comes from the pseudoinverse and residual variance; ADF uses the level-coefficient t statistic. `_deterministic._scores` correctly uses RSS/n and the prescribed parameter penalties. Normalized rank checking after fitting does not repair raw-fit accuracy.

Saved exact binary64 designs were reused, including the ±1e-12 p=1 cases. Their raw ct condition is about 1.16e14, versus normalized 3.82; all columns remain above the current raw pseudoinverse cutoff. The Mac raw RSS excess remains **8.64784e-6**, inside its unchanged 1e-5 test budget; Colab's documented **2.51061e-4** exceeds that budget. Interpreter migration does not explain or remedy this solver/build sensitivity.

New 80/120-digit Decimal coefficient/covariance references agree at binary64 output precision across nine saved fits. Against the 120-digit reference:

| Diagnostic | Scaled pivoted QR | Scaled SVD |
|---|---:|---:|
| Maximum relative RSS error | 2.22e-16 | 2.22e-16 |
| Maximum covariance error, standardized by sqrt(Cii*Cjj) | 1.06e-15 | 8.64e-16 |
| Maximum relative standard-error error | 4.44e-16 | 4.44e-16 |

QR t statistics differ from the reference by at most 5.33e-15. The current raw fit's standardized covariance error reaches **3.34e-3** on the scaled ct case: RSS accuracy alone understates its inferential error. Scaled QR also agrees with current ordinary-scale results on 27 saved p=0 fits: maximum t difference 1.11e-13; all checked decisions unchanged. Raw Statsmodels QR succeeds on these extreme fixtures, but its covariance route forms/inverts R'R; simply switching a method flag is not a complete integration/robustness design.

Single-thread warmed medians, microseconds/fit (7 batches ×100; not simulation timing):

| Design | Current backend + rank screen | Raw SVD arithmetic | Scaled QR diagnostic | Scaled SVD diagnostic |
|---|---:|---:|---:|---:|
| n=118, p=1 ct, ordinary | 102.2 | 57.1 | 85.0 | 60.1 |
| Same, extreme scale | 100.3 | 57.3 | 84.3 | 59.7 |
| n=500, p=0 ct, ordinary | 109.6 | 65.0 | 92.3 | 78.6 |

All include a common normalized rank screen and coefficient/RSS/covariance work; the backend additionally creates its full public object. QR's rank-screen SVD plus QR factorization costs more than scaled SVD here. These diagnostic timings establish no full-API speedup or cross-platform winner.

Near-collinear full-rank probes remain ill-conditioned after scaling: at normalized condition ~3.43e12, both solvers have unstable large coefficients and differing residual arithmetic. Both reject the tested rank-deficient/zero-column cases under the existing rank threshold. Scaling cannot cure intrinsic identification problems; no reduced-rank statistical fallback is justified.

## E. Implementation proposal — not applied

Recommend a **shared, column-equilibrated full-rank fitting layer**. Scaled SVD is the stronger initial default candidate on present accuracy/cost evidence; scaled pivoted QR remains a viable alternative and independent validation reference. Require Linux/Colab covariance and full-result-object validation before selecting/authorizing the production solver.

For positive column scales D and optional positive response scale c, let Z=X D^-1 and v=y/c. With pivoted QR, Z[:,p]=QR; solve g[p]=R^-1 Q'v and map **beta=c D^-1 g**. Map covariance as **Cov(beta)=sigma² D^-1 P(R^-1 R^-T)P' D^-1**, where sigma²=RSS/(n-k), undoing pivot order. Scaled SVD substitutes V diag(1/s²)V' for the middle matrix. Use triangular solves/SVD factors, not a binary64 inverse of X'X. This is column scaling only: no centering, detrending, term removal or sample change. [SciPy's QR convention](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr.html) specifies the pivot mapping.

Proposed source scope:

- New private full-rank fitting helper for coefficients, fitted values, residuals, normalized covariance, covariance, standard errors, t values and rank/condition diagnostics.
- `_adapters.py`: ADF branch of `_estimate` and `_check_regression`. Preserve raw named design, sample, lags and public regression behavior. Recompute statistic **and matching p-value** from the stable fit while retaining the same ARCH 8.0.0 calibration/critical-value definitions. Never attach corrected RSS to an otherwise inaccurate ARCH result. Do not silently change DFGLS/PP/KPSS or the separate lag selector.
- `_enders.py`: use the same stable layer for unrestricted, restricted `_joint` and `_auxiliary` fits; preserve restrictions, df, tables, branches and tests.
- `_deterministic.py`: retain `_scores`, parameter counts, common sample, exact tie rule and strict decisions. Strengthen consistency checks if needed; changing IC formulas is not the remedy.

Preserve the present normalized-rank rejection threshold and conditioning warnings. Use SVD to diagnose difficult designs; do not silently discard columns, alter df or accept a pseudoinverse reduced-rank model. Any changed acceptance/warning policy needs explicit review.

For RSS, compare direct raw-unit residuals against scaled/projection calculations. They agree within 2.22e-16 relatively on the saved nine fits, but near-singular probes show that cancellation still matters. Scale residual norms to avoid overflow/underflow, verify representable positive RSS and variance, and reject out-of-range fits rather than clipping or inventing residual variance. IC uses RSS/n; coefficient covariance uses RSS/(n-k). If response normalization is used, restore original RSS/log-likelihood units, including the common score shift. No compensation or RSS-only repair can fix inaccurate coefficients.

Keep public signatures and regression properties coherent: named params, raw model arrays, residuals/fitted values, df, cov_params, bse, tvalues, likelihood, prediction and restriction tests. Existing objects expose Statsmodels behavior; preserving that contract needs a reviewed adapter/result construction, not undocumented field patching or a new look-alike object. Existing backend parity tests may expose intentional numerical changes; preserve current tolerances and review any contract conflict rather than weakening tests.

Required future tests: exact saved ±scale arrays; high-precision coefficient/covariance references; RSS/IC scale identities and ranks; pivot reversal; p=0/positive-p/ctt; full/auxiliary/restricted fits and Wald/RSS agreement; true/near rank deficiency; zero residuals; overflow/underflow; ordinary backend parity; near-tie and critical-boundary decisions; regression-object properties; source/wheel installation; full Linux/Mac/Colab suites. Preserve the RNG requirements separately.

**Decision:** adopt this qualified 3.13 environment for the next authorized development work; preserve 3.12. Approve neither a production QR replacement nor full Monte Carlo execution from this investigation alone.

Evidence is wholly outside existing docs. Commands, freezes, fixture comparisons, numerical results and preservation inventories are adjacent. Two diagnostic-only setup issues are recorded: the old installed 3.12 wheel required the historical explicit source path, and an early JSON round-trip rounded nanosecond timestamps. The retained raw audit capture proves that the latter caused no actual file-time change; `preservation-verified.json` is authoritative over the preliminary preservation output.
