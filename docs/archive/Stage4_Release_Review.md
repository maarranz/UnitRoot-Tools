# Stage 4 compatibility and release review

27 October 2026. **Compatibility and numerical validation pass in the six
available macOS matrix environments.** Clean artifact installation checks are
recorded separately below. No econometric functionality, public API, formula,
calibration, lag/bandwidth convention, numerical tolerance or package version
was changed. No commit, push, tag or publication was performed.

## Compatibility matrix

Every successful row uses a fresh virtual environment with system site-packages
disabled and released ARCH **8.0.0**. Test tools are pytest 9.1.1 and pytest-cov
7.1.0. The host is macOS 26.7.1 arm64. These are tests on one host, not evidence
for every OS, BLAS, dependency combination or Python patch release.

| Python | Dependency stack | Install mode | Full suite | Numerical audits / examples / pip check |
|---|---|---|---:|---|
| 3.11.18 | current | repository | 1,164 pass | pass |
| 3.11.18 | lower direct bounds | repository | 1,164 pass | pass |
| 3.12.16 | current | repository | 1,164 pass | pass |
| 3.12.16 | lower direct bounds | repository | 1,164 pass | pass |
| 3.13.18 | current | repository | 1,164 pass | pass |
| 3.14.10 | current | editable | 1,164 pass | pass |
| 3.13 / 3.14 | exact lower direct bounds | installation unavailable | not run | NumPy 1.26.4 has no usable wheel for these interpreters |
| Linux / Colab | either stack | unavailable | not run | no Linux host, Docker or Podman available |

Current stack: NumPy **2.4.4**, SciPy **1.17.1**, Statsmodels **0.15.0**,
Pandas **3.1.0**. Lower-bound stack: NumPy **1.26.4**, SciPy **1.11.4**,
Statsmodels **0.14.4**, with Pandas **2.2.3**. The latter tests all three declared
direct numerical lower bounds simultaneously; it is not a test of the oldest
possible transitive dependency set. Complete resolved versions are saved per
row. ARCH remains pinned at 8.0.0 in both stacks.

The 3.13/3.14 lower-stack attempts used binary-only dependency installation;
source compilation was not attempted. This is an unavailable combination,
not a failure of UnitRoot-Tools numerical code. No claim is made that exact old
bounds work on every newer interpreter. Windows, macOS x86_64, free-threaded
Python and other BLAS/OS combinations were not tested. Runtime downloads were
kept in the task workspace and used directly; inability to create optional
global executable links did not affect their use.

## Installation and distribution checks

Repository installation succeeds in five matrix environments; editable
installation succeeds in Python 3.14. Distribution builders use the standard
setuptools backend. Both a wheel and a source distribution build, including a
wheel rebuilt from the source distribution by `python -m build`.

Two additional fresh environments check **actual pip installation** of the
built artifacts: Python 3.14/current from the wheel and Python 3.12/current
from the source distribution. Both run all **1,164 tests**, public-import/
all-five-function smoke checks, installed-document checks and `pip check`.
These are additional runs of the same tests, not additional independent test
cases. Runtime imports resolve to installed site-packages, not the checkout.
The eight successful suite runs total 9,312 passing test executions.

The compatibility runner copies tests and fixtures outside the repository and
uses a separate pytest configuration without `pythonpath=src`. Thus ordinary
installation tests cannot accidentally pass by importing checkout source.
The editable row intentionally imports the checkout through its installed
editable mapping. Artifact tests use the tests and fixtures extracted from the
source distribution, again without the checkout's pytest path override.

Archive inspection verifies:

- All eight runtime modules match repository source bytes.
- Wheel metadata: name `unitroot-tools`, unchanged version `0.1.0.dev0`, Python
  `>=3.11`, MIT license expression and exact runtime dependency `arch==8.0.0`.
- MIT license bytes and license metadata are present.
- README, all Markdown documents and the three executable examples are shipped
  under `share/doc/unitroot-tools`, retaining the docs/examples directory layout.
- Source distribution includes source, tests, reviewed numerical fixtures,
  validation scripts, documentation and license. Machine-specific Stage 4 build
  evidence is retained in the repository rather than recursively archived.

The required packaging correction was to add documentation/example data files
and an explicit `MANIFEST.in`; the prior wheel contained runtime modules and
license but omitted the standalone teaching documents. No runtime dependency
or dependency bound changed. Artifact names, sizes and SHA-256 digests, content
lists and metadata checks are in `validation/stage4/distributions.json`.
Build logs and metadata/rendering checks are retained beside that report.

## Consolidated numerical regression

Each of the six matrix rows executes all of the following, without rewriting
historical results:

| Evidence | Checks per environment | Outcome |
|---|---:|---|
| Full Stage 3 suite | 1,164 tests | all pass; no skips |
| Stage 2 independent audit | 1,188 historical comparisons, 840 expected equalities | historical counts/provenance confirmed; no expected failures |
| Stage 2 independent reconstruction | 1,620 MIC candidate rows, 264 selector outputs, 120 public-ARCH compositions | all satisfy existing checks |
| Stage 3B replay | 312 model pairs, including 240 automatic choices | all selected orders, samples, statistics and inference agree |
| Stage 3C replay | 1,944 scalar rows, including 1,800 expected equalities | no expected failures; package-specific differences retained |
| Executable examples | original ADF/DF-GLS, PP/KPSS and new teaching walkthrough | all run |

The full suite preserves distinct MAIC/MBIC versus conventional AIC/BIC scores,
common-sample OLS/Perron–Qu conventions, explicit/fixed zero orders, automatic
selection metadata, PP tau/rho and deterministic specifications, Schwert PP
bandwidth/caps, public Hobijn KPSS behavior, and fixed covariance bandwidths.
Right-tail KPSS and left-tail ADF/DF-GLS/PP decisions, exact equality and
approximate-p-value boundary notes remain tested. All four tests retain ARCH
8.0.0 calibration. No R/Statsmodels critical values are substituted.

Current-stack numerical output matches the saved Stage 3B/3C ARCH fields
exactly on all four Python versions. Lower-bound ADF/DF-GLS maximum statistic
difference is **2.7924329515371937e−12**, p-value difference
**1.5816514764566136e−13**, and critical-value difference zero. Its largest
statistic difference is fixed ADF, `rw_negma`, `ctt`, p=0, nobs=199
(−11.75225690105483 versus −11.752256901057622).

Lower-bound PP/KPSS maximum saved-ARCH difference is
**7.958078640513122e−12**, for PP rho on `rw_negma`, ct, b=12
(−258.3201102065688 versus −258.32011020656086). Wrapper output equals fresh
public ARCH within each environment exactly. Selected augmentation orders and
bandwidths remain unchanged; independent reconstruction checks also pass.
The saved NumPy configurations show lower-stack OpenBLAS versus current-stack
Apple Accelerate, alongside the dependency version changes. These findings
are consistent with floating-point linear-algebra rounding; they do not isolate
one library as its cause. Both worst-case rows and configuration logs are saved.
No numerical discrepancy required a formula change or tolerance relaxation;
the reviewed absolute **1e−8** comparison tolerance remains unchanged.

The documented urca PP finite-sample differences and Statsmodels KPSS
calibration differences remain visible, rather than being counted as equality
failures. R was not rerun: comparisons use the preserved Stage 2 R evidence.
These six simulated vectors are regression evidence, not size/power validation.

Coverage with branch measurement is **99.50–99.54%** across the matrix. The
minor percentage/statement-count difference comes from interpreter-specific
coverage analysis of unchanged source; the same two defensive lines and one
branch remain uncovered. No existing tests or runtime modules were edited.
SHA-256 preservation checks cover all existing source, tests, examples and
all eleven historical numerical result files. The package version is unchanged
in both runtime and build metadata.

## Teaching documentation

The improved README links to `User_Guide.md` and gives installation and a
runnable walkthrough. The guide covers all five functions, deterministic
terms, fixed/automatic augmentation orders, all four candidate criteria,
separate covariance bandwidths, null hypotheses/rejection tails, p-values,
critical-value decisions, regression access and candidate diagnostics.

All ordinary Python blocks in the guide execute sequentially. All three
example scripts run across the six matrix rows. The new
`examples/teaching_walkthrough.py` can be adapted to Colab by uploading the
locally built wheel; the guide does not assume an unpublished GitHub/PyPI
package exists. Notebook installation/restart and Linux/Colab behavior remain
untested. The guide explicitly warns that failure to reject proves neither a
unit root nor stationarity, and does not frame these tests as automatic model
selection or a mechanical voting rule.

## Recommended version and remaining decisions

**Recommend `0.1.0rc1` as the first review/release candidate after instructor
approval**, followed by `0.1.0` only after the chosen release acceptance work.
The actual version remains **0.1.0.dev0**; no version file or tag was changed.

Before a broadly advertised first release:

1. Add and run Linux CI (including a clean wheel installation); test the intended
   Colab workflow if it is advertised. Resolve any observed differences without
   changing approved methodology or automatically widening tolerances.
2. Agree the supported platform/stack policy. The successful exact matrix is
   evidence for these combinations, not every version permitted by open ranges.
   Exact old numerical bounds cannot be promised on newer Python versions.
3. Review this guide, release notes and the proposed candidate version; approve
   any eventual version change, commit/tag and publication separately.
4. If the release promises regenerated cross-language evidence, rerun the frozen
   R study in a separately recorded environment. Otherwise clearly retain its
   existing saved-result provenance and earlier paper-access limitations.

No numerical or software compatibility failure blocks review of this stage
on the tested Mac environments. Linux/Colab validation is the principal remaining
platform acceptance gap. No new inference method or econometric scope is
proposed. Stop here for instructor review.

## Evidence and reproduction

See `validation/stage4/matrix.json` for commands, versions, import paths, clean
environment settings, test counts and coverage; per-environment subdirectories
contain frozen dependencies, XML test results, numerical replays and logs.
`validation/stage4/Reproduce.md` describes commands that avoid source shadowing.
`validation/validate_stage4.py` replays saved ADF/DF-GLS cases;
`validation/inspect_distributions.py` verifies local build contents.
`validation/stage4/preservation.json` records original-file hashes and the exact
allowed packaging change. No upstream source was modified.
