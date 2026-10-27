# Stage 4B — 0.1.0rc1 candidate review

27 October 2026. **Candidate prepared locally for instructor review; not
published.** No commit, push or tag was made. The production change is limited
to the authorized version update from 0.1.0.dev0 to **0.1.0rc1**; no test
statistic, formula, calibration, public API behavior, bandwidth convention or
historical numerical evidence was changed.

## Platform status

Linux and Google Colab execution are **unavailable and unverified** in this
session. The host is macOS 26.7.1 arm64, without Docker/Podman or a connected
Linux executor/Colab runtime. No platform pass is inferred from package
metadata, a pure-Python wheel, or local execution of validation code.

Prepared for later execution:

- `.github/workflows/platform-validation.yml`: manually triggered Linux matrix,
  current stack on Python 3.11–3.14 plus direct lower bounds on 3.11/3.12.
  It builds/inspects distributions, installs the wheel into an isolated
  environment, runs the complete tests, saved audits and teaching examples,
  and retains logs, versions and numerical differences. It has no repository
  write, tagging or publishing step. It has **not been dispatched**; using it
  remotely first requires separately authorized delivery of repository files.
- `validation/Platform_Validation.ipynb`: unexecuted Colab notebook. Upload the
  candidate wheel, source archive and SHA256SUMS, then run cells in order.
  It verifies checksums, creates an isolated environment within the runtime,
  runs the same checks and downloads the evidence. It avoids modifying Colab's
  base kernel; it does not establish arbitrary base-kernel dependency support.
- `validation/run_platform_validation.py`: shared installed-candidate runner,
  locally exercised. It copies tests/fixtures outside checkout source,
  removes PYTHONPATH, checks installed imports and package/backend versions,
  records dependency/BLAS metadata and retains failures without overwriting
  historical results. A local pass is not a Linux/Colab pass.

The workflow structure follows the official
[checkout](https://github.com/actions/checkout),
[setup-python](https://github.com/actions/setup-python) and
[upload-artifact](https://github.com/actions/upload-artifact) documentation.
Its YAML, shell blocks and matrix are checked locally; GitHub-hosted execution
remains pending. Notebook Python cells are compiled locally and their common
validation runner is exercised, but Colab upload/download integration is not.

## Exact local candidate matrix

All rows use fresh virtual environments with system site-packages disabled,
released ARCH **8.0.0**, pytest **9.1.1** and pytest-cov **7.1.0**, and install
the built candidate wheel. Each runs the full suite, three saved numerical
audit/replay jobs, all three teaching scripts and pip check.

| Platform | Python | Stack | Full suite | Numerical audits / examples / pip check |
|---|---|---|---:|---|
| macOS 26.7.1 arm64 | 3.11.18 | current | 1,164 pass | pass |
| macOS 26.7.1 arm64 | 3.12.16 | current | 1,164 pass | pass |
| macOS 26.7.1 arm64 | 3.13.18 | current | 1,164 pass | pass |
| macOS 26.7.1 arm64 | 3.14.10 | current | 1,164 pass | pass |
| macOS 26.7.1 arm64 | 3.11.18 | direct lower bounds | 1,164 pass | pass |
| macOS 26.7.1 arm64 | 3.12.16 | direct lower bounds | 1,164 pass | pass |

Current stack: NumPy **2.4.4**, SciPy **1.17.1**, Statsmodels **0.15.0**,
Pandas **3.1.0**. Lower stack: NumPy **1.26.4**, SciPy **1.11.4**,
Statsmodels **0.14.4**, Pandas **2.2.3**. These exact combinations are tested;
all permitted dependency combinations are not. The Stage 4 inability to install
NumPy 1.26.4 wheels on Python 3.13/3.14 remains a limitation, not a newly passed
case. Windows, x86_64 and other platforms remain untested.

A separate Python 3.12/current environment installs the candidate source
archive and runs the same complete checks. Candidate tests thus total seven
runs of the same 1,164 cases (8,148 passing executions), not 8,148 independent
experiments. Prior Stage 4 matrix/artifact checks remain preserved separately.

## Numerical outcomes and preservation

Each candidate environment confirms:

- Stage 2's historical **1,188 comparisons / 840 expected equalities**;
  independent reconstruction of **1,620 candidate rows**, **264 selector
  outputs** and **120 released-ARCH composition checks**.
- Stage 3B's **312 model pairs**, including **240 automatic selections**;
  all four distinct methods, fixed/automatic orders and effective samples.
- Stage 3C's **1,944 scalar comparisons / 1,800 expected equalities**;
  preserved PP/KPSS bandwidths, statistics, inference and package differences.
- All ARCH 8.0.0 calibrations, different nulls/tails and strict critical-value
  equality rules through the unchanged complete test suite.

Current-stack saved ARCH comparisons match exactly. Lower-stack maximum
saved-reference statistic differences remain **2.7924329515371937e−12** for
ADF/DF-GLS and **7.958078640513122e−12** for PP/KPSS; critical values match.
Fresh public-ARCH results match within each environment exactly. These are the
same investigated stack/BLAS rounding differences as Stage 4; no tolerance
was relaxed from the reviewed absolute 1e−8, no lag/bandwidth choice changed,
and no formula was altered. Stored urca PP and Statsmodels KPSS differences
remain expected, not erased. R was not rerun; prior provenance limitations
and the lack of size/power validation remain.

Preservation hashes cover the original Stage 2–4 numerical results and logs,
all existing tests/examples and runtime modules. Only `_version.py` changes
among runtime modules. Existing calibration source, selector and API code are
byte-identical. Historical reports remain historical; their development
version statements are not rewritten as candidate evidence.

## Packaging, installation and documentation

Runtime/build versions agree at **0.1.0rc1**. Runtime dependencies and minimum
Python declaration are unchanged. Wheel and source archive build successfully;
archive inspection verifies module bytes, MIT license/expression, ARCH pin,
README, complete Markdown documentation, examples, notebook and validation
materials. The sdist retains tests and reviewed fixtures; machine-specific
Stage 4/4B run evidence remains in the repository, avoiding recursive artifacts.
Twine metadata/description checks pass without uploading anything.

The archive inspector now compares metadata with pyproject.toml rather than
hard-coding the historical development version. This is a packaging-boundary
update, not a numerical test change. README/user-guide installation commands
name the actual candidate artifact. They are checked with clean candidate
wheel/source installs, pip check, public imports/all five functions and
installed-document access. Ordinary Python documentation examples execute;
notebook magic/Colab integration remains unexecuted.

Concise release notes: `docs/Release_Notes_0.1.0rc1.md`.
Local artifacts are `unitroot_tools-0.1.0rc1-py3-none-any.whl` and
`unitroot_tools-0.1.0rc1.tar.gz`; SHA256SUMS identifies their exact bytes. These
are local review artifacts, not an announcement of a public release.

## Readiness and remaining blockers

**Ready for instructor review as an unpublished candidate on the tested Mac
matrix.** No local numerical, API or packaging blocker remains. Linux and
Colab platform acceptance remains open: run the prepared files and inspect
all failures/differences before claiming those platforms passed or advertising
broad compatibility. Further dependency combinations and clean notebook-kernel
integration require their own evidence if promised.

No stable release, repository write, tag or publication is authorized by this
candidate preparation. Stop here for review and separate authorization.

## Evidence and reproduction

`validation/stage4b/` contains candidate distribution inspection, preservation
hashes, dependency lists, installed module paths, XML test/coverage reports,
numerical replays, build/check logs and per-environment `status.json` files.
The six wheel rows are labeled `311-current`, `312-current`, `313-current`,
`314-current`, `311-lower`, `312-lower`; `312-sdist` is the additional source
installation. Use the shared runner with an isolated candidate environment:

```sh
python validation/run_platform_validation.py \
  --source /path/to/extracted/source --out /tmp/unitroot-rc-results
```

Run from Python without `-O`, and retain outputs even if a check fails. The
original `validation/stage4/` evidence and Stage 2/3 result files are unchanged.
