# Experiment A — hosted Colab qualification result

30 October 2026 · **NO-GO for full Colab execution. Bounded hosted qualification failed.** This report supersedes the earlier access/prerequisite status reports. Full execution always requires separate user authorization.

## Completed and preserved

Authenticated hosted Google Compute Engine runtime, Drive mount, original archive/source/fixture integrity checks and the prepared bounded validation suite were executed. Handoff and support ZIP SHA256 matched unchanged notebook constants; six original handoff entries, 181 source entries, 181 Linux source-reference entries and 27 diagnostic entries passed hosted verification. Original static integrity review is preserved alongside this report.

Environment: CPython 3.13.18, Linux x86_64/glibc 2.39, NumPy 2.1.3 (OpenBLAS 0.3.27), SciPy 1.16.3, statsmodels 0.15.0, pandas 2.2.3, pytest 8.4.2. User-authorized arch 8.0.0 installation used the official PyPI binary wheel with --no-deps; other installed distribution versions were checked unchanged. Environment provenance, numerical build configuration, thread settings, source hashes and pip freeze are in the evidence bundle. No Mac dependency or repository changes occurred.

| Prepared check | Hosted result |
|---|---|
| Environment and dependency freeze | Commands passed |
| Package tests | Failed during collection: 563 items found, 2 import errors; no complete package pass |
| Research tests | 24 passed, 2 failed: order/chunk/resume invariance and serial/parallel exact parity |
| Independent validator | Exit 0 |
| Enders validator | Exit 1: missing installed unitroot-tools distribution metadata |
| Deterministic validator | Exit 0 |
| Original portable fixtures | Exit 1: validation_100/warnings list length mismatch |
| Exact numerical diagnosis | Exit 0; raw-fit numerical limitation remains |
| Operational validator | Exit 1 at resumed-versus-reference exact row equality |
| Validation all-grid | Complete: 100 DGPs × 4 replications = 400 datasets |
| Strict all-grid Linux comparison | Failed: warning-list differences in rows 0 and 4 |

The all-grid command took 27.44 seconds including interpreter startup; this is a small validation measurement, not a full-run forecast. Separate research and operational tests also performed their original bounded validation workloads. No pilot or full main simulation ran. Guard tests did not authorize main computation.

The all-grid comparator checked 43,591 floats with unchanged 1e-8 tolerance; max absolute difference was 8.178946409032051e-11. It reports exact innovation hashes, level-input hashes and discrete paths. Nevertheless all_passed=false because warning lists differ. No warning fields or exact requirements were filtered, modified or waived.

## Failure interpretation and numerical limitation

Package collection errors originate from Matplotlib calling deprecated pyparsing oneOf, promoted to errors under the unchanged test warning policy. Replay also records Matplotlib/pyparsing deprecations. These are demonstrated environment problems. The precise causes of both research equality failures and operational equality failure have not been isolated; the displayed row differences were truncated. They must not be dismissed as harmless warning differences without complete comparison. Enders failed at metadata collection, so its command cannot be called a passed validator.

At ±1e-12 scale with trend, the exact diagnostic reports raw RSS relative error 0.00025106143581843376 and raw-versus-QR t difference -0.002780858693755661. Normalized QR RSS error reports zero at binary64 output precision. Checked rank orders and Enders paths agree. Exit 0 does not certify numerical robustness or resolve the raw solver limitation.

## Strict unresolved Linux–Mac gate

The existing replicate-2 innovation discrepancy for DGP f37e47c1c96adc9ce4005f8fc992f1eda1b065f853468e10d5aa651aa38dd033 remains unresolved. Colab matching Linux innovations provides additional platform evidence, not a causal explanation or waiver. No original arrays, expectations, tolerances or algorithms were changed.

## Checkpoint boundary and evidence

The notebook archived all validation results and then failed its required-check assertion. The checkpoint stage cell was not executed. No runtime was terminated and no fresh-runtime restore or deterministic checkpoint continuation was demonstrated. Hosted snapshot corruption/locking/environment rejection and cloud durability therefore remain unqualified; operational validation itself failed. Drive read-back of a log bundle is not checkpoint recovery evidence. Runtime remained connected at this stopping boundary.

Drive evidence: My Drive/UnitRoot-Tools/qualification/prepared-qualification-01/validation-1793383951472643089.zip, with SHA256 sidecar. SHA256: 5f919ba9d253dc74cc63a71194c3c6997ff72bdf6694f04fabbc85697152946e. Downloaded local copy Hosted_Validation_Evidence.zip matches that SHA256 and passes ZIP CRC verification. It contains full logs, provenance, diagnostic output, validation shards/manifests and operational artifacts. Installation evidence is a separate arch-setup ZIP in the same Drive folder.

All local outputs remain outside the Git repository and M1 Pro evidence tree. Original scientific cells/source were preserved; supplemental hosted cells only installed the explicitly authorized dependency, inspected prerequisites/results and downloaded evidence.

## Work required before reconsideration

Prepare and review a compatible hosted dependency environment and genuine package installation metadata, preserving scientific source and strict warning/hash rules. Diagnose complete differing records for research/operations, then rerun only checks justified by those fixes under a separately identified environment. Pass all strict gates before the prepared stage/delete/fresh-runtime-restore exercise. Resolve the Linux–Mac innovation mismatch and separately review the raw-fit limitation. Do not launch full execution without explicit authorization.
