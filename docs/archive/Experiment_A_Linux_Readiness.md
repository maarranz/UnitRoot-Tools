# Experiment A — Manjaro preparation and readiness

29 October 2026. **Preparation complete; Linux main execution is not cleared.** Four unchanged package tests fail in the existing numerical environment. The research runner, independent arithmetic, dry run and operational checks pass. Resolve the package validation discrepancy and obtain explicit main-execution authorization before using the main commands below. No full main experiment or repeat pilot was executed.

## Repository and scientific scope

Repository: `/home/maarranz/Projects/UnitRoot-Tools`. Initial tracked and untracked status was clean on `main`. HEAD, locally recorded origin/main and the actual GitHub `refs/heads/main` were all `37fa6efab21a5f332c062c263007411cdb46623e`. Remote verification used `GIT_SSH_COMMAND='ssh -F /dev/null' git ls-remote origin refs/heads/main` because the system SSH configuration reported bad ownership/permissions for `/etc/ssh/ssh_config.d/20-systemd-ssh-proxy.conf`. No SSH configuration was changed. No repository AGENTS.md or applicable ancestor AGENTS.md was found.

The committed continuation notes, research protocol, validation report, pilot summary and research implementation governed this preparation. Historical documents retain their historical status statements. This report supersedes their “main disabled” description only for the prepared runner; it does not authorize execution.

The exact catalogue remains 100 laws, 153 aliases, 17 nulls and 83 alternatives; T counts 22/43/35 at T=100/250/500. Gaussian AR(1), y0=0, no burn-in, sigma=1, T+1 levels, p=0, both intercept-continuation I and fixed-mean M families remain unchanged. All seven methods and strict critical-value decisions at 1/5/10% are unchanged. Main is exactly 2,000 replications per law: 200,000 datasets and 4,200,000 decisions. The existing master seed and RNG construction are unchanged; `main` is a distinct phase in the hashed stream key. No pilot streams are reused as main streams.

The original pilot and its external Mac evidence have not been edited, copied over, or rerun. They were not available locally for a new cross-platform raw-pilot audit. The original pilot identity remains `d8265584c70f2b6f4e13dd98dfe3759c993a35a06ae6168d1a9bae5a69c130b7`.

## Changed files

| File | Change |
|---|---|
| `research/experiment_a/catalogue.py` | Adds the `main` phase name; leaves all laws, seed construction and old phase keys unchanged. |
| `research/experiment_a/runner.py` | Main opt-in, exact main count/full-catalogue guard, shared coordinator/worker leases and fail-closed stale-lock recovery, numerical-build/thread provenance, correct session start/finish UTC, task CPU timing, verified chunk iterator. |
| `research/experiment_a/report.py` | Reads one DGP's histories at a time; phase-aware summary metadata. Statistical formulas unchanged. |
| `research/experiment_a/benchmark.py` | Records summed task wall time and actual worker CPU time separately from elapsed run time. |
| `research/experiment_a/tests/test_research.py` | Checks main-stream separation in place of the obsolete unsupported-phase assertion. |
| `research/experiment_a/tests/test_main_preparation.py` | Main configuration without computation, live/foreign/legacy lock refusal, evidence-preserving recovery, surviving-process lease, ordering/chunk/worker/resume invariance. |
| `research/experiment_a/validate_operations.py` | Reproducible small real SIGTERM/SIGKILL/resume checks in validation phase. |
| `research/experiment_a/README.md` and this report | Current commands, limitations and recovery/transfer runbook. |

No changes to `src/`, existing package tests, dependency declarations, public API, handout, DGP parameters or method rules. No commit, push, tag or release. Changes remain local and uncommitted. A patch and complete source snapshot accompany this report; do not assume GitHub or the Mac has these changes.

## Linux environment and validation

Existing interpreter: `/home/maarranz/miniforge3/envs/unitroot_validation/bin/python`.

- Manjaro Linux 6.18.49-1, x86_64, glibc 2.44; Intel i7-4770, four physical cores/eight logical CPUs.
- Python 3.13.18; ARCH 8.0.0; NumPy 2.5.4; SciPy 1.18.1; Statsmodels 0.15.0; pandas 3.1.0.
- 15 GiB RAM, approximately 9 GiB available during inspection, no swap; approximately 566 GiB free on the evidence filesystem.
- One numerical thread per process for OpenBLAS, OMP, MKL, Accelerate and NumExpr; no OPENBLAS_CORETYPE override in accepted measurements. No dependency installs or upgrades.
- `environment.json`, `pip-freeze.txt`, CPU/memory/disk records and source checksums are in the evidence bundle. The manifest records package/research bytes, import path, versions, OS/architecture, build-configuration hash and numerical thread settings. This is a conservative resume identity, not a guarantee of cross-machine numerical equality.

| Check | Result |
|---|---|
| Existing package suite | **1,363 passed; 4 failed**. All four failures are the pre-existing ±1e−12 scale tests for AIC/BIC. |
| Final research suite | **26 passed**, including actual spawned workers, recovery and main-configuration guard. |
| Independent QR/tree fixtures | 12 vectors; passed at unchanged 1e−8 tolerance. Largest discrepancy 3.41e−13. |
| Existing independent Enders validator | Passed. Historical published rounding caveat retained. |
| Existing independent deterministic-selector validator | Passed; largest reported discrepancy about 8.53e−14. |
| Serial/parallel benchmark | 180 full records exactly equal at 1/2/4 workers. |
| Alternate chunk size | Same 180 records exactly equal for chunk sizes 5 and 20. |
| Reversed DGP order, changed worker count and resume | Exact full-record equality in research tests. |
| Real SIGTERM/resume | Partial checkpoints retained; final 80 records equal uninterrupted reference. |
| Real SIGKILL/recovery/resume | Stale lock archived; final 80 records equal uninterrupted reference. |
| Lock safeguards | Live coordinator and surviving worker refused; foreign host/owner, legacy lock and ambiguous liveness fail closed. Linux PID namespace/boot are checked. |
| Corrupt checkpoint and identity mismatch | Rejected by research tests; orphan files without commit markers are recomputed. |
| Small all-grid dry run | 100×4 validation-phase records; 8,400 decisions; zero failures/warnings/nesting violations. |
| Independent dry-run audit | All 400 input/innovation hashes and decisions checked; 7,456 history nodes reconstructed; maximum history error 6.05e−13. |
| Original-code preservation | Exact catalogue and old RNG outputs agree; all seven decompressed summary CSVs exactly match the original report implementation on the same dry run. |
| Portable fixtures | Nine saved vectors/records across T and validation/pilot/main keys replay successfully on Linux; Mac replay is supplied but has not been executed. |

The package discrepancy is concrete: the scaled/unscaled RSS ratio was **1.000063259311376**, outside the existing 1±1e−5 tolerance. It occurs on extreme ±1e−12 inputs, not the sigma=1 Experiment A dry run. Haswell, Sandybridge and Nehalem OpenBLAS overrides did not resolve the four tests. Other existing environments lacked required statistical packages. These observations do not justify loosening tests, modifying the public econometric code, or silently declaring a full green gate. The cause across numerical builds remains unresolved. Initial failed logs are preserved, including a development-only provenance warning subsequently fixed by requesting configuration dictionaries rather than formatted output.

## Measured resources and conditional runtime

Accepted measurements ran after package tests finished, with no simultaneous validation workload launched by this task. Benchmarks are small, so startup, CPU scheduling and branch composition affect extrapolation.

| Workers | 180-row elapsed wall time | Speedup | Sum of worker CPU time | Sum of process RSS peaks* |
|---:|---:|---:|---:|---:|
| 1 | 8.00 s | 1.00× | 7.85 s | 150 MiB |
| 2 | 5.57 s | 1.43× | 8.87 s | 468 MiB |
| 4 | 4.06 s | 1.97× | 10.98 s | 811 MiB |

*Sum of individual high-water marks, not simultaneous RAM use; shared pages and coordinator accumulation may inflate it. The largest individual process in the dry run was about 150 MiB.

The distinct all-grid dry run took **7.964 seconds end-to-end wall time** for 400 rows with four workers and chunks of four. Summed shard elapsed time was **22.001 seconds**; summed shard CPU time **21.809 seconds**. End-to-end descendant CPU, including imports/coordinator/startup, was **26.987 seconds**. None of these sums is elapsed wall runtime.

Direct all-grid wall scaling gives `7.964 × (200000/400) = 3982 seconds`, about **66.4 minutes**. The nine-cell four-worker benchmark gives about **75.3 minutes**; chunk-20 measurement (4.590 s for 180 rows) gives about **85 minutes**, with poor final-wave balance in only nine tasks. These are conditional projections, not a main run or confidence interval. Plan approximately **65–90 minutes computation, with a 90–120 minute operational window** for startup, disk/CPU variation and summary/audit work. The all-grid measurement weights every DGP equally; four replications cannot precisely capture all branch-cost variation.

Recommend **four workers, one numerical thread each**, and **chunk size 20** as an explicitly measured bounded-memory starting setting. Two workers provide a lower-load alternative; worker count can change on resume. Freeze chunk size at launch because it is part of checkpoint identity. Chunk 50 remains the existing default, but the commands below explicitly use 20.

The 400-row dry run occupied 687,242 bytes including manifest/checkpoint overhead. Naive scaling is about 344 MB; reserve at least **1 GiB per run plus separate copies/backups** rather than treating that small-sample estimate as exact. The dry summary used about **110 MiB maximum process RSS** and **4.20 seconds**. Summary now retains one DGP (2,000 main histories) at a time rather than all 200,000, with unchanged fixed-size summary tables. Full-main memory has not been measured; allow at least **2 GiB available RAM** as an engineering allowance and monitor actual use. The available host RAM exceeds the small-run measurements comfortably.

## Linux preflight and commands

These commands refer to the revised local source. The main commands are supplied for later review/authorization and **were not executed**. Do not run Python with `-O`: validation uses assertions.

```bash
export R=/home/maarranz/Projects/UnitRoot-Tools
export P=/home/maarranz/miniforge3/envs/unitroot_validation/bin/python
export E=/home/maarranz/Documents/Codex/Experiment_A_preflight_2026-10-29
export PYTHONPATH="$R/src:$R" PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
unset OPENBLAS_CORETYPE
cd "$R"
mkdir -p "$E"
"$P" -m pytest tests research/experiment_a/tests -p no:cacheprovider
# STOP if this fails. The recorded Linux environment currently has four failures.
"$P" -m research.experiment_a.validate --report "$E/independent.json"
"$P" validation/validate_enders.py --report "$E/enders.json"
"$P" validation/validate_deterministic.py --report "$E/ic.json"
# Each of the following output directories must be new.
"$P" -m research.experiment_a.validate_operations --output "$E/operations"
"$P" -m research.experiment_a.benchmark --output "$E/benchmark"
```

**Start, only after a resolved validation gate and explicit main authorization:**

```bash
export MAIN=/home/maarranz/Documents/Codex/Experiment_A_main_2026-10-29
# First launch only: refuse an existing output directory.
test ! -e "$MAIN" || { echo 'Choose a new main directory; do not overwrite evidence.'; exit 1; }
"$P" -m research.experiment_a.runner --output "$MAIN" \
  --phase main --reps 2000 --allow-main --workers 4 --chunk-size 20
```

**Resume the same run**, with the same exports, source bytes, environment and chunk size:

```bash
"$P" -m research.experiment_a.runner --output "$MAIN" \
  --phase main --reps 2000 --allow-main --workers 4 --chunk-size 20
```

This command verifies all committed shards and skips only verified completed ranges. A completed run is a no-op apart from another session record. You may use `--workers 2` or `--workers 1` on resume. Do not change phase, count, catalogue, chunk size, source, import path or numerical environment. Do not edit the manifest to bypass a mismatch.

For a graceful stop, from another terminal on the same host and PID namespace, inspect the lock and process owner/command, then send `kill -TERM <coordinator-pid>`. Wait for the current bounded worker wave and normal exit. A shell PID from inside an isolated task may not be meaningful in another PID namespace: never blindly signal a PID copied from a different namespace.

After a complete authorized run:

```bash
"$P" -m research.experiment_a.report --source "$MAIN" --output "${MAIN}_summary"
```

## Recovery runbook, independent of this task

1. Keep the original directory and all manifest/shard/lock bytes. Do not delete `coordinator.lock`, `coordinator.lease`, `.done.json` markers or orphan shards manually. The persistent lease file must never be replaced or unlinked: its inode coordinates writers.
2. If a process is still running, prefer graceful SIGTERM and wait. Do not “recover” a live job. A hard-killed coordinator can leave working children; their shared leases block recovery while they write.
3. On the **same host, user, boot and Linux PID namespace**, with the revised source, run:

   ```bash
   "$P" -m research.experiment_a.runner --output "$MAIN" --recover-lock
   ```

   Recovery exits without running any simulation. It requires a recognized lock format, matching host/UID/namespace, an absent owner PID and an exclusive OS lease. A live or reused PID, inaccessible PID, surviving writer, foreign/legacy lock or ambiguity causes refusal. On success it renames the exact lock into `recovery/` and writes a timestamped checksum/audit record. Queued workers carrying the old token cannot write after recovery.
4. Then use the identical resume command above. Checksummed committed chunks are preserved; uncommitted partial files are overwritten only by deterministic recomputation of their own task. If a committed checksum fails, stop and investigate or restore an independently verified backup copy; do not delete a marker merely to hide corruption.
5. If the host, boot or PID namespace changed, or the lock is from the old pilot runner, the recovery command intentionally refuses. Preserve a frozen copy and obtain independent evidence that the original coordinator and workers have stopped. The supplied tool does not override this ambiguity. Start a separately identified run on a validated environment if necessary; do not forcibly remove a foreign lock.

Use a local filesystem with working POSIX advisory locks. Network/shared-drive lock semantics and sudden power-loss durability have not been validated. Atomic renames protect normal process interruption; no new fsync/power-loss guarantee is claimed. Copy checkpoints only after the original job stops, or take a consistent filesystem snapshot, and verify copied file hashes before using them. Do not run two writable copies as one study.

Reproduce the operational checks with `python -m research.experiment_a.validate_operations --output NEW_DIRECTORY`; it uses only two approved cells, 40 validation replications each and two signal scenarios, never main or pilot. Research tests also check surviving-worker refusal and ambiguous locks. All scripts are preserved in the source snapshot.

## M1 Pro fallback and transfer

Mac repository: `/Users/maarranz/GitHub/UnitRoot-Tools`. No Mac access or new Mac execution is claimed. Preserve the original pilot evidence under `/Users/maarranz/Documents/Codex/2026-10-28/continue-the-unfinished-theoretical-audit-of/` and its original ZIPs unchanged.

Transfer `Experiment_A_Preparation.patch`, `Experiment_A_Preparation_Source.zip`, `Experiment_A_Linux_Evidence.zip`, this report and `SHA256SUMS` through the user's chosen channel. No changes were pushed. On the Mac, verify the transferred SHA-256 checksums, inspect the patch, and require a clean checkout at the base commit before applying:

```bash
cd /Users/maarranz/GitHub/UnitRoot-Tools
git status --short
git rev-parse HEAD
# Require clean status and HEAD 37fa6efab21a5f332c062c263007411cdb46623e.
# Set BUNDLE to the actual directory containing the transferred artifacts.
export BUNDLE="$HOME/Downloads/Experiment_A_Preparation"
(cd "$BUNDLE" && shasum -a 256 -c SHA256SUMS)
git apply --check "$BUNDLE/Experiment_A_Preparation.patch"
git apply "$BUNDLE/Experiment_A_Preparation.patch"
```

If status or base differs, do not reset, clean or overwrite the Mac repository. Review/merge the patch separately, or extract the source snapshot into a new directory and validate that snapshot. The source archive includes the base commit, all tracked source/docs and these local changes; it contains no `.git`, environments or pilot evidence. Its `SOURCE_SHA256SUMS` identifies the exact snapshot files. Do not assume a package with the same version string contains the required code.

The historically validated Mac interpreter was:

```bash
export P=/Users/maarranz/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/stage4b/env312-current/bin/python
export R=/Users/maarranz/GitHub/UnitRoot-Tools
export PYTHONPATH="$R/src:$R" PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
unset OPENBLAS_CORETYPE
```

Verify that interpreter still exists; its historical stack was Python 3.12.16, ARCH 8.0.0, NumPy 2.4.4, SciPy 1.17.1, Statsmodels 0.15.0 and pandas 3.1.0. Rerun all package/research tests, independent validators and the small benchmark/operations checks with new external Mac output paths. Do not upgrade the pilot environment in place. Extract the Linux evidence archive into a new directory, then replay the provided saved vectors and stream keys:

```bash
# EVIDENCE is the extracted evidence/ directory, containing portable-fixtures.*.
export EVIDENCE="$BUNDLE/linux-evidence/evidence"
"$P" "$EVIDENCE/replay_fixtures.py" --evidence "$EVIDENCE" \
  --report "$BUNDLE/mac-portable-replay.json"
```

The replay requires exact innovation/input hashes and discrete decisions/routes; floating statistics must agree within the unchanged 1e−8 absolute tolerance. It records maximum differences and the Mac provenance. A mismatch requires investigation, not editing the fixtures or tolerances. The fixtures include complete histories and nine small vectors, not a main-study execution.

After a green Mac gate and explicit main authorization, set `MAIN=/Users/maarranz/Documents/Codex/Experiment_A_main_M1_2026-10-29` and use the same main start/resume command with `--phase main --reps 2000 --allow-main --workers 4 --chunk-size 20`. Rerun the benchmark to choose workers on that machine. Historical M1 Pro pilot wall time was 93.09 s for 20,000 rows, giving a rough 15.5-minute tenfold projection; it is not a measurement of this revised main runner and should not replace the new benchmark.

**Checkpoint compatibility:** JSON/gzip/NPZ evidence is portable for inspection, audit and summarization. In-place Linux→Mac resume is intentionally rejected by platform, architecture, path, numerical-build and source identity. Matching dependency version strings alone do not make a Linux checkpoint resumable on the Mac. Transfer a stopped, checksummed copy for inspection; preserve its manifest unchanged. A Mac fallback starts a separately identified full run from replicate zero with the same fixed main stream keys, then becomes the chosen complete run. Do not mix or pool partial Linux and Mac runs as independent replications. The original pilot must be inspected using its frozen source/evidence, never re-labelled as main or updated to claim it ran with the new code.

## Remaining risks and review decision

- The four package test failures are the current Linux launch blocker. Investigating or preparing another numerical environment is follow-up work; no dependencies or econometric rules were changed here.
- Cross-platform numerical identity is unverified until the supplied replay and full Mac validation run there. Linux tests do not certify the Mac fallback.
- Runtime, storage and memory projections come from small runs; full-main execution and full-main summary memory were intentionally not measured.
- Recovery is conservative across hosts, reboots and PID namespaces. It refuses ambiguity; there is no unsafe force option. Filesystem power-loss and network-drive behavior remain uncertified.
- Pilot null rejection inflation remains a scientific finding, not a bug fixed by this preparation. Main retains raw methods, pointwise Monte Carlo intervals and paired comparisons. No size-adjusted power or uniform-size guarantee is introduced.

The deliverable is a reviewable prepared runner with transparent validation evidence, **not permission to launch**. No main output directory was created by the preparation commands. A unit test constructed a main manifest with computation replaced by an immediate exception; zero main replications were computed by the runner. Nine portable fixtures included three isolated main-phase keys solely to validate stream separation and replay.
