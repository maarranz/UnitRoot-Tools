# Experiment A research runner

This directory is separate from `src/unitroot_tools` and the public package API.
It implements only the approved 100-cell AR(1), zero-start, fixed-p=0 core in
`docs/Monte_Carlo_Research_Protocol.md`. All 153 parameter labels are retained as
aliases of 100 exact laws (22/43/35 cells at T=100/250/500; 17 null, 83 alternative).
Initial-condition/sign blocks, B, automatic lags and calibration are absent.
Main configuration is prepared, but execution requires separate authorization and
`--phase main --reps 2000 --allow-main`; see
[`Experiment_A_Linux_Readiness.md`](../../docs/archive/Experiment_A_Linux_Readiness.md)
for the Linux validation blocker, measured resources, commands and Mac fallback.

Run from a checkout with the validated numerical dependencies already available:

```sh
export PYTHONPATH="$PWD/src:$PWD"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m pytest research/experiment_a/tests -p no:cacheprovider
python -m research.experiment_a.validate --report /external/path/validation.json
python -m research.experiment_a.benchmark --output /external/path/benchmark
# Only after all validation passes:
python -m research.experiment_a.runner --output /external/path/pilot --workers 4
python -m research.experiment_a.report --source /external/path/pilot --output /external/path/summary
```

The fixed public master seed is hexadecimal `20261028b6a84fe0e1a0000000000001`.
See `catalogue.py` for byte-order/domain/phase rules. Each replicate is directly
addressable and makes one fixed-size PCG64DXSM standard-normal draw. Worker IDs,
output paths and task completion order do not enter the stream. Main uses its own phase key, exactly 2,000 replications and all 100 ordered cells.
Without `--allow-main` it is rejected. Validation/benchmark remain limited to
1–200 replications; pilot still requires exactly 200.
Defaults are all 100 cells and chunks of 50. `--cell-id` is an exact catalogue ID
filter for validation/benchmarks; `--max-tasks` provides a bounded checkpoint stop.

The runner calls unchanged public ADF, IC and Enders methods. It does not cache
or substitute regression fits. All alpha levels use critical-value decisions;
IC's downstream fixed-ADF rejection is an explicit research pipeline. Original
and default Conservative are compared on matched valid outputs. For the pilot's
first 20 replications per cell, permitted Conservative is checked against the
complete Original history; it is not another reported performance method.

Per-task gzip JSONL stores deterministic records and scalar Enders histories;
NPZ stores numeric decision/selection arrays with failures encoded -1, not as
nonrejection. JSON commit markers/checksums are written last. A manifest records
code/package bytes, versions, catalogue, completed ranges and timing. Resume
verifies every committed shard and rejects mismatched code/environment/identity;
missing markers mean a task must be recomputed, even if an orphan file exists.
SIGTERM requests a graceful stop after the current bounded worker wave. A lock
prevents concurrent coordinators. Shared OS leases protect the coordinator and
each writing worker. `--output PATH --recover-lock` archives a stale local lock
only after ownership/host/PID-namespace/liveness checks and an exclusive lease.
It refuses live workers, legacy/foreign locks and ambiguity. Never manually
remove the lock or the persistent `coordinator.lease` inode. See the recovery
runbook in the readiness report. `validate_operations` reproduces bounded real
SIGTERM and SIGKILL/resume exercises using validation streams only.

Cross-platform launch uses the same code/configuration and existing dependencies.
Replay saved inputs and record numerical margins rather than promising bitwise
floating-point identity across BLAS builds. Cross-environment resumption is
intentionally rejected to avoid silently mixing numerical builds in one run;
use a separately identified run and replay comparison instead. Linux validation
and limitations are recorded in the readiness report; Mac revalidation and
Colab/Kaggle execution are not implied. Frozen identity now also includes the
numerical build configuration and thread settings; old pilot manifests must
not be rewritten or resumed with revised source.

`report.py` exports compressed CSV rejection/selection rates, terminal branches,
full path counts, all 21 paired contrasts, score margins and timings, plus JSON
failure/warning lists. Individual rates use Wilson intervals; paired differences
use paired MCSE and conservative discordance bounds (Wilson for the nested
Original-minus-Conservative event). Failed methods retain attempted denominators
and failure-sensitivity bounds; no redraw, sample change or lag reduction occurs.
Pilot estimates are exploratory and do not establish selective inference or
uniform size. All output directories should be external to the repository.

Summary aggregation now loads one DGP at a time through verified chunks; the
statistical formulas and exported table values are unchanged. New sessions
record both actual start and finish UTC times. Historical pilot timestamps
remain untouched. Use local Linux/macOS storage with working advisory locks;
network filesystems and power-loss durability have not been certified.
