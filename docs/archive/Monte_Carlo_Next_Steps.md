# Monte Carlo continuation notes

28 October 2026 · UnitRoot-Tools · end-of-day preservation

## Checkpoint and current state

- Repository: `/Users/maarranz/GitHub/UnitRoot-Tools`.
- Branch/remote: `main`, `origin` = `git@github.com:maarranz/UnitRoot-Tools.git`.
- Package version: **0.1.0rc1**, unchanged.
- Git commit **before this preservation operation**:
  `b6a84fec7ca28a1802cdd722c163f624d16dca24`.
- The preservation commit is the commit introducing this document. Identify its
  hash with `git log -1 -- docs/Monte_Carlo_Next_Steps.md`; the pre-operation hash
  above remains the package checkpoint used for the pilot.

The deterministic AIC/BIC selector, its corrected rank handling and the existing
Enders algorithms are already in the pre-operation checkpoint. The older final
review's rank blocker was subsequently corrected; that report is historical.
Today's additions preserve the research runner, theoretical reviews and compact
pilot summaries. No public econometric algorithm, dependency or API is changed
by this preservation operation. No main experiment has been launched.

## Experiment A: approved and completed scope

The governing design is [Monte_Carlo_Research_Protocol.md](../Monte_Carlo_Research_Protocol.md),
with **only its 100-cell A core** implemented and piloted:

- Gaussian AR(1) DGPs; zero-start initialization; no burn-in; innovation SD 1.
- Effective ADF sample sizes **T=100,250,500**; **fixed p=0**; T+1 input levels.
- Zero, local-to-zero and fixed drift nulls.
- Both alternative families: I retains the regression intercept; M fixes the
  equilibrium mean through a=(1−rho)m, m=0,1,5. Their interpretation differs.
- Exact parameter deduplication: **100 laws, 153 aliases; 17 null and 83 alternative
  cells**, with 22/43/35 cells at the three sample sizes.
- Seven procedures: prespecified ADF n/c/ct, Original Enders, default Conservative
  Enders, global AIC→fixed-specification ADF, and BIC→fixed-specification ADF.
- Nominal levels 1%,5%,10%; primary level 5%.
- Pilot: **200 replications per DGP, 20,000 cell-replications, 420,000 decisions**.

Permitted Conservative is an equivalence diagnostic, not an eighth method.
There is no separate verified Hamilton/Enders algorithm. No size-adjusted
calibration, Experiment B, automatic lag selection, sign/initial-condition study,
or additional trend experiment has been implemented or run.

## Validation and principal empirical findings

See [Experiment_A_Validation_Report.md](../Experiment_A_Validation_Report.md) and
[Experiment_A_Pilot_Summary.md](../Experiment_A_Pilot_Summary.md).

- Existing package suite: **1,367 passing tests**; research suite: **22 passing**.
- Independent QR/decision reconstruction agreed within about **1.2e−13** against
  a fixed 1e−8 tolerance. The independent Enders report matches the historical
  Mac baseline; its largest numerical discrepancy remains about 3.12e−12.
- Serial/parallel records match. Actual SIGTERM/resume reproduces an uninterrupted
  reference. The complete pilot audit regenerated all input/innovation hashes
  and checked all 420,000 decisions and 120,000 scalar Enders histories.
- **Zero numerical failures, zero recorded warnings, zero nesting violations**;
  **6,000 permitted-Conservative/Original history checks** passed.
- At nominal 5%, observed null rejection frequencies across the 17 null laws were
  **8.5–15.5% for Original Enders**, **3–11.5% for default Conservative**, and
  **3–14.5% for each IC→ADF pipeline**. Prespecified ct ranged **3–7%**.
- These are exploratory per-cell pilot estimates, not uniform-size guarantees
  or a reason to tune the DGP grid/rules. Higher raw power cannot establish
  superior size-adjusted power when null rejection rates differ.
- Default Conservative's rejection-set subset relation to Original held on all
  matched outputs. That property does not establish nominal size control.

R=200 is deliberately exploratory: at probability .05, expected MCSE is about
1.54 percentage points. The summary tables contain pointwise Wilson intervals
and paired uncertainties. Do not pool distinct laws/aliases/levels as independent
replications or interpret retained ct/c as proof of true deterministic terms.

## Files preserved in Git and how to use them

Research implementation, tests and usage instructions are under
[`research/experiment_a`](../../research/experiment_a/README.md), outside the installed
public package. It calls unchanged public methods without substituting a cached
regression engine. The current runner deliberately rejects `main` phase and
replication counts above 200; pilot requires exactly 200. Do not bypass these
safeguards without separate authorization and a reviewed runner revision.

Compact results in this directory:

| File | Contents |
|---|---|
| `rejection_rates.csv.gz` | 2,100 cell/method/alpha rows; filter `is_null=True` for null frequencies; includes counts, MCSE and intervals |
| `paired_differences.csv.gz` | 6,300 paired comparisons; difference is method_a minus method_b; primary 5% contrasts flagged |
| `cell_timings.csv.gz` | 100 measured cell timing summaries |
| `summary.json` | Completion counts, failures, warnings and run identity |

All files were checked against the original external evidence and decompressed
successfully before preservation. Additional selection, branch, path and IC-margin
tables are in the external archive below. Raw checkpoints and large archives are
intentionally not committed.

## Original external evidence: preserve unchanged

The original evidence directory on this Mac is:

```text
/Users/maarranz/Documents/Codex/2026-10-28/continue-the-unfinished-theoretical-audit-of/work/experiment-a-pilot/
```

Important relative paths there:

- `pilot/manifest.json`, `pilot/catalogue.json`, `pilot/shards/`: frozen execution
  identity, exact aliases, checksummed raw JSONL/NPZ checkpoints and commit markers.
- `summary/`: complete portable summaries, including selection/branches/margins.
- `benchmark-accepted/benchmark.json`: accepted 1/2/4-worker measurements.
- `validation-gate.json`, `independent-validation.json`,
  `independent-pilot-audit.json`, `interruption-validation.json`: validation evidence.
- `environment.json`, `pilot-command.json`: runtime and exact launch settings.

User-facing external deliverables are in:

```text
/Users/maarranz/Documents/Codex/2026-10-28/continue-the-unfinished-theoretical-audit-of/outputs/
```

Keep these original files unchanged:

- `Experiment_A_Pilot_Evidence.zip` — frozen research source, raw pilot checkpoints,
  summaries, benchmark and validation evidence (about 16 MB).
- `Experiment_A_Pilot_Summary_ChatGPT_Upload.zip` — compact upload bundle, with
  decompressed CSVs, JSON timing records and reports (430,483 bytes).
- `Experiment_A_Runtime_Estimates.md` — measured timings, memory and explicitly
  conditional workload projections.

Pilot run identity:
`d8265584c70f2b6f4e13dd98dfe3759c993a35a06ae6168d1a9bae5a69c130b7`.
Public master seed: `20261028b6a84fe0e1a0000000000001`.
Streams are independently addressable by law, phase, replicate and component;
workers/output paths do not change data. Main must use a distinct approved phase,
not reuse pilot streams. Version number alone is insufficient to identify the
package: older same-version wheels predate the features; verify commit and hashes.

## Measured resources and known pre-main issues

On native M1 Pro/macOS, four workers completed the pilot in **93.09 seconds**.
Direct A-only scaling suggests roughly **15.5 minutes** for the 100-cell main
experiment, plus aggregation/I/O and platform variation; this is a projection,
not an executed main run. The full proposed 300-cell study is not implemented;
its runtime estimates assume unmeasured B costs and are not cloud guarantees.

The summary process used about 736 MiB RSS at 20,000 records and materializes
histories in memory. Review streaming aggregation before larger runs. The
checkpoint computation is chunk-bounded. No GPU benefit is assumed.

Execution/replay occurred on this Mac only. Linux/Colab/Kaggle reproducibility
remains to be verified. Cross-environment resume is intentionally rejected rather
than silently mixing numerical builds. A hard kill can leave a coordinator lock;
confirm old workers stopped before manual recovery. Graceful SIGTERM/resume was
validated.

Minor recorded bookkeeping issue: the frozen pilot manifest's session field
`started_utc` contains its completion timestamp; the monotonic elapsed duration
is correct. Treat that field accordingly. Rename/add timestamps in a separately
reviewed runner revision, preserving the original executed snapshot and hashes.
Do not edit old manifests to pretend a different source version produced them.

## Next task: Google Colab readiness, then authorization

**Prepare a Google Colab notebook; do not launch main yet.** The notebook should
be a thin launcher for the same reviewed research code and statistical rules.

1. Pin the package/research checkpoint and compatible dependency versions; verify
   import provenance and numerical fixtures. The Mac baseline used Python 3.12.16,
   ARCH 8.0.0, NumPy 2.4.4, SciPy 1.17.1, Statsmodels 0.15.0 and pandas 3.1.0.
   Do not assume Colab defaults match or upgrade this repository environment.
2. Replay designated existing input/stream keys; compare data hashes, statistics,
   decisions and complete Enders histories. Record numerical-build differences
   and margins instead of promising bitwise identity across platforms.
3. Use validation/benchmark-only runs to verify serial/parallel behavior,
   interruption/resume, checksum rejection and Google Drive checkpoint persistence
   across a runtime restart. Confirm durable storage paths and recovery procedure.
4. Measure Colab runtime, memory, output size and worker efficiency. Review the
   summary-memory and timestamp issues without changing econometric algorithms
   or simulation rules. Keep historical evidence unchanged.
5. Present the notebook, reproducibility/persistence checks and cost estimate.
   **Request explicit authorization before launching the 2,000-replication main
   experiment** or extending the runner to enable it.

Planned main A scope, subject to that authorization: 100 unchanged laws ×2,000
replications = **200,000 cell-replications and 4.2 million method/alpha decisions**,
with distinct main streams. Do not add B, automatic lags, size adjustment, trends
or other sensitivity blocks by implication. Their specifications/budgets require
separate review. Do not tune methods or remove cells based on pilot outcomes.
