# Experiment A implementation and independent validation report

28 October 2026 · UnitRoot-Tools · Mac M1 Pro · research-only runner

## Outcome

The approved 100-cell A-core runner is implemented separately under `research/`.
All acceptance checks passed **before the pilot launch**, and the completed pilot
contains exactly 20,000 unique data keys and 420,000 decisions. Independent
post-run audit regenerated every input/innovation hash and checked all decisions,
IC formulas/winners/gaps, 120,000 Enders histories and exported counts. There are
no recorded numerical failures/warnings or nesting/equivalence violations.

Package source, APIs, metadata and dependencies are unchanged. No main run,
Experiment B, automatic lags, size adjustment, additional trend study, commit,
push or tag occurred. HEAD remains `b6a84fec7ca28a1802cdd722c163f624d16dca24`.

## Files added

```text
research/__init__.py
research/experiment_a/README.md
research/experiment_a/__init__.py
research/experiment_a/benchmark.py
research/experiment_a/catalogue.py
research/experiment_a/engine.py
research/experiment_a/report.py
research/experiment_a/runner.py
research/experiment_a/tests/test_reporting.py
research/experiment_a/tests/test_research.py
research/experiment_a/validate.py
```

The root research namespace is outside the installed `src/unitroot_tools` package.
No existing source/test/doc file was edited. Before/after SHA-256 checks preserve
all 727 existing non-Git/non-bytecode files, including unrelated
untracked reports and the protocol. The public API remains the scientific
reference; this implementation does not insert a fit cache or replace the OLS
engine.

## Scope and reproducibility

Exactly 100 laws (17 null/83 alternative; T counts 22/43/35) represent 153 original
aliases. Rational roots and square-free radical drift descriptors establish exact
law identity before float64 resolution. A-law aliases do not affect streams.
The published master seed is `20261028b6a84fe0e1a0000000000001`; PCG64DXSM/SeedSequence
uses fixed little-endian words and a SHA-256 domain containing RNG version,
canonical cell ID, phase, replicate and component. One T-sized normal draw and
zero-start recursion produce T+1 levels. All final fits use p=0 and T response
rows `[1:T+1]`. Main phase/counts above 200 are rejected; pilot requires exactly
200. Validation/benchmark/pilot streams are distinct.

Public calls preserve the existing rejection rules and critical values. IC's
selected fixed-specification rejection is the protocol's research pipeline,
without adding any selector-level p-value. Complete scalar histories and failure
records are stored; backend objects are not serialized. All seven procedures
receive common data per replication.

Checkpoints use deterministic gzip JSONL and portable numeric NPZ, with SHA-256
commit markers written last. Resume verifies frozen code/environment identity
and every committed shard; unfinished orphan files are not completion evidence.
A lock prevents concurrent coordinators. Operational worker count does not enter
random streams or scientific records.

## Tests and independent numerical checks

| Check | Result |
|---|---|
| Existing complete package suite | 1,367 passed; 0 failed/errors/skipped |
| Research acceptance suite | 22 passed; 0 failed/errors/skipped |
| Independent QR/tree reconstruction | 12 fixed vectors across T=100/250/500, all three α; passed |
| Existing independent Enders validator | 18 regression cases; report exactly matches earlier Mac baseline |
| Existing independent IC validator | 18 distinct fits /12 selections; passed |
| Real SIGTERM/resume | Stopped after 4 checkpointed rows; final 16 rows exactly equal uninterrupted reference |
| Serial/parallel benchmark | All 180 records identical at workers 1,2,4 |
| Full-pilot streaming audit | 20,000 unique keys, 420,000 decisions, 120,000 histories checked |
| Permitted Conservative versus Original | 6,000 complete-history comparisons passed, not another performance method |
| Conservative subset invariant | 60,000 matched α comparisons; 0 reverse discordances |
| Complete record counts and summaries | All 100×200 keys; all exported rejection/selection/paired counts agree |

One additional reporting test initially failed: the Wilson lower endpoint at
zero events rounded to 1.7347e−18 instead of exactly zero. It was corrected by
explicit endpoint handling; the same test then passed. No existing package test
or tolerance was weakened. The initial failure log is retained. The first 21
runner tests and numerical/routing checks passed before the preliminary benchmark;
the accepted benchmark and pilot followed the complete 22-test acceptance gate.

The independent validator builds response/design rows explicitly, solves QR,
computes covariance/standard errors and phi through both RSS and Wald formulas,
sums Gaussian densities and reconstructs routing independently. Tau calibration
is checked against a separate call to released public ARCH; this verifies correct
calibration use, not independent derivation of the response surfaces. Phi thresholds
were separately transcribed from Enders Table B; Student thresholds use independently
computed df. The frozen absolute tolerance was 1e−8 and was not relaxed.

| New independent quantity | Maximum absolute discrepancy |
|---|---:|
| coefficient | 4.496403249731884e-15 |
| ssr | 5.684341886080801e-14 |
| statistic | 6.183942247162122e-14 |
| phi | 1.13686837721616e-13 |
| auxiliary | 1.06581410364015e-14 |
| aic | 1.172395514004165e-13 |
| bic | 1.172395514004165e-13 |
| log_likelihood | 1.13686837721616e-13 |
| critical_value | 0 |

The existing Enders validator's largest discrepancy is about 3.12e−12, also below
its unchanged 1e−8 tolerance. Its known published A2.3 intercept-t rounding mismatch
(2.646058967 versus 2.64) remains documented and unchanged; it is not a new runner
failure. Fixed curvature vectors used only for numerical branch validation are
not additional simulation DGPs or a trend experiment.

The post-run audit independently rebuilt the seed layout and recurrence, checked
all 153 catalogue aliases, verified strict rejection inequalities and sample/df
metadata at every recorded node, and checked phi/table and Student tail/df mappings.
The stored IC arithmetic audit had maximum error 0.
It checked every actual history; controlled edge/equality/rank cases also remain
covered by the unchanged package suite. Pilot branch rarity is not treated as
proof that an unobserved branch is impossible.

## Portable outputs and commands

The evidence archive includes the exact research source, all raw checkpoint data,
NPZ masks/decisions, manifests/checksums, catalogue/aliases, benchmark, validation,
summary tables and audit/report scripts. The completed run identity is
`d8265584c70f2b6f4e13dd98dfe3759c993a35a06ae6168d1a9bae5a69c130b7`. `summary.json` reports all row/decision counts.
Raw failed-method masks would be −1 rather than nonrejection; no failures occurred.

Interpreter used:
`/Users/maarranz/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/stage4b/env312-current/bin/python`.
Python 3.12.16, macOS 26.7.1 arm64; dependencies:
`{"arch": "8.0.0", "numpy": "2.4.4", "pandas": "3.1.0", "scipy": "1.17.1", "statsmodels": "0.15.0"}`.
Numerical threads were one per worker. No dependencies were installed or upgraded.

From the repository, with P set to that interpreter, R to the repository and W to
an external evidence directory:

```sh
export PYTHONPATH="$R/src:$R" PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
"$P" -m pytest tests -p no:cacheprovider
"$P" -m pytest research/experiment_a/tests -p no:cacheprovider
"$P" -m research.experiment_a.validate --report "$W/independent-validation.json"
"$P" validation/validate_enders.py --report "$W/enders-validation.json"
"$P" validation/validate_deterministic.py --report "$W/ic-validation.json"
"$P" -m research.experiment_a.benchmark --output "$W/benchmark"
# Only following acceptance:
"$P" -m research.experiment_a.runner --output "$W/pilot" --phase pilot --reps 200 --workers 4 --chunk-size 50
"$P" -m research.experiment_a.report --source "$W/pilot" --output "$W/summary"
```

`pilot-command.json`, the validation gate and interruption script record exact
executed arguments and environment overrides. Use a fresh output directory or
resume the exact frozen identity; do not overwrite another study.

## Remaining review items before main

1. Adaptive null rejection inflation is substantive pilot evidence. Preserve raw
   rules, distinguish raw from size-adjusted power, and review interpretation
   before authorizing R=2,000. No rule was tuned to pilot outcomes.
2. R=200 and pointwise intervals do not establish uniform size or precise power
   rankings. Secondary 1% results and rare branches are especially uncertain.
3. Portability was designed/tested through spawn processes and portable files,
   but execution occurred on this Mac only. Linux/Colab/Kaggle replay acceptance
   remains necessary; cross-environment resume is intentionally rejected rather
   than silently mixing numerical builds.
4. SIGTERM and bounded interruption/resume were tested. After a hard kill, the
   coordinator lock requires confirming all old workers stopped before manual
   recovery; surviving-worker crash recovery is not automatic.
5. The current summary loads histories into memory. Review streaming aggregation
   for a larger main/full study; measured resource figures are in the runtime report.
6. Minor bookkeeping limitation: the current manifest session field `started_utc`
   is populated at completion, while elapsed time uses a correct monotonic clock.
   Treat it as a completion timestamp; the acceptance gate supplies a prior UTC
   launch record. Rename/add timestamps in a separately versioned pre-main runner
   revision, retaining this executed snapshot and its hashes. This affects neither
   streams, estimates nor measured durations. The validated/pilot source was not
   changed after freezing the gate.

No remaining statistical-arithmetic, routing or numerical failure was found.
Promotion requires review and fresh authorization; the runner cannot execute the
main experiment as delivered.

## Final repository status

```text
## main...origin/main
?? docs/Deterministic_AIC_BIC_Correction_Commit_Report_2026-10-28.md
?? docs/Deterministic_AIC_BIC_Final_Review_2026-10-28.md
?? docs/Deterministic_AIC_BIC_Implementation_Validation_Report_2026-10-28.md
?? docs/Enders_Mac_M1_Pro_Validation_Report_2026-10-28.md
?? docs/Hamilton_Enders_Theoretical_Review.md
?? docs/Monte_Carlo_Research_Protocol.md
?? docs/Theoretical_Audit_Sequential_Strategies.md
?? research/
```

No commit, push, tag, API/version/dependency change or public econometric edit.
