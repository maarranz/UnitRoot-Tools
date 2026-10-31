# Experiment A — prepared Mac validation and Colab qualification

30 October 2026 · M1 Pro · preparation only

**Recommendation: GO for bounded Colab qualification only. NO-GO for full execution.** The original Linux preparation is integrated and passes Mac tests. One exact innovation-hash discrepancy in broader Linux–Mac replay remains unresolved, and hosted Colab numerical validation and fresh-runtime Drive recovery have not been performed. Extreme-scale raw-fit inaccuracy also remains; passing Mac package tests does not repair it or authorize restricted-domain acceptance.

## 1. Handoff integrity and integration

The ZIP was found at `/Users/maarranz/GitHub/UnitRoot-Tools/UnitRoot_Experiment_A_Linux_Handoff.zip`, not the expected external `handoff/` directory. It was the only untracked item. After filesystem access was granted, it was relocated unchanged to `/Users/maarranz/Documents/ChatGPT/UnitRoot-Tools/handoff/UnitRoot_Experiment_A_Linux_Handoff.zip`. SHA256 before and after: `625bff8da7ed57c9c088c7cd076cda45758bd825d4a27e44e459c71d24a41537`.

Extracted originals reside under `handoff/original-linux-handoff/`. All six original file checksums pass; outer and three nested ZIP integrity checks pass. The prepared source's 181 files match `SOURCE_SHA256SUMS`, the Linux evidence source manifest matches those same 181 files, and all 27 entries in the numerical diagnosis manifest pass. Evidence ZIP bytes are covered by the handoff checksum; that archive has a source manifest, not a separate complete per-evidence-file manifest.

A historical `validation/stage4b/SHA256SUMS` refers to release wheel/sdist builds omitted from the source snapshot. Its own bytes are verified by the source manifest, but those historical binaries are not verified. An attempted comparison against existing Mac `dist/` files found different hashes: those are separate build artifacts and cannot authenticate the unbundled historical Linux builds. This is recorded as an unavailable historical-artifact check, not a checksum failure of supplied source or diagnostic fixtures.

Git HEAD and recorded `origin/main` both equal `37fa6efab21a5f332c062c263007411cdb46623e`. After moving the known handoff archive, status was clean before integration. No fetch was needed beyond the successful prior review. `git apply --check` passed on the actual repository. The patch was first applied to an isolated baseline export: every resulting file exactly matched the original source snapshot. It was then applied to the Mac repository and all 181 source hashes matched again. Final `git diff --check` passes.

The nine intended changed/new files are the research README, benchmark, catalogue phase support, report, runner, two research test files, operational validator and Linux readiness document. Public source, package tests, dependencies, critical values, deterministic rules, exact DGP catalogue and old RNG key construction are unchanged. Changes remain uncommitted. No push, tag, pull, reset, clean, merge or rebase occurred.

## 2. Updated Mac validation

Existing interpreter: `/Users/maarranz/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/stage4b/env312-current/bin/python`. Actual Python 3.12.16; NumPy 2.4.4/Accelerate; SciPy 1.17.1; statsmodels 0.15.0; ARCH 8.0.0; pandas 3.1.0. PYTHONPATH was repo/src:repo, bytecode/user-site disabled, numerical thread controls set to 1, and OPENBLAS_CORETYPE unset. No dependencies installed/upgraded. Full source/build/thread provenance is in `environment.log` and replay reports.

| Check | Prepared Mac result |
|---|---|
| Complete package tests | **1,367 passed** |
| Research tests | **26 passed**, including main configuration guard with computation intercepted, worker/chunk/order/resume invariance, corruption/identity rejection, foreign/live/legacy lock safeguards and surviving-process lease refusal |
| Independent arithmetic/tree validator | 12 cases passed at unchanged 1e-8 tolerance |
| Existing Enders and deterministic validators | Both passed |
| Original Linux portable fixtures | 9 passed; exact innovation/input hashes and discrete records; 1,071 floating comparisons, max absolute difference 1.99e-13 |
| Actual SIGTERM and SIGKILL | Two-worker interrupted jobs resumed with one worker; both final 80-record sets exactly equal uninterrupted reference. Live recovery refused; SIGKILL stale lock archived using lease-aware recovery, not manually deleted |
| Small all-grid validation run | 100 DGPs × 4 replications = 400 datasets/8,400 decisions; zero warnings, failures or nesting violations |
| Independent all-grid audit | Input/innovation hashes checked against current Mac generation, 7,456 history nodes reconstructed; max error 5.65e-13; catalogue and old phase streams unchanged; all seven summary CSVs match original reporting code on this same new dry run |
| Larger Linux–Mac record comparison | **FAILED strict exact innovation parity:** one hash mismatch; all other fields meet exact/discrete or unchanged 1e-8 numerical requirements |

No pilot was rerun and no main workload executed. Portable replay uses the supplied isolated main-phase fixture keys only. Main guard tests construct configuration without computing main records.

## 3. Exact innovation discrepancy

The failed key is validation replicate **2**, DGP `f37e47c1c96adc9ce4005f8fc992f1eda1b065f853468e10d5aa651aa38dd033`: T=500, rho=.95, intercept=0. Across all 400 records and 43,591 float comparisons, this is the only exact-field mismatch:

- Linux innovation hash: `4381d297820c9e18827f3ac699f3007d8a56f99d72e602c488740963049ddfc5`
- Mac innovation hash: `44f1ffd79c4b73f67e60ea8fce1d338a1f881a971b083ec0ee0f53402c2379dc`
- Generated level-input hash is identical on both: `eea9b470f1a8688d1177d2bde39e1d2b2d0beb543f3ceb7b392d1a9adf4a3397`

All decisions, paths and other exact fields match; maximum numerical difference is 3.41e-13. The strict comparison still exits with failure, and no hash requirement or tolerance was weakened. `innovation-mismatch-Mac-only.npz` preserves this Mac's actual innovations/levels for investigation. The original Linux innovations for this key are not included as an array in the handoff, so the differing elements and cause cannot be isolated from the saved hashes alone. NumPy versions/platform builds differ; no causal attribution is established. Original Linux evidence was neither regenerated nor edited.

## 4. Exact numerical diagnosis

Saved `inputs.npz` and all nine `design_*.npz` were loaded unchanged. Designs constructed by the current public selector from those saved levels match the archived arrays exactly. The original Decimal solver was reused on these arrays at 80/120 digits; the 120-digit RSS references agree with the saved Linux references within 1e-70 relative error. The diagnostic compares Statsmodels raw SVD/pseudoinverse, raw/normalized SciPy SVD/QR alternatives, coefficients, RSS, t statistics, AIC/BIC rankings, fixed critical decisions and complete Enders routes. No production fitting code was replaced.

For scaled ct at ±1e-12, raw RSS relative excess is **8.64784e-6 on Mac**, versus **6.32593e-5 on Linux**. Mac normalized QR RSS agrees with the high-precision result at binary64 reporting precision. Mac raw versus normalized level t is approximately **-11.498044 versus -11.496230**; AIC differs by about **0.00102044**. The Mac package test budget of 1e-5 is met, but the unbalanced raw fit still loses accuracy. All checked fixed-ADF decisions, full AIC/BIC rank orders and Original/Conservative Enders paths remain identical. These vectors are not near their decision boundaries.

A read-only Mac re-audit of the **400 original saved Linux validation observations** reconstructs matching level-input hashes, 1,200 fits and 7,456 history nodes. Largest raw condition number: 1306.92; normalized: 127.87. Backend versus normalized QR RSS relative difference ≤6.66e-16, t difference ≤6.07e-13, IC difference ≤2.21e-13. Zero discrete disagreements. Minimum observed IC, fixed-ADF and reached-Enders margins are about 9.83e-4, 9.98e-4 and 5.66e-4.

Thus the discrepancy reproduces at extreme scale, while **no impact is detected in these checked ordinary-scale Experiment A observations**. This is finite evidence, not a guarantee over all 200,000 future datasets or a resolution of the separate innovation-byte discrepancy.

## 5. Prepared-runner Mac resources

Benchmarks ran sequentially after tests/operations. Nine representative cells ×20 rows per worker setting; full records identical across 1/2/4 workers.

| Workers | 180-row wall seconds | Datasets/second | Linear 200,000-row projection |
|---:|---:|---:|---:|
| 1 | 2.77 | 65.0 | 51.3 minutes |
| 2 | 2.55 | 70.5 | 47.3 minutes |
| 4 | 1.92 | 93.6 | 35.6 minutes |

The 400-row all-grid run took 3.54 s including interpreter startup (runner session 2.88 s), projecting 29.5 minutes. Plan approximately **30–55 computation minutes**, provisionally, plus separate reporting/auditing/snapshot time; a 60–90-minute operational allowance is prudent rather than a measured bound. Tiny samples, different chunks, sustained thermals and method branch mix can change this. No Colab speed estimate is asserted.

Four-worker process-peak sum was ~726 MiB, not simultaneous physical RAM. Checkpoint serialization/write overhead was ~0.067 summed worker-wall seconds over 4.47 task seconds (~1.5%) in that benchmark, and ~0.166/7.98 (~2.1%) in the all-grid run. These sums are not end-to-end wall time. All-grid evidence occupied 687,706 bytes, projecting ~344 MB; reserve at least 1 GiB per complete run plus separate snapshots/backups. Summary took ~1.79 s on 400 rows. Full-main memory/output/report costs remain unmeasured; allow at least 2 GiB available RAM and monitor. Drive transfer cost remains unmeasured.

## 6. Colab deliverable and remaining qualification

The old baseline notebook is preserved. New notebook: `colab-qualification/prepared-01/Colab_Prepared_Qualification.ipynb`. Its code cells compile. It uses the original handoff ZIP plus a checksummed support ZIP, verifies source and numerical fixtures, records numerical provenance, runs bounded tests/operations/exact replay/all-grid hash checks, and provides separate **stage** and **fresh-runtime restore** modes. It performs no installations and includes no main command. Missing dependencies or any strict gate failure stop promotion and preserve logs.

Computation stays on local Colab storage. Snapshot creation holds an exclusive OS lease, requires no coordinator lock, verifies shards, copies immutable archives to Drive, and records archive/internal hashes. Restore verifies hashes and current environment/source provenance before extraction/resume. Local helper tests passed: partial 2-row checkpoint restored and completed to 12 identical records, live writer lease rejected, corruption and changed environment rejected before extraction. This local POSIX test does **not** certify Drive/cloud durability.

No authenticated hosted Colab run, dependency setup, Drive mount, fresh-runtime restore or cloud benchmark was performed. The notebook requires those steps explicitly. It retains exact original Linux innovation expectations, so the known larger replay discrepancy cannot be silently waived. Before full execution: explain/resolve that discrepancy, review the raw-fit limitation without changing algorithms/tolerances in this assignment, pass hosted package/numerical/operational checks, demonstrate durable fresh-runtime restoration, and obtain separate full-run authorization. Colab remains preferred for Python convenience; the Mac remains appropriate for development/R/Julia.

Original Mac pilot evidence and all previous external reports were preserved; new outputs are confined to the existing external evidence tree and requested handoff directory. Final source hashes, Git status and checksummed evidence accompany this report.
