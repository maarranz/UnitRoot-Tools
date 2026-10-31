# Experiment A: Linux–Mac innovation discrepancy

30 October 2026 · numerical investigation only

**Finding:** the discrepancy is explained by one-ULP variation in the scalar `log1p` calculation used by NumPy's Gaussian-tail sampler. A controlled independent calculation, using identical raw random integers, reproduces **both complete recorded innovation hashes**. Binary64 AR(1) rounding then makes the complete level series identical. This is a real innovation-byte difference, not a hashing/serialization error.

**Gate recommendation: keep exact cross-platform innovation reproducibility blocked.** Explaining a mismatch does not satisfy the existing equality requirement. No correction, tolerance change or exception has been implemented. No main/pilot run, integration repeat, dependency change, repository edit or interaction with the separate Colab Work thread occurred.

## Provenance and scope

DGP `f37e47c1c96adc9ce4005f8fc992f1eda1b065f853468e10d5aa651aa38dd033`, validation replicate 2, T=500, rho=.95, intercept=0. Repository HEAD is `37fa6efab21a5f332c062c263007411cdb46623e` plus the already integrated preparation. Catalogue SHA256 is `92e4700e9ddcef9108c600f6a2edd18dc0316e8e91f0000d3ed1526197c42db3`; engine SHA256 is `283438c3060a116d675d1771f89f86a97a55b4106fe239f116fde5235eb5bd90`. Full initial/final repository inventories and interpreter/package/build/thread provenance are preserved alongside this report. All pre-existing repository file bytes are unchanged. An additional root `Colab_Prepared_Qualification.ipynb` appeared through other activity during this investigation; it was left untouched.

Mac: existing CPython 3.12.16, NumPy 2.4.4/Accelerate, SciPy 1.17.1, statsmodels .15.0, ARCH 8.0.0, pandas 3.1.0; arm64 macOS 26.7.1. Original Linux metadata: CPython 3.13.18, NumPy 2.5.4, SciPy 1.18.1, x86-64/glibc 2.44. No NumPy environment was installed: the controlled elementary-function experiment was more direct and smaller than a second package environment.

The original Linux innovation array/state/raw stream were not archived for this key. All new arrays are explicitly labelled diagnostic/Mac/inferred, never original Linux evidence. Original handoff hashes and source/diagnostic manifests remain intact.

## RNG and hash audit

The project uses `Generator(PCG64DXSM(SeedSequence(entropy)))`, then one `standard_normal(500, dtype=float64)` call. The canonical JSON key is `["experiment-a-rng-v1", DGP_ID, "validation", 2, "innovations"]`. Its SHA256 is decoded as eight little-endian uint32 words, appended to four least-significant-first uint32 words of master seed `20261028b6a84fe0e1a0000000000001`. No worker ID, scheduling order, output path or OS entropy enters this key.

`engine.evaluate` generates **one** innovation array, computes levels from it, hashes that same innovation array and passes those same levels to evaluation. `array_hash` converts to little-endian float64 and hashes C-order bytes; it does not round values. Shape is (500,), dtype float64, contiguous. Native, big-endian, Fortran-layout, strided-view and copied representations give the same normalized hash.

Five same-process and fifteen spawned-process replays across 1/2/4 workers reproduced identical innovations, initial/final states and raw integers; reversed task order and unrelated RNG calls had no effect. An independent Python integer implementation of PCG64DXSM matches all 1,024 saved raw uint64 values. The Gaussian interpreter consumes 520 uint64 values; its independently advanced final state matches NumPy exactly. Raw-stream SHA256: `c49a8ea2def21220acc4177e32f2ca6fef276917dda7880e488819aa48213e53`.

Actual historical Linux state equality is **not directly observed**. Identical archived key/source construction predicts the same seed; the same-bits controlled experiment fully reproduces its recorded output hash. Direct historical state/raw-bit confirmation remains available through the supplied one-key probe, without claiming it has run there.

NumPy guarantees PCG64DXSM's fixed-seed **integer** stream, while Generator's broader reproducibility conditions include the same build/environment/machine. The Gaussian guarantee cannot be inferred from the integer guarantee. See [PCG64DXSM documentation](https://numpy.org/doc/stable/reference/random/bit_generators/pcg64dxsm.html) and [compatibility policy](https://numpy.org/doc/stable/reference/random/compatibility.html).

## Mechanism and byte-level result

The tagged NumPy 2.4.4 and 2.5.4 scalar normal routines and Ziggurat tables are byte-identical; their PCG header is also identical. The normal routine uses the Ziggurat method, with scalar `log1p` in tail draws. NumPy maps `npy_log1p` to the C math-library function. See the saved, checksummed [normal sampler source](https://github.com/numpy/numpy/blob/v2.5.4/numpy/random/src/distributions/distributions.c) and [math header](https://github.com/numpy/numpy/blob/v2.5.4/numpy/_core/include/numpy/npy_math.h).

Only innovation index **236** (zero-based) takes the relevant tail path. Its raw-word offset is 244. The first uniform argument is `u = 0x1.4ae116fcad1ddp-1`:

| Quantity | Mac/system math | Generic glibc controlled calculation |
|---|---|---|
| `log1p(-u)` | `-0x1.0a066ffd3348dp+0` | `-0x1.0a066ffd3348ep+0` |
| innovation z[236] | `-0x1.f821c7c81f557p+1` | `-0x1.f821c7c81f558p+1` |
| innovation SHA256 | `44f1ffd79c4b73f67e60ea8fce1d338a1f881a971b083ec0ee0f53402c2379dc` | `4381d297820c9e18827f3ac699f3007d8a56f99d72e602c488740963049ddfc5` |

All other 499 innovations agree. Difference: -4.440892098500626e-16, exactly one innovation ULP. Tail acceptance and raw consumption are unchanged. Mac system math, NumPy's log1p ufunc and a 100-digit Decimal correctly-rounded log1p calculation give the Mac result. A locally compiled generic glibc log1p implementation gives the Linux recorded hash, with FP contraction both disabled and enabled. Fusing the final tail arithmetic does not change this example.

This confirms a sufficient platform-math mechanism without invoking a changed Gaussian algorithm. The controlled glibc source was downloaded from current master and retained by SHA256; the glibc-2.44 source tag was unavailable. It was compiled with Apple clang 21.0.0 using a minimal diagnostic wrapper, **not** the actual Linux binary. The original Linux symbol dispatch/compiler build was not directly inspected; do not equate the experiment with a historical Linux rerun.

At output time 237, both variants round to level `-0x1.0a30ebda06570p+3`. Its ULP is 1.7763568394002505e-15. The glibc-style exact sum lies precisely at the lower rounding midpoint and ties-to-even selects this level; the Mac sum also rounds to it. Since subsequent innovations match, all subsequent levels match. Common level SHA256: `eea9b470f1a8688d1177d2bde39e1d2b2d0beb543f3ceb7b392d1a9adf4a3397`. This is ordinary non-injective floating-point addition, not a separately generated or rounded hash representation.

## Competing explanations and downstream checks

| Explanation | Assessment |
|---|---|
| Platform-dependent Gaussian math | **Confirmed mechanism:** substituting generic glibc log1p alone reproduces the complete Linux hash from the same raw bits. Exact historical binary path remains unobserved. |
| NumPy-version algorithm change | **Rejected for the inspected normal routine/tables/PCG core.** Version numbers alone do not establish causation; installed binary/build differences are not fully isolated. |
| Hashing, byte order, dtype, memory layout | **Rejected** by explicit normalization/replay checks and independently reproduced hashes. |
| Different arrays used for innovation/level hashes | **Rejected** by engine source and byte-exact single-observation reconstruction. |
| Incorrect deterministic evidence record | **No error demonstrated; rejected as a necessary explanation.** An ordinary library calculation naturally reproduces the archived hash. This is not a blanket certification of every Linux record. |
| RNG-key mismatch | **Rejected as necessary explanation:** same key/raw stream reproduces both outputs. Direct historical Linux state remains unobserved. |
| Worker/order/global-RNG dependence | **Rejected** by narrow spawned-process and order tests. |

Independent recurrence reproduces the Mac levels exactly. Independent QR confirms p=0 regression inputs, coefficients, RSS, ADF statistics, AIC/BIC rankings, every reached Original/Conservative Enders step and all 21 final decisions. Both IC rank orders are n,c,ct. Largest independent arithmetic discrepancy is 1.14e-13; comparison with the original Linux observation differs by at most 2.28e-13 in checked floating fields. The known innovation hash is excluded **only from this downstream impact comparison**, not from the strict reproducibility gate, which remains failed.

## Implications and corrective action for review

This observation's actual simulated levels and econometric decisions are reproducible byte-for-byte/discretely; its innovation bytes are not. No defect in the unit-root testing algorithms or intended stochastic DGP is demonstrated. One case cannot establish full statistical equivalence, future bytewise portability or universal numerical safety. The separate extreme-scale OLS issue is unchanged. There is no basis here to invalidate the existing pilot; its distinct streams/evidence remain frozen, and its previous scientific cautions remain.

Original Linux evidence remains internally consistent and useful for its recorded environment. Preserve its original hashes and label the new calculations as independent diagnostic evidence. Do not replace its missing original innovation array with a candidate.

For review, propose either (a) an explicitly versioned deterministic Gaussian transform with controlled elementary functions/arithmetic and independent cross-platform validation, or (b) a canonical, checksummed reference innovation artifact with an approved generation/replay protocol. Merely pinning NumPy or dropping innovation hashes is insufficient. A new transform/RNG version can change exact streams; assess versioning and preserve old pilot streams rather than silently relabelling them. Any implementation/policy change needs authorization. A one-key original-Linux probe would close the remaining historical state/log1p attribution gap but would not itself fix portability.

**Strict gate remains NO-GO until an authorized reproducibility solution passes unchanged exact requirements.** Procedures, commands and small fixtures are saved here; no correction has been applied.
