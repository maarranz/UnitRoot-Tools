# Experiment A core: exploratory pilot summary

28 October 2026 · UnitRoot-Tools · commit b6a84fe · Gaussian AR(1), p=0

## Completion and interpretation

Completed **100 distinct configurations × 200 replications = 20,000 cell-replications**,
covering 17 unit-root null and 83 alternative cells and all 153 original aliases.
The seven procedures at 1%,5%,10% produced **420,000 decisions**. There were
**zero numerical failures, zero recorded warnings, zero Conservative/Original
nesting violations, and zero mismatches in 6,000 permitted-Conservative equivalence
checks**. No 2,000-replication main experiment or other experiment was run.

The pilot shows appreciable null rejection inflation for several adaptive pipelines
in some cells. This is not a software failure to repair by changing the rules:
independent arithmetic/routing validation and the all-record audit agree. It is
an exploratory finding about the defined complete procedures. Default Conservative
rejects less often than Original but is not thereby calibrated. IC→ADF means the
explicit downstream selected fixed-specification test; it is not an IC p-value.

All tables use equal denominators of 200 per cell. Individual confidence intervals
are Wilson 95% intervals; they are pointwise Monte Carlo uncertainty intervals,
not simultaneous tests of all cells or confidence statements about uniform size.
At probability .05 the anticipated MCSE at R=200 is .0154 (1.54 percentage points).
At .01 only about two rejections are expected under nominal behavior; secondary-level
comparisons are especially imprecise. Full per-cell results at every level are in
the evidence archive's compressed CSV tables.

## 1. Exact DGPs and scientific scope

Null: y_t=y_(t−1)+a+z_t. Family I: y_t=ρ y_(t−1)+a+z_t, retaining the specified
intercept index. Family M: a=(1−ρ)m for m∈{0,1,5}. Innovations are iid standard
normal, y0=0, no burn-in; T=100,250,500 effective response observations and T+1
supplied levels. Stable-root alternatives are conditionally zero-start processes,
not claimed to be strictly stationary from observation zero.

Exact deduplication yields T-specific cell counts 22,43,35 and null counts 5,6,6.
All aliases are retained, including the local/fixed drift coincidence at T=100,
root coincidences, and I/M overlap. One finite-sample law is never counted twice
within an estimated rate. Family-specific views can legitimately show the same
cell under two labels; they must not be pooled as independent replicates.

The pilot does not examine true nonzero difference-equation trends, non-Gaussian
errors, breaks, short-run AR/MA errors, automatic lags, alternative initial states,
negative-drift blocks or size-adjusted calibration. In this core the difference-
equation trend coefficient is zero; retaining ct does not prove a genuine trend.

## 2. Null rejection frequencies at nominal 5%

These are rejection frequencies under the declared stochastic-unit-root DGPs.
For n/c tests outside their maintained deterministic calibration, annotate the
comparison as misspecification/conditional-inference sensitivity rather than a
universal nominal-size check. Drifted levels do not imply a difference-equation trend.

| T | a/σ | ADF n | ADF c | ADF ct | Original | Conservative | AIC→ADF | BIC→ADF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 0 | 5.5% | 6.0% | 6.0% | 14.5% | 9.5% | 12.0% | 12.5% |
| 100 | 0.05 | 6.5% | 4.5% | 6.5% | 14.5% | 9.0% | 12.5% | 14.0% |
| 100 | 0.1 | 4.0% | 1.5% | 5.0% | 11.0% | 6.5% | 9.0% | 9.5% |
| 100 | 0.2 | 1.0% | 4.0% | 4.5% | 13.5% | 7.5% | 8.0% | 8.0% |
| 100 | 0.5 | 0.0% | 0.0% | 3.0% | 13.5% | 3.0% | 3.0% | 3.0% |
| 250 | 0 | 6.5% | 3.0% | 3.5% | 12.5% | 6.0% | 10.5% | 11.5% |
| 250 | 0.0316228 | 7.0% | 7.0% | 4.0% | 15.5% | 10.0% | 13.5% | 14.5% |
| 250 | 0.0632456 | 3.5% | 5.0% | 7.0% | 15.0% | 11.5% | 14.5% | 14.5% |
| 250 | 0.1 | 1.5% | 1.5% | 3.0% | 8.5% | 4.5% | 5.5% | 6.0% |
| 250 | 0.126491 | 1.0% | 1.5% | 7.0% | 13.0% | 9.0% | 9.0% | 9.5% |
| 250 | 0.5 | 0.0% | 1.0% | 4.0% | 9.5% | 5.0% | 5.0% | 5.0% |
| 500 | 0 | 7.0% | 6.5% | 3.5% | 14.5% | 7.5% | 11.5% | 13.5% |
| 500 | 0.0223607 | 2.0% | 3.0% | 4.0% | 9.5% | 7.0% | 7.5% | 7.0% |
| 500 | 0.0447214 | 3.0% | 2.5% | 5.5% | 12.0% | 8.0% | 9.5% | 9.5% |
| 500 | 0.0894427 | 0.0% | 2.0% | 5.5% | 13.0% | 7.5% | 7.0% | 7.0% |
| 500 | 0.1 | 0.0% | 2.5% | 6.0% | 13.5% | 8.5% | 8.0% | 8.0% |
| 500 | 0.5 | 0.0% | 0.0% | 5.0% | 13.0% | 6.0% | 5.0% | 5.0% |

Ranges across these 17 different laws are descriptive extrema, not pooled size:

| Method | Observed 5% null-rate range |
|---|---:|
| ADF n | 0.0%–7.0% |
| ADF c | 0.0%–7.0% |
| ADF ct | 3.0%–7.0% |
| Original | 8.5%–15.5% |
| Conservative | 3.0%–11.5% |
| AIC→ADF | 3.0%–14.5% |
| BIC→ADF | 3.0%–14.5% |

A single-cell uncertainty example, T=100 and zero drift:

| Method | Rejections | Rate | MCSE | Wilson 95% interval |
|---|---:|---:|---:|---|
| ADF n | 11/200 | 5.5% | 1.6% | [3.1%, 9.6%] |
| ADF c | 12/200 | 6.0% | 1.7% | [3.5%, 10.2%] |
| ADF ct | 12/200 | 6.0% | 1.7% | [3.5%, 10.2%] |
| Original | 29/200 | 14.5% | 2.5% | [10.3%, 20.0%] |
| Conservative | 19/200 | 9.5% | 2.1% | [6.2%, 14.4%] |
| AIC→ADF | 24/200 | 12.0% | 2.3% | [8.2%, 17.2%] |
| BIC→ADF | 25/200 | 12.5% | 2.3% | [8.6%, 17.8%] |

The repeated elevated adaptive rates warrant methodological review before main
promotion. They do not establish their limits as T grows, reveal uniform size,
or justify selecting whichever method happened to perform best in this pilot.

## 3. Raw power under alternatives

Zero-intercept alternatives (the overlap of family I and M with m=0):

| T | ρ | ADF n | ADF c | ADF ct | Original | Conservative | AIC→ADF | BIC→ADF |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 0.900 | 79.0% | 29.5% | 16.5% | 80.5% | 33.5% | 74.5% | 80.5% |
| 100 | 0.950 | 30.5% | 13.5% | 8.5% | 36.5% | 16.0% | 33.0% | 36.5% |
| 100 | 0.990 | 9.5% | 6.0% | 7.0% | 17.0% | 10.5% | 14.0% | 17.0% |
| 250 | 0.900 | 100.0% | 96.5% | 85.5% | 100.0% | 97.0% | 99.5% | 100.0% |
| 250 | 0.950 | 91.0% | 46.5% | 25.5% | 92.5% | 48.5% | 84.5% | 92.5% |
| 250 | 0.960 | 76.0% | 27.0% | 16.0% | 78.0% | 31.0% | 72.5% | 77.5% |
| 250 | 0.980 | 34.0% | 7.5% | 9.0% | 38.0% | 12.5% | 33.0% | 36.0% |
| 250 | 0.990 | 12.0% | 7.5% | 3.5% | 17.5% | 9.5% | 15.0% | 17.0% |
| 500 | 0.900 | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% | 100.0% |
| 500 | 0.950 | 100.0% | 96.0% | 80.5% | 100.0% | 96.5% | 99.5% | 100.0% |
| 500 | 0.980 | 77.5% | 35.5% | 17.0% | 78.5% | 38.5% | 70.5% | 77.5% |
| 500 | 0.990 | 31.0% | 9.0% | 9.0% | 36.0% | 13.5% | 30.0% | 34.5% |

A contrast of parameterizations along ρ=1−5/T:

| T | Family label | ρ | Intercept a | Equilibrium mean | ADF n | ADF c | ADF ct | Original | Conservative | AIC→ADF | BIC→ADF |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | M: m=5 | 0.950 | 0.25 | 5 | 4.0% | 20.5% | 9.0% | 24.5% | 22.5% | 23.0% | 24.5% |
| 100 | I: a=.1 | 0.950 | 0.1 | 2 | 22.0% | 11.0% | 10.0% | 32.5% | 16.0% | 28.0% | 32.5% |
| 250 | M: m=5 | 0.980 | 0.1 | 5 | 11.5% | 9.0% | 8.5% | 23.5% | 13.5% | 18.5% | 22.5% |
| 250 | I: a=.1 | 0.980 | 0.1 | 5 | 11.5% | 9.0% | 8.5% | 23.5% | 13.5% | 18.5% | 22.5% |
| 500 | M: m=5 | 0.990 | 0.05 | 5 | 15.5% | 13.0% | 12.0% | 27.0% | 17.0% | 26.0% | 25.0% |
| 500 | I: a=.1 | 0.990 | 0.1 | 10 | 7.0% | 11.0% | 9.0% | 18.0% | 15.0% | 15.5% | 17.0% |

At T=250 the two displayed family labels identify the same law; they are one
configuration, not independent evidence. Across T, fixed mean and fixed intercept
do not generate interchangeable alternatives. These are **raw** power figures;
higher rejection under alternatives cannot be interpreted as better size-adjusted
power when the procedures have different null rejection rates. No size adjustment
was implemented or estimated. The full 83-cell power matrix and all α levels
are exported in `summary/rejection_rates.csv.gz`, with aliases and Wilson intervals.

## 4. Retained deterministic specifications

Illustrative 5% profiles. Entries are **counts n/c/ct out of 200**, not probabilities
that a deterministic component is true. The complete 6,300-row selection table
covers every method, cell, α and specification.

| T | Null a | Original | Conservative | AIC | BIC |
|---|---:|---|---|---|---|
| 100 | 0 | 177/10/13 | 0/187/13 | 51/54/95 | 137/38/25 |
| 100 | 0.5 | 3/191/6 | 0/194/6 | 0/25/175 | 23/82/95 |
| 250 | 0 | 186/7/7 | 0/193/7 | 52/57/91 | 161/23/16 |
| 250 | 0.5 | 0/192/8 | 0/192/8 | 0/28/172 | 1/115/84 |
| 500 | 0 | 183/10/7 | 0/193/7 | 61/47/92 | 171/22/7 |
| 500 | 0.5 | 0/188/12 | 0/188/12 | 0/35/165 | 0/137/63 |

IC selection is α-independent; its downstream rejection is not. Enders selection
can depend on α because it includes stopping on unit-root rejection. Even a tau
rejection while retaining ct does not establish a nonzero deterministic coefficient.

## 5. Enders terminal branches

Per-cell frequencies, MCSE and Wilson intervals are in `enders_terminals.csv.gz`;
full executed-path counts are in `enders_paths.csv.gz`. Only observed branches/paths
have rows: an absent branch for a valid cell has count zero. The following are
**design-mix descriptive counts** across the 17 null cells at 5% (3,400 observations
per variant), not a single DGP's rejection probability or a pooled binomial interval.

| Variant | Retained spec | Terminal test | Reject | Policy stop | Count | Share of null design mix |
|---|---|---|---|---|---:|---:|
| original | c | conditional_t | False | False | 700 | 20.588% |
| original | c | conditional_t | True | False | 105 | 3.088% |
| original | c | tau2 | True | False | 78 | 2.294% |
| original | ct | conditional_t | True | False | 8 | 0.235% |
| original | ct | tau3 | True | False | 166 | 4.882% |
| original | n | tau1 | False | False | 2267 | 66.676% |
| original | n | tau1 | True | False | 76 | 2.235% |
| conservative | c | tau2 | False | True | 3148 | 92.588% |
| conservative | c | tau2 | True | False | 78 | 2.294% |
| conservative | ct | conditional_t | True | False | 8 | 0.235% |
| conservative | ct | tau3 | True | False | 166 | 4.882% |

The c conditional branch in Original is bypassed by default Conservative. The
trend-stage conditional branch is shared. Ancestor and deterministic tests condition
later decisions, so conventional limits under maintained DGP assumptions do not
calibrate the adaptive rejection event. Branch counts explain which routes occurred;
they do not alone identify a causal explanation for all observed inflation.

## 6. Paired differences and uncertainty

All methods use the same data within a replication. The all-pairs CSV contains
6,300 comparisons (100 cells ×3 levels ×21 pairs), paired MCSE, discordant counts,
interval method and all-attempted failure bounds. The six prespecified primary
5% contrasts are flagged and oriented according to the protocol.

Original minus Conservative is a nested extra-rejection event. The reverse
discordant count was zero in every one of the 60,000 matched method/α comparisons.
Its null-cell 5% differences are:

| T | a/σ | Original-only rejections | Difference | Paired MCSE | Wilson 95% interval |
|---|---:|---:|---:|---:|---|
| 100 | 0 | 10/200 | 5.0% | 1.5% | [2.7%, 9.0%] |
| 100 | 0.05 | 11/200 | 5.5% | 1.6% | [3.1%, 9.6%] |
| 100 | 0.1 | 9/200 | 4.5% | 1.5% | [2.4%, 8.3%] |
| 100 | 0.2 | 12/200 | 6.0% | 1.7% | [3.5%, 10.2%] |
| 100 | 0.5 | 21/200 | 10.5% | 2.2% | [7.0%, 15.5%] |
| 250 | 0 | 13/200 | 6.5% | 1.7% | [3.8%, 10.8%] |
| 250 | 0.0316228 | 11/200 | 5.5% | 1.6% | [3.1%, 9.6%] |
| 250 | 0.0632456 | 7/200 | 3.5% | 1.3% | [1.7%, 7.0%] |
| 250 | 0.1 | 8/200 | 4.0% | 1.4% | [2.0%, 7.7%] |
| 250 | 0.126491 | 8/200 | 4.0% | 1.4% | [2.0%, 7.7%] |
| 250 | 0.5 | 9/200 | 4.5% | 1.5% | [2.4%, 8.3%] |
| 500 | 0 | 14/200 | 7.0% | 1.8% | [4.2%, 11.4%] |
| 500 | 0.0223607 | 5/200 | 2.5% | 1.1% | [1.1%, 5.7%] |
| 500 | 0.0447214 | 8/200 | 4.0% | 1.4% | [2.0%, 7.7%] |
| 500 | 0.0894427 | 11/200 | 5.5% | 1.6% | [3.1%, 9.6%] |
| 500 | 0.1 | 10/200 | 5.0% | 1.5% | [2.7%, 9.0%] |
| 500 | 0.5 | 14/200 | 7.0% | 1.8% | [4.2%, 11.4%] |

For non-nested pairs, reported intervals use conservative simultaneous exact
bounds for the two discordant probabilities. No independent-binomial variance
formula is substituted for the paired variance. With no failures, valid-only and
attempted-denominator rates coincide. The implementation nevertheless preserves
failure masks and lower/upper sensitivity bounds rather than coding failure as
nonrejection.

## 7. IC margins

Near-tie diagnostic is prespecified gap<.1 and never changes the exact comparator.
The median range is across cell-specific medians; counts below .1 are workload-
weighted descriptive diagnostics, not a pooled inferential estimate.

| Criterion | DGP block | Valid selections | Exact ties | Gap<.1 | Range of per-cell median gaps |
|---|---|---:|---:|---:|---:|
| AIC | null | 3400 | 0 | 86 | 1.435–2.875 |
| AIC | alternative | 16600 | 0 | 575 | 1.112–1.882 |
| BIC | null | 3400 | 0 | 66 | 1.639–4.047 |
| BIC | alternative | 16600 | 0 | 209 | 1.671–5.723 |

Every stored winner, score formula and gap was checked independently during the
streaming audit. Retain margins when comparing numerical builds: very close
comparisons may change under floating-point differences without changing the
specified exact tie rule.

## 8. Promotion review

- Review the sizable adaptive null rejection rates and conventional-t branch
  interpretation before authorizing main. Do not change the grid/rules because
  a method has attractive power or unattractive size.
- R=200 is exploratory; rates near .05 have roughly ±3 percentage-point normal
  half-width. Per-cell Wilson intervals are provided, but not simultaneous claims.
- Fixed-mean versus intercept-continuation families and zero-start transients
  must remain separate interpretations. Broader claims need the deferred designs.
- No platform comparison on Linux/Colab/Kaggle was executed here; those platforms
  should replay the fixed numerical fixtures and stream hashes before main.
- Review runtime/memory in the separate estimates. Main is deliberately disabled
  in this runner; it requires separate authorization and a reviewed phase/count
  extension. Any later optimization must preserve complete public-method parity.

The observed empirical behavior is not a reason to alter the public econometric
implementation in this assignment. Software validation passed; methodological
interpretation and main-study authorization remain the user's next review step.
