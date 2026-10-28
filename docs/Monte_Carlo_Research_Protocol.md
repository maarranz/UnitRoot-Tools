# Monte Carlo research protocol: deterministic selection and ADF inference

UnitRoot-Tools · 28 October 2026 · Protocol proposal v1.0

**Planning only. No simulation runner has been implemented and no simulations have been run.** The user-specified study requirements are distinguished below from proposed methodological choices that require approval before implementation. This protocol does not authorize new package functions, dependency changes, commits or releases.

## 1. Objectives, estimands and hypotheses

Primary question: how does data-dependent deterministic-term selection change unit-root rejection probabilities, especially when drift is small relative to innovation noise?

Estimate, for each completely specified DGP, sample size, lag policy and method: its unit-root rejection probability, deterministic-specification frequencies, terminal-branch frequencies, numerical-failure rate, and paired differences from comparison methods. Separate raw power from any separately calibrated power. A selected regression or a component p-value is not itself a calibrated testing procedure.

Predeclared hypotheses to investigate, not assumed findings:

- Drift of order T⁻¹ᐟ² produces different finite-sample behavior from zero drift and fixed nonzero drift; pointwise conventional-limit arguments need not approximate the local regime well.
- Adaptive procedures may alter rejection rates through both model selection and the calibration used after selection. Extra rejection can represent size distortion rather than improved power.
- Default Conservative Enders has a rejection set contained in Original Enders for matched inputs/configuration. This exact logical property supplies a software invariant, not a claim of nominal size control or better size-adjusted power.
- IC margins and Enders branch frequencies help explain where rejection rates change. Small margins are diagnostics, not model probabilities.
- Correctly specified short-run dynamics and fixed-order approximations should not be pooled; residual serial correlation can otherwise be misattributed to deterministic selection.

The protocol follows an aims/DGPs/estimands/methods/performance-measures structure and reports Monte Carlo uncertainty explicitly, consistent with [Morris, White and Crowther (2019), Statistics in Medicine 38, 2074–2102](https://pmc.ncbi.nlm.nih.gov/articles/PMC6492164/).

### 1.1 Fixed requirements versus proposed choices

| Status | Decision |
|---|---|
| User-specified | Seven methods below; no eighth Hamilton/Enders method or permitted-Conservative method |
| User-specified | A: Gaussian AR(1), p=0; B: fixed short-run dynamics; C: automatic lags deferred |
| User-specified | Effective ADF T=100,250,500; pilot R=200; main R=2,000; primary α=.05, secondary .01/.10 |
| User-supplied candidate grid, adopted with exact deduplication | Six drift labels and five root labels in §4 |
| Proposed, approval required | Mean/intercept alternative families; zero-start primary design; 8 symmetry and 16 initial-condition cells |
| Proposed, approval required | Exact AR dynamics and invertible MA approximation grids in §5; long-run noise normalization |
| Proposed, approval required | Optional pointwise size calibration, its additional training budget and research-only threshold rule in §8 |
| Proposed, approval required | RNG scheme, pilot gates, storage/checkpoint design and staged workload |

The recommended first executable milestone is **Experiment A core: 100 distinct cells**. The complete proposed A/B package has **300 analysis cells**. The larger package is a costed proposal, not an implicit instruction to run it now.

## 2. Implementation checkpoint and exact seven-method comparison

Read-only inspection confirms branch `main`, commit `b6a84fec7ca28a1802cdd722c163f624d16dca24`, version `0.1.0rc1`. This includes the corrected normalized-rank handling in the experimental selector. The tracked working tree was clean; existing untracked reports/audits were preserved. The historical source/wheel/sdist validation recorded 1,367 passes after that correction; this protocol does not repeat those tests or turn their success into statistical size evidence.

Implementation references: [ADF adapter](../src/unitroot_tools/_adapters.py), [fixed-test result and rejection rule](../src/unitroot_tools/results.py), [Enders implementation](../src/unitroot_tools/_enders.py), [Enders documentation](Enders.md), [deterministic selector](../src/unitroot_tools/_deterministic.py), [IC methodology](Deterministic_Selection.md), and [Hamilton/Enders theoretical review](Hamilton_Enders_Theoretical_Review.md). The theoretical review supplies no independently verified new Hamilton selection tree.

| Method ID | Operation at fixed common p | Raw rejection rule at α | Output retained |
|---|---|---|---|
| `adf_n` | `adf(y, trend="n", lags=p)` | τn < ARCH critical value for n | statistic, critical values, p, sample, versions |
| `adf_c` | Same with c | τc < ARCH critical value for c | Same |
| `adf_ct` | Same with ct | τct < ARCH critical value for ct | Same |
| `enders_original` | `enders(y, variant="original", lags=p, alpha=α)` | Existing complete adaptive tree | unit-root decision, retained specification, full executed history and terminal test |
| `enders_conservative` | Same with conservative and default no-constant permission false | Same trend subtree; at c tau nonrejection stop with nonrejection | Same plus policy-stop indicator |
| `det_aic` | `select_deterministic(y, lags=p, criterion="aic")` | **Research pipeline:** selected fixed-specification ADF statistic < its ARCH critical value | all candidate scores, winner, gap, selected fixed ADF result |
| `det_bic` | Same with bic | Same pipeline with BIC winner | Same |

All comparisons are strict: equality does not reject. Use critical values, not approximate p-values, to define decisions. The IC selectors themselves expose no p-value/reject method. The study explicitly defines a downstream testing pipeline using `selected_fixed_spec_adf.reject(α)`; it does not add a selector-level inference API.

**Fixed-specification ADF p-value reported after selection.** Such a p-value ignores deterministic selection and any upstream lag selection. It is not a selective p-value conditional on winning and establishes no overall size control. No automatic lag selection occurs in A/B, but this qualification remains essential for later C.

### 2.1 Enders routing and calibration, frozen for raw comparisons

Write the regression as Δy=a+γt+b y[-1]+ΣβjΔy[-j]+e. At ct, test b=0 with nonstandard τ3. If nonrejected, test (b,γ)=(0,0) with φ3. A joint nonrejection moves to c; rejection leads to the b=0 auxiliary regression and a two-sided trend t test. An auxiliary nonrejection also moves to c; rejection leads to a terminal conventional left-tail test of b in the original ct regression.

At c, τ2 rejection terminates. Original otherwise tests (b,a)=(0,0) with φ1, then, if needed, an auxiliary intercept test imposing b=0. Nonrejection at either selection gate moves to n and terminal τ1; auxiliary rejection leads to terminal conventional b testing in c. Default Conservative stops immediately at c τ2 nonrejection, bypassing all these later branches. Permitted Conservative equals Original and is only an equivalence diagnostic.

- Tau critical values: unchanged ARCH 8.0.0 ADF calibration; MacKinnon approximate p-values are descriptive.
- Phi: `(RSS_R−RSS_U)/2 / (RSS_U/df_U)`; nonstandard Enders Table B, never conventional F. Grid row is largest Tgrid≤effective T, capped at 500. Here T=100,250,500 exactly matches rows. At 5%, φ1=(4.71,4.63,4.61), φ3=(6.49,6.34,6.30).
- Auxiliary t: two-sided finite-Student thresholds using auxiliary residual df. Conditional b test: left-tail finite-Student threshold using unrestricted residual df. Preserve these conventions even though their theoretical motivation is conditional/asymptotic, not exact finite-sample Student inference.
- Joint or deterministic-coefficient rejection alone never means unit-root rejection. No φ2 occurs.

The manual's deterministic-component argument does not establish the assumptions merely because a pretest rejects. Fixed nonzero drift without a separately fitted trend and a drifted random walk fitted with ct have different asymptotic experiments. The resulting vulnerability is part of what this study measures, not a rationale to modify Enders during the experiment.

### 2.2 IC conventions, frozen

Always compare unrestricted n/c/ct at one supplied p and identical response rows. For deterministic dimension d=0,1,2, q=1+p+d coefficients, k=q+1 including variance, and n=T:

AIC = T log(RSS/T)+2k; BIC = T log(RSS/T)+log(T)k.

Scores omit only the common Gaussian constant. Both are stored. Use the global unrounded minimum; exact ties prefer n, then c, then ct. No near-tie tolerance or winner refit. One failed candidate invalidates that selector comparison; never choose among survivors.

## 3. Time, initialization, noise and sample conventions

**T means effective final ADF regression observations**, not supplied levels, total generated history, or a selected-lag comparison sample. With fixed p, supply N=T+p+1 levels y0,…,y(T+p). Each method uses response-level indices p+1,…,T+p, i.e. zero-based slice `[p+1:N]`. No values are dropped and no frequency is inferred. The ct regressor is the existing backend's estimation-sample sequence 1,…,T.

A primary cell generates all required initial/extra observations under its stated DGP and then conditions on the first p+1 supplied levels for estimation. Report both N and T. Every method in a cell uses identical y and p; no per-specification lag choices are allowed.

Experiment A primary initialization: y0=0 and **no burn-in**, under null and alternatives. Thus stable-root alternatives are zero-start conditional AR processes, not processes claimed to be strictly stationary at the first observation. Exact stationary initialization is an explicit sensitivity block. Burn-in of a near-unit-root process would need sample-size-dependent justification and would change the experiment; it is not silently added.

Set σ=1 in A. Multiplying every coefficient, state and innovation by the same positive unit scale gives no additional substantive DGP in exact arithmetic. Use a small fixed-fixture software scale check in preflight, not another main simulation factor. Initial level shifts are different: including n destroys general translation invariance.

For B normalize the **long-run standard deviation of the difference innovation under the null** to ω=1, so drift/noise contrasts remain comparable across dynamics. Also report one-step innovation SD and short-run coefficients. Distinguish this normalization from A's σ, where σ=ω. The formulas in §5 define it exactly.

## 4. Experiment A: AR(1), fixed p=0

### 4.1 Unit-root null

For t=1,…,T:

\[
y_t=y_{t-1}+\mu_T+\sigma z_t,\qquad z_t\overset{iid}{\sim}N(0,1),\quad y_0=0.
\]

Use drift labels δT=μT/σ in

\[
\mathcal D_T=\{0,\,1/(2\sqrt T),\,1/\sqrt T,\,2/\sqrt T,\,1/10,\,1/2\}.
\]

This parameterization keeps the local and fixed regimes distinct: under local δ=c/√T the deterministic level displacement and cumulative innovation scale are both order √T. Under fixed nonzero δ, deterministic displacement is order T. A nonzero drift in levels does **not** imply a nonzero time coefficient in the difference equation.

At T=100, 1/√T=0.1: simulate that law once, retaining both `local_c1` and `fixed_0.1` labels. There are 5,6,6 distinct null laws at T=100,250,500, respectively: **17 core null cells**.

### 4.2 Alternative family I: retain the regression intercept

\[
y_t=\rho y_{t-1}+a_T+\sigma z_t,\quad a_T/\sigma\in\mathcal D_T,\quad y_0=0,
\]

with candidate roots

\[
\mathcal R_T=\{9/10,19/20,99/100,1-5/T,1-10/T\}.
\]

Here aT is the same intercept index as μT under the matching unit-root null. For ρ<1 its equilibrium mean is aT/(1−ρ). It is **not** a mean-preserving alternative. For ρ=1−c/T, fixed a yields an equilibrium mean of order T; local a of order T⁻¹ᐟ² yields one of order √T. These are meaningful conditional-start alternatives but can look strongly trending over the observed sample. Do not call their power differences evidence about fixed-mean stationary alternatives.

The deduplicated root sets are:

| T | Distinct roots | Root aliases retained |
|---|---|---|
| 100 | .90, .95, .99 | 1−10/T=.90; 1−5/T=.95 |
| 250 | .90, .95, .96, .98, .99 | No coincidence |
| 500 | .90, .95, .98, .99 | 1−5/T=.99 |

A fixed .99 and a local-root expression can coincide at one T without defining the same sequence as T varies. Preserve both labels while simulating the common finite-sample law once.

### 4.3 Alternative family M: preserve a fixed mean

Proposed additional family:

\[
y_t=m+\rho(y_{t-1}-m)+\sigma z_t,\quad m/\sigma\in\{0,1,5\},\quad\rho\in\mathcal R_T,\quad y_0=0.
\]

The intercept is (1−ρ)m. Under local-to-unity roots and fixed m this is order T⁻¹, a different experiment from local drift T⁻¹ᐟ². At the unit-root boundary the intercept vanishes and m is unidentified; with y0 fixed at zero the matching null is the zero-drift null. This family supplies the coherent fixed-mean comparison missing from an intercept-only grid. Its zero initialization intentionally includes a transient toward m.

Deduplicate identical `(T,ρ,a/σ,initial-law)` configurations across I and M while retaining all family aliases. For example, at T=100 and ρ=.90, m=1 gives a=.1 and m=5 gives a=.5, already present in I. Do not give those laws double weight in pooled summaries.

### 4.4 Complete core matrix and counts

| T | Nulls | I alternatives | M alternatives | I∩M alternatives | Distinct alternatives | Total A core |
|---|---:|---:|---:|---:|---:|---:|
| 100 | 5 | 15 | 9 | 7 | 17 | 22 |
| 250 | 6 | 30 | 15 | 8 | 37 | 43 |
| 500 | 6 | 24 | 12 | 7 | 29 | 35 |
| Total | **17** | 69 | 36 | 22 | **83** | **100** |

The matrix is the Cartesian products defined above, followed by exact law deduplication; it is not an unspecified sampling of interesting cells. Use rational roots and reduced algebraic drift descriptors in the configuration catalogue, not broad floating-point “near-duplicate” tolerances. Expose excluded aliases and their canonical identifiers.

### 4.5 Limited negative-drift symmetry block: 8 cells

For each T, add null δ∈{−1/√T,−.1}, deduplicating their coincidence at T=100: 5 cells. Add I alternatives with ρ=.95 and a/σ=−1/√T at each T: 3 cells.

Generate these as exact sign reflections of their positive counterparts, including innovations and initial state. This is a paired symmetry diagnostic, not independent evidence from extra random samples. All seven rules should be invariant to sign reversal in exact arithmetic: level statistics and phi values are invariant and deterministic t tests are two-sided. Report any discrepancy with input hashes and margins; do not silently average it away. Each reflected cell still has R independent replication indices internally.

### 4.6 Initial-condition block: 16 cells

Two proposed sensitivities, separate from core estimates:

1. **Large deterministic initial level:** for each T and null δ∈{0,.1}, replace y0=0 with y0=5σ√T. This adds 6 null cells and investigates an initial state on the unit-root stochastic scale.
2. **Exactly stationary initialization and matched null:** for m=σ and each distinct ρ∈{.95,1−5/T}, draw y0=m+σ Zinit/√(1−ρ²), independently of future innovations. Simulate the stable AR(1) alternative with this law; simulate its zero-drift unit-root comparator with the **same initial distribution descriptor**. This adds 5 alternative and 5 null cells after the T=100 root coincidence is removed. The ρ label in the matched null identifies the initial variance, not its unit autoregressive coefficient.

No burn-in is needed for the exactly stationary alternative. Comparing it with a zero-initialized null would confound initial conditions with calibration; the matched null avoids that mistake. These cells are sensitivity analyses, not pooled replicates of the core. A now totals **124 cells: 33 null, 91 alternative**.

## 5. Experiment B: short-run dynamics, fixed lags

Define \(\mathcal H_T=\{0,1/\sqrt T,.1\}\) and \(\mathcal Q_T=\{1,.95,1-5/T\}\), deduplicating both at T=100 as applicable. There are 4,9,9 distinct `(drift-index,root)` pairs at T=100,250,500. For alternatives, the drift index is a continuation parameter; only under ρ=1 is it the mean increment.

### 5.1 B-AR: exact finite ADF representation, 44 cells

For ψ∈{−.5,.5}, let

\[
(1-\rho L)(1-\psi L)y_t=(1-\psi)\mu_T+(1-\psi)\omega z_t,
\quad \mu_T/\omega\in\mathcal H_T,\quad\rho\in\mathcal Q_T,\quad\omega=1.
\]

Initialize y(−1)=y0=0; no burn-in. Generate y1,…,y(T+1); supply y0,…,y(T+1). The exact ADF form is

\[
\Delta y_t=(1-\psi)\mu_T+(\rho-1)(1-\psi)y_{t-1}
 +\rho\psi\Delta y_{t-1}+(1-\psi)\omega z_t.
\]

Use **p=1** for every method. The level coefficient is (ρ−1)(1−ψ), not ρ−1, and the lagged-difference coefficient is ρψ. Under the null, increments follow a stable AR(1) with mean μT and long-run innovation SD ω. Their initial increment is zero rather than drawn from stationarity; the stable transient is explicit. Under alternatives the equilibrium mean remains μT/(1−ρ), matching family I's continuation interpretation.

Other roots are stable because |ψ|<1; when ρ<1 both factors are stable. Counts: **8,18,18**, total **44**, of which 16 are null and 28 alternative. This is a correctly specified finite-order ADF experiment; fixed p does not remove the separate deterministic-model misspecification of individual methods.

### 5.2 B-MA: invertible errors, finite-order approximation, 132 cells

For θ∈{−.5,.5}, define

\[
y_t=\rho y_{t-1}+\mu_T+v_t,\qquad
v_t=\frac{\omega}{1+\theta}(z_t+\theta z_{t-1}),
\quad\mu_T/\omega\in\mathcal H_T,\quad\rho\in\mathcal Q_T,\quad\omega=1.
\]

Set y0=0, draw z0 independently with the future standard normals, and use **p∈{1,4,8}** as three fixed truncation settings. No nonzero finite p exactly whitens a general MA(1); these are approximation experiments. Invertibility holds and the null long-run SD is ω. One-step innovation SD is ω/(1+θ), so negative and positive θ differ substantially in their short-run variability even with matched long-run scale.

Use one underlying level path through T+8 for all three p settings and take the appropriate prefixes of N=T+p+1 levels. Each has T final regression rows, but their response windows start/end at different times. Do not claim row-for-row matching **across p**; row matching holds across the seven methods **within a cell**. Save all three p settings as paired analysis cells sharing a generator identifier. A common-maximum-lag response window would be another design, not a silent replacement.

Counts: **24,54,54**, total **132**, of which 48 are null and 84 alternative. There are 44 underlying `(T,θ,μ,ρ)` generator configurations, each analyzed at three p settings.

### 5.3 Diagnostics and interpretation

Record residual autocorrelations at a prespecified set of lags and a descriptive whiteness diagnostic, without using them to change p or discard outcomes. No rule assumes p=8 is sufficient for every MA configuration. Negative-MA size distortion must be reported as a dynamics/truncation interaction, not solely a failure of deterministic selection.

B adds **176 analysis cells: 64 null, 112 alternative**. ψ=0 or θ=0 would overlap the A-type innovation model, so they are not added as redundant dynamic factors; a positive-p overaugmentation baseline could be a separately approved extension.

## 6. Deferred experiments and limits of the initial grid

### 6.1 C: automatic lags, not included in A/B counts

Reserve MAIC/MBIC/AIC/BIC lag selection for a later protocol amendment. Specify K, common selection sample, detrending, once-per-series versus per-specification selection, final sample and failure rules before that experiment. Do not mix existing automatic defaults into fixed-p calls. A fair extension could select p once under ct for all seven methods, but that is a **proposed** common-lag experiment, not the same as each test choosing p independently.

### 6.2 Linear level trends and genuine difference trends are separate extensions

Core A/B have no genuine nonzero γt term in their difference-equation DGP. They answer the drift-near-zero question and assess unnecessary ct inclusion; they cannot establish performance when ct is necessary for a stationary linear-trend alternative.

A coherent **linear-level-trend** extension is y_t=m+κT t+x_t with x_t=ρx_(t−1)+σzt. Its ADF coefficients at p=0 are a=(1−ρ)m+ρκT, γ=(1−ρ)κT and b=ρ−1. Under the unit-root null it has constant drift κT and γ=0; under stable alternatives γ may be nonzero. State initialization of x, not just y, must be fixed.

A **genuine difference-trend** extension is Δy_t=aT+γT t+b y_(t−1)+…+e_t with b=0 and γT≠0. It permits quadratic level mean growth. To put its accumulated trend on the stochastic √T scale at p=0 requires γT/σ of order T⁻³ᐟ², not T⁻¹ᐟ². Fixed difference trends are a different asymptotic regime again.

Recommend defining separate grids and budgets for these extensions before making broad claims about all deterministic-trend choices. They are **not** silently included in the 300 cells, and an intercept-significance result is not evidence of a difference-equation trend. Breaks, heteroskedasticity and non-Gaussian errors are also outside the initial Gaussian study.

## 7. Replications, outputs and Monte Carlo uncertainty

Pilot R=200 and main R=2,000 are separate stages with independent random streams. Pilot findings guide debugging and a documented protocol amendment, not discretionary deletion of unfavorable cells. Freeze the final matrix and decision rules before main execution.

### 7.1 Rejection rates and precision

For R independent replication indices within a cell, let I(r,m,α) be a method's rejection indicator. With no failures:

\[
\widehat\pi=R^{-1}\sum I_r,\qquad
MCSE(\widehat\pi)=\sqrt{\widehat\pi(1-\widehat\pi)/R}.
\]

Use Wilson 95% intervals for individual proportions; normal half-widths below describe anticipated precision, not the interval rule at boundaries.

| True probability | Pilot MCSE, R=200 | Approx. 95% half-width | Main MCSE, R=2,000 | Approx. 95% half-width |
|---|---:|---:|---:|---:|
| .01 | .00704 | .01379 | .00222 | .00436 |
| .05 | .01541 | .03021 | .00487 | .00955 |
| .10 | .02121 | .04158 | .00671 | .01315 |
| .50 | .03536 | .06930 | .01118 | .02191 |

At 5%, main uncertainty is about ±0.96 percentage points. The pilot is for implementation and large-effect diagnosis, not fine ranking. The 1% main level expects only about 20 rejections under nominal behavior and is exploratory unless more precision is approved.

Additional replications require a decision-relevant target. For a prespecified MCSE e, use R≥π(1−π)/e²: at π=.05, e=.0025 needs about 7,600 replications; a 95% half-width of .005 needs about 7,300. Recommend an independent fixed-size confirmatory tranche only when a scientifically relevant contrast is unresolved at R=2,000. Report that confirmation tranche independently; do not naively pool an outcome-selected first stage with it under fixed-R intervals. Do not stop adaptively when a desired significance result appears, or treat repeated interim Wilson intervals as a confidence sequence.

### 7.2 Size versus null rejection under misspecification

Report every method's null rejection rate against the same DGP, but annotate whether its maintained deterministic/dynamic assumptions are satisfied. A no-constant test applied to a drifted DGP is an omitted-deterministic-term comparison, not a clean nominal-size verification of its maintained calibration. Similarly, finite MA truncation may leave misspecified errors. “Empirical size” in a global comparison means rejection under the declared stochastic-unit-root DGP family; it is not automatic validation of each component's narrower null assumptions.

Under ρ<1 report raw power conditional on the stated family and initialization. Fixed ρ=.90/.95/.99 and ρT=1−c/T answer different questions. Show T-specific values and their aliases rather than a misleading pooled root label.

### 7.3 Required diagnostics

- Selection frequencies over n/c/ct for IC and Enders, separately by α where routing depends on α. Fixed ADFs have their prespecified model frequency one.
- Enders terminal frequencies: τ3, conditional-ct, τ2 rejection, conservative c-policy nonrejection, conditional-c, τ1 rejection/nonrejection; also full path identifiers. Path labels include the decision, not just a reused `conditional_t` name.
- Conditional rejection frequencies by selected model or terminal branch, with branch counts. These are descriptive selection diagnostics, **not** selective p-value calibration.
- IC candidate scores and differences; winning gap, relative and absolute gap summaries, exact ties and a declared descriptive near-tie bin (e.g. gap<.1). Binning never changes selection. Store both AIC/BIC margins and selected specifications once per series.
- Numerical warning/failure counts, exception categories, candidate/stage, df, condition diagnostics where available, and cross-platform replay discrepancies.
- For core estimating equations, γ=0; classifying a selected ct as proven “trend present” is wrong. Under a unit root, fixed mean m is unidentified. Use “retained specification” and “unnecessary included term relative to the stated equation,” with those qualifications.

### 7.4 Paired contrasts and multiplicity

All seven methods use common data within each replication. For methods A/B, d_r=I(r,A)−I(r,B), estimate Δ=mean(d) and MCSE(Δ)=sd(d)/√R. Retain the discordant counts (A=1,B=0) and (A=0,B=1); independent-binomial SEs discard useful pairing.

Use a paired normal interval when discordant counts support it. At rare/zero discordance, give a conservative interval from simultaneous 97.5% exact-binomial intervals for the two discordant probabilities and subtract their bounds; the union bound does not require them to be independent. For Original−Conservative, d∈{0,1} on valid matched outputs, so an interval for the extra-rejection proportion is directly appropriate.

Primary 5% contrasts: Original−Conservative; AIC−BIC; Original−ADFct; Conservative−ADFct; AIC−ADFct; BIC−ADFct. Other comparisons and 1%/10% results are secondary. Default intervals are **pointwise Monte Carlo intervals**; hundreds of cells do not support simultaneous claims from unadjusted intervals. Prespecify any confirmatory contrast family and adjustment before using those results to make such claims. Do not treat methods, α levels, reflected signs or shared MA-p analyses as independent replications.

### 7.5 Numerical failures

Never redraw a failed replication or silently change p, drop a candidate, rescale the data, or repair an estimator mid-study. Keep attempted R fixed, record successes S and failures F, and report the valid-only rejection fraction K/S clearly as conditional on successful computation. The all-attempted rejection fraction is unresolved within [K/R,(K+F)/R] absent a failure-decision policy; report these bounds rather than counting failure as nonrejection.

For paired contrasts, use common-valid replications and disclose their count and selection. Also average per-replication lower/upper possible differences over all attempts to obtain failure-sensitivity bounds. Preserve successful methods' results when another method fails. Invariant failures are a reason to pause/diagnose the run, not to edit the historical outcomes. Conditioning warnings are logged distinctly from raised estimation failures.

## 8. Optional size-adjusted power: calibrate the whole rule

**Raw results are primary. Size-adjusted power is optional and requires separate approval.** A final selected-model ADF quantile alone cannot calibrate Original Enders, because earlier rejections and conditional conventional branches are also part of the rejection event.

### 8.1 Proposed research-only rejection-evidence statistic

At a fixed nominal α, leave all phi and auxiliary selection gates at their existing α thresholds. Define each potential unit-root node's margin as e_j=c_j(α)−t_j, positive precisely on strict rejection. This includes **all** tau and conditional-t nodes that can occur along the deterministic-gate path, not only the node finally reached in the raw execution.

For Enders, construct that path counterfactually assuming each preceding tau node did not reject: ct tau is included; if the trend gates retain ct, include conditional ct and end; otherwise include c tau. Conservative then ends. Original continues according to the constant gates, ending at conditional c or n tau. Gate results depend on data and α but not on a common shift of unit-root critical values. Define Eα=max_j e_j over this path. On valid potential fits, the existing raw rejection decision is exactly `Eα>0`: an earlier tau rejection already makes the union true.

For fixed ADF, Eα is its one margin. For IC, Eα uses the selected fixed-model margin after executing global selection. Calibration samples must repeat that selection. No probability is manufactured from these margins.

The calibrated research variant rejects `Eα>h`, equivalently shifting **every potential unit-root critical value** to c_j−h while preserving deterministic-gate thresholds. Replay the resulting tree if a calibrated terminal path is reported. Calibrating only its observed final t ratio, or changing only its final ADF critical value, is not equivalent. If a counterfactual required fit fails even when raw execution stopped earlier, record a calibration-specific failure; do not assert equivalence on an undefined statistic. Verify Eα>0 against full public-API histories before using this construction.

These are **new supplementary research calibrations**, not altered raw versions of the seven methods and not a package API change. Separately calibrating Original/Conservative can destroy their raw rejection-set nesting; that is not a contradiction.

### 8.2 Independent calibration and evaluation

For each chosen null cell, T, p, α and method, generate a separate calibration sample of M datasets using a disjoint phase stream. Sort its E values and use threshold h=E_(k), k=ceil((M+1)(1−α)), with strict exceedance (and +∞ if k>M). Under exchangeability and complete valid observations this gives an unconditional calibration/evaluation bound (M+1−k)/(M+1)≤α; ties make it conservative. It is not a conditional guarantee for a particular realized training set or a proof over nuisance parameters.

Freeze h before evaluating independent main null and alternative samples. The existing R=2,000 main samples may serve as evaluation because their streams are disjoint from calibration, provided the calibration specification was fixed without optimizing it on those outcomes. Report independent null holdout performance alongside power and distinguish evaluation MCSE conditional on h from additional calibration uncertainty. Order-statistic/binomial uncertainty or independent calibration batches can assess the latter; none is implemented here.

Proposed limited calibration subset: A's zero-start nulls with δ∈{0,1/√T,.1}, deduplicated to **8 cells**. If calibration uncertainty materially affects power comparisons, propose M=10,000 per cell: 80,000 additional training datasets. At 5% its rough probability-scale training uncertainty is .00218 versus main evaluation MCSE .00487. This extra budget is not part of the approved pilot/main R and must be authorized explicitly. Training can serve all seven methods and α levels through paired data, but Eα and its cutoff are method/α-specific.

Reference-null mapping must be declared:

- Intercept-continuation family I: match its aT to the same null μT, T and initial law.
- Fixed-mean family M with zero start: use the zero-drift boundary null, where (1−ρ)m vanishes.
- If one alternative has both family aliases, it can have two scientifically different reference-null calibrations. Label `(alternative ID, reference-null ID, rule)` explicitly; do not present one as universal “size-adjusted power.”
- Initial-condition, AR and MA cells require their own matching calibration laws; do not borrow the 8 A cutoffs. Report no adjusted power where calibration is unavailable.

### 8.3 Pointwise versus uniform claims

Pointwise calibration holds, at best, at a specified null parameter and initialization under its maintained numerical conditions. A finite-grid maximum rejection rate is a useful summary, not uniform size control over continuous drift, local sequences, dynamics and initial conditions. Mixture calibration controls a specified mixture, not each component. A common worst-grid threshold would need its own finite-grid uncertainty assessment and still would not prove continuum control. No such claim is part of this protocol.

## 9. Platform-independent reproducibility architecture

Proposed layers, all outside the econometric package: immutable configuration catalogue; DGP generator; seven-method adapter; replicate evaluator; atomic checkpoint writer; deterministic aggregator; reports. A Linux command-line entry point and thin Colab/Kaggle notebooks should invoke the same configuration and runner. Notebook cells must not contain different statistical rules. No runner is written in this assignment.

### 9.1 Provenance and configuration identity

Pin the package to commit b6a84fe or a subsequently reviewed immutable checkpoint, plus wheel SHA-256 and public import path. Version 0.1.0rc1 alone is insufficient: historical artifacts under that version predate features. Freeze a compatible Linux environment lock and record Python, ARCH 8.0.0, NumPy, SciPy, Statsmodels, pandas, OS, architecture and BLAS/LAPACK. The earlier Mac stack is reference evidence, not a declaration that Linux/Colab/Kaggle has been tested.

Canonical sorted JSON identifies the full law, T, p, initialization, scale, stage and method conventions. Store both exact parameter descriptors and resolved float64 values. Use a SHA-256 configuration ID and retain aliases. Do not use Python's session-dependent `hash()`. Output paths and worker counts are operational settings, not parameters that change the law or seed.

### 9.2 RNG streams and parallel scheduling

Propose an explicit NumPy `Generator` with PCG64DXSM and `SeedSequence`, not an implicit default RNG or worker-ID seed. Freeze a published 128-bit master seed before execution. Form fixed-length unsigned 32-bit words from that seed plus a SHA-256 digest of `(study RNG version, generator ID, phase, replicate ID, component label)` with declared byte order. Use separate component labels for initial states and innovations, so adding an initial-state diagnostic cannot shift the innovation draws.

Each replicate is directly addressable; worker count/order, retries and resume position cannot change its input. Pilot, main, calibration, benchmark and confirmation have distinct phase identifiers. MA generator IDs exclude analysis p and generate the declared maximum-length stream once; their p-specific analysis IDs retain p. Reflected sign cases reuse the designated positive generator and transform the entire state/path.

Draw float64 standard normals in one fixed-size call per component, with fixed argument sizes for each generator ID. Do not assume changing chunk sizes of RNG calls preserves distribution-method streams. Save selected input/innovation hashes and fixture arrays to verify cross-platform behavior. NumPy's documented stream-compatibility conditions are stricter than merely sharing a seed; different builds or machines can differ. Therefore require identical methods/configurations and replay checks, not an unsupported promise of bitwise identical fitted statistics across platforms. [NumPy parallel RNG](https://numpy.org/doc/stable/reference/random/parallel.html), [PCG64DXSM](https://numpy.org/doc/stable/reference/random/bit_generators/pcg64dxsm.html), [compatibility policy](https://numpy.org/doc/stable/reference/random/compatibility.html).

Use process-level parallelism with a fixed small worker pool; limit BLAS/OpenMP numerical threads per worker to one to avoid oversubscription. Avoid forking a shared initialized RNG. CPU OLS/Python overhead is the likely bottleneck; GPU availability is not a reason to change the design or expect acceleration.

### 9.3 Checkpointing and outputs

A work unit is a deterministic `(analysis cell, phase, replication range)`; begin with chunks of 50, then tune only the operational chunk size to roughly a minute of measured work. Store completed shards atomically, with row ranges, config/code/environment hashes and content checksums. A manifest lists pending/completed tasks; resumption skips only verified complete shards. A repeated task must reproduce the same keys and never append duplicate rows. Aggregate in canonical key order and use integer counts before floating summaries.

Use portable JSON manifests, numeric NPZ shards with `allow_pickle=False`, and compressed CSV exports/summary tables. No Parquet dependency is required. Keep a compact wide decision array indexed by replication, α and seven methods; diagnostics include three tau statistics/critical values, IC scores/winners/gaps, Enders executed path/terminal metadata, failure masks and warnings. Store full histories/design arrays only for designated pilot/replay cases and failures, with a regenerable seed key for every observation. Never serialize mutable fitted objects as essential study data.

The 300-cell main design contains 600,000 cell-replicates and 12.6 million method/α decisions. Wide binary/CSV storage avoids a separate verbose row for each decision. Measure actual shard sizes in the benchmark before committing to a storage budget; keep remote storage transfers outside timed estimation.

For Colab/Kaggle, install/attach the same reviewed artifacts and lock, verify import provenance, and use notebook-specific writable/output mount settings only. Internet availability, session termination and resource quotas affect scheduling/checkpoint frequency, not T, R, lags or decision rules. Test serial, parallel and interrupted/resumed execution on the same pilot keys. Carry final manifests/artifacts to durable storage before a session ends.

## 10. Computational workload and runtime benchmark

### 10.1 Complete proposed A/B matrix

| Block | T=100 | T=250 | T=500 | Total cells | Null | Alternative |
|---|---:|---:|---:|---:|---:|---:|
| A core | 22 | 43 | 35 | 100 | 17 | 83 |
| A sign checks | 2 | 3 | 3 | 8 | 5 | 3 |
| A initial conditions | 4 | 6 | 6 | 16 | 11 | 5 |
| B exact AR | 8 | 18 | 18 | 44 | 16 | 28 |
| B MA × fixed p | 24 | 54 | 54 | 132 | 48 | 84 |
| Total | **60** | **124** | **116** | **300** | **97** | **203** |

There are 212 structural law/initialization configurations before counting the three MA p analyses separately; the eight sign checks are transformations rather than fresh independent generators. With the recommended pairing, at most 204 independent generator keys per replication index are needed. These distinctions reduce generation/storage duplication, not the 300 analysis cells or their per-cell R.

| Scope | Pilot cell-replicates | Main cell-replicates | Main method/α decisions |
|---|---:|---:|---:|
| A core only | 20,000 | 200,000 | 4.2 million |
| A plus checks | 24,800 | 248,000 | 5.208 million |
| B | 35,200 | 352,000 | 7.392 million |
| All proposed A/B | **60,000** | **600,000** | **12.6 million** |

Pilot plus main are 660,000 analysis cell-replicates, not 660,000 independent observations for a pooled statistical estimator. Optional trends, calibration, confirmation and later C are excluded.

### 10.2 Model-evaluation accounting

Reference execution per cell-replicate fits the three fixed ADFs once for all α (3 fits), calls each IC selector once (6 more fits), and calls two Enders variants at three α (6 Enders executions). Thus there are **11 public API calls**, implementing seven scientific methods.

Worst-case regression count for each α: Original ≤3 ADF +2 auxiliary fits; default Conservative ≤2 ADF +1 auxiliary fit. Across all α plus the 9 fixed/IC fits: **≤24 ADF fits +9 auxiliary OLS fits =33 full regression fits**, plus ≤9 restricted RSS least-squares solves for phi. That is ≤42 least-squares evaluations per cell-replicate, excluding backend rank/deterministic checks and other overhead. Early termination often reduces the total.

The all-block main budget is ≤19.8 million full regression fits plus ≤5.4 million restricted solves; pilot adds ≤1.98 million fits and .54 million restricted solves. Count actual calls/fits in the pilot. The first-20-replicates-per-cell permitted-Conservative equivalence check adds at most 126,000 least-squares evaluations in the full pilot; it is not an eighth performance method and need not run in main.

A later research-only fit cache could compute three ADF fits, two auxiliaries and two restrictions once and replay all α/trees/IC rankings, at ≤7 distinct regression/solve operations per cell-replicate (≤4.2 million in main). It is permitted as a **future optimization proposal only after complete statistic/history/decision parity with public calls**, including failures, warnings and counterfactual auxiliary handling. Do not refactor package code or change numerical engines for speed.

### 10.3 Small runtime benchmark, proposed but not run

Use 20 replicates in each of six representative cells: A null T100/δ0; A null T500/local c1; A fixed-mean alternative T250/ρ.95/m5; B-AR T500/ψ−.5/null/local c1; B-MA T100/θ−.5/null/p8; B-MA T500/θ+.5/ρ.95/p8. Use benchmark-only streams, two warm-up evaluations per cell, and record generation, regression, summary and I/O times, peak memory, branch counts and shard sizes.

Compare 1,2 and at most 4 CPU workers where available; time both normal and long Enders paths. Extrapolate weighted per-cell timings rather than scaling a short-path median across all cells. If measured average time is t seconds per cell-replicate, serial main time is approximately 600,000t seconds; t=.1 would imply 16.7 hours, an **illustration, not a measurement**. Measure parallel efficiency and checkpoint/cloud costs before approving the main budget. GPU benchmarking is unnecessary initially.

## 11. Pilot validation and stage gates

All items below are future work requiring implementation approval; none was executed for this protocol.

1. **Preflight:** verify package commit/artifact/import path; existing test suite and independent Enders/IC validation; reject wrong same-version packages. Record environment lock. Inspect independent source equations rather than trusting generated labels.
2. **Generator validation:** use fixed innovation vectors and algebraic recursion identities, not Monte Carlo agreement alone. Verify AR(1) expectations/variances, B-AR factorization and ADF coefficients, MA covariance/long-run scale, sample slices, initial laws and exact aliases. Check that no automatic lag path is reachable in A/B.
3. **Decision reconstruction:** independently reconstruct OLS/QR statistics, RSS restrictions, IC normalization, critical lookups and Enders histories on designated pilot vectors at all α. Enforce per-replicate Conservative⊆Original and permitted-Conservative parity; strict equality tests use controlled statistics.
4. **Symmetry/units:** use the sign-reflection block and fixed scale fixtures. Preserve the ordinary-scale validation tolerances; extreme units can have existing backend rounding and require separately declared bounds, not a relaxed global tolerance.
5. **Execution reproducibility:** serial/parallel/resume keys yield identical input hashes on the pinned build. Replay saved vectors on Linux/Colab/Kaggle, comparing statistics within recorded numerical tolerances and recording near-boundary decisions. Do not promise exact bit identity across BLAS builds.
6. **Pilot R=200:** inspect failures, path/selection frequencies, zero/rare branches, residual dependence in MA, workload and file sizes. Missing rare branches at R200 is not proof they cannot occur; controlled branch fixtures provide routing coverage. Monte Carlo estimates are exploratory.
7. **Promotion gate:** require no unexplained identity, recursion, import, nesting or parity failures. Investigate every failure; a nonzero numerical failure rate can itself be an outcome only after its cause and reporting policy are accepted. Amend/freeze the protocol explicitly before main; do not prune cells based on apparent power.
8. **Main R=2,000:** new independent streams, locked rules, checkpoints and predetermined reports. Select any confirmation/calibration follow-up by scientific precision needs under the rules in §§7–8, not by favorable outcomes.

No pilot criterion says nominal size must be correct: discovering a reproducible size distortion is a substantive result, not automatically a software defect.

## 12. Reporting, approvals and recommended order

Primary report panels: null rejection against drift/noise regime at each T; raw power against actual ρ separately for I and M; selection and Enders terminal paths; paired 5% contrasts with MC uncertainty; IC-gap distributions; dynamics/truncation and initialization sensitivities; numerical-failure sensitivity bounds. Any size-adjusted panels name their reference-null cell, training/evaluation counts and research rule. Never average duplicated laws, initialization variants or α levels as though they were replications of one estimand.

Before implementation, resolve these proposed choices:

1. **Alternative/initialization contract:** approve both I and M families, m/σ={0,1,5}, zero-start primary processes, and the separately matched stationary-initialization block; or revise and recost explicitly.
2. **Staged scope:** authorize the 100-cell A core pilot first, then its 24 checks, then the 176 B cells with the exact factors/normalization above. The 300-cell package is recommended in stages, not one unreviewed automatic run.
3. **Size-adjustment scope:** recommend raw results first. Decide whether to approve the whole-rule evidence construction, its parity tests and the optional 8×10,000 independent calibration budget. Do not substitute an easier final-stage adjustment.
4. **Trend scope:** approve deferral of linear-level/difference-trend extensions for the initial drift study, or require a separately specified matrix before broader deterministic-trend claims.
5. **Execution/precision contract:** approve the frozen RNG/configuration scheme, primary contrast family, confirmation policy and benchmark-informed resource budget. Choose/publish the master seed and lock before data generation, without trial-and-error seed selection.

Recommended next milestone is implementation of a separately reviewed research runner and generator/decision validation for **A core**, followed by its R=200 pilot. It is not a package feature or release task. No separate Hamilton/Enders selector, Monte Carlo execution, dependency installation, public API change, bootstrap, commit or push is part of the present assignment.
