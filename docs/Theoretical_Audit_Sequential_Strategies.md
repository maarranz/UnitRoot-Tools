# Theoretical audit of sequential deterministic-specification strategies

28 October 2026 — research assessment updated following the full assignment. Original-source verification remains incomplete for Perron (1988) and the final journal version of DJS; these limits constrain implementation recommendations.

Repository inspected read-only: `/Users/maarranz/GitHub/UnitRoot-Tools`.
HEAD: `8c4f1a17cf12207dba3c051be6f1de7943df895c`.

## 1. Completed versus outstanding

At the start of this continuation, HEAD was reconfirmed as `8c4f1a17cf12207dba3c051be6f1de7943df895c`. Working-tree status was `?? docs/Theoretical_Audit_Sequential_Strategies.md`; that user-supplied untracked copy was byte-identical to the first-pass report. It was preserved. Only this chat's output report was updated. Commit `329468c` adds the Enders implementation, documentation, tests and numerical-validation script relative to `44505ec`. Its successor `8c4f1a1` changes only the CSV SHA-256 recorded in `validation/enders/source.json`.

Inspected directly: `docs/Enders.md`, `src/unitroot_tools/_enders.py`, relevant test structure, `validation/validate_enders.py`, ADF calibration labels, and the Stage 4B candidate report. No tests were executed in this audit. The earlier 1,287-pass Linux statement is recovered conversation history; it is not a newly verified M1 result. The saved 1,164-case Stage 4B matrix predates Enders.

The full theoretical-audit assignment has now been supplied. The continuation examined DJS's original Banco de España working paper, West's full original article, Perron's own later account of his 1988 sequence, Ryan's final published paper and a pinned revision of its R archive. Relevant PDF equations/flowcharts were visually checked. Sections 2–4 retain the completed first-pass analysis; sections 5 onward contain the new source assessment, operational trees, architecture assessment and implementation order.

Research conclusions: no separately established Hamilton/West sequential algorithm was verified; DJS's original working-paper sequence is less operationally specific than Enders; Ryan's published tree and code require separate qualifications. Perron's 1988 full text could not be obtained, and the user also could not locate it. Its exact tree is therefore not certified. Independent numerical validation on the current Mac checkpoint and all future implementation, simulation and release work remain separate tasks. Bootstrap is deferred. No repository files were changed and no tests, implementations, commits or simulations were run.

## 2. Model and interpretation

Use

\[
\Delta y_t=a+\gamma t+b y_{t-1}+\sum_{j=1}^{p}\beta_j\Delta y_{t-j}+\varepsilon_t.
\]

`ct` retains both deterministic terms, `c` sets gamma to zero, and `n` sets both to zero. Unit-root branches test `b=0` against `b<0`. Calling that alternative stationary additionally requires stability of the remaining dynamics; a negative coefficient alone is insufficient for arbitrary autoregressions.

Deterministic coefficients in this difference equation must not be confused with deterministic coefficients in a decomposition of levels. For p=0, under b=0, summation yields an a*t component and a gamma*t*(t+1)/2 component in levels. Thus gamma nonzero permits quadratic mean growth under the unit-root null. Under a stable alternative, the same equation ordinarily yields a linear deterministic mean. Lag augmentation changes coefficients but not this basic distinction when the difference dynamics are stable.

An auxiliary regression imposed under b=0 is useful for investigating that null model. Its significance cannot establish that the maintained null is true, establish an economically appropriate model, or restore unconditional validity after selection.

## 3. Enders: source and branch verification

Reference: Walter Enders, *Supplementary Manual to accompany Applied Econometric Time Series*, fourth edition, §4.4, PDF pages 63–66. Enders explicitly describes his tree as a modification of DJS. Accordingly, “Original Enders” identifies this package's reference tree, not verified identity with DJS. [Author's manual](https://www.time-series.net/assets/docs/SupplementaryManual_all.117125921.pdf).

The numbered tree was checked against the code:

| Stage | Null/statistic and calibration | Stop or continuation |
|---|---|---|
| ct | b=0, tau3, nonstandard DF | Reject: stop; otherwise phi3 |
| ct joint | b=gamma=0, phi3, nonstandard joint | Nonreject: c; reject: auxiliary |
| ct auxiliary, b=0 | gamma=0, conventional two-sided t | Nonreject: c; reject: conditional b test, stop |
| c | b=0, tau2, nonstandard DF | Reject: stop; otherwise phi1 |
| c joint | b=a=0, phi1, nonstandard joint | Nonreject: n; reject: auxiliary |
| c auxiliary, b=0 | a=0, conventional two-sided t | Nonreject: n; reject: conditional b test, stop |
| n | b=0, tau1, nonstandard DF | Always stop |

The published conditional tests require maintained DGP assumptions; preselection supplies no guarantee. The manual itself warns about cascading inference. No phi2 branch is required by this reference tree.

### Package-specific choices verified in code

These are implementation observations, not additional claims about the publication:

- Lag order is chosen once under ct, then held fixed. Nested regressions retain the same response sample and lagged differences. Automatic selection uses its own maximum-lag comparison sample before the final fit.
- Tau critical values come through the package ADF adapter's ARCH 8.0.0 calibration.
- Phi statistics use `(RSS_R-RSS_U)/2` divided by `RSS_U/df_U`. This is an F-shaped statistic; the code correctly applies nonstandard phi critical values.
- Phi lookup uses the largest grid size no greater than effective n, with a cap at 500; n below 25 is rejected. This is a finite-grid convention, not established exact calibration for augmented or selected-lag regressions.
- Auxiliary tests use finite residual-df Student t, two-sided. Terminal conditional tests reuse the unrestricted level-coefficient t statistic with a one-sided finite-df Student critical value. Neither exact finite-sample Student inference nor uniquely prescribed tail choice follows merely from an asymptotic conventional limit.
- Strict inequalities determine rejection. Results do not advertise a post-selection p-value.

The independent-validation script reconstructs regressions by QR and checks nested statistics through a Wald calculation. This provides numerical independence from the fitting routine. Its reference comparison and calibration assumptions still require theoretical scrutiny. Controlled branch tests verify routing, not size or power.

## 4. Conservative Enders: what is genuinely different

Direct code comparison establishes that the trend stage is identical. With the default `allow_no_constant=False`, reaching c and failing to reject tau2 stops immediately. Both phi1 and the conditional constant-only t branch are bypassed. The modification therefore does more than prohibit selection of n: it also removes an additional route to rejection while retaining c.

With `allow_no_constant=True`, routing, statistics and decisions coincide with Original Enders for the same inputs and settings. Treat this as an equivalence check, not another distinct research candidate.

There is also a useful finite-sample logical result. For matched data/configuration, default Conservative can reject only on a path where Original also rejects: their common earlier paths are identical, and Conservative stops with nonrejection exactly where Original may continue. Consequently its rejection set is a subset of Original's. Its rejection probability cannot exceed Original's under any fixed DGP, including alternatives. This establishes a mechanical tradeoff, not nominal size control, superior size-adjusted power, or general statistical superiority. Trend-stage conditional inference remains unchanged.

## 5. Source record and level of verification

| Reference | Material actually examined | What it establishes / limitation |
|---|---|---|
| Dolado, Jenkinson and Sosvilla-Rivero (DJS), 1990 | Original Banco de España Working Paper 9005, *Cointegration and Unit Roots: A Survey*, printed pp. 12–15 and notes 7–8, printed pp. 38–39 | Authors' original sequence and conditional-distribution argument. This is a working-paper version, not a verified transcription of the final journal text. |
| Perron, 1988 | Publisher metadata and abstract; author bibliography | Exact publication identified, but original full text unavailable. No exact 1988 tree is certified. |
| Perron, 1992 | *Racines unitaires en macroéconomie: le cas d'une variable*, L'Actualité économique 68(1–2), 325–356, especially pp. 335–338 | Author's own later explanation explicitly attributes a sequence to his 1988 article. It supports partial reconstruction, not every missing branch. |
| Hamilton, 1994 | *Time Series Analysis*, §17.4, Cases 3–4 and discussion, pp. 494–501, in the previously inspected text reproduction | Different maintained DGP/regression cases, not a fully specified automatic selector. A clean page-image copy was not obtained in this continuation. |
| West, 1988 | Full original Econometrica 56(6), 1397–1417; introduction, Assumptions 1–6, Theorems 2.1/3.1 and covariance discussion; relevant theorem page visually checked | Nonzero-drift asymptotic normality with specified stationary regressors, regularity and consistent covariance estimation. |
| Ryan, Haslbeck and Waldorp, 2025 | Final published article, Multivariate Behavioral Research 60(3), 556–588; §§3.3–3.4, Figure 5, Figure 6, Appendices D–E | Published sequence and stated simulation design now examined directly. Figure 5 visually checked. |
| Ryan code archive | Revision `844d4fc04e81df194dc03e90ff0f9569d449039e`; `adf_flow.R`, `analysis_adfsim.R`, `aux_functions.R`, and figure-script references to saved results | Static implementation and simulation-call inspection. No R execution or regeneration of saved simulation results. |

Primary-source links: [DJS original working paper](https://www.bde.es/f/webbde/SES/Secciones/Publicaciones/PublicacionesSeriadas/DocumentosTrabajo/90/Fich/dt_9005e.pdf); [DJS final journal bibliographic record](https://doi.org/10.1111/j.1467-6419.1990.tb00088.x); [Perron 1988](https://doi.org/10.1016/0165-1889(88)90043-7); [Perron 1992 full article](https://www.erudit.org/fr/revues/ae/1992-v68-n1-2-ae2733/602070ar.pdf); [Hamilton text reproduction](https://djvu.online/file/mozWnvm1yoLAH); [West original article copy](https://www3.nd.edu/~nmark/Climate/West_Econometrica1913104.pdf); [Ryan final published PDF, with repository cover](https://pure.uva.nl/ws/files/306761281/Non-Stationarity_in_Time-Series_Analysis.pdf); [pinned Ryan archive](https://github.com/ryanoisin/NonStationaryTS/tree/844d4fc04e81df194dc03e90ff0f9569d449039e).

The earlier inability to read Ryan's paper is resolved. The earlier Perron/DJS provenance question is only partly resolved: DJS's original working paper is available; Perron's 1988 full text and the final DJS journal text remain inaccessible through the attempted public routes. Secondary reconstructions are not substituted for originals.

## 6. DJS and Perron: reconstruct only what the sources specify

### 6.1 DJS original working-paper sequence

DJS write their general DF regression as equation (9), `Δy = β0 + β1*t + γ*y[-1] + e`; translate their γ to this report's b and β1 to this report's gamma. Printed p. 14 (PDF page 15) gives the sequence below. Printed p. 15 extends the same strategy to sufficiently augmented regressions with white-noise errors. DJS discuss the usual DF t statistics and joint statistics before the sequence, but do not explicitly insert phi gates into its step-by-step paragraph. Their notes also warn that drift may have to be large for the deterministic component to dominate in finite samples.

| Node | Regression / hypothesis | Statistic and reference | Documented action |
|---|---|---|---|
| D1 | ct; H0 b=0 versus b<0 | Level-coefficient t ratio, DF trend calibration (their tau_tau) | Reject: stop; otherwise D2 |
| D2 | Trend significance **under the unit-root null** | The paragraph does not uniquely specify restricted auxiliary versus unrestricted deterministic-coefficient statistic, its critical values, or significance level | Significant: D3; insignificant: drop trend, D4 |
| D3 | Retained ct; H0 b=0 | Standardized normal reference for the level-coefficient statistic | Retest; a terminal decision is the natural reading, but detailed tail/finite-sample/stopping conventions are not enumerated |
| D4 | c; H0 b=0 versus b<0 | DF intercept calibration (their tau_mu) | Reject: stop; otherwise test deterministic significance under the null “and so on” |
| D5 | Remaining constant/no-constant progression | Only an analogical continuation is supplied | A complete c/n subtree would be an operational completion, not a literal transcription |

**Finding:** DJS establish the historical general-to-specific rationale and conditional normal retesting. The inspected original version does not uniquely specify a new production algorithm with all hypotheses, critical values and stopping rules. In particular, we cannot infer that its paragraph mandates Enders's phi3 → restricted auxiliary → conditional test sequence. Conversely, its brevity does not prove that a particular later formalization is inconsistent with its intent.

Do not publish an independent `djs` entry point by filling these gaps silently. A no-phi, restricted-auxiliary sequence would be genuinely different from current Enders, but it would be a **DJS-inspired operational variant**, with our restrictions, tails, significance levels, lag policy and stopping rules explicitly documented. Normal critical values would be asymptotic conventions; finite-df Student critical values would be a further convention. Such a variant remains an experimental candidate rather than a newly verified established test.

### 6.2 Perron's original attribution and later primary account

The 1988 target is **Pierre Perron, “Trends and random walks in macroeconomic time series: Further evidence from a new approach,” Journal of Economic Dynamics and Control 12(2–3), 297–332**. It is not the Phillips–Perron coauthored Biometrika article and not Perron's 1989 break test. Public publisher, author, working-paper and repository routes did not yield the full 1988 text; the user also could not locate it.

There is nevertheless new primary evidence. Perron's own 1992 article, pp. 335–338, says individual deterministic-coefficient t tests can depend on nuisance parameters, describes a sequence attributed to his 1988 article, and calls for a joint test. Its documented fragment is:

| Node | Regression / hypothesis | Reference / action actually stated |
|---|---|---|
| P1 | General model with constant and trend; unit-root null | Start here; if rejected, stop |
| P2 | After P1 nonrejection, investigate redundant trend using a joint null of a unit-root process **without drift** | Joint likelihood-ratio-type test, rather than an individual deterministic t test |
| P3 | After rejection of that joint null | The 1992 paragraph says a regression without a trend can be used to test the unit root |
| Remaining branches | Exact restricted model, statistic normalization, critical values, joint nonrejection route and c/n stopping rules | Not supplied completely by this passage; not certified |

The P3 wording was checked on the page image; it is not an OCR inference. We do not reverse its rejection direction to make it resemble another textbook tree. Nor do we identify “without drift” with phi3 by name: in the report's unrestricted ct difference regression, a pure driftless random walk imposes `a=gamma=b=0`, whereas the implemented phi3 imposes only `gamma=b=0`. Initial-condition and model parameterization details also matter. The missing original text is needed to settle the actual statistic and sequence.

**Finding:** there is primary support for a Perron sequence with joint testing, and primary support for DJS's conditional-retesting proposal. There is insufficient evidence to certify their complete trees as identical, or to certify two fully distinct ready-to-code published procedures. Treating “Perron,” “DJS,” and “Enders” as three independently specified algorithms would currently overstate the evidence. The same-author later account is evidence of provenance, not a substitute for the missing 1988 operational specification.

## 7. Hamilton and West: conditional inference, not another selection tree

Hamilton's Cases 3 and 4 distinguish two experiments. In `y_t = a + rho*y[t−1] + e_t`, a fixed genuinely nonzero a under rho=1 makes the level regressor trend-dominated and permits asymptotically conventional inference under the maintained conditions. Adding a separate time trend changes the result: for a random walk with or without drift, the unit-root statistic has the nonstandard trend-case distribution. Hamilton recommends a specification plausible under both null and alternative; the inspected treatment does not specify the proposed automatic “Hamilton/West” cascade.

West provides stronger original-source verification. His p. 1397 unit-root example explicitly allows a constant and lagged differences, **excluding a time trend from the stationary-regressor block**. Assumption 1(b), p. 1400, requires nonzero `E Δy`. Theorem 2.1, p. 1402, establishes joint normality with differing convergence rates; Theorem 3.1, p. 1405, extends the analysis to additional stationary regressors/instruments. The covariance discussion on pp. 1403–1406 requires a consistent covariance estimate; ordinary computer-output standard errors are not a universal replacement for the appropriate long-run covariance.

| Conditional benchmark | Maintained conditions | Terminal inference and stopping |
|---|---|---|
| c regression with genuinely nonzero drift | Stable difference dynamics; fixed nonzero long-run mean of differences; appropriate errors/covariance and regularity | H0 b=0 versus b<0; asymptotic normal left-tail statistic; stop after one test |
| ct regression under a random walk with drift | Trend coefficient in the difference equation is zero; separate time trend retained in regression | Nonstandard ct DF inference; stop after one test |
| ct regression under genuine nonzero difference-equation trend | A different null allowing quadratic mean growth in levels | Conventional-limit motivation comes from polynomial-dominance arguments; West's stationary-regressor theorem alone does not establish this expanded case |

Our inference from these results is that automatic significance-based switching requires its own assessment. A fixed nonzero drift result is pointwise; it does not supply uniform calibration as drift approaches zero. A c-only conventional test also fails as a general substitute when trend stationarity is a plausible alternative: omitting the necessary trend can destroy power. Perron himself explains this in his 1992 Rules 3–4.

**Recommendation:** no public `hamilton_west` strategy. A prespecified conditional normal benchmark is theoretically meaningful, but must state its maintained null and alternative and must not silently relabel an ordinary ADF test. It adds no automatic deterministic-selection algorithm. Preselection is not proof of the required drift or trend.

## 8. Ryan et al.: published tree versus archive implementation

Reference: Ryan, O., Haslbeck, J. M. B., and Waldorp, L. J. (2025), “Non-Stationarity in Time-Series Analysis: Modeling Stochastic and Deterministic Trends,” **Multivariate Behavioral Research 60(3), 556–588**, DOI 10.1080/00273171.2024.2436413. The article is an accessible tutorial; §§3.3–3.4 explicitly restrict the illustrated experiments to DF/AR(1), explaining ADF as an extension. Therefore the absence of augmentation in `adf_flow.R` is a scope limit, not by itself a paper/code contradiction.

### 8.1 Published decision tree

Figure 5 on printed p. 570 labels its models in the order ct, c, n. This differs from the R object's numerical naming: `m3` is ct, `m2` is c, `m1` is n.

1. In the current model, test the level coefficient for a unit root. Section 3.3 requires the relevant nonstandard DF reference.
2. On rejection, terminate the **unit-root decision** with rejection.
3. On nonrejection, if no simpler model exists, terminate with nonrejection.
4. Otherwise test simplification: `gamma=0` in ct or `a=0` in c.
5. If simplification is not rejected, adopt the simpler model and return to step 1.
6. If simplification is rejected, the diagram calls for a **normal-reference** unit-root retest in the retained model, then terminates on either outcome.

The diagram does not specify an auxiliary regression or phi gate. It does not fully tabulate deterministic-pretest calibration or the tail convention. Section 3.4 additionally discusses further simplification after a unit-root decision for descriptive model classification. The R function's immediate returns do not implement that additional classification pass.

### 8.2 Pinned code tree and calibration

The inspected archive revision is `844d4fc04e81df194dc03e90ff0f9569d449039e`. [Pinned function with line numbers](https://github.com/ryanoisin/NonStationaryTS/blob/844d4fc04e81df194dc03e90ff0f9569d449039e/adf_flow.R).

| Branch | Code statistic / calibration | Action |
|---|---|---|
| ct initial | Level t ratio; interpolated `table3` DF p approximation | Reject: stop; otherwise deterministic pretest |
| ct simplification | Absolute unrestricted trend t ratio; separate `table_m3t` | Significant: conventional retest; otherwise c |
| ct retained | `summary(m3)` level p-value: ordinary **two-sided finite-df Student t** | Stop on either outcome |
| c initial | Level t ratio; interpolated `table2` DF p approximation | Reject: stop; otherwise deterministic pretest |
| c simplification | Absolute unrestricted intercept t ratio; `table_m2c` | Significant: conventional retest; otherwise n |
| c retained | Ordinary two-sided finite-df Student level p-value | Stop on either outcome |
| n terminal | Ordinary two-sided finite-df Student level p-value; `table1` is defined but unused here | Stop on either outcome |

For the deterministic tables, the labels 1%, 2.5%, 5%, 10% index right-tail thresholds for absolute statistics in the function. Neither a comment citing Dickey–Fuller nor a numerical table establishes uniform validity over drift/initial-condition nuisance parameters. That calibration provenance needs independent verification before any reuse.

### 8.3 Separate implementation defects, methodological concerns and unaffected claims

| Finding | Classification and consequence |
|---|---|
| Final n branch uses ordinary Student p-values, lines 223–235 | Code/calibration defect relative to the unit-root problem and §3.3's nonstandard-reference requirement. The defined n-case DF table is bypassed. |
| Figure 5 says normal retest; code uses two-sided finite-df Student p-values, lines 165–178 and 205–220 | Paper/code calibration mismatch. Changing to a one-sided test would also change the implemented rejection set; it must be an explicit alternative definition, not a silent repair. |
| Deterministic coefficient significance licenses the normal retest in Figure 5 and code | A methodological issue already present in the published tree, not solely a coding mistake. Significance does not establish the maintained DGP needed for conventional inference. |
| IC log-likelihood evaluates `xt1` rather than response `yt`, line 119 | Definite code defect in the optional IC modes. Residual-df sigma is used instead of ML `sqrt(RSS/n)`; BIC uses `log(N−1)` although N is response length. |
| `adf_translate`, lines 258 onward | Code classification defect: with a stable AR root, c gives a constant mean and ct gives a linear mean, not the printed linear/quadratic labels. It only translates the selected regression; it does not perform the post-decision simplification described in the text. |
| Paper's drift explanation on p. 570 associates growing levels with a nonzero slope in the ct DF model | Methodological explanation needs correction. A random walk with constant drift has `E Δy=a`, so the true coefficient of t in the difference equation is zero. Growth in levels does not by itself justify retaining that coefficient or switching its unit-root calibration. |
| Figure 6 and Appendix E performance summaries | Reported outcomes of a specified implementation/design, not whole-procedure size guarantees or evidence for the current Python Enders tree. |

The last explanatory criticism can be checked directly: for p=0 under the drifted-unit-root null, `y[t−1]=y0+a(t−1)+S[t−1]`. Projection on the ct regression's constant and trend removes the deterministic part, leaving the accumulated innovations. Drift therefore does not turn the ct level statistic into a conventional normal statistic. An unrestricted deterministic pretest can nevertheless react to nuisance-dependent estimates; that is a selection/calibration problem, not proof of a true difference-equation trend.

### 8.4 What the code defects do and do not establish about the published study

The pinned `analysis_adfsim.R` sets `alpha=0.01`, `modSel="pval"`, and `iter=500`. The call in `aux_functions.R`, lines 54–60, uses `aTSA::adf.test(..., nlag=1)` for the known-model comparison and calls `adf_flow` for selection. The figure script reads the saved simulation objects. Therefore the optional IC bug does **not** explain the reported p-value-based selection curves. The terminal n and retained-model calibration issues are in the called path; the translator is not called by this simulation function.

No new simulations were run, so the effect sizes of correcting these defects are not established here. The study's caution about relying mechanically on trend classification remains reasonable, but its reported error rates cannot validate a corrected algorithm or this package's different Enders implementation. Distinguish design/code correspondence from numerical reproduction: only the former was audited.

## 9. Comparative assessment and candidate definitions

| Candidate | Distinctness and theoretical status | Implementation recommendation |
|---|---|---|
| Prespecified ADF n/c/ct | Established component tests under their respective maintained conditions; no deterministic preselection | Ready using existing API; provide matched fixed-lag/sample benchmarks and economic specification guidance |
| Existing Enders manual tree | Published operational cascade with conditional-inference qualifications; package lag/sample/calibration conventions explicit | Preserve implementation and API; whole-procedure calibration remains unestablished |
| Default Conservative Enders | Our deletion of later c-stage rejection paths | Already implemented; experimental, with a guaranteed rejection-set subset but no guarantee of nominal size |
| Conservative allowing n | Same tree as Original Enders | Do not count as a distinct candidate |
| DJS-inspired restricted-auxiliary tree without phi gates | Distinct if explicitly defined; historical motivation, incomplete original operational specification | Research candidate; settle operational choices and label as our variant before coding |
| Unrestricted deterministic-pretest tree inspired by Ryan | Different from Enders and from restricted-auxiliary selection | Hold for methodological research: pretest nuisance calibration and terminal inference unresolved |
| Perron 1988 tree | Primary provenance confirmed through author's later account; exact full tree unverified | Hold original-named implementation; no placeholder alias |
| “Hamilton/West” | Conditional asymptotic theory, no separate verified selection algorithm | No named selector; optional prespecified conditional-normal research benchmark only |
| Global AIC/BIC deterministic ranking | Well-defined model-comparison exercise; not a validated selected unit-root test | Ready to specify/code as an explicitly experimental selector, with conditional fixed-model inference clearly separated |
| Calibrated union of rejections | Established alternative family requiring joint calibration | Valuable later benchmark; exact published member/calibration must be selected and independently checked first; bootstrap variants remain outside initial release |

“Ready” means the statistical operation can be specified without inventing a missing named published procedure. It does not authorize implementation or assert nominal size after selection.

### 9.1 A concrete IC candidate, without Ryan's hybrid stopping rule

For an initial, isolated candidate, prespecify p (or select it once under ct by an explicitly identified package convention), fit n/c/ct on exactly the same response rows and the same lagged differences, and calculate

\[
\widehat\sigma_d^2=RSS_d/n,\qquad
IC_d=n\log(RSS_d/n)+c_n k_d,
\]

where `c_n=2` for AIC or `log(n)` for BIC and `k_d=1+p+d+1` counts the level coefficient, p augmentation coefficients, d deterministic coefficients and the variance parameter. Common parameter counts can be omitted only consistently. Resolve exact ties by the smaller deterministic dimension, record every candidate score, and choose the global minimum once. Report selected-model ADF inference as **conditional on that specification**, together with all fixed-model results; no selected-procedure p-value or nominal-size claim.

This is not Ryan's IC mode: his tree first screens with a unit-root test and makes sequential adjacent IC comparisons, potentially stopping before the global winner. It is also not justified to compare the package's existing MAIC/MBIC scores across n/c/ct: they are lag-selection criteria on separately detrended series, with omitted constants and a purpose different from deterministic-model likelihood comparison. Joint lag/deterministic selection would be another candidate requiring a common maximum-lag sample and an explicit final-refit convention; defer it initially.

### 9.2 Reference context without adding unexamined methods to the release

DJS, Perron and Enders provide historical sequential rationale, not an exhaustive modern comparison. Ayat and Burridge (2000), Journal of Econometrics 95(1), 71–96, develop a pseudo-GLS/BIC sequence allowing quadratic trends; this audit has not verified its full tree and does not mark it implementation-ready. [Primary publisher record](https://doi.org/10.1016/S0304-4076(99)00030-5).

Harvey, Leybourne and Taylor (2009), Econometric Theory 25(3), 587–636, support a calibrated union approach under deterministic/initial-condition uncertainty. Rejecting when any unadjusted component alpha test rejects is not a justified nominal-alpha union. [Primary publication](https://www.cambridge.org/core/product/identifier/S026646660809018X/type/journal_article). This is a benchmark family to examine next, not a reason to introduce bootstrap into the initial release.

## 10. Architecture assessment: preserve the validated path

Inspected `src/unitroot_tools/_enders.py`, `_adapters.py`, `selection.py`, and `results.py`. The existing design already separates fixed-lag estimation/calibration, lag selection and the Enders decision sequence. It is adequate for the current release candidates. No broad refactor is required by this audit.

1. **Preserve Enders behavior.** Its public result types, fitting calls, phi tables, strict comparisons, lag policy and stopping decisions should remain unchanged. New algorithms must not be introduced as aliases or undocumented switches inside `enders()`.
2. **Add a separate selector only if approved.** A deterministic AIC/BIC candidate can initially call `adf(..., lags=p)` for each specification. For a common p, these fits have matching response boundaries; validate that invariant explicitly. Use the regression's raw-response RSS and actual regressor count, not `select_lags()` scores. Joint selection over different p cannot assume those final-fit samples are identical.
3. **Keep selection and inference as separate records.** Record strategy identifier/version, source, resolved configuration, candidates, selection scores, sample, lag policy and terminal decision. A retained deterministic specification is not a calibration proof. Existing Enders history already records much of the necessary provenance; avoid changing its schema merely for stylistic uniformity.
4. **Generalize only when a verified second tree needs it.** Reusable named-column restriction/auxiliary helpers may eventually be useful. Extract them only with numerical and full-history parity checks for Enders. A general decision-tree framework is unnecessary now.
5. **New calibration belongs in an explicit layer.** Any future phi2 or deterministic-pretest table requires independently checked null restrictions, nuisance conditions, source edition, sample indexing and tail. Never map every F-shaped statistic to an ordinary F distribution or every coefficient t ratio to Student t.
6. **No speculative performance work.** Optional regression storage, a different numerical core and bootstrap machinery are not prerequisites for the selected initial candidates. Consider lightweight storage later if measured simulation needs justify it; preserve current defaults and evidence.

Meaningful future validation would compare fitted quantities with an independent calculation and exercise cases where candidate trees give different paths, including retained-model nonrejections and boundary equality. Branch coverage is software evidence; it does not answer the inference question. Those checks are recommendations only and were not executed here.

## 11. Unresolved limitations and proposed implementation order

The research assessment now supports a concrete scope decision, but **original-source certification of a full Perron (1988) tree remains unfinished**. DJS is verified from its original working-paper version only. No exact published significance-level policy or complete c/n continuation was recovered from its abbreviated sequence. Ryan's deterministic-pretest tables need independent nuisance/calibration verification; its corrected rejection probabilities require later simulation. Hamilton's clean page-image verification remains a bibliographic limitation, while West's full original theorem was directly examined.

Proceed in this order, with implementation separately approved:

1. **Use and independently validate what already exists.** Fixed-specification ADF and the two distinct implemented Enders variants form the initial comparison set. Verify the present Mac checkpoint; do not repeat the old 1,287-pass claim as a new result.
2. **First optional addition: deterministic AIC/BIC ranking.** The common-p, common-sample global selector in §9.1 is operationally ready as an experimental model-selection feature. Keep its output distinct from a calibrated unit-root procedure. If the initial release should expose only established inference, defer it instead of overstating its status.
3. **Research before a second sequential implementation.** Resolve the DJS-inspired branch conventions and deterministic-pretest calibration. Define one complete candidate, with separate unit-root and selection significance levels if desired, then implement only that approved specification. Do not multiply DJS/Ryan/Hamilton names for the same tree. A prespecified conditional-normal benchmark is lower priority and should not be marketed as automatic selection.
4. **Defer the Perron-labelled strategy until original verification.** Obtain the 1988 procedure pages, identify the precise joint null/statistic and all branches, and compare the resulting rejection set with Enders/DJS before deciding whether another implementation is warranted. The current evidence does not justify a release commitment to this label.
5. **Select a calibrated non-bootstrap union benchmark for research comparison.** Verify its exact published definition and calibration before coding or using an external reference. This may improve the research design without expanding the initial public API.
6. **After strategy definitions and numerical independence checks, authorize Monte Carlo assessment.** Evaluate overall rejection probabilities under specified null families, raw and size-adjusted power, branch frequencies and sensitivity to local drift/trend, initial conditions, sample size and lag policy. Only then finalize release claims. Bootstrap remains subsequent work.

The recommendation is to release a smaller set of transparently defined procedures rather than promise three additional named strategies whose operational or inferential foundations are incomplete. This audit authorizes no code change and establishes no new size/power or platform-validation result.
