# Hamilton/Enders sequential unit-root procedure: theoretical review

28 October 2026 · UnitRoot-Tools · research only

## Finding and recommendation

**Do not implement another selector called Hamilton/Enders on the evidence currently available.** Two different things are being conflated:

1. Hamilton supplies distribution theory and specification guidance under particular maintained models. The inspected passages do not define a complete automatic deterministic-selection cascade.
2. Enders supplies an operational cascade. Its numbered instructions correspond to the package's **Original Enders**, subject to the package's explicit lag, table, tail and finite-df conventions. Giving that same tree a Hamilton/Enders name would not create a distinct method.

If “Hamilton/Enders” means a different no-phi or unrestricted deterministic-pretest algorithm, it remains **incompletely specified** until its exact source and every branch/calibration are supplied. We should not manufacture those missing steps from Hamilton's asymptotic cases.

The default Conservative Enders is genuinely different by policy, but is our modification, not a separately verified Hamilton procedure. Allowing the no-constant branches in Conservative reproduces Original Enders.

## 1. Scope and repository checkpoint

Repository inspected read-only: `/Users/maarranz/GitHub/UnitRoot-Tools`.
Branch `main`; HEAD `b6a84fec7ca28a1802cdd722c163f624d16dca24`.
Tracked working tree was clean; five existing untracked audit/validation documents were present and preserved. Reviewed `docs/Theoretical_Audit_Sequential_Strategies.md`, `docs/Enders.md`, `_enders.py`, the ADF adapter and saved validation provenance. The earlier audit describes checkpoint `8c4f1a1`; that historical checkpoint is not the current HEAD.

This assignment adds only this review in the chat's outputs folder. No package files, tests, dependencies or APIs were changed. No tests, estimation experiments, simulations, commits, pushes or releases were performed. Saved numerical evidence is identified as historical, not a new run.

## 2. Source record: what was actually verified

| Source | Edition and locator | Verification and limit |
|---|---|---|
| James D. Hamilton, *Time Series Analysis*, Princeton University Press, 1994 | §17.4, especially Cases 3–4, printed pp. 495–501; Table 17.1, pp. 502–503; §17.7 begins p. 516 | Original book text consulted through an OCR reproduction, including equations 17.4.38–17.4.55 and the p. 501 specification paragraph. Publisher/JSTOR metadata confirms the 1994 book. A clean page-image copy was **not** obtained; an attempted PDF retrieval timed out. OCR is not treated as image verification. |
| Walter Enders, *Supplementary Manual to accompany Applied Econometric Time Series*, fourth edition | §4.4, PDF physical pages 63–64; GDP discussion pp. 64–65; Figure A2.1 p. 66; Tables A/B pp. 122–123 | Downloaded from the author's site. Numbered steps, figure, equations and tables examined directly; relevant pages rendered and visually inspected. The manual is an original author source. Its title page does not supply a publication year; none is invented. |
| Enders, *Applied Econometric Time Series*, fourth edition, Wiley | Chapter 4, *Models with Trend* | Publisher/companion-site metadata verified. Wiley lists the paperback ISBN 9781118808566 as October 2014; many citations use 2015. The main book's procedure pages were **not** directly inspected, so no unverified main-book page numbers are assigned. Manual PDF page numbers must not be represented as textbook pages. |
| Sims, Stock and Watson (1990), *Inference in Linear Time Series Models with Some Unit Roots*, Econometrica 58(1), 113–144 | Theorem 1 and discussion pp. 123–124; Theorem 2, §4; §5, Examples 1–2, pp. 127–129 | Full original article obtained from Watson's university page. The nonzero-drift example and its change on adding a time regressor were inspected, including page images 128–129. These support conditional distribution results, not an automatic model-selection rule. |
| Kenneth D. West (1988), *Asymptotic Normality, When Regressors Have a Unit Root*, Econometrica 56(6), 1397–1417 | p. 1397 example; Assumption 1(b), p. 1400; Theorem 2.1, p. 1402; extension in §3 | Original PDF retained from the earlier audit was re-examined. Nonzero mean differences and covariance assumptions are explicit. A new web retrieval failed; the retained original, not a secondary summary, supports the discussion. |
| Dolado, Jenkinson and Sosvilla-Rivero (1990), *Cointegration and Unit Roots: A Survey*, Banco de España WP 9005 | Printed pp. 14–15, PDF pp. 15–16 | Original retained working paper re-examined. Enders attributes his modification to DJS. The final journal article is bibliographically identified as *Journal of Economic Surveys* 4(3), 249–273, DOI 10.1111/j.1467-6419.1990.tb00088.x, but its full text was not verified here. A repository record was found; direct access returned 403. |

Sources: [Hamilton publication record](https://www.jstor.org/stable/j.ctv14jx6sm), [Hamilton text reproduction](https://djvu.online/file/mozWnvm1yoLAH), [Enders author-hosted manual](https://www.time-series.net/assets/docs/SupplementaryManual_all.117125921.pdf), [Wiley edition record](https://www.wiley.com/en-ie/applied-econometric-time-series-4th-edition-p-9781118808566), [Wiley Chapter 4 companion record](https://bcs.wiley.com/he-bcs/Books?action=chapter&bcsId=9288&chapterId=105761&itemId=1118808568), [Sims–Stock–Watson original](https://www.princeton.edu/~mwatson/papers/Sims_Stock_Watson_Ecta_1990.pdf), [West original copy](https://www3.nd.edu/~nmark/Climate/West_Econometrica1913104.pdf), [DJS original working paper](https://www.bde.es/f/webbde/SES/Secciones/Publicaciones/PublicacionesSeriadas/DocumentosTrabajo/90/Fich/dt_9005e.pdf), [DJS journal repository metadata](https://e-archivo.uc3m.es/entities/publication/48ab5dd8-4531-4cdd-b025-abf57ae435f2).

Underlying research includes Dickey–Fuller distribution theory and Fuller tables; Enders Table A explicitly credits Fuller (1996). Dickey and Fuller (1981), *Likelihood Ratio Statistics for Autoregressive Time Series with a Unit Root*, Econometrica 49(4), 1057–1072, is relevant background for joint statistics. Its full article was **not** successfully read here; an earlier locally saved file with a PDF suffix is HTML, not evidence of verification. Bibliography is corroborated by [Dickey's own retrospective](https://garfield.library.upenn.edu/classics1993/A1993MC56300001.pdf). The table entries used below were checked in Enders directly, not inferred from an unread paper.

Secondary descriptions found during searching were not used to supply missing branches. No verified original passage names a separate “Hamilton/Enders” algorithm.

## 3. Common notation and statistical experiments

Use the package parameterization

\[
\Delta y_t=a+\gamma t+b y_{t-1}+\sum_{j=1}^{p}\beta_j\Delta y_{t-j}+u_t.
\]

The unit-root hypothesis is b=0 versus b<0. Hamilton's levels coefficient is ρ, so b=ρ−1; Enders's manual calls the level coefficient γ, which is **our b**, and calls the difference-equation trend coefficient a₂, which is **our γ**. Confusing these symbols would reverse the meaning of several restrictions.

`ct` estimates a, γ and b freely; `c` imposes γ=0; `n` imposes a=γ=0. A rejection of b=0 is a unit-root rejection; a stationary interpretation additionally requires the remaining autoregressive roots to be stable. Neither a selected model code nor an individual pretest proves the true deterministic structure.

For a fixed p and N input levels, the implementation uses n=N−p−1 response observations `[p+1:N]`. Let L contain the p matched lagged differences. Define:

- Uct: regress Δy on `[y[-1], L, 1, t]`, q=p+3.
- Uc: regress Δy on `[y[-1], L, 1]`, q=p+2.
- Un: regress Δy on `[y[-1], L]`, q=p+1.
- Act: impose b=0, regress Δy on `[L, 1, t]`.
- Ac: impose b=γ=0, regress Δy on `[L, 1]`.

These matrix definitions are an explicit operational translation and code comparison; they are not an assertion that Hamilton prescribed this package's sample or lag policy.

## 4. What Hamilton does and does not establish

Hamilton's Case 3 uses a c-only regression under a genuinely nonzero fixed drift and obtains conventional asymptotic inference. Case 4 adds a time regressor while the true process remains a random walk with or without drift: its unit-root statistic is nonstandard. A joint test of ρ=1 and the trend coefficient zero is also nonstandard. At p. 501 he recommends a model plausible under both hypotheses, rather than enumerating a deterministic-pretest cascade. Section 17.7 extends the regression analysis to augmentation. [Hamilton, §17.4 and §17.7](https://djvu.online/file/mozWnvm1yoLAH).

The following is a **case-selection map of that theory**, not a newly discovered sequential algorithm:

| Maintained experiment | Regression | Reference for level-coefficient statistic | What remains unspecified as an automatic algorithm |
|---|---|---|---|
| Driftless unit root, no deterministic regressors | n | Nonstandard DF no-constant | No rule deciding from data that this is the correct experiment |
| Driftless unit root, intercept estimated | c | Nonstandard DF intercept case | No complete ct-to-c-to-n route |
| Fixed genuinely nonzero drift, no separate trend estimated | c | Conventional asymptotic normal | No source-certified pretest event establishing nonzero drift |
| Random walk with arbitrary drift, constant and time trend estimated | ct | Nonstandard DF trend case | No license to replace it with normal inference because levels visibly grow |

Thus a Hamilton-only reconstruction cannot honestly provide the requested initial model, deterministic tests, switch thresholds and stopping rules. Those elements are absent from the inspected passages, not silently filled in below. The complete tree in §5 is Enders's operational presentation.

### Our explanation of why the distinction matters

For p=0 and b=γ=0,

\[
y_{t-1}=y_0+a(t-1)+\sum_{s=1}^{t-1}u_s.
\]

In a c-only regression, fixed a≠0 gives a deterministic O(T) regressor component dominating the Oₚ(√T) innovation accumulation. In ct, partialling out `[1,t]` removes that deterministic component. The remaining integrated component explains why drift does not make the ct unit-root statistic conventional. This algebra is our explanatory derivation of the distinction, not an extra branch attributed to Hamilton.

Sims–Stock–Watson's §5 Example 2 makes the same essential contrast for an AR(2): conventional asymptotic coefficient behavior without the additional trend, nonstandard behavior when it is included. Their canonical-regressor results require inspecting the actual restrictions; a broad slogan that the presence of *any* deterministic regressor makes *all* inference conventional is unsafe. [SSW, pp. 123–129](https://www.princeton.edu/~mwatson/papers/Sims_Stock_Watson_Ecta_1990.pdf).

For augmentation with stable difference dynamics and γ=0, the long-run difference mean is a/(1−Σβⱼ), not generally a. West requires this mean to be nonzero and appropriate regularity/covariance estimation; his stationary-regressor block excludes a deterministic time trend. This does not justify a ct cascade under arbitrary drift. [West, pp. 1397–1405](https://www3.nd.edu/~nmark/Climate/West_Econometrica1913104.pdf).

Our further inference: fixed-nonzero-drift results are pointwise. For local drift of order T⁻¹ᐟ², deterministic and accumulated stochastic terms have comparable scale. Significance-based selection does not establish uniform validity near zero drift, or validity conditional on having passed that selection event.

## 5. Complete operational tree: Enders and the current implementation

Enders's §4.4 numbered steps specify tau testing, joint phi gates, b=0 auxiliary deterministic tests, conditional retests and simplification. He identifies the method as a modification of DJS. [Enders, PDF pp. 63–64](https://www.time-series.net/assets/docs/SupplementaryManual_all.117125921.pdf).

The table is a complete formalization checked against those steps and `_enders.py`. R means reject the stated null; NR means fail to reject. Start at E1. Sample and p remain common. No joint or auxiliary rejection is itself a unit-root rejection.

| Node | Fit, null and statistic | Reference / rejection direction | R branch | NR branch |
|---|---|---|---|---|
| E1 | Uct; b=0; τ₃= b̂/se(b̂) | Nonstandard DF ct, left | Stop: reject unit root; retained ct | E2 |
| E2 | Uct versus Ac; b=γ=0; φ₃ | Nonstandard joint phi, right | E3 | E5 |
| E3 | Act; γ=0; auxiliary trend t ratio | Conventional deterministic-coefficient reference, two-sided in package | E4 | E5 |
| E4 | Original Uct; b=0; same level t ratio as E1 | Conditional conventional reference, left in package | Stop: reject unit root, retained ct | Stop: do not reject, retained ct |
| E5 | Uc; b=0; τ₂ | Nonstandard DF c, left | Stop: reject unit root, retained c | E6 |
| E6 | Uc versus regression on L only; b=a=0; φ₁ | Nonstandard joint phi, right | E7 | E9 |
| E7 | Ac; a=0; auxiliary intercept t ratio | Conventional deterministic-coefficient reference, two-sided in package | E8 | E9 |
| E8 | Original Uc; b=0; same level t ratio as E5 | Conditional conventional reference, left in package | Stop: reject unit root, retained c | Stop: do not reject, retained c |
| E9 | Un; b=0; τ₁ | Nonstandard DF n, left | Stop: reject unit root, retained n | Stop: do not reject, retained n |

At E3, only b is constrained; γ remains estimated. At E7, b=0 is imposed and γ was already removed. E4/E8 return to the unrestricted level-coefficient regressions; they do **not** test for a unit root in the restricted auxiliary regression. A repeated fit is unnecessary because the statistic was already calculated at E1/E5.

### Statistic definitions and calibration ledger

With unrestricted fit U, two restrictions and a matched restricted fit R,

\[
\phi=\frac{(RSS_R-RSS_U)/2}{RSS_U/(n-q_U)}.
\]

For φ₃, remove `[y[-1],t]`, retaining intercept and L. For φ₁, remove `[y[-1],1]`, retaining L. At p=0 the latter restricted RSS is Δy′Δy. This F-shaped algebra does **not** permit conventional F quantiles. No φ₂ is used.

| Statistic | Actual package threshold at nominal α | Source status / interpretation |
|---|---|---|
| τ₁, τ₂, τ₃ | ARCH 8.0.0 ADF critical values for n/c/ct; strict `<` | Nonstandard fixed-specification component calibration; not a sequential p-value |
| φ₁, φ₃ | Enders Table B; strict `>` | Nonstandard joint calibration; finite grid rule is a package convention |
| Auxiliary γ t | `abs(t) > t_(1−α/2, n−p−2)` | Conventional reference after imposing b=0; finite Student df and two-sided choice explicitly operational |
| Auxiliary a t | `abs(t) > t_(1−α/2, n−p−1)` | Same distinction; imposing a unit root is a maintained restriction |
| Conditional b t in ct | `t < t_(α, n−p−3)` | Reuses Uct statistic; conventional-limit motivation is conditional on the maintained DGP |
| Conditional b t in c | `t < t_(α, n−p−2)` | Reuses Uc statistic; retaining an intercept by policy is insufficient justification |

The manual specifies conventional t inference, but does not uniquely prescribe the package's tail and finite-df conventions. Normal limits do not establish exact finite-sample Student distributions in an autoregression. An asymptotic-normal version would use zα and z₁₋α/₂ instead; at 5% these are approximately −1.645 and 1.960. That would be an explicit calibration variant, not evidence of a separate named Hamilton tree.

A useful source check is the 5% T=250 row: Enders Table A gives τ₁=−1.95, τ₂=−2.88, τ₃=−3.43; Table B gives φ₁=4.63 and φ₃=6.34. The package retains ARCH tau surfaces rather than rounded Table A values. It uses the largest phi grid T≤n, capped at 500, and refuses n<25. These numerical choices are not Hamilton prescriptions. [Enders Tables A/B, PDF pp. 122–123](https://www.time-series.net/assets/docs/SupplementaryManual_all.117125921.pdf).

### Flowchart cross-reference and stopping ambiguity

Figure A2.1, PDF p. 66, was visually checked. At the trend conditional branch it says “Go back to Step 1” and then instructs a t-distribution retest. At the intercept conditional branch it says “Go to Step 2” before another t retest. A literal restart would not yield a coherent new algorithm: unchanged data would repeat the same earlier tests, and numbered Step 2 is the trend joint test. Numbered prose instead specifies return to the relevant unrestricted regression and terminal conventional inference. This review follows that explicit prose, rather than inventing looping/re-estimation rules or silently treating the figure's cross-reference as a second procedure. An author correction could clarify the figure's numbering. [Enders §4.4 and Figure A2.1](https://www.time-series.net/assets/docs/SupplementaryManual_all.117125921.pdf).

### Terminal interpretation

A conditional ct nonrejection corresponds to a maintained unit-root model with a genuine trend in the **difference equation**, potentially quadratic mean growth in levels. With p=0, summation gives a t + γt(t+1)/2. A conditional c nonrejection corresponds to a maintained unit root with drift. These are descriptions of maintained models, not conclusions proved by sequential nonrejection.

A nonzero trend in levels does not imply γ≠0 in the difference equation. Nor does selecting ct prove trend stationarity, selecting c prove drift, or selecting n prove a zero deterministic component. The alternative must satisfy the remaining stability conditions.

## 6. Branch-by-branch comparison with our algorithms

This table records code behavior at the inspected commit, independently of any proposed Hamilton name.

| Event in §5 | Original Enders | Default Conservative | Conservative with allow_no_constant=True |
|---|---|---|---|
| Entry | ct; choose p once or use supplied p | Identical | Identical |
| E1 rejection | Stop, unit-root rejection | Identical | Identical |
| E1 NR → E2 | Evaluate φ₃ | Identical | Identical |
| E2 NR, or E3 NR | Move to c | Identical | Identical |
| E2 R and E3 R | E4; terminal conditional ct test | Identical | Identical |
| E5 rejection | Stop, unit-root rejection | Identical | Identical |
| E5 nonrejection | E6 | **Stop immediately: retain c, no unit-root rejection** | Identical to Original |
| E6 NR, or E7 NR | Move to n, E9 | Unreachable | Identical to Original |
| E6 R and E7 R | E8; terminal conditional c test | Unreachable | Identical to Original |
| E9 | Terminal n DF decision | Unreachable | Identical to Original |

Thus default Conservative removes both the n paths and an additional c-only conditional rejection path. It is not merely a prohibition on returning n. For matched configuration its rejection set is a subset of Original's, because every rejection before E5 is shared and Conservative adds no later rejection. This is a logical property of the code, not proof of nominal size or superior size-adjusted power.

The existing once-under-ct MAIC/MBIC/AIC/BIC lag policy, common final sample, finite-grid rule and equality handling are explicit package choices. The manual's fixed-lag GDP illustration does not establish any per-stage automatic lag selection. They must remain part of a reproducible strategy definition.

## 7. Are there genuinely distinct candidates?

| Meaning assigned to “Hamilton/Enders” | Assessment | Recommendation |
|---|---|---|
| The nine-node Enders cascade above | Operationally established presentation already implemented as Original | No duplicate function, alias or Monte Carlo row |
| Hamilton's maintained nonzero-drift case | A conditional inference experiment, without deterministic selection | Possible separately specified research benchmark, not a selector |
| Same cascade with normal rather than finite-Student conventional thresholds | Different finite-sample rejection set, same branch structure | Label calibration sensitivity explicitly if later authorized; do not attribute an unverified algorithm |
| Delete φ gates and test deterministic terms directly | Different tree if fully specified | Source-incomplete; do not infer restricted/unrestricted fit, threshold or stop from Hamilton |
| Test unrestricted deterministic coefficients using ordinary t/F because a constant is present | Neither the Enders tree nor generally justified by the verified theory | Do not implement under the proposed name |
| Simplify deterministic terms after a unit-root rejection | Additional post-decision specification task | No complete original procedure verified here; define separately if needed |

DJS's inspected working-paper paragraph motivates general-to-specific testing and conditional normal retesting but does not uniquely enumerate all the auxiliary, phi-gate and final stopping conventions. Its “and so on” progression cannot certify an additional full c/n tree. [DJS WP 9005, printed pp. 14–15](https://www.bde.es/f/webbde/SES/Secciones/Publicaciones/PublicacionesSeriadas/DocumentosTrabajo/90/Fich/dt_9005e.pdf).

No verified passage resolves the unnamed method into a genuinely new ready-to-code selector. This is a scoped finding about the inspected sources, not proof that no edition, lecture note or other original publication could contain another procedure.

## 8. Implementation specification and independent validation, conditional on future approval

The established operational object is already available. Its specification is §5 plus the package conventions in §6; preserve `enders()` and its history/calibration metadata. No architecture change or new public strategy is warranted by this review.

Before considering another algorithm, require a source-backed contract naming: initial model; all hypotheses; restricted versus unrestricted deterministic tests; critical-value table and edition; separate selection/unit-root levels if any; tails and equality; terminal nonrejection rules; lag/sample policy; and maintained DGP/covariance assumptions. Identify choices absent from the source as our research conventions. A new entry must differ in rejection set or substantive output, not merely its name.

If a later contract is approved, independent numerical validation should:

1. Build ct/c/n and b=0 auxiliary designs by explicit time indexing; solve QR independently of the production fit and compare coefficients, residuals, RSS, standard errors and t ratios.
2. Reconstruct φ by restricted/unrestricted RSS and independently by a two-restriction Wald calculation. Inspect retained columns, unrestricted residual df and named restrictions; validate critical lookup independently of implementation constants.
3. Exercise every R/NR edge, including terminal conditional **nonrejection**, exact threshold equality, p=0 and positive p, and identical response samples. Test that unit-root testing returns to the unrestricted fit rather than its auxiliary.
4. Prove full-history parity with Original for any purported alternative presentation. For an actual new tree, use controlled statistics/fixtures on paths where its decision differs; do not count duplicated outcomes as independent algorithms.
5. Reuse the GDP fixture as a historical regression check. The saved evidence distinguishes Original rejection from default Conservative nonrejection and documents a printed intercept-t rounding discrepancy. Neither the example nor software parity validates nominal size.
6. Verify calibration separately from arithmetic. Fixed-nonzero-drift normality cannot be “validated” just by matching ordinary regression software. A maintained normal benchmark needs a source-specific covariance estimator and null; no deterministic-selection step should be implied.

These are plans only. No algorithm or tests were written and no validation run was performed for this review.

## 9. Implications for the planned Monte Carlo study

**Do not add a Hamilton/Enders-labelled row until it has a distinct, complete definition.** Begin with the two distinct existing Enders trees, prespecified ADF specifications, and the experimental AIC/BIC selectors. Treat Conservative with permission as a parity check, not a third independent Enders algorithm. IC selection followed by fixed-specification inference remains an uncalibrated complete testing pipeline; specify that pipeline explicitly for research.

A later protocol should separate:

- b=0, γ=0 with zero, small/local and fixed nonzero long-run drift;
- b=0, genuine γ≠0, allowing quadratic level mean growth, as a separately labelled expanded null;
- stationary and local-to-unity alternatives around zero/nonzero means or linear deterministic trends, with other roots stable;
- fixed common p first, then any separately specified automatic-lag policy.

Use the same innovations/data across methods within replications for paired comparison. Record unconditional unit-root rejection, selection and branch frequencies, terminal calibration, conditional descriptive rejection rates, numerical failures, sample/lag choices and Monte Carlo uncertainty. Report raw and size-adjusted power with independent calibration/evaluation samples if adjustment is used. Do not discard rare branches or failed fits without reporting them.

Include sensitivity to normal versus finite-Student conventional thresholds and phi table-row policy as **calibration variants**, not new historical names. A prespecified nonzero-drift normal benchmark may explain pointwise behavior but must not serve as a universal oracle at zero/local drift. Local drift sequences and initial conditions must be fixed in the protocol before simulation.

The most consequential theoretical question is the size of the **entire adaptive rejection rule**, especially near deterministic-coefficient boundaries. Passing branches changes the joint distribution of subsequent statistics; asymptotic validity under maintained parameters is not validity conditional on the selection event. Enders himself warns about cascading inference in the GDP discussion. [Enders, PDF pp. 64–65](https://www.time-series.net/assets/docs/SupplementaryManual_all.117125921.pdf).

No simulations or bootstrap work are authorized or performed here.

## 10. Outstanding questions and disposition

1. **Identity of the unnamed procedure:** no exact external page or algorithm supplied beyond the Hamilton/Enders label. Another source could change the classification; a synthesized cascade cannot be treated as verified publication.
2. **Hamilton page images:** p. 501 and Cases 3–4 are verified through reproduced original text, not clean images. Obtain a legitimate page-image copy before claiming bibliographic certification beyond that level.
3. **Enders figure numbering:** reconcile the conditional-branch cross-references with the explicit numbered prose; do not implement literal loops. Main-book treatment in the relevant edition remains unverified here.
4. **Conventional calibration:** tails, finite versus asymptotic thresholds, covariance assumptions and near-zero deterministic parameters remain methodological choices/limitations, not solved by significance pretests.
5. **Historical provenance:** Enders identifies his tree as a DJS modification; final DJS journal operational details were not read here. Do not assert exact identity with every DJS/Perron presentation.
6. **Whole-procedure performance:** neither Hamilton's conditional theory nor Enders's published example establishes uniform size, empirical size for our conventions, or selective p-values.

**Disposition:** preserve Original and Conservative Enders, add no Hamilton/Enders implementation now, and carry the explicit branch/calibration distinctions into the separately approved Monte Carlo protocol. The review is complete at the attainable source-verification level; the listed limitations remain visible rather than being filled by assumptions.
