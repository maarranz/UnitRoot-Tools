# Sequential Enders procedures

```python
from unitroot_tools import enders

original = enders(y, lags=1)
conservative = enders(y, variant="conservative", max_lags=12)
permitted = enders(y, variant="conservative", allow_no_constant=True, max_lags=12)
print(conservative.summary())
for step in original.history:
    print(step.test, step.statistic, step.critical_value, step.reject, step.next_action)
```

## API

```python
enders(y, *, variant="original", alpha=0.05, lags=None,
       max_lags=None, method="maic", allow_no_constant=False)
```

`alpha` accepts only 0.01, 0.05 and 0.10. Data must be finite, real,
one-dimensional and equally spaced; gaps are not dropped. A minimum of 25
**effective regression observations** is required because the phi table starts
at 25. Automatic lag selection that leaves fewer observations raises an error;
reduce `max_lags` explicitly rather than silently changing the selected order.

`method` accepts MAIC, MBIC, AIC and BIC (case-insensitive). Selection runs once
under `ct`; the selected augmentation order is then fixed. Every regression
uses the same response slice `[p+1:N]`, the same lagged differences and, when
present, the ARCH trend index `1,...,N-p-1`. Selection itself retains the
existing maximum-lag common sample; it need not equal the final sample.

With fixed `lags`, the default `method="maic"` and `None` are inactive. An
explicit nondefault method or `max_lags` raises an error. This resolves the
proposed signature's default/fixed-lag conflict without changing `adf()`.
For `variant="original"`, `allow_no_constant` is ignored and effective permission
is always true. Both the requested and effective flags are recorded.

## Reference and implementation conventions

Primary reference: Walter Enders, *Supplementary Manual to accompany Applied
Econometric Time Series*, fourth edition,
[Section 4.4, PDF pages 63–66; Table B, PDF page 123](https://www.time-series.net/assets/docs/SupplementaryManual_all.117125921.pdf).
The numbered text verifies both auxiliary regressions and subsequent
conventional t-tests. The GDP illustration also demonstrates why the
resulting sequential conclusion can be economically implausible.

The common-lag choice is our explicit implementation convention, compatible
with the fixed-lag illustration; the reference does not prescribe MAIC or a
per-stage automatic-lag algorithm. We do not describe this MAIC extension as
part of Enders's original publication.

Three calibrations remain separate:

- **tau1/tau2/tau3:** existing `adf()` statistics and ARCH 8.0.0 critical values;
  reject strictly in the left tail. These response surfaces need not equal
  rounded textbook tables.
- **phi1/phi3:** nested-regression RSS statistics with two restrictions and
  unrestricted residual degrees of freedom. Reject strictly in the right
  tail using Enders Table B; never use ordinary F quantiles. The chosen finite
  table row is the largest T not exceeding effective regression n, capped at
  T=500. This is an explicit finite-grid approximation, not interpolation or
  an exact finite-sample calibration. Below n=25 we refuse inference. Every
  phi history entry records the table row, restricted RSS and retained columns.
- **Conventional tests:** auxiliary deterministic coefficients use two-sided
  Student-t quantiles with their own residual degrees of freedom. The terminal
  conditional level-coefficient test uses the original unrestricted ADF
  statistic and its residual degrees of freedom, in the left tail, matching
  the assignment's alternative b<0. The manual specifies a t-distribution but
  not explicitly its tail. The finite-df and tail conventions are therefore
  exposed here rather than attributed to a uniquely specified published rule.

The inspected [urca source](https://github.com/cran/urca/blob/master/R/ur-df.R)
contains different T=250 phi3 entries at 5% and 10% (6.49 and 5.47), whereas
Enders's primary table gives 6.34 and 5.39. We transcribed the primary table
independently; no urca source code was copied. Our sample-size row rule also
differs from urca's rule. No phi2 enters either tree.

## Decisions and interpretation

Original starts with `ct`. A tau rejection terminates immediately. Otherwise
a phi rejection leads to the b=0 auxiliary regression. A significant auxiliary
deterministic coefficient leads to the conditional t-test and termination;
otherwise the procedure moves to the simpler specification. At `n`, tau1 is
terminal. Conservative follows the same trend stage, but by default stops at
a constant-only tau2 nonrejection. Setting `allow_no_constant=True` enables
the remaining original branches; it does not force selection of `n`.

**Conservative Enders is a proposed modification, not an established or
independently size-validated test.** Retaining an intercept by policy never
licenses ordinary t inference. A joint rejection is never a unit-root
rejection; an auxiliary rejection is never a unit-root rejection either.

The conventional branches presume that the retained deterministic component
really belongs in the difference-equation DGP, with adequate augmentation and
regular errors. Statistical preselection does not establish those premises.
Under b=0, a nonzero constant implies drift in levels; a nonzero difference-
equation trend can imply a quadratic level path. Do not interpret a conditional
nonrejection as proof of these features or of a unit root. Neither procedure
is robust by construction to breaks, outliers or arbitrary volatility changes.

Economic priors, plots and residual diagnostics should guide whether a
sequence starting at `ct` is appropriate at all. Do not choose `n` because it
produces a smaller p-value. The nominal alpha is per step, not a guarantee of
the whole procedure's size. Regressing the level on deterministic terms under
stationarity can be a descriptive rule of thumb; it is not an automatic stage
of these procedures.

## Results

`EndersResult` exposes `variant`, `alpha`, effective `allow_no_constant`,
`allow_no_constant_requested`, selected `trend`, `lags`, `lag_method`,
`lag_selection`, `reject_unit_root`, `terminal_test`, `statistic`,
`critical_value`, sample boundaries, dependency `versions`, `notes` and
chronological `history`. No post-selection p-value is advertised.

Each `EndersStep` records specification, named test, hypothesis, signed
statistic, critical value, tail, rejection, calibration, branch/termination
reason and fitted regression. Joint steps additionally retain restricted
RSS/columns, q and table sample size; conventional steps retain reference df.
For two-sided tests compare the absolute statistic with the positive critical
value. Scalar records are frozen; exposed Statsmodels regression objects are
not deeply immutable. History includes only steps actually executed.

## Reproduction and subsequent research

Run `python -m pytest` and
`python validation/validate_enders.py --report /tmp/enders-validation.json`.
The GDP fixture is the DATE/RGDP subset of the author's
[RGDP workbook](https://time-series.net/assets/docs/rgdp.350347.xls);
provenance and hashes are in `validation/enders/source.json`.
No workbook-reading dependency is required for normal installation or tests.

A subsequent, separately authorized Monte Carlo study should compare original,
conservative (both flag settings), and fixed-specification ADF with matched
lag/sample conventions. Vary sample size, local-to-unity alternatives, zero and
small drift/trend coefficients, and short-run MA dynamics. Report size, raw
power, size-adjusted power, branch frequencies and specification frequencies,
with simulation uncertainty. Include sensitivity to table-row and conventional
critical-value conventions. Do not describe this numerical implementation
validation as that study.
