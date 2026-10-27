# UnitRoot-Tools

A small teaching API for four unit-root/stationarity tests and transparent
lag selection. This is an unpublished release candidate; the version is
`0.1.0rc1`. No PyPI or GitHub release is assumed.

## Install and learn

From a local checkout, in an activated Python virtual environment:

```sh
python -m pip install .
python -m pip check
python examples/teaching_walkthrough.py
```

Use `python -m pip install -e '.[test]'` for development, or install a supplied wheel with
`python -m pip install /path/to/unitroot_tools-0.1.0rc1-py3-none-any.whl`. See the [user guide](docs/User_Guide.md) for complete
examples, candidate diagnostics and Google Colab adaptation, and the
[candidate review](docs/Stage4B_Release_Candidate.md) for tested compatibility.

## Tests and interpretation

ADF, DF-GLS and PP test a unit-root null; KPSS tests stationarity around the
stated deterministic terms. Their rejection tails differ. Failure to reject
proves neither stationarity nor a unit root. These tests do not automatically
select a model or determine whether to difference a series.


Stage 3C provides independent lag selection and public `adf()`, `dfgls()`,
`pp()` and `kpss()` functions backed by **released ARCH 8.0.0**.

```python
import numpy as np
from unitroot_tools import adf, dfgls, select_lags

y = np.random.default_rng(20261027).standard_normal(240).cumsum()
fixed = adf(y, trend="c", lags=4)          # explicit zero is also supported
print(fixed.summary())
print(fixed.reject(0.05))                 # uses ARCH's 5% critical value

automatic = dfgls(y, trend="ct", max_lags=12)  # default MAIC
print(automatic.lags, automatic.nobs)
print(automatic.lag_selection.selection_nobs)
print(automatic.regression.params)       # public Statsmodels regression object

conventional = adf(y, method="bic", max_lags=12)
selection = select_lags(y, max_lags=12, method="mbic")
for row in selection.candidates:
    print(row.lags, row.aic, row.bic, row.maic, row.mbic)
```

Run the short executable examples with `python examples/adf_dfgls.py` and
`python examples/pp_kpss.py`.
Install first with `python -m pip install -e '.[test]'`.

`maic`, `mbic`, `aic` and `bic` use complete-series OLS detrending and a
common maximum-lag selection sample. Conventional scores omit the modified
tau penalty. These two-step conventions need not match one-step ADF autolag
from another library. MAIC is the automatic default for both tests.

ADF supports `n`, `c`, `ct`; `ctt` requires fixed lags. DF-GLS supports `c`
and `ct`. Fixed `lags` cannot be combined with `method` or `max_lags`.
Automatic maximum orders may be capped for final-regression feasibility;
this emits a warning and records the effective maximum in the result.
The selection sample and final regression sample are reported separately.

The backend receives original validated levels and an explicit integer lag.
DF-GLS performs its GLS detrending exactly once inside ARCH. Its returned
regression uses GLS-detrended data, while the selection regressions use
OLS-detrended data. The exposed regression is a mutable backend object;
its ordinary coefficient p-values are not the unit-root p-value.

Results and metadata are frozen; critical-value and version maps are read-only
copies. `reject()` uses the strict left-tail comparison for ADF/DF-GLS/PP and
right-tail comparison for KPSS at 1%, 5% or 10%.
Exact equality does not reject. Approximate p-values remain visible but never
override the critical-value decision; failure to reject does not establish H0.
ARCH calibration is preserved, including its dedicated DFGLS surfaces.

Missing observations are rejected, never dropped. Data must be equally spaced
for the model; this vector interface does not infer temporal frequency.
Singular/deterministic regressions and unsuitable numerical scales fail with
`EstimationError`. Unexpected programming errors are not silently swallowed.

Run tests with `python -m pytest --cov=unitroot_tools --cov-branch`.
Python 3.11+ is the declared range; see the [Stage 4 matrix](docs/Stage4_Release_Review.md)
for environments actually tested. Compatibility is not inferred from metadata.

The selector remains independent NumPy/SciPy code. Adapters use public ARCH
APIs only, check the active version, and need no R or development snapshot.
See the [governing specification](docs/Implementation_Specification.md).

PP supports `trend="n"`, `"c"`, `"ct"` and `test_type="tau"` (default) or
`"rho"`. KPSS supports `"c"` (default) and `"ct"`. Both accept `bandwidth=0`
and other fixed nonnegative integers; they do not accept augmentation `lags`,
`method` or `max_lags` arguments.

```python
from unitroot_tools import pp, kpss

p = pp(y, trend="ct", bandwidth=4, test_type="rho")
k = kpss(y, trend="c")  # public ARCH Hobijn-style automatic bandwidth
print(p.summary())
print(k.summary())
print(k.bandwidth, k.bandwidth_method, k.reject(0.05))
assert p.lags is None and k.regression is None
```

With N input observations, PP's automatic bandwidth is
`min(ceil(12*(N/100)**0.25), N-2)`, with a warning if capped. PP requires
`bandwidth < N-1`; KPSS requires `bandwidth < N`. KPSS delegates its automatic
choice to public ARCH and rejects invalid resolved bandwidths instead of
silently capping them. Very short samples can meet the length check yet fail
ARCH's automatic calculation; specify a feasible fixed bandwidth or use more
observations. PP uses N-1 regression observations; KPSS uses all N levels.
Bandwidth metadata is separate from augmentation-lag metadata.

PP tests the unit-root null and exposes ARCH's public HAC regression result.
KPSS tests stationarity around an intercept or a linear trend; rejection does
not by itself prove a unit root. KPSS exposes no public regression. Approximate
p-values and endpoint notes retain ARCH's own calibration, which can differ
from R urca and Statsmodels. Only ARCH's exact historical KPSS default-change
deprecation notice is suppressed; other warnings remain visible.
