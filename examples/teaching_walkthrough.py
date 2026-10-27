"""Reproducible teaching example; paste ordinary Python cells into Colab."""
import numpy as np
from unitroot_tools import adf, dfgls, pp, kpss, select_lags

y = np.random.default_rng(20261027).standard_normal(240).cumsum()
results = (
    adf(y, trend="c", lags=4),
    dfgls(y, trend="ct", lags=0),
    adf(y, trend="c", max_lags=12),
    dfgls(y, trend="ct", method="mbic", max_lags=12),
    pp(y, trend="c", bandwidth=4),
    pp(y, trend="ct", test_type="rho"),
    kpss(y, trend="c", bandwidth=0),
    kpss(y, trend="ct"),
)
for result in results:
    print(result.summary())
    print("Reject stated null at 5%:", result.reject(.05))
    print()
for method in ("maic", "mbic", "aic", "bic"):
    selection = select_lags(y, trend="c", method=method, max_lags=12)
    print(method, "selected:", selection.selected_lags,
          "common sample:", selection.selection_nobs)
    for row in selection.candidates:
        print(row.lags, row.aic, row.bic, row.maic, row.mbic)
print("ADF regression:", results[0].regression.params)
print("DF-GLS regression:", results[1].regression.params)
print("PP covariance:", results[4].regression.cov_type)
assert results[-1].regression is None
