"""Run with: python examples/adf_dfgls.py (after installing unitroot-tools)."""
import numpy as np

from unitroot_tools import adf, dfgls


def main():
    y = np.random.default_rng(20261027).standard_normal(240).cumsum()
    fixed = adf(y, trend="c", lags=4)
    automatic = dfgls(y, trend="ct", max_lags=12)  # default MAIC
    conventional = adf(y, trend="c", method="bic", max_lags=12)
    for result in (fixed, automatic, conventional):
        print(result.summary())
        print("Reject at 1%, 5%, 10%:", [result.reject(a) for a in (.01, .05, .10)])
        print()
    # Public regression access is available for inspection. Its ordinary
    # coefficient p-values are not the unit-root p-value returned above.
    print("DF-GLS regression coefficients:", automatic.regression.params)


if __name__ == "__main__":
    main()
