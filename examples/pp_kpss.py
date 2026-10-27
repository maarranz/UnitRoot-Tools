"""Run after installing UnitRoot-Tools and its released ARCH dependency."""
import numpy as np
from unitroot_tools import pp, kpss

series = np.random.default_rng(20261027).standard_normal(240).cumsum()
for result in (
    pp(series),
    pp(series, trend="ct", bandwidth=4, test_type="rho"),
    kpss(series, bandwidth=0),
    kpss(series, trend="ct"),
):
    print(result.summary())
    print("Bandwidth:", result.bandwidth, result.bandwidth_method)
    print("Reject stated null at 5%:", result.reject(.05))
    print()
