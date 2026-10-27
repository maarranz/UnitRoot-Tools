from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
import math

import numpy as np
import pytest
from numpy.testing import assert_allclose
from statsmodels.regression.linear_model import OLS

import unitroot_tools
from unitroot_tools import EstimationError, SpecificationWarning, select_lags
from unitroot_tools.selection import _choose

Y = np.array([1., 2., 1., 3., 2., 5.])


def test_hand_calculated_fraction_example():
    # N=6, K=1, response [-1,2,-1,3], levels [2,1,3,2].
    # X'X at p=1 is [[18,5],[5,7]], X'y=[3,-8].
    expected = [(Fraction(1,6), Fraction(29,2), Fraction(4,29)),
                (Fraction(61,101), Fraction(60,101), Fraction(22326,505))]
    result = select_lags(Y, trend="n", max_lags=1)
    for row, (b, ssr, tau) in zip(result.candidates, expected, strict=True):
        assert row.nobs == 4
        assert_allclose([row.level_coefficient, row.ssr, row.variance, row.level_sum_squares, row.tau],
                        [float(b), float(ssr), float(ssr)/4, 18, float(tau)], rtol=1e-13)
        assert_allclose(row.maic-row.aic, 2*float(tau)/4)
        assert_allclose(row.mbic-row.bic, math.log(4)*float(tau)/4)
    assert result.selected_lags == 0
    assert select_lags(Y, trend="n", max_lags=1, method="aic").selected_lags == 1
    assert select_lags(Y, trend="n", max_lags=1, method="bic").selected_lags == 1
    assert select_lags(Y, trend="n", max_lags=1, method="mbic").selected_lags == 0


@pytest.mark.parametrize("trend", ["n", "c", "ct"])
def test_conventional_criteria_against_matched_public_ols(trend):
    N, K = 80, 5
    y = np.random.default_rng(43).normal(size=N).cumsum()
    result = select_lags(y, trend=trend, max_lags=K, method="bic")
    if trend == "n":
        v = y
    else:
        z = np.ones((N,1)) if trend == "c" else np.column_stack((np.ones(N), np.arange(1,N+1)))
        v = OLS(y, z).fit().resid
    dy = np.diff(v); n = N-1-K
    for row in result.candidates:
        x = np.column_stack([v[K:-1]]+[dy[K-j:-j] for j in range(1,row.lags+1)])
        fit = OLS(dy[K:], x).fit()
        # Statsmodels counts p+1 fitted coefficients, and Gaussian constants.
        common = math.log(2*math.pi)+1
        assert_allclose(row.aic, fit.aic/n-common-2/n, rtol=0, atol=1e-12)
        assert_allclose(row.bic, fit.bic/n-common-math.log(n)/n, rtol=0, atol=1e-12)
        assert_allclose(row.ssr, fit.ssr, rtol=0, atol=1e-10)
        assert_allclose(row.variance, fit.ssr/n, rtol=0, atol=1e-12)


@pytest.mark.parametrize("method", ["maic", "mbic", "aic", "bic"])
@pytest.mark.parametrize("trend", ["n", "c", "ct"])
@pytest.mark.parametrize("factor", [1e-100, -7., 1e100])
def test_units_and_sign_invariance(method, trend, factor):
    y = np.random.default_rng(25).normal(size=60).cumsum()
    a = select_lags(y, trend=trend, max_lags=4, method=method)
    b = select_lags(y*factor, trend=trend, max_lags=4, method=method)
    assert a.selected_lags == b.selected_lags
    for x,z in zip(a.candidates,b.candidates, strict=True):
        assert_allclose([z.level_coefficient,z.tau], [x.level_coefficient,x.tau], rtol=1e-12, atol=1e-12)
        assert_allclose(z.variance/(factor**2), x.variance, rtol=2e-13)
        assert_allclose(z.criterion_value-x.criterion_value, 2*math.log(abs(factor)), rtol=0, atol=2e-12)


@pytest.mark.parametrize("trend", ["c", "ct"])
@pytest.mark.parametrize("method", ["maic", "mbic", "aic", "bic"])
def test_deterministic_shift_invariance(trend, method):
    y = np.random.default_rng(7).normal(size=70)
    added = np.full(70, 37.)
    if trend == "ct":added += .5*np.arange(1,71)
    a = select_lags(y, trend=trend, max_lags=5, method=method)
    b = select_lags(y+added, trend=trend, max_lags=5, method=method)
    assert a.selected_lags == b.selected_lags
    assert_allclose([x.criterion_value for x in a.candidates], [x.criterion_value for x in b.candidates], rtol=0,atol=1e-10)


def test_zero_order_and_metadata():
    r = select_lags(Y, trend="n", max_lags=np.int64(0), method="BIC")
    assert r.max_lags_requested == r.max_lags == r.selected_lags == 0
    assert r.uncapped_rule_value is None and r.runner_up_gap is None
    assert r.selection_nobs == 5 and (r.sample_start,r.sample_stop)==(1,6)
    assert r.method=="bic" and r.selection_detrending=="none"
    assert r.criterion_normalization=="per_observation_omit_common_constants"
    assert r.notes==() and r.input_name is None


def test_defaults_caps_and_warning():
    with pytest.warns(SpecificationWarning,match="capped"):
        r = select_lags(Y, trend="c")
    assert r.max_lags==1 and r.max_lags_requested is None
    assert r.uncapped_rule_value==math.ceil(12*(6/100)**.25)
    assert "capped" in r.notes[0] and r.selection_detrending=="ols"
    y=np.random.default_rng(1).normal(size=200)
    r=select_lags(y)
    assert r.max_lags==r.uncapped_rule_value==15
    assert r.notes==() and r.method=="maic"


@pytest.mark.parametrize("trend,N", [("n",3),("c",3),("ct",4)])
def test_shortest_sample(trend,N):
    y=np.array([1.,4.,2.,8.])[:N]
    r=select_lags(y,trend=trend,max_lags=0)
    assert r.selection_nobs==N-1
    with pytest.raises(ValueError,match="observations"):
        select_lags(y[:-1],trend=trend,max_lags=0)


def test_immutability_and_input_copy():
    y=Y.copy();r=select_lags(y,trend="n",max_lags=1)
    with pytest.raises(FrozenInstanceError):r.selected_lags=99
    with pytest.raises(FrozenInstanceError):r.candidates[0].tau=0
    with pytest.raises(TypeError):r.candidates[0]=r.candidates[1]
    y[:]=0
    assert r.candidates[0].level_sum_squares==pytest.approx(18)


def test_exact_tie_and_near_tie_without_tolerance():
    # Isolate ranking contract; do not rely on accidental floating-point ties.
    row=select_lags(Y,trend="n",max_lags=1).candidates[0]
    a=replace(row,lags=0,criterion_value=1.)
    b=replace(row,lags=1,criterion_value=1.)
    assert _choose((b,a))==(a,0.)
    b=replace(b,criterion_value=np.nextafter(1.,0.))
    winner,gap=_choose((a,b))
    assert winner.lags==1 and gap>0


def test_pandas_series_name_without_mutation():
    import pandas as pd
    y=pd.Series(Y,name="inflation");r=select_lags(y,max_lags=1,trend="C")
    assert r.input_name=="inflation" and r.trend=="c"
    assert_allclose(y.to_numpy(),Y)


def test_stage_exports_and_public_backend_boundary():
    import ast
    from pathlib import Path
    assert all(hasattr(unitroot_tools,x) for x in ("adf","dfgls","pp","kpss","UnitRootResult"))
    assert set(unitroot_tools.__all__)=={"select_lags","adf","dfgls","pp","kpss","UnitRootResult","LagSelectionResult","CandidateDiagnostics","EstimationError","SpecificationWarning"}
    for path in Path(unitroot_tools.__file__).parent.glob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node,ast.Import):
                for n in node.names:
                    if n.name.startswith("arch"):
                        assert path.name=="_adapters.py" and n.name=="arch"
            if isinstance(node,ast.ImportFrom) and (node.module or "").startswith("arch"):
                assert path.name=="_adapters.py" and node.module=="arch.unitroot"
                assert {n.name for n in node.names} <= {"ADF","DFGLS","PhillipsPerron","KPSS"}
