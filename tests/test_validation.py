import numpy as np
import pytest

from unitroot_tools import EstimationError, SpecificationWarning, select_lags


@pytest.mark.parametrize("y", [1., [[1],[2],[3]], [[1,2],[3,4]], True, [1,True,3],
                                    [True,False,True], ["1","2","3"], [1+1j,2,3],
                                    np.array([1,2,3],dtype=object), [None,2,3], [[1],[2,3]]])
def test_invalid_types(y):
    with pytest.raises(TypeError):select_lags(y,max_lags=0)


@pytest.mark.parametrize("y", [[], [1,np.nan,2], [1,np.inf,2], [1,-np.inf,2]])
def test_invalid_numeric_data(y):
    with pytest.raises(ValueError):select_lags(y,max_lags=0)


@pytest.mark.parametrize("name,value,error", [("method","t-stat",ValueError),("method","unknown",ValueError),
    ("method",None,TypeError),("trend","ctt",ValueError),("trend",2,TypeError),
    ("max_lags",True,TypeError),("max_lags",np.bool_(False),TypeError),
    ("max_lags",1.0,TypeError),("max_lags",-1,ValueError),("max_lags",100,ValueError)])
def test_invalid_options(name,value,error):
    kwargs={"max_lags":0,name:value}
    with pytest.raises(error):select_lags([1.,4.,2.,3.,8.],**kwargs)


@pytest.mark.parametrize("trend,y", [("n",np.zeros(8)),("n",np.ones(8)),("c",np.ones(8)),
                                      ("ct",2+.4*np.arange(10)),("n",2.**np.arange(10))])
def test_deterministic_and_zero_variance(trend,y):
    with pytest.raises(EstimationError):select_lags(y,trend=trend,max_lags=0)


def test_zero_candidate_regressor_and_singular_search():
    # N=5,K=1: matched levels at indices 1:4 are zero.
    with pytest.raises(EstimationError,match="range|rank"):
        select_lags([1.,0.,0.,0.,2.],trend="n",max_lags=1)
    # All lagged levels are geometric, but the last response breaks prediction:
    # p=0 has positive SSR; at p=1 its two regressor columns are proportional.
    y=2.**np.arange(10);y[-1]+=1
    with pytest.raises(EstimationError,match="rank"):
        select_lags(y,trend="n",max_lags=1)


def test_poor_condition_is_warned_and_recorded():
    y=1.1**np.arange(20)+np.random.default_rng(4).normal(scale=1e-9,size=20)
    with pytest.warns(SpecificationWarning,match="conditioned"):
        r=select_lags(y,trend="n",max_lags=1)
    assert any("conditioned" in note for note in r.notes)


@pytest.mark.parametrize("scale", [1e-250,1e250])
def test_unrepresentable_original_unit_diagnostics(scale):
    with pytest.raises(EstimationError,match="rescale"):
        select_lags(np.array([1.,4.,2.,3.,8.])*scale,trend="n",max_lags=0)


def test_nearly_deterministic_but_resolvable():
    y=np.arange(30.)+np.random.default_rng(9).normal(scale=1e-8,size=30)
    r=select_lags(y,trend="ct",max_lags=0)
    assert r.candidates[0].variance>0


def test_sample_cap_does_not_skip_unestimable_candidate():
    # The sample-length cap is necessary, not sufficient, for a valid search.
    y=np.arange(1.,9.)
    with pytest.warns(SpecificationWarning,match="capped"):
        with pytest.raises(EstimationError):select_lags(y,trend="c")


def test_backend_independence(monkeypatch):
    import builtins
    original=builtins.__import__
    def guarded(name,*args,**kwargs):
        if name.startswith("arch"):raise AssertionError("ARCH must not be called")
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,"__import__",guarded)
    assert select_lags([1.,4.,2.,3.,8.],trend="n",max_lags=0).selected_lags==0


def test_linalg_failure_preserves_cause(monkeypatch):
    from unitroot_tools import selection
    failure=np.linalg.LinAlgError("SVD did not converge")
    def fail(*args, **kwargs):raise failure
    monkeypatch.setattr(selection.linalg,"lstsq",fail)
    with pytest.raises(EstimationError,match="SVD least squares failed") as caught:
        select_lags([1.,4.,2.,3.,8.],trend="n",max_lags=0)
    assert caught.value.__cause__ is failure


def test_lag_column_zero_is_infeasible_without_skipping():
    # Only the final response moves; every lagged difference is exactly zero.
    y=np.ones(10);y[-1]=2
    with pytest.raises(EstimationError,match="zero regressor column"):
        select_lags(y,trend="n",max_lags=1)
