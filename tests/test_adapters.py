"""Released public ARCH integration and preserved Stage 2 evidence."""
import math
from types import SimpleNamespace

import arch
from arch.unitroot import ADF, DFGLS
import numpy as np
import pytest
from numpy.testing import assert_allclose

from unitroot_tools import EstimationError, SpecificationWarning, adf, dfgls, select_lags

FUNCTIONS = [(adf, ADF, "ADF", ["n","c","ct","ctt"]), (dfgls, DFGLS, "DFGLS", ["c","ct"])]
FIXED_CASES = [(f, cls, test, tr) for f,cls,test,trends in FUNCTIONS for tr in trends]
AUTO_CASES = [x for x in FIXED_CASES if x[-1]!="ctt"]
SERIES=["rw","rw_negma","ar07","trend_ar08","hetero_rw","mic_boundary"]


def check_inference(result, expected):
    assert_allclose(result.statistic, expected.stat, rtol=0, atol=1e-12)
    assert_allclose(result.pvalue, expected.pvalue, rtol=0, atol=1e-12)
    assert result.nobs==int(expected.nobs)
    assert dict(result.critical_values)=={.01:expected.critical_values["1%"],.05:expected.critical_values["5%"],.10:expected.critical_values["10%"]}
    assert result.regression.nobs==result.nobs
    for alpha in (.01,.05,.10):assert result.reject(alpha)==bool(expected.stat < result.critical_values[alpha])


@pytest.mark.parametrize("name", SERIES)
@pytest.mark.parametrize("f,cls,test,trend", FIXED_CASES)
@pytest.mark.parametrize("p", [0,4])
def test_fixed_reference_evidence(stage2,name,f,cls,test,trend,p):
    y=stage2[0][name]
    result=f(y,trend=trend,lags=p)
    check_inference(result,cls(y,trend=trend,lags=p))
    assert result.lags==p and result.lag_selection is None
    assert (result.input_nobs,result.nobs,result.sample_start,result.sample_stop)==(len(y),len(y)-1-p,p+1,len(y))
    # Read the saved release rows directly, including their calibration.
    from conftest import RESULTS
    import csv
    with (RESULTS/'python_release.csv').open() as file:
        saved=[r for r in csv.DictReader(file) if (r['series'],r['backend'],r['test'],r['trend'],r['lags'])==
               (name,'ARCH',test,trend,str(p))][0]
    assert_allclose([result.statistic,result.pvalue,*result.critical_values.values()],
                    [float(saved[k]) for k in ['statistic','pvalue','cv1','cv5','cv10']],rtol=0,atol=1e-8)


@pytest.mark.parametrize("name", SERIES)
@pytest.mark.parametrize("f,cls,test,trend", AUTO_CASES)
@pytest.mark.parametrize("method", ["maic","mbic","aic","bic"])
@pytest.mark.parametrize("K", [4,12])
def test_automatic_public_fixed_lag_composition(stage2,name,f,cls,test,trend,method,K):
    y=stage2[0][name]
    selection=select_lags(y,trend=trend,max_lags=K,method=method)
    result=f(y,trend=trend,max_lags=K,method=method)
    assert result.lag_selection==selection and result.lags==selection.selected_lags
    assert result.nobs==len(y)-1-result.lags
    check_inference(result,cls(y,trend=trend,lags=selection.selected_lags))
    assert 'OLS' in result.regression_description
    assert result.backend=='arch' and result.backend_version=='8.0.0'
    assert result.versions['arch']=='8.0.0'
    assert result.bandwidth is None and result.test_type is None


@pytest.mark.parametrize("f,cls,test,trend", AUTO_CASES)
def test_automatic_lag_zero_explicit_K_and_default_maic(f,cls,test,trend):
    y=np.random.default_rng(40).normal(size=50)
    a=f(y,trend=trend,max_lags=0)
    assert a.lags==0 and a.lag_selection.method=='maic'
    assert a.lag_selection.runner_up_gap is None
    check_inference(a,cls(y,trend=trend,lags=0))
    b=f(y,trend=trend,method='MAIC',max_lags=0)
    assert a.statistic==b.statistic


@pytest.mark.parametrize("f,cls,test,trend", FIXED_CASES)
def test_shortest_fixed_samples(f,cls,test,trend):
    d={'n':0,'c':1,'ct':2,'ctt':3}[trend]
    y=np.array([1.,4.,2.,8.,3.,6.])[:d+3]
    result=f(y,trend=trend,lags=np.int64(0))
    check_inference(result,cls(y,trend=trend,lags=0))
    with pytest.raises(ValueError,match='observations'):f(y[:-1],trend=trend,lags=0)


@pytest.mark.parametrize("f,cls,test,trend", AUTO_CASES)
def test_default_caps_preserve_request_metadata(f,cls,test,trend):
    y=np.array([1.,4.,2.,8.,3.,6.,1.,7.])
    d={'n':0,'c':1,'ct':2}[trend]
    bound=min((len(y)-3)//2,(len(y)-3-d)//2 if test=='ADF' else len(y)-3-d)
    with pytest.warns(SpecificationWarning,match='capped') as caught:
        result=f(y,trend=trend)
    assert len(caught)==1
    s=result.lag_selection
    assert s.max_lags==bound and s.max_lags_requested is None
    assert s.uncapped_rule_value==math.ceil(12*(len(y)/100)**.25)
    assert s.selection_nobs==len(y)-1-bound and s.notes==result.notes
    check_inference(result,cls(y,trend=trend,lags=s.selected_lags))
    with pytest.raises(ValueError,match='maximum'):
        f(y,trend=trend,max_lags=bound+1)


@pytest.mark.parametrize("f", [adf,dfgls])
@pytest.mark.parametrize("kwargs", [dict(lags=0,method='maic'),dict(lags=1,max_lags=3),
    dict(lags=0,method='aic',max_lags=3),dict(lags=True),dict(lags=0.),dict(lags=-1),
    dict(lags=100),dict(max_lags=1.),dict(max_lags=False),dict(max_lags=-1),
    dict(method='t-stat'),dict(method=True),dict(trend='t'),dict(lags=0,trend=2)])
def test_invalid_options(f,kwargs):
    with pytest.raises((TypeError,ValueError)):f(np.random.default_rng(1).normal(size=20),**kwargs)


@pytest.mark.parametrize("f", [adf,dfgls])
@pytest.mark.parametrize("y", [[],[1,np.nan,2,3,4],[1,np.inf,2,3,4],[[1],[4],[2],[8],[3]],
    [True,False,True,False,True],[1.,True,3.,2.,4.],['1','4','2','8','3'],[1+2j,4,2,8,3],
    np.array([1,4,2,8,3],dtype=object),1.])
def test_rejected_inputs(f,y):
    with pytest.raises((TypeError,ValueError)):f(y,lags=0)


@pytest.mark.parametrize("f,trend,y",[(adf,'c',np.ones(20)),(dfgls,'c',np.ones(20)),
    (adf,'ct',np.arange(20.)),(dfgls,'ct',np.arange(20.)),
    (adf,'ctt',np.arange(20.)**2),(adf,'n',2.**np.arange(20)),
    (adf,'c',np.arange(20.)),(adf,'n',np.r_[np.ones(19),2.])])
def test_degenerate_fixed_models(f,trend,y):
    p=1 if np.array_equal(y,np.r_[np.ones(19),2.]) else 0
    with pytest.raises(EstimationError):f(y,trend=trend,lags=p)


@pytest.mark.parametrize("f", [adf,dfgls])
def test_large_p_leaves_no_degrees_of_freedom(f):
    with pytest.raises(ValueError,match='degrees'):
        f(np.random.default_rng(1).normal(size=10),lags=4)


def test_fixed_ctt_only_and_dfgls_no_n():
    y=np.random.default_rng(1).normal(size=30)
    for method in [None,'maic','mbic','aic','bic']:
        with pytest.raises(ValueError,match='fixed-lag'):adf(y,trend='ctt',method=method)
    for tr in ['n','ctt']:
        with pytest.raises(ValueError,match='trend'):dfgls(y,trend=tr,lags=0)


@pytest.mark.parametrize("f,cls,test,trend", FIXED_CASES)
def test_fixed_scaling_and_copy(f,cls,test,trend):
    y=np.random.default_rng(33).normal(size=70).cumsum()
    original=y.copy()
    r=f(y,trend=trend,lags=2)
    scaled=f(y*-3,trend=trend,lags=2)
    assert_allclose(r.statistic,scaled.statistic,rtol=0,atol=1e-10)
    assert_allclose(y,original)
    y[:]=0
    assert np.any(r.regression.model.endog!=0)


def test_dfgls_original_levels_once_and_integer_lags(monkeypatch):
    import arch.unitroot as public
    calls=[];real=public.DFGLS
    def spy(y,**kwargs):
        calls.append((np.array(y,copy=True),kwargs.copy()))
        return real(y,**kwargs)
    monkeypatch.setattr(public,'DFGLS',spy)
    y=np.random.default_rng(8).normal(size=90).cumsum()+np.arange(90)*.2
    r=dfgls(y,trend='ct',method='mbic',max_lags=7)
    assert len(calls)==1
    assert_allclose(calls[0][0],y,rtol=0,atol=0)
    assert calls[0][1]==dict(trend='ct',lags=r.lags)
    check_inference(r,real(y,trend='ct',lags=r.lags))



def test_reject_unreleased_or_different_backend_before_constructor(monkeypatch):
    monkeypatch.setattr(arch,'__version__','8.0.1.dev0')
    with pytest.raises(RuntimeError,match='released 8.0.0'):
        adf([1.,4.,2.,8.,3.],lags=0)


@pytest.mark.parametrize("f", [adf,dfgls])
def test_eager_backend_failure_and_exception_chain(monkeypatch,f):
    import arch.unitroot as public
    def fail(*args,**kwargs):raise RuntimeError('backend failed')
    monkeypatch.setattr(public,'ADF' if f is adf else 'DFGLS',fail)
    with pytest.raises(EstimationError,match='ARCH estimation failed') as caught:
        f([1.,4.,2.,8.,3.],lags=0)
    assert isinstance(caught.value.__cause__,RuntimeError)


def test_unexpected_programming_error_is_not_swallowed(monkeypatch):
    import arch.unitroot as public
    def fail(*args,**kwargs):raise KeyError('programming error')
    monkeypatch.setattr(public,'ADF',fail)
    with pytest.raises(KeyError):adf([1.,4.,2.,8.,3.],lags=0)


def test_nonfinite_inference_and_wrong_sample_rejected(monkeypatch):
    import arch.unitroot as public
    y=np.random.default_rng(7).normal(size=30)
    real=ADF(y,lags=0)
    def fake(**overrides):
        values=dict(stat=real.stat,pvalue=real.pvalue,nobs=real.nobs,lags=0,
                    critical_values=real.critical_values,regression=real.regression)
        values.update(overrides)
        return lambda *args,**kwargs:SimpleNamespace(**values)
    for override in [dict(stat=float('nan')),dict(pvalue=1.2),dict(pvalue=-.1),
                     dict(critical_values={'1%':float('inf'),'5%':-2.,'10%':-1.}),
                     dict(nobs=28),dict(lags=1)]:
        monkeypatch.setattr(public,'ADF',fake(**override))
        with pytest.raises(EstimationError):adf(y,lags=0)


def test_bad_numeric_scale_rejected():
    with pytest.raises(EstimationError,match='scale|non-finite'):
        adf(np.array([1.,4.,2.,8.,3.,6.,1.])*1e250,lags=0)


def test_poor_condition_warning_is_returned_as_a_note():
    y=1.1**np.arange(20)+np.random.default_rng(4).normal(scale=1e-9,size=20)
    with pytest.warns(SpecificationWarning,match='conditioned'):
        r=adf(y,trend='n',lags=1)
    assert any('conditioned' in note for note in r.notes)
    assert np.isfinite(r.statistic)


def test_regression_guard_rejects_invalid_backend_state(monkeypatch):
    import arch.unitroot as public
    y=np.random.default_rng(20).normal(size=20)
    real=ADF(y,lags=0)
    reg=real.regression
    valid=dict(model=SimpleNamespace(exog=reg.model.exog,endog=reg.model.endog),
               resid=reg.resid,df_resid=reg.df_resid,ssr=reg.ssr,scale=reg.scale,nobs=reg.nobs)
    bad_states=[
        dict(model=SimpleNamespace(exog=np.ones((2,2)),endog=reg.model.endog)),
        dict(df_resid=0),
        dict(resid=np.full(19,float('nan'))),
        dict(model=SimpleNamespace(exog=np.column_stack([np.zeros(19),np.ones(19)]),endog=reg.model.endog)),
        dict(model=SimpleNamespace(exog=np.ones((19,2)),endog=reg.model.endog)),
        dict(model=SimpleNamespace(exog=reg.model.exog,endog=np.zeros(19))),
        dict(resid=np.zeros(19)),
        dict(ssr=0),
    ]
    for overrides in bad_states:
        fake_reg=SimpleNamespace(**(valid|overrides))
        backend=SimpleNamespace(stat=real.stat,pvalue=real.pvalue,nobs=real.nobs,lags=0,
                                critical_values=real.critical_values,regression=fake_reg)
        monkeypatch.setattr(public,'ADF',lambda *args,_b=backend,**kwargs:_b)
        with pytest.raises(EstimationError):adf(y,lags=0)


def test_runtime_warning_is_not_hidden(monkeypatch):
    import warnings
    import arch.unitroot as public
    def fail(*args,**kwargs):warnings.warn('invalid variance',RuntimeWarning)
    monkeypatch.setattr(public,'ADF',fail)
    with pytest.raises(EstimationError) as caught:adf([1.,4.,2.,8.,3.],lags=0)
    assert isinstance(caught.value.__cause__,RuntimeWarning)


@pytest.mark.parametrize('name',SERIES)
@pytest.mark.parametrize('trend',['c','ct'])
@pytest.mark.parametrize('p',[0,4])
def test_saved_urca_statistics_do_not_replace_arch_calibration(stage2,name,trend,p):
    r=dfgls(stage2[0][name],trend=trend,lags=p)
    saved=[row for row in stage2[3] if (row['series'],row['test'],row['trend'],row['lags'],row['variant'])==
           (name,'DFGLS',trend,str(p),'fixed')][0]
    assert_allclose(r.statistic,float(saved['statistic']),rtol=0,atol=1e-8)
    assert any(abs(r.critical_values[a]-float(saved[key]))>1e-6 for a,key in [(0.01,'cv1'),(0.05,'cv5'),(0.10,'cv10')])


@pytest.mark.parametrize('f,cls,test,trend',AUTO_CASES)
def test_uncapped_default_matches_independent_default(f,cls,test,trend):
    y=np.random.default_rng(35).normal(size=200)
    r=f(y,trend=trend)
    selection=select_lags(y,trend=trend)
    assert r.lag_selection==selection
    assert r.lag_selection.max_lags_requested is None
    assert r.lag_selection.max_lags==r.lag_selection.uncapped_rule_value==15
    check_inference(r,cls(y,trend=trend,lags=selection.selected_lags))


@pytest.mark.parametrize('f,cls,test,trend',AUTO_CASES)
def test_automatic_shortest_feasible_sample(f,cls,test,trend):
    d={'n':0,'c':1,'ct':2}[trend]
    y=np.array([1.,4.,2.,8.,3.])[:d+3]
    with pytest.warns(SpecificationWarning,match='capped'):
        r=f(y,trend=trend)
    assert r.lags==0 and r.lag_selection.max_lags==0
    assert r.lag_selection.selection_nobs==r.nobs==len(y)-1
    check_inference(r,cls(y,trend=trend,lags=0))
