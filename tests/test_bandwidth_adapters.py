"""PP/KPSS public-ARCH integration and reviewed peer-package differences."""
from dataclasses import replace
import math
from types import SimpleNamespace
import warnings

import arch
from arch.unitroot import PhillipsPerron, KPSS
import numpy as np
from numpy.testing import assert_allclose
import pytest

from unitroot_tools import EstimationError, SpecificationWarning, pp, kpss

SERIES=['rw','rw_negma','ar07','trend_ar08','hetero_rw','mic_boundary']
PP_CASES=[(tr,typ) for tr in ['n','c','ct'] for typ in ['tau','rho']]


def native_kpss(y,**kwargs):
    with warnings.catch_warnings():
        # Reference-call noise only. Production suppression is tested separately.
        warnings.simplefilter('ignore',DeprecationWarning)
        obj=KPSS(y,**kwargs)
        _=obj.stat
        return obj


def saved_rows(filename):
    import csv
    from conftest import RESULTS
    with (RESULTS/filename).open() as f:return list(csv.DictReader(f))


def check_result(result,native):
    assert result.statistic==native.stat and result.pvalue==native.pvalue
    assert dict(result.critical_values)=={.01:native.critical_values['1%'],.05:native.critical_values['5%'],.10:native.critical_values['10%']}
    assert result.bandwidth==native.lags and result.nobs==native.nobs
    assert result.lags is None and result.lag_selection is None
    for alpha in [.01,.05,.10]:
        c=result.critical_values[alpha]
        assert result.reject(alpha)==bool(result.statistic>c if result.test=='KPSS' else result.statistic<c)


@pytest.mark.parametrize('name',SERIES)
@pytest.mark.parametrize('trend,typ',PP_CASES)
@pytest.mark.parametrize('b',[0,4,12])
def test_fixed_pp_references(stage2,name,trend,typ,b):
    y=stage2[0][name];r=pp(y,trend=trend,bandwidth=b,test_type=typ)
    native=PhillipsPerron(y,trend=trend,lags=b,test_type=typ)
    check_result(r,native)
    assert (r.input_nobs,r.nobs,r.sample_start,r.sample_stop)==(len(y),len(y)-1,1,len(y))
    assert r.bandwidth_method=='fixed' and r.bandwidth_requested==b and r.bandwidth_uncapped is None
    assert r.test_type==typ and r.tail=='left'
    assert r.regression.cov_type=='HAC' and r.regression.cov_kwds['maxlags']==b
    assert 'HAC/Newey-West' in r.regression_description
    assert ('ADF-t' if typ=='tau' else 'ADF-z') in r.calibration
    assert 'unit root' in r.null_hypothesis
    if trend!='n' and b!=0:
        saved=next(row for row in saved_rows('python_release.csv') if
                   (row['series'],row['backend'],row['test'],row['trend'],row['lags'])==
                   (name,'ARCH','PP_'+typ,trend,str(b)))
        assert_allclose([r.statistic,r.pvalue,*r.critical_values.values()],
                        [float(saved[k]) for k in ['statistic','pvalue','cv1','cv5','cv10']],rtol=0,atol=1e-8)


@pytest.mark.parametrize('name',SERIES)
@pytest.mark.parametrize('trend',['c','ct'])
@pytest.mark.parametrize('b',[0,4,12])
def test_fixed_kpss_references(stage2,name,trend,b):
    y=stage2[0][name];r=kpss(y,trend=trend,bandwidth=b)
    check_result(r,native_kpss(y,trend=trend,lags=b))
    assert (r.input_nobs,r.nobs,r.sample_start,r.sample_stop)==(len(y),len(y),0,len(y))
    assert r.bandwidth_method=='fixed' and r.bandwidth_requested==b and r.bandwidth_uncapped is None
    assert r.regression is None and r.test_type is None and r.tail=='right'
    assert 'stationary around' in r.null_hypothesis and 'unit root' not in r.null_hypothesis
    assert 'not stationary' in r.alternative_hypothesis
    assert 'simulated' in r.calibration and 'interpolation' in r.pvalue_provenance
    if b!=0:
        saved=next(row for row in saved_rows('python_release.csv') if
                   (row['series'],row['backend'],row['test'],row['trend'],row['lags'])==
                   (name,'ARCH','KPSS',trend,str(b)))
        assert_allclose([r.statistic,r.pvalue,*r.critical_values.values()],
                        [float(saved[k]) for k in ['statistic','pvalue','cv1','cv5','cv10']],rtol=0,atol=1e-8)


@pytest.mark.parametrize('name',SERIES)
@pytest.mark.parametrize('trend,typ',PP_CASES)
def test_auto_pp_schwert_uncapped(stage2,name,trend,typ):
    y=stage2[0][name];rule=math.ceil(12*(len(y)/100)**.25)
    r=pp(y,trend=trend,test_type=typ)
    assert r.bandwidth==rule and r.bandwidth_uncapped==rule
    assert r.bandwidth_requested is None and r.bandwidth_method=='schwert_capped'
    check_result(r,PhillipsPerron(y,trend=trend,lags=rule,test_type=typ))


@pytest.mark.parametrize('name',SERIES)
@pytest.mark.parametrize('trend',['c','ct'])
def test_auto_kpss_hobijn_public_bandwidth(stage2,name,trend):
    y=stage2[0][name]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        r=kpss(y,trend=trend)
    assert not caught
    check_result(r,native_kpss(y,trend=trend,lags=None))
    assert r.bandwidth_requested is None and r.bandwidth_uncapped is None
    assert r.bandwidth_method=='hobijn' and 0<=r.bandwidth<len(y)


@pytest.mark.parametrize('trend,typ',PP_CASES)
def test_pp_shortest_samples_and_automatic_cap(trend,typ):
    d={'n':0,'c':1,'ct':2}[trend];y=np.array([1.,4.,2.,8.,3.])[:d+3]
    for b in [0,len(y)-2]:check_result(pp(y,trend=trend,bandwidth=b,test_type=typ),PhillipsPerron(y,trend=trend,lags=b,test_type=typ))
    with pytest.warns(SpecificationWarning,match='capped') as caught:
        r=pp(y,trend=trend,test_type=typ)
    assert len(caught)==1 and r.bandwidth==len(y)-2
    assert r.bandwidth_uncapped>r.bandwidth and 'capped' in r.notes[0]
    check_result(r,PhillipsPerron(y,trend=trend,lags=len(y)-2,test_type=typ))
    with pytest.raises(ValueError,match='observations'):pp(y[:-1],trend=trend,bandwidth=0,test_type=typ)
    with pytest.raises(ValueError,match='less than'):pp(y,trend=trend,bandwidth=len(y)-1,test_type=typ)


@pytest.mark.parametrize('trend',['c','ct'])
def test_kpss_shortest_fixed_and_invalid_auto(trend):
    d={'c':1,'ct':2}[trend];y=np.array([1.,4.,2.,8.])[:d+2]
    for b in [0,len(y)-1]:check_result(kpss(y,trend=trend,bandwidth=b),native_kpss(y,trend=trend,lags=b))
    with pytest.raises(ValueError,match='observations'):kpss(y[:-1],trend=trend,bandwidth=0)
    with pytest.raises(ValueError,match='less than'):kpss(y,trend=trend,bandwidth=len(y))
    # These short samples pass length feasibility but native Hobijn fails.
    with pytest.raises(RuntimeError):native_kpss(y,trend=trend,lags=None)
    with pytest.raises(EstimationError) as caught:kpss(y,trend=trend)
    assert isinstance(caught.value.__cause__,RuntimeError)


@pytest.mark.parametrize('f,trend',[ (pp,'n'),(pp,'c'),(pp,'ct'),(kpss,'c'),(kpss,'ct') ])
def test_largest_fixed_bandwidth(f,trend):
    y=np.random.default_rng(16).normal(size=40)
    b=len(y)-2 if f is pp else len(y)-1
    r=f(y,trend=trend,bandwidth=np.int64(b))
    native=PhillipsPerron(y,trend=trend,lags=b) if f is pp else native_kpss(y,trend=trend,lags=b)
    check_result(r,native)
    with pytest.raises(ValueError):f(y,trend=trend,bandwidth=b+1)


@pytest.mark.parametrize('f',[pp,kpss])
@pytest.mark.parametrize('b',[True,np.bool_(False),0.,'auto',-1,1000])
def test_bandwidth_input_errors(f,b):
    with pytest.raises((TypeError,ValueError)):f([1.,4.,2.,8.,3.,5.],bandwidth=b)


@pytest.mark.parametrize('f',[pp,kpss])
@pytest.mark.parametrize('y',[[],1.,[[1],[4],[2],[8],[3]],[1,np.nan,2,3,4],[1,np.inf,2,3,4],
                             [True,False,True,False,True],[1.,True,2.,8.,3.],['1','4','2','8','3'],
                             np.array([1,4,2,8,3],dtype=object),[1+1j,4,2,8,3]])
def test_shared_malformed_inputs(f,y):
    with pytest.raises((TypeError,ValueError)):f(y,bandwidth=0)


@pytest.mark.parametrize('f',[pp,kpss])
@pytest.mark.parametrize('kwargs',[{'lags':0},{'method':'maic'},{'max_lags':4},{'trend':'ctt'},{'trend':None}])
def test_options_cannot_conflate_bandwidth_and_augmentation(f,kwargs):
    with pytest.raises((TypeError,ValueError)):f([1.,4.,2.,8.,3.,5.],**kwargs)


def test_test_type_and_trend_labels():
    y=np.random.default_rng(7).normal(size=30)
    assert pp(y,trend='CT',test_type='RHO',bandwidth=0).test_type=='rho'
    assert kpss(y,trend='C',bandwidth=0).trend=='c'
    for v in ['t',None,1]:
        with pytest.raises((TypeError,ValueError)):pp(y,test_type=v,bandwidth=0)
    with pytest.raises(ValueError):kpss(y,trend='n',bandwidth=0)


@pytest.mark.parametrize('f,trend,y',[(pp,'c',np.ones(20)),(pp,'n',np.ones(20)),
    (pp,'ct',np.arange(20.)),(kpss,'c',np.ones(20)),(kpss,'ct',np.arange(20.)),
    (pp,'n',2.**np.arange(20))])
def test_degenerate_and_zero_variance_models(f,trend,y):
    with pytest.raises(EstimationError):f(y,trend=trend,bandwidth=0)


@pytest.mark.parametrize('f',[pp,kpss])
def test_bad_numeric_scale(f):
    with pytest.raises(EstimationError):f(np.array([1.,4.,2.,8.,3.,5.])*1e250,bandwidth=0)


@pytest.mark.parametrize('trend',['c','ct'])
@pytest.mark.parametrize('typ',['tau','rho'])
@pytest.mark.parametrize('b',[4,12])
@pytest.mark.parametrize('name',SERIES)
def test_saved_pp_finite_sample_difference_from_urca(stage2,name,trend,typ,b):
    r=pp(stage2[0][name],trend=trend,bandwidth=b,test_type=typ)
    peer=next(row for row in stage2[3] if (row['series'],row['test'],row['trend'],row['lags'],row['variant'])==
              (name,'PP_'+typ,trend,str(b),'fixed'))
    own=next(row for row in saved_rows('python_release.csv') if (row['series'],row['backend'],row['test'],row['trend'],row['lags'])==
             (name,'ARCH','PP_'+typ,trend,str(b)))
    delta=r.statistic-float(peer['statistic'])
    assert abs(delta)>1e-8
    assert_allclose(delta,float(own['statistic'])-float(peer['statistic']),rtol=0,atol=1e-8)


@pytest.mark.parametrize('name',SERIES)
@pytest.mark.parametrize('trend',['c','ct'])
@pytest.mark.parametrize('b',[4,12])
def test_kpss_statsmodels_calibration_difference(stage2,name,trend,b):
    r=kpss(stage2[0][name],trend=trend,bandwidth=b)
    peer=next(row for row in saved_rows('python_release.csv') if (row['series'],row['backend'],row['test'],row['trend'],row['lags'])==
              (name,'Statsmodels','KPSS',trend,str(b)))
    assert_allclose(r.statistic,float(peer['statistic']),rtol=0,atol=1e-8)
    assert r.pvalue!=float(peer['pvalue'])
    assert all(r.critical_values[a]!=float(peer[k]) for a,k in [(.01,'cv1'),(.05,'cv5'),(.10,'cv10')])


@pytest.mark.parametrize('b',[0,4,12])
@pytest.mark.parametrize('typ',['tau','rho'])
def test_pp_no_deterministics_independent_formula(b,typ):
    y=np.random.default_rng(31).normal(size=80).cumsum()
    x,z=y[:-1],y[1:];n=len(x);rho=x@z/(x@x);u=z-rho*x
    s2=u@u/(n-1);gamma=u@u/n;se=np.sqrt(s2/(x@x))
    lam=gamma+sum(2*(1-j/(b+1))*(u[j:]@u[:-j])/n for j in range(1,b+1))
    value=np.sqrt(gamma/lam)*(rho-1)/se-.5*(lam-gamma)/np.sqrt(lam)*(n*se/np.sqrt(s2)) if typ=='tau' else n*(rho-1)-.5*n*n*se**2/s2*(lam-gamma)
    assert_allclose(pp(y,trend='n',bandwidth=b,test_type=typ).statistic,value,rtol=0,atol=1e-8)


def test_no_augmentation_selector_is_called(monkeypatch):
    import unitroot_tools.selection as selection
    import unitroot_tools._adapters as adapters
    def forbidden(*args,**kwargs):raise AssertionError('augmentation selector called')
    monkeypatch.setattr(selection,'select_lags',forbidden)
    monkeypatch.setattr(adapters,'select_lags',forbidden)
    y=np.random.default_rng(12).normal(size=50)
    assert pp(y,bandwidth=0).lags is None
    assert kpss(y,bandwidth=0).lags is None


def test_public_bandwidth_handoff_and_no_kpss_regression(monkeypatch):
    import arch.unitroot as public
    y=np.random.default_rng(13).normal(size=60)
    real_pp,real_kpss=public.PhillipsPerron,public.KPSS
    calls=[]
    def pp_spy(values,**kw):
        calls.append(('PP',np.array(values,copy=True),kw.copy()))
        return real_pp(values,**kw)
    def kpss_spy(values,**kw):
        calls.append(('KPSS',np.array(values,copy=True),kw.copy()))
        return real_kpss(values,**kw)
    monkeypatch.setattr(public,'PhillipsPerron',pp_spy)
    monkeypatch.setattr(public,'KPSS',kpss_spy)
    a=pp(y,bandwidth=0,test_type='rho')
    b=kpss(y)
    assert calls[0][2]=={'trend':'c','lags':0,'test_type':'rho'}
    assert calls[1][2]=={'trend':'c','lags':None}
    for _,values,_ in calls:assert_allclose(values,y,rtol=0,atol=0)
    assert a.regression is not None and b.regression is None
    native=native_kpss(y,lags=4)
    class NoRegression:
        stat=native.stat;pvalue=native.pvalue;nobs=native.nobs;lags=native.lags;critical_values=native.critical_values
        @property
        def regression(self):raise AssertionError('KPSS regression accessed')
    monkeypatch.setattr(public,'KPSS',lambda *a,**kw:NoRegression())
    assert kpss(y,bandwidth=4).regression is None


def test_only_exact_legacy_kpss_warning_is_suppressed(monkeypatch):
    import arch.unitroot as public
    real=public.KPSS
    def warn(values,**kwargs):
        warnings.warn('Different KPSS deprecation',DeprecationWarning)
        warnings.warn('Unrelated user warning',UserWarning)
        return real(values,**kwargs)
    monkeypatch.setattr(public,'KPSS',warn)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        kpss(np.random.default_rng(10).normal(size=80))
    assert [str(w.message) for w in caught]==['Different KPSS deprecation','Unrelated user warning']


@pytest.mark.parametrize('f,public_name',[(pp,'PhillipsPerron'),(kpss,'KPSS')])
def test_numerical_backend_errors_and_programming_errors(monkeypatch,f,public_name):
    import arch.unitroot as public
    y=np.random.default_rng(12).normal(size=50)
    def fail(*args,**kwargs):raise RuntimeError('numerical failure')
    monkeypatch.setattr(public,public_name,fail)
    with pytest.raises(EstimationError) as caught:f(y,bandwidth=0)
    assert isinstance(caught.value.__cause__,RuntimeError)
    def unexpected(*args,**kwargs):raise KeyError('programming error')
    monkeypatch.setattr(public,public_name,unexpected)
    with pytest.raises(KeyError):f(y,bandwidth=0)


def test_numeric_warnings_are_chained(monkeypatch):
    import arch.unitroot as public
    def bad(*args,**kwargs):warnings.warn('invalid covariance',RuntimeWarning)
    monkeypatch.setattr(public,'KPSS',bad)
    with pytest.raises(EstimationError) as caught:kpss(np.arange(20.)**2,bandwidth=0)
    assert isinstance(caught.value.__cause__,RuntimeWarning)


def test_kpss_invalid_public_bandwidth_is_rejected_not_capped(monkeypatch):
    import arch.unitroot as public
    y=np.random.default_rng(14).normal(size=30);native=native_kpss(y,lags=0)
    for b in [len(y),-1,True,1.5]:
        obj=SimpleNamespace(stat=native.stat,pvalue=native.pvalue,nobs=len(y),lags=b,critical_values=native.critical_values)
        monkeypatch.setattr(public,'KPSS',lambda *args,_b=obj,**kw:_b)
        with pytest.raises(EstimationError,match='specify a fixed bandwidth'):kpss(y)


@pytest.mark.parametrize('test',['PP','KPSS'])
def test_invalid_backend_inference_or_sample(monkeypatch,test):
    import arch.unitroot as public
    y=np.random.default_rng(15).normal(size=30)
    f=pp if test=='PP' else kpss
    native=PhillipsPerron(y,lags=4) if test=='PP' else native_kpss(y,lags=4)
    values=dict(stat=native.stat,pvalue=native.pvalue,nobs=native.nobs,lags=4,critical_values=native.critical_values)
    if test=='PP':values|=dict(regression=native.regression,test_type='tau')
    changes=[{'stat':float('nan')},{'pvalue':1.2},{'nobs':1},{'lags':3},
             {'critical_values':{'1%':float('inf'),'5%':-2.,'10%':-1.}}]
    if test=='KPSS':changes += [{'stat':0.},{'stat':-1.},{'pvalue':0.}]
    else:changes += [{'test_type':'rho'}]
    for change in changes:
        fake=SimpleNamespace(**(values|change))
        monkeypatch.setattr(public,'PhillipsPerron' if test=='PP' else 'KPSS',lambda *a,_b=fake,**kw:_b)
        with pytest.raises(EstimationError):f(y,bandwidth=4)


@pytest.mark.parametrize('f',[pp,kpss])
def test_released_backend_version_required(monkeypatch,f):
    monkeypatch.setattr(arch,'__version__','8.1.0')
    with pytest.raises(RuntimeError,match='released 8.0.0'):f([1.,4.,2.,8.,3.,5.],bandwidth=0)


@pytest.mark.parametrize('alpha',[.01,.05,.10])
def test_kpss_right_tail_boundary_and_metadata(alpha):
    r=kpss(np.random.default_rng(20).normal(size=80),bandwidth=4)
    cv=r.critical_values[alpha]
    assert replace(r,statistic=cv,pvalue=.0001).reject(alpha) is False
    assert replace(r,statistic=np.nextafter(cv,np.inf),pvalue=1.).reject(alpha) is True
    assert replace(r,statistic=np.nextafter(cv,-np.inf),pvalue=.0001).reject(alpha) is False
    with pytest.raises(ValueError,match='right-tail'):replace(r,tail='left')
    text=r.summary()
    assert 'KPSS stationarity test' in text and 'Long-run-variance bandwidth b: 4 (fixed)' in text
    assert 'Augmentation order' not in text and 'Lag selection:' not in text
    assert 'does not establish the null' in text and 'accept H0' not in text
    assert r.lags is None and r.lag_selection is None


@pytest.mark.parametrize('typ',['tau','rho'])
def test_pp_left_tail_boundary_and_summary(typ):
    r=pp(np.random.default_rng(21).normal(size=80),bandwidth=0,test_type=typ)
    for a,c in r.critical_values.items():
        assert replace(r,statistic=c).reject(a) is False
        assert replace(r,statistic=np.nextafter(c,-np.inf)).reject(a) is True
    text=r.summary()
    assert f'Statistic type: {typ}' in text and 'bandwidth b: 0 (fixed)' in text
    assert 'HAC/Newey-West' in text and 'Augmentation order' not in text


@pytest.mark.parametrize('test',['PP','KPSS'])
@pytest.mark.parametrize('endpoint',['lower','upper'])
def test_endpoint_approximation_notes(monkeypatch,test,endpoint):
    import arch.unitroot as public
    y=np.random.default_rng(22).normal(size=50)
    native=PhillipsPerron(y,lags=0) if test=='PP' else native_kpss(y,lags=0)
    p=1. if endpoint=='upper' else 0. if test=='PP' else .0001
    values=dict(stat=native.stat,pvalue=p,nobs=native.nobs,lags=0,critical_values=native.critical_values)
    if test=='PP':values|=dict(regression=native.regression,test_type='tau')
    obj=SimpleNamespace(**values)
    monkeypatch.setattr(public,'PhillipsPerron' if test=='PP' else 'KPSS',lambda *a,**kw:obj)
    r=(pp if test=='PP' else kpss)(y,bandwidth=0)
    assert r.pvalue_boundary==endpoint+'_endpoint'
    assert ('simulated-table' if test=='KPSS' else 'response-surface') in r.notes[0]
    assert 'not an exact tail probability' in r.notes[0]
    assert r.pvalue==p


def test_input_name_and_copy_and_immutable_bandwidth():
    import pandas as pd
    from dataclasses import FrozenInstanceError
    y=pd.Series(np.random.default_rng(23).normal(size=60),name='inflation')
    a=pp(y,bandwidth=0);b=kpss(y,bandwidth=0)
    assert a.input_name==b.input_name=='inflation'
    assert 'Series: inflation' in a.summary() and 'Series: inflation' in b.summary()
    original=a.regression.model.endog.copy();y[:]=0
    assert_allclose(a.regression.model.endog,original)
    with pytest.raises(FrozenInstanceError):b.bandwidth=2


@pytest.mark.parametrize('name',SERIES)
@pytest.mark.parametrize('trend',['c','ct'])
def test_fresh_statsmodels_kpss_statistic_and_distinct_critical_values(stage2,name,trend):
    from statsmodels.tsa.stattools import kpss as sm_kpss
    from statsmodels.tools.sm_exceptions import InterpolationWarning
    y=stage2[0][name]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',InterpolationWarning)
        warnings.filterwarnings('ignore',message='^kpss currently returns a plain tuple.*',category=FutureWarning)
        statistic,p,b,critical=sm_kpss(y,regression=trend,nlags=0)
    r=kpss(y,trend=trend,bandwidth=0)
    assert_allclose(r.statistic,statistic,rtol=0,atol=1e-8)
    assert any(r.critical_values[a]!=critical[k] for a,k in [(.01,'1%'),(.05,'5%'),(.10,'10%')])
