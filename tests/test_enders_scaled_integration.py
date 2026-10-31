"""Independent Enders regression accuracy and unchanged sequential contracts."""
from dataclasses import replace
import importlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from statsmodels.regression.linear_model import OLSResults, RegressionResultsWrapper

from unitroot_tools import adf,enders,EstimationError,SpecificationWarning
from unitroot_tools._enders import _joint,_auxiliary,_enders_regression

module=importlib.import_module('unitroot_tools._enders')
ROOT=Path(__file__).resolve().parents[1]
CASES=json.loads((ROOT/'validation/enders/scaled_references.json').read_text())['cases']
INPUTS=json.loads((ROOT/'validation/scaled_svd/adf_inputs.json').read_text())['inputs']


def check_fit(fit,reference):
    assert isinstance(fit,RegressionResultsWrapper) and isinstance(fit._results,OLSResults)
    cov=np.array(reference['covariance']);se=np.sqrt(cov.diagonal())
    assert list(fit.params.index)==reference['names']
    assert np.max(abs(np.asarray(fit.params)-reference['params'])/se)<1e-12
    assert np.max(abs(np.asarray(fit.cov_params())-cov)/np.sqrt(np.outer(cov.diagonal(),cov.diagonal())))<1e-12
    assert_allclose(fit.ssr,float(reference['rss']),rtol=1e-12,atol=0)
    assert_allclose(fit.tvalues,reference['tvalues'],rtol=0,atol=1e-12)
    assert fit.df_resid==reference['df_resid']
    assert_array_equal(fit.resid,fit.model.endog-fit.predict())


@pytest.mark.parametrize('case',CASES,ids=lambda c:f"{c['key']}-{c['trend']}-p{c['lags']}")
def test_independent_exact_designs_and_high_precision(case,monkeypatch):
    y=np.array([float.fromhex(v) for v in INPUTS[case['key']]])
    trend,p=case['trend'],case['lags'];delta=np.diff(y);response=delta[p:]
    columns=[y[p:-1]]+[delta[p-j:len(delta)-j] for j in range(1,p+1)]+[np.ones(len(response))]
    names=['Level.L1']+[f'Diff.L{j}' for j in range(1,p+1)]+['const']
    if trend=='ct':columns.append(np.arange(1,len(response)+1,dtype=float));names.append('trend')
    x=np.column_stack(columns);r=adf(y,trend=trend,lags=p)
    assert_array_equal(r.regression.model.exog,x);assert_array_equal(r.regression.model.endog,response)
    check_fit(r.regression,case['fits']['unrestricted'])
    calls=[];real=module._enders_regression
    def spy(y,design,nobs,context):
        fit=real(y,design,nobs,context);calls.append((design,fit,context));return fit
    monkeypatch.setattr(module,'_enders_regression',spy)
    phi,rss_r,kept=_joint(r)
    removed=('Level.L1','trend' if trend=='ct' else 'const')
    indices=[j for j,name in enumerate(names) if name not in removed]
    assert kept==tuple(names[j] for j in indices)
    assert_allclose(rss_r,float(case['fits']['restricted']['rss']),rtol=1e-12,atol=0)
    assert_allclose(phi,case['phi'],rtol=0,atol=1e-12)
    assert phi==max(0.,(rss_r-r.regression.ssr)/2/(r.regression.ssr/r.regression.df_resid))
    if indices:
        assert_array_equal(calls[0][0],x[:,indices]);assert_array_equal(calls[0][1].model.endog,response)
        check_fit(calls[0][1],case['fits']['restricted'])
    else:assert calls==[]
    aux,stat=_auxiliary(r)
    assert_array_equal(aux.model.exog,x[:,1:]);assert_array_equal(aux.model.endog,response)
    check_fit(aux,case['fits']['auxiliary'])
    assert stat==aux.tvalues.iloc[-1]


MODES=['reject','drop','aux_drop','conditional_reject','conditional_nonreject']
@pytest.mark.parametrize('a,b,n_reject,flag',list(itertools.product(MODES,MODES,[False,True],[False,True])))
def test_all_conservative_paths_and_policy(a,b,n_reject,flag,monkeypatch):
    y=np.array([float.fromhex(v) for v in INPUTS['y']])
    actual={tr:adf(y,trend=tr,lags=1) for tr in ('ct','c','n')}
    modes={'ct':a,'c':b,'n':'reject' if n_reject else 'drop'};calls=[]
    def fake_adf(y,*,trend,lags,**kw):
        calls.append('fit_'+trend);mode=modes[trend]
        return replace(actual[trend],statistic=-20. if mode=='reject' else -2. if mode=='conditional_reject' else 0.)
    def fake_joint(r):
        calls.append('joint_'+r.trend);return (0. if modes[r.trend]=='drop' else 100.),1.,()
    def fake_aux(r):
        calls.append('aux_'+r.trend);fit,_=_auxiliary(actual[r.trend]);return fit,0. if modes[r.trend]=='aux_drop' else 100.
    monkeypatch.setattr(module,'adf',fake_adf);monkeypatch.setattr(module,'_joint',fake_joint);monkeypatch.setattr(module,'_auxiliary',fake_aux)
    r=enders(y,lags=1,variant='conservative',allow_no_constant=flag)
    expected_calls=[];expected_tests=[]
    for tr,tau,phi,term in [('ct','tau3','phi3','trend'),('c','tau2','phi1','intercept'),('n','tau1',None,None)]:
        mode=modes[tr];expected_calls.append('fit_'+tr);expected_tests.append(tau)
        if mode=='reject' or tr=='n' or (tr=='c' and not flag):terminal=tr;rejected=mode=='reject';break
        expected_calls.append('joint_'+tr);expected_tests.append(phi)
        if mode=='drop':continue
        expected_calls.append('aux_'+tr);expected_tests.append('auxiliary_'+term)
        if mode=='aux_drop':continue
        expected_tests.append('conditional_t');terminal=tr;rejected=mode=='conditional_reject';break
    assert calls==expected_calls and [h.test for h in r.history]==expected_tests
    assert r.trend==terminal and r.reject_unit_root==rejected
    assert r.terminal_test==expected_tests[-1] and r.lag_method=='fixed' and r.lag_selection is None
    assert r.allow_no_constant==flag and r.allow_no_constant_requested==flag
    assert (r.nobs,r.sample_start,r.sample_stop)==(len(y)-2,2,len(y))


@pytest.mark.parametrize('kind',['zero_column','duplicate','zero_response','perfect','nonfinite','no_df'])
def test_enders_fit_fails_closed(kind):
    t=np.arange(30.);y=np.sin(.7*t);x=np.c_[np.ones(30),t]
    if kind=='zero_column':x[:,1]=0
    if kind=='duplicate':x[:,1]=1
    if kind=='zero_response':y[:]=0
    if kind=='perfect':y=x[:,0]+x[:,1]
    if kind=='nonfinite':x[3,1]=np.inf
    if kind=='no_df':x=np.eye(30)
    with pytest.raises(EstimationError):_enders_regression(y,pd.DataFrame(x),len(y),'Enders test regression')


def test_full_rank_near_singular_warning_not_reduced_rank():
    t=np.linspace(-1,1,118);x=pd.DataFrame({'a':np.ones(118),'b':1+1e-9*t})
    with pytest.warns(SpecificationWarning,match='conditioned'):
        r=_enders_regression(np.sin(3*t)+.1*np.cos(7*t),x,len(t),'Enders test regression')
    assert r.df_resid==116 and np.isfinite(r.tvalues).all()
    x['b']=1+1e-14*t
    with pytest.raises(EstimationError,match='rank-deficient'):_enders_regression(np.sin(3*t),x,len(t),'Enders test regression')


@pytest.mark.parametrize('scale',[1e-250,1e250])
@pytest.mark.parametrize('variant',['original','conservative'])
def test_public_extreme_squared_unit_failure(scale,variant):
    y=np.array([float.fromhex(v) for v in INPUTS['y']])*scale
    with pytest.raises(EstimationError,match='range|scale|variance'):enders(y,lags=1,variant=variant)


def test_no_raw_ols_fallback(monkeypatch):
    from statsmodels.regression.linear_model import OLS
    def fail(*a,**kw):raise AssertionError('raw OLS.fit must not be used')
    monkeypatch.setattr(OLS,'fit',fail)
    y=np.array([float.fromhex(v) for v in INPUTS['y']]);r=adf(y,trend='ct',lags=1)
    assert np.isfinite(_joint(r)[0]) and np.isfinite(_auxiliary(r)[1])


def test_enders_fit_error_context_and_chain(monkeypatch):
    def fail(*a,**kw):raise RuntimeError('numeric failure')
    monkeypatch.setattr(module,'scaled_ols',fail)
    with pytest.raises(EstimationError,match='Enders test regression') as error:
        _enders_regression(np.sin(np.arange(30.)),pd.DataFrame({'const':np.ones(30)}),30,'Enders test regression')
    assert isinstance(error.value.__cause__,RuntimeError)


@pytest.mark.parametrize('alpha',[.01,.05,.10])
@pytest.mark.parametrize('variant',['original','conservative'])
@pytest.mark.parametrize('boundary',['tau','phi','auxiliary','conditional'])
def test_exact_significance_boundaries(alpha,variant,boundary,monkeypatch):
    from scipy.stats import t
    y=np.array([float.fromhex(v) for v in INPUTS['y']]);actual={tr:adf(y,trend=tr,lags=1) for tr in ('ct','c','n')}
    def fake_adf(y,*,trend,**kw):
        r=actual[trend];value=0.
        if trend=='ct':
            if boundary=='tau':value=r.critical_values[alpha]
            elif boundary=='conditional':value=float(t.ppf(alpha,r.regression.df_resid))
        return replace(r,statistic=value)
    def fake_joint(r):
        cv,_=module._phi_critical('phi3' if r.trend=='ct' else 'phi1',r.nobs,alpha)
        return (cv if boundary=='phi' else cv+100.),1.,()
    def fake_aux(r):
        fit,_=_auxiliary(actual[r.trend]);cv=float(t.ppf(1-alpha/2,fit.df_resid))
        return fit,cv if boundary=='auxiliary' else cv+100.
    monkeypatch.setattr(module,'adf',fake_adf);monkeypatch.setattr(module,'_joint',fake_joint);monkeypatch.setattr(module,'_auxiliary',fake_aux)
    result=enders(y,lags=1,alpha=alpha,variant=variant)
    name={'tau':'tau3','phi':'phi3','auxiliary':'auxiliary_trend','conditional':'conditional_t'}[boundary]
    step=next(h for h in result.history if h.test==name)
    assert step.statistic==step.critical_value and not step.reject
