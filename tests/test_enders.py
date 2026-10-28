"""Numerical references, published GDP replication, and controlled branch tests."""
import csv
from dataclasses import replace
import hashlib
import importlib
import itertools
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from scipy.stats import t

from unitroot_tools import adf, enders, EstimationError
from unitroot_tools._enders import _joint, _auxiliary, _phi_critical

module = importlib.import_module('unitroot_tools._enders')
REFERENCE = Path(__file__).resolve().parents[1]/'validation/enders'


def qr_fit(x, y):
    q, r = np.linalg.qr(x, mode='reduced')
    coef = np.linalg.solve(r, q.T@y)
    residual = y-x@coef
    rss = residual@residual
    inv = np.linalg.solve(r, np.eye(len(r)))
    cov = rss/(len(y)-x.shape[1])*(inv@inv.T)
    return coef, np.sqrt(np.diag(cov)), rss, cov


@pytest.mark.parametrize('trend,p', list(itertools.product(['ct', 'c'], [0, 1, 4])))
def test_independent_regressions_and_joint_restrictions(trend, p):
    y = np.random.default_rng(123).normal(size=160).cumsum()+.03*np.arange(160)
    d = np.diff(y); response=d[p:]
    cols=[y[p:-1]]+[d[p-j:-j] for j in range(1,p+1)]+[np.ones(len(response))]
    if trend=='ct': cols.append(np.arange(1,len(response)+1))
    x=np.column_stack(cols)
    coef,se,rss,cov=qr_fit(x,response)
    result=adf(y,trend=trend,lags=p); reg=result.regression
    np.testing.assert_allclose(reg.params,coef,rtol=1e-10,atol=1e-12)
    np.testing.assert_allclose(reg.bse,se,rtol=1e-10,atol=1e-12)
    assert reg.ssr==pytest.approx(rss,rel=1e-12)
    np.testing.assert_allclose(reg.tvalues,coef/se,rtol=1e-10,atol=1e-12)
    indices=[0,x.shape[1]-1] # level and trend/intercept
    restriction=np.eye(x.shape[1])[indices]
    wald=float((restriction@coef)@np.linalg.solve(restriction@cov@restriction.T,restriction@coef)/2)
    value,restricted_rss,names=_joint(result)
    assert value==pytest.approx(wald,rel=1e-10,abs=1e-12)
    keep=[j for j in range(x.shape[1]) if j not in indices]
    rr=qr_fit(x[:,keep],response)[2] if keep else response@response
    assert restricted_rss==pytest.approx(rr,rel=1e-12)
    assert value==pytest.approx((rr-rss)/2/(rss/(len(response)-x.shape[1])),rel=1e-10,abs=1e-12)
    aux,stat=_auxiliary(result); ac,ase,arss,_=qr_fit(x[:,1:],response)
    np.testing.assert_allclose(aux.params,ac,rtol=1e-10,atol=1e-12)
    np.testing.assert_allclose(aux.bse,ase,rtol=1e-10,atol=1e-12)
    assert aux.ssr==pytest.approx(arss,rel=1e-12)
    assert stat==pytest.approx(ac[-1]/ase[-1],rel=1e-10,abs=1e-12)
    np.testing.assert_array_equal(aux.model.endog,reg.model.endog)


def test_published_gdp_and_conservative_difference():
    meta=json.loads((REFERENCE/'source.json').read_text())
    assert hashlib.sha256((REFERENCE/'rgdp.csv').read_bytes()).hexdigest()==meta['csv_sha256']
    with (REFERENCE/'rgdp.csv').open() as f:y=np.log([float(row['rgdp']) for row in csv.DictReader(f)])
    r=enders(y,lags=1)
    assert [h.test for h in r.history]==['tau3','phi3','tau2','phi1','auxiliary_intercept','conditional_t']
    # Published values rounded to two decimals; 0.005 accommodates rounding only.
    for step,published in zip(r.history,[-1.49,2.97,-2.05,25.49,6.80,-2.05]):
        assert abs(step.statistic-published)<=.005
    assert r.history[1].critical_value==6.34
    assert r.history[3].critical_value==4.63
    assert r.history[1].table_sample_size==250
    assert r.trend=='c' and r.reject_unit_root and r.nobs==262
    # ARCH's trend starts at 1 on the estimation sample. Original time is t=3,...,264.
    reg=r.history[0].regression
    adjusted_intercept=reg.params['const']-2*reg.params['trend']
    assert abs(adjusted_intercept-.1248)<.00005
    for name,expected in [('Level.L1',-.0156),('Diff.L1',.3663),('trend',.0001)]:
        assert abs(reg.params[name]-expected)<.00005
    constant=r.history[2].regression
    for name,expected in [('Level.L1',-.0019),('Diff.L1',.3539),('const',.0215)]:
        assert abs(constant.params[name]-expected)<.00005
    aux=r.history[4].regression
    np.testing.assert_allclose(aux.params,[.3706,.0049],atol=.00005,rtol=0)
    conservative=enders(y,lags=1,variant='conservative')
    assert not conservative.reject_unit_root and conservative.terminal_test=='tau2'
    assert 'policy' in conservative.history[-1].next_action
    permitted=enders(y,lags=1,variant='conservative',allow_no_constant=True)
    assert [h.test for h in permitted.history]==[h.test for h in r.history]
    assert permitted.reject_unit_root
    assert 'proposed modification' in conservative.summary()


# Cases are controlled inference inputs, not size/power experiments. Numerical
# correctness of each regression is tested independently above.
@pytest.mark.parametrize('a,b,n_reject',list(itertools.product(
    ['reject','drop','aux_drop','conditional_reject','conditional_nonreject'],
    ['reject','drop','aux_drop','conditional_reject','conditional_nonreject'],[False,True])))
def test_all_original_paths_and_immediate_stopping(a,b,n_reject):
    y=np.random.default_rng(77).normal(size=90).cumsum()
    actual={tr:adf(y,trend=tr,lags=1) for tr in ['ct','c','n']}
    modes={'ct':a,'c':b,'n':'reject' if n_reject else 'drop'}
    calls=[]
    def fake_adf(y,*,trend,lags=None,**kw):
        calls.append('fit_'+trend)
        mode=modes[trend]
        value=-20. if mode=='reject' else -2. if mode=='conditional_reject' else 0.
        return replace(actual[trend],statistic=value)
    def fake_joint(r):
        calls.append('joint_'+r.trend)
        return (0. if modes[r.trend]=='drop' else 100.),1.,()
    def fake_aux(r):
        calls.append('aux_'+r.trend)
        fit,_=_auxiliary(actual[r.trend])
        return fit,0. if modes[r.trend]=='aux_drop' else 100.
    with patch.object(module,'adf',side_effect=fake_adf),patch.object(module,'_joint',side_effect=fake_joint),patch.object(module,'_auxiliary',side_effect=fake_aux):
        r=enders(y,lags=1)
    expected_calls=[];expected_tests=[]
    for tr,mode,tau,phi,term in [('ct',a,'tau3','phi3','trend'),('c',b,'tau2','phi1','intercept'),('n',modes['n'],'tau1',None,None)]:
        expected_calls.append('fit_'+tr);expected_tests.append(tau)
        if mode=='reject' or tr=='n':
            terminal=tr;reject=mode=='reject';break
        expected_calls.append('joint_'+tr);expected_tests.append(phi)
        if mode=='drop':continue
        expected_calls.append('aux_'+tr);expected_tests.append('auxiliary_'+term)
        if mode=='aux_drop':continue
        expected_tests.append('conditional_t');terminal=tr;reject=mode=='conditional_reject';break
    assert calls==expected_calls
    assert [h.test for h in r.history]==expected_tests
    assert r.trend==terminal and r.reject_unit_root==reject
    assert r.terminal_test==expected_tests[-1]
    assert all(h.null_hypothesis and h.inference and h.next_action for h in r.history)
    for h in r.history:
        assert h.regression.nobs==r.nobs
        if h.test.startswith('auxiliary_'):
            assert h.critical_value==t.ppf(.975,h.regression.df_resid)
        if h.test=='conditional_t':
            assert h.critical_value==t.ppf(.05,h.regression.df_resid)


@pytest.mark.parametrize('flag',[False,True])
def test_original_flag_irrelevant(flag):
    y=np.random.default_rng(5).normal(size=200).cumsum()
    r=enders(y,lags=1,allow_no_constant=flag)
    assert r.trend=='n' and r.allow_no_constant
    assert r.allow_no_constant_requested==flag


@pytest.mark.parametrize('method',['maic','mbic','aic','bic'])
def test_select_once_and_preserve_sample(method):
    y=np.random.default_rng(5).normal(size=200).cumsum()
    with patch.object(module,'adf',wraps=adf) as fit:
        r=enders(y,method=method,max_lags=8)
    initial=adf(y,trend='ct',method=method,max_lags=8)
    assert r.lags==initial.lags and r.lag_method==method
    assert r.lag_selection==initial.lag_selection
    assert fit.call_args_list[0].kwargs==dict(trend='ct',max_lags=8,method=method)
    for call in fit.call_args_list[1:]:
        assert call.kwargs==dict(trend=call.kwargs['trend'],lags=r.lags)
    for h in r.history:
        assert h.regression.nobs==r.nobs
        np.testing.assert_array_equal(h.regression.model.endog,np.diff(y)[r.lags:])


@pytest.mark.parametrize('n,row',[(25,25),(49,25),(50,50),(99,50),(100,100),(249,100),(250,250),(499,250),(500,500),(10000,500)])
@pytest.mark.parametrize('alpha',[.01,.05,.10])
def test_phi_table_boundaries(n,row,alpha):
    # Independently transcribed primary Table B, 10%, 5%, 1% order.
    phi1={25:[4.12,5.18,7.88],50:[3.94,4.86,7.06],100:[3.86,4.71,6.70],250:[3.81,4.63,6.52],500:[3.79,4.61,6.47]}
    phi3={25:[5.91,7.24,10.61],50:[5.61,6.73,9.31],100:[5.47,6.49,8.73],250:[5.39,6.34,8.43],500:[5.36,6.30,8.34]}
    for name,table in [('phi1',phi1),('phi3',phi3)]:
        assert _phi_critical(name,n,alpha)==(table[row][[.10,.05,.01].index(alpha)],row)


@pytest.mark.parametrize('kw,exc',[
    ({'variant':'unknown'},ValueError),({'variant':1},TypeError),
    ({'alpha':True},TypeError),({'alpha':np.bool_(True)},TypeError),({'alpha':'0.05'},TypeError),
    ({'alpha':.025},ValueError),({'alpha':np.nan},ValueError),
    ({'allow_no_constant':1},TypeError),({'lags':True},TypeError),({'lags':-1},ValueError),
    ({'lags':1,'max_lags':4},ValueError),({'lags':1,'method':'bic'},ValueError),
    ({'method':'unknown'},ValueError),({'max_lags':-1},ValueError),({'lags':99},ValueError),
])
def test_invalid_options(kw,exc):
    with pytest.raises(exc):enders(np.random.default_rng(1).normal(size=80),**kw)


@pytest.mark.parametrize('y,exc',[(np.arange(25.),ValueError),(np.ones(40),EstimationError),
    (np.arange(40.),EstimationError),([1.,np.nan]*20,ValueError),
    (np.ones((40,1)),TypeError),([True]*40,TypeError)])
def test_invalid_series(y,exc):
    with pytest.raises(exc):enders(y,lags=0)


@pytest.mark.parametrize('alpha',[.01,.05,.1])
def test_level_and_tail_metadata(alpha):
    y=np.random.default_rng(5).normal(size=200).cumsum()
    r=enders(y,lags=1,alpha=alpha)
    for step in r.history:
        if step.tail=='left':expected=step.statistic<step.critical_value
        elif step.tail=='right':expected=step.statistic>step.critical_value
        else:expected=abs(step.statistic)>step.critical_value
        assert step.reject==expected
    assert 'Nonrejection does not establish' in r.summary()


def test_conservative_never_evaluates_phi1_by_default():
    y=np.random.default_rng(5).normal(size=200).cumsum()
    with patch.object(module,'_joint',wraps=_joint) as joint:
        result=enders(y,lags=1,variant='conservative')
    assert result.terminal_test=='tau2' and not result.reject_unit_root
    assert [call.args[0].trend for call in joint.call_args_list]==['ct']


def test_selected_lag_with_too_few_effective_observations():
    y=np.random.default_rng(33).normal(size=30).cumsum()
    candidate=adf(y,trend='ct',lags=6)
    with patch.object(module,'adf',return_value=candidate):
        with pytest.raises(ValueError,match='fewer than 25'):
            enders(y,max_lags=8)
    with pytest.raises(ValueError,match='25'):
        _phi_critical('phi3',24,.05)


@pytest.mark.parametrize('boundary',['tau','phi','auxiliary','conditional'])
def test_equality_does_not_reject(boundary):
    y=np.random.default_rng(7).normal(size=100).cumsum()
    actual={tr:adf(y,trend=tr,lags=1) for tr in ['ct','c','n']}
    def fake_adf(y,*,trend,**kw):
        r=actual[trend]
        value=0.
        if trend=='ct':
            if boundary=='tau':value=r.critical_values[.05]
            elif boundary=='conditional':value=float(t.ppf(.05,r.regression.df_resid))
        return replace(r,statistic=value)
    def fake_joint(r):
        cv,_=_phi_critical('phi3' if r.trend=='ct' else 'phi1',r.nobs,.05)
        return (cv if boundary=='phi' else 100.),1.,()
    def fake_aux(r):
        fit,_=_auxiliary(actual[r.trend])
        value=float(t.ppf(.975,fit.df_resid)) if boundary=='auxiliary' else 100.
        return fit,value
    with patch.object(module,'adf',side_effect=fake_adf),patch.object(module,'_joint',side_effect=fake_joint),patch.object(module,'_auxiliary',side_effect=fake_aux):
        r=enders(y,lags=1)
    name={'tau':'tau3','phi':'phi3','auxiliary':'auxiliary_trend','conditional':'conditional_t'}[boundary]
    entry=next(h for h in r.history if h.test==name)
    assert entry.statistic==entry.critical_value
    assert not entry.reject
