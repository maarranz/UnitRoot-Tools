"""Experimental fixed-lag deterministic IC selection and interpretation contracts."""
from dataclasses import FrozenInstanceError, replace
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import unitroot_tools as u
from unitroot_tools import _deterministic as module


def fixture(offset=.1, slope=0.):
    t=np.arange(120,dtype=float)
    return offset+slope*t+np.sin(.8*t)+.6*np.cos(1.71*t)+.15*np.sin(.17*t*t)


@pytest.mark.parametrize('p',[0,1,3])
@pytest.mark.parametrize('criterion',['aic','bic','AIC','BIC'])
def test_independent_scores_and_sample(p,criterion):
    y=fixture(); original=y.copy()
    result=u.select_deterministic(y,lags=p,criterion=criterion)
    assert result.criterion==criterion.lower()
    assert (result.input_nobs,result.nobs,result.sample_start,result.sample_stop)==(120,119-p,p+1,120)
    assert [c.specification for c in result.candidates]==['n','c','ct']
    response=np.array([y[t]-y[t-1] for t in range(p+1,len(y))])
    for d,c in enumerate(result.candidates):
        x=np.array([[y[t-1]]+[y[t-j]-y[t-j-1] for j in range(1,p+1)]
                    +([1.] if d else [])+([t-p] if d==2 else []) for t in range(p+1,len(y))])
        beta=np.linalg.lstsq(x,response,rcond=None)[0]
        residual=response-x@beta
        rss=residual@residual; n=len(response); k=x.shape[1]+1
        np.testing.assert_allclose(c.adf_result.regression.params,beta,rtol=0,atol=1e-11)
        np.testing.assert_array_equal(c.adf_result.regression.model.exog,x)
        np.testing.assert_array_equal(c.adf_result.regression.model.endog,response)
        assert c.ssr==pytest.approx(rss,abs=1e-11)
        assert c.variance_mle==pytest.approx(rss/n,abs=1e-12)
        assert c.log_likelihood==pytest.approx(-n/2*(np.log(2*np.pi)+1+np.log(rss/n)),abs=1e-10)
        assert c.aic_score==pytest.approx(n*np.log(rss/n)+2*k,abs=1e-10)
        assert c.bic_score==pytest.approx(n*np.log(rss/n)+np.log(n)*k,abs=1e-10)
        assert (c.n_coefficients,c.n_parameters,c.df_resid)==(1+p+d,2+p+d,n-1-p-d)
        fixed=u.adf(y,trend=c.specification,lags=p)
        assert (fixed.statistic,fixed.pvalue,fixed.critical_values,fixed.calibration)==(
            c.adf_result.statistic,c.adf_result.pvalue,c.adf_result.critical_values,c.adf_result.calibration)
    assert result.selected_fixed_spec_adf is result.selected_candidate.adf_result
    scores=[getattr(c,result.criterion+'_score') for c in result.candidates]
    ordered=sorted(scores)
    assert result.criterion_value==min(scores)
    assert result.runner_up_gap==ordered[1]-ordered[0]
    assert dict(result.score_differences)=={c.specification:s-min(scores) for c,s in zip(result.candidates,scores)}
    np.testing.assert_array_equal(y,original)


def test_real_disagreement_and_trend():
    assert u.select_deterministic(fixture(),lags=1,criterion='aic').selected_specification=='c'
    assert u.select_deterministic(fixture(),lags=1).selected_specification=='n'
    for criterion in ('aic','bic'):
        assert u.select_deterministic(fixture(slope=.03),lags=1,criterion=criterion).selected_specification=='ct'


def test_analytic_fixture():
    rows=[module._scores(rss,100,k) for rss,k in zip(
        [100,100*np.exp(-.03),100*np.exp(-.035)],[3,4,5])]
    np.testing.assert_allclose([r[2] for r in rows],[6,5,6.5],rtol=0,atol=1e-12)
    np.testing.assert_allclose([r[3] for r in rows],[3*np.log(100),4*np.log(100)-3,5*np.log(100)-3.5],rtol=0,atol=1e-12)
    for k,row in zip([3,4,5],rows):
        assert row[2]==pytest.approx(-2*row[1]+2*k-100*(1+np.log(2*np.pi)),abs=1e-12)


@pytest.mark.parametrize('criterion',['aic','bic'])
def test_exact_and_near_ties(criterion):
    base=u.select_deterministic(fixture(),lags=1)
    tied=tuple(replace(c,aic_score=1.,bic_score=1.) for c in base.candidates)
    # Order independence; a genuine exact tie favors the smallest dimension.
    assert module._winner(tied[::-1],criterion).specification=='n'
    near=list(tied); near[2]=replace(near[2],**{criterion+'_score':np.nextafter(1.,0.)})
    assert module._winner(tuple(near),criterion).specification=='ct'
    result=replace(base,criterion=criterion,selected_specification='n',candidates=tied)
    assert result.runner_up_gap==0 and set(result.score_differences.values())=={0.}


@pytest.mark.parametrize('scale',[-3.,.001,100.])
def test_scale_equivariance(scale):
    a=u.select_deterministic(fixture(),lags=1)
    b=u.select_deterministic(fixture()*scale,lags=1)
    assert a.selected_specification==b.selected_specification
    for ca,cb in zip(a.candidates,b.candidates):
        for key in ('aic_score','bic_score'):
            assert getattr(cb,key)-getattr(ca,key)==pytest.approx(2*a.nobs*np.log(abs(scale)),abs=1e-9)


@pytest.mark.parametrize('lags',[None,True,np.bool_(False),1.0,'1',-1])
def test_invalid_lags(lags):
    with pytest.raises((TypeError,ValueError),match='lags'):
        u.select_deterministic(fixture(),lags=lags)


@pytest.mark.parametrize('criterion',['maic','mbic','invalid',None,1])
def test_invalid_criterion(criterion):
    with pytest.raises((TypeError,ValueError),match='criterion'):
        u.select_deterministic(fixture(),lags=1,criterion=criterion)


@pytest.mark.parametrize('kwargs',[{}, {'lags':1,'max_lags':2},{'lags':1,'method':'bic'},
                                    {'lags':1,'alpha':.05},{'lags':1,'candidates':['c','ct']}])
def test_public_signature(kwargs):
    with pytest.raises(TypeError): u.select_deterministic(fixture(),**kwargs)


@pytest.mark.parametrize('y',[[],[1.,np.nan,2.],[1.,np.inf,2.],np.ones((10,1)),[1,True,2],
                              [1+1j,2],['1','2'],np.array([1,2],dtype=object)])
def test_invalid_series(y):
    with pytest.raises((TypeError,ValueError)): u.select_deterministic(y,lags=0)


@pytest.mark.parametrize('p',[0,1,3])
def test_insufficient_sample(p):
    with pytest.raises(ValueError,match='candidate ct: insufficient sample'):
        u.select_deterministic(fixture()[:2*p+4],lags=p)


@pytest.mark.parametrize('y',[np.ones(30),np.zeros(30),np.arange(30.),2.**np.arange(30)])
def test_degenerate(y):
    with pytest.raises((u.EstimationError,ValueError),match='Deterministic candidate [nct]+:'):
        u.select_deterministic(y,lags=0)


@pytest.mark.parametrize('bad_spec',['n','c','ct'])
@pytest.mark.parametrize('exception',[TypeError,ValueError,u.EstimationError])
def test_failure_identifies_candidate(monkeypatch,bad_spec,exception):
    actual=module.adf; calls=[]
    def broken(y,*,trend,lags):
        calls.append(trend)
        if trend==bad_spec: raise exception('controlled failure')
        return actual(y,trend=trend,lags=lags)
    monkeypatch.setattr(module,'adf',broken)
    with pytest.raises(exception,match=f'candidate {bad_spec}: controlled failure'):
        u.select_deterministic(fixture(),lags=1)
    assert calls==list(('n','c','ct')[:('n','c','ct').index(bad_spec)+1])


@pytest.mark.parametrize('corruption',['metadata','response','names','values','rank','df','inference'])
def test_alignment_guard(monkeypatch,corruption):
    actual=module.adf
    def broken(y,*,trend,lags):
        r=actual(y,trend=trend,lags=lags)
        if trend!='c': return r
        if corruption=='metadata': return replace(r,sample_start=r.sample_start+1)
        if corruption=='inference': return replace(r,pvalue=float('nan'))
        reg=r.regression
        x=reg.model.exog.copy(); response=reg.model.endog.copy();names=list(reg.model.exog_names)
        df=reg.df_resid
        if corruption=='response': response[0]+=1
        if corruption=='names': names.reverse()
        if corruption=='values': x[0,0]+=1
        if corruption=='rank': x[:,1]=x[:,0]
        if corruption=='df': df+=1
        fake=SimpleNamespace(model=SimpleNamespace(exog=x,endog=response,exog_names=names),nobs=reg.nobs,df_resid=df)
        return replace(r,regression=fake)
    monkeypatch.setattr(module,'adf',broken)
    with pytest.raises(u.EstimationError,match='candidate c:'):
        u.select_deterministic(fixture(),lags=1)


@pytest.mark.parametrize('ssr',[0.,-1.,float('nan'),float('inf'),np.nextafter(0.,1.)])
def test_score_variance_guard(ssr):
    with pytest.raises(u.EstimationError,match='variance'): module._scores(ssr,100,3)


def test_provenance_summary_and_frozen_records():
    y=pd.Series(fixture(),name='economic series')
    r=u.select_deterministic(y,lags=np.int64(1))
    assert r.input_name=='economic series' and r.experimental
    assert r.lag_policy=='supplied_common_integer'
    assert r.criterion_normalization=='total_omit_common_gaussian_constant_include_variance'
    assert r.tie_break=='exact_score_then_smaller_deterministic_dimension'
    assert r.inference_status=='fixed_specification_calibration_after_selection'
    assert dict(r.versions)==dict(r.selected_fixed_spec_adf.versions)
    for attribute in ('pvalue','reject','reject_unit_root'): assert not hasattr(r,attribute)
    summary=r.summary()
    for text in ('Experimental','fixed p=1','Fixed-specification ADF p-value reported after selection',
                 'upstream lag','not a selective p-value','nominal size','Calibration:', 'Including n'):
        assert text in summary
    assert 'Decision at 5%' not in summary
    with pytest.raises(FrozenInstanceError): r.lags=2
    with pytest.raises(FrozenInstanceError): r.candidates[0].ssr=1.
    with pytest.raises(TypeError): r.versions['arch']='bad'
    with pytest.raises(TypeError): r.score_differences['n']=1.
    assert {'select_deterministic','DeterministicCandidate','DeterministicSelectionResult'}<=set(u.__all__)


def test_only_three_fixed_fits_and_warning_propagation(monkeypatch):
    import warnings
    actual=module.adf; calls=[]
    def tracked(y,*,trend,lags):
        calls.append((trend,lags))
        if trend=='c': warnings.warn('inherited conditioning warning',u.SpecificationWarning)
        return actual(y,trend=trend,lags=lags)
    monkeypatch.setattr(module,'adf',tracked)
    with pytest.warns(u.SpecificationWarning,match='inherited conditioning warning'):
        r=u.select_deterministic(fixture(),lags=1)
    assert calls==[('n',1),('c',1),('ct',1)]
    assert r.selected_fixed_spec_adf is r.selected_candidate.adf_result


@pytest.mark.parametrize('p',[0,1,2])
def test_algebraic_sample_boundary(p):
    # Only feasibility is asserted; this is not a recommended inferential sample.
    y=fixture()[:2*p+5]
    r=u.select_deterministic(y,lags=p)
    assert r.candidates[-1].df_resid==1


@pytest.mark.parametrize('scale',[1e-12,-1e-12])
@pytest.mark.parametrize('criterion',['aic','bic'])
def test_small_scale_preserves_valid_fits_scores_and_ranking(scale,criterion):
    y=fixture()
    original=u.select_deterministic(y,lags=1,criterion=criterion)
    scaled=u.select_deterministic(scale*y,lags=1,criterion=criterion)
    assert original.selected_specification==scaled.selected_specification
    n=original.nobs
    # ARCH/Statsmodels solves the raw design: at these extreme units the ct
    # residual RSS differs by ~8.7e-6 relatively despite normalized full rank.
    # A 1e-5 RSS budget implies |n log(RSS_scaled/(s² RSS))| <=
    # -n log(1-1e-5), about .00118 score units here; two scores bound a gap.
    # Ordinary-scale tolerances in the existing tests remain unchanged.
    rss_rtol=1e-5
    score_atol=-n*np.log1p(-rss_rtol)+1e-9
    for a,b in zip(original.candidates,scaled.candidates):
        fixed=u.adf(scale*y,trend=b.specification,lags=1)
        assert fixed.regression.df_resid==n-b.n_coefficients
        assert fixed.regression.ssr>0 and np.isfinite(fixed.regression.ssr)
        assert fixed.statistic==b.adf_result.statistic
        x=fixed.regression.model.exog
        normalized=x/np.linalg.norm(x,axis=0)
        assert np.linalg.matrix_rank(normalized)==b.n_coefficients
        assert np.linalg.cond(normalized)<4
        assert b.ssr/(scale*scale*a.ssr)==pytest.approx(1.,rel=0,abs=rss_rtol)
        for key in ('aic_score','bic_score'):
            assert getattr(b,key)-getattr(a,key)==pytest.approx(
                2*n*np.log(abs(scale)),rel=0,abs=score_atol)
    for key in ('aic_score','bic_score'):
        order=lambda r: [c.specification for c in sorted(r.candidates,key=lambda c:getattr(c,key))]
        assert order(original)==order(scaled)
    assert scaled.runner_up_gap==pytest.approx(original.runner_up_gap,rel=0,abs=2*score_atol)
