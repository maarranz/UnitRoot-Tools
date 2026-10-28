"""Independent QR, source-table and decision-tree reconstruction for Experiment A."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import stats
from arch.unitroot import ADF
from unitroot_tools import adf
from .engine import ALPHAS, METHODS, SPECS, evaluate_series
from .catalogue import catalogue, innovations, generate

# Independently transcribed Enders manual Table B (T rows; 1%,5%,10%).
PHI = {100: {'ct':(8.73,6.49,5.47),'c':(6.70,4.71,3.86)},
       250: {'ct':(8.43,6.34,5.39),'c':(6.52,4.63,3.81)},
       500: {'ct':(8.34,6.30,5.36),'c':(6.47,4.61,3.79)}}


def qr(X,Y):
    q,r=np.linalg.qr(X,mode='reduced')
    b=np.linalg.solve(r,q.T@Y)
    residual=Y-X@b
    rss=float(residual@residual)
    df=len(Y)-X.shape[1]
    inv=np.linalg.solve(r,np.eye(len(r)))
    cov=rss/df*(inv@inv.T)
    return dict(b=b,rss=rss,df=df,cov=cov,t=b/np.sqrt(np.diag(cov)),X=X,Y=Y)


def reconstruct(y):
    T=len(y)-1
    Y=np.array([y[t]-y[t-1] for t in range(1,T+1)])
    level=np.array([y[t-1] for t in range(1,T+1)])
    refs={}
    for d,spec in enumerate(SPECS):
        columns=[level]
        if d:columns.append(np.ones(T))
        if d==2:columns.append(np.arange(1,T+1,dtype=float))
        ref=qr(np.column_stack(columns),Y)
        # Calibration is intentionally the released public ARCH calibration,
        # independently requested, not a new theoretical calibration claim.
        backend=ADF(y,lags=0,trend=spec)
        ref['critical']={a:float(backend.critical_values[k]) for a,k in zip(ALPHAS,('1%','5%','10%'))}
        k=d+2
        sigma2=ref['rss']/T
        ref['aic']=T*np.log(sigma2)+2*k
        ref['bic']=T*np.log(sigma2)+np.log(T)*k
        ref['llf']=float(stats.norm.logpdf(Y-ref['X']@ref['b'],scale=np.sqrt(sigma2)).sum())
        if spec!='n':
            aux=qr(ref['X'][:,1:],Y)
            restricted=qr(np.ones((T,1)),Y)['rss'] if spec=='ct' else float(Y@Y)
            ref['aux']=aux
            ref['phi']=(restricted-ref['rss'])/2/(ref['rss']/ref['df'])
            indices=[0,ref['X'].shape[1]-1]
            b=ref['b'][indices];cov=ref['cov'][np.ix_(indices,indices)]
            ref['phi_wald']=float(b@np.linalg.solve(cov,b)/2)
        refs[spec]=ref
    return refs


def expected_history(refs,T,alpha,conservative):
    ai=ALPHAS.index(alpha)
    history=[]
    def step(spec,test,value,critical,tail):
        reject=bool(abs(value)>critical if tail=='two-sided' else value>critical if tail=='right' else value<critical)
        history.append(dict(specification=spec,test=test,statistic=float(value),critical_value=float(critical),reject=reject,tail=tail))
        return reject
    for spec,tau,phi in [('ct','tau3','phi3'),('c','tau2','phi1'),('n','tau1',None)]:
        ref=refs[spec]
        rejected=step(spec,tau,ref['t'][0],ref['critical'][alpha],'left')
        if rejected or spec=='n' or (conservative and spec=='c'):
            return history
        if not step(spec,phi,ref['phi'],PHI[T][spec][ai],'right'):
            continue
        aux=ref['aux']
        name='auxiliary_trend' if spec=='ct' else 'auxiliary_intercept'
        if not step(spec,name,aux['t'][-1],stats.t.ppf(1-alpha/2,aux['df']),'two-sided'):
            continue
        step(spec,'conditional_t',ref['t'][0],stats.t.ppf(alpha,ref['df']),'left')
        return history
    raise AssertionError('unreachable reference tree')


def validate():
    maxima={k:0. for k in ['coefficient','ssr','statistic','phi','auxiliary','aic','bic','log_likelihood','critical_value']}
    tol=1e-8
    def check(key,a,b):
        error=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
        maxima[key]=max(maxima[key],error)
        assert error<=tol,(key,error)
    seen=set();cases=[]
    for T in (100,250,500):
        cells=[c for c in catalogue() if c['law']['T']==T]
        examples=[next(c for c in cells if c['is_null'] and c['intercept']==0),
                  next(c for c in cells if c['is_null'] and c['intercept']==.5),
                  next(c for c in cells if not c['is_null'] and c['rho']==.9 and c['intercept']==0)]
        fixtures=[(c['id'],generate(c,innovations(c,19,'validation'))) for c in examples]
        t=np.arange(T+1,dtype=float)
        fixtures.append(('fixed_curvature',.002*t*t+np.cumsum(np.sin(1.23*t)+.2*np.cos(.31*t*t))))
        for name,y in fixtures:
            refs=reconstruct(y)
            actual=evaluate_series(y,cell_id=name,equivalence=True)
            assert not actual['failures']
            assert actual['equivalence_checks']==3
            for spec,ref in refs.items():
                r=adf(y,trend=spec,lags=0)
                assert (r.lags,r.nobs,r.sample_start,r.sample_stop)==(0,T,1,T+1)
                np.testing.assert_array_equal(r.regression.model.exog,ref['X'])
                np.testing.assert_array_equal(r.regression.model.endog,ref['Y'])
                check('coefficient',r.regression.params,ref['b'])
                check('ssr',r.regression.ssr,ref['rss'])
                check('statistic',r.statistic,ref['t'][0])
                for alpha in ALPHAS:
                    check('critical_value',r.critical_values[alpha],ref['critical'][alpha])
                    assert actual['decisions'][ALPHAS.index(alpha)][METHODS.index('adf_'+spec)]==int(ref['t'][0]<ref['critical'][alpha])
                if spec!='n':check('phi',ref['phi'],ref['phi_wald'])
            for criterion in ['aic','bic']:
                ordered=sorted(SPECS,key=lambda s:(refs[s][criterion],SPECS.index(s)))
                assert actual['ic'][criterion]['selected']==ordered[0]
                assert abs(actual['ic'][criterion]['gap']-(refs[ordered[1]][criterion]-refs[ordered[0]][criterion]))<=tol
                for c in actual['ic'][criterion]['scores']:
                    ref=refs[c['specification']]
                    check('aic',c['aic'],ref['aic']);check('bic',c['bic'],ref['bic'])
                    check('log_likelihood',c['log_likelihood'],ref['llf'])
                    assert c['k']==SPECS.index(c['specification'])+2
                for i,alpha in enumerate(ALPHAS):
                    assert actual['decisions'][i][METHODS.index('det_'+criterion)]==int(refs[ordered[0]]['t'][0]<refs[ordered[0]]['critical'][alpha])
            for alpha in ALPHAS:
                for variant in ['original','conservative']:
                    expected=expected_history(refs,T,alpha,variant=='conservative')
                    observed=actual['enders'][str(alpha)][variant]
                    assert len(expected)==len(observed['steps'])
                    for e,o in zip(expected,observed['steps']):
                        assert (e['specification'],e['test'],e['tail'],e['reject'])==(o['specification'],o['test'],o['tail'],o['reject'])
                        check('auxiliary' if e['test'].startswith('auxiliary') else 'phi' if e['test'].startswith('phi') else 'statistic',o['statistic'],e['statistic'])
                        check('critical_value',o['critical_value'],e['critical_value'])
                        if e['test'].startswith('auxiliary'):
                            assert o['reference_df']==T-(2 if e['specification']=='ct' else 1)
                        if e['test']=='conditional_t':
                            assert o['reference_df']==T-(3 if e['specification']=='ct' else 2)
                        seen.add((e['specification'],e['test'],e['reject']))
                    assert observed['selected']==expected[-1]['specification']
                    assert observed['reject']==expected[-1]['reject']
                    assert observed['terminal']==expected[-1]['test']
            cases.append(dict(T=T,fixture=name,failures=0))
    return dict(all_passed=True,cases=cases,absolute_tolerance=tol,max_absolute_errors=maxima,
                history_edges_seen=sorted(seen),calibration='Public released ARCH for tau; independent Enders Table B transcription; Student quantiles checked with correct independent df',
                scope='Fixed-vector numerical/decision validation; not pilot rejection rates or size evidence')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report',type=Path,required=True)
    args=p.parse_args()
    if not __debug__:raise RuntimeError('Assertions required')
    report=validate();args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['cases','history_edges_seen']},indent=2))


if __name__=='__main__':main()
