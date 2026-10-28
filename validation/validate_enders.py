"""Independent QR regression checks and Enders's rounded published GDP example.

Run against an installed package or with PYTHONPATH=src. Writes only --report.
No simulation of size/power and no changes to historical evidence.
"""
import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform

import numpy as np
from unitroot_tools import adf, enders
from unitroot_tools._enders import _joint, _auxiliary

ROOT = Path(__file__).resolve().parent/'enders'


def qr(x, y):
    q,r=np.linalg.qr(x,mode='reduced')
    b=np.linalg.solve(r,q.T@y)
    e=y-x@b
    ssr=float(e@e)
    inverse=np.linalg.solve(r,np.eye(len(r)))
    covariance=ssr/(len(y)-len(r))*(inverse@inverse.T)
    return b,np.sqrt(covariance.diagonal()),ssr,covariance


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    if not __debug__:raise RuntimeError('Assertions required; do not use -O')
    metadata=json.loads((ROOT/'source.json').read_text())
    assert hashlib.sha256((ROOT/'rgdp.csv').read_bytes()).hexdigest()==metadata['csv_sha256']
    with (ROOT/'rgdp.csv').open() as f: levels=np.array([float(r['rgdp']) for r in csv.DictReader(f)])
    y=np.log(levels)
    # Controlled deterministic vectors are regression checks, not MC trials.
    vectors=[y, np.sin(np.arange(120)*.31)+.001*np.arange(120)**2+.03*np.random.default_rng(71).normal(size=120),
             np.cumsum(np.random.default_rng(72).normal(size=180))]
    maxima=dict(coefficient=0.,standard_error=0.,ssr=0.,t_statistic=0.,phi=0.,aux_coefficient=0.,aux_standard_error=0.,aux_ssr=0.,aux_t=0.)
    cases=0
    def check(name,actual,reference):
        error=float(np.max(np.abs(np.asarray(actual)-np.asarray(reference))))
        maxima[name]=max(maxima[name],error)
        assert error<1e-8,(name,error)
    for v in vectors:
        for trend in ['ct','c']:
            for p in [0,1,4]:
                d=np.diff(v); response=d[p:]
                columns=[v[p:-1]]+[d[p-j:-j] for j in range(1,p+1)]+[np.ones(len(response))]
                if trend=='ct':columns.append(np.arange(1,len(response)+1))
                x=np.column_stack(columns)
                b,se,ssr,cov=qr(x,response)
                a=adf(v,trend=trend,lags=p);reg=a.regression
                for name,actual,reference in [('coefficient',reg.params,b),('standard_error',reg.bse,se),('ssr',reg.ssr,ssr),('t_statistic',reg.tvalues,b/se)]:check(name,actual,reference)
                restrictions=np.eye(x.shape[1])[[0,-1]]
                contrast=restrictions@b
                wald=float(contrast@np.linalg.solve(restrictions@cov@restrictions.T,contrast)/2)
                phi,rss_r,_=_joint(a)
                check('phi',phi,wald)
                keep=[i for i in range(x.shape[1]) if i not in [0,x.shape[1]-1]]
                reference_rss=qr(x[:,keep],response)[2] if keep else response@response
                assert abs(rss_r-reference_rss)<1e-8
                check('phi',phi,(reference_rss-ssr)/2/(ssr/(len(response)-x.shape[1])))
                aux,stat=_auxiliary(a); ab,ase,arss,_=qr(x[:,1:],response)
                check('aux_coefficient',aux.params,ab);check('aux_standard_error',aux.bse,ase)
                check('aux_ssr',aux.ssr,arss);check('aux_t',stat,ab[-1]/ase[-1])
                cases+=1
    original=enders(y,lags=1)
    conservative=enders(y,lags=1,variant='conservative')
    permitted=enders(y,lags=1,variant='conservative',allow_no_constant=True)
    published=[-1.49,2.97,-2.05,25.49,6.80,-2.05]
    rows=[]
    for h,reference in zip(original.history,published):
        assert abs(h.statistic-reference)<=.005
        rows.append(dict(test=h.test,statistic=h.statistic,published_statistic=reference,
                         absolute_difference=abs(h.statistic-reference),critical_value=h.critical_value,
                         reject=h.reject,next_action=h.next_action))
    assert len(rows)==6 and original.reject_unit_root and permitted.reject_unit_root
    assert original.trend==conservative.trend==permitted.trend=='c'
    assert not conservative.reject_unit_root and conservative.terminal_test=='tau2'
    assert original.history[1].critical_value==6.34 and original.history[3].critical_value==4.63
    general=original.history[0].regression
    coefficient_rows=[]
    groups=[('A2.2',general,{'const':.1248,'trend':.0001,'Level.L1':-.0156,'Diff.L1':.3663}),
            ('A2.3',original.history[2].regression,{'const':.0215,'Level.L1':-.0019,'Diff.L1':.3539}),
            ('A2.1',original.history[4].regression,{'const':.0049,'Diff.L1':.3706})]
    for equation,reg,expected in groups:
        for term,ref in expected.items():
            value=float(reg.params[term])
            if equation=='A2.2' and term=='const':value-=2*float(reg.params['trend'])
            assert abs(value-ref)<.00005
            coefficient_rows.append(dict(equation=equation,term=term,value=value,published=ref))
    t_rows=[]
    for equation,reg,expected in [
        ('A2.2',general,{'const':1.58,'trend':1.31,'Level.L1':-1.49,'Diff.L1':6.26}),
        ('A2.3',original.history[2].regression,{'const':2.64,'Level.L1':-2.05,'Diff.L1':6.12}),
        ('A2.1',original.history[4].regression,{'const':6.80,'Diff.L1':6.44})]:
        values=np.asarray(reg.tvalues)
        if equation=='A2.2':
            # Express the constant on the publication's original t=3,...,264 index.
            transform=np.eye(len(reg.params))
            names=list(reg.model.exog_names)
            transform[names.index('const'),names.index('trend')]=-2
            beta=transform@np.asarray(reg.params)
            cov=transform@np.asarray(reg.cov_params())@transform.T
            values=beta/np.sqrt(np.diag(cov))
        for term,reference in expected.items():
            value=float(values[list(reg.model.exog_names).index(term)])
            difference=abs(value-reference)
            match=difference<=.005
            t_rows.append(dict(equation=equation,term=term,value=value,published=reference,
                               absolute_difference=difference,within_rounding=match))
            if (equation,term)!=('A2.3','const'):
                assert match,(equation,term,difference)
    # Retain the printed intercept-t discrepancy explicitly, without relaxing tolerance.
    discrepancies=[row for row in t_rows if not row['within_rounding']]
    assert [(row['equation'],row['term']) for row in discrepancies]==[('A2.3','const')]
    report=dict(platform=platform.platform(),python=platform.python_version(),
        versions={x:importlib.metadata.version(x) for x in ['unitroot-tools','arch','numpy','scipy','statsmodels','pandas']},
        regression_cases=cases,absolute_tolerance=1e-8,max_absolute_errors=maxima,
        published_gdp=rows,published_coefficients=coefficient_rows,
        published_t_statistics=t_rows,published_t_discrepancies=discrepancies,
        published_t_statistics_all_within_rounding=not discrepancies,
        gdp_nobs=original.nobs,published_rounding_tolerance=.005,
        original_reject=original.reject_unit_root,conservative_reject=conservative.reject_unit_root,
        conservative_permitted_reject=permitted.reject_unit_root,
        limitations=['Source text mentions 244 and 262 observations; 264 quarterly levels produce 262 regression observations.',
                    'Published tau tables are rounded; adapter retains ARCH response surfaces.',
                    'Conditional left tail follows assignment b<0; reference says t-distribution without explicit tail.',
                    'No Monte Carlo size/power study; these checks establish numerical and branch correctness.'],
        all_passed=True)
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
