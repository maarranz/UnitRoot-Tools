"""Independent QR checks for experimental deterministic AIC/BIC selection.

Run with PYTHONPATH=src, or against an isolated installed package:
    python validation/validate_deterministic.py --report /tmp/deterministic.json
Fixed analytic regression fixtures only; no Monte Carlo size/power experiments.
"""
import argparse
import json
from pathlib import Path
import platform

import numpy as np
from scipy.stats import norm
from unitroot_tools import adf, select_deterministic


def fixture(offset=0.1, slope=0.0):
    t = np.arange(120, dtype=float)
    return offset + slope*t + np.sin(.8*t) + .6*np.cos(1.71*t) + .15*np.sin(.17*t*t)


def reconstruct(y, p, d):
    # Construct observations by explicit time indices, independently of selector slicing.
    times = range(p+1, len(y))
    response = np.array([y[t]-y[t-1] for t in times])
    rows = []
    for i, t in enumerate(times):
        row = [y[t-1]] + [y[t-j]-y[t-j-1] for j in range(1, p+1)]
        if d >= 1:
            row.append(1.)
        if d == 2:
            row.append(float(i+1))
        rows.append(row)
    x = np.array(rows)
    q, r = np.linalg.qr(x, mode='reduced')
    coefficients = np.linalg.solve(r, q.T@response)
    residual = response-x@coefficients
    rss = float(residual@residual)
    n = len(response)
    k = x.shape[1]+1
    variance = rss/n
    # Sum Gaussian log densities instead of using selector's closed form.
    llf = float(norm.logpdf(residual, scale=np.sqrt(variance)).sum())
    common = n*(1+np.log(2*np.pi))
    aic = -2*llf+2*k-common
    bic = -2*llf+np.log(n)*k-common
    return dict(x=x, response=response, coefficients=coefficients, ssr=rss,
                variance_mle=variance, log_likelihood=llf, aic_score=aic, bic_score=bic,
                n_parameters=k, n_coefficients=x.shape[1], df_resid=n-x.shape[1])


def validate():
    maxima = {k: 0. for k in ('coefficients','ssr','variance_mle','log_likelihood',
                              'aic_score','bic_score','runner_up_gap','score_differences')}
    # Absolute thresholds fixed in advance, including small-sample p=0 and positive lags.
    tolerances = {k: 1e-9 for k in maxima}
    cases = []
    for name, y in [('disagreement', fixture()), ('trend', fixture(slope=.03))]:
        for p in (0, 1, 3):
            expected = {spec: reconstruct(y, p, d) for d, spec in enumerate(('n','c','ct'))}
            for criterion in ('aic','bic'):
                result = select_deterministic(y, lags=p, criterion=criterion)
                assert (result.sample_start,result.sample_stop,result.nobs)==(p+1,len(y),len(y)-p-1)
                scores = {s: row[criterion+'_score'] for s,row in expected.items()}
                ordering = sorted(scores, key=lambda s: (scores[s], ('n','c','ct').index(s)))
                assert result.selected_specification == ordering[0]
                def compare(key, actual, reference):
                    error = float(np.max(np.abs(np.asarray(actual)-np.asarray(reference))))
                    maxima[key] = max(maxima[key], error)
                    assert error <= tolerances[key], (key, error, tolerances[key])
                compare('runner_up_gap',result.runner_up_gap,scores[ordering[1]]-scores[ordering[0]])
                for c in result.candidates:
                    ref = expected[c.specification]
                    reg = c.adf_result.regression
                    assert np.array_equal(reg.model.exog,ref['x'])
                    assert np.array_equal(reg.model.endog,ref['response'])
                    compare('coefficients',reg.params,ref['coefficients'])
                    for key in ('ssr','variance_mle','log_likelihood','aic_score','bic_score'):
                        compare(key,getattr(c,key),ref[key])
                    for key in ('n_parameters','n_coefficients','df_resid'):
                        assert getattr(c,key)==ref[key]
                    compare('score_differences',result.score_differences[c.specification],scores[c.specification]-scores[ordering[0]])
                    fixed = adf(y,trend=c.specification,lags=p)
                    assert c.adf_result.statistic == fixed.statistic
                    assert c.adf_result.pvalue == fixed.pvalue
                    assert c.adf_result.critical_values == fixed.critical_values
                    assert c.adf_result.calibration == fixed.calibration
                assert result.selected_fixed_spec_adf is result.selected_candidate.adf_result
                cases.append(dict(fixture=name,lags=p,criterion=criterion,selected=result.selected_specification,
                                  runner_up_gap=result.runner_up_gap,scores=scores))
    a = next(c for c in cases if c['fixture']=='disagreement' and c['lags']==1 and c['criterion']=='aic')
    b = next(c for c in cases if c['fixture']=='disagreement' and c['lags']==1 and c['criterion']=='bic')
    assert (a['selected'],b['selected'])==('c','n')
    assert all(c['selected']=='ct' for c in cases if c['fixture']=='trend')
    # Section 9.1 is an arithmetic fixture, not three estimated regressions.
    n=100
    rss=np.array([100.,100*np.exp(-.03),100*np.exp(-.035)])
    k=np.array([3,4,5])
    aic=n*np.log(rss/n)+2*k
    bic=n*np.log(rss/n)+np.log(n)*k
    np.testing.assert_allclose(aic,[6,5,6.5],rtol=0,atol=1e-12)
    assert np.argmin(aic)==1 and np.argmin(bic)==0
    return dict(all_passed=True,python=platform.python_version(),architecture=platform.machine(),
                cases=cases,regressions=18,selection_calls=12,max_absolute_differences=maxima,
                absolute_tolerances=tolerances,arithmetic_fixture=dict(n=n,rss=rss.tolist(),k=k.tolist(),
                aic=aic.tolist(),bic=bic.tolist()),
                scope='Numerical implementation validation only; no nominal-size or power claim.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    if not __debug__:
        raise RuntimeError('Assertions required; do not use -O')
    result=validate()
    args.report.write_text(json.dumps(result,indent=2)+'\n')
    print('Independent deterministic-selection QR validation passed (12 selections, 18 distinct fits).')
    print(json.dumps(result['max_absolute_differences'],indent=2))


if __name__=='__main__':
    main()
