"""Independent fixed-lag IC verification on saved binary64 inputs; no RNG.

Uses 80/120-digit Decimal elimination from the independent Enders validator,
then Decimal logarithms for Gaussian scores. Production fitting is used only
for the actual values under test. Never generates or overwrites input fixtures.
"""
from decimal import Decimal as D, localcontext
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
from unitroot_tools import select_deterministic

ROOT = Path(__file__).resolve().parents[1]
loader = importlib.util.spec_from_file_location('independent_enders', ROOT/'validation/enders/verify_scaled_references.py')
independent = importlib.util.module_from_spec(loader)
loader.loader.exec_module(independent)


def verify_case(key, p, criterion):
    input_path = ROOT/'validation/scaled_svd/adf_inputs.json'
    assert hashlib.sha256(input_path.read_bytes()).hexdigest() == '0a5ef87c88547da466d37b1368e195037fe010de7a32a75e9234d8f2524ae8d4'
    inputs = json.loads(input_path.read_text())['inputs']
    y = np.array([float.fromhex(v) for v in inputs[key]])
    actual = select_deterministic(y, lags=p, criterion=criterion)
    n = len(y)-1-p
    response = np.array([y[t]-y[t-1] for t in range(p+1,len(y))])
    rows = []
    for d, candidate in enumerate(actual.candidates):
        x = np.array([[y[t-1]]+[y[t-j]-y[t-j-1] for j in range(1,p+1)]
                      +([1.] if d else [])+([float(t-p)] if d==2 else [])
                      for t in range(p+1,len(y))])
        refs = [independent.reference(x,response,precision) for precision in (80,120)]
        assert float(refs[0]['rss']) == float(refs[1]['rss'])
        expected = []
        for precision, ref in zip((80,120),refs):
            with localcontext() as ctx:
                ctx.prec = precision
                variance = D(ref['rss'])/D(n)
                base = D(n)*variance.ln()
                # 120+ digits of pi, independent of binary64 math.pi.
                pi = D('3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446')
                common = D(n)*(1+(2*pi).ln())
                k = x.shape[1]+1
                expected.append(dict(ssr=float(ref['rss']), variance_mle=float(variance),
                                     log_likelihood=float(-(base+common)/2),
                                     aic_score=float(base+2*k), bic_score=float(base+D(n).ln()*k)))
        assert expected[0] == expected[1]
        ref = expected[1]
        reg = candidate.adf_result.regression
        assert np.array_equal(reg.model.exog,x) and np.array_equal(reg.model.endog,response)
        assert (candidate.n_coefficients,candidate.n_parameters,candidate.df_resid)==(1+p+d,2+p+d,n-1-p-d)
        assert (actual.sample_start,actual.sample_stop,actual.nobs)==(p+1,len(y),n)
        assert candidate.ssr == float(reg.ssr)
        assert candidate.variance_mle == candidate.ssr/n
        assert reg.scale == candidate.ssr/candidate.df_resid
        common = n*(1+np.log(2*np.pi))
        np.testing.assert_allclose(candidate.log_likelihood,reg.llf,rtol=0,atol=2e-12)
        for score in ('aic','bic'):
            # Statsmodels default excludes sigma²: request dk_params=1 explicitly.
            np.testing.assert_allclose(getattr(candidate,score+'_score'),
                                       reg.info_criteria(score,dk_params=1)-common,rtol=0,atol=4e-12)
        errors = {name:abs(getattr(candidate,name)-value) for name,value in ref.items()}
        for name in ('ssr','variance_mle'):
            errors[name] /= abs(ref[name])
            assert errors[name] <= 1e-12,(key,p,name,errors[name])
        for name in ('log_likelihood','aic_score','bic_score'):
            assert errors[name] <= 4e-12,(key,p,name,errors[name])
        rows.append(dict(specification=candidate.specification,reference=ref,errors=errors))
    order = sorted(rows,key=lambda row:(row['reference'][criterion+'_score'],('n','c','ct').index(row['specification'])))
    assert actual.selected_specification == order[0]['specification']
    scores = {row['specification']:row['reference'][criterion+'_score'] for row in rows}
    np.testing.assert_allclose(actual.runner_up_gap,scores[order[1]['specification']]-scores[order[0]['specification']],rtol=0,atol=8e-12)
    for spec, difference in actual.score_differences.items():
        np.testing.assert_allclose(difference,scores[spec]-scores[actual.selected_specification],rtol=0,atol=8e-12)
    assert actual.selected_fixed_spec_adf is actual.selected_candidate.adf_result
    return dict(key=key,lags=p,criterion=criterion,selected=actual.selected_specification,rows=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report',type=Path,required=True);args=p.parse_args()
    cases=[verify_case(key,lag,criterion) for key in ('y','positive','negative') for lag in (0,1,4) for criterion in ('aic','bic')]
    maxima={name:max(row['errors'][name] for case in cases for row in case['rows']) for name in cases[0]['rows'][0]['errors']}
    with args.report.open('x') as f:json.dump(dict(cases=cases,max_errors=maxima,distinct_fits=27,selection_calls=18,precision_80_120_binary64_equal=True),f,indent=2);f.write('\n')
    print(json.dumps(maxima,indent=2))

if __name__=='__main__':main()
