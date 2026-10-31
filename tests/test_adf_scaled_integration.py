"""Stable ADF accuracy, econometric definition and separate calibration contracts."""
import json
import math
from pathlib import Path

import numpy as np
import pytest
from arch.unitroot import ADF
from arch.unitroot.unitroot import mackinnonp as arch_mackinnonp
from arch.unitroot.unitroot import mackinnoncrit as arch_mackinnoncrit
from numpy.testing import assert_allclose, assert_array_equal
from scipy import stats
from statsmodels.tsa.adfvalues import mackinnonp, mackinnoncrit
from statsmodels.regression.linear_model import OLSResults, RegressionResultsWrapper

from unitroot_tools import EstimationError, adf
from unitroot_tools._adapters import _adf_calibration
from unitroot_tools._regression import scaled_fit

V = Path(__file__).resolve().parents[1] / 'validation/scaled_svd'
INPUTS = json.loads((V/'adf_inputs.json').read_text())['inputs']
REFERENCES = json.loads((V/'exact_references.json').read_text())['fits']
CALIBRATION = json.loads((V/'adf_calibration_references.json').read_text())['comparisons']


# Frozen values cross numerical-library/platform boundaries. Reserve 128 binary64
# eps per CDF evaluation; two such evaluations differ by at most 2*b/(1-b)
# relative to the stored one, provided each meets that accuracy budget. External
# 80/120-digit Decimal CDF analysis checks this budget independently of the small
# Mac/Colab pairwise difference. This is a reviewed accuracy target, not a proof
# for every libm implementation. See docs/Numerical_Validation.md.
_CALIBRATION_EVAL_RTOL = 128 * np.finfo(np.float64).eps
_FROZEN_PVALUE_RTOL = 2 * _CALIBRATION_EVAL_RTOL / (1 - _CALIBRATION_EVAL_RTOL)


def _assert_frozen_pvalue(actual, expected):
    assert math.isfinite(actual) and 0 <= actual <= 1
    assert math.isfinite(expected) and 0 <= expected <= 1
    if expected in (0., 1.):
        assert actual == expected  # Exact endpoint/saturation contract.
    else:
        assert 0 < actual < 1
        # Divide the error rather than multiplying a tiny expected value by
        # rtol: no absolute floor, including for positive subnormal p-values.
        assert abs(actual - expected) / expected <= _FROZEN_PVALUE_RTOL
    # P-value classifications are a separate compatibility requirement. The
    # package's actual reject() decision uses critical values, tested elsewhere.
    for alpha in (.01, .05, .10):
        assert (actual < alpha) == (expected < alpha)


def hexarray(values):
    return np.array([[float.fromhex(v) for v in row] for row in values])


@pytest.mark.parametrize('trend',['n','c','ct','ctt'])
@pytest.mark.parametrize('p',[0,4])
def test_exact_design_response_names_sample_and_df(stage2,trend,p):
    y=stage2[0]['rw_negma']; saved=y.copy()
    r=adf(y,trend=trend,lags=p); raw=ADF(y,trend=trend,lags=p)
    assert_array_equal(r.regression.model.exog,raw.regression.model.exog)
    assert_array_equal(r.regression.model.endog,raw.regression.model.endog)
    assert r.regression.model.exog_names==raw.regression.model.exog_names
    assert_array_equal(y,saved)
    assert (r.input_nobs,r.nobs,r.sample_start,r.sample_stop)==(len(y),len(y)-p-1,p+1,len(y))
    assert r.regression.df_resid==raw.regression.df_resid
    assert r.regression.df_model==raw.regression.df_model
    assert list(r.critical_values)==[.01,.05,.10]
    assert r.critical_values[.01]<r.critical_values[.05]<r.critical_values[.10]
    fit=scaled_fit(r.regression.model.exog,r.regression.model.endog)
    assert_array_equal(r.regression.params,fit.params)
    assert_allclose(r.regression.ssr,fit.ssr,rtol=0,atol=0)
    assert_allclose(r.regression.cov_params(),fit.covariance,rtol=0,atol=0)
    assert_allclose(r.regression.tvalues,fit.tvalues,rtol=0,atol=0)
    assert r.statistic==r.regression.tvalues.iloc[0]
    assert r.notes==() and 'column-scaled SVD' in r.regression_description
    assert 'public Statsmodels' in r.pvalue_provenance
    assert 'public Statsmodels' in r.critical_values_provenance


@pytest.mark.parametrize('case',REFERENCES,ids=lambda c:c['name'])
def test_public_adf_exact_binary64_high_precision(case):
    scale,spec=case['name'].removeprefix('design_').removesuffix('.npz').split('_')
    key={'1.0':'y','1e-12':'positive','-1e-12':'negative'}[scale]
    y=np.array([float.fromhex(v) for v in INPUTS[key]])
    r=adf(y,trend=spec,lags=1); reg=r.regression
    assert_array_equal(reg.model.exog,hexarray(case['X']))
    assert_array_equal(reg.model.endog,[float.fromhex(v) for v in case['Y']])
    beta=np.array([float.fromhex(v) for v in case['beta']]);cov=hexarray(case['covariance']);se=np.sqrt(np.diag(cov))
    assert np.max(np.abs(np.asarray(reg.params)-beta)/se)<1e-12
    assert np.max(np.abs(reg.cov_params()-cov)/np.sqrt(np.outer(np.diag(cov),np.diag(cov))))<1e-12
    assert_allclose(reg.ssr,float(case['rss']),rtol=1e-12,atol=0)
    assert_allclose(reg.bse,se,rtol=1e-12,atol=0)
    assert_allclose(reg.tvalues,beta/se,rtol=0,atol=1e-12)
    assert_allclose(r.statistic,beta[0]/se[0],rtol=0,atol=1e-12)
    assert r.statistic==reg.tvalues.iloc[0]


@pytest.mark.parametrize('trend',['n','c','ct','ctt'])
def test_public_calibration_against_frozen_arch_evaluations(trend):
    pchecks=cvchecks=0
    for row in CALIBRATION:
        if row['spec']!=trend: continue
        if row['kind']=='pvalue':
            p,_=_adf_calibration(row['statistic'],trend,118)
            _assert_frozen_pvalue(p,row['arch']);pchecks+=1
        else:
            nobs=np.inf if row['nobs']=='infinity' else row['nobs']
            _,critical=_adf_calibration(-3.,trend,nobs)
            assert_array_equal(list(critical.values()),row['arch']);cvchecks+=1
    assert pchecks>=160 and cvchecks==10


@pytest.mark.parametrize('trend',['n','c','ct','ctt'])
def test_public_regression_object_compatibility(trend):
    r=adf(np.array([float.fromhex(v) for v in INPUTS['positive']]),trend=trend,lags=1)
    reg=r.regression
    assert isinstance(reg,RegressionResultsWrapper) and isinstance(reg._results,OLSResults)
    assert list(reg.params.index)==reg.model.exog_names
    assert_array_equal(reg.fittedvalues,reg.predict())
    assert_array_equal(reg.resid,reg.model.endog-reg.predict())
    assert_allclose(reg.cov_params(),reg.scale*reg.normalized_cov_params,rtol=0,atol=0)
    assert_allclose(reg.bse,np.sqrt(np.diag(reg.cov_params())),rtol=0,atol=0)
    assert_allclose(reg.pvalues,2*stats.t.sf(abs(reg.tvalues),reg.df_resid),rtol=1e-12)
    restriction=np.eye(len(reg.params))[0]
    assert_allclose(reg.t_test(restriction).tvalue,r.statistic,rtol=1e-12)
    assert_allclose(reg.f_test(restriction).fvalue,r.statistic**2,rtol=1e-12)
    assert np.isfinite(reg.llf) and np.isfinite(np.asarray(reg.conf_int())).all()


def test_no_unused_raw_adf_constructor(monkeypatch):
    import arch.unitroot as public
    def fail(*args,**kwargs): raise AssertionError('raw ADF must not be fitted')
    monkeypatch.setattr(public,'ADF',fail)
    r=adf([1.,4.,2.,8.,3.,6.],lags=0)
    assert np.isfinite(r.statistic)


@pytest.mark.parametrize('scale',[1e-250,1e250])
def test_public_adf_squared_unit_range_failures(scale):
    y=np.array([float.fromhex(v) for v in INPUTS['y']])*scale
    with pytest.raises(EstimationError,match='range|scale|variance'):
        adf(y,trend='ct',lags=1)


@pytest.mark.parametrize('trend',['n','c','ct','ctt'])
def test_calibration_exact_current_runtime_upstream_compatibility(trend):
    # ARCH/Statsmodels share surfaces/formulas: compatibility, not independent
    # mathematical verification. Frozen and high-precision tests remain separate.
    for row in CALIBRATION:
        if row['spec'] != trend:
            continue
        if row['kind'] == 'pvalue':
            p,_ = _adf_calibration(row['statistic'],trend,118)
            assert p == mackinnonp(row['statistic'],regression=trend,N=1)
            assert p == arch_mackinnonp(row['statistic'],regression=trend)
        else:
            nobs = np.inf if row['nobs'] == 'infinity' else row['nobs']
            _,critical = _adf_calibration(-3.,trend,nobs)
            assert_array_equal(list(critical.values()),mackinnoncrit(N=1,regression=trend,nobs=nobs))
            assert_array_equal(list(critical.values()),arch_mackinnoncrit(regression=trend,nobs=nobs))


@pytest.mark.parametrize('expected',[2.0221237233452862e-30,1e-300,1e-200,1e-30,.2,.9])
@pytest.mark.parametrize('direction',[0.,1.])
def test_frozen_pvalue_accepts_last_bit_platform_differences(expected,direction):
    actual = expected
    for _ in range(8):
        actual = float(np.nextafter(actual,direction))
    _assert_frozen_pvalue(actual,expected)


@pytest.mark.parametrize('expected',[2.0221237233452862e-30,1e-300,1e-100,.02,.2,.8])
@pytest.mark.parametrize('relative_error',[-1e-10,1e-10])
def test_frozen_pvalue_rejects_material_relative_calibration_error(expected,relative_error):
    with pytest.raises(AssertionError):
        _assert_frozen_pvalue(expected*(1+relative_error),expected)


@pytest.mark.parametrize('expected',[1e-300,.2])
def test_frozen_pvalue_enforces_reviewed_relative_budget(expected):
    _assert_frozen_pvalue(expected*(1+_FROZEN_PVALUE_RTOL/2),expected)
    with pytest.raises(AssertionError):
        _assert_frozen_pvalue(expected*(1+2*_FROZEN_PVALUE_RTOL),expected)


@pytest.mark.parametrize('endpoint',[0.,1.])
def test_frozen_pvalue_preserves_exact_endpoints(endpoint):
    _assert_frozen_pvalue(endpoint,endpoint)
    neighbor = float(np.nextafter(endpoint,1-endpoint))
    with pytest.raises(AssertionError):
        _assert_frozen_pvalue(neighbor,endpoint)
    with pytest.raises(AssertionError):
        _assert_frozen_pvalue(endpoint,neighbor)


@pytest.mark.parametrize('invalid',[np.nan,np.inf,-np.inf,-.1,1.1])
def test_frozen_pvalue_rejects_invalid_probabilities(invalid):
    with pytest.raises(AssertionError):
        _assert_frozen_pvalue(invalid,.2)
    with pytest.raises(AssertionError):
        _assert_frozen_pvalue(.2,invalid)


@pytest.mark.parametrize('alpha',[.01,.05,.10])
@pytest.mark.parametrize('below_reference',[False,True])
def test_frozen_pvalue_rejects_threshold_crossings_despite_numeric_closeness(alpha,below_reference):
    below = float(np.nextafter(alpha,0.))
    actual,expected = (alpha,below) if below_reference else (below,alpha)
    assert abs(actual-expected)/expected < _FROZEN_PVALUE_RTOL
    with pytest.raises(AssertionError):
        _assert_frozen_pvalue(actual,expected)


@pytest.mark.parametrize('alpha',[.01,.05,.10])
@pytest.mark.parametrize('direction',[0.,1.])
def test_frozen_pvalue_accepts_same_side_threshold_neighbors(alpha,direction):
    expected = float(np.nextafter(alpha,direction))
    actual = float(np.nextafter(expected,direction))
    _assert_frozen_pvalue(actual,expected)


@pytest.mark.parametrize('expected',[float(np.nextafter(0.,1.)),1e-320])
def test_positive_subnormal_pvalues_have_no_vacuous_absolute_floor(expected):
    _assert_frozen_pvalue(expected,expected)
    with pytest.raises(AssertionError):
        _assert_frozen_pvalue(0.,expected)
    with pytest.raises(AssertionError):
        _assert_frozen_pvalue(float(np.nextafter(expected,1.)),expected)
