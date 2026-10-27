from dataclasses import FrozenInstanceError, replace
import numpy as np
import pytest

from unitroot_tools import adf, dfgls


@pytest.fixture
def result():
    return adf(np.random.default_rng(11).normal(size=60).cumsum(),lags=1)


@pytest.mark.parametrize('alpha',[.01,.05,.10,np.float64(.05)])
def test_strict_left_tail_at_critical_boundary(result,alpha):
    critical=result.critical_values[alpha]
    at=replace(result,statistic=critical,pvalue=0.)
    assert at.reject(alpha) is False
    assert replace(at,statistic=np.nextafter(critical,-np.inf),pvalue=1.).reject(alpha) is True
    assert replace(at,statistic=np.nextafter(critical,np.inf),pvalue=0.).reject(alpha) is False


@pytest.mark.parametrize('alpha',[0.,1.,.025,float('nan'),float('inf'),-.05])
def test_invalid_alpha_values(result,alpha):
    with pytest.raises(ValueError):result.reject(alpha)


@pytest.mark.parametrize('alpha',[True,False,np.bool_(True),'5%',None,[.05]])
def test_invalid_alpha_types(result,alpha):
    with pytest.raises(TypeError):result.reject(alpha)


def test_frozen_results_and_copied_maps(result):
    with pytest.raises(FrozenInstanceError):result.statistic=0.
    with pytest.raises(TypeError):result.critical_values[.05]=-99.
    with pytest.raises(TypeError):result.versions['arch']='dev'
    source=dict(result.critical_values);versions=dict(result.versions)
    copied=replace(result,critical_values=source,versions=versions,notes=['test'])
    source[.05]=-99.;versions['arch']='dev'
    assert copied.critical_values[.05]==result.critical_values[.05]
    assert copied.versions['arch']=='8.0.0' and copied.notes==('test',)
    # Regression is a documented public backend object, not deeply frozen.
    assert copied.regression is result.regression


def test_provenance_and_regression_labels(result):
    assert result.calibration=='ARCH 8.0.0 ADF-t calibration'
    assert 'MacKinnon' in result.pvalue_provenance
    assert 'no-constant exception' in result.critical_values_provenance
    assert result.pvalue_status=='approximate' and result.tail=='left'
    assert result.pvalue_boundary is None
    assert set(result.versions)=={'unitroot_tools','arch','numpy','scipy','statsmodels'}
    r=dfgls(np.random.default_rng(11).normal(size=60),lags=1,trend='ct')
    assert 'dedicated DFGLS' in r.calibration
    assert 'simulation' in r.pvalue_provenance and 'simulation' in r.critical_values_provenance
    assert 'GLS-detrended' in r.regression_description


def test_summary_teaches_failure_and_shows_samples(result):
    r=replace(result,statistic=0.,pvalue=1.)
    text=r.summary()
    assert 'fail to reject H0' in text and 'accept H0' not in text
    assert 'Failure to reject does not establish the null.' in text
    assert 'fixed order' in text and 'Final observations: 58 / 60' in text
    assert '1%:' in text and '5%:' in text and '10%:' in text
    assert 'approximate p-value' in text and 'H0:' in text
    r=replace(result,statistic=-100.,pvalue=1.)
    assert 'Decision at 5%: reject H0' in r.summary()


def test_automatic_summary_and_series_name():
    import pandas as pd
    y=pd.Series(np.random.default_rng(13).normal(size=80),name='inflation')
    r=dfgls(y,max_lags=5,method='bic')
    assert r.input_name=='inflation' and r.lag_selection.input_name=='inflation'
    text=r.summary()
    assert 'Series: inflation' in text and 'DF-GLS' in text
    assert 'BIC' in text and 'Selection observations: 74' in text
    assert 'full-series ols detrending' in text


def test_approximation_endpoints_are_not_claimed_exact(monkeypatch):
    from types import SimpleNamespace
    import arch.unitroot as public
    y=np.random.default_rng(9).normal(size=40)
    native=public.ADF(y,lags=0)
    for p,endpoint in [(0.,'lower_endpoint'),(1.,'upper_endpoint')]:
        wrapped=SimpleNamespace(stat=native.stat,pvalue=p,nobs=native.nobs,lags=0,
                                critical_values=native.critical_values,regression=native.regression)
        monkeypatch.setattr(public,'ADF',lambda *args,_wrapped=wrapped,**kwargs:_wrapped)
        r=adf(y,lags=0)
        assert r.pvalue_boundary==endpoint
        assert 'not an exact tail probability' in r.summary()
        assert r.reject()==bool(native.stat < native.critical_values['5%'])


def test_unsupported_tail_is_not_a_kpss_implementation(result):
    with pytest.raises(ValueError,match='left-tail'):replace(result,tail='right')


def test_noncanonical_float_alpha_is_value_error_not_key_error(result):
    with pytest.raises(ValueError):result.reject(np.float32(.05))
