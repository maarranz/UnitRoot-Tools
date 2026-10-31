"""Exact-array/high-precision checks for the private full-rank OLS layer."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy import stats
from statsmodels.regression.linear_model import OLSResults, RegressionResultsWrapper

from unitroot_tools import EstimationError
from unitroot_tools._regression import scaled_fit, scaled_ols

REFERENCE = json.loads((Path(__file__).resolve().parents[1] /
                       'validation/scaled_svd/exact_references.json').read_text())['fits']


def floats(a):
    return np.array([[float.fromhex(v) for v in row] for row in a])


@pytest.mark.parametrize('case', REFERENCE, ids=lambda c: c['name'])
def test_saved_binary64_high_precision_coefficients_and_covariance(case):
    x = floats(case['X']); y = np.array([float.fromhex(v) for v in case['Y']])
    beta = np.array([float.fromhex(v) for v in case['beta']]); cov = floats(case['covariance'])
    r = scaled_fit(x, y)
    se = np.sqrt(np.diag(cov))
    assert np.max(np.abs(r.params - beta) / se) < 1e-12
    assert np.max(np.abs(r.covariance - cov) / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))) < 1e-12
    assert_allclose(r.ssr, float(case['rss']), rtol=1e-12, atol=0)
    assert_allclose(r.bse, se, rtol=1e-12, atol=0)
    assert_allclose(r.tvalues, beta / se, rtol=0, atol=1e-12)
    assert r.df_resid == len(y) - x.shape[1]
    assert_array_equal(r.fittedvalues, x @ r.params)
    assert_array_equal(r.resid, y - r.fittedvalues)


def test_genuine_named_statsmodels_result_contract():
    t = np.arange(1., 81.)
    x = pd.DataFrame({'level': np.sin(t * .73), 'const': np.ones(80), 'trend': t})
    y = 2 * x['level'].to_numpy() + .2 + .01 * t + np.cos(t * .31)
    r = scaled_ols(y, x)
    assert isinstance(r, RegressionResultsWrapper) and isinstance(r._results, OLSResults)
    assert list(r.params.index) == list(x.columns)
    assert_array_equal(r.model.exog, x.to_numpy())
    assert_array_equal(r.model.endog, y)
    assert_array_equal(r.resid, y - r.predict(x))
    assert_allclose(r.cov_params(), r.scale * r.normalized_cov_params, rtol=0, atol=0)
    assert_allclose(r.bse, np.sqrt(np.diag(r.cov_params())), rtol=0, atol=0)
    assert_allclose(r.tvalues, r.params / r.bse, rtol=0, atol=0)
    assert_allclose(r.pvalues, 2 * stats.t.sf(abs(r.tvalues), r.df_resid), rtol=1e-12)
    assert_allclose(r.t_test([1, 0, 0]).tvalue, r.tvalues.iloc[0], rtol=1e-12)
    assert_allclose(r.f_test([1, 0, 0]).fvalue, r.tvalues.iloc[0] ** 2, rtol=1e-12)
    assert np.isfinite(r.llf) and np.isfinite(np.asarray(r.conf_int())).all()
    assert r.df_resid == 77 and r.nobs == 80


@pytest.mark.parametrize('x,y', [
    (np.ones(4), np.ones(4)), (np.ones((4, 2)), np.ones(3)),
    (np.ones((2, 2)), np.ones(2)), (np.c_[np.ones(10), np.zeros(10)], np.arange(10.)),
    (np.ones((10, 2)), np.arange(10.)), (np.ones((10, 1)), np.zeros(10)),
    (np.full((10, 1), np.nan), np.arange(10.)),
    (np.arange(1., 11.)[:, None], np.arange(1., 11.)),
])
def test_invalid_designs_and_zero_variance_fail_closed(x, y):
    with pytest.raises(EstimationError):
        scaled_fit(x, y)


@pytest.mark.parametrize('scale', [1e-250, 1e250])
def test_squared_unit_underflow_and_overflow_are_explicit(scale):
    t = np.arange(30.)
    with pytest.raises(EstimationError, match='range|variance'):
        scaled_fit(np.c_[np.sin(.8 * t), np.ones(30)], scale * np.cos(.33 * t))


def test_existing_normalized_rank_threshold_not_a_reduced_rank_fallback():
    t = np.linspace(-1, 1, 118)
    y = np.sin(3 * t) + .1 * np.cos(7 * t)
    r = scaled_fit(np.c_[np.ones(118), 1 + 1e-9 * t], y)
    assert r.condition > 1 / np.sqrt(np.finfo(float).eps) and r.df_resid == 116
    with pytest.raises(EstimationError, match='rank-deficient'):
        scaled_fit(np.c_[np.ones(118), 1 + 1e-14 * t], y)
