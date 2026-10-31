"""Private full-rank OLS with column equilibration and original-unit results."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import linalg

from .errors import EstimationError

_EPS = np.finfo(np.float64).eps


@dataclass(frozen=True)
class ScaledFit:
    params: np.ndarray
    fittedvalues: np.ndarray
    resid: np.ndarray
    ssr: float
    scale: float
    normalized_cov_params: np.ndarray
    covariance: np.ndarray
    bse: np.ndarray
    tvalues: np.ndarray
    pinv: np.ndarray
    df_resid: int
    condition: float


def residual_ssr(residual, context="OLS") -> float:
    """Check squared-unit range before using the ordinary OLS residual dot product."""
    r = np.asarray(residual, dtype=float)
    if r.ndim != 1 or not np.isfinite(r).all():
        raise EstimationError(f"{context}: non-finite residuals")
    maximum = float(np.max(np.abs(r))) if r.size else 0.
    if maximum == 0:
        raise EstimationError(f"{context}: zero residual variance")
    try:
        with np.errstate(all="raise"):
            length = maximum * float(linalg.norm(r / maximum))
            safe = length * length
            rss = float(r @ r)
        if (not np.isfinite((safe, rss)).all() or safe <= 0 or rss <= 0
                or abs(rss / safe - 1) > _EPS * max(1, len(r))):
            raise FloatingPointError("squared residual units out of range")
    except (FloatingPointError, OverflowError) as exc:
        raise EstimationError(f"{context}: residual variance out of numeric range; rescale y") from exc
    return rss


def scaled_fit(x, y, context="OLS") -> ScaledFit:
    """Fit all columns, retaining the adapter's normalized-SVD rank threshold.

    X = Z D, y = c v, beta = c D^-1 V S^-1 U' v.
    The normalized covariance is D^-1 V S^-2 V' D^-1, never (X'X)^-1.
    No centering, column removal, reduced-rank fallback or sample change occurs.
    """
    if np.iscomplexobj(x) or np.iscomplexobj(y):
        raise EstimationError(f"{context}: real regression arrays required")
    try:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
    except (TypeError, ValueError) as exc:
        raise EstimationError(f"{context}: numeric regression arrays required") from exc
    if x.ndim != 2 or y.ndim != 1 or x.shape[0] != len(y):
        raise EstimationError(f"{context}: regression dimensions do not match")
    n, k = x.shape
    if k < 1 or n <= k:
        raise EstimationError(f"{context}: positive residual degrees of freedom required")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise EstimationError(f"{context}: non-finite regression arrays")
    maxima = np.max(np.abs(x), axis=0)
    if np.any(maxima == 0):
        raise EstimationError(f"{context}: zero regressor column (rank deficient)")
    response_scale = float(np.max(np.abs(y)))
    if response_scale == 0:
        raise EstimationError(f"{context}: zero response and residual variance")
    try:
        with np.errstate(all="raise"):
            z = x / maxima
            lengths = linalg.norm(z, axis=0)
            z /= lengths
            u, singular, vt = np.linalg.svd(z, full_matrices=False)
            if singular[-1] <= _EPS * max(x.shape) * singular[0]:
                raise EstimationError(f"{context}: rank-deficient final regression")
            right = vt.T / singular
            coefficient = right @ (u.T @ (y / response_scale))
            beta = (coefficient * response_scale) / lengths / maxima
            mapped_right = right / lengths[:, None] / maxima[:, None]
            normalized_covariance = mapped_right @ mapped_right.T
            pinv = mapped_right @ u.T
            fitted = x @ beta
            residual = y - fitted
            if linalg.norm(residual / response_scale) <= _EPS * max(x.shape) * linalg.norm(y / response_scale):
                raise EstimationError(f"{context}: zero residual variance (numerically deterministic)")
            rss = residual_ssr(residual, context)
            variance = rss / (n - k)
            covariance = variance * normalized_covariance
            se = np.sqrt(np.diag(covariance))
            tvalues = beta / se
        if (variance <= 0 or np.any(se <= 0)
                or not all(np.isfinite(a).all() for a in
                           (beta, fitted, residual, normalized_covariance, pinv, covariance, se, tvalues))):
            raise FloatingPointError("non-finite coefficients or covariance")
    except (ValueError, linalg.LinAlgError, FloatingPointError, OverflowError) as exc:
        raise EstimationError(f"{context}: scaled SVD fit out of numeric range; rescale y") from exc
    return ScaledFit(beta, fitted, residual, rss, variance, normalized_covariance,
                     covariance, se, tvalues, pinv, n - k,
                     float(singular[0] / singular[-1]))


def scaled_ols(y, x, context="OLS"):
    """Construct a genuine Statsmodels OLS result on the original named design.

    OLSResults is initialized from a single coherent fit, rather than editing an
    already-fitted raw regression or replacing it with a look-alike object.
    """
    from statsmodels.regression.linear_model import OLS, OLSResults, RegressionResultsWrapper

    model = OLS(y, x)
    fit = scaled_fit(model.exog, model.endog, context)
    model.pinv_wexog = fit.pinv
    model.normalized_cov_params = fit.normalized_cov_params
    model.rank = model.exog.shape[1]
    # Public condition_number/eigenvals describe the raw model, not Z.
    model.wexog_singular_values = np.linalg.svd(model.exog, compute_uv=False)
    result = OLSResults(model, fit.params, normalized_cov_params=fit.normalized_cov_params)
    return RegressionResultsWrapper(result)
