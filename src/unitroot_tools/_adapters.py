"""Public fixed-lag ARCH adapters; own selection always precedes estimation."""
from __future__ import annotations

from dataclasses import replace
import math
import warnings

import numpy as np
from numpy.typing import ArrayLike
from scipy import linalg

from ._validation import label, optional_order, series
from ._version import __version__
from .errors import EstimationError, SpecificationWarning
from .results import LagSelectionResult, UnitRootResult
from .selection import _fit, select_lags

_EPS = np.finfo(np.float64).eps
_NUMERICAL = (ValueError, RuntimeError, linalg.LinAlgError, FloatingPointError,
              OverflowError, ZeroDivisionError, RuntimeWarning)


def _backend():
    # Lazy imports preserve the independent selector's backend-free execution.
    import arch
    if arch.__version__ != "8.0.0":
        raise RuntimeError(f"Unsupported ARCH version {arch.__version__}; released 8.0.0 is required")
    from arch.unitroot import ADF, DFGLS, KPSS, PhillipsPerron
    return {"ADF": ADF, "DFGLS": DFGLS, "PP": PhillipsPerron, "KPSS": KPSS}, arch.__version__


def _check_levels(data: np.ndarray, trend: str, context: str) -> None:
    """Validate deterministic degeneracy without GLS-transforming backend data."""
    scale = float(np.max(np.abs(data)))
    if scale == 0 or np.all(data == data[0]):
        raise EstimationError(f"{context}: constant series")
    if trend != "n":
        t = np.arange(1, len(data)+1) / len(data)
        d = {"c": 1, "ct": 2, "ctt": 3}[trend]
        design = np.column_stack([t**j for j in range(d)])
        normalized = data / scale
        _, residual, _ = _fit(design, normalized, context + ": deterministic check")
        if linalg.norm(residual) <= _EPS * max(design.shape) * linalg.norm(normalized):
            raise EstimationError(f"{context}: complete deterministic fit")


def _selection(y, N, d, test, trend, requested, method) -> LagSelectionResult:
    bound = min((N-3)//2, (N-3-d)//2 if test == "ADF" else N-3-d)
    rule = math.ceil(12 * (N/100)**.25) if requested is None else None
    K = min(rule, bound) if requested is None else requested
    if K > bound:
        raise ValueError(f"{test} trend={trend}: max_lags={K} infeasible for N={N}; maximum is {bound}")
    cap_notes = ()
    if rule is not None and K < rule:
        note = f"{test} trend={trend}: default max_lags capped from {rule} to {K} for final-regression feasibility (N={N})"
        warnings.warn(note, SpecificationWarning, stacklevel=3)
        cap_notes = (note,)
    # Explicit effective K bypasses only the selector's default cap calculation;
    # candidate regressions, detrending, penalties and diagnostic values are unchanged.
    selected = select_lags(y, trend=trend, max_lags=K, method=method)
    return replace(selected, max_lags_requested=requested, uncapped_rule_value=rule,
                   notes=cap_notes + selected.notes)


def _check_regression(regression, nobs: int, columns: int, context: str) -> tuple[str, ...]:
    """Check the actual public backend design and residuals after eager fitting."""
    x = np.asarray(regression.model.exog, dtype=float)
    y = np.asarray(regression.model.endog, dtype=float).ravel()
    residual = np.asarray(regression.resid, dtype=float).ravel()
    if x.shape != (nobs, columns) or len(y) != nobs or len(residual) != nobs:
        raise EstimationError(f"{context}: unexpected backend regression dimensions")
    if regression.df_resid != nobs-columns:
        raise EstimationError(f"{context}: backend used a reduced-rank final regression; rescale y")
    if not all(np.isfinite(a).all() for a in (x, y, residual)):
        raise EstimationError(f"{context}: non-finite backend regression; rescale y")
    maxima = np.max(np.abs(x), axis=0)
    if np.any(maxima == 0):
        raise EstimationError(f"{context}: zero regressor column (rank deficient)")
    scaled = x / maxima
    scaled /= linalg.norm(scaled, axis=0)
    singular = linalg.svdvals(scaled)
    if singular[-1] <= _EPS * max(x.shape) * singular[0]:
        raise EstimationError(f"{context}: rank-deficient final regression")
    response_scale = float(np.max(np.abs(y)))
    if response_scale == 0:
        raise EstimationError(f"{context}: zero response and residual variance")
    if linalg.norm(residual/response_scale) <= _EPS * max(x.shape) * linalg.norm(y/response_scale):
        raise EstimationError(f"{context}: zero residual variance (numerically deterministic)")
    if not np.isfinite((regression.ssr, regression.scale)).all() or regression.ssr <= 0 or regression.scale <= 0:
        raise EstimationError(f"{context}: residual variance is out of numeric range; rescale y")
    condition = float(singular[0] / singular[-1])
    if condition > 1 / math.sqrt(_EPS):
        note = f"{context}: poorly conditioned full-rank final design (condition={condition:.3g})"
        warnings.warn(note, SpecificationWarning, stacklevel=3)
        return (note,)
    return ()


def _estimate(y, *, test, trend, lags, max_lags, method) -> UnitRootResult:
    trend = label(trend, "trend", ("n", "c", "ct", "ctt") if test == "ADF" else ("c", "ct"))
    p = optional_order(lags, "lags")
    requested = optional_order(max_lags, "max_lags")
    if p is not None and (method is not None or requested is not None):
        raise ValueError("Fixed lags cannot be combined with method or max_lags")
    chosen_method = label(method, "method", ("maic", "mbic", "aic", "bic")) if method is not None else "maic"
    if p is None and trend == "ctt":
        raise ValueError("ctt is supported for fixed-lag ADF only; specify lags")
    data, name = series(y)
    N = len(data)
    d = {"n": 0, "c": 1, "ct": 2, "ctt": 3}[trend]
    if N < d+3:
        raise ValueError(f"{test} trend={trend}: at least {d+3} observations required")
    if p is not None and (N < 3+d+p or N-1-p <= 1+p+(d if test == "ADF" else 0)):
        raise ValueError(f"{test} trend={trend}: lags={p} infeasible for N={N}; positive residual degrees of freedom required")
    context = f"{test} trend={trend}"
    _check_levels(data, trend, context)
    selected = None
    if p is None:
        selected = _selection(y, N, d, test, trend, requested, chosen_method)
        p = selected.selected_lags
    classes, backend_version = _backend()
    context += f" p={p}"
    expected_nobs = N-1-p
    # The original validated levels go to ARCH. No GLS or OLS selection residual
    # is substituted, and lags is always an integer, including zero.
    from statsmodels.tools.sm_exceptions import SingularMatrixWarning
    try:
        with warnings.catch_warnings(), np.errstate(over="raise", divide="raise", invalid="raise"):
            warnings.simplefilter("error", RuntimeWarning)
            warnings.simplefilter("error", SingularMatrixWarning)
            backend = classes[test](data, lags=p, trend=trend)
            stat, pvalue = float(backend.stat), float(backend.pvalue)
            nobs, used_lags = backend.nobs, backend.lags
            critical = {alpha: float(backend.critical_values[key]) for alpha, key in
                        ((.01, "1%"), (.05, "5%"), (.10, "10%"))}
            regression = backend.regression
            regression_notes = _check_regression(regression, expected_nobs, 1+p+(d if test == "ADF" else 0), context)
    except EstimationError:
        raise
    except (_NUMERICAL + (SingularMatrixWarning,)) as exc:
        raise EstimationError(f"{context}: released ARCH estimation failed; check data scale and specification") from exc
    if nobs != expected_nobs or used_lags != p or regression.nobs != expected_nobs:
        raise EstimationError(f"{context}: backend sample/order does not match the requested fixed-lag model")
    if not np.isfinite((stat, pvalue, *critical.values())).all() or not 0 <= pvalue <= 1:
        raise EstimationError(f"{context}: backend inference is non-finite or out of range")
    notes = (selected.notes if selected else ()) + regression_notes
    boundary = None
    if pvalue in (0., 1.):
        boundary = "lower_endpoint" if pvalue == 0 else "upper_endpoint"
        notes += (f"Approximate p-value is at the {boundary.replace('_', ' ')}; not an exact tail probability.",)
    if test == "ADF":
        calibration = "ARCH 8.0.0 ADF-t calibration"
        p_provenance = "MacKinnon asymptotic p-value response-surface approximation"
        cv_provenance = "ARCH 8.0.0 sample-size ADF-t critical-value surfaces; historical no-constant exception"
        description = "OLS ADF difference regression on original levels with requested deterministic terms"
    else:
        calibration = "ARCH 8.0.0 dedicated DFGLS calibration"
        p_provenance = "ARCH 8.0.0 dedicated DFGLS simulation p-value response surfaces"
        cv_provenance = "ARCH 8.0.0 dedicated DFGLS simulation critical-value response surfaces"
        description = "OLS no-deterministic difference regression on ARCH GLS-detrended levels"
    import scipy
    import statsmodels
    return UnitRootResult(
        test=test, statistic=stat, pvalue=pvalue, critical_values=critical, trend=trend,
        input_nobs=N, nobs=int(nobs), sample_start=p+1, sample_stop=N, lags=p,
        lag_selection=selected, regression=regression, regression_description=description,
        calibration=calibration, pvalue_provenance=p_provenance, critical_values_provenance=cv_provenance,
        versions={"unitroot_tools": __version__, "arch": backend_version, "numpy": np.__version__,
                  "scipy": scipy.__version__, "statsmodels": statsmodels.__version__},
        null_hypothesis="The stochastic component has a unit root.",
        alternative_hypothesis="The stochastic component is stationary around the stated deterministic terms.",
        notes=notes, input_name=name, backend_version=backend_version, pvalue_boundary=boundary,
    )


def adf(
    y: ArrayLike, *, trend: str = "c", lags: int | None = None,
    max_lags: int | None = None, method: str | None = None
) -> UnitRootResult:
    """ADF via released ARCH 8.0.0, at a fixed or independently selected lag.

    Fixed lags (including zero) prohibit method/max_lags. Otherwise default MAIC
    or explicit MBIC/AIC/BIC uses the approved common-sample OLS selector.
    Trends n/c/ct are supported; ctt requires fixed lags. Final inference is
    ARCH's, with strict left-tail critical-value decisions at 1%, 5% and 10%.
    """
    return _estimate(y, test="ADF", trend=trend, lags=lags, max_lags=max_lags, method=method)


def dfgls(
    y: ArrayLike, *, trend: str = "c", lags: int | None = None,
    max_lags: int | None = None, method: str | None = None
) -> UnitRootResult:
    """DF-GLS via released ARCH 8.0.0 with exactly one backend GLS detrending.

    Independent selection uses full-series OLS detrending (default MAIC;
    explicit MBIC/AIC/BIC). The original levels and an explicit integer lag are
    passed to ARCH. Trends c/ct only; fixed lags prohibit method/max_lags.
    """
    return _estimate(y, test="DFGLS", trend=trend, lags=lags, max_lags=max_lags, method=method)
