"""Independent common-sample OLS lag selection; no backend selectors."""
from __future__ import annotations

import math
import warnings

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import linalg

from ._validation import label, optional_order, series
from .errors import EstimationError, SpecificationWarning
from .results import CandidateDiagnostics, LagSelectionResult

_EPS = np.finfo(np.float64).eps


def _fit(
    x: NDArray[np.float64], y: NDArray[np.float64], context: str
) -> tuple[NDArray[np.float64], NDArray[np.float64], float]:
    """Column-scaled SVD solve with explicit residuals and relative rank test."""
    norms = linalg.norm(x, axis=0)
    if np.any(norms == 0):
        raise EstimationError(f"{context}: zero regressor column (rank deficient)")
    scaled = x / norms
    cutoff = _EPS * max(scaled.shape)
    try:
        coefficient, _, rank, singular = linalg.lstsq(
            scaled, y, cond=cutoff, lapack_driver="gelsd", check_finite=False
        )
    except (ValueError, linalg.LinAlgError) as exc:
        raise EstimationError(f"{context}: SVD least squares failed") from exc
    if rank != x.shape[1]:
        raise EstimationError(f"{context}: rank-deficient candidate design")
    beta = coefficient / norms
    residual = y - scaled @ coefficient
    return beta, residual, float(singular[0] / singular[-1])


def _detrend(y: NDArray[np.float64], trend: str) -> tuple[NDArray[np.float64], float]:
    scale = float(np.max(np.abs(y)))
    if scale == 0 or np.all(y == y[0]):
        raise EstimationError(f"selection trend={trend}: constant series")
    v = y / scale
    if trend != "n":
        d = np.ones((len(y), 1))
        if trend == "ct":
            d = np.column_stack((d, np.arange(1, len(y) + 1) / len(y)))
        _, residual, _ = _fit(d, v, f"selection trend={trend}: OLS detrending")
        if linalg.norm(residual) <= _EPS * max(d.shape) * linalg.norm(v):
            raise EstimationError(f"selection trend={trend}: complete deterministic fit")
        v = residual
    residual_scale = float(np.max(np.abs(v)))
    return v / residual_scale, scale * residual_scale


def _restore(square: float, scale: float, context: str) -> float:
    """Use logs to avoid overflow/underflow in a premature scale**2."""
    try:
        value = math.exp(math.log(square) + 2 * math.log(scale))
    except (ValueError, OverflowError) as exc:
        raise EstimationError(f"{context}: squared-unit diagnostics out of range; rescale y") from exc
    if not math.isfinite(value) or value <= 0:
        raise EstimationError(f"{context}: squared-unit diagnostics out of range; rescale y")
    return value


def _choose(candidates: tuple[CandidateDiagnostics, ...]) -> tuple[CandidateDiagnostics, float | None]:
    ordered = sorted(candidates, key=lambda r: (r.criterion_value, r.lags))
    gap = ordered[1].criterion_value - ordered[0].criterion_value if len(ordered) > 1 else None
    return ordered[0], gap


def select_lags(
    y: ArrayLike, *, trend: str = "c", max_lags: int | None = None, method: str = "maic"
) -> LagSelectionResult:
    """Select 0..K on one OLS-detrended sample using MAIC/MBIC/AIC/BIC.

    Parameters
    ----------
    y : one-dimensional real numeric input
        Equally spaced observations. Missing values are rejected, never dropped.
    trend : {"n", "c", "ct"}
        Deterministics removed by OLS from the complete series before trimming.
    max_lags : int or None
        Common-sample maximum K. None uses ceil(12*(N/100)**.25), capped at
        floor((N-3)/2) with a SpecificationWarning; explicit invalid K is an error.
    method : {"maic", "mbic", "aic", "bic"}
        MAIC is default. Conventional scores omit tau; all four scores are
        reported in every candidate. Labels are case insensitive.

    Notes
    -----
    Candidate response indices are K+1..N-1 in zero-based level indexing,
    so n=N-1-K for every candidate. Variance is SSR/n, never SSR/df_resid.
    Sum squared lagged levels over all those rows. Scores use natural logs,
    per observation, omitting Gaussian/common parameter-count constants.
    Exact ties select the smallest lag; no near-tie tolerance is used.
    Final unit-root tests and their calibrations are outside this function.
    """
    trend = label(trend, "trend", ("n", "c", "ct"))
    method = label(method, "method", ("maic", "mbic", "aic", "bic"))
    requested = optional_order(max_lags, "max_lags")
    data, name = series(y)
    N = len(data)
    d = {"n": 0, "c": 1, "ct": 2}[trend]
    if N < max(3, d + 2):
        raise ValueError(f"selection trend={trend}: at least {max(3, d + 2)} observations required")
    bound = (N - 3) // 2
    rule = math.ceil(12 * (N / 100) ** .25) if requested is None else None
    K = min(rule, bound) if requested is None else requested
    if K > bound:
        raise ValueError(f"max_lags={K} is infeasible for N={N}; maximum is {bound}")
    notes: list[str] = []
    if rule is not None and K < rule:
        notes.append(f"Default max_lags capped from {rule} to {K} for N={N}")
        warnings.warn(notes[-1], SpecificationWarning, stacklevel=2)
    v, scale = _detrend(data, trend)
    dy = np.diff(v)
    lhs = dy[K:]
    level = v[K:-1]
    n = len(lhs)
    level_ss = float(level @ level)
    if level_ss == 0:
        raise EstimationError(f"selection trend={trend} K={K}: zero matched level sum (rank deficient)")
    restored_ss = _restore(level_ss, scale, f"selection trend={trend} K={K}")
    log_scale = 2 * math.log(scale)
    rows = []
    for p in range(K + 1):
        context = f"selection trend={trend} K={K} p={p}"
        x = np.column_stack([level] + [dy[K-j:len(dy)-j] for j in range(1, p + 1)])
        beta, residual, condition = _fit(x, lhs, context)
        if linalg.norm(residual) <= _EPS * max(x.shape) * linalg.norm(lhs):
            raise EstimationError(f"{context}: zero residual variance (numerically deterministic)")
        ssr = float(residual @ residual)
        variance = ssr / n
        tau = float(beta[0]**2 * level_ss / variance)
        log_var = math.log(variance) + log_scale
        aic = log_var + 2 * p / n
        bic = log_var + math.log(n) * p / n
        maic = aic + 2 * tau / n
        mbic = bic + math.log(n) * tau / n
        if not np.isfinite((beta[0], tau, aic, bic, maic, mbic)).all():
            raise EstimationError(f"{context}: non-finite coefficient or criterion; rescale y")
        if condition > 1 / math.sqrt(_EPS):
            notes.append(f"{context}: poorly conditioned full-rank design (condition={condition:.3g})")
            warnings.warn(notes[-1], SpecificationWarning, stacklevel=2)
        score = {"aic": aic, "bic": bic, "maic": maic, "mbic": mbic}[method]
        rows.append(CandidateDiagnostics(
            lags=p, nobs=n, df_resid=n-p-1, level_coefficient=float(beta[0]),
            ssr=_restore(ssr, scale, context), variance=_restore(variance, scale, context),
            level_sum_squares=restored_ss, tau=tau, aic=aic, bic=bic,
            maic=maic, mbic=mbic, criterion_value=score,
        ))
    candidates = tuple(rows)
    winner, gap = _choose(candidates)
    return LagSelectionResult(
        method=method, selection_detrending="none" if trend == "n" else "ols", trend=trend,
        input_nobs=N, max_lags_requested=requested, max_lags=K, uncapped_rule_value=rule,
        selected_lags=winner.lags, selection_nobs=n, sample_start=K+1, sample_stop=N,
        criterion_value=winner.criterion_value, runner_up_gap=gap, candidates=candidates,
        notes=tuple(notes), input_name=name,
    )
