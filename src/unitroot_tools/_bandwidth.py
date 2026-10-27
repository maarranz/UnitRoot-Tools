"""Covariance-bandwidth adapters, separate from augmentation-lag selection."""
from __future__ import annotations

import math
import re
import warnings

import numpy as np
from numpy.typing import ArrayLike

from ._adapters import _backend, _check_levels, _check_regression, _NUMERICAL
from ._validation import label, optional_order, series
from ._version import __version__
from .errors import EstimationError, SpecificationWarning
from .results import UnitRootResult

_KPSS_DEFAULT_WARNING = (
    "Lag selection has changed to use a data-dependent method. To use the "
    "old method that only depends on time, set lags=-1"
)


def _estimate_bandwidth(y, *, test, trend, bandwidth, test_type=None) -> UnitRootResult:
    is_pp = test == "PP"
    trend = label(trend, "trend", ("n", "c", "ct") if is_pp else ("c", "ct"))
    if is_pp:
        test_type = label(test_type, "test_type", ("tau", "rho"))
    requested = optional_order(bandwidth, "bandwidth")
    data, name = series(y)
    N = len(data)
    d = {"n": 0, "c": 1, "ct": 2}[trend]
    minimum = d + (3 if is_pp else 2)
    if N < minimum:
        raise ValueError(f"{test} trend={trend}: at least {minimum} observations required")
    m = N-1 if is_pp else N
    if requested is not None and requested >= m:
        raise ValueError(f"{test}: bandwidth must be less than covariance sample size {m}")
    context = f"{test} trend={trend}"
    _check_levels(data, trend, context)
    notes = ()
    rule = None
    effective = requested
    if requested is not None:
        bandwidth_method = "fixed"
    elif is_pp:
        bandwidth_method = "schwert_capped"
        rule = math.ceil(12 * (N/100)**.25)
        effective = min(rule, N-2)
        if effective < rule:
            note = f"PP: default bandwidth capped from {rule} to {effective} for N={N}"
            warnings.warn(note, SpecificationWarning, stacklevel=3)
            notes = (note,)
    else:
        bandwidth_method = "hobijn"
    classes, backend_version = _backend()
    from statsmodels.tools.sm_exceptions import SingularMatrixWarning
    try:
        with warnings.catch_warnings(), np.errstate(over="raise", divide="raise", invalid="raise"):
            warnings.simplefilter("error", RuntimeWarning)
            warnings.simplefilter("error", SingularMatrixWarning)
            if not is_pp and requested is None:
                warnings.filterwarnings("ignore", message="^" + re.escape(_KPSS_DEFAULT_WARNING) + "$",
                                        category=DeprecationWarning)
            kwargs = dict(trend=trend, lags=effective)
            if is_pp:
                kwargs["test_type"] = test_type
            backend = classes[test](data, **kwargs)
            stat, pvalue = float(backend.stat), float(backend.pvalue)
            used = backend.lags
            nobs = backend.nobs
            critical = {a: float(backend.critical_values[k]) for a, k in
                        ((.01, "1%"), (.05, "5%"), (.10, "10%"))}
            regression = backend.regression if is_pp else None
            if is_pp:
                notes += _check_regression(regression, m, d+1, context)
    except EstimationError:
        raise
    except (_NUMERICAL + (SingularMatrixWarning,)) as exc:
        raise EstimationError(f"{context}: released ARCH estimation failed; check data scale and bandwidth") from exc
    if isinstance(used, (bool, np.bool_)) or not isinstance(used, (int, np.integer)) or not 0 <= used < m:
        raise EstimationError(f"{context}: invalid backend bandwidth {used}; specify a fixed bandwidth in [0, {m-1}]")
    if nobs != m or (effective is not None and used != effective):
        raise EstimationError(f"{context}: backend sample/bandwidth does not match the requested model")
    if is_pp and (regression.nobs != m or backend.test_type != test_type):
        raise EstimationError(f"{context}: unexpected backend regression sample or PP statistic type")
    if not np.isfinite((stat, pvalue, *critical.values())).all() or not 0 <= pvalue <= 1:
        raise EstimationError(f"{context}: backend inference is non-finite or out of range")
    if not is_pp and (stat <= 0 or pvalue < .0001):
        raise EstimationError(f"{context}: nonpositive statistic/long-run variance or invalid table probability")
    boundary = None
    lower = 0. if is_pp else .0001
    if pvalue in (lower, 1.):
        boundary = "lower_endpoint" if pvalue == lower else "upper_endpoint"
        source = "response-surface" if is_pp else "simulated-table"
        notes += (f"Approximate {test} p-value is at the {source} {boundary.replace('_', ' ')} ({pvalue:g}); not an exact tail probability.",)
    if is_pp:
        dist = "ADF-t" if test_type == "tau" else "ADF-z"
        calibration = f"ARCH 8.0.0 PP {test_type} ({dist}) calibration"
        p_provenance = f"ARCH 8.0.0 {dist} p-value response-surface approximation"
        cv_provenance = f"ARCH 8.0.0 sample-size {dist} critical-value response surfaces"
        description = f"OLS level AR(1) regression with HAC/Newey-West covariance (bandwidth {used})"
        null = "The stochastic component has a unit root."
        alternative = "The stochastic component is stationary around the stated deterministic terms."
    else:
        calibration = "ARCH 8.0.0 KPSS simulated-quantile interpolation"
        p_provenance = "ARCH 8.0.0 simulated KPSS quantile table; linear interpolation with endpoint saturation"
        cv_provenance = "ARCH 8.0.0 simulated KPSS quantile table (100 million replications, series length 2000)"
        description = "None: ARCH KPSS exposes no public regression result"
        null = "The series is stationary around " + ("an intercept." if trend == "c" else "a linear time trend.")
        alternative = "The series is not stationary around the stated deterministic terms."
    import scipy
    import statsmodels
    return UnitRootResult(
        test=test, statistic=stat, pvalue=pvalue, critical_values=critical, trend=trend,
        input_nobs=N, nobs=int(nobs), sample_start=1 if is_pp else 0, sample_stop=N,
        lags=None, lag_selection=None, regression=regression, regression_description=description,
        calibration=calibration, pvalue_provenance=p_provenance, critical_values_provenance=cv_provenance,
        versions={"unitroot_tools": __version__, "arch": backend_version, "numpy": np.__version__,
                  "scipy": scipy.__version__, "statsmodels": statsmodels.__version__},
        null_hypothesis=null, alternative_hypothesis=alternative, notes=notes,
        input_name=name, backend_version=backend_version, tail="left" if is_pp else "right",
        pvalue_boundary=boundary, bandwidth=int(used), bandwidth_method=bandwidth_method,
        test_type=test_type, bandwidth_requested=requested, bandwidth_uncapped=rule,
    )


def pp(
    y: ArrayLike, *, trend: str = "c", bandwidth: int | None = None, test_type: str = "tau"
) -> UnitRootResult:
    """Phillips-Perron via public ARCH 8.0.0, with a covariance bandwidth.

    Trends n/c/ct and statistic types tau/rho are supported. None uses
    min(ceil(12*(N/100)**.25), N-2), warning if capped. Explicit b, including
    zero, requires 0 <= b < N-1. This is never augmentation-lag selection.
    The public regression uses HAC covariance; inference is PP calibration.
    """
    return _estimate_bandwidth(y, test="PP", trend=trend, bandwidth=bandwidth, test_type=test_type)


def kpss(
    y: ArrayLike, *, trend: str = "c", bandwidth: int | None = None
) -> UnitRootResult:
    """KPSS stationarity test via public ARCH 8.0.0 (right-tail rejection).

    Trends c/ct only. None delegates Hobijn-style automatic bandwidth to
    public ARCH, rejecting any resolved b outside 0 <= b < N. Fixed zero is
    supported. The calibration is ARCH's simulated table, with boundary notes.
    No public regression result is available; regression is None.
    """
    return _estimate_bandwidth(y, test="KPSS", trend=trend, bandwidth=bandwidth)
