"""Sequential Enders rules; inference conventions are in docs/Enders.md."""
from dataclasses import dataclass
from numbers import Real
from types import MappingProxyType

import numpy as np
from scipy import linalg, stats

from ._adapters import adf, _check_regression
from ._validation import label, optional_order, series
from .errors import EstimationError

# Enders, Supplementary Manual (4th ed.), Table B, PDF page 123.
# Columns: 1%, 5%, 10%. Transcribed from the primary table, not urca code.
_PHI = {
    'phi1': ((7.88, 5.18, 4.12), (7.06, 4.86, 3.94), (6.70, 4.71, 3.86),
             (6.52, 4.63, 3.81), (6.47, 4.61, 3.79)),
    'phi3': ((10.61, 7.24, 5.91), (9.31, 6.73, 5.61), (8.73, 6.49, 5.47),
             (8.43, 6.34, 5.39), (8.34, 6.30, 5.36)),
}
_GRID = (25, 50, 100, 250, 500)


def _phi_critical(test, nobs, alpha):
    if nobs < 25:
        raise ValueError('Enders requires at least 25 effective observations for phi tables')
    row = max(i for i, n in enumerate(_GRID) if n <= nobs)
    return _PHI[test][row][(.01, .05, .10).index(alpha)], _GRID[row]


@dataclass(frozen=True, slots=True)
class EndersStep:
    """One actual test; policy termination is recorded as its next_action.

    Regression objects are exposed for inspection, not deeply immutable.
    Signed t statistics use tail='two-sided' for deterministic coefficients.
    """
    specification: str
    test: str
    null_hypothesis: str
    statistic: float
    critical_value: float
    reject: bool
    tail: str
    inference: str
    next_action: str
    regression: object
    restricted_ssr: float | None = None
    restricted_columns: tuple[str, ...] | None = None
    restrictions: int | None = None
    table_sample_size: int | None = None
    reference_df: float | None = None


@dataclass(frozen=True, slots=True)
class EndersResult:
    variant: str
    alpha: float
    allow_no_constant: bool
    allow_no_constant_requested: bool
    trend: str
    lags: int
    lag_method: str
    lag_selection: object
    reject_unit_root: bool
    terminal_test: str
    statistic: float
    critical_value: float
    history: tuple[EndersStep, ...]
    input_nobs: int
    nobs: int
    sample_start: int
    sample_stop: int
    versions: object
    notes: tuple[str, ...]

    def __post_init__(self):
        object.__setattr__(self, 'versions', MappingProxyType(dict(self.versions)))

    def summary(self):
        title = 'Original Enders' if self.variant == 'original' else 'Conservative Enders (proposed modification)'
        lines = [title, f'Nominal alpha: {self.alpha:g}; common lag order: {self.lags} ({self.lag_method})',
                 f'No-constant selection permitted: {self.allow_no_constant}',
                 f'Selected specification: {self.trend}; observations: {self.nobs}; response slice [{self.sample_start}:{self.sample_stop}]']
        for i, h in enumerate(self.history, 1):
            lines.append(f'{i}. {h.specification} {h.test}: H0 {h.null_hypothesis}; statistic={h.statistic:.6g}, '
                         f'critical={h.critical_value:.6g} ({h.tail}); '
                         f'{"reject" if h.reject else "do not reject"}; {h.next_action}')
        lines += [f'Final unit-root decision: {"reject" if self.reject_unit_root else "do not reject"} H0 using {self.terminal_test}.',
                  'Nonrejection does not establish a unit root. Joint and deterministic tests are not unit-root decisions.',
                  'Sequential pretesting does not guarantee overall nominal size.']
        lines.extend('Note: '+n for n in self.notes)
        return '\n'.join(lines)


def _joint(result):
    reg = result.regression
    names = list(reg.model.exog_names)
    removed = ('Level.L1', 'trend') if result.trend == 'ct' else ('Level.L1', 'const')
    # ARCH names are checked rather than silently relying on column positions.
    if any(name not in names for name in removed):
        raise EstimationError('Unexpected ARCH ADF regressor names')
    keep = [i for i, name in enumerate(names) if name not in removed]
    x = np.asarray(reg.model.exog)[:, keep]
    y = np.asarray(reg.model.endog)
    residual = y - x @ linalg.lstsq(x, y)[0] if keep else y
    rss = float(residual @ residual)
    value = (rss-float(reg.ssr))/2 / (float(reg.ssr)/reg.df_resid)
    if not np.isfinite(value) or value < -1e-8:
        raise EstimationError('Invalid nested-regression joint statistic')
    return max(0., float(value)), rss, tuple(names[i] for i in keep)


def _auxiliary(result):
    reg = result.regression
    names = list(reg.model.exog_names)
    keep = [i for i, name in enumerate(names) if name != 'Level.L1']
    if len(keep) != len(names)-1:
        raise EstimationError('Unexpected ARCH ADF level regressor')
    from statsmodels.regression.linear_model import OLS
    from pandas import DataFrame
    design = DataFrame(np.asarray(reg.model.exog)[:, keep], columns=[names[i] for i in keep])
    fit = OLS(np.asarray(reg.model.endog), design).fit()
    _check_regression(fit, result.nobs, len(keep), 'Enders auxiliary regression')
    term = 'trend' if result.trend == 'ct' else 'const'
    j = [names[i] for i in keep].index(term)
    return fit, float(np.asarray(fit.tvalues)[j])


def enders(y, *, variant='original', alpha=.05, lags=None, max_lags=None,
           method='maic', allow_no_constant=False):
    """Sequential Enders procedure, with a complete chronological audit trail.

    See docs/Enders.md for primary references, conditional Student-t inference,
    finite phi-table lookup, and the common-lag convention. MAIC selects once
    under ct. For fixed lags, default 'maic' (or None) is inactive; other methods
    or max_lags conflict. Original ignores allow_no_constant and always permits
    n. Conservative defaults to retaining c on tau2 nonrejection.
    """
    variant = label(variant, 'variant', ('original', 'conservative'))
    if isinstance(alpha, (bool, np.bool_)) or not isinstance(alpha, Real):
        raise TypeError('alpha must be 0.01, 0.05 or 0.10')
    alpha = float(alpha)
    if alpha not in (.01, .05, .10):
        raise ValueError('alpha must be 0.01, 0.05 or 0.10')
    if not isinstance(allow_no_constant, (bool, np.bool_)):
        raise TypeError('allow_no_constant must be boolean')
    allowed = variant == 'original' or bool(allow_no_constant)
    p = optional_order(lags, 'lags')
    max_lags = optional_order(max_lags, 'max_lags')
    method = 'maic' if method is None else label(method, 'method', ('maic', 'mbic', 'aic', 'bic'))
    if p is not None and (max_lags is not None or method != 'maic'):
        raise ValueError('Fixed lags cannot be combined with max_lags or a nondefault method')
    data, _ = series(y)
    if len(data)-1-(p or 0) < 25:
        raise ValueError('Enders requires at least 25 effective observations')
    initial = adf(data, trend='ct', lags=p) if p is not None else adf(data, trend='ct', max_lags=max_lags, method=method)
    if initial.nobs < 25:
        raise ValueError('Selected lag order leaves fewer than 25 effective observations; specify a smaller max_lags')
    p = initial.lags
    history = []
    notes = [*initial.notes, 'Lag order selected once under ct and retained; no per-stage reselection.',
             'Phi tables use the largest finite tabulated T <= effective n (500 for n >= 500); no interpolation.',
             'ADF tau critical values retain ARCH calibration; phi uses Enders Table B.',
             'Conditional Student-t inference follows the published sequence, not a size guarantee.']
    if variant == 'original':
        notes.append('allow_no_constant is ignored for original; no-constant selection is always permitted.')
    else:
        notes.append('Conservative Enders is a proposed modification, not an independently validated procedure.')

    def finish(result):
        last = history[-1]
        return EndersResult(variant, alpha, allowed, bool(allow_no_constant), result.trend, p,
                            initial.lag_selection.method if initial.lag_selection else 'fixed',
                            initial.lag_selection, last.reject, last.test, last.statistic, last.critical_value,
                            tuple(history), len(data), result.nobs, p+1, len(data), initial.versions, tuple(notes))

    for trend, tau, phi, term, lower in [('ct', 'tau3', 'phi3', 'trend', 'c'), ('c', 'tau2', 'phi1', 'intercept', 'n'), ('n', 'tau1', None, None, None)]:
        result = initial if trend == 'ct' else adf(data, trend=trend, lags=p)
        if (result.nobs, result.sample_start, result.sample_stop) != (initial.nobs, initial.sample_start, initial.sample_stop):
            raise EstimationError('Sequential regression samples do not match')
        rejected = result.reject(alpha)
        policy_stop = trend == 'c' and not allowed and not rejected
        stop = rejected or trend == 'n' or policy_stop
        action = ('terminate: retain constant by policy with tau2 nonrejection' if policy_stop else
                  'terminate: unit-root decision' if stop else 'evaluate '+phi)
        history.append(EndersStep(trend, tau, 'b=0 against b<0', result.statistic,
                       result.critical_values[alpha], rejected, 'left', result.calibration, action, result.regression))
        if stop:
            return finish(result)
        value, rss, columns = _joint(result)
        critical, table_n = _phi_critical(phi, result.nobs, alpha)
        rejected = value > critical
        history.append(EndersStep(trend, phi, 'gamma=b=0' if trend == 'ct' else 'a=b=0',
                       value, critical, rejected, 'right', 'Enders Table B nonstandard joint distribution',
                       'estimate auxiliary '+term if rejected else 'proceed to '+lower,
                       result.regression, rss, columns, 2, table_n))
        if not rejected:
            continue
        aux, value = _auxiliary(result)
        critical = float(stats.t.ppf(1-alpha/2, aux.df_resid))
        rejected = abs(value) > critical
        history.append(EndersStep(trend, 'auxiliary_'+term, 'gamma=0' if trend == 'ct' else 'a=0',
                       value, critical, rejected, 'two-sided', 'Student t in b=0 auxiliary regression',
                       'conditional unit-root test' if rejected else 'proceed to '+lower,
                       aux, reference_df=float(aux.df_resid)))
        if not rejected:
            continue
        critical = float(stats.t.ppf(alpha, result.regression.df_resid))
        history.append(EndersStep(trend, 'conditional_t', 'b=0 against b<0', result.statistic,
                       critical, result.statistic < critical, 'left',
                       'Conditional Student t in unrestricted ADF regression', 'terminate: conditional unit-root decision',
                       result.regression, reference_df=float(result.regression.df_resid)))
        return finish(result)
    raise AssertionError('Unreachable Enders branch')
