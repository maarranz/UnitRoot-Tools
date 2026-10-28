"""Reference evaluation using unchanged public APIs; no fit-cache substitution."""
import warnings
import numpy as np
from unitroot_tools import adf, enders, select_deterministic, EstimationError
from .catalogue import generate, innovations, array_hash

METHODS = ('adf_n', 'adf_c', 'adf_ct', 'enders_original', 'enders_conservative', 'det_aic', 'det_bic')
ALPHAS = (.01, .05, .10)
SPECS = ('n', 'c', 'ct')
EXPECTED_FAILURES = (EstimationError, TypeError, ValueError, FloatingPointError, np.linalg.LinAlgError)


class InvariantError(RuntimeError):
    pass


def enders_record(r):
    steps = []
    for h in r.history:
        steps.append(dict(specification=h.specification, test=h.test, null=h.null_hypothesis,
                          statistic=h.statistic, critical_value=h.critical_value, reject=h.reject,
                          tail=h.tail, inference=h.inference, next_action=h.next_action,
                          reference_df=h.reference_df, table_sample_size=h.table_sample_size,
                          restricted_ssr=h.restricted_ssr, restricted_columns=h.restricted_columns,
                          restrictions=h.restrictions))
    return dict(selected=r.trend, reject=r.reject_unit_root, terminal=r.terminal_test,
                statistic=r.statistic, critical_value=r.critical_value, lags=r.lags,
                nobs=r.nobs, sample_start=r.sample_start, sample_stop=r.sample_stop, steps=steps)


def evaluate(cell, replicate, phase='pilot', equivalence=False):
    z = innovations(cell, replicate, phase)
    y = generate(cell, z)
    return evaluate_series(y, cell['id'], replicate, phase, array_hash(z), equivalence)


def evaluate_series(y, cell_id='fixture', replicate=0, phase='validation', innovation_hash=None, equivalence=False):
    T = len(y)-1
    failures, warning_rows = [], []
    decisions = [[-1]*len(METHODS) for _ in ALPHAS]
    selected = [[-1]*len(METHODS) for _ in ALPHAS]
    row = dict(cell_id=cell_id, replicate=replicate, phase=phase, input_hash=array_hash(y),
               innovation_hash=innovation_hash, T=T, adf={}, ic={}, enders={},
               decisions=decisions, selected=selected, failures=failures, warnings=warning_rows,
               equivalence_checks=0)

    def call(method, alpha, fn):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            try:
                value = fn()
            except EXPECTED_FAILURES as exc:
                failures.append(dict(method=method, alpha=alpha, exception=type(exc).__name__, message=str(exc)))
                value = None
        warning_rows.extend(dict(method=method, alpha=alpha, category=w.category.__name__, message=str(w.message)) for w in caught)
        return value

    def check(r):
        if (r.lags, r.nobs, r.sample_start, r.sample_stop) != (0, T, 1, T+1):
            raise InvariantError('Fixed p=0/common-sample invariant failed')
        if hasattr(r, 'lag_selection') and r.lag_selection is not None:
            raise InvariantError('Automatic lag selection reached')

    def record(method, spec, result):
        j = METHODS.index(method)
        for i, alpha in enumerate(ALPHAS):
            decisions[i][j] = int(result.reject(alpha))
            selected[i][j] = SPECS.index(spec)

    for spec in SPECS:
        method = 'adf_'+spec
        r = call(method, None, lambda: adf(y, trend=spec, lags=0))
        if r is not None:
            check(r)
            row['adf'][spec] = dict(statistic=r.statistic, critical_values={str(a):v for a,v in r.critical_values.items()},
                                    pvalue=r.pvalue, ssr=float(r.regression.ssr), df_resid=float(r.regression.df_resid),
                                    coefficients=np.asarray(r.regression.params).tolist(), calibration=r.calibration)
            record(method, spec, r)
    for criterion in ('aic', 'bic'):
        method = 'det_'+criterion
        r = call(method, None, lambda: select_deterministic(y, lags=0, criterion=criterion))
        if r is not None:
            check(r)
            if r.nobs != T or tuple(c.specification for c in r.candidates) != SPECS:
                raise InvariantError('IC candidate sample/order invariant failed')
            for c in r.candidates:
                check(c.adf_result)
                base = row['adf'].get(c.specification)
                if base and (c.adf_result.statistic != base['statistic'] or c.ssr != base['ssr']):
                    raise InvariantError('IC/fixed-ADF fitted results differ')
            row['ic'][criterion] = dict(selected=r.selected_specification, gap=r.runner_up_gap,
                                        scores=[dict(specification=c.specification, aic=c.aic_score, bic=c.bic_score,
                                                     ssr=c.ssr, k=c.n_parameters, log_likelihood=c.log_likelihood)
                                                for c in r.candidates])
            record(method, r.selected_specification, r.selected_fixed_spec_adf)
    for i, alpha in enumerate(ALPHAS):
        records = {}
        for variant in ('original', 'conservative'):
            method = 'enders_'+variant
            r = call(method, alpha, lambda: enders(y, variant=variant, lags=0, alpha=alpha))
            if r is not None:
                check(r)
                records[variant] = enders_record(r)
                j = METHODS.index(method)
                decisions[i][j] = int(r.reject_unit_root)
                selected[i][j] = SPECS.index(r.trend)
        row['enders'][str(alpha)] = records
        if 'original' in records and 'conservative' in records:
            if records['conservative']['reject'] and not records['original']['reject']:
                raise InvariantError('Conservative rejection is outside Original rejection set')
        if equivalence and 'original' in records:
            permitted = call('equivalence', alpha, lambda: enders(y, variant='conservative', allow_no_constant=True, lags=0, alpha=alpha))
            if permitted is not None:
                if enders_record(permitted) != records['original']:
                    raise InvariantError('Permitted Conservative/Original history mismatch')
                row['equivalence_checks'] += 1
    if failures:
        row['input_levels_on_failure'] = np.asarray(y).tolist()
    return row
