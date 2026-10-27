#!/usr/bin/env python3
"""Stage 2 research harness, not UnitRoot-Tools implementation.

Requires release arch 8.0.0, statsmodels 0.15.0, NumPy/SciPy/Pandas and R
urca/bootUR/jsonlite. A separate installed ARCH development snapshot is mandatory.
Nothing is installed or downloaded automatically. See docs/Numerical_Validation.md.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import warnings

import numpy as np

ARCH_COMMIT = '09a6be5f1b490346d8ba5d3b5850991734aaaf3a'
# Filled from the immutable upstream implementation inspected in Stage 1.
ARCH_SOURCE_SHA256 = 'b6be3547403ab880d425c780e1358c86eb9073a7b903ad8975650b45ad445fc9'
FIELDS = ['series', 'backend', 'test', 'trend', 'lags', 'variant', 'statistic',
          'pvalue', 'nobs', 'cv1', 'cv5', 'cv10', 'estimate', 'reported_lags', 'criterion_value']


def write_csv(path, rows, fields=None):
    if not rows:
        raise ValueError(f'No rows for {path}')
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def read_csv(path):
    with Path(path).open(newline='') as f:
        return list(csv.DictReader(f))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def deterministics(n, trend):
    t = np.arange(1, n + 1, dtype=float)
    if trend == 'n':
        return np.empty((n, 0))
    return np.column_stack([t**j for j in range({'c': 1, 'ct': 2, 'ctt': 3}[trend])])


def detrend(y, trend, gls=False):
    z = deterministics(len(y), trend)
    if z.shape[1] == 0:
        return y.copy()
    if gls:
        a = 1 - (7 if trend == 'c' else 13.5) / len(y)
        qy, qz = y.copy(), z.copy()
        qy[1:] = y[1:] - a * y[:-1]
        qz[1:] = z[1:] - a * z[:-1]
        beta = np.linalg.lstsq(qz, qy, rcond=None)[0]
    else:
        beta = np.linalg.lstsq(z, y, rcond=None)[0]
    return y - z @ beta


def mic_curve(y, trend, K, variant):
    """Independent least-squares diagnostic; no ARCH or R selector is called.

    PQ_full: full-series OLS detrending and full matched lag-level sum.
    PQ_bootUR: same regression, penalty drops first matched lag-level value.
    NP_GLS: original-style GLS detrending for selection, full sum.
    one_step_ARCH / one_step_bootUR: deterministic terms fitted jointly.
    """
    one = variant.startswith('one_step')
    v = y.copy() if one else detrend(y, trend, gls=(variant == 'NP_GLS'))
    dy = np.diff(v)
    lhs = dy[K:]
    level = v[K:-1]
    z = deterministics(len(lhs), trend) if one else np.empty((len(lhs), 0))
    residual_level = level - z @ np.linalg.lstsq(z, level, rcond=None)[0]
    ss = (level[1:] @ level[1:]) if variant.endswith('bootUR') else residual_level @ residual_level
    n = len(lhs)
    out = []
    for p in range(K + 1):
        x = np.column_stack([level] + [dy[K-j:len(dy)-j] for j in range(1, p + 1)] + [z])
        beta = np.linalg.lstsq(x, lhs, rcond=None)[0]
        e = lhs - x @ beta
        var = e @ e / n
        tau = beta[0]**2 * ss / var
        out.append(dict(p=p, n=n, beta=float(beta[0]), sigma2=float(var),
                        level_ss=float(ss), tau=float(tau),
                        MAIC=float(np.log(var) + 2 * (tau + p) / n),
                        MBIC=float(np.log(var) + np.log(n) * (tau + p) / n)))
    return out


def simulated_series():
    # PCG64 named explicitly; same floating-point CSV values are read in R.
    rng = np.random.Generator(np.random.PCG64(20261027))
    e = rng.standard_normal((5, 501))
    ar = np.zeros(501)
    for t in range(1, 501):
        ar[t] = .7 * ar[t-1] + e[2, t]
    tr = np.zeros(501)
    for t in range(1, 501):
        tr[t] = .8 * tr[t-1] + e[3, t]
    data = {
        'rw': np.cumsum(e[0, 1:])[-200:],
        'rw_negma': np.cumsum(e[1, 1:] - .8 * e[1, :-1])[-200:],
        'ar07': ar[-200:],
        'trend_ar08': tr[-200:] + 2 + .03 * np.arange(1, 201),
        'hetero_rw': np.cumsum(np.r_[e[4, 1:101], 3 * e[4, 101:201]]),
    }
    # Deterministic diagnostic search, not a random size/power experiment. Save the
    # first seed in 0..999 exposing bootUR's omitted penalty observation.
    witness = None
    for seed in range(1000):
        u = np.random.Generator(np.random.PCG64(seed)).standard_normal(65)
        y = np.cumsum(u[1:] - .8 * u[:-1])
        a, b = mic_curve(y, 'c', 4, 'PQ_full'), mic_curve(y, 'c', 4, 'PQ_bootUR')
        if np.argmin([r['MAIC'] for r in a]) != np.argmin([r['MAIC'] for r in b]):
            data['mic_boundary'] = y
            witness = seed
            break
    if witness is None:
        raise RuntimeError('Diagnostic boundary search found no witness')
    return data, witness


def load_series(path):
    rows = read_csv(path)
    return {k: np.array([float(r[k]) for r in rows if r[k] != '']) for k in rows[0]}


def package_row(id, backend, test, tr, lag, variant, f):
    cv = f.critical_values
    return dict(zip(FIELDS, [id, backend, test, tr, lag, variant, float(f.stat),
        float(f.pvalue), int(f.nobs), float(cv['1%']), float(cv['5%']), float(cv['10%']), '']))


def worker(args):
    import arch
    import arch.unitroot.unitroot as au
    import statsmodels
    from arch.unitroot import ADF, DFGLS, KPSS, PhillipsPerron
    from statsmodels.tsa.stattools import adfuller, kpss
    from statsmodels.tools.sm_exceptions import InterpolationWarning
    data = load_series(args.series)
    rows = []
    if args.worker == 'development':
        if digest(au.__file__) != ARCH_SOURCE_SHA256:
            raise RuntimeError('Development unitroot.py hash does not match audited revision')
        for id, y in data.items():
            for tr in ('n', 'c', 'ct'):
                for K in (4, 12):
                    for lm in (False, True):
                        variant = f'K{K}_' + ('low' if lm else 'full')
                        f = ADF(y, trend=tr, max_lags=K, method='maic', low_memory=lm)
                        row = package_row(id, 'ARCH_dev', 'MAIC_ADF', tr, f.lags, variant, f)
                        row['criterion_value'] = float(f._ic_best)  # Audit diagnostic only.
                        rows.append(row)
                        if tr != 'n':
                            f = DFGLS(y, trend=tr, max_lags=K, method='maic', low_memory=lm)
                            row = package_row(id, 'ARCH_dev', 'MAIC_DFGLS', tr, f.lags, variant, f)
                            value, selected = au._df_select_lags(detrend(y, tr), 'n', K, 'maic', low_memory=lm)
                            if selected != f.lags:
                                raise AssertionError('DFGLS selector/class disagreement')
                            row['criterion_value'] = float(value)
                            rows.append(row)
    else:
        if arch.__version__ != '8.0.0' or statsmodels.__version__ != '0.15.0':
            raise RuntimeError('Release worker requires arch==8.0.0 and statsmodels==0.15.0')
        for id, y in data.items():
            for tr in ('n', 'c', 'ct', 'ctt'):
                for p in (0, 4):
                    f = ADF(y, lags=p, trend=tr)
                    rows.append(package_row(id, 'ARCH', 'ADF', tr, p, 'fixed', f))
                    s = adfuller(y, maxlag=p, regression=tr, autolag=None, result_object=True)
                    rows.append(dict(zip(FIELDS, [id, 'Statsmodels', 'ADF', tr, p, 'fixed',
                        s.statistic, s.pvalue, s.nobs, s.critical_values['1%'],
                        s.critical_values['5%'], s.critical_values['10%'], ''])))
            for tr in ('c', 'ct'):
                for b in (4, 12):
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', InterpolationWarning)  # Table clipping is reported.
                        f = KPSS(y, lags=b, trend=tr)
                        s = kpss(y, regression=tr, nlags=b, result_object=True)
                    rows.append(package_row(id, 'ARCH', 'KPSS', tr, b, 'fixed', f))
                    rows.append(dict(zip(FIELDS, [id, 'Statsmodels', 'KPSS', tr, b, 'fixed',
                        s.statistic, s.pvalue, len(y), s.critical_values['1%'],
                        s.critical_values['5%'], s.critical_values['10%'], ''])))
                    for typ in ('tau', 'rho'):
                        f = PhillipsPerron(y, lags=b, trend=tr, test_type=typ)
                        rows.append(package_row(id, 'ARCH', 'PP_' + typ, tr, b, 'fixed', f))
                for p in (0, 4):
                    f = DFGLS(y, lags=p, trend=tr)
                    rows.append(package_row(id, 'ARCH', 'DFGLS', tr, p, 'fixed', f))
    write_csv(args.result, rows, FIELDS)
    metadata = dict(python=sys.version.split()[0], platform=platform.platform(),
        packages={p: importlib.metadata.version(p) for p in
                  ('numpy', 'scipy', 'pandas', 'arch', 'statsmodels', 'patsy')},
        all_distributions={d.metadata['Name']: d.version for d in importlib.metadata.distributions()},
        arch_unitroot_sha256=digest(au.__file__), arch_runtime_version=arch.__version__)
    # --target metadata may resolve a different installed distribution; runtime
    # import version and source digest above identify the development code.
    Path(str(args.result) + '.metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')


def pp_reference(y, trend, bandwidth, use_current):
    """Reconcile finite-sample PP correction using lagged vs current levels."""
    lhs, level = y[1:], y[:-1]
    z = deterministics(len(lhs), trend)
    x = np.column_stack([level, z])
    b = np.linalg.lstsq(x, lhs, rcond=None)[0]
    u = lhs - x @ b
    n, k = x.shape
    gamma = u @ u / n
    lam = gamma + sum(2 * (1-j/(bandwidth+1)) * (u[j:] @ u[:-j]) / n
                      for j in range(1, bandwidth+1))
    se = np.sqrt((u @ u / (n-k)) * np.linalg.inv(x.T @ x)[0, 0])
    v = lhs if use_current else level
    vr = v - z @ np.linalg.lstsq(z, v, rcond=None)[0]
    M = vr @ vr / n**2
    if use_current and trend == 'ct':
        # urca's explicit trend-case M has this finite-sample factor.
        M *= 1 - n**-2
    tau = np.sqrt(gamma/lam) * (b[0]-1) / se - .5 * (lam-gamma)/np.sqrt(lam*M)
    rho = n * (b[0]-1) - .5 * (lam-gamma)/M
    return {'tau': tau, 'rho': rho}


def summarize(out, data):
    raw = sum([read_csv(out / f) for f in ('python_release.csv', 'python_development.csv', 'r_results.csv')], [])
    comps = []
    def one(backend, id, test, tr, lag=None, variant='fixed'):
        r = [r for r in raw if r['backend'] == backend and r['series'] == id and
             r['test'] == test and r['trend'] == tr and r['variant'] == variant and
             (lag is None or int(float(r['lags'])) == lag)]
        if len(r) != 1:
            raise ValueError((backend, id, test, tr, lag, variant, len(r)))
        return r[0]
    def cmp(group, id, tr, lag, field, a, b, expected=False):
        a, b = float(a), float(b)
        delta = abs(a-b)
        comps.append(dict(group=group, series=id, trend=tr, lags=lag, field=field,
                          a=a, b=b, abs_difference=delta, expected_equal=expected,
                          within_1e_8=bool(delta <= 1e-8)))
    for id, y in data.items():
        for test, trends, lags in [('ADF', ('n','c','ct','ctt'), (0,4)), ('KPSS', ('c','ct'), (4,12)),
                                    ('DFGLS', ('c','ct'), (0,4))]:
            other = 'R' if test == 'DFGLS' else 'Statsmodels'
            for tr in trends:
                for p in lags:
                    a, b = one('ARCH', id, test, tr, p), one(other, id, test, tr, p)
                    for field in ('statistic','cv1','cv5','cv10'):
                        cmp(test, id, tr, p, field, a[field], b[field], field == 'statistic' or test == 'ADF')
                    if test != 'DFGLS':
                        cmp(test, id, tr, p, 'pvalue', a['pvalue'], b['pvalue'], test == 'ADF')
        for tr in ('c','ct'):
            for bw in (4,12):
                for typ in ('tau','rho'):
                    a, b = one('ARCH', id, 'PP_'+typ, tr, bw), one('R', id, 'PP_'+typ, tr, bw)
                    cmp('PP_'+typ, id, tr, bw, 'statistic', a['statistic'], b['statistic'])
                    for backend, row, current in [('ARCH',a,False), ('R',b,True)]:
                        cmp('PP_reconciled_'+backend, id, tr, bw, typ, row['statistic'],
                            pp_reference(y, tr, bw, current)[typ], True)
    candidate_rows = []
    for id, y in data.items():
        for tr in ('n','c','ct'):
            for K in (4,12):
                curves = {v: mic_curve(y,tr,K,v) for v in
                          ('PQ_full','PQ_bootUR','NP_GLS','one_step_ARCH','one_step_bootUR')}
                for v, curve in curves.items():
                    for row in curve:
                        candidate_rows.append(dict(series=id,trend=tr,K=K,variant=v,**row))
                for crit in ('MAIC','MBIC'):
                    for step, ref in [('two_step','PQ_bootUR'), ('one_step','one_step_bootUR')]:
                        r = one('R',id,crit,tr,variant=f'{step}_K{K}_unscaled')
                        chosen = int(np.argmin([x[crit] for x in curves[ref]]))
                        cmp('bootUR_'+crit+'_oracle_'+step,id,tr,K,'selected_lag',r['lags'],chosen,True)
                        if step == 'one_step':
                            dc = {'n':0, 'c':1, 'ct':2}[tr]
                            cmp('bootUR_'+crit+'_raw_lag_reporting',id,tr,K,'reported_lag',
                                r['reported_lags'],chosen+dc,True)
                    if crit == 'MAIC':
                        for test, ref in [('MAIC_ADF','one_step_ARCH'),('MAIC_DFGLS','PQ_full')]:
                            if test.endswith('DFGLS') and tr == 'n':
                                continue
                            a = one('ARCH_dev',id,test,tr,variant=f'K{K}_full')
                            b = one('ARCH_dev',id,test,tr,variant=f'K{K}_low')
                            chosen = int(np.argmin([x[crit] for x in curves[ref]]))
                            cmp('ARCH_dev_'+test+'_oracle',id,tr,K,'selected_lag',a['lags'],chosen,True)
                            cmp('ARCH_dev_'+test+'_memory',id,tr,K,'selected_lag',a['lags'],b['lags'],True)
                            for route, native in [('full', a), ('low', b)]:
                                cmp('ARCH_dev_'+test+'_criterion_'+route,id,tr,K,'minimum_MAIC',
                                    native['criterion_value'], min(x[crit] for x in curves[ref]), True)
                            r = one('R',id,crit,tr,variant=f'{"two_step" if test.endswith("DFGLS") else "one_step"}_K{K}_unscaled')
                            cmp('ARCH_bootUR_'+test,id,tr,K,'selected_lag',a['lags'],r['lags'])
                    if crit == 'MBIC':
                        # No native ARCH MBIC; compare recommended convention to bootUR.
                        chosen = int(np.argmin([x[crit] for x in curves['PQ_full']]))
                        r = one('R',id,crit,tr,variant=f'two_step_K{K}_unscaled')
                        cmp('MBIC_PQ_full_vs_bootUR',id,tr,K,'selected_lag',chosen,r['lags'])
                u = one('R',id,'MAIC',tr,variant=f'two_step_K{K}_unscaled')
                s = one('R',id,'MAIC',tr,variant=f'two_step_K{K}_scaled')
                cmp('bootUR_rescaling_sensitivity',id,tr,K,'selected_lag',u['lags'],s['lags'])
    write_csv(out / 'mic_candidates.csv', candidate_rows)
    write_csv(out / 'comparisons.csv', comps)
    groups = {}
    for group in sorted({r['group'] for r in comps}):
        r = [r for r in comps if r['group'] == group]
        groups[group] = dict(comparisons=len(r),max_abs_difference=max(x['abs_difference'] for x in r),
                             differences_over_1e_8=sum(not x['within_1e_8'] for x in r),
                             expected_failures=sum(x['expected_equal'] and not x['within_1e_8'] for x in r))
    (out / 'summary.json').write_text(json.dumps(groups, indent=2)+'\n')
    failures = [r for r in comps if r['expected_equal'] and not r['within_1e_8']]
    print(json.dumps(groups, indent=2))
    if failures:
        raise AssertionError(f'{len(failures)} unexplained diagnostic mismatches; inspect comparisons.csv')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parent/'results')
    parser.add_argument('--arch-dev-path', type=Path)
    parser.add_argument('--rscript', default='Rscript')
    parser.add_argument('--r-library', default='')
    parser.add_argument('--worker', choices=['release','development'], help=argparse.SUPPRESS)
    parser.add_argument('--series', type=Path, help=argparse.SUPPRESS)
    parser.add_argument('--result', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if args.arch_dev_path is None:
        parser.error('--arch-dev-path is required (installed audited snapshot, outside this repository)')
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    data, witness = simulated_series()
    write_csv(out/'series.csv', [{k: (format(v[i], '.17g') if i<len(v) else '') for k,v in data.items()}
                                for i in range(max(map(len,data.values())))])
    script = str(Path(__file__).resolve())
    # Child imports are isolated: only the development worker gets the override.
    env = os.environ.copy(); env.pop('PYTHONPATH',None)
    for mode in ('release','development'):
        child_env = env.copy()
        if mode == 'development':
            child_env['PYTHONPATH'] = str(args.arch_dev_path.resolve())
        subprocess.run([sys.executable, script, '--worker', mode, '--series', str(out/'series.csv'),
                        '--result', str(out/f'python_{mode}.csv')], env=child_env, check=True)
    subprocess.run([args.rscript, '--vanilla', str(Path(__file__).with_suffix('.R').resolve()),
                    str(out/'series.csv'), str(out/'r_results.csv'), str(out/'r_metadata.json'),
                    args.r_library], check=True)
    manifest = dict(arch_commit=ARCH_COMMIT, arch_source_sha256=ARCH_SOURCE_SHA256,
                    generator='NumPy Generator(PCG64)', seed=20261027,
                    mic_boundary_search=dict(seeds='0..999 in ascending order', selected_seed=witness,
                                             T=64, theta=-.8, K=4, trend='c'),
                    tolerance=dict(absolute=1e-8, selected_lags='exact'),
                    lengths={k:len(v) for k,v in data.items()},
                    file_sha256={name:digest(out/name) for name in
                                 ('series.csv','python_release.csv','python_development.csv','r_results.csv')},
                    script_sha256={Path(script).name:digest(script),
                                   'cross_validate.R':digest(Path(__file__).with_suffix('.R'))})
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    summarize(out, data)


if __name__ == '__main__':
    main()
