"""Reproduce Stage 3C comparisons without changing reviewed Stage 2 fixtures.

Run from any directory with this package and its test dependencies installed:
    python validation/validate_stage3c.py
Writes stage3c_comparisons.csv and stage3c_numerics.json beside this script.
R comparisons use saved Stage 2 output, not a newly executed R installation.
"""
import csv
import hashlib
import json
import platform
from pathlib import Path
import warnings

import arch
from arch.unitroot import KPSS, PhillipsPerron
import numpy as np
import scipy
import statsmodels
from unitroot_tools import kpss, pp

ROOT = Path(__file__).resolve().parent
FIXTURES = ROOT / 'results'
TOL = 1e-8


def read(name):
    with (FIXTURES / name).open(newline='') as f:
        return list(csv.DictReader(f))


def main():
    assert arch.__version__ == '8.0.0'
    manifest = json.loads((FIXTURES / 'manifest.json').read_text())
    columns = read('series.csv')
    vectors = {name: np.array([row[name] for row in columns[:n]], dtype=float)
               for name, n in manifest['lengths'].items()}
    saved = read('python_release.csv')
    r_rows = read('r_results.csv')
    rows = []
    def compare(group, name, test, trend, b, field, actual, reference, expected):
        difference = abs(float(actual) - float(reference))
        rows.append(dict(group=group, series=name, test=test, trend=trend,
                         bandwidth=b, field=field, actual=float(actual), reference=float(reference),
                         absolute_difference=difference, expected_equal=expected,
                         within_tolerance=difference <= TOL))
        if expected:
            assert difference <= TOL, rows[-1]
    for name, y in vectors.items():
        for test, adapter, native_class, trends, types in [
            ('PP', pp, PhillipsPerron, ('n', 'c', 'ct'), ('tau', 'rho')),
            ('KPSS', kpss, KPSS, ('c', 'ct'), (None,))]:
            for trend in trends:
                for typ in types:
                    label = test + ('_' + typ if typ else '')
                    for b in (0, 4, 12, None):
                        kw = dict(trend=trend)
                        if typ:
                            kw['test_type'] = typ
                        result = adapter(y, bandwidth=b, **kw)
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore', DeprecationWarning)
                            # PP's capped rule is tested separately on short samples.
                            native = native_class(y, lags=b, **kw)
                            _ = native.stat
                        actual = dict(statistic=result.statistic, pvalue=result.pvalue,
                                      nobs=result.nobs, bandwidth=result.bandwidth,
                                      cv1=result.critical_values[.01], cv5=result.critical_values[.05],
                                      cv10=result.critical_values[.10])
                        reference = dict(statistic=native.stat, pvalue=native.pvalue,
                                         nobs=native.nobs, bandwidth=native.lags,
                                         cv1=native.critical_values['1%'], cv5=native.critical_values['5%'],
                                         cv10=native.critical_values['10%'])
                        for field in actual:
                            compare('fresh_ARCH', name, label, trend, b, field,
                                    actual[field], reference[field], True)
                        if b not in (4, 12) or trend == 'n':
                            continue
                        def match(row):
                            return (row['series'], row['test'], row['trend'], row['lags']) == (name, label, trend, str(b))
                        frozen = next(row for row in saved if row['backend'] == 'ARCH' and match(row))
                        for field in ('statistic', 'pvalue', 'nobs', 'cv1', 'cv5', 'cv10'):
                            compare('saved_ARCH', name, label, trend, b, field, actual[field], frozen[field], True)
                        if test == 'PP':
                            peer = next(row for row in r_rows if match(row))
                            compare('saved_urca', name, label, trend, b, 'statistic', actual['statistic'], peer['statistic'], False)
                        else:
                            peer = next(row for row in saved if row['backend'] == 'Statsmodels' and match(row))
                            for field in ('statistic', 'pvalue', 'cv1', 'cv5', 'cv10'):
                                compare('saved_Statsmodels', name, label, trend, b, field, actual[field], peer[field], field == 'statistic')
    with (ROOT / 'stage3c_comparisons.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    groups = {}
    for row in rows:
        key = row['group'] + ':' + row['test'] + ':' + row['field']
        item = groups.setdefault(key, dict(comparisons=0, max_absolute_difference=0., expected_equal=row['expected_equal'], outside_tolerance=0))
        item['comparisons'] += 1
        item['max_absolute_difference'] = max(item['max_absolute_difference'], row['absolute_difference'])
        item['outside_tolerance'] += int(not row['within_tolerance'])
    metadata = dict(python=platform.python_version(), platform=platform.platform(),
                    versions=dict(arch=arch.__version__, numpy=np.__version__, scipy=scipy.__version__, statsmodels=statsmodels.__version__),
                    absolute_tolerance=TOL, comparisons=len(rows),
                    expected_equality_checks=sum(row['expected_equal'] for row in rows),
                    expected_equality_failures=sum(row['expected_equal'] and not row['within_tolerance'] for row in rows),
                    fresh_ARCH_model_pairs=192, saved_ARCH_model_pairs=72,
                    saved_urca_model_pairs=48, saved_Statsmodels_model_pairs=24,
                    fixture_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(FIXTURES.iterdir()) if p.is_file()},
                    groups=groups,
                    limitations=['Saved R references only; no live R regeneration.', 'Six reviewed series are regression evidence, not a size/power study.', 'Python 3.11, lower dependency bounds and other platforms untested.'])
    (ROOT / 'stage3c_numerics.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps({key: metadata[key] for key in ('comparisons', 'expected_equality_checks', 'expected_equality_failures')}))


if __name__ == '__main__':
    main()
