"""Replay all 312 saved Stage 3B model cases against the installed package.

The complete suite and audit_saved_results.py additionally verify selector
candidate equations and the historical Stage 2 comparison ledger. This script
never rewrites fixtures. Use --report to save environment-specific results.
"""
import argparse
import csv
import json
from pathlib import Path
import platform

import numpy as np
from arch.unitroot import ADF, DFGLS
from unitroot_tools import adf, dfgls, select_lags


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    def read(path):
        with path.open(newline='') as f:
            return list(csv.DictReader(f))
    columns = read(root/'results/series.csv')
    vectors = {name: np.array([float(row[name]) for row in columns if row[name]]) for name in columns[0]}
    rows = read(root/'stage3b_comparisons.csv')
    maxima = dict(saved_statistic=0., saved_pvalue=0., saved_critical_value=0., public_statistic=0., public_pvalue=0., public_critical_value=0.)
    selections = 0
    for row in rows:
        adapter, native = (adf, ADF) if row['test'] == 'ADF' else (dfgls, DFGLS)
        y = vectors[row['series']]
        if row['mode'] == 'fixed':
            result = adapter(y, trend=row['trend'], lags=int(row['p']))
        else:
            result = adapter(y, trend=row['trend'], method=row['method'], max_lags=int(row['K']))
            selected = select_lags(y, trend=row['trend'], method=row['method'], max_lags=int(row['K']))
            assert result.lags == selected.selected_lags == int(row['p'])
            assert result.lag_selection.selection_nobs == int(row['selection_nobs'])
            selections += 1
        assert result.nobs == int(row['nobs']) == len(y)-1-result.lags
        assert result.tail == 'left' and result.versions['arch'] == '8.0.0'
        backend = native(y, trend=row['trend'], lags=result.lags)
        maxima['saved_statistic'] = max(maxima['saved_statistic'], abs(result.statistic-float(row['statistic'])))
        maxima['saved_pvalue'] = max(maxima['saved_pvalue'], abs(result.pvalue-float(row['pvalue'])))
        maxima['public_statistic'] = max(maxima['public_statistic'], abs(result.statistic-backend.stat))
        maxima['public_pvalue'] = max(maxima['public_pvalue'], abs(result.pvalue-backend.pvalue))
        for alpha, field, key in ((.01,'cv1','1%'),(.05,'cv5','5%'),(.1,'cv10','10%')):
            maxima['saved_critical_value'] = max(maxima['saved_critical_value'], abs(result.critical_values[alpha]-float(row[field])))
            maxima['public_critical_value'] = max(maxima['public_critical_value'], abs(result.critical_values[alpha]-backend.critical_values[key]))
            assert result.reject(alpha) == (result.statistic < result.critical_values[alpha])
    assert len(rows) == 312 and selections == 240
    assert max(maxima.values()) <= 1e-8
    args.report.write_text(json.dumps(dict(python=platform.python_version(), versions=dict(result.versions), model_pairs=len(rows), automatic_selections=selections, absolute_tolerance=1e-8, max_absolute_differences=maxima), indent=2)+'\n')
    print(json.dumps(maxima))


if __name__ == '__main__':
    main()
