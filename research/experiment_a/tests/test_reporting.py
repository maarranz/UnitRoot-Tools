import csv
import gzip
import numpy as np
import pytest
from research.experiment_a import report,runner
from research.experiment_a.catalogue import catalogue


def test_intervals_and_complete_exports(tmp_path):
    assert report.wilson(0,200)[0]==0
    assert 0<report.wilson(0,200)[1]<.02
    assert report.exact_interval(0,200)[0]==0
    assert report.exact_interval(200,200)[1]==1
    out=tmp_path/'run';summary=tmp_path/'summary'
    runner.run(out,catalogue()[:1],reps=4,phase='validation',chunk_size=2)
    result=report.summarize(out,summary)
    _,rows=runner.load_rows(out)
    assert result['cell_replications']==4 and result['method_alpha_decisions']==84
    with gzip.open(summary/'rejection_rates.csv.gz','rt') as f:
        values=list(csv.DictReader(f))
    assert len(values)==21
    D=np.array([r['decisions'] for r in rows])
    for row in values:
        i=report.ALPHAS.index(float(row['alpha']));j=report.METHODS.index(row['method'])
        assert int(row['rejections'])==int(D[:,i,j].sum())
        assert float(row['rejection_rate'])==D[:,i,j].mean()
        assert int(row['failures'])==0
    with gzip.open(summary/'paired_differences.csv.gz','rt') as f:
        pairs=list(csv.DictReader(f))
    assert len(pairs)==63
    for p in pairs:
        i=report.ALPHAS.index(float(p['alpha']));a=report.METHODS.index(p['method_a']);b=report.METHODS.index(p['method_b'])
        d=D[:,i,a]-D[:,i,b]
        assert float(p['difference'])==d.mean()
        assert float(p['mcse'])==pytest.approx(d.std(ddof=1)/2)
        assert float(p['ci95_low'])<=float(p['ci95_high'])
    primary=[p for p in pairs if p['primary_5pct']=='True']
    assert len(primary)==6
    assert any(p['method_a']=='det_aic' and p['method_b']=='adf_ct' for p in primary)
