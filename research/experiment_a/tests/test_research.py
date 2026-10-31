"""Catalogue/DGP/decision/storage validation independent of the study outcomes."""
from fractions import Fraction as F
import json
from pathlib import Path
import numpy as np
import pytest
from research.experiment_a.catalogue import catalogue,exact,resolved,digest,generate,innovations
from research.experiment_a import engine,runner


@pytest.mark.parametrize('T,total,nulls',[(100,22,5),(250,43,6),(500,35,6)])
def test_exact_catalogue(T,total,nulls):
    cells=[c for c in catalogue() if c['law']['T']==T]
    assert len(cells)==total and sum(c['is_null'] for c in cells)==nulls
    assert sum(len(c['aliases']) for c in cells)==51
    assert len({c['id'] for c in cells})==total
    for c in cells:
        assert digest(c['law'])==c['id']
        assert c['law']['p']==0 and c['law']['y0']==[0,1]
        for alias in c['aliases']:
            if alias['family']=='M':
                assert c['law']['intercept']==exact((1-F(*c['law']['rho']))*alias['mean'])
    # Rational zero/intercept aliases are identical laws, not near-duplicates.
    if T==100:
        c=next(c for c in cells if c['is_null'] and c['intercept']==.1)
        assert {a['drift'] for a in c['aliases']}=={'local_c1','fixed_0.1'}


@pytest.mark.parametrize('T',[100,250,500])
@pytest.mark.parametrize('rho',[.9,.95,1.])
def test_recursion_closed_form(T,rho):
    c=next(c for c in catalogue() if c['law']['T']==T and c['rho']==rho and c['intercept']>0)
    z=np.sin(np.arange(T)*.7)+np.cos(np.arange(T)*1.31)
    y=generate(c,z);r=c['rho'];a=c['intercept']
    expected=np.array([sum(r**(t-j)*(a+z[j-1]) for j in range(1,t+1)) for t in range(T+1)])
    np.testing.assert_allclose(y,expected,rtol=0,atol=1e-10)
    np.testing.assert_allclose(y[1:]-r*y[:-1]-a,z,rtol=0,atol=1e-12)
    assert y[0]==0 and len(y)==T+1


def test_stream_addressability():
    a,b=catalogue()[:2]
    z=innovations(a,7,'validation')
    assert np.array_equal(z,innovations(a,7,'validation'))
    assert not np.array_equal(z,innovations(a,8,'validation'))
    assert not np.array_equal(z,innovations(a,7,'benchmark'))
    assert not np.array_equal(z,innovations(b,7,'validation'))
    assert not np.array_equal(z,innovations(a,7,'main'))
    assert not np.array_equal(innovations(a,7,'pilot'),innovations(a,7,'main'))


@pytest.mark.parametrize('phase,reps',[('main',2000),('pilot',2000),('pilot',20),('validation',201)])
def test_scope_guard(tmp_path,phase,reps):
    with pytest.raises(ValueError):runner.run(tmp_path,catalogue()[:1],reps=reps,phase=phase)


def test_all_calls_fixed_and_failed_candidate_recorded(monkeypatch):
    c=catalogue()[0];y=generate(c,np.sin(np.arange(100)*1.1))
    actual=engine.adf;calls=[]
    def wrapped(y,*,trend,lags):
        calls.append((trend,lags))
        if trend=='n':raise ValueError('controlled candidate failure')
        return actual(y,trend=trend,lags=lags)
    monkeypatch.setattr(engine,'adf',wrapped)
    row=engine.evaluate_series(y)
    assert calls==[('n',0),('c',0),('ct',0)]
    assert [v[0] for v in row['decisions']]==[-1,-1,-1]
    assert row['failures'][0]['method']=='adf_n'
    assert 'input_levels_on_failure' in row


def test_repeated_and_checkpoint_resume(tmp_path):
    cells=catalogue()[:2]
    out=tmp_path/'resume';ref=tmp_path/'reference'
    first=runner.run(out,cells,reps=4,phase='validation',chunk_size=2,max_tasks=1)
    assert not first['complete'] and first['completed_rows']==2
    # A partial uncommitted task file must not count as a completed shard.
    (out/'shards'/'uncommitted.tmp').write_text('incomplete')
    resumed=runner.run(out,cells,reps=4,phase='validation',chunk_size=2)
    runner.run(ref,cells,reps=4,phase='validation',chunk_size=2)
    assert resumed['complete']
    _,a=runner.load_rows(out);_,b=runner.load_rows(ref)
    assert a==b and len(a)==8
    noop=runner.run(out,cells,reps=4,phase='validation',chunk_size=2)
    assert noop['sessions'][-1]['new_rows']==0
    with pytest.raises(ValueError,match='identity'):
        runner.run(out,cells,reps=5,phase='validation',chunk_size=2)
    done=next((out/'shards').glob('*.jsonl.gz'))
    done.write_bytes(done.read_bytes()+b'corrupt')
    with pytest.raises(ValueError,match='checksum'):
        runner.run(out,cells,reps=4,phase='validation',chunk_size=2)


def test_serial_parallel_parity(tmp_path):
    cells=catalogue()[:2]
    runner.run(tmp_path/'serial',cells,reps=3,phase='validation',chunk_size=1,workers=1)
    runner.run(tmp_path/'parallel',cells,reps=3,phase='validation',chunk_size=1,workers=2)
    assert runner.load_rows(tmp_path/'serial')[1]==runner.load_rows(tmp_path/'parallel')[1]


def test_immutable_grid_and_seven_methods():
    assert len(catalogue())==100 and sum(c['is_null'] for c in catalogue())==17
    assert len(engine.METHODS)==7 and engine.ALPHAS==(.01,.05,.1)
    assert len({digest(c['law']) for c in catalogue()})==100
