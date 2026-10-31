"""Operational main guard and fail-closed local recovery (no main computation)."""
import json
import os
import platform
import subprocess
import sys
from pathlib import Path
import pytest
from research.experiment_a import runner
from research.experiment_a.catalogue import catalogue


def dead_pid():
    child=subprocess.Popen([sys.executable,'-c','pass'])
    child.wait()
    return child.pid


def lock_value(pid):
    return dict(lock_version=1,process_namespace=runner.process_namespace(),host=platform.node(),uid=os.getuid(),pid=pid,token='fixture')


def test_main_config_without_computation(tmp_path,monkeypatch):
    def stop(task,token):
        raise RuntimeError('configuration accepted; computation deliberately stopped')
    monkeypatch.setattr(runner,'locked_chunk',stop)
    with pytest.raises(RuntimeError,match='deliberately stopped'):
        runner.run(tmp_path,reps=2000,phase='main',allow_main=True)
    m=json.loads((tmp_path/'manifest.json').read_text())
    assert len(m['catalogue'])==100 and len(m['tasks'])==4000
    assert sum(t['stop']-t['start'] for t in m['tasks'])==200000
    assert m['phase']=='main' and not list((tmp_path/'shards').iterdir())
    assert not (tmp_path/'coordinator.lock').exists()
    for kwargs in [dict(reps=200),dict(cells=catalogue()[:1]),dict(cells=catalogue()[::-1])]:
        params=dict(reps=2000,phase='main',allow_main=True);params.update(kwargs)
        with pytest.raises(ValueError):runner.run(tmp_path/'invalid',**params)


def test_recovery_preserves_evidence_and_rejects_ambiguous(tmp_path):
    lock=tmp_path/'coordinator.lock'
    for value in [lock_value(os.getpid()),dict(lock_value(dead_pid()),host='foreign-host'),
                  dict(lock_value(dead_pid()),uid=-1),dict(pid=dead_pid(),host=platform.node())]:
        raw=json.dumps(value);lock.write_text(raw)
        with pytest.raises(RuntimeError):runner.recover_lock(tmp_path)
        assert lock.read_text()==raw
    raw=json.dumps(lock_value(dead_pid()));lock.write_text(raw)
    with runner.lease(tmp_path):
        with pytest.raises(RuntimeError,match='lease is live'):runner.recover_lock(tmp_path)
    archived=Path(runner.recover_lock(tmp_path))
    assert archived.read_text()==raw and not lock.exists()
    assert archived.with_suffix('.audit.json').exists()
    with pytest.raises(RuntimeError,match='expired coordinator'):
        runner.locked_chunk((str(tmp_path/'shards'),catalogue()[0],0,1,'validation'),'fixture')


def test_order_chunk_and_resume_invariance(tmp_path):
    cells=catalogue()[:2]
    runner.run(tmp_path/'a',cells,4,'validation',1,2)
    runner.run(tmp_path/'b',cells[::-1],4,'validation',2,1,max_tasks=3)
    runner.run(tmp_path/'b',cells[::-1],4,'validation',1,1)
    key=lambda row:(row['cell_id'],row['replicate'])
    assert sorted(runner.load_rows(tmp_path/'a')[1],key=key)==sorted(runner.load_rows(tmp_path/'b')[1],key=key)
    m=json.loads((tmp_path/'b'/'manifest.json').read_text())
    assert all(s['started_utc']<=s['finished_utc'] for s in m['sessions'])


def test_surviving_process_lease_blocks_recovery(tmp_path):
    lock=tmp_path/'coordinator.lock'
    lock.write_text(json.dumps(lock_value(dead_pid())))
    child=subprocess.Popen([sys.executable,'-c',
        'import sys; from research.experiment_a.runner import lease; '
        'h=lease(sys.argv[1]); print("ready",flush=True); sys.stdin.read()',str(tmp_path)],
        stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
    try:
        assert child.stdout.readline().strip()=='ready'
        with pytest.raises(RuntimeError,match='lease is live'):runner.recover_lock(tmp_path)
        assert lock.exists()
    finally:
        child.communicate('',timeout=10)
    assert Path(runner.recover_lock(tmp_path)).exists()
