"""Atomic, resumable Experiment A reference runner. Main execution requires an explicit opt-in."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import scipy
import fcntl
import uuid
import gzip
import hashlib
import importlib.metadata
import json
import multiprocessing
import os
from pathlib import Path
import platform
import resource
import signal
import sys
import time

import numpy as np
import unitroot_tools
from .catalogue import MASTER_SEED, PHASES, RNG_VERSION, canonical, catalogue, digest
from .engine import ALPHAS, METHODS, evaluate


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n')
    os.replace(temporary, path)


def rss_bytes():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == 'darwin' else value*1024)


def provenance():
    here = Path(__file__).resolve().parent
    package = Path(unitroot_tools.__file__).resolve().parent
    build = dict(numpy=np.show_config(mode='dicts'), scipy=scipy.show_config(mode='dicts'))
    return dict(numerical_build_sha256=digest(build),
                numerical_threads={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','OPENBLAS_CORETYPE')},
                rng=RNG_VERSION, master_seed_hex=f'{MASTER_SEED:032x}',
                methods=list(METHODS), alphas=list(ALPHAS), python=platform.python_version(),
                implementation=platform.python_implementation(), platform=platform.platform(),
                architecture=platform.machine(), module_path=str(package),
                versions={n:importlib.metadata.version(n) for n in ('arch','numpy','scipy','statsmodels','pandas')},
                unitroot_version=unitroot_tools.__version__,
                package_sha256={p.name:sha(p) for p in sorted(package.glob('*.py'))},
                runner_sha256={n:sha(here/n) for n in ('__init__.py','catalogue.py','engine.py','runner.py')})


def task_id(cell, start, stop):
    return f"{cell['id']}_{start:06d}_{stop:06d}"


def compute_chunk(task):
    root, cell, start, stop, phase = task
    root = Path(root)
    identifier = task_id(cell, start, stop)
    start_time = time.perf_counter()
    start_cpu = time.process_time()
    rows = [evaluate(cell, rep, phase, equivalence=(phase == 'pilot' and rep < 20))
            for rep in range(start, stop)]
    computed = time.perf_counter()
    raw = b''.join((canonical(row)+'\n').encode() for row in rows)
    gz = root/(identifier+'.jsonl.gz')
    tmp = gz.with_suffix('.tmp')
    with tmp.open('wb') as f:
        with gzip.GzipFile(filename='', mode='wb', fileobj=f, mtime=0, compresslevel=6) as z:
            z.write(raw)
    os.replace(tmp, gz)
    npz = root/(identifier+'.npz')
    tmp = npz.with_suffix('.tmp')
    with tmp.open('wb') as f:
        np.savez_compressed(f, replicate=np.arange(start, stop, dtype=np.int32),
                            decisions=np.array([r['decisions'] for r in rows], dtype=np.int8),
                            selected=np.array([r['selected'] for r in rows], dtype=np.int8))
    os.replace(tmp, npz)
    done = dict(task_id=identifier, cell_id=cell['id'], start=start, stop=stop, rows=len(rows),
                hashes={gz.name:sha(gz), npz.name:sha(npz)}, input_hashes=[r['input_hash'] for r in rows],
                compute_seconds=computed-start_time, cpu_seconds=time.process_time()-start_cpu, elapsed_seconds=time.perf_counter()-start_time,
                process_id=os.getpid(), worker_peak_rss_bytes=rss_bytes(),
                failure_records=sum(len(r['failures']) for r in rows),
                warning_records=sum(len(r['warnings']) for r in rows),
                equivalence_checks=sum(r['equivalence_checks'] for r in rows))
    atomic_json(root/(identifier+'.done.json'), done)  # Commit marker is always last.
    return done


def verified_done(root, task):
    root = Path(root)
    _, cell, start, stop, _ = task
    identifier = task_id(cell, start, stop)
    path = root/(identifier+'.done.json')
    if not path.exists():
        return None
    value = json.loads(path.read_text())
    if (value['task_id'],value['cell_id'],value['start'],value['stop'],value['rows']) != (identifier,cell['id'],start,stop,stop-start):
        raise ValueError('Corrupt checkpoint metadata: '+identifier)
    expected = {identifier+'.jsonl.gz', identifier+'.npz'}
    if set(value['hashes']) != expected or any(sha(root/n) != h for n,h in value['hashes'].items()):
        raise ValueError('Checkpoint checksum mismatch: '+identifier)
    return value


def iter_chunks(root, manifest):
    root = Path(root)
    for t in manifest['tasks']:
        cell = next(c for c in manifest['catalogue'] if c['id'] == t['cell_id'])
        task = (str(root/'shards'),cell,t['start'],t['stop'],manifest['phase'])
        done = verified_done(root/'shards',task)
        if done is None:
            continue
    # Iterate in declared cell/range order, not filesystem/completion order.
        identifier = done['task_id']
        with gzip.open(root/'shards'/(identifier+'.jsonl.gz'),'rt') as f:
            chunk = [json.loads(line) for line in f]
        if [(r['cell_id'],r['replicate']) for r in chunk] != [(cell['id'],i) for i in range(t['start'],t['stop'])]:
            raise ValueError('Checkpoint row keys mismatch')
        with np.load(root/'shards'/(identifier+'.npz'),allow_pickle=False) as arrays:
            if not (np.array_equal(arrays['replicate'],np.arange(t['start'],t['stop']))
                    and np.array_equal(arrays['decisions'],[r['decisions'] for r in chunk])
                    and np.array_equal(arrays['selected'],[r['selected'] for r in chunk])):
                raise ValueError('Portable checkpoint arrays disagree with JSON records')
        yield chunk


def load_rows(root):
    manifest = json.loads((Path(root)/'manifest.json').read_text())
    return manifest, [row for chunk in iter_chunks(root, manifest) for row in chunk]


def lease(output, exclusive=False):
    # Never unlink this inode. Workers retain a shared lease while writing.
    handle = (Path(output)/'coordinator.lease').open('a+')
    try:
        fcntl.flock(handle, (fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH) | fcntl.LOCK_NB)
    except OSError:
        handle.close()
        raise RuntimeError('Coordinator or worker lease is live; recovery refused')
    return handle


def locked_chunk(task, token):
    output = Path(task[0]).parent
    with lease(output):
        lock = output/'coordinator.lock'
        if not lock.exists() or json.loads(lock.read_text()).get('token') != token:
            raise RuntimeError('Worker belongs to an expired coordinator; write refused')
        return compute_chunk(task)


def process_namespace():
    # PID liveness is meaningful only in the same Linux PID namespace and boot.
    if sys.platform == 'linux':
        return dict(boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                    pid_namespace=os.readlink('/proc/self/ns/pid'))
    return None


def recover_lock(output):
    """Fail closed for live owners/workers, legacy locks and foreign hosts."""
    output = Path(output).resolve()
    with lease(output, exclusive=True):
        lock = output/'coordinator.lock'
        value = json.loads(lock.read_text())
        if (value.get('lock_version') != 1 or value.get('host') != platform.node()
                or value.get('uid') != os.getuid() or lock.stat().st_uid != os.getuid()
                or value.get('process_namespace') != process_namespace()):
            raise RuntimeError('Legacy, foreign-host, foreign-namespace or foreign-owner lock: recovery refused')
        pid = value.get('pid')
        if not isinstance(pid, int) or isinstance(pid, bool) or pid <= 0:
            raise RuntimeError('Invalid owner PID: recovery refused')
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            pass
        except PermissionError as exc:
            raise RuntimeError('Owner liveness is ambiguous: recovery refused') from exc
        else:
            raise RuntimeError('Owner PID is live (or reused): recovery refused')
        evidence = output/'recovery'
        evidence.mkdir(exist_ok=True)
        archived = evidence/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'-'+uuid.uuid4().hex+'.lock.json')
        # Rename preserves the exact original bytes; no checkpoints are deleted.
        os.rename(lock, archived)
        atomic_json(archived.with_suffix('.audit.json'), dict(recovered_utc=datetime.now(timezone.utc).isoformat(),
                    host=platform.node(),uid=os.getuid(),recovery_pid=os.getpid(),
                    original_sha256=sha(archived),reason='Owner absent and exclusive local worker lease acquired'))
        return str(archived)


def run(output, cells=None, reps=200, phase='pilot', workers=1, chunk_size=50, max_tasks=None, allow_main=False):
    if phase not in PHASES or not isinstance(reps,int) or isinstance(reps,bool):
        raise ValueError('Invalid phase or replication count')
    if phase == 'main':
        if reps != 2000 or not allow_main:
            raise ValueError('Main requires exactly 2000 replications and explicit --allow-main execution opt-in')
    elif not 1 <= reps <= 200:
        raise ValueError('Validation/benchmark/pilot require 1..200 replications')
    if phase == 'pilot' and reps != 200:
        raise ValueError('The pilot has exactly 200 replications per cell')
    if workers < 1 or chunk_size < 1 or (max_tasks is not None and max_tasks < 1):
        raise ValueError('Worker, chunk and task limits must be positive')
    full = catalogue()
    cells = full if cells is None else cells
    if not cells or len({c['id'] for c in cells}) != len(cells):
        raise ValueError('Nonempty distinct catalogue cells required')
    official = {c['id']:c for c in full}
    if any(c['id'] not in official or c != official[c['id']] for c in cells):
        raise ValueError('Only exact approved Experiment A core cells are supported')
    if phase == 'main' and cells != full:
        raise ValueError('Main requires the full ordered 100-cell catalogue')
    output = Path(output).resolve()
    output.mkdir(parents=True,exist_ok=True)
    shards = output/'shards'
    shards.mkdir(exist_ok=True)
    software = provenance()
    identity = dict(protocol='experiment-a-core-v1', catalogue=cells, reps=reps, phase=phase,
                    chunk_size=chunk_size, provenance=software)
    identity_hash = digest(identity)
    lock = output/'coordinator.lock'
    guard = lease(output)
    token = uuid.uuid4().hex
    try:
        fd = os.open(lock,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    except OSError:
        guard.close()
        raise RuntimeError('Coordinator lock exists or cannot be created; use --recover-lock only for a stale owner')
    with os.fdopen(fd,'w') as f:
        f.write(canonical(dict(lock_version=1,process_namespace=process_namespace(),pid=os.getpid(),uid=os.getuid(),host=platform.node(),identity=identity_hash,token=token)))
    old_handler = None
    stopped = False
    def on_term(signum,frame):
        nonlocal stopped
        stopped = True  # Finish only the current bounded wave; leave valid checkpoints.
    try:
        if (output/'manifest.json').exists():
            manifest = json.loads((output/'manifest.json').read_text())
            if manifest['identity_hash'] != identity_hash:
                raise ValueError('Resume identity/code/environment differs from the frozen run')
        else:
            unknown = list(shards.iterdir())
            if unknown:
                raise ValueError('Nonempty shards directory without a run manifest')
            tasks = [dict(cell_id=c['id'],start=i,stop=min(i+chunk_size,reps)) for c in cells for i in range(0,reps,chunk_size)]
            manifest = dict(identity,identity_hash=identity_hash,tasks=tasks,complete=False,sessions=[])
            atomic_json(output/'catalogue.json',cells)
            atomic_json(output/'manifest.json',manifest)
        task_list = [(str(shards),official[t['cell_id']],t['start'],t['stop'],phase) for t in manifest['tasks']]
        done = [verified_done(shards,t) for t in task_list]
        expected_markers = {task_id(t[1],t[2],t[3])+'.done.json' for t in task_list}
        if any(p.name not in expected_markers for p in shards.glob('*.done.json')):
            raise ValueError('Unexpected committed task outside the declared run')
        pending = [t for t,d in zip(task_list,done) if d is None]
        if max_tasks is not None:
            pending = pending[:max_tasks]
        started_utc = datetime.now(timezone.utc).isoformat()
        start_time = time.perf_counter()
        completed_new = []
        old_handler = signal.signal(signal.SIGTERM,on_term)
        if workers == 1:
            for task in pending:
                if stopped:
                    break
                completed_new.append(locked_chunk(task, token))
                print(f'completed {sum(d is not None for d in done)+len(completed_new)}/{len(task_list)} tasks',flush=True)
        else:
            with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
                for start in range(0,len(pending),workers):
                    if stopped:
                        break
                    futures = [pool.submit(locked_chunk,t,token) for t in pending[start:start+workers]]
                    for future in as_completed(futures):
                        completed_new.append(future.result())
                    print(f'completed {sum(d is not None for d in done)+len(completed_new)}/{len(task_list)} tasks',flush=True)
        final = [verified_done(shards,t) for t in task_list]
        manifest['complete'] = all(d is not None for d in final)
        manifest['completed_tasks'] = sum(d is not None for d in final)
        manifest['completed_rows'] = sum(d['rows'] for d in final if d is not None)
        session = dict(started_utc=started_utc,finished_utc=datetime.now(timezone.utc).isoformat(),workers=workers,
                       new_tasks=len(completed_new),new_rows=sum(d['rows'] for d in completed_new),
                       wall_seconds=time.perf_counter()-start_time,interrupted=stopped,
                       coordinator_peak_rss_bytes=rss_bytes())
        manifest['sessions'].append(session)
        atomic_json(output/'manifest.json',manifest)
        return manifest
    finally:
        if old_handler is not None:
            signal.signal(signal.SIGTERM,old_handler)
        if lock.exists() and json.loads(lock.read_text()).get('token') == token:
            lock.unlink()
        guard.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--phase',choices=PHASES,default='pilot')
    p.add_argument('--reps',type=int,default=200)
    p.add_argument('--workers',type=int,default=1)
    p.add_argument('--chunk-size',type=int,default=50)
    p.add_argument('--max-tasks',type=int)
    p.add_argument('--allow-main',action='store_true',help='Explicit execution opt-in; use only after main authorization')
    p.add_argument('--recover-lock',action='store_true',help='Archive a provably stale local lock and exit')
    p.add_argument('--cell-id',action='append')
    args = p.parse_args()
    if args.recover_lock:
        print(recover_lock(args.output))
        return
    cells = catalogue()
    if args.cell_id:
        ids = set(args.cell_id)
        cells = [c for c in cells if c['id'] in ids]
        if {c['id'] for c in cells} != ids:
            p.error('Every cell ID must match the exact catalogue')
    m = run(args.output,cells,args.reps,args.phase,args.workers,args.chunk_size,args.max_tasks,args.allow_main)
    print(canonical(dict(complete=m['complete'],rows=m['completed_rows'],identity=m['identity_hash'])))


if __name__ == '__main__':
    main()
