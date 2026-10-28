"""Atomic, resumable Experiment A reference runner. Main experiment is disabled."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
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
    return dict(rng=RNG_VERSION, master_seed_hex=f'{MASTER_SEED:032x}',
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
                compute_seconds=computed-start_time, elapsed_seconds=time.perf_counter()-start_time,
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


def load_rows(root):
    root = Path(root)
    manifest = json.loads((root/'manifest.json').read_text())
    rows = []
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
        rows.extend(chunk)
    return manifest, rows


def run(output, cells=None, reps=200, phase='pilot', workers=1, chunk_size=50, max_tasks=None):
    if phase not in PHASES or not isinstance(reps,int) or isinstance(reps,bool) or not 1 <= reps <= 200:
        raise ValueError('Only validation/benchmark/pilot with 1..200 replications is authorized')
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
    output = Path(output).resolve()
    output.mkdir(parents=True,exist_ok=True)
    shards = output/'shards'
    shards.mkdir(exist_ok=True)
    software = provenance()
    identity = dict(protocol='experiment-a-core-v1', catalogue=cells, reps=reps, phase=phase,
                    chunk_size=chunk_size, provenance=software)
    identity_hash = digest(identity)
    lock = output/'coordinator.lock'
    try:
        fd = os.open(lock,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    except FileExistsError as exc:
        raise RuntimeError('Coordinator lock exists; confirm old coordinator/workers stopped before explicit manual recovery') from exc
    with os.fdopen(fd,'w') as f:
        f.write(canonical(dict(pid=os.getpid(),host=platform.node(),identity=identity_hash)))
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
        start_time = time.perf_counter()
        completed_new = []
        old_handler = signal.signal(signal.SIGTERM,on_term)
        if workers == 1:
            for task in pending:
                if stopped:
                    break
                completed_new.append(compute_chunk(task))
                print(f'completed {sum(d is not None for d in done)+len(completed_new)}/{len(task_list)} tasks',flush=True)
        else:
            with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
                for start in range(0,len(pending),workers):
                    if stopped:
                        break
                    futures = [pool.submit(compute_chunk,t) for t in pending[start:start+workers]]
                    for future in as_completed(futures):
                        completed_new.append(future.result())
                    print(f'completed {sum(d is not None for d in done)+len(completed_new)}/{len(task_list)} tasks',flush=True)
        final = [verified_done(shards,t) for t in task_list]
        manifest['complete'] = all(d is not None for d in final)
        manifest['completed_tasks'] = sum(d is not None for d in final)
        manifest['completed_rows'] = sum(d['rows'] for d in final if d is not None)
        session = dict(started_utc=datetime.now(timezone.utc).isoformat(),workers=workers,
                       new_tasks=len(completed_new),new_rows=sum(d['rows'] for d in completed_new),
                       wall_seconds=time.perf_counter()-start_time,interrupted=stopped,
                       coordinator_peak_rss_bytes=rss_bytes())
        manifest['sessions'].append(session)
        atomic_json(output/'manifest.json',manifest)
        return manifest
    finally:
        if old_handler is not None:
            signal.signal(signal.SIGTERM,old_handler)
        lock.unlink(missing_ok=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--phase',choices=PHASES,default='pilot')
    p.add_argument('--reps',type=int,default=200)
    p.add_argument('--workers',type=int,default=1)
    p.add_argument('--chunk-size',type=int,default=50)
    p.add_argument('--max-tasks',type=int)
    p.add_argument('--cell-id',action='append')
    args = p.parse_args()
    cells = catalogue()
    if args.cell_id:
        ids = set(args.cell_id)
        cells = [c for c in cells if c['id'] in ids]
        if {c['id'] for c in cells} != ids:
            p.error('Every cell ID must match the exact catalogue')
    m = run(args.output,cells,args.reps,args.phase,args.workers,args.chunk_size,args.max_tasks)
    print(canonical(dict(complete=m['complete'],rows=m['completed_rows'],identity=m['identity_hash'])))


if __name__ == '__main__':
    main()
