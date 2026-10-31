"""Experiment-A-only throughput/memory/output benchmark; no B or main execution."""
import argparse
from collections import Counter
import json
from pathlib import Path
import time

from .catalogue import catalogue
from .engine import evaluate
from .runner import run,load_rows


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    cells=[]
    for T in (100,250,500):
        group=[c for c in catalogue() if c['law']['T']==T]
        cells.extend([next(c for c in group if c['is_null'] and c['intercept']==0),
                      next(c for c in group if c['is_null'] and c['intercept']==.5),
                      next(c for c in group if not c['is_null'] and c['rho']==.9 and c['intercept']==0)])
    for c in cells:evaluate(c,1000,'benchmark')  # Warm-up is outside timed rows.
    records=[];reference=None
    for workers in (1,2,4):
        output=args.output/f'workers-{workers}'
        start=time.perf_counter()
        manifest=run(output,cells,reps=20,phase='benchmark',workers=workers,chunk_size=5)
        total=time.perf_counter()-start
        _,rows=load_rows(output)
        if reference is None:reference=rows
        assert rows==reference,'Serial/parallel benchmark records differ'
        done=[json.loads(f.read_text()) for f in (output/'shards').glob('*.done.json')]
        worker_peaks={}
        for d in done:worker_peaks[d['process_id']]=max(worker_peaks.get(d['process_id'],0),d['worker_peak_rss_bytes'])
        parent_peak=manifest['sessions'][-1]['coordinator_peak_rss_bytes']
        aggregate=(sum(worker_peaks.values())+parent_peak) if workers>1 else max([parent_peak,*worker_peaks.values()])
        branches=Counter((alpha,variant,len(value['steps']),value['selected'],value['terminal'])
                         for row in rows for alpha,variants in row['enders'].items() for variant,value in variants.items())
        records.append(dict(workers=workers,rows=len(rows),wall_seconds=total,rows_per_second=len(rows)/total,
                            coordinator_peak_rss_bytes=parent_peak,worker_peak_rss_bytes=worker_peaks,
                            summed_task_wall_seconds=sum(d['elapsed_seconds'] for d in done),
                            summed_compute_wall_seconds=sum(d['compute_seconds'] for d in done),
                            summed_worker_cpu_seconds=sum(d['cpu_seconds'] for d in done),
                            sum_process_peaks_upper_estimate_bytes=aggregate,
                            output_bytes=sum(f.stat().st_size for f in output.rglob('*') if f.is_file()),
                            branch_counts=[dict(alpha=k[0],variant=k[1],history_length=k[2],specification=k[3],terminal=k[4],count=v) for k,v in sorted(branches.items())]))
    serial=records[0]['wall_seconds']
    for r in records:
        r['speedup']=serial/r['wall_seconds'];r['parallel_efficiency']=r['speedup']/r['workers']
    report=dict(complete=True,cells=[c['id'] for c in cells],replications_per_cell=20,records=records,
                data_and_decisions_identical_across_workers=True,
                memory_note='Sum of process-specific high-water marks, not simultaneous RSS; shared pages may be counted repeatedly.',
                scope='9 A-core cells, including null zero/high drift and rho=.90 alternatives at each T. No B generation or main experiment.')
    (args.output/'benchmark.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ('workers','wall_seconds','rows_per_second','speedup','parallel_efficiency')} for r in records],indent=2))


if __name__=='__main__':main()
