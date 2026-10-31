"""Small validation-only SIGTERM/SIGKILL recovery exercise, never main or pilot."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from .catalogue import catalogue
from .runner import load_rows, recover_lock, run


def wait_for(predicate, timeout=30):
    deadline=time.monotonic()+timeout
    while not predicate():
        if time.monotonic()>deadline:
            raise TimeoutError('Operational validation timed out')
        time.sleep(.02)


def validate(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    cells=[c for c in catalogue() if c['law']['T']==500 and c['is_null'] and c['intercept'] in (0,.5)]
    reference=output/'reference'
    run(reference,cells,40,'validation',2,2)
    reference_rows=load_rows(reference)[1]
    results=[]
    for name,sig in [('sigterm',signal.SIGTERM),('sigkill',signal.SIGKILL)]:
        target=output/name
        command=[sys.executable,'-m','research.experiment_a.runner','--output',str(target),
                 '--phase','validation','--reps','40','--workers','2','--chunk-size','2']
        for cell in cells:command.extend(['--cell-id',cell['id']])
        with (output/(name+'.log')).open('w') as log:
            process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
            try:
                wait_for(lambda: bool(list((target/'shards').glob('*.done.json'))) or process.poll() is not None)
                assert process.poll() is None, 'Run finished before interruption test'
                try:
                    recover_lock(target)
                except RuntimeError:
                    pass
                else:
                    raise AssertionError('Live coordinator recovery incorrectly accepted')
                process.send_signal(sig)
                code=process.wait(timeout=30)
            finally:
                if process.poll() is None:
                    process.kill();process.wait()
        recovered=None
        if sig==signal.SIGKILL:
            deadline=time.monotonic()+30
            while True:
                try:
                    recovered=recover_lock(target)
                    break
                except RuntimeError as exc:
                    if 'lease is live' not in str(exc) or time.monotonic()>deadline:raise
                    time.sleep(.02)
            assert code==-signal.SIGKILL
        else:
            assert code==0 and not (target/'coordinator.lock').exists()
            assert json.loads((target/'manifest.json').read_text())['sessions'][-1]['interrupted']
        committed=len(list((target/'shards').glob('*.done.json')))
        assert 0<committed<40
        # Worker-count change on resume is deliberately exercised.
        run(target,cells,40,'validation',1,2)
        assert load_rows(target)[1]==reference_rows
        results.append(dict(signal=name,committed_before_resume=committed,recovery_archive=recovered,
                            rows=80,exact_record_parity=True,live_recovery_refused=True))
    report=dict(all_passed=True,phase='validation',results=results)
    (output/'operations.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if not __debug__:raise RuntimeError('Assertions required')
    print(json.dumps(validate(args.output),indent=2))


if __name__=='__main__':main()
