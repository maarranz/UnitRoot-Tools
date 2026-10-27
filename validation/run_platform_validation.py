"""Validate an installed candidate without source-path shadowing.

Run using the environment containing the installed candidate and test tools:
python validation/run_platform_validation.py --source /path/to/source --out /tmp/rc-results
No original results are overwritten. Linux/Colab status requires actual execution there.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--expected-version', default='0.1.0rc1')
    args = parser.parse_args()
    source, out = args.source.resolve(), args.out.resolve()
    if out == source or source in out.parents:
        parser.error('--out must be outside the source tree')
    out.mkdir(parents=True, exist_ok=True)
    suite = out/'suite'
    suite.mkdir(exist_ok=True)
    fixtures = source/'validation/results'
    before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in fixtures.iterdir() if p.is_file()}
    for name in ('tests', 'validation', 'examples'):
        shutil.copytree(source/name, suite/name, dirs_exist_ok=True,
                        ignore=shutil.ignore_patterns('__pycache__', 'stage4', 'stage4b'))
    config = suite/'pytest.ini'
    config.write_text('[pytest]\ntestpaths = tests\naddopts = --strict-markers\nfilterwarnings = error\n')
    env = dict(os.environ, PYTHONNOUSERSITE='1')
    env.pop('PYTHONPATH', None)
    smoke = '''import importlib.metadata as m,json,pathlib,platform,sys,numpy as np,unitroot_tools as u
expected=sys.argv[1]
assert u.__version__==m.version('unitroot-tools')==expected
assert m.version('arch')=='8.0.0'
assert 'site-packages' in u.__file__,u.__file__
y=np.random.default_rng(20261027).normal(size=240).cumsum()
for f in (u.adf,u.dfgls,u.pp,u.kpss):
 r=f(y);assert r.versions['arch']=='8.0.0';assert np.isfinite(r.statistic)
assert u.select_lags(y,max_lags=12).method=='maic'
print(json.dumps(dict(python=sys.version,platform=platform.platform(),module=u.__file__,versions={x:m.version(x) for x in ('unitroot-tools','arch','numpy','scipy','statsmodels','pandas','pytest','pytest-cov')},venv_config=pathlib.Path(sys.prefix,'pyvenv.cfg').read_text() if pathlib.Path(sys.prefix,'pyvenv.cfg').exists() else None),indent=2))
'''
    commands = [
        ('installed-smoke', [sys.executable, '-c', smoke, args.expected_version]),
        ('freeze', [sys.executable, '-m', 'pip', 'freeze']),
        ('pip-check', [sys.executable, '-m', 'pip', 'check']),
        ('pytest', [sys.executable, '-m', 'pytest', '-c', str(config), '--junitxml='+str(out/'pytest.xml'),
                    '--cov=unitroot_tools', '--cov-branch', '--cov-report=json:'+str(out/'coverage.json')]),
        ('stage2-audit', [sys.executable, str(suite/'validation/audit_saved_results.py'), '--report', str(out/'stage2_audit.json')]),
        ('stage3b-replay', [sys.executable, str(suite/'validation/validate_stage4.py'), '--report', str(out/'stage3b_replay.json')]),
        ('stage3c-replay', [sys.executable, str(suite/'validation/validate_stage3c.py')]),
        *[(name, [sys.executable, str(suite/'examples'/name)]) for name in ('adf_dfgls.py','pp_kpss.py','teaching_walkthrough.py')],
        ('numpy-config', [sys.executable, '-c', 'import numpy; numpy.show_config()']),
    ]
    results = {}
    for name, command in commands:
        with (out/(name+'.log')).open('w') as log:
            proc = subprocess.run(command, cwd=suite, env=env, stdout=log, stderr=subprocess.STDOUT)
        results[name] = dict(returncode=proc.returncode, command=command)
    for name in ('stage3c_comparisons.csv','stage3c_numerics.json'):
        p = suite/'validation'/name
        if p.exists():
            shutil.copy2(p, out/name)
    after = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in fixtures.iterdir() if p.is_file()}
    report = dict(expected_version=args.expected_version, python=platform.python_version(), platform=platform.platform(),
                  commands=results, original_fixtures_unchanged=before == after, fixture_sha256=after,
                  all_passed=before == after and all(r['returncode']==0 for r in results.values()))
    (out/'status.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({key: report[key] for key in ('platform','python','all_passed','original_fixtures_unchanged')}))
    return 0 if report['all_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
