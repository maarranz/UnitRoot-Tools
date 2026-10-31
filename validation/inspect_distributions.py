"""Inspect local release artifacts without publishing or changing versions."""
import argparse
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path
import tarfile
import tomllib
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dist', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    wheels = list(args.dist.glob('*.whl'))
    sources = list(args.dist.glob('*.tar.gz'))
    assert len(wheels) == len(sources) == 1
    with zipfile.ZipFile(wheels[0]) as wheel:
        names = wheel.namelist()
        assert wheel.testzip() is None
        for module in (repo/'src/unitroot_tools').glob('*.py'):
            assert wheel.read('unitroot_tools/'+module.name) == module.read_bytes()
        metadata_name = next(n for n in names if n.endswith('.dist-info/METADATA'))
        metadata = BytesParser().parsebytes(wheel.read(metadata_name))
        assert metadata['Name'] == 'unitroot-tools'
        assert metadata['Version'] == tomllib.loads((repo/'pyproject.toml').read_text())['project']['version']
        assert metadata['Requires-Python'] == '>=3.11'
        assert metadata['License-Expression'] == 'MIT'
        assert 'arch==8.0.0' in metadata.get_all('Requires-Dist')
        license_name = next(n for n in names if n.endswith('.dist-info/licenses/LICENSE'))
        assert wheel.read(license_name) == (repo/'LICENSE').read_bytes()
        for relative in ['README.md', *[str(p.relative_to(repo)) for p in (repo/'docs').rglob('*.md')],
                         *[str(p.relative_to(repo)) for p in (repo/'examples').glob('*.py')],
                         *[str(p.relative_to(repo)) for folder in ['src/unitroot_tools', 'validation', 'validation/results', 'validation/enders', 'validation/scaled_svd', 'validation/experiment_a/pilot-summary', 'validation/evidence', 'research', 'research/experiment_a', 'research/experiment_a/tests'] for p in (repo/folder).glob('*') if p.is_file() and p.suffix in ('.md', '.py', '.R', '.csv', '.json', '.xml', '.ipynb', '.gz', '.zip') and p.name != 'results.zip']]:
            name = next(n for n in names if n.endswith('/share/doc/unitroot-tools/'+relative))
            assert wheel.read(name) == (repo/relative).read_bytes()
    with tarfile.open(sources[0]) as source:
        members = source.getnames()
        prefix = members[0].split('/')[0]
        required = ['LICENSE', 'README.md', 'MANIFEST.in', 'pyproject.toml',
                    'docs/User_Guide.md', 'docs/archive/Stage4_Release_Review.md',
                    'examples/teaching_walkthrough.py', 'validation/validate_stage4.py',
                    'docs/archive/Stage4B_Release_Candidate.md', 'docs/archive/Release_Notes_0.1.0rc1.md',
                    'validation/Platform_Validation.ipynb', 'validation/run_platform_validation.py',
                    '.github/workflows/platform-validation.yml']
        required += [str(p.relative_to(repo)) for folder in ['src/unitroot_tools','tests','validation/results']
                     for p in (repo/folder).glob('*') if p.is_file() and p.suffix in ('.py','.csv','.json')]
        required += [str(p.relative_to(repo)) for folder in ['docs', 'research', 'validation/enders', 'validation/scaled_svd', 'validation/experiment_a/pilot-summary', 'validation/evidence'] for p in (repo/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix in ('.md', '.py', '.json', '.csv', '.gz', '.zip')]
        for name in required:
            assert source.extractfile(prefix+'/'+name).read() == (repo/name).read_bytes(), name
    report = dict(version=metadata['Version'], requires_python=metadata['Requires-Python'],
                  license_expression=metadata['License-Expression'], wheel_members=names,
                  source_members=members, runtime_modules_byte_identical=True,
                  documentation_license_examples_present=True, reviewed_fixtures_present_in_sdist=True,
                  artifacts={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [wheels[0],sources[0]]})
    args.report.write_text(json.dumps(report,indent=2)+'\n')
    print('Distribution contents, metadata, license and source bytes verified.')


if __name__ == '__main__':
    main()
