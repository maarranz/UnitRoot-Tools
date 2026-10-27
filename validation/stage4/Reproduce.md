# Reproduce the Stage 4 checks

The checked-in evidence describes macOS arm64 runs; executing these commands
on another host is new evidence, not automatic confirmation of that host.
Runtime downloads, environments and build artifacts in the original run were
kept outside the repository. No global Python links were required.

Use the desired Python 3.11, 3.12, 3.13 or 3.14 interpreter to create a new
virtual environment. Do not enable system site-packages. The tested patch
versions and complete resolved package lists are in `matrix.json` and each
row's `freeze.log`. Released ARCH 8.0.0 is required in every environment.

Current numerical stack:

```sh
python -m pip install arch==8.0.0 numpy==2.4.4 scipy==1.17.1 statsmodels==0.15.0 pandas==3.1.0
```

Exact direct lower-bound stack, exercised on Python 3.11 and 3.12:

```sh
python -m pip install --only-binary=:all: arch==8.0.0 numpy==1.26.4 scipy==1.11.4 statsmodels==0.14.4 pandas==2.2.3
```

The binary-only lower-stack attempts fail on 3.13/3.14 because NumPy 1.26.4 has
no usable wheel. No source-build success is claimed. Test tools were pytest
9.1.1 and pytest-cov 7.1.0; build tooling included setuptools 84.2.1 and build
1.6.1. Ordinary repository installs used the same setuptools build backend
through `uv pip install --no-build-isolation --no-deps /path/to/repo`;
editable installation added `-e`. Artifact installations used actual
`python -m pip install --no-deps --no-build-isolation /path/to/artifact`.
Install the runtime dependencies first when using `--no-deps`.

## Test installed code rather than accidentally importing checkout source

After installation, use a test working directory outside the checkout. From
the repository, with the new environment activated:

```sh
repo_path="$PWD"
check_dir="$(mktemp -d)"
cp -R tests validation examples "$check_dir/"
cat > "$check_dir/pytest.ini" <<'EOF'
[pytest]
testpaths = tests
addopts = --strict-markers
filterwarnings = error
EOF
cd "$check_dir"
python -c 'import unitroot_tools; print(unitroot_tools.__file__)'
python -m pytest -c "$check_dir/pytest.ini" --junitxml=pytest.xml \
  --cov=unitroot_tools --cov-branch --cov-report=json:coverage.json
python validation/audit_saved_results.py --report stage2_audit.json
python validation/validate_stage4.py --report stage3b_replay.json
python validation/validate_stage3c.py
python examples/adf_dfgls.py
python examples/pp_kpss.py
python examples/teaching_walkthrough.py
python -m pip check
```

The copied Stage 3C harness writes only its copy of the Stage 3C numerical
report; original historical artifacts remain intact. An ordinary installed
package should resolve to the environment's site-packages. An editable install
intentionally resolves to checkout source. The original checkout's pytest
configuration includes `pythonpath=src`, so running it unchanged cannot by
itself verify an ordinary installed wheel.

## Build and inspect distributions

From the repository, after installing build/setuptools (and optionally twine):

```sh
python -m build --no-isolation --outdir /path/to/artifacts
python validation/inspect_distributions.py --dist /path/to/artifacts \
  --report /path/to/distributions.json
python -m twine check /path/to/artifacts/*
```

The standard build command builds an sdist and then a wheel from that sdist.
For independent install tests, create two more fresh environments, install the
current dependencies and then install the wheel or sdist. Run the full tests
from the extracted source distribution with the external pytest configuration
above, and check all five functions plus installed documents. The actual Stage
4 artifact test commands and results are in `314-wheel/` and `312-sdist/`.

No network service, repository publication, tag, or new production method is
needed. A Linux rerun and intended Colab workflow should precede broad platform
claims. These scripts use assertions: run without Python `-O`.
