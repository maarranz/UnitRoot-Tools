# Validation evidence in Git

The canonical Stage 2 input/result files are in `results/`; their manifests
and reviewed numerical tolerances are unchanged. Historical Stage 3 comparison
CSVs and structured numerical reports are retained. Stage 4/4B structured
matrix, numerical, distribution and preservation records retain tested
versions, results and provenance. Per-environment `dependencies.txt` files
record resolved dependencies without local installation URLs; install the
candidate separately rather than treating these as published-package locks.

Local copied `suite/` sandboxes, verbose logs, coverage details, test XML,
build outputs and the duplicate `results.zip` are excluded from Git. They
remain on the validation machine; generated evidence has not been deleted.
Absolute paths inside retained structured reports are historical provenance,
not portable execution paths or requirements to recreate those directories.
Older reports describe the full local evidence, including excluded logs/XML.

Use `run_platform_validation.py` with a separately installed candidate and an
output directory outside the source tree to regenerate complete evidence.
See `stage4/Reproduce.md` for installed-package testing without source-path
shadowing, and `Platform_Validation.ipynb` / the manually triggered workflow
for pending Colab/Linux validation. Preparing or committing those files does
not establish a platform pass. No tolerance is relaxed on a new platform.

Relocated frozen pilot summaries are in `experiment_a/pilot-summary/`; the
original external pilot remains unchanged. Historical hosted evidence and its
checksums are in [evidence/README.md](evidence/README.md). These older-source
records do not constitute current-code Linux/Colab qualification. See the
[current documentation index](../docs/README.md) for scope and historical reports.
