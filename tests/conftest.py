"""Read reviewed fixtures without importing the research harness or any backend."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "validation" / "results"


@pytest.fixture(scope="session")
def stage2():
    manifest = json.loads((RESULTS / "manifest.json").read_text())
    # Fail visibly if the reviewed input/result fixture provenance changes.
    for name in ("series.csv", "python_development.csv", "r_results.csv"):
        assert hashlib.sha256((RESULTS / name).read_bytes()).hexdigest() == manifest["file_sha256"][name]
    def read(name):
        with (RESULTS / name).open(newline="") as f:
            return list(csv.DictReader(f))
    columns = read("series.csv")
    data = {}
    for name in columns[0]:
        values = [r[name] for r in columns]
        length = manifest["lengths"][name]
        assert all(values[:length]) and not any(values[length:])
        data[name] = np.array(values[:length], dtype=float)
    candidates = read("mic_candidates.csv")
    assert hashlib.sha256((RESULTS / "mic_candidates.csv").read_bytes()).hexdigest() == "b50fe30d96de6c3b5688a8347ceb9bf9de872a6d47163736575c214081bbd1e7"
    return data, candidates, read("python_development.csv"), read("r_results.csv")
