# Deterministic AIC/BIC rank correction, validation and commit

28 October 2026 · UnitRoot-Tools · MacBook Pro M1 Pro

The blocker is resolved. Commit `b6a84fec7ca28a1802cdd722c163f624d16dca24` was pushed successfully to existing
`origin/main` by an ordinary fast-forward push from `8c4f1a1`. Local HEAD and
origin/main agree. Version remains 0.1.0rc1; no tag or release was created.

## Exact correction and documentation

Removed only `np.linalg.matrix_rank(expected) != expected.shape[1]` from the
selector's `_check_alignment`. Added a comment explaining reliance on the
unchanged ADF adapter's normalized full-rank safeguard. All separate response,
column-value/name, sample, lag-order and residual-df checks remain. No existing
ADF adapter code, numerical tolerance, Enders code, dependency or package metadata
changed. Existing invalid-rank/degenerate-data tests still pass.

The approved design document now identifies its historical status, links to
`Deterministic_Selection.md`, removes machine-specific links, corrects obsolete
implementation assertions and uses “Fixed-specification ADF p-value reported
after selection” in its example. Its theoretical design is preserved.

Four regression cases were added: positive/negative 1e−12 scaling crossed with
AIC/BIC on the documented 120-observation fixture at p=1. They verify all three
ADF fits, normalized rank/conditioning, successful selection, unchanged winners,
full AIC/BIC ordering, the common 2n log|s| score shift, and runner-up gaps.
No existing test or tolerance was weakened.

## Validation

| Path | Suite | Passed | Failures | Errors | Skipped |
|---|---|---:|---:|---:|---:|
| source | full | 1367 | 0 | 0 | 0 |
| source | focused | 80 | 0 | 0 | 0 |
| source | enders | 123 | 0 | 0 | 0 |
| wheel | full | 1367 | 0 | 0 | 0 |
| wheel | focused | 80 | 0 | 0 | 0 |
| wheel | enders | 123 | 0 | 0 | 0 |
| sdist | full | 1367 | 0 | 0 | 0 |
| sdist | focused | 80 | 0 | 0 | 0 |
| sdist | enders | 123 | 0 | 0 | 0 |

All final test runs emitted zero uncaptured warnings, with warnings configured
as errors. The complete suite is the original 1,287 tests plus 80 selector tests.
Focused deterministic tests and independent validation passed, followed by the
complete suite and focused Enders checks on each path. Fresh wheel and sdist were
built externally from HEAD plus exactly the nine approved paths; installed checks
proved their imports came from isolated target directories with identical runtime
bytes. Both installed dependency checks and distribution inspection passed.
The current documentation example executed successfully on all three paths.

The standalone independent deterministic validator retained **1e−9 absolute**
tolerances. Its maximum observed ordinary-scale discrepancy remains
1.4210854715202e-13; reports agree exactly across
source/wheel/sdist. It covers 18 distinct fits and 12 selection calls per path.

A separate independent normalized-QR diagnostic confirms the small-scale fix.
The normalized ct design has rank 4 and condition approximately 3.82. The
unchanged backend's raw-unit regression produces a maximum relative RSS error
of 8.647838938280472e-06 against normalized QR. Consequently its
maximum score-shift discrepancy is 0.00102044058348838,
which is distinct from the ordinary-scale validator result above.

The new extreme-scale test explicitly allows a relative RSS error of 1e−5.
For positive RSS this implies an absolute score bound of
−n log(1−1e−5)+1e−9 = 0.001180006900039334 at n=118,
and twice that bound for differences/gaps. The independent normalized-QR check
verifies this RSS budget directly. These tolerances describe existing backend
roundoff under extreme units; the prior 1e−9 ordinary-scale tests and the ADF
adapter's rank tolerances were unchanged. Winners and rankings match exactly.
All three distributions give identical small-scale diagnostic reports.

The independent Enders validator and additional 84-case variant-equivalence
checks were repeated on all three paths. Their complete JSON reports match the
saved earlier Mac baseline exactly. The published GDP example therefore remains:
Original rejects, default Conservative does not, and permitted Conservative agrees
with Original. The previously recorded printed A2.3 intercept t discrepancy
(2.646058967 versus 2.64, difference about .00606 exceeding .005 rounding)
remains documented; no tolerance was changed to conceal it.

Nonfatal diagnostics: the existing manifest finds no validation XML files;
byte-compilation is intentionally disabled in external builds; NumPy configuration
display requests optional pyyaml; system_profiler reports an hw.cpufamily lookup
error while reporting hardware successfully. There are no remaining validation
failures. These diagnostics are separate from test warnings.

## Commit scope and preservation

Exact committed paths:

```text
README.md
docs/Deterministic_AIC_BIC_Design.md
docs/Deterministic_Selection.md
docs/User_Guide.md
src/unitroot_tools/__init__.py
src/unitroot_tools/_deterministic.py
tests/test_deterministic.py
tests/test_selection.py
validation/validate_deterministic.py
```

The staged path set and every staged byte were checked against the validated
snapshot; `git diff --cached --check` passed before commit. The only files changed
during this correction were `_deterministic.py`, `test_deterministic.py`, and the
approved design document. Before/after checks covered 724 pre-existing non-Git,
non-bytecode files. All unrelated material and the earlier implementation files
outside those three paths were preserved. No historical ZIP/report or broader
audit was staged. Git's separately protected index initially blocked staging;
after the requested internal write access was granted, staging and commit passed.

Commit message: `Add experimental deterministic AIC/BIC selection for ADF`.
Push result: `8c4f1a1..b6a84fe main -> main` to
`git@github.com:maarranz/UnitRoot-Tools.git`. No force push.

Final Git status:

```text
## main...origin/main
?? docs/Deterministic_AIC_BIC_Final_Review_2026-10-28.md
?? docs/Deterministic_AIC_BIC_Implementation_Validation_Report_2026-10-28.md
?? docs/Enders_Mac_M1_Pro_Validation_Report_2026-10-28.md
?? docs/Theoretical_Audit_Sequential_Strategies.md
```

The tracked working tree and index are clean, with four intentionally untracked
report/audit documents retained. There is no remaining blocking finding.

## Reproduction and limitations

The adjacent evidence archive contains exact commands/environment overrides,
source/wheel/sdist logs and JUnit, independent JSON reports, baseline equivalence
checks, staged patch, snapshot/preservation hashes, distribution inspection,
validation scripts and fresh artifacts. Hardware identifiers are redacted.
The interpreter is the previously validated Python 3.12.16 arm64 environment:
`/Users/maarranz/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/stage4b/env312-current/bin/python`.
Dependencies remain ARCH 8.0.0, NumPy 2.4.4, SciPy 1.17.1, Statsmodels 0.15.0,
pandas 3.1.0 and pytest 9.1.1 on macOS 26.7.1/M1 Pro.

Artifacts are validation-only builds under the unchanged version:

```json
{
  "unitroot_tools-0.1.0rc1-py3-none-any.whl": {
    "bytes": 127233,
    "sha256": "992b1efc80ba1eed803eea4450ed148d0c6b926155a634525c74bc6500c73764"
  },
  "unitroot_tools-0.1.0rc1.tar.gz": {
    "bytes": 280518,
    "sha256": "dd2234f3e91fb9baddf0f256117c1e8141acd9c766501ce9da7398ac009a8186"
  }
}
```

The feature remains experimental deterministic specification selection. Its
fixed-specification ADF inference does not account for deterministic or upstream
lag selection and does not establish overall nominal size. Backend rounding can
affect almost-ties, particularly under extreme measurement units. No Monte Carlo,
bootstrap, automatic-lag extension or new strategy was implemented or run.
