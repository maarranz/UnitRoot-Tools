# Deterministic AIC/BIC final review — stopped before commit

28 October 2026 · UnitRoot-Tools · MacBook Pro M1 Pro

## Outcome

One blocking numerical issue was confirmed. Following the assignment's stop rule,
no corrections, staging, commit or push were performed. Existing files and unrelated
untracked material were preserved.

## Repository and baseline

- Branch: main.
- Local HEAD and live origin/main: `8c4f1a17cf12207dba3c051be6f1de7943df895c`.
- All eight implementation paths in the prior validation snapshot match their saved SHA-256 hashes exactly.
- The repository's untracked implementation report matches the delivered report byte-for-byte.
- The approved design is additionally untracked and listed for the proposed nine-file commit.
- The index was empty and remains empty. Four tracked implementation files are modified; five approved paths remain untracked. Three unrelated untracked reports/audit documents also remain untracked.

## Blocking finding: scale-dependent duplicate rank check

Location: `src/unitroot_tools/_deterministic.py:158`, in `_check_alignment`.

After the existing ADF adapter has validated full rank using normalized columns,
the selector checks `np.linalg.matrix_rank(expected)` on unscaled columns. With
small measurement units, the level/difference columns become tiny relative to
the unchanged constant/trend columns. NumPy's raw absolute rank threshold then
classifies a well-conditioned, column-normalized design as deficient.

Reproducer: use the documented fixture

```python
import numpy as np
from unitroot_tools import adf, select_deterministic

t = np.arange(120, dtype=float)
y = .1 + np.sin(.8*t) + .6*np.cos(1.71*t) + .15*np.sin(.17*t*t)
select_deterministic(y, lags=1)  # selects n
for spec in ("n", "c", "ct"):
    adf(1e-12*y, trend=spec, lags=1)  # all succeed
select_deterministic(1e-12*y, lags=1)  # raises for ct
```

Observed ct design: 118 observations, four columns, residual df 114. The existing
adapter accepts it. Raw numerical rank is 2, while normalized rank is 4 and the
normalized condition number is approximately 3.8218463344. RSS is positive and
finite, approximately 3.9974846954e-23. The selector raises:

```text
EstimationError: Deterministic candidate ct: common-sample response or regressor alignment/rank mismatch
```

This creates an unnecessary rejection after valid fixed-lag fits and undermines
the intended scale-equivariance convention. The existing scale tests cover
−3, .001 and 100 and do not expose this case. The earlier 1,363 passing tests
remain valid evidence for their tested inputs, but do not resolve this finding.

Recommended smallest correction, subject to discussion: rely on the adapter's
existing normalized full-rank safeguard while retaining exact alignment and df
checks, or normalize columns before the additional rank check. Add a focused
small-scale regression test that confirms successful selection and compares
independent scoring/routing at an appropriately stated tolerance. Do not change
existing econometric adapters or relax rank tolerances to mask the issue.

## Nonblocking documentation observations

The implementation, README, User Guide and current methodology page correctly
describe experimental specification selection, use “Fixed-specification ADF
p-value reported after selection,” and distinguish the fixed calibration from
selective inference and overall size control. External lag-selection uncertainty
and economic prespecification are explicit.

The approved design document still says “design only” and “No procedure, tests or
package changes have been implemented,” contains prospective API wording and an
older p-value example label, and links to local output files. Before committing
it as current project documentation, clarify its historical/approved-design status,
point to the implementation page, and align its example terminology. These are
small documentation corrections; none were made after finding the blocker.

## Other source-review conclusions

The reviewed fixed-p response/lag alignment, parameter counts including variance,
RSS/n likelihood, total reduced IC conventions, global unrounded minimization,
exact dimension-based tie rule, all-or-error behavior, three-fit/no-refit flow,
nested UnitRootResult compatibility and absence of selector-level unit-root
inference agree with the approved contract. Frozen scalars/read-only maps and
mutable nested backend objects follow the documented conventions. No stylistic
refactor is recommended.

## Checks performed

Repository/status/remotes inspection; exact snapshot/report byte comparison;
source, exports, tests, independent validator and documentation review; live
`git ls-remote origin refs/heads/main`; targeted numerical reproducer with raw
and normalized ranks/condition numbers. No full suite or distribution rebuild
was repeated. No package dependencies or environment settings were modified.

The reproducer used the previously validated Python 3.12.16 environment with
`PYTHONPATH=/Users/maarranz/GitHub/UnitRoot-Tools/src`,
`PYTHONDONTWRITEBYTECODE=1` and `PYTHONNOUSERSITE=1`.

## Commit, push and status

Files committed: none. New commit hash: none. Push: not attempted. No force push,
tag or publication. Remote main had not moved at inspection. Final tracked and
untracked status is unchanged from the start of this review.

Next step: discuss/authorize resolution of the numerical blocker, then make the
small correction and affected tests, update the design's historical wording,
review the exact nine-file staged diff and resume the authorized commit/push
only when the blocking finding is resolved.
