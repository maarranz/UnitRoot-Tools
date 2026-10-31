# Experiment A — Colab qualification, Stage 2

30 October 2026 · **NO-GO for full execution. Gate A fails in the minimally corrected hosted environment.** No checkpoint promotion or fresh-runtime recovery is claimed.

## Outcome

The original warning and metadata failures are diagnosed, with complete original records preserved. An isolated Colab environment corrects the import incompatibility and genuinely installs the unchanged project. Full package testing now collects and executes all 1,367 tests: **1,363 pass; four extreme-scale numerical tests fail**. Independent, Enders and deterministic validators pass. Gate A remains failed, so Stage 2 Gates B–D and the conditional stage/delete/fresh-runtime restore exercise were not executed.

This is a result for the tested environment, not proof that every possible Colab environment must fail. No algorithms, tolerances, expectations or dependency constraints were changed. No full simulation, pilot, Mac integration, Mac dependency installation, repository edit, commit or push occurred.

## Evidence examined before revalidation

Verified prior Hosted_Validation_Evidence.zip SHA256: `5f919ba9d253dc74cc63a71194c3c6997ff72bdf6694f04fabbc85697152946e`. Read Prepared_Validation_Report.md, Hosted_Final_Report.md and the M1 Pro Innovation_Discrepancy_Report.md. Original Linux evidence was never regenerated or edited. The original hosted runtime was still alive: its saved pytest temporary directories were copied before any environment correction, allowing full record comparisons rather than reliance on truncated pytest output.

Complete records, warning lists, manifests, shards and field differences are in the new evidence ZIP. A separate Field_Comparison_Audit.json makes the comparisons reviewable. Comparisons use every field, including hashes, ADF coefficients/statistics, AIC/BIC candidates, Enders paths, decisions, failure/status fields and serialization content. No warning field was removed from a passing qualification comparison.

| Original failing comparison | First difference, after sorting by cell ID and replicate | All record differences |
|---|---|---|
| Order/chunk/resume: 8 vs 8 rows | /0/warnings: 0 vs 14 entries | Only warnings in sorted rows 0, 1, 4 |
| Serial/parallel: 6 vs 6 rows | /3/warnings: 0 vs 14 entries | Only warnings in sorted rows 3, 4 |
| Operational SIGTERM resumed/reference: 80 vs 80 rows | /4/warnings: 0 vs 14 entries | Only warnings in sorted row 4 |
| Portable validation_100, cold vs warm diagnostic | /warnings: 14 vs 0 entries | Only warnings; all other fields exactly equal |
| Original all-grid rows 0 and 4 | Warning lists differ (0 vs 14) | Also ordinary floating differences within the original tolerance; no innovation/input hash or discrete mismatch |

The research and recovery comparisons establish exact equality of all non-warning record fields, not merely numerical tolerance agreement. RNG identity uses the unchanged cell ID, phase, replicate, master seed and key construction. Records do not archive the raw RNG state; no historical raw-state equality is claimed.

Index conventions matter: sorting by cell ID places the two all-grid warning-bearing records at 192 and 376; the original runner load order calls them rows 0 and 4. Both full records in original order are preserved. In research_order the catalogue order/chunk configuration differs intentionally, so manifest identity, task partition and session metadata differ; both manifests are complete. Serial/parallel and operational pairs retain identical run identity, with expected session/timing differences. Manifest differences are archived separately from observation differences.

The original portable failure did not save its actual record. The cold/warm record is explicitly a new two-evaluation diagnostic on the unchanged archived validation_100 vector in the original dependency environment. It does not masquerade as a recovered original record.

## Warning mechanism and correction

Matplotlib 3.10.0 calls legacy pyparsing APIs. With pyparsing 3.3.3, a cold econometric evaluation records 14 PyparsingDeprecationWarning entries on its first adf_n call: oneOf, parseString, resetCache and enablePackrat. A repeated evaluation in the same process records zero. Independent cold Matplotlib import captures 13 warnings with filenames/line numbers in matplotlib/_fontconfig_pattern.py; the complete econometric import adds the packrat initialization warning. Source tracebacks and original records distinguish these compatibility warnings from econometric estimation warnings. No econometric/numerical-warning divergence was observed in the compared records.

Spawned workers import independently, so worker count, task order, chunk partition and a resumed coordinator determine which observation receives initialization warnings. The unchanged package warning-as-error policy instead turns the cold import warning into a collection error. This explains the observed order/process dependence without invoking RNG or fit differences.

Official upstream documentation states that pyparsing 3.3 begins emitting deprecation warnings for these legacy API names: https://pyparsing-docs.readthedocs.io/en/latest/whats_new_in_3_3.html . The targeted correction pins pyparsing 3.2.5 only inside an isolated environment. It does not add warning filters or pre-import modules merely to hide warnings. A strict `-W error` arch/Matplotlib import then succeeds.

A before/after comparison of the complete cold validation_100 record changes only warnings (14 to zero). Its warm record is exactly unchanged, and corrected cold/warm records are exactly equal. This is targeted evidence, not a substitute for Gates B–D or a 400-row revalidation in the corrected environment.

## Reproducible hosted environment and metadata

The isolated environment inherits existing numerical binaries through --system-site-packages. Versions: CPython 3.13.18; Linux x86_64/glibc 2.39; NumPy 2.1.3; SciPy 1.16.3; statsmodels 0.15.0; pandas 2.2.3; ARCH 8.0.0; pytest 8.4.2; Matplotlib 3.10.0; pyparsing 3.2.5. Exact freeze and build configuration are archived. The numerical-build fingerprint remains `c9a799bbd799a8d28f09b26085779899faf43d64b4e01ab5b0cee1b8458db9f2`, identical to Stage 1.

Colab ensurepip failed, so the environment was created with --without-pip and used inherited pip within the isolated prefix. Installation logs confirm the global pyparsing installation was not uninstalled. The unchanged project was copied to a separate build directory, built with --no-deps --no-build-isolation, and installed as a real wheel with --no-deps. No source snapshot or project constraints were edited. Wheel SHA256: `970812d59b1e35a9c1aeff8b17c7796f56ef77227944ca7caf97b57bb166be7d`. Its installed Python source files were checked byte-for-byte against the original source. Genuine distribution metadata reports unitroot-tools 0.1.0rc1, and the Enders validator now completes.

Commands, wheel, PyPI artifact hash/install reports, executed notebook cells, logs and before/after evidence are included. This overlay is a reproducible recipe conditional on the same base numerical builds; it has not been recreated in a fresh hosted runtime and is not certified durable or fully qualified. A diagnostic pip check also reports a pre-existing inherited IPython requirement: jedi is missing. That problem was preserved rather than installing unrelated packages; it is an additional environment limitation.

## Gate A results and stop condition

Hosted source integrity recheck: original 181-file source manifest, 181-file Linux source-reference manifest and 27 diagnostic fixture entries pass. Installed package source matches original bytes; provenance and dependency freeze complete. Package testing reports four failures in `test_small_scale_preserves_valid_fits_scores_and_ranking`, covering AIC/BIC and scales +1e-12/-1e-12. In each case the RSS scale-invariance ratio is **1.0002510614358182**, against the unchanged absolute budget **1e-5** around one. Relative excess is about **2.510614×10^-4**, approximately **25.1 times** the permitted budget.

These are numerical assertion failures after successful collection, not warning or missing-metadata failures. They reproduce the Stage 1 extreme-scale diagnosis with identical numerical binaries. The code reaches the RSS assertion after matching the selected specification; that does not waive the RSS requirement or establish that subsequent assertions passed.

Independent arithmetic/tree validation passes 12 cases at its original 1e-8 tolerance. Enders validation passes with maximum reported auxiliary-t arithmetic error about 2.23e-12; deterministic validation passes its original checks. Their ordinary-scale success does not override the four package failures.

Gate A therefore fails. Gates B–D, corrected-environment 400-row replay, checkpoint staging, runtime termination, genuine fresh-runtime restore and durable deterministic continuation remain **NOT RUN** in Stage 2. Original Stage 1 failures remain preserved. No snapshot was promoted. Changing estimator/scaling/tolerances or treating the four failures as exclusions would exceed this assignment.

## Linux–Mac innovation gate

The Mac report supplies a sufficient mechanism: platform log1p rounding in one Gaussian-tail draw, reproducing both saved innovation hashes from the same underlying random integers. That explains the mismatch but does not satisfy unchanged exact cross-platform equality. Stage 1 Colab matched all 400 original Linux innovation/input hashes; no Stage 2 400-row replay was run after Gate A failed. The strict Linux–Mac gate remains blocked, and the extreme-scale raw-fit issue is separate.

## Deliverables and decision

New Drive directory: `My Drive/UnitRoot-Tools/qualification/prepared-qualification-stage2-01/`. Prior failed qualification remains untouched. `Stage2_Qualification_Evidence.zip` SHA256: **e1d33f3f62d81031bdc8d685cfa5553adf31c61b0ba3a20578dbfdbfd4d99316**. Drive copy read-back matched; the downloaded local ZIP matches that hash, passes CRC and verifies all **370 internal file hashes**. This proves bundle integrity, not runtime-loss recovery.

The browser repeatedly reported JavaScript output-loading errors. Reload restored notebook control; the final evidence was downloaded through Colab's Files pane when the output-frame download did not complete. No cookie/security settings were altered. The hosted runtime was left connected for diagnosis.

**NO-GO for full Colab execution.** Review the four scientific numerical failures and exact cross-platform innovation policy before authorizing any implementation change. Any future environment candidate must preserve the same strict tests and pass Gates A–D before checkpoint promotion. Full Monte Carlo execution always requires separate authorization.
