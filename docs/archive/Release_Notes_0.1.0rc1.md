# UnitRoot-Tools 0.1.0rc1 — candidate release notes

**Unpublished candidate for instructor review.** This stage changes the version
and release/validation materials, with no econometric or API changes.

- Five teaching functions: `select_lags`, `adf`, `dfgls`, `pp`, `kpss`.
- Independent common-sample Perron–Qu OLS selection; MAIC is the default,
  with distinct MBIC and conventional AIC/BIC options.
- Released ARCH 8.0.0 final estimation and calibration, public regression access
  where available, immutable metadata and critical-value-based decisions.
- Separate PP/KPSS covariance bandwidth handling and clear stationarity versus
  unit-root nulls. No new tests or inference methods are introduced.
- Wheel/source documentation, examples, license, manual Linux workflow and
  unexecuted Colab validation notebook are included for review.

Validated locally on macOS arm64 with Python 3.11–3.14/current stack and
3.11/3.12/direct lower bounds. Exact versions and candidate rerun outcomes are
in `Stage4B_Release_Candidate.md`. Linux and Google Colab are **not verified**;
the supplied workflow/notebook must actually run before claiming those passes.
Exact older bounds lack usable NumPy wheels on Python 3.13/3.14. R evidence is
preserved from the reviewed Stage 2 study, not newly regenerated.

Installation uses the supplied `unitroot_tools-0.1.0rc1-py3-none-any.whl` or
`unitroot_tools-0.1.0rc1.tar.gz`; no public package release is assumed. Match
SHA-256 checksums before using the artifacts. No commit, push, tag or publication
has been authorized by this preparation step.
