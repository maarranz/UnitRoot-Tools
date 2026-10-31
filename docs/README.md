# Documentation index

## Current package and methods

- [User guide](User_Guide.md): APIs, assumptions, diagnostics and interpretation.
- [MAIC/MBIC methodology](MIC_Methodology.md): adopted lag-selection definitions.
- [Enders procedures](Enders.md): Original/Conservative decision trees and caveats.
- [Deterministic AIC/BIC selection](Deterministic_Selection.md): common-sample scores and unadjusted post-selection ADF inference.

ADF uses column-scaled SVD OLS and public Statsmodels MacKinnon calibration.
Enders restricted/auxiliary fits use the stable fitting layer; AIC/BIC reuses
integrated ADF fits. DFGLS, PP and KPSS were not changed by that integration.
The package remains 0.1.0rc1; this index announces no release.

## Current numerical evidence

- [Python 3.13 qualification and numerical investigation](Python313_Numerical_Investigation_Report.md).
- [Calibration-interface investigation](ARCH_Comparison_ADF_Design_Report.md).
- [ADF integration](ADF_Scaled_SVD_Integration_Report.md).
- [Enders integration](Enders_Scaled_SVD_Integration_Report.md).
- [AIC/BIC verification](AIC_BIC_Verification_Report.md).
- [Original cross-package numerical validation](Numerical_Validation.md).

The latest numerical checkpoint recorded 1,615 passing source and installed-wheel
tests on the qualified Mac Python 3.13 stack. Historical matrices do not qualify
new code on other platforms. Exact fixtures and independent validators remain in
[`validation/`](../validation/README.md), including the Enders and scaled-SVD subdirectories.

## Experiment A research

- [Scientific protocol](Monte_Carlo_Research_Protocol.md), preserved as a dated planning record.
- [Research runner and recovery conventions](../research/experiment_a/README.md).
- [Original runner validation](Experiment_A_Validation_Report.md).
- [Original pilot interpretation](Experiment_A_Pilot_Summary.md).
- [Frozen pilot summary files](../validation/experiment_a/pilot-summary/README.md).
- [Theoretical audit](Theoretical_Audit_Sequential_Strategies.md) and [Hamilton/Enders review](Hamilton_Enders_Theoretical_Review.md).

The original pilot and its source identity are immutable historical evidence.
Do not resume it with changed source or present old results as new-code results.

## Historical development

Dated designs, reviews, failed attempts and versioned release notes are preserved
under [`archive/`](archive/README.md). They describe their original checkpoints;
their prospective instructions are not new execution authorization. See that
index for relocation/link provenance and original checksums.

## Outstanding Linux and Colab qualification

Current-code Linux/Colab execution and fresh-runtime persistent-checkpoint
qualification remain separate work. Preserve the [innovation-byte discrepancy](Innovation_Discrepancy_Report.md)
and the historical [hosted failure](archive/Hosted_Final_Report.md) and
[Stage 2 failure](archive/Stage2_Qualification_Report.md). Their unchanged evidence
archives are in [`validation/evidence/`](../validation/evidence/README.md).
The old transfer bundle is archived externally, and the retained root notebook
is an older-source preparation artifact, not current qualification. No main
simulation, platform pass or full-run authorization follows from this cleanup.
