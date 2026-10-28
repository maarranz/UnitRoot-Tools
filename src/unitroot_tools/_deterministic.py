"""Experimental fixed-lag deterministic selection; no post-selection calibration."""
from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from ._adapters import adf
from ._validation import label, optional_order, series
from .errors import EstimationError
from .results import UnitRootResult


_INFERENCE_NOTE = (
    "Fixed-specification ADF p-value reported after selection. "
    "Ordinary ADF calibration ignores deterministic selection and any upstream lag "
    "selection; this is not a selective p-value conditional on winning and does "
    "not establish nominal size control for the complete procedure."
)


@dataclass(frozen=True, slots=True)
class DeterministicCandidate:
    """One unrestricted fixed-lag ADF fit and total Gaussian IC scores."""

    specification: str
    deterministic_dimension: int
    n_coefficients: int
    n_parameters: int
    df_resid: int
    ssr: float
    variance_mle: float
    log_likelihood: float
    aic_score: float
    bic_score: float
    adf_result: UnitRootResult


@dataclass(frozen=True, slots=True)
class DeterministicSelectionResult:
    """Frozen scalar record; nested backend regression objects remain mutable.

    Response boundaries are zero-based level indices with exclusive stop.
    No selector-level p-value or unit-root decision is provided.
    """

    criterion: str
    selected_specification: str
    lags: int
    input_nobs: int
    nobs: int
    sample_start: int
    sample_stop: int
    candidates: tuple[DeterministicCandidate, ...]
    versions: Mapping[str, str]
    notes: tuple[str, ...]
    input_name: str | None = None
    lag_policy: str = "supplied_common_integer"
    criterion_normalization: str = "total_omit_common_gaussian_constant_include_variance"
    tie_break: str = "exact_score_then_smaller_deterministic_dimension"
    inference_status: str = "fixed_specification_calibration_after_selection"
    experimental: bool = True

    def __post_init__(self) -> None:
        object.__setattr__(self, "versions", MappingProxyType(dict(self.versions)))
        object.__setattr__(self, "notes", tuple(self.notes))
        object.__setattr__(self, "candidates", tuple(self.candidates))

    @property
    def selected_candidate(self) -> DeterministicCandidate:
        return next(c for c in self.candidates if c.specification == self.selected_specification)

    @property
    def selected_fixed_spec_adf(self) -> UnitRootResult:
        """Existing ADF result, with fixed-specification calibration after selection."""
        return self.selected_candidate.adf_result

    @property
    def criterion_value(self) -> float:
        return getattr(self.selected_candidate, self.criterion + "_score")

    @property
    def score_differences(self) -> Mapping[str, float]:
        """Requested-criterion scores minus the winning score, without rounding."""
        return MappingProxyType({c.specification: getattr(c, self.criterion + "_score")
                                 - self.criterion_value for c in self.candidates})

    @property
    def runner_up_gap(self) -> float:
        ordered = sorted(getattr(c, self.criterion + "_score") for c in self.candidates)
        return ordered[1] - ordered[0]

    def summary(self) -> str:
        """Display selection and qualified fixed-specification inference together."""
        selected = self.selected_fixed_spec_adf
        lines = [f"Experimental deterministic selection: {self.criterion.upper()}, fixed p={self.lags}"]
        if self.input_name is not None:
            lines.append(f"Series: {self.input_name}")
        lines += [f"Common observations: {self.nobs} / {self.input_nobs}; response slice [{self.sample_start}:{self.sample_stop}]",
                  "Scores omit the common Gaussian constant and include the variance parameter.",
                  "spec  coefficients  parameters  AIC  BIC  selected-criterion difference"]
        lines += [f"{c.specification:4}  {c.n_coefficients}  {c.n_parameters}  {c.aic_score:.12g}  {c.bic_score:.12g}  {self.score_differences[c.specification]:.12g}"
                  for c in self.candidates]
        lines += [f"Selected: {self.selected_specification}; runner-up gap: {self.runner_up_gap:.12g}",
                  "Fixed-specification ADF calibration for the selected model",
                  f"ADF statistic: {selected.statistic:.12g}",
                  f"Fixed-specification ADF p-value reported after selection: {selected.pvalue:.12g}",
                  _INFERENCE_NOTE,
                  "Critical values: " + "; ".join(f"{100*a:g}%: {selected.critical_values[a]:.12g}" for a in (0.01, 0.05, 0.10)),
                  f"Calibration: {selected.calibration}",
                  f"P-value provenance: {selected.pvalue_provenance}",
                  f"Critical-value provenance: {selected.critical_values_provenance}"]
        lines.extend(f"Note: {note}" for note in self.notes if note != _INFERENCE_NOTE)
        return "\n".join(lines)


def _scores(ssr: float, nobs: int, n_parameters: int) -> tuple[float, float, float, float]:
    variance = ssr / nobs
    if not math.isfinite(variance) or variance <= 0:
        raise EstimationError("positive finite ML residual variance is required")
    base = nobs * math.log(variance)
    llf = -0.5 * (base + nobs * (1 + math.log(2 * math.pi)))
    aic = base + 2 * n_parameters
    bic = base + math.log(nobs) * n_parameters
    if not all(math.isfinite(v) for v in (llf, aic, bic)):
        raise EstimationError("nonfinite Gaussian likelihood or information criterion")
    return variance, llf, aic, bic


def _winner(candidates: tuple[DeterministicCandidate, ...], criterion: str) -> DeterministicCandidate:
    return min(candidates, key=lambda c: (getattr(c, criterion + "_score"), c.deterministic_dimension))


def _check_alignment(result: UnitRootResult, data: np.ndarray, p: int, spec: str) -> None:
    """Check actual response and named columns against the common raw sample."""
    n = len(data) - 1 - p
    d = {"n": 0, "c": 1, "ct": 2}[spec]
    reg = result.regression
    if (result.test != "ADF" or result.trend != spec or result.lags != p
            or result.lag_selection is not None or result.input_nobs != len(data)
            or result.nobs != n or result.sample_start != p + 1 or result.sample_stop != len(data)
            or reg is None):
        raise EstimationError("fixed-lag common-sample metadata mismatch")
    dy = np.diff(data)
    columns = [data[p:-1]] + [dy[p-j:len(data)-1-j] for j in range(1, p+1)]
    names = ["Level.L1"] + [f"Diff.L{j}" for j in range(1, p+1)]
    if d:
        columns.append(np.ones(n))
        names.append("const")
    if d == 2:
        columns.append(np.arange(1, n+1, dtype=float))
        names.append("trend")
    expected = np.column_stack(columns)
    # The ADF adapter already checks rank using normalized columns. Rechecking
    # unscaled columns would make acceptance depend on the series units.
    if (list(reg.model.exog_names) != names or not np.array_equal(reg.model.endog, dy[p:])
            or not np.array_equal(reg.model.exog, expected)
            or reg.nobs != n or reg.df_resid != n - expected.shape[1]):
        raise EstimationError("common-sample response or regressor alignment/rank mismatch")
    if not all(math.isfinite(float(v)) for v in (result.statistic, result.pvalue, *result.critical_values.values())):
        raise EstimationError("nonfinite fixed-specification ADF inference")


def select_deterministic(y: object, *, lags: int, criterion: str = "bic") -> DeterministicSelectionResult:
    """Select n/c/ct by total Gaussian AIC/BIC at one supplied common lag order.

    Experimental specification selection, not a calibrated selection-and-testing
    procedure. Fixed-specification ADF p-value reported after selection.
    The input must be equally spaced; calendar spacing is the caller's responsibility.
    """
    criterion = label(criterion, "criterion", ("aic", "bic"))
    if lags is None:
        raise TypeError("lags must be a nonnegative integer; None is invalid")
    p = optional_order(lags, "lags")
    data, input_name = series(y)
    if len(data) < 2 * p + 5:
        raise ValueError("Deterministic candidate ct: insufficient sample for positive residual degrees of freedom")
    candidates = []
    for d, spec in enumerate(("n", "c", "ct")):
        try:
            result = adf(data, trend=spec, lags=p)
            _check_alignment(result, data, p, spec)
            q = 1 + p + d
            ssr = float(result.regression.ssr)
            variance, llf, aic, bic = _scores(ssr, result.nobs, q + 1)
            candidates.append(DeterministicCandidate(spec, d, q, q+1, result.nobs-q,
                                                     ssr, variance, llf, aic, bic, result))
        except (TypeError, ValueError, EstimationError) as exc:
            raise type(exc)(f"Deterministic candidate {spec}: {exc}") from exc
    candidates = tuple(candidates)
    winner = _winner(candidates, criterion)
    return DeterministicSelectionResult(
        criterion=criterion, selected_specification=winner.specification, lags=p,
        input_nobs=len(data), nobs=len(data)-1-p, sample_start=p+1, sample_stop=len(data),
        candidates=candidates, versions=winner.adf_result.versions, input_name=input_name,
        notes=(_INFERENCE_NOTE,
               "Including n may be inappropriate for economic series. Selection does not prove "
               "the presence or absence of drift/trend; use prespecified ADF when theory gives strong guidance.",
               "Equally spaced observations are required; calendar spacing is not checked."),
    )
