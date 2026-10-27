"""Immutable scalar diagnostics for common-sample lag selection."""
from dataclasses import dataclass
from numbers import Real
from types import MappingProxyType
from typing import Mapping


@dataclass(frozen=True, slots=True)
class CandidateDiagnostics:
    """One regression on the selection sample, in original squared units."""

    lags: int
    nobs: int
    df_resid: int
    level_coefficient: float
    ssr: float
    variance: float
    level_sum_squares: float
    tau: float
    aic: float
    bic: float
    maic: float
    mbic: float
    criterion_value: float


@dataclass(frozen=True, slots=True)
class LagSelectionResult:
    """Resolved method, sample boundaries and a tuple of frozen candidates.

    Sample boundaries are zero-based response-level indices with exclusive
    stop. The selection sample is not the future final-test regression sample.
    """

    method: str
    selection_detrending: str
    trend: str
    input_nobs: int
    max_lags_requested: int | None
    max_lags: int
    uncapped_rule_value: int | None
    selected_lags: int
    selection_nobs: int
    sample_start: int
    sample_stop: int
    criterion_value: float
    runner_up_gap: float | None
    candidates: tuple[CandidateDiagnostics, ...]
    notes: tuple[str, ...]
    input_name: str | None = None
    criterion_normalization: str = "per_observation_omit_common_constants"


# The maps below are frozen copies; the exposed backend regression remains a
# Statsmodels object and is intentionally not advertised as deeply immutable.


@dataclass(frozen=True, slots=True)
class UnitRootResult:
    """Inference from public released-ARCH tests with explicit sample metadata.

    Decisions use critical values, not ordinary regression p-values. A frozen
    record does not make its exposed Statsmodels regression object immutable.
    """

    test: str
    statistic: float
    pvalue: float
    critical_values: Mapping[float, float]
    trend: str
    input_nobs: int
    nobs: int
    sample_start: int
    sample_stop: int
    lags: int | None
    lag_selection: LagSelectionResult | None
    regression: object | None
    regression_description: str
    calibration: str
    pvalue_provenance: str
    critical_values_provenance: str
    versions: Mapping[str, str]
    null_hypothesis: str
    alternative_hypothesis: str
    notes: tuple[str, ...] = ()
    input_name: str | None = None
    backend: str = "arch"
    backend_version: str = "8.0.0"
    tail: str = "left"
    pvalue_status: str = "approximate"
    pvalue_boundary: str | None = None
    bandwidth: int | None = None
    bandwidth_method: str | None = None
    test_type: str | None = None
    bandwidth_requested: int | None = None
    bandwidth_uncapped: int | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "critical_values", MappingProxyType(dict(self.critical_values)))
        object.__setattr__(self, "versions", MappingProxyType(dict(self.versions)))
        object.__setattr__(self, "notes", tuple(self.notes))
        expected = "right" if self.test == "KPSS" else "left"
        if self.tail != expected:
            raise ValueError(f"{self.test} requires {expected}-tail decisions")

    def reject(self, alpha: float = 0.05) -> bool:
        """Reject in the specified tail at the backend 1%, 5% or 10% critical value.

        Equality does not reject; the approximate p-value never overrides this
        decision. Unsupported alpha values are not interpolated.
        """
        if isinstance(alpha, bool) or not isinstance(alpha, Real):
            raise TypeError("alpha must be a real number: 0.01, 0.05 or 0.10")
        alpha = float(alpha)
        if alpha not in (0.01, 0.05, 0.10):
            raise ValueError("alpha must be 0.01, 0.05 or 0.10")
        critical = self.critical_values[alpha]
        return bool(self.statistic > critical) if self.tail == "right" else bool(self.statistic < critical)

    def summary(self) -> str:
        """Teaching-oriented display, with the nominal 5% decision highlighted."""
        trends = {"n": "no deterministic terms", "c": "intercept", "ct": "intercept and linear trend",
                  "ctt": "intercept, linear and quadratic trends"}
        selection = self.lag_selection
        title = "DF-GLS" if self.test == "DFGLS" else self.test
        kind = "stationarity" if self.test == "KPSS" else "unit-root"
        lines = [f"{title} {kind} test (ARCH {self.backend_version})"]
        if self.input_name is not None:
            lines.append(f"Series: {self.input_name}")
        lines += [f"H0: {self.null_hypothesis}", f"Alternative: {self.alternative_hypothesis}",
                  f"Deterministic terms: {self.trend} ({trends[self.trend]})"]
        if self.bandwidth is None:
            lines.append(f"Augmentation order p: {self.lags}")
        else:
            lines.append(f"Long-run-variance bandwidth b: {self.bandwidth} ({self.bandwidth_method})")
            lines.append(f"Requested bandwidth: {self.bandwidth_requested}; uncapped rule: {self.bandwidth_uncapped}")
        if self.test_type is not None:
            lines.append(f"Statistic type: {self.test_type}")
        lines.append(f"Final observations: {self.nobs} / {self.input_nobs}; response slice [{self.sample_start}:{self.sample_stop}]")
        if self.bandwidth is not None:
            lines.append("Covariance bandwidth; no augmentation-lag selection")
        elif selection is None:
            lines.append("Lag selection: fixed order")
        else:
            lines += [f"Lag selection: {selection.method.upper()}, full-series {selection.selection_detrending} detrending; K={selection.max_lags}",
                      f"Selection observations: {selection.selection_nobs}; response slice [{selection.sample_start}:{selection.sample_stop}]"]
        lines += [f"Statistic: {self.statistic:.6g}; approximate p-value: {self.pvalue:.6g}",
                  "Critical values: " + "; ".join(f"{100*a:g}%: {self.critical_values[a]:.6g}" for a in (0.01, 0.05, 0.10)),
                  "Decision at 5%: " + ("reject H0" if self.reject() else "fail to reject H0"),
                  "Failure to reject does not establish the null.",
                  f"Calibration: {self.calibration}", f"Regression: {self.regression_description}"]
        lines.extend(f"Note: {note}" for note in self.notes)
        return "\n".join(lines)
