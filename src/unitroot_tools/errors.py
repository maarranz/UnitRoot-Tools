"""Public numerical errors and transparent specification warnings."""


class EstimationError(RuntimeError):
    """The requested regression/search is numerically infeasible."""


class SpecificationWarning(UserWarning):
    """A default was capped or an estimable design is poorly conditioned."""
