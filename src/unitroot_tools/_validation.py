"""Input checks shared by the selector and future adapters."""
from numbers import Integral

import numpy as np
from numpy.typing import NDArray


def label(value: str, name: str, choices: tuple[str, ...]) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a string")
    normalized = value.lower()
    if normalized not in choices:
        raise ValueError(f"{name} must be one of {', '.join(choices)}")
    return normalized


def optional_order(value: int | None, name: str) -> int | None:
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an integer or None; booleans are invalid")
    value = int(value)
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def series(y: object) -> tuple[NDArray[np.float64], str | None]:
    # NumPy would silently convert mixed lists such as [1, True, 2] to integers.
    if isinstance(y, (list, tuple)) and any(isinstance(x, (bool, np.bool_)) for x in y):
        raise TypeError("y must contain real numbers, not booleans")
    try:
        a = np.asarray(y)
    except (TypeError, ValueError) as exc:
        raise TypeError("y must be a one-dimensional real numeric series") from exc
    if a.ndim != 1:
        raise TypeError("y must be one-dimensional; column arrays are not accepted")
    if a.size == 0:
        raise ValueError("y must not be empty")
    if a.dtype.kind not in "iuf":
        raise TypeError("y must be real numeric data, not boolean, complex, string or object data")
    with np.errstate(over="ignore", invalid="ignore"):
        a = np.array(a, dtype=np.float64, copy=True)
    if not np.isfinite(a).all():
        raise ValueError("y must be finite after float64 conversion; gaps are not dropped")
    name = getattr(y, "name", None)
    return a, None if name is None else str(name)
