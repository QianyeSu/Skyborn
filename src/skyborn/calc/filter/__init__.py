"""
Lanczos filter module for Skyborn.

Provides Lanczos filtering capabilities for 1D, 2D, and multi-dimensional data.
Includes Fortran 90 accelerated implementations for high performance.

References
----------
Duchon, C. E. (1979). Lanczos Filtering in One and Two Dimensions.
    Journal of Applied Meteorology, 18(8), 1016-1022.
"""

from __future__ import annotations

__all__ = [
    "lanczos_filter",
    "lanczos_lowpass",
    "lanczos_highpass",
    "lanczos_bandpass",
    "lanczos_weights",
]

from importlib import import_module
from typing import Any

_LAZY_EXPORTS = {
    "lanczos_filter": ("skyborn.calc.filter.lanczos", "lanczos_filter"),
    "lanczos_lowpass": ("skyborn.calc.filter.lanczos", "lanczos_lowpass"),
    "lanczos_highpass": ("skyborn.calc.filter.lanczos", "lanczos_highpass"),
    "lanczos_bandpass": ("skyborn.calc.filter.lanczos", "lanczos_bandpass"),
    "lanczos_weights": ("skyborn.calc.filter.lanczos", "lanczos_weights"),
}


def __getattr__(name: str) -> Any:
    """Lazy load filter functions on first access."""
    if name in _LAZY_EXPORTS:
        module_name, object_name = _LAZY_EXPORTS[name]
        value = getattr(import_module(module_name), object_name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Include lazy exports in dir() output."""
    return sorted(set(globals()) | set(__all__))
