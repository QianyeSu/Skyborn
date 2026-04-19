"""Public Mann-Kendall package exports."""

from .core import (
    mann_kendall_multidim,
    mann_kendall_test,
    mann_kendall_xarray,
    mk_multidim,
    mk_test,
    trend_analysis,
)
from .partial import (
    partial_mann_kendall_multidim,
    partial_mann_kendall_test,
    partial_mann_kendall_xarray,
    partial_multidim,
    partial_test,
)

__all__ = [
    "mann_kendall_test",
    "mann_kendall_multidim",
    "mann_kendall_xarray",
    "partial_mann_kendall_test",
    "partial_mann_kendall_multidim",
    "partial_mann_kendall_xarray",
    "partial_test",
    "partial_multidim",
    "trend_analysis",
    "mk_test",
    "mk_multidim",
]
