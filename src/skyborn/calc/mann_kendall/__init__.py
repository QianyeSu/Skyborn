"""Public Mann-Kendall package exports."""

from .core import (
    mann_kendall_multidim,
    mann_kendall_test,
    mann_kendall_xarray,
    mk_multidim,
    mk_test,
    trend_analysis,
)

__all__ = [
    "mann_kendall_test",
    "mann_kendall_multidim",
    "mann_kendall_xarray",
    "trend_analysis",
    "mk_test",
    "mk_multidim",
]
