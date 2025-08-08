"""
Calculation module for Skyborn package.

This module contains various calculation functions including:
- Statistical calculations and linear regression
- Emergent constraint methods
- PDF calculations and analysis
- Mann-Kendall trend analysis
"""

from .calculations import (
    linear_regression,
    convert_longitude_range,
    pearson_correlation,
    spearman_correlation,
    kendall_correlation,
    calculate_potential_temperature,
)

from .emergent_constraints import (
    # New improved function names
    gaussian_pdf,
    emergent_constraint_posterior,
    emergent_constraint_prior,
    # Legacy function names for backward compatibility
    calc_GAUSSIAN_PDF,
    calc_PDF_EC,
    find_std_from_PDF,
    calc_PDF_EC_PRIOR,
)

from .mann_kendall import (
    mann_kendall_test,
    mann_kendall_multidim_numpy,
    trend_analysis,
)

# Try to import xarray-specific functions if xarray is available
try:
    from .mann_kendall import mann_kendall_xarray

    _has_xarray_mk = True
except ImportError:
    _has_xarray_mk = False

__all__ = [
    # From calculations.py
    "linear_regression",
    "convert_longitude_range",
    "pearson_correlation",
    "spearman_correlation",
    "kendall_correlation",
    "calculate_potential_temperature",
    # From emergent_constraints.py - New names
    "gaussian_pdf",
    "emergent_constraint_posterior",
    "emergent_constraint_prior",
    # Legacy names for backward compatibility
    "calc_GAUSSIAN_PDF",
    "calc_PDF_EC",
    "find_std_from_PDF",
    "calc_PDF_EC_PRIOR",
    # From mann_kendall.py
    "mann_kendall_test",
    "mann_kendall_multidim_numpy",
    "trend_analysis",
]

# Add xarray function to __all__ if available
if _has_xarray_mk:
    __all__.append("mann_kendall_xarray")
