"""
Calculation module for Skyborn package.

This module contains various calculation functions including:
- Statistical calculations and linear regression
- Emergent constraint methods
- PDF calculations and analysis
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
]
