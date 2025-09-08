"""
Calculation module for Skyborn package.

This module contains various calculation functions including:
- Statistical calculations and linear regression
- Emergent constraint methods
- PDF calculations and analysis
- Mann-Kendall trend analysis
- WMO tropopause calculation (with Fortran extensions)
- Geostrophic wind calculation (with SIMD-optimized Fortran extensions)
- Tropical cyclone potential intensity calculation (with Fortran extensions)
"""

from .calculations import (
    calculate_potential_temperature,
    convert_longitude_range,
    kendall_correlation,
    linear_regression,
    pearson_correlation,
    spatial_correlation,
    spearman_correlation,
)
from .emergent_constraints import (  # New improved function names; Legacy function names for backward compatibility
    calc_GAUSSIAN_PDF,
    calc_PDF_EC,
    calc_PDF_EC_PRIOR,
    emergent_constraint_posterior,
    emergent_constraint_prior,
    find_std_from_PDF,
    gaussian_pdf,
)

# Import geostrophic wind calculation (high-level interface)
from .geostrophic import (
    GeostrophicWind,
    geostrophic_speed,
    geostrophic_uv,
    geostrophic_wind,
)

# Import tropical cyclone potential intensity calculation
from .GPI import potential_intensity
from .mann_kendall import mk_multidim  # alias
from .mann_kendall import mk_test  # alias
from .mann_kendall import (
    mann_kendall_multidim,
    mann_kendall_test,
    mann_kendall_xarray,
    trend_analysis,
)

# Import tropopause calculation (requires compiled extensions)
from .troposphere import trop_wmo, trop_wmo_profile
