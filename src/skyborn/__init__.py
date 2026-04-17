# skyborn/__init__.py

# Import calculation functions
# Import submodules
from . import ROF, calc, gridfill, interp, plot, spharm, windspharm

# Expose calc submodules at package top-level for convenience
from .calc import potential_intensity  # GPI module's main function
from .calc import (  # Mann-Kendall functions; New emergent constraint function names; Legacy names for backward compatibility
    calc_GAUSSIAN_PDF,
    calc_PDF_EC,
    calc_PDF_EC_PRIOR,
    calculate_potential_temperature,
    convert_longitude_range,
    emergent_constraint_posterior,
    emergent_constraint_prior,
    find_std_from_PDF,
    gaussian_pdf,
    geostrophic,
    kendall_correlation,
    linear_regression,
    mann_kendall_multidim,
    mann_kendall_test,
    mann_kendall_xarray,
    mk_multidim,
    mk_test,
    pearson_correlation,
    spearman_correlation,
    trend_analysis,
    troposphere,
)
from .causality import granger_causality, liang_causality
from .gradients import (
    calculate_gradient,
    calculate_meridional_gradient,
    calculate_vertical_gradient,
    calculate_zonal_gradient,
)

# Import key gridfill functions for convenient access
from .gridfill import fill as gridfill_fill

# Expose gridfill functions at top level with clear names
fill = gridfill_fill


__version__ = "0.3.20"
