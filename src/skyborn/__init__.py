# skyborn/__init__.py

# Import calculation functions
from .calc import (
    linear_regression,
    convert_longitude_range,
    pearson_correlation,
    spearman_correlation,
    kendall_correlation,
    calculate_potential_temperature,
    # New emergent constraint function names
    gaussian_pdf,
    emergent_constraint_posterior,
    emergent_constraint_prior,
    # Legacy names for backward compatibility
    calc_GAUSSIAN_PDF,
    calc_PDF_EC,
    find_std_from_PDF,
    calc_PDF_EC_PRIOR,
)

from .gradients import (
    calculate_gradient,
    calculate_meridional_gradient,
    calculate_zonal_gradient,
    calculate_vertical_gradient,
)

from .causality import liang_causality, granger_causality

# Import conversion functions for easy access
from .conversion import (
    convert_grib_to_nc,
    convert_grib_to_nc_simple,
    batch_convert_grib_to_nc,
    grib2nc,
    grib_to_netcdf,
)

# Import submodules
from . import plot
from . import interp
from . import ROF
from . import conversion
from . import calc
from . import spharm
from . import windspharm
from . import gridfill

# Import key gridfill functions for convenient access
# (Users can use either skyborn.fill or skyborn.gridfill.fill)
from .gridfill import fill as gridfill_fill, fill_cube as gridfill_fill_cube

# Expose gridfill functions at top level with clear names
fill = gridfill_fill
fill_cube = gridfill_fill_cube

# Note: Enhanced gridfill module now includes modern xarray interface:
# from skyborn.gridfill.xarray import fill  # For xarray DataArrays

__version__ = "0.3.10"  # Updated to version 0.3.10
