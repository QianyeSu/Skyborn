# skyborn/__init__.py
from .calculations import convert_longitude_range, linear_regression

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

from . import plot
from . import interp
from . import ROF
from . import conversion

__version__ = "0.3.6"  # Updated to version 0.3.6
