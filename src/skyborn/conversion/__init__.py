"""
Data format conversion utilities for the skyborn library.

This module provides functions for converting between different atmospheric
data formats, including GRIB to NetCDF conversion using eccodes.

Author: Qianye Su
Email: suqianye2000@gmail.com
"""

import subprocess  # For testing patches

from .grib_to_netcdf import (
    convert_grib_to_nc,
    convert_grib_to_nc_simple,
    batch_convert_grib_to_nc,
    grib2nc,
    grib_to_netcdf,
    GribToNetCDFError,
    # Private functions for testing
    _check_grib_to_netcdf_available,
    _validate_grib_files,
    _build_grib_to_netcdf_command,
)

__all__ = [
    "convert_grib_to_nc",
    "convert_grib_to_nc_simple",
    "batch_convert_grib_to_nc",
    "grib2nc",
    "grib_to_netcdf",
    "GribToNetCDFError",
    # Private functions (for testing)
    "_check_grib_to_netcdf_available",
    "_validate_grib_files",
    "_build_grib_to_netcdf_command",
]
