"""
WMO tropopause calculation for xarray DataArrays.

This module provides functions to calculate tropopause properties from xarray
DataArrays using the WMO definition. It automatically detects spatial coordinates
and preserves coordinate information and metadata throughout the computation process.

*** DATA REQUIREMENTS ***
- CRITICAL: This function requires ISOBARIC (constant pressure level) data
- Temperature data must be provided on constant pressure levels
- Pressure levels must be sorted in ASCENDING order (low pressure/high altitude to high pressure/low altitude)
- For model level data, first interpolate to pressure levels before using this function

Main Functions:
    trop_wmo : Calculate WMO tropopause properties for xarray DataArray

Examples:
    >>> import xarray as xr
    >>> import numpy as np
    >>> from skyborn.calc.troposphere.xarray import trop_wmo
    >>>
    >>> # Load isobaric atmospheric data
    >>> ds = xr.open_dataset('era5_pressure_levels.nc')  # Already on isobaric levels
    >>> result = trop_wmo(ds.temperature)  # Auto-generates pressure from level coordinate
    >>> print(result.pressure.attrs)  # Original attributes preserved
"""

from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

__all__ = ["trop_wmo"]


import numpy as np
import numpy.ma as ma
import xarray as xr

from . import tropopause

# Type aliases
DataArray = xr.DataArray


def _detect_atmospheric_dimensions(
    data_array: DataArray,
) -> Tuple[Optional[int], Optional[int], int, Optional[int]]:
    """
    Auto-detect dimension indices for atmospheric data in xarray DataArray.
    Supports 2D, 3D, and 4D data.

    Parameters
    ----------
    data_array : xr.DataArray
        Atmospheric data to analyze

    Returns
    -------
    xdim : int, optional
        Longitude dimension index (None for 2D data without longitude)
    ydim : int, optional
        Latitude dimension index (None for 2D data without latitude)
    levdim : int
        Vertical level dimension index
    timedim : int, optional
        Time dimension index (None if not found)

    Raises
    ------
    ValueError
        If level dimension cannot be identified
    """
    dims = data_array.dims

    # Common dimension name patterns
    lon_names = {
        "lon",
        "longitude",
        "x",
        "X",
        "LON",
        "XLON",
        "LONS",
        "LONG",
        "LONGITUDE",
    }
    lat_names = {"lat", "latitude", "y", "Y", "LAT", "YLAT", "LATS", "LATI", "LATITUDE"}
    lev_names = {
        "level",
        "lev",
        "plev",
        "pressure",
        "pressure_level",
        "z",
        "Z",
        "LEV",
        "PRES",
        "LEVEL",
        "PLEVEL",
        "PRES_LEVEL",
        "HEIGHT",
        "height",
        "altitude",
        "ALTITUDE",
        "depth",
        "DEPTH",
        "isobaric",
        "ISOBARIC",
        "model_level",
        "PRESSURE_LEVEL",
        "PRESSURE",
    }
    time_names = {"time", "t", "T", "year", "month", "yr", "mn", "season"}

    xdim = ydim = levdim = timedim = None

    for i, dim_name in enumerate(dims):
        dim_lower = dim_name.lower()

        if any(name in dim_lower for name in lon_names):
            xdim = i
        elif any(name in dim_lower for name in lat_names):
            ydim = i
        elif any(name in dim_lower for name in lev_names):
            levdim = i
        elif any(name in dim_lower for name in time_names):
            timedim = i

    # Level dimension is required
    if levdim is None:
        raise ValueError(
            f"Could not auto-detect level dimension. Found dims: {dims}. "
            f"Please specify levdim explicitly."
        )

    # Check dimension requirements based on data dimensionality
    ndim = data_array.ndim
    if ndim == 1:
        # 1D profile data - only need level dimension
        pass
    elif ndim == 2:
        # 2D data - need at least one spatial dimension (lat or lon)
        if xdim is None and ydim is None:
            raise ValueError(
                f"For 2D data, need at least one spatial dimension (lat or lon). "
                f"Found dims: {dims}"
            )
    elif ndim >= 3:
        # 3D+ data - need both lat and lon dimensions
        if xdim is None or ydim is None:
            raise ValueError(
                f"For {ndim}D data, need both lat and lon dimensions. Found dims: {dims}. "
                f"Please specify xdim and ydim explicitly."
            )

    return xdim, ydim, levdim, timedim


def trop_wmo(
    temperature: DataArray,
    pressure: Optional[DataArray] = None,
    xdim: Optional[Union[int, str]] = None,
    ydim: Optional[Union[int, str]] = None,
    levdim: Optional[Union[int, str]] = None,
    timedim: Optional[Union[int, str]] = None,
    pressure_unit: str = "hPa",
    lapse_criterion: float = 2.0,
    missing_value: float = -999.0,
    keep_attrs: bool = True,
    auto_sort_levels: bool = True,
) -> xr.Dataset:
    """
    Calculate WMO tropopause properties for xarray DataArrays.

    This function processes gridded atmospheric data to find the tropopause
    following the WMO definition (lapse rate < 2 K/km). It automatically
    detects coordinate dimensions and preserves all metadata.

    *** DESIGNED FOR ISOBARIC DATA ***

    Parameters
    ----------
    temperature : xarray.DataArray
        Atmospheric temperature data [K] on isobaric (constant pressure) levels.
        Can be 1D profile, 2D, 3D, or 4D array.
        Must have level coordinate (for pressure generation) if pressure is not provided.
    pressure : xarray.DataArray, optional
        Atmospheric pressure data [hPa or Pa] on isobaric levels. If None, will be
        automatically generated from temperature's level coordinate.
        **CRITICAL**: Pressure levels MUST be sorted in ASCENDING order
        (from low pressure/high altitude to high pressure/low altitude).
        This is required by the underlying WMO tropopause algorithm.
    xdim : int or str, optional
        Longitude dimension index/name. Auto-detected if None. Not required for 1D profiles.
    ydim : int or str, optional
        Latitude dimension index/name. Auto-detected if None. Not required for 1D profiles.
    levdim : int or str, optional
        Vertical level dimension index/name. Auto-detected if None.
    timedim : int or str, optional
        Time dimension index/name. Auto-detected if None.
    pressure_unit : str, default 'hPa'
        Unit of pressure data ('hPa' or 'Pa')
    lapse_criterion : float, default 2.0
        WMO lapse rate criterion [K/km] for tropopause definition
    missing_value : float, default -999.0
        Value to use for missing data
    keep_attrs : bool, default True
        Preserve input DataArray attributes in output
    auto_sort_levels : bool, default True
        Automatically sort pressure levels in ascending order along the level dimension

    Returns
    -------
    xarray.Dataset
        Dataset containing tropopause properties:

        For multi-dimensional data (2D, 3D, 4D):
        - 'pressure': Tropopause pressure [hPa] with spatial/temporal coordinates
        - 'height': Tropopause height [m] with spatial/temporal coordinates
        - 'level_index': Tropopause level index (0-based)
        - 'lapse_rate': Tropopause lapse rate [K/km]
        - 'success': Success flag for each grid point

        For 1D profile data:
        - Same variables but as scalar values (0D DataArrays)

    Examples
    --------
    **1D Profile Analysis (Isobaric Data):**

    >>> import xarray as xr
    >>> import numpy as np
    >>> from skyborn.calc.troposphere.xarray import trop_wmo
    >>>
    >>> # Create 1D isobaric profile (ascending pressure order)
    >>> temp_profile = xr.DataArray(
    ...     [210, 230, 250, 270, 280, 288],  # Temperature decreasing with altitude
    ...     dims=['level'],
    ...     coords={'level': [100, 300, 500, 700, 850, 1000]}  # hPa - ASCENDING order
    ... )
    >>> result = trop_wmo(temp_profile)
    >>> print(f"Tropopause: {float(result.pressure)} hPa, {float(result.height)} m")

    **Simplified Interface (Auto-pressure generation from isobaric levels):**

    >>> # Load ERA5 pressure level data (already isobaric)
    >>> ds = xr.open_dataset('era5_pressure_levels.nc')  # Has 'level' coordinate in hPa
    >>> result = trop_wmo(ds.temperature)  # Pressure auto-generated from level coordinate
    >>> print(f"Tropopause pressure shape: {result.pressure.shape}")

    **2D Spatial Analysis (Isobaric Cross-sections):**

    >>> # Analyze latitude or longitude cross-sections on isobaric levels
    >>> temp_2d = ds.temperature.isel(time=0, lon=0)  # (level, lat) - isobaric levels
    >>> result = trop_wmo(temp_2d)
    >>> # Result shape: (lat,)

    **Advanced usage with explicit isobaric pressure:**

    >>> # Ensure pressure levels are in ascending order
    >>> result = trop_wmo(
    ...     temperature_data,  # On isobaric levels
    ...     pressure=pressure_data,  # Corresponding isobaric pressure levels
    ...     xdim='longitude', ydim='latitude', levdim='level',
    ...     lapse_criterion=2.5  # Custom WMO criterion
    ... )

    **4D Time Series (Isobaric Data):**

    >>> # Multi-year isobaric climate data
    >>> result = trop_wmo(temperature_4d)  # (time, level, lat, lon) - isobaric levels
    >>> # Result preserves time dimension: (time, lat, lon)
    >>> seasonal_mean = result.height.groupby('time.season').mean()

    Notes
    -----
    - This function is optimized for **ISOBARIC data** (constant pressure levels).
    - For model level data, first interpolate to pressure levels before using this function.
    - Requires compiled Fortran extensions. Install with: pip install skyborn[fortran]
    - The underlying algorithm follows the WMO (1957) tropopause definition.

    The function automatically:
    - Detects spatial and temporal coordinates using metadata
    - Handles missing values (NaN or masked arrays)
    - Preserves all coordinate information and attributes
    - Works with multi-dimensional isobaric data

    See Also
    --------
    skyborn.calc.troposphere.tropopause.trop_wmo : Lower-level function for numpy arrays
    """
    # Validate input types
    if not isinstance(temperature, xr.DataArray):
        raise TypeError(
            f"temperature must be xarray.DataArray, got {type(temperature).__name__}"
        )

    if pressure is not None and not isinstance(pressure, xr.DataArray):
        raise TypeError(
            f"pressure must be xarray.DataArray or None, got {type(pressure).__name__}"
        )

    # Generate pressure from temperature level coordinate if not provided
    if pressure is None:
        # Auto-detect dimensions from temperature first
        xdim_auto, ydim_auto, levdim_auto, timedim_auto = (
            _detect_atmospheric_dimensions(temperature)
        )

        # Use detected level dimension
        level_dim_name = list(temperature.dims)[levdim_auto]

        if level_dim_name not in temperature.coords:
            raise ValueError(
                f"Cannot generate pressure: temperature must have '{level_dim_name}' coordinate "
                f"when pressure is not provided. Available coordinates: {list(temperature.coords.keys())}"
            )

        # Get level coordinate values (assumed to be pressure levels)
        level_coord = temperature.coords[level_dim_name]
        pressure_levels = level_coord.values

        # Create 1D pressure array (optimized for isobaric data)
        pressure = xr.DataArray(
            pressure_levels,
            dims=[level_dim_name],
            coords={level_dim_name: level_coord},
            attrs={
                "units": pressure_unit,
                "long_name": "Atmospheric pressure levels",
                "description": f"1D isobaric pressure levels from {level_dim_name} coordinate",
            },
        )

        print(
            f"Generated pressure from level coordinate '{level_dim_name}' "
            f"with {len(pressure_levels)} levels ({pressure_levels.min():.1f}-{pressure_levels.max():.1f} {pressure_unit})"
        )

    # Auto-detect dimensions if not provided
    # Use temperature for dimension detection since pressure might be 1D
    if (
        levdim is None
        or (temperature.ndim >= 3 and (xdim is None or ydim is None))
        or (temperature.ndim == 2 and xdim is None and ydim is None)
    ):
        xdim_auto, ydim_auto, levdim_auto, timedim_auto = (
            _detect_atmospheric_dimensions(temperature)
        )
        if xdim is None:
            xdim = xdim_auto
        if ydim is None:
            ydim = ydim_auto
        if levdim is None:
            levdim = levdim_auto
        if timedim is None:
            timedim = timedim_auto

    # Store coordinate information for output (use temperature since pressure might be 1D)
    original_coords = temperature.coords
    original_dims = temperature.dims

    # Convert dimension names to indices if needed (based on temperature dimensions)
    if isinstance(xdim, str):
        xdim = list(temperature.dims).index(xdim) if xdim in temperature.dims else None
    if isinstance(ydim, str):
        ydim = list(temperature.dims).index(ydim) if ydim in temperature.dims else None
    if isinstance(levdim, str):
        levdim = (
            list(temperature.dims).index(levdim) if levdim in temperature.dims else None
        )
    if isinstance(timedim, str):
        timedim = (
            list(temperature.dims).index(timedim)
            if timedim in temperature.dims
            else None
        )

    # Sort pressure levels if requested
    if auto_sort_levels:
        # Get level dimension name for sorting (from temperature)
        level_dim_name = list(temperature.dims)[levdim]

        # For 1D pressure, check if it's sorted
        if pressure.ndim == 1:
            # Check 1D pressure values directly
            if not np.all(pressure.values[:-1] <= pressure.values[1:]):
                # Sort both pressure and temperature
                sort_indices = np.argsort(pressure.values)
                pressure = pressure.isel({pressure.dims[0]: sort_indices})
                temperature = temperature.isel({level_dim_name: sort_indices})
        else:
            # Multi-dimensional pressure - use existing logic
            # Check if sorting is needed by examining pressure coordinate
            if level_dim_name in pressure.coords:
                level_coord = pressure.coords[level_dim_name]
                level_values = level_coord.values

                # Check if already sorted in ascending order
                if not np.all(level_values[:-1] <= level_values[1:]):
                    # Sort both pressure and temperature by level coordinate
                    pressure = pressure.sortby(level_dim_name)
                    temperature = temperature.sortby(level_dim_name)
            else:
                # No level coordinate, check actual pressure values
                # Take a sample pressure profile to check ordering
                sample_indices = {
                    dim: 0 for dim in pressure.dims if dim != level_dim_name
                }
                sample_profile = pressure.isel(**sample_indices)

                if not np.all(sample_profile.values[:-1] <= sample_profile.values[1:]):
                    # Need to sort - create sorting index
                    sort_indices = np.argsort(sample_profile.values)

                    # Apply sorting to both arrays
                    pressure = pressure.isel({level_dim_name: sort_indices})
                    temperature = temperature.isel({level_dim_name: sort_indices})

    # Extract numpy arrays
    pressure_data = pressure.values
    temperature_data = temperature.values

    # Handle different dimensionalities
    if temperature.ndim == 1:
        # 1D profile - use dedicated profile function
        result = tropopause.trop_wmo_profile(
            temperature_data,
            pressure_data,
            pressure_unit=pressure_unit,
            lapse_criterion=lapse_criterion,
            missing_value=missing_value,
        )

        # Convert scalar results to 0D arrays for consistency
        for key in result:
            if not isinstance(result[key], np.ndarray):
                result[key] = np.array(result[key])
    else:
        # Multi-dimensional data - use grid function with optimized 1D pressure
        # Use -1 as a placeholder for missing dimensions
        xdim_arg = xdim if xdim is not None else -1
        ydim_arg = ydim if ydim is not None else -1
        timedim_arg = timedim if timedim is not None else -1

        # Call the core tropopause calculation function with 1D pressure optimization
        result = tropopause.trop_wmo(
            temperature_data,
            pressure_data,  # This will be 1D for isobaric data
            xdim=xdim_arg,
            ydim=ydim_arg,
            levdim=levdim,
            timedim=timedim_arg,
            pressure_unit=pressure_unit,
            lapse_criterion=lapse_criterion,
            missing_value=missing_value,
            check_pressure_order=not auto_sort_levels,  # Skip check if we already sorted
        )

    # Create output Dataset with proper coordinates
    if temperature.ndim == 1:
        # For 1D profiles, output is scalar
        output_dims = []
        output_coords = {}
    else:
        # For multi-dimensional data, remove level dimension
        output_dims = list(temperature.dims)
        level_dim_name = output_dims.pop(levdim)

        # Create coordinates for output (excluding level)
        output_coords = {}
        for dim_name in output_dims:
            if dim_name in temperature.coords:
                output_coords[dim_name] = temperature.coords[dim_name]

    # Create DataArrays for each output variable
    data_vars = {}

    # Pressure variable
    data_vars["pressure"] = xr.DataArray(
        result["pressure"],
        dims=output_dims,
        coords=output_coords,
        attrs={
            "long_name": "Tropopause pressure",
            "units": "hPa",
            "standard_name": "tropopause_air_pressure",
            "description": "WMO tropopause pressure calculated using lapse rate criterion",
        },
    )

    # Height variable
    data_vars["height"] = xr.DataArray(
        result["height"],
        dims=output_dims,
        coords=output_coords,
        attrs={
            "long_name": "Tropopause height",
            "units": "m",
            "standard_name": "tropopause_altitude",
            "description": "WMO tropopause height above sea level",
        },
    )

    # Level index variable
    data_vars["level_index"] = xr.DataArray(
        result["level_index"],
        dims=output_dims,
        coords=output_coords,
        attrs={
            "long_name": "Tropopause level index",
            "units": "1",
            "description": "Zero-based index of tropopause level in input array",
        },
    )

    # Lapse rate variable
    data_vars["lapse_rate"] = xr.DataArray(
        result["lapse_rate"],
        dims=output_dims,
        coords=output_coords,
        attrs={
            "long_name": "Tropopause lapse rate",
            "units": "K km-1",
            "standard_name": "air_temperature_lapse_rate",
            "description": f"Temperature lapse rate at tropopause (threshold: {lapse_criterion} K/km)",
        },
    )

    # Success flag variable
    data_vars["success"] = xr.DataArray(
        result["success"],
        dims=output_dims,
        coords=output_coords,
        attrs={
            "long_name": "Tropopause calculation success flag",
            "description": "True where tropopause was successfully identified",
        },
    )

    # Create Dataset
    ds = xr.Dataset(data_vars)

    # Add global attributes
    ds.attrs = {
        "title": "WMO tropopause calculation results",
        "description": "Tropopause properties calculated using WMO definition",
        "lapse_rate_criterion_K_per_km": lapse_criterion,
        "pressure_unit": pressure_unit,
        "missing_value": missing_value,
        "method": "WMO tropopause definition (lapse rate criterion)",
        "software": "skyborn atmospheric calculation package",
    }

    # Preserve original attributes if requested
    if keep_attrs:
        if hasattr(pressure, "attrs") and pressure.attrs:
            ds.attrs.update(
                {f"source_pressure_{k}": v for k, v in pressure.attrs.items()}
            )
        if hasattr(temperature, "attrs") and temperature.attrs:
            ds.attrs.update(
                {f"source_temperature_{k}": v for k, v in temperature.attrs.items()}
            )

    return ds
