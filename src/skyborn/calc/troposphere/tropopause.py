"""
WMO tropopause calculation for numpy arrays and masked arrays.

This module provides functions to calculate tropopause properties from gridded
atmospheric data using the WMO definition (lapse rate < 2 K/km). It handles
multi-dimensional arrays with dimension reordering and supports both pressure
and height output.

Key Features:
    - WMO tropopause definition implementation
    - Multi-dimensional array support (3D and 4D)
    - Dimension reordering for arbitrary input layouts
    - Pressure and height calculation
    - Integration with compiled Fortran routines
    - Support for both 1D and multi-dimensional pressure inputs

Mathematical Background:
    The WMO tropopause is defined as the lowest level at which the lapse rate
    decreases to 2 K/km or less, provided also the average lapse rate between
    this level and all higher levels within 2 km doesn't exceed 2 K/km.

Data Requirements:
    - **CRITICAL**: This function requires isobaric (constant pressure) level data.
    - Pressure levels must be sorted in ascending order (low pressure/high altitude
      to high pressure/low altitude).
    - Temperature data must be provided on the same pressure levels.

Examples:
    >>> import numpy as np
    >>> import numpy.ma as ma
    >>> from skyborn.calc.tropopause import trop_wmo
    >>>
    >>> # Create test data (time, level, lat, lon) - isobaric data
    >>> pressure_levels = np.array([10, 20, 50, 100, 200, 300, 500, 700, 850, 1000])  # hPa
    >>> temperature = 300 - np.random.rand(12, 10, 180, 360) * 80  # K
    >>>
    >>> # Calculate tropopause properties
    >>> result = trop_wmo(temperature, pressure_levels,
    ...                   xdim=3, ydim=2, levdim=1, timedim=0)
    >>> print(f"Tropopause pressure shape: {result['pressure'].shape}")
    >>> print(f"Tropopause height shape: {result['height'].shape}")
"""

from __future__ import absolute_import, print_function

import warnings
from typing import Dict, Optional, Tuple, Union

import numpy as np
import numpy.ma as ma

from . import tropopause_height

__all__ = ["trop_wmo", "trop_wmo_profile"]


def _order_dims_for_fortran(
    grid: np.ndarray, xdim: int, ydim: int, levdim: int, timedim: Optional[int] = None
) -> Tuple[np.ndarray, Dict]:
    """
    Reorder array dimensions for Fortran routines.
    Handles 2D, 3D, and 4D data.

    Parameters
    ----------
    grid : np.ndarray
        Input array to reorder
    xdim : int
        Longitude dimension index (-1 if not present)
    ydim : int
        Latitude dimension index (-1 if not present)
    levdim : int
        Vertical level dimension index
    timedim : int, optional
        Time dimension index for 4D data (-1 or None if not present)

    Returns
    -------
    grid : np.ndarray
        Reordered array
    info : dict
        Metadata for dimension restoration
    """
    original_shape = grid.shape
    original_ndim = grid.ndim

    # Build new dimension order based on what dimensions are present
    new_order = []

    # For 2D data, handle (level, lat) or (level, lon) cases
    if original_ndim == 2:
        if ydim >= 0:  # (level, lat) case
            new_order = [ydim, levdim]  # -> (lat, level)
        elif xdim >= 0:  # (level, lon) case
            new_order = [xdim, levdim]  # -> (lon, level)
        else:
            raise ValueError("2D data must have either lat or lon dimension")

    # For 3D data and above, use standard order: (lat, lon, level, [time])
    else:
        if ydim >= 0:
            new_order.append(ydim)
        if xdim >= 0:
            new_order.append(xdim)
        if levdim >= 0:
            new_order.append(levdim)
        if timedim is not None and timedim >= 0:
            new_order.append(timedim)

    # Add any remaining dimensions
    remaining_dims = [i for i in range(original_ndim) if i not in new_order]
    new_order.extend(remaining_dims)

    # Transpose to new order
    grid = np.transpose(grid, new_order)

    # Store restoration info
    info = {
        "original_shape": original_shape,
        "original_ndim": original_ndim,
        "new_order": new_order,
        "inverse_order": np.argsort(new_order),
    }

    return grid, info


def _restore_dims_from_fortran(
    grid: np.ndarray, info: Dict, remove_level_dim: bool = True
) -> np.ndarray:
    """
    Restore original dimension order after Fortran processing.

    Parameters
    ----------
    grid : np.ndarray
        Processed array from Fortran
    info : dict
        Metadata from _order_dims_for_fortran
    remove_level_dim : bool
        If True, accounts for removal of level dimension in output

    Returns
    -------
    np.ndarray
        Array restored to original dimensional structure
    """
    if remove_level_dim:
        # Output arrays don't have level dimension, adjust inverse order
        # Remove level dimension index and shift higher indices down
        inverse_order = []
        levdim_pos = 2  # Level is at position 2 in (lat, lon, level, [time])

        for orig_pos in info["inverse_order"]:
            if orig_pos == levdim_pos:
                continue  # Skip level dimension
            elif orig_pos > levdim_pos:
                inverse_order.append(orig_pos - 1)  # Shift down
            else:
                inverse_order.append(orig_pos)

        inverse_order = np.array(inverse_order)
    else:
        inverse_order = info["inverse_order"]

    # Transpose back to original order
    restored = np.transpose(grid, inverse_order)
    return restored


def trop_wmo(
    temperature: Union[np.ndarray, ma.MaskedArray],
    pressure: Union[np.ndarray, ma.MaskedArray, list],
    xdim: int,
    ydim: int,
    levdim: int,
    timedim: Optional[int] = None,
    pressure_unit: str = "hPa",
    lapse_criterion: float = 2.0,
    missing_value: float = -999.0,
    check_pressure_order: bool = True,
) -> Dict[str, np.ndarray]:
    """
    Calculate WMO tropopause properties for multi-dimensional atmospheric data.

    This function processes gridded atmospheric data to find the tropopause
    following the WMO definition (lapse rate < 2 K/km). Designed for
    **isobaric (constant pressure level) data**.

    Parameters
    ----------
    temperature : array_like
        Atmospheric temperature data [K] on isobaric levels. Must be ordered
        to correspond with the pressure levels.
    pressure : array_like
        Atmospheric pressure data [hPa or Pa]. Can be either:

        - **1D array**: Pressure levels (recommended for isobaric data).
          Length must match the level dimension of temperature.
        - **Multi-dimensional array**: Same shape as temperature (legacy support).

        **CRITICAL**: Pressure levels must be sorted in ascending order
        (from low pressure/high altitude to high pressure/low altitude).
        This is required by the underlying WMO tropopause algorithm.
    xdim : int
        Longitude dimension index
    ydim : int
        Latitude dimension index
    levdim : int
        Vertical level dimension index
    timedim : int, optional
        Time dimension index for 4D data
    pressure_unit : str, default 'hPa'
        Unit of pressure data ('hPa' or 'Pa')
    lapse_criterion : float, default 2.0
        WMO lapse rate criterion [K/km] for tropopause definition
    missing_value : float, default -999.0
        Value to use for missing data
    check_pressure_order : bool, default True
        Whether to check and enforce ascending pressure order along the level dimension

    Returns
    -------
    dict
        Dictionary containing:
        - 'pressure': Tropopause pressure [hPa]
        - 'height': Tropopause height [m]
        - 'level_index': Tropopause level index (0-based)
        - 'lapse_rate': Tropopause lapse rate [K/km]
        - 'success': Success flag for each grid point

    Examples
    --------
    Calculate tropopause for 4D isobaric data (time, level, lat, lon):

    >>> # Define isobaric pressure levels (ascending order: high altitude to surface)
    >>> pressure_levels = np.array([10, 20, 50, 100, 200, 300, 500, 700, 850, 1000])  # hPa
    >>> temperature = 300 - np.random.rand(12, 10, 180, 360) * 80  # K
    >>> result = trop_wmo(temperature, pressure_levels,
    ...                   xdim=3, ydim=2, levdim=1, timedim=0)
    >>> print(f"Pressure shape: {result['pressure'].shape}")  # (12, 180, 360)
    >>> print(f"Height shape: {result['height'].shape}")      # (12, 180, 360)

    Calculate tropopause for 3D isobaric data (level, lat, lon):

    >>> pressure_levels = np.array([50, 100, 200, 300, 500, 700, 850, 1000])  # hPa
    >>> temperature_3d = 300 - np.random.rand(8, 180, 360) * 80  # K
    >>> result = trop_wmo(temperature_3d, pressure_levels,
    ...                   xdim=2, ydim=1, levdim=0)
    >>> print(f"Pressure shape: {result['pressure'].shape}")  # (180, 360)

    Notes
    -----
    - This function is optimized for **isobaric data** (constant pressure levels).
    - For model level data, first interpolate to pressure levels.
    - Requires compiled Fortran extensions. Install with: pip install skyborn[fortran]
    - The underlying algorithm follows the WMO (1957) tropopause definition.
    """

    # Convert inputs to arrays
    pressure = np.asarray(pressure, dtype=np.float32)
    temperature = np.asarray(temperature, dtype=np.float32)

    # Check if pressure length matches temperature level dimension
    temp_level_size = temperature.shape[levdim]
    if len(pressure) != temp_level_size:
        raise ValueError(
            f"Pressure length ({len(pressure)}) must match temperature "
            f"level dimension size ({temp_level_size}). "
            f"Pressure shape: {pressure.shape}, Temperature shape: {temperature.shape}, levdim: {levdim}"
        )

    # Use temperature dimensions for validation
    ndim = temperature.ndim
    dims_to_check = []

    # Only check dimensions that are not -1 (placeholder for missing)
    if xdim >= 0:
        dims_to_check.append(xdim)
    if ydim >= 0:
        dims_to_check.append(ydim)
    if levdim >= 0:
        dims_to_check.append(levdim)
    if timedim is not None and timedim >= 0:
        dims_to_check.append(timedim)

    if dims_to_check and (max(dims_to_check) >= ndim or min(dims_to_check) < 0):
        raise ValueError(
            f"Dimension indices must be valid for {ndim}D temperature array"
        )

    if len(set(dims_to_check)) != len(dims_to_check):
        raise ValueError("Dimension indices must be unique")

    # Handle masked arrays
    if ma.is_masked(pressure):
        pressure = ma.filled(pressure, missing_value)
    if ma.is_masked(temperature):
        temperature = ma.filled(temperature, missing_value)

    # Reorder dimensions to (lat, lon, level, [time])
    # 1D pressure doesn't need reordering, use temperature dimensions for ordering
    temperature_ordered, dim_info = _order_dims_for_fortran(
        temperature, xdim, ydim, levdim, timedim
    )
    pressure_ordered = pressure  # Keep 1D pressure as is

    # Check and sort pressure levels if needed
    if check_pressure_order:
        # Check 1D pressure array directly
        if not np.all(pressure_ordered[:-1] <= pressure_ordered[1:]):
            # Need to sort - get sorting indices
            sort_indices = np.argsort(pressure_ordered)

            # Apply sorting to both pressure and temperature
            pressure_ordered = pressure_ordered[sort_indices]
            # Apply same sorting to temperature along level dimension
            # Find the level dimension position in the reordered temperature array
            # For 2D: (level, lat) -> (lat, level), so level is now axis 1
            # For 3D: (level, lat, lon) -> (lat, lon, level), so level is now axis 2
            # For 4D: (time, level, lat, lon) -> (lat, lon, level, time), so level is now axis 2
            if temperature_ordered.ndim == 2:
                level_axis = 1  # (spatial, level)
            elif temperature_ordered.ndim == 3:
                level_axis = 2  # (lat, lon, level)
            else:  # 4D
                level_axis = 2  # (lat, lon, level, time)

            temperature_ordered = np.take(
                temperature_ordered, sort_indices, axis=level_axis
            )

    # Get dimensions based on temperature data shape (since pressure is now 1D)
    shape = temperature_ordered.shape

    # Get number of levels from pressure (which is always 1D now)
    nlev = len(pressure_ordered)

    if len(shape) == 2:  # 2D data (spatial_dim, level)
        nspatial = shape[0]
        nlat = nlon = nspatial  # Use same value for missing dimension
        ntime = None
        is_4d = False
        is_2d = True
    elif len(shape) == 3:  # 3D data (lat, lon, level)
        nlat, nlon, _ = shape  # level dimension from pressure
        ntime = None
        is_4d = False
        is_2d = False
    elif len(shape) == 4:  # 4D data (lat, lon, level, time)
        nlat, nlon, _, ntime = shape  # level dimension from pressure
        is_4d = True
        is_2d = False
    else:
        raise ValueError(f"Unsupported temperature data shape: {shape}")

    nlevm = nlev + 1

    # Set pressure unit flag (0 for hPa, 1 for Pa)
    punit = 0 if pressure_unit.lower() == "hpa" else 1

    # Use the 1D pressure array (already sorted if needed)
    pressure_levels_1d = pressure_ordered

    # Call appropriate Fortran routine based on dimensionality
    if is_2d:
        # Use 2D routine with 1D pressure
        ptrop_hpa, htrop_m, itrop, lapse_rate, success = (
            tropopause_height.tropopause_grid_2d(
                nlevm,
                pressure_levels_1d,
                temperature_ordered,
                missing_value,
                lapse_criterion,
                punit,
            )
        )

    elif is_4d:
        # Use 4D routine with 1D pressure
        ptrop_hpa, htrop_m, itrop, lapse_rate, success = (
            tropopause_height.tropopause_grid_4d(
                nlevm,
                pressure_levels_1d,
                temperature_ordered,
                missing_value,
                lapse_criterion,
                punit,
            )
        )
    else:  # 3D case
        # Use 3D routine with 1D pressure
        ptrop_hpa, htrop_m, itrop, lapse_rate, success = (
            tropopause_height.tropopause_grid_3d(
                nlevm,
                pressure_levels_1d,
                temperature_ordered,
                missing_value,
                lapse_criterion,
                punit,
            )
        )

    # Restore original dimension order (removing level dimension)
    if not is_2d:  # For 2D data, results are already in correct shape
        ptrop_hpa = _restore_dims_from_fortran(
            ptrop_hpa, dim_info, remove_level_dim=True
        )
        htrop_m = _restore_dims_from_fortran(htrop_m, dim_info, remove_level_dim=True)
        itrop = _restore_dims_from_fortran(itrop, dim_info, remove_level_dim=True)
        lapse_rate = _restore_dims_from_fortran(
            lapse_rate, dim_info, remove_level_dim=True
        )
        success = _restore_dims_from_fortran(success, dim_info, remove_level_dim=True)

    # Create output dictionary
    result = {
        "pressure": ptrop_hpa,
        "height": htrop_m,
        "level_index": itrop,
        "lapse_rate": lapse_rate,
        "success": success,
    }

    return result


def trop_wmo_profile(
    temperature: Union[np.ndarray, ma.MaskedArray],
    pressure: Union[np.ndarray, ma.MaskedArray, list],
    pressure_unit: str = "hPa",
    lapse_criterion: float = 2.0,
    missing_value: float = -999.0,
) -> Dict[str, Union[float, int, bool]]:
    """
    Calculate WMO tropopause properties for a single vertical profile.

    This function processes a single atmospheric vertical profile to find
    the tropopause following the WMO definition (lapse rate < 2 K/km).
    Optimized for profile data analysis with **isobaric level data**.

    Parameters
    ----------
    temperature : array_like
        Atmospheric temperature profile [K] on isobaric levels. 1D array.
    pressure : array_like
        Atmospheric pressure profile [hPa or Pa]. Must have same length as temperature.
        **CRITICAL**: Pressure levels must be sorted in ascending order
        (from low pressure/high altitude to high pressure/low altitude).
        This is required by the underlying WMO tropopause algorithm.
    pressure_unit : str, default 'hPa'
        Unit of pressure data ('hPa' or 'Pa')
    lapse_criterion : float, default 2.0
        WMO lapse rate criterion [K/km] for tropopause definition
    missing_value : float, default -999.0
        Value to use for missing data

    Returns
    -------
    dict
        Dictionary containing:
        - 'pressure': Tropopause pressure [hPa] (float)
        - 'height': Tropopause height [m] (float)
        - 'level_index': Tropopause level index (int, 0-based)
        - 'lapse_rate': Tropopause lapse rate [K/km] (float)
        - 'success': Success flag (bool)

    Examples
    --------
    Basic usage with isobaric profile data:

    >>> import numpy as np
    >>> from skyborn.calc.troposphere.tropopause import trop_wmo_profile
    >>>
    >>> # Create atmospheric profile on isobaric levels (ascending order)
    >>> pressure = np.array([10, 20, 50, 100, 200, 300, 500, 700, 850, 1000])  # hPa
    >>> temperature = np.array([190, 200, 210, 220, 230, 250, 270, 280, 285, 288])  # K
    >>>
    >>> result = trop_wmo_profile(temperature, pressure)
    >>> print(f"Tropopause pressure: {result['pressure']:.1f} hPa")
    >>> print(f"Tropopause height: {result['height']:.0f} m")

    Notes
    -----
    - This function is optimized for **isobaric data** (constant pressure levels).
    - For model level data, first interpolate to pressure levels.
    - Uses dedicated Fortran routine for single profile analysis.
    - For gridded data, use trop_wmo() instead.
    """
    # Validate inputs
    pressure = np.asarray(pressure, dtype=np.float32)
    temperature = np.asarray(temperature, dtype=np.float32)

    if pressure.ndim != 1 or temperature.ndim != 1:
        raise ValueError("Profile inputs must be 1D arrays")

    if len(pressure) != len(temperature):
        raise ValueError("Pressure and temperature profiles must have same length")

    nlev = len(pressure)
    nlevm = nlev + 1

    # Handle masked arrays
    if ma.is_masked(pressure):
        pressure = ma.filled(pressure, missing_value)
    if ma.is_masked(temperature):
        temperature = ma.filled(temperature, missing_value)

    # Set pressure unit flag (0 for hPa, 1 for Pa)
    punit = 0 if pressure_unit.lower() == "hpa" else 1

    # Use dedicated 1D profile function
    ptrop_hpa, htrop_m, itrop, lapse_rate, success = (
        tropopause_height.tropopause_profile_1d(
            nlevm, pressure, temperature, missing_value, lapse_criterion, punit
        )
    )

    # Create output dictionary
    result = {
        "pressure": ptrop_hpa,
        "height": htrop_m,
        "level_index": itrop,
        "lapse_rate": lapse_rate,
        "success": success,
    }

    return result
