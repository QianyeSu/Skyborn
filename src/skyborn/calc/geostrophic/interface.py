"""
High-level Python interface for geostrophic wind calculations.

This module provides user-friendly interfaces for calculating geostrophic winds
from geopotential height fields, with support for multi-dimensional data using
the windspharm data preparation utilities.

The interface handles automatic data reshaping, dimension reordering, and
integration with the optimized Fortran backend.
"""

import warnings
from typing import Optional, Tuple, Union

import numpy as np

# Import windspharm tools for data preparation
from skyborn.windspharm.tools import prep_data, recover_data

# Import the compiled Fortran functions (only 2D and 3D needed)
from . import geostrophicwind as _geostrophic_module

z2geouv = _geostrophic_module.z2geouv  # For 2D data
# For 3D data (handles combined dimensions)
z2geouv_3d = _geostrophic_module.z2geouv_3d


def _is_longitude_cyclic(glon: np.ndarray, tolerance: float = 1.0) -> bool:
    """
    Determine if longitude data is cyclic by checking if it spans 360 degrees.

    Works with different grid resolutions by using adaptive tolerance based
    on grid spacing.

    Parameters
    ----------
    glon : ndarray
        Longitude coordinates in degrees
    tolerance : float
        Base tolerance for cyclicity check (default: 1.0 degrees)

    Returns
    -------
    bool
        True if longitude appears to be cyclic (spans ~360°)
    """
    if len(glon) < 3:  # Need at least 3 points for meaningful cyclicity
        return False

    lon_range = glon[-1] - glon[0]
    dlon = np.mean(np.diff(glon))

    # Adaptive tolerance based on grid spacing
    # For coarse grids (large dlon), use larger tolerance
    adaptive_tolerance = max(tolerance, abs(dlon) * 0.5)

    # Check if the range plus one grid spacing is approximately 360°
    expected_range = lon_range + dlon
    is_cyclic_360 = abs(expected_range - 360.0) < adaptive_tolerance

    # Also check if range is already ~360° (for grids that include both 0 and 360)
    is_already_360 = abs(lon_range - 360.0) < adaptive_tolerance

    return is_cyclic_360 or is_already_360


def _ensure_south_to_north(
    z: np.ndarray, glat: np.ndarray, dim_order: str
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Ensure latitude dimension is ordered south-to-north as required by Fortran code.

    Parameters
    ----------
    z : ndarray
        Geopotential height data
    glat : ndarray
        Latitude coordinates
    dim_order : str
        Dimension order string

    Returns
    -------
    z_ordered : ndarray
        Data with latitude ordered south-to-north
    glat_ordered : ndarray
        Latitude coordinates ordered south-to-north
    """
    # Check if latitude needs to be reversed (north-to-south -> south-to-north)
    if glat[0] > glat[-1]:  # Currently north-to-south, need to reverse
        glat_ordered = glat[::-1].copy()

        # Find latitude axis in the data
        lat_axis = dim_order.lower().find("y")
        if lat_axis == -1:
            raise ValueError("Latitude dimension 'y' not found in dim_order")

        # Reverse the latitude dimension in the data
        z_ordered = np.flip(z, axis=lat_axis)
    else:
        # Already south-to-north or same latitude
        z_ordered = z
        glat_ordered = glat

    return z_ordered, glat_ordered


def geostrophic_wind(
    z: np.ndarray,
    glon: np.ndarray,
    glat: np.ndarray,
    dim_order: str,
    missing_value: float = -999.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate geostrophic wind components from geopotential height.

    This function can handle various input shapes by using windspharm's prep_data
    and recover_data utilities to reshape data for batch processing.

    Parameters
    ----------
    z : ndarray
        Geopotential height data [gpm]. Can be 2D, 3D, or 4D.
        Must contain latitude ('y') and longitude ('x') dimensions.
    glon : ndarray, shape (nlon,)
        Longitude coordinates in degrees
    glat : ndarray, shape (nlat,)
        Latitude coordinates in degrees (south to north)
    dim_order : str
        String specifying dimension order using:
        - 'x' for longitude
        - 'y' for latitude
        - 't' for time
        - 'z' for level
        Example: 'tzyx' for (time, level, lat, lon)
    missing_value : float, optional
        Missing value identifier (default: -999.0)

    Returns
    -------
    ug : ndarray
        Zonal geostrophic wind component [m/s] (same shape as input z)
    vg : ndarray
        Meridional geostrophic wind component [m/s] (same shape as input z)

    Examples
    --------
    # 2D case: single time/level
    >>> z2d = np.random.randn(73, 144)  # (lat, lon)
    >>> ug, vg = geostrophic_wind(z2d, glon, glat, 'yx')

    # 3D case: multiple times
    >>> z3d = np.random.randn(73, 144, 12)  # (lat, lon, time)
    >>> ug, vg = geostrophic_wind(z3d, glon, glat, 'yxt')

    # 4D case: multiple levels and times
    >>> z4d = np.random.randn(12, 17, 73, 144)  # (time, level, lat, lon)
    >>> ug, vg = geostrophic_wind(z4d, glon, glat, 'tzyx')

    # Alternative 4D ordering
    >>> z4d_alt = np.random.randn(73, 144, 17, 12)  # (lat, lon, level, time)
    >>> ug, vg = geostrophic_wind(z4d_alt, glon, glat, 'yxzt')
    """
    # Auto-detect longitude cyclicity
    cyclic = _is_longitude_cyclic(glon)
    iopt = 1 if cyclic else 0

    # Ensure latitude is ordered south-to-north as required by Fortran code
    z, glat = _ensure_south_to_north(z, glat, dim_order)

    # Handle multi-dimensional data using windspharm tools
    if len(z.shape) > 2:
        return _geostrophic_wind_multidim(z, glon, glat, dim_order, missing_value, iopt)
    else:
        # Direct 2D geostrophic wind calculation (ug, vg components)
        return _calc_geostrophic_2d(z, glon, glat, missing_value, iopt)


def _geostrophic_wind_multidim(
    z: np.ndarray,
    glon: np.ndarray,
    glat: np.ndarray,
    dim_order: str,
    missing_value: float,
    iopt: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Internal function for multi-dimensional geostrophic wind calculation.
    Uses windspharm prep_data/recover_data for dimension handling.
    """
    # Step 1: Prepare data (reshape to (nlat, nlon, combined_other_dims))
    prepared_z, recovery_info = prep_data(z, dim_order)
    nlat, nlon, n_combined = prepared_z.shape

    # Ensure coordinate arrays are compatible
    if len(glat) != nlat:
        raise ValueError(
            f"Latitude array length ({len(glat)}) doesn't match data ({nlat})"
        )
    if len(glon) != nlon:
        raise ValueError(
            f"Longitude array length ({len(glon)}) doesn't match data ({nlon})"
        )

    # Step 2: Use 3D function for all cases (handles any n_combined size)
    prepared_z = np.asarray(prepared_z, dtype=np.float32)
    ug_prepared, vg_prepared = z2geouv_3d(
        prepared_z, zmsg=missing_value, glon=glon, glat=glat, iopt=iopt
    )

    # Step 3: Recover original shape and dimension order
    ug_final = recover_data(ug_prepared, recovery_info)
    vg_final = recover_data(vg_prepared, recovery_info)

    return ug_final, vg_final


def _calc_geostrophic_2d(z, glon, glat, missing_value, iopt):
    """Calculate 2D geostrophic wind components using Fortran backend."""
    return z2geouv(
        np.asarray(z, dtype=np.float32),
        zmsg=missing_value,
        glon=glon,
        glat=glat,
        iopt=iopt,
    )


class GeostrophicWind:
    """
    Class-based interface for geostrophic wind calculations.

    This class provides a high-level interface similar to windspharm's VectorWind,
    allowing for easy calculation of various geostrophic wind quantities.

    Parameters
    ----------
    z : ndarray
        Geopotential height data [gpm]
    glon : ndarray
        Longitude coordinates [degrees]
    glat : ndarray
        Latitude coordinates [degrees] (south to north)
    dim_order : str
        Dimension ordering specification
    missing_value : float, optional
        Missing value identifier (default: -999.0)

    Examples
    --------
    >>> # Create GeostrophicWind instance (longitude cyclicity auto-detected)
    >>> gw = GeostrophicWind(z, glon, glat, 'tzyx')
    >>>
    >>> # Get wind components
    >>> ug, vg = gw.uv_components()
    >>>
    >>> # Calculate derived quantities
    >>> speed = gw.speed()
    >>>
    >>> # Access original data
    >>> z_orig = gw.geopotential_height
    """

    def __init__(
        self,
        z: np.ndarray,
        glon: np.ndarray,
        glat: np.ndarray,
        dim_order: str,
        missing_value: float = -999.0,
    ):
        self._z_original = np.asarray(z)
        self._glon = np.asarray(glon)
        self._glat = np.asarray(glat)
        self._dim_order = dim_order
        self._missing_value = missing_value
        # Calculate winds on initialization
        self._ug, self._vg = geostrophic_wind(z, glon, glat, dim_order, missing_value)

    @property
    def geopotential_height(self) -> np.ndarray:
        """Original geopotential height data."""
        return self._z_original

    @property
    def longitude(self) -> np.ndarray:
        """Longitude coordinates."""
        return self._glon

    @property
    def latitude(self) -> np.ndarray:
        """Latitude coordinates."""
        return self._glat

    def uv_components(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Return zonal and meridional wind components.

        Returns
        -------
        ug : ndarray
            Zonal (eastward) geostrophic wind component [m/s]
        vg : ndarray
            Meridional (northward) geostrophic wind component [m/s]
        """
        return self._ug, self._vg

    def speed(self) -> np.ndarray:
        """
        Calculate geostrophic wind speed.

        Returns
        -------
        speed : ndarray
            Geostrophic wind speed [m/s]
        """
        # Handle missing values properly
        ug_valid = np.where(self._ug == self._missing_value, 0, self._ug)
        vg_valid = np.where(self._vg == self._missing_value, 0, self._vg)

        speed = np.sqrt(ug_valid**2 + vg_valid**2)

        # Restore missing values where either component was missing
        missing_mask = (self._ug == self._missing_value) | (
            self._vg == self._missing_value
        )
        speed = np.where(missing_mask, self._missing_value, speed)

        return speed


# Convenience functions matching windspharm naming conventions


def geostrophic_uv(z, glon, glat, dim_order, **kwargs):
    """
    Calculate geostrophic wind components directly.

    This function calculates geostrophic wind components (ug, vg) from
    geopotential height fields. Uses the same implementation as GeostrophicWind
    class for consistency.

    Parameters
    ----------
    z : ndarray
        Geopotential height data [gpm]. Can be 2D, 3D, or 4D.
        Must contain latitude ('y') and longitude ('x') dimensions.
    glon : ndarray, shape (nlon,)
        Longitude coordinates in degrees
    glat : ndarray, shape (nlat,)
        Latitude coordinates in degrees (automatically ordered south-to-north)
    dim_order : str
        String specifying dimension order using:
        - 'x' for longitude
        - 'y' for latitude
        - 't' for time
        - 'z' for level
        Example: 'tzyx' for (time, level, lat, lon)
    missing_value : float, optional
        Missing value identifier (default: -999.0)

    Returns
    -------
    ug : ndarray
        Zonal geostrophic wind component [m/s] (same shape as input z)
    vg : ndarray
        Meridional geostrophic wind component [m/s] (same shape as input z)

    Notes
    -----
    - Longitude cyclicity is automatically detected
    - Latitude ordering is automatically ensured to be south-to-north
    - Uses optimized Fortran backend with SIMD optimization
    """
    gw = GeostrophicWind(z, glon, glat, dim_order, **kwargs)
    return gw.uv_components()


def geostrophic_speed(z, glon, glat, dim_order, **kwargs):
    """
    Calculate geostrophic wind speed directly.

    This function calculates geostrophic wind speed from geopotential height
    fields. Uses the GeostrophicWind class internally for consistent results.

    Parameters
    ----------
    z : ndarray
        Geopotential height data [gpm]. Can be 2D, 3D, or 4D.
        Must contain latitude ('y') and longitude ('x') dimensions.
    glon : ndarray, shape (nlon,)
        Longitude coordinates in degrees
    glat : ndarray, shape (nlat,)
        Latitude coordinates in degrees (automatically ordered south-to-north)
    dim_order : str
        String specifying dimension order using:
        - 'x' for longitude
        - 'y' for latitude
        - 't' for time
        - 'z' for level
        Example: 'tzyx' for (time, level, lat, lon)
    missing_value : float, optional
        Missing value identifier (default: -999.0)

    Returns
    -------
    speed : ndarray
        Geostrophic wind speed [m/s] (same shape as input z)

    Notes
    -----
    - Longitude cyclicity is automatically detected
    - Latitude ordering is automatically ensured to be south-to-north
    - Speed calculated as sqrt(ug² + vg²)
    - Missing values are properly handled
    """
    gw = GeostrophicWind(z, glon, glat, dim_order, **kwargs)
    return gw.speed()
