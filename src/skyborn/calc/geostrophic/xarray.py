"""
Geostrophic wind calculation for xarray DataArrays.

This module provides functions to calculate geostrophic wind components from
geopotential height fields using xarray DataArrays. It automatically detects
spatial coordinates and preserves coordinate information and metadata throughout
the computation process.

Main Functions:
    geostrophic_wind : Calculate geostrophic wind components for xarray DataArray
    GeostrophicWind : Class-based interface for xarray DataArrays

Examples:
    >>> import xarray as xr
    >>> import numpy as np
    >>> from skyborn.calc.geostrophic.xarray import geostrophic_wind
    >>>
    >>> # Load geopotential height data
    >>> z = xr.open_dataarray('geopotential_500hPa.nc')
    >>> result = geostrophic_wind(z)  # Auto-detects coordinates
    >>> print(result.ug.attrs)  # Original attributes preserved
"""

from __future__ import annotations

import warnings
from typing import Any, Dict, List, Optional, Tuple, Union

__all__ = ["geostrophic_wind", "GeostrophicWind"]

import numpy as np
import xarray as xr

from . import interface

# Type aliases
DataArray = xr.DataArray


def _detect_spatial_dimensions(
    data_array: DataArray,
) -> Tuple[Optional[int], Optional[int]]:
    """
    Auto-detect longitude and latitude dimension indices in xarray DataArray.

    Parameters
    ----------
    data_array : xr.DataArray
        Input data to analyze

    Returns
    -------
    xdim : int, optional
        Longitude dimension index (None if not found)
    ydim : int, optional
        Latitude dimension index (None if not found)

    Raises
    ------
    ValueError
        If both longitude and latitude dimensions cannot be identified
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

    xdim = ydim = None

    for i, dim_name in enumerate(dims):
        dim_lower = dim_name.lower()

        if any(name.lower() in dim_lower for name in lon_names):
            xdim = i
        elif any(name.lower() in dim_lower for name in lat_names):
            ydim = i

    # Both longitude and latitude are required for geostrophic wind calculation
    if xdim is None or ydim is None:
        raise ValueError(
            f"Could not auto-detect both longitude and latitude dimensions. "
            f"Found dims: {dims}. Expected longitude and latitude coordinates."
        )

    return xdim, ydim


def _extract_coordinates(data_array: DataArray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract longitude and latitude coordinate arrays from xarray DataArray.

    Parameters
    ----------
    data_array : xr.DataArray
        Input data containing coordinate information

    Returns
    -------
    glon : np.ndarray
        Longitude coordinates in degrees
    glat : np.ndarray
        Latitude coordinates in degrees

    Raises
    ------
    ValueError
        If longitude or latitude coordinates are not found
    """
    dims = data_array.dims
    coords = data_array.coords

    # Find longitude and latitude coordinate names
    lon_coord = lat_coord = None

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

    for coord_name in coords:
        coord_lower = coord_name.lower()
        if any(name.lower() in coord_lower for name in lon_names):
            lon_coord = coord_name
        elif any(name.lower() in coord_lower for name in lat_names):
            lat_coord = coord_name

    if lon_coord is None or lat_coord is None:
        raise ValueError(
            f"Could not find longitude and latitude coordinates. "
            f"Available coordinates: {list(coords.keys())}"
        )

    # Extract coordinate values
    glon = coords[lon_coord].values
    glat = coords[lat_coord].values

    return glon, glat


def _create_dim_order_string(data_array: DataArray, xdim: int, ydim: int) -> str:
    """
    Create dimension order string for interface functions.

    Parameters
    ----------
    data_array : xr.DataArray
        Input data array
    xdim : int
        Longitude dimension index
    ydim : int
        Latitude dimension index

    Returns
    -------
    dim_order : str
        Dimension order string (e.g., 'tzyx', 'yx', etc.)
    """
    dims = list(data_array.dims)
    dim_order = [""] * len(dims)

    # Set longitude and latitude
    dim_order[xdim] = "x"
    dim_order[ydim] = "y"

    # Common patterns for other dimensions
    time_names = {"time", "t", "T", "year", "month", "yr", "mn", "season"}
    level_names = {
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
        "height",
        "altitude",
        "isobaric",
    }

    # Fill in other dimensions
    for i, dim_name in enumerate(dims):
        if dim_order[i] == "":  # Not yet assigned
            dim_lower = dim_name.lower()
            if any(name.lower() in dim_lower for name in time_names):
                dim_order[i] = "t"
            elif any(name.lower() in dim_lower for name in level_names):
                dim_order[i] = "z"
            else:
                # Default to 't' for unrecognized dimensions
                dim_order[i] = "t"

    return "".join(dim_order)


def geostrophic_wind(
    z: DataArray,
    missing_value: float = -999.0,
    keep_attrs: bool = True,
) -> xr.Dataset:
    """
    Calculate geostrophic wind components for xarray DataArrays.

    This function processes geopotential height data to calculate geostrophic
    wind components (ug, vg). It automatically detects coordinate dimensions
    and preserves all metadata.

    Parameters
    ----------
    z : xarray.DataArray
        Geopotential height data [gpm]. Can be 2D, 3D, or 4D.
        Must contain longitude and latitude dimensions with coordinate information.
    missing_value : float, optional
        Missing value identifier (default: -999.0)
    keep_attrs : bool, optional
        Preserve input DataArray attributes in output (default: True)

    Returns
    -------
    xarray.Dataset
        Dataset containing geostrophic wind components:
        - 'ug': Zonal geostrophic wind component [m/s] with spatial/temporal coordinates
        - 'vg': Meridional geostrophic wind component [m/s] with spatial/temporal coordinates

    Examples
    --------
    **2D Geopotential Height Analysis:**

    >>> import xarray as xr
    >>> import numpy as np
    >>> from skyborn.calc.geostrophic.xarray import geostrophic_wind
    >>>
    >>> # Load 500 hPa geopotential height
    >>> z = xr.open_dataarray('z500_era5.nc')  # Shape: (lat, lon)
    >>> result = geostrophic_wind(z)
    >>> print(f"Wind components: ug{result.ug.shape}, vg{result.vg.shape}")

    **3D Time Series Analysis:**

    >>> # Multi-time geopotential height data
    >>> z_3d = xr.open_dataarray('z500_timeseries.nc')  # Shape: (time, lat, lon)
    >>> result = geostrophic_wind(z_3d)
    >>> # Result preserves time dimension: (time, lat, lon)
    >>> monthly_mean = result.ug.groupby('time.month').mean()

    **4D Multi-level Analysis:**

    >>> # Multi-level, multi-time data
    >>> z_4d = xr.open_dataarray('z_multilevel.nc')  # Shape: (time, level, lat, lon)
    >>> result = geostrophic_wind(z_4d)
    >>> # Result shape: (time, level, lat, lon)
    >>> surface_winds = result.sel(level=1000)  # 1000 hPa level

    **Simplified Interface (No coordinate specification needed):**

    >>> # Automatic coordinate detection
    >>> result = geostrophic_wind(z_data)  # Longitude cyclicity auto-detected
    >>> print(f"Longitude cyclicity auto-detected: {result.attrs['longitude_cyclic']}")
    >>> print(f"Latitude ordering: {result.attrs['latitude_ordering']}")

    Notes
    -----
    - Longitude cyclicity is automatically detected from coordinate spacing
    - Latitude ordering is automatically ensured to be south-to-north as required
    - Requires compiled Fortran extensions for optimal performance
    - All coordinate information and attributes are preserved

    The function automatically:
    - Detects longitude and latitude coordinates using metadata
    - Handles missing values (NaN or masked arrays)
    - Preserves all coordinate information and attributes
    - Works with multi-dimensional data of any supported shape

    See Also
    --------
    skyborn.calc.geostrophic.interface.geostrophic_wind : Lower-level function for numpy arrays
    GeostrophicWind : Class-based interface for xarray DataArrays
    """
    # Validate input type
    if not isinstance(z, xr.DataArray):
        raise TypeError(f"z must be xarray.DataArray, got {type(z).__name__}")

    # Auto-detect spatial dimensions
    xdim, ydim = _detect_spatial_dimensions(z)

    # Extract coordinate arrays
    glon, glat = _extract_coordinates(z)

    # Create dimension order string
    dim_order = _create_dim_order_string(z, xdim, ydim)

    # Store original coordinate information
    original_coords = z.coords
    original_dims = z.dims

    # Extract numpy array
    z_data = z.values

    # Call the core geostrophic calculation function
    ug_data, vg_data = interface.geostrophic_wind(
        z_data, glon, glat, dim_order, missing_value=missing_value
    )

    # Create output coordinates (same as input)
    output_coords = {}
    for dim_name in z.dims:
        if dim_name in z.coords:
            output_coords[dim_name] = z.coords[dim_name]

    # Create DataArrays for wind components
    ug = xr.DataArray(
        ug_data,
        dims=z.dims,
        coords=output_coords,
        attrs={
            "long_name": "Zonal geostrophic wind component",
            "units": "m s-1",
            "standard_name": "eastward_geostrophic_wind",
            "description": "Zonal (eastward) component of geostrophic wind calculated from geopotential height",
        },
    )

    vg = xr.DataArray(
        vg_data,
        dims=z.dims,
        coords=output_coords,
        attrs={
            "long_name": "Meridional geostrophic wind component",
            "units": "m s-1",
            "standard_name": "northward_geostrophic_wind",
            "description": "Meridional (northward) component of geostrophic wind calculated from geopotential height",
        },
    )

    # Create Dataset
    ds = xr.Dataset({"ug": ug, "vg": vg})

    # Add global attributes
    ds.attrs = {
        "title": "Geostrophic wind calculation results",
        "description": "Geostrophic wind components calculated from geopotential height",
        "longitude_cyclic": interface._is_longitude_cyclic(glon),
        "latitude_ordering": "south_to_north",
        "missing_value": missing_value,
        "method": "Finite difference approximation with geostrophic balance",
        "software": "skyborn atmospheric calculation package",
        "equations": "ug = -(g/f)*dZ/dy, vg = (g/f)*dZ/dx",
    }

    # Preserve original attributes if requested
    if keep_attrs and hasattr(z, "attrs") and z.attrs:
        ds.attrs.update({f"source_geopotential_{k}": v for k, v in z.attrs.items()})

    return ds


class GeostrophicWind:
    """
    Class-based geostrophic wind analysis using xarray DataArrays.

    This class provides a high-level interface for geostrophic wind calculations
    that preserves xarray coordinate information and metadata. It wraps the standard
    interface implementation while maintaining CF-compliant attributes.

    Parameters
    ----------
    z : xarray.DataArray
        Geopotential height data [gpm]. Must contain longitude and latitude
        dimensions with appropriate coordinate information.
    missing_value : float, optional
        Missing value identifier (default: -999.0)
    keep_attrs : bool, optional
        Preserve input DataArray attributes in output (default: True)

    Attributes
    ----------
    _z_original : xarray.DataArray
        Original geopotential height data
    _result : xarray.Dataset
        Computed geostrophic wind components
    _glon : np.ndarray
        Longitude coordinates
    _glat : np.ndarray
        Latitude coordinates

    Examples
    --------
    >>> import xarray as xr
    >>> from skyborn.calc.geostrophic.xarray import GeostrophicWind
    >>>
    >>> # Load geopotential height
    >>> z = xr.open_dataarray('z500.nc')
    >>>
    >>> # Create GeostrophicWind instance
    >>> gw = GeostrophicWind(z)
    >>>
    >>> # Get wind components with preserved metadata
    >>> ug, vg = gw.uv_components()
    >>> print(ug.attrs)  # CF-compliant attributes
    >>>
    >>> # Calculate derived quantities
    >>> speed = gw.speed()
    >>> print(f"Max wind speed: {float(speed.max()):.1f} m/s")
    >>>
    >>> # Access original data
    >>> z_orig = gw.geopotential_height
    """

    def __init__(
        self,
        z: DataArray,
        missing_value: float = -999.0,
        keep_attrs: bool = True,
    ):
        """Initialize GeostrophicWind with xarray DataArray."""
        self._z_original = z
        self._missing_value = missing_value
        self._keep_attrs = keep_attrs

        # Extract coordinates for later use
        self._glon, self._glat = _extract_coordinates(z)

        # Calculate geostrophic wind components
        self._result = geostrophic_wind(
            z, missing_value=missing_value, keep_attrs=keep_attrs
        )

    @property
    def geopotential_height(self) -> DataArray:
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

    def uv_components(self) -> Tuple[DataArray, DataArray]:
        """
        Return zonal and meridional wind components.

        Returns
        -------
        ug : xarray.DataArray
            Zonal (eastward) geostrophic wind component [m/s]
        vg : xarray.DataArray
            Meridional (northward) geostrophic wind component [m/s]
        """
        return self._result.ug, self._result.vg

    def speed(self) -> DataArray:
        """
        Calculate geostrophic wind speed.

        Returns
        -------
        speed : xarray.DataArray
            Geostrophic wind speed [m/s]
        """
        ug, vg = self.uv_components()

        # Handle missing values properly
        ug_valid = ug.where(ug != self._missing_value, 0)
        vg_valid = vg.where(vg != self._missing_value, 0)

        # Calculate speed
        speed = np.hypot(ug_valid, vg_valid)

        # Restore missing values where either component was missing
        missing_mask = (ug == self._missing_value) | (vg == self._missing_value)
        speed = speed.where(~missing_mask, self._missing_value)

        # Set attributes
        speed.attrs = {
            "long_name": "Geostrophic wind speed",
            "units": "m s-1",
            "standard_name": "geostrophic_wind_speed",
            "description": "Speed of geostrophic wind calculated as sqrt(ug² + vg²)",
        }

        return speed
