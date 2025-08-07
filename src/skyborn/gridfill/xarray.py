"""
Gridfill procedures for missing value interpolation with xarray interface.

This module provides functions to fill missing values in xarray DataArrays using
iterative relaxation methods to solve Poisson's equation. It preserves coordinate
information and metadata throughout the computation process.

Main Functions:
    fill : Fill missing values in xarray DataArray

Examples:
    >>> import xarray as xr
    >>> import numpy as np
    >>> from skyborn.gridfill.xarray import fill
    >>>
    >>> # Load data with missing values
    >>> data = xr.open_dataarray('temperature_with_gaps.nc')
    >>>
    >>> # Fill missing values preserving metadata
    >>> filled_data = fill(data, eps=1e-4)
    >>> print(filled_data.attrs)  # Original attributes preserved
"""

from __future__ import annotations
from typing import Optional, Tuple, Union, Any, Dict, List
import warnings

try:
    import xarray as xr
except ImportError:
    try:
        import xray as xr  # Fallback for older versions
    except ImportError:
        raise ImportError(
            "xarray is required for the xarray interface. "
            "Install with: pip install xarray"
        )

import numpy as np
import numpy.ma as ma

from . import gridfill

# Type aliases for better readability
DataArray = xr.DataArray


def _find_spatial_coordinates(
    data: DataArray,
) -> Tuple[str, str, int, int]:
    """
    Find spatial coordinate dimensions in xarray DataArray.

    This function automatically detects latitude/longitude or y/x coordinate
    dimensions in the DataArray using common naming conventions and coordinate
    attributes.

    Parameters
    ----------
    data : xarray.DataArray
        Input DataArray to analyze

    Returns
    -------
    y_name : str
        Name of the y-coordinate dimension (latitude)
    x_name : str
        Name of the x-coordinate dimension (longitude)
    y_dim : int
        Index of y-coordinate dimension
    x_dim : int
        Index of x-coordinate dimension

    Raises
    ------
    ValueError
        If spatial coordinates cannot be identified

    Notes
    -----
    Detection priority:
    1. Standard names: 'latitude'/'longitude'
    2. Axis attributes: axis='Y'/'X'
    3. Common dimension names: 'lat'/'lon', 'y'/'x'
    4. Unit attributes: 'degrees_north'/'degrees_east'
    """
    x_name = None
    y_name = None

    # Try to find coordinates by standard_name first
    for name, coord in data.coords.items():
        if hasattr(coord, "standard_name"):
            if coord.standard_name == "latitude":
                y_name = name
            elif coord.standard_name == "longitude":
                x_name = name

    # Try to find by axis attribute
    if x_name is None or y_name is None:
        for name, coord in data.coords.items():
            if hasattr(coord, "axis"):
                if coord.axis == "Y" and y_name is None:
                    y_name = name
                elif coord.axis == "X" and x_name is None:
                    x_name = name

    # Try to find by common dimension names
    if x_name is None or y_name is None:
        for name in data.dims:
            name_lower = name.lower()
            if y_name is None and name_lower in ["lat", "latitude", "y"]:
                y_name = name
            elif x_name is None and name_lower in ["lon", "lng", "longitude", "x"]:
                x_name = name

    # Try to find by units
    if x_name is None or y_name is None:
        for name, coord in data.coords.items():
            if hasattr(coord, "units"):
                units = getattr(coord, "units", "")
                if y_name is None and "degrees_north" in str(units):
                    y_name = name
                elif x_name is None and "degrees_east" in str(units):
                    x_name = name

    # Validate that we found both coordinates
    if x_name is None or y_name is None:
        available_dims = list(data.dims)
        available_coords = list(data.coords.keys())
        raise ValueError(
            f"Could not identify spatial coordinates automatically. "
            f"Available dimensions: {available_dims}, "
            f"Available coordinates: {available_coords}. "
            f"Please ensure your data has recognizable latitude/longitude "
            f"coordinates with appropriate metadata (standard_name, axis, or units attributes)."
        )

    # Get dimension indices
    try:
        y_dim = data.dims.index(y_name)
        x_dim = data.dims.index(x_name)
    except ValueError as e:
        raise ValueError(f"Coordinate dimension not found in data dimensions: {e}")

    return y_name, x_name, y_dim, x_dim


def _detect_cyclic_longitude(lon_coord: xr.DataArray) -> bool:
    """
    Detect if longitude coordinate is cyclic (wraps around).

    Parameters
    ----------
    lon_coord : xarray.DataArray
        Longitude coordinate to analyze

    Returns
    -------
    bool
        True if longitude appears to be cyclic (global coverage)

    Notes
    -----
    Detection criteria:
    1. Check for 'circular' attribute (iris convention)
    2. Check if span is approximately 360 degrees
    3. Check if values span from approximately -180 to 180 or 0 to 360
    """
    # Check for explicit circular attribute
    if hasattr(lon_coord, "circular") and lon_coord.circular:
        return True

    # Get longitude values
    lon_vals = lon_coord.values

    # Check if span covers ~360 degrees
    lon_span = np.max(lon_vals) - np.min(lon_vals)
    if np.abs(lon_span - 360.0) < 10.0:  # Allow some tolerance
        return True

    # Check common global longitude ranges
    lon_min, lon_max = np.min(lon_vals), np.max(lon_vals)

    # Check for 0-360 range
    if np.abs(lon_min) < 10.0 and np.abs(lon_max - 360.0) < 10.0:
        return True

    # Check for -180 to 180 range
    if np.abs(lon_min + 180.0) < 10.0 and np.abs(lon_max - 180.0) < 10.0:
        return True

    return False


def fill(
    data: DataArray,
    eps: float,
    x_dim: Optional[str] = None,
    y_dim: Optional[str] = None,
    relax: float = 0.6,
    itermax: int = 100,
    initzonal: bool = False,
    cyclic: Optional[bool] = None,
    verbose: bool = False,
    keep_attrs: bool = True,
) -> DataArray:
    """
    Fill missing values in xarray DataArray using Poisson equation solver.

    This function fills missing values (NaN or masked values) in gridded data by
    solving Poisson's equation (∇²φ = 0) using an iterative relaxation scheme.
    The method provides smooth interpolation while preserving coordinate information
    and metadata from the input DataArray.

    Parameters
    ----------
    data : xarray.DataArray
        Input DataArray containing data with missing values to fill. Missing
        values can be NaN or masked values.
    eps : float
        Convergence tolerance. Iteration stops when the maximum residual
        falls below this threshold.
    x_dim : str, optional
        Name of the x-coordinate dimension (longitude). If None, will be
        detected automatically using coordinate metadata.
    y_dim : str, optional
        Name of the y-coordinate dimension (latitude). If None, will be
        detected automatically using coordinate metadata.
    relax : float, default 0.6
        Relaxation parameter for the iterative scheme. Must be in range
        (0, 1). Values between 0.45-0.6 typically work well.
    itermax : int, default 100
        Maximum number of iterations.
    initzonal : bool, default False
        Initialization method for missing values:
        - False: Initialize with zeros
        - True: Initialize with zonal (x-direction) mean
    cyclic : bool, optional
        Whether the x-coordinate is cyclic (e.g., longitude wrapping).
        If None, will be detected automatically for longitude coordinates.
    verbose : bool, default False
        Print convergence information for each slice.
    keep_attrs : bool, default True
        Preserve input DataArray attributes in output.

    Returns
    -------
    filled_data : xarray.DataArray
        DataArray with missing values filled, preserving coordinates and
        optionally attributes from input.

    Raises
    ------
    ValueError
        If spatial coordinates cannot be identified or are invalid
    TypeError
        If input is not an xarray DataArray

    Warnings
    --------
    Issues warning if algorithm fails to converge on any slices

    Notes
    -----
    The algorithm solves:
    ∇²φ = (∂²φ/∂x²) + (∂²φ/∂y²) = 0

    using a finite difference relaxation scheme. The method automatically:
    - Detects spatial coordinates using metadata
    - Handles cyclic longitude boundaries for global data
    - Preserves all coordinate information and attributes
    - Works with multi-dimensional data (time series, levels, etc.)

    For missing value detection, both NaN values and xarray/numpy masked
    arrays are supported.

    Examples
    --------
    Basic usage with automatic coordinate detection:

    >>> import xarray as xr
    >>> import numpy as np
    >>> from skyborn.gridfill.xarray import fill
    >>>
    >>> # Load data with missing values
    >>> data = xr.open_dataarray('sst_with_gaps.nc')
    >>>
    >>> # Fill missing values
    >>> filled = fill(data, eps=1e-4)
    >>> print(f"Original shape: {data.shape}")
    >>> print(f"Filled shape: {filled.shape}")
    >>> print(f"Attributes preserved: {filled.attrs == data.attrs}")

    Advanced usage with explicit parameters:

    >>> # Create test data with gaps
    >>> lons = np.linspace(0, 360, 72, endpoint=False)
    >>> lats = np.linspace(-90, 90, 36)
    >>> time = pd.date_range('2020-01-01', periods=12, freq='M')
    >>>
    >>> # Create DataArray with metadata
    >>> data = xr.DataArray(
    ...     np.random.rand(12, 36, 72),
    ...     coords={'time': time, 'lat': lats, 'lon': lons},
    ...     dims=['time', 'lat', 'lon'],
    ...     attrs={'units': 'K', 'long_name': 'temperature'}
    ... )
    >>>
    >>> # Add some missing values
    >>> data = data.where(np.random.rand(*data.shape) > 0.1)
    >>>
    >>> # Fill with custom settings
    >>> filled = fill(
    ...     data,
    ...     eps=1e-5,
    ...     relax=0.55,
    ...     initzonal=True,
    ...     verbose=True
    ... )

    Working with specific coordinate dimensions:

    >>> # Explicitly specify coordinate dimensions
    >>> filled = fill(data, eps=1e-4, x_dim='longitude', y_dim='latitude')

    See Also
    --------
    skyborn.gridfill.fill : Lower-level function for numpy arrays
    skyborn.gridfill.fill_cube : Fill iris cubes
    """
    # Validate input type
    if not isinstance(data, xr.DataArray):
        raise TypeError(f"data must be xarray.DataArray, got {type(data).__name__}")

    # Find spatial coordinates if not provided
    if x_dim is None or y_dim is None:
        y_name, x_name, y_dim_idx, x_dim_idx = _find_spatial_coordinates(data)
        if x_dim is None:
            x_dim = x_name
        if y_dim is None:
            y_dim = y_name
    else:
        # Validate provided dimension names
        if x_dim not in data.dims:
            raise ValueError(
                f"x_dim '{x_dim}' not found in data dimensions: {list(data.dims)}"
            )
        if y_dim not in data.dims:
            raise ValueError(
                f"y_dim '{y_dim}' not found in data dimensions: {list(data.dims)}"
            )
        x_dim_idx = data.dims.index(x_dim)
        y_dim_idx = data.dims.index(y_dim)

    # Detect cyclic boundary if not specified
    if cyclic is None:
        x_coord = data.coords[x_dim]
        cyclic = _detect_cyclic_longitude(x_coord)

        if verbose:
            print(f"Auto-detected cyclic={cyclic} for coordinate '{x_dim}'")

    # Convert to masked array for gridfill processing
    data_values = data.values

    # Handle different types of missing values
    if hasattr(data_values, "mask"):
        # Already a masked array
        masked_data = data_values
    else:
        # Create mask from NaN values
        mask = np.isnan(data_values)
        if not np.any(mask):
            # No missing values to fill
            warnings.warn("No missing values found in input data")
            return data if not keep_attrs else data.copy()
        masked_data = ma.array(data_values, mask=mask)

    # Call the core gridfill function
    filled_values, converged = gridfill.fill(
        masked_data,
        xdim=x_dim_idx,
        ydim=y_dim_idx,
        eps=eps,
        relax=relax,
        itermax=itermax,
        initzonal=initzonal,
        cyclic=cyclic,
        verbose=verbose,
    )

    # Check convergence and issue warnings
    not_converged = np.logical_not(converged)
    if np.any(not_converged):
        warnings.warn(
            f"gridfill did not converge on {not_converged.sum()} out of "
            f"{not_converged.size} slices. Consider increasing itermax or "
            f"relaxing eps tolerance."
        )

    # Create output DataArray preserving coordinates
    filled_data = xr.DataArray(
        filled_values,
        coords=data.coords,
        dims=data.dims,
        name=data.name,
    )

    # Preserve attributes if requested
    if keep_attrs:
        filled_data.attrs.update(data.attrs)
        # Add processing history
        if "history" not in filled_data.attrs:
            filled_data.attrs["history"] = ""
        filled_data.attrs[
            "history"
        ] += f"; Filled missing values using gridfill (eps={eps})"

    return filled_data


def fill_multiple(
    datasets: List[DataArray],
    eps: float,
    x_dim: Optional[str] = None,
    y_dim: Optional[str] = None,
    **kwargs,
) -> List[DataArray]:
    """
    Fill missing values in multiple DataArrays with consistent parameters.

    This convenience function applies the same gridfill parameters to multiple
    DataArrays, ensuring consistent processing across related datasets.

    Parameters
    ----------
    datasets : list of xarray.DataArray
        List of DataArrays to process
    eps : float
        Convergence tolerance for all datasets
    x_dim : str, optional
        X-coordinate dimension name (applied to all)
    y_dim : str, optional
        Y-coordinate dimension name (applied to all)
    **kwargs
        Additional parameters passed to fill()

    Returns
    -------
    list of xarray.DataArray
        List of filled DataArrays in same order as input

    Examples
    --------
    >>> from skyborn.gridfill.xarray import fill_multiple
    >>>
    >>> # Fill multiple related variables
    >>> temp_filled, humid_filled = fill_multiple(
    ...     [temperature_data, humidity_data],
    ...     eps=1e-4,
    ...     verbose=True
    ... )
    """
    return [
        fill(data, eps=eps, x_dim=x_dim, y_dim=y_dim, **kwargs) for data in datasets
    ]


def validate_grid_coverage(
    data: DataArray,
    x_dim: Optional[str] = None,
    y_dim: Optional[str] = None,
    min_coverage: float = 0.1,
) -> Dict[str, Any]:
    """
    Validate grid data coverage and suitability for gridfill.

    This function analyzes the input data to determine if it's suitable
    for gap filling and provides diagnostic information.

    Parameters
    ----------
    data : xarray.DataArray
        Input data to analyze
    x_dim : str, optional
        X-coordinate dimension name
    y_dim : str, optional
        Y-coordinate dimension name
    min_coverage : float, default 0.1
        Minimum fraction of valid data required (0.0 to 1.0)

    Returns
    -------
    dict
        Dictionary containing validation results:
        - 'valid': bool, whether data is suitable for filling
        - 'coverage': float, fraction of valid data points
        - 'total_points': int, total number of grid points
        - 'missing_points': int, number of missing points
        - 'messages': list, diagnostic messages

    Examples
    --------
    >>> from skyborn.gridfill.xarray import validate_grid_coverage
    >>>
    >>> # Check data quality before filling
    >>> validation = validate_grid_coverage(data, min_coverage=0.2)
    >>> if validation['valid']:
    ...     filled = fill(data, eps=1e-4)
    ... else:
    ...     print("Data quality issues:", validation['messages'])
    """
    # Find spatial coordinates if not provided
    if x_dim is None or y_dim is None:
        y_name, x_name, y_dim_idx, x_dim_idx = _find_spatial_coordinates(data)
        if x_dim is None:
            x_dim = x_name
        if y_dim is None:
            y_dim = y_name

    messages = []

    # Calculate coverage statistics
    data_values = data.values
    if hasattr(data_values, "mask"):
        missing_mask = data_values.mask
    else:
        missing_mask = np.isnan(data_values)

    total_points = missing_mask.size
    missing_points = np.sum(missing_mask)
    valid_points = total_points - missing_points
    coverage = valid_points / total_points if total_points > 0 else 0.0

    # Validate coverage
    valid = True
    if coverage < min_coverage:
        valid = False
        messages.append(
            f"Insufficient data coverage: {coverage:.1%} < {min_coverage:.1%}"
        )

    # Check for completely empty slices
    if data.ndim > 2:
        # For multi-dimensional data, check each 2D slice
        other_dims = [dim for dim in data.dims if dim not in [x_dim, y_dim]]
        if other_dims:
            slice_coverages = []
            for coords in data.groupby(other_dims[0]):
                slice_data = coords[1]
                slice_missing = (
                    np.isnan(slice_data.values)
                    if not hasattr(slice_data.values, "mask")
                    else slice_data.values.mask
                )
                slice_coverage = 1.0 - (np.sum(slice_missing) / slice_missing.size)
                slice_coverages.append(slice_coverage)

            min_slice_coverage = np.min(slice_coverages)
            if min_slice_coverage < min_coverage:
                messages.append(
                    f"Some slices have insufficient coverage (min: {min_slice_coverage:.1%})"
                )

    # Check coordinate regularity
    x_coord = data.coords[x_dim]
    y_coord = data.coords[y_dim]

    # Check for regular spacing
    if len(x_coord) > 1:
        x_diffs = np.diff(x_coord.values)
        if not np.allclose(x_diffs, x_diffs[0], rtol=1e-3):
            messages.append("X-coordinate spacing is not regular")

    if len(y_coord) > 1:
        y_diffs = np.diff(y_coord.values)
        if not np.allclose(y_diffs, y_diffs[0], rtol=1e-3):
            messages.append("Y-coordinate spacing is not regular")

    return {
        "valid": valid,
        "coverage": coverage,
        "total_points": total_points,
        "missing_points": missing_points,
        "messages": messages,
    }
