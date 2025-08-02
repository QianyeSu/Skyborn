"""
Common functionality shared across windspharm interfaces.

This module provides utility functions that are used by multiple windspharm
interfaces (standard, xarray, etc.) for grid type detection, dimension
ordering, and data reshaping operations.

Main Functions:
    get_apiorder: Calculate dimension reordering for API compatibility
    inspect_gridtype: Determine grid type from latitude coordinates
    to3d: Reshape array to 3D format for spherical harmonic analysis
"""

from __future__ import annotations
from typing import Tuple, List, Union, Literal
import numpy as np
from skyborn.spharm import gaussian_lats_wts

# Type aliases for better readability
ArrayLike = Union[np.ndarray, np.ma.MaskedArray]
GridType = Literal["regular", "gaussian"]


def get_apiorder(
    ndim: int, latitude_dim: int, longitude_dim: int
) -> Tuple[List[int], List[int]]:
    """
    Calculate dimension ordering for API-compatible array transposition.

    This function determines the proper dimension ordering to move latitude
    and longitude dimensions to the first two positions, which is required
    by the spherical harmonic analysis routines.

    Parameters
    ----------
    ndim : int
        Total number of dimensions in the array
    latitude_dim : int
        Index of the latitude dimension
    longitude_dim : int
        Index of the longitude dimension

    Returns
    -------
    apiorder : list of int
        List of dimension indices in the order required for API compatibility.
        Latitude and longitude dimensions are moved to positions 0 and 1.
    reorder : list of int
        Inverse mapping to restore original dimension order after processing

    Examples
    --------
    >>> # For 4D array (time, level, lat, lon)
    >>> apiorder, reorder = get_apiorder(4, 2, 3)
    >>> print(apiorder)  # [2, 3, 0, 1] - lat, lon, time, level
    >>> print(reorder)   # [2, 3, 0, 1] - inverse mapping
    """
    apiorder = list(range(ndim))
    apiorder.remove(latitude_dim)
    apiorder.remove(longitude_dim)
    apiorder.insert(0, latitude_dim)
    apiorder.insert(1, longitude_dim)
    reorder = [apiorder.index(i) for i in range(ndim)]
    return apiorder, reorder


def inspect_gridtype(latitudes: ArrayLike) -> GridType:
    """
    Determine grid type by examining latitude coordinate values.

    Analyzes the spacing and distribution of latitude points to determine
    whether they correspond to a regular (equally-spaced) or Gaussian grid.
    This is essential for selecting the appropriate spherical harmonic
    analysis method.

    Parameters
    ----------
    latitudes : array_like
        Array of latitude coordinate values in degrees. Should be ordered
        from north to south (90° to -90°).

    Returns
    -------
    gridtype : {'regular', 'gaussian'}
        Detected grid type:
        - 'regular': Equally-spaced latitude grid
        - 'gaussian': Gaussian latitude grid (optimal for spectral methods)

    Raises
    ------
    ValueError
        If the grid type cannot be determined or if latitudes don't match
        either regular or Gaussian patterns within tolerance

    Examples
    --------
    >>> # Regular grid with 73 points including poles
    >>> lats_regular = np.linspace(90, -90, 73)
    >>> gridtype = inspect_gridtype(lats_regular)
    >>> print(gridtype)  # 'regular'
    >>>
    >>> # Gaussian grid
    >>> from skyborn.spharm import gaussian_lats_wts
    >>> lats_gauss, _ = gaussian_lats_wts(64)
    >>> gridtype = inspect_gridtype(lats_gauss)
    >>> print(gridtype)  # 'gaussian'
    """
    # Convert to numpy array for consistent handling
    latitudes = np.asarray(latitudes)

    # Define tolerance for floating-point comparisons
    # This must be much smaller than typical grid spacings
    tolerance = 5e-4

    # Get the number of latitude points
    nlat = len(latitudes)
    if nlat < 2:
        raise ValueError(
            f"Need at least 2 latitude points for grid type detection, got {nlat}"
        )

    # Check if latitudes are equally spaced
    diffs = np.abs(np.diff(latitudes))
    equally_spaced = (np.abs(diffs - diffs[0]) < tolerance).all()

    if not equally_spaced:
        # Not equally-spaced - check if they match Gaussian latitudes
        try:
            gauss_reference, _ = gaussian_lats_wts(nlat)
            difference = np.abs(latitudes - gauss_reference)
            if (difference > tolerance).any():
                raise ValueError(
                    f"Latitudes are neither equally-spaced nor Gaussian. "
                    f"Maximum difference from Gaussian: {difference.max():.6f} degrees, "
                    f"tolerance: {tolerance:.6f} degrees"
                )
            gridtype = "gaussian"
        except Exception as e:
            raise ValueError(
                f"Failed to generate Gaussian reference latitudes for {nlat} points: {e}"
            ) from e
    else:
        # Equally-spaced - verify they match global regular grid pattern
        if nlat % 2:
            # Odd number of latitudes includes both poles
            equal_reference = np.linspace(90, -90, nlat)
        else:
            # Even number of latitudes excludes poles
            delta_latitude = 180.0 / nlat
            equal_reference = np.linspace(
                90 - 0.5 * delta_latitude, -90 + 0.5 * delta_latitude, nlat
            )

        difference = np.abs(latitudes - equal_reference)
        if (difference > tolerance).any():
            raise ValueError(
                f"Equally-spaced latitudes don't match global regular grid pattern. "
                f"Maximum difference: {difference.max():.6f} degrees, "
                f"tolerance: {tolerance:.6f} degrees. "
                f"Grid may be regional or use non-standard spacing."
            )
        gridtype = "regular"
    return gridtype


def to3d(array: ArrayLike) -> np.ndarray:
    """
    Reshape array to 3D format for spherical harmonic analysis.

    Converts input arrays to the (nlat, nlon, nfields) format required
    by spherical harmonic routines, where any additional dimensions beyond
    latitude and longitude are flattened into a single 'fields' dimension.

    Parameters
    ----------
    array : array_like
        Input array with shape (nlat, nlon, ...) where latitude and longitude
        are the first two dimensions

    Returns
    -------
    ndarray
        Reshaped array with shape (nlat, nlon, nfields) where nfields is the
        product of all dimensions beyond the first two

    Examples
    --------
    >>> # 4D array: (nlat=73, nlon=144, ntime=12, nlevel=17)
    >>> data_4d = np.random.randn(73, 144, 12, 17)
    >>> data_3d = to3d(data_4d)
    >>> print(data_3d.shape)  # (73, 144, 204)
    >>>
    >>> # Already 3D array remains unchanged
    >>> data_3d_input = np.random.randn(73, 144, 10)
    >>> result = to3d(data_3d_input)
    >>> print(result.shape)  # (73, 144, 10)
    """
    array = np.asarray(array)
    if array.ndim < 2:
        raise ValueError(
            f"Array must have at least 2 dimensions (lat, lon), got {array.ndim}D"
        )

    # Calculate new shape: (nlat, nlon, product_of_remaining_dims)
    new_shape = array.shape[:2] + (np.prod(array.shape[2:], dtype=int),)
    return array.reshape(new_shape)
