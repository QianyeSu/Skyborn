"""
Data management tools for windspharm.VectorWind and spharm.Spharmt.

This module provides utilities for preparing data for use with spherical harmonic
wind analysis, including dimension reordering, shape management, and coordinate
system handling.

Main Functions:
    prep_data: Prepare data for VectorWind input
    recover_data: Restore original data shape/order
    get_recovery: Create recovery function for multiple arrays
    reverse_latdim: Reverse latitude dimension order
    order_latdim: Ensure north-to-south latitude ordering
"""

from __future__ import annotations
from typing import Tuple, Dict, Any, List, Callable, Union
import numpy as np

# Type aliases
ArrayLike = Union[np.ndarray, np.ma.MaskedArray]
DimOrder = str
RecoveryInfo = Dict[str, Any]


def _order_dims(data: np.ndarray, dim_order: DimOrder) -> Tuple[np.ndarray, DimOrder]:
    """
    Reorder array dimensions to place lat/lon first.

    Internal function to reorder dimensions so that latitude ('y') and
    longitude ('x') are the first two dimensions.

    Parameters
    ----------
    data : ndarray
        Input data array
    dim_order : str
        String describing dimension order using 'x' for longitude, 'y' for latitude

    Returns
    -------
    reordered_data : ndarray
        Data with dimensions reordered
    new_order : str
        New dimension order string

    Raises
    ------
    ValueError
        If 'x' or 'y' not found in dim_order
    """
    if "x" not in dim_order or "y" not in dim_order:
        raise ValueError(
            "Dimension order must contain 'x' (longitude) and 'y' (latitude). "
            f"Got: '{dim_order}'"
        )

    # Find positions of longitude and latitude
    lon_pos = dim_order.lower().find("x")
    lat_pos = dim_order.lower().find("y")

    # Move longitude to front
    data = np.rollaxis(data, lon_pos)
    if lat_pos < lon_pos:
        lat_pos += 1

    # Move latitude to front (after longitude)
    data = np.rollaxis(data, lat_pos)

    # Create new dimension order string
    out_order = dim_order.replace("x", "").replace("y", "")
    out_order = "yx" + out_order

    return data, out_order


def _reshape_for_spharm(data: np.ndarray) -> Tuple[np.ndarray, Tuple[int, ...]]:
    """
    Reshape data for spherical harmonic analysis.

    Reshapes to (nlat, nlon, other) where 'other' combines all remaining dimensions.

    Parameters
    ----------
    data : ndarray
        Input data array with lat/lon as first two dimensions

    Returns
    -------
    reshaped : ndarray
        Reshaped array
    original_shape : tuple
        Original shape for recovery
    """
    new_shape = data.shape[:2] + (np.prod(data.shape[2:], dtype=int),)
    return data.reshape(new_shape), data.shape


def prep_data(data: ArrayLike, dim_order: DimOrder) -> Tuple[np.ndarray, RecoveryInfo]:
    """
    Prepare data for input to VectorWind or spherical harmonic transforms.

    This function reorders and reshapes input data to the format required by
    VectorWind and spherical harmonic analysis routines. It returns the prepared
    data along with information needed to recover the original format.

    Parameters
    ----------
    data : array_like
        Input data array with at least 2 dimensions. Must contain latitude
        and longitude dimensions.
    dim_order : str
        String specifying dimension order using 'x' for longitude and 'y' for
        latitude. Other characters represent additional dimensions.

    Returns
    -------
    prepared_data : ndarray
        Data reshaped to (latitude, longitude, other) format
    recovery_info : dict
        Dictionary containing information needed to recover original format.
        Contains keys: 'intermediate_shape', 'intermediate_order', 'original_order'

    Raises
    ------
    ValueError
        If dim_order doesn't contain both 'x' and 'y'

    See Also
    --------
    recover_data : Recover original data format
    get_recovery : Create recovery function for multiple arrays

    Examples
    --------
    >>> import numpy as np
    >>> # Prepare 4D data: (time=12, level=17, lat=73, lon=144)
    >>> data = np.random.randn(12, 17, 73, 144)
    >>> pdata, info = prep_data(data, 'tzyx')
    >>> print(pdata.shape)  # (73, 144, 204)
    >>>
    >>> # Prepare data with arbitrary dimension labels
    >>> data = np.random.randn(144, 16, 73, 21)
    >>> pdata, info = prep_data(data, 'xayb')
    """
    # Validate input
    data = np.asarray(data)
    if data.ndim < 2:
        raise ValueError(
            f"Data must have at least 2 dimensions, got {data.ndim}D array"
        )

    # Reorder dimensions to put lat/lon first
    prepared_data, intermediate_order = _order_dims(data, dim_order)

    # Reshape for spherical harmonic analysis
    prepared_data, intermediate_shape = _reshape_for_spharm(prepared_data)

    # Create recovery information
    recovery_info = {
        "intermediate_shape": intermediate_shape,
        "intermediate_order": intermediate_order,
        "original_order": dim_order,
    }

    return prepared_data, recovery_info


def recover_data(prepared_data: np.ndarray, recovery_info: RecoveryInfo) -> np.ndarray:
    """
    Recover original shape and dimension order of processed data.

    This function reverses the operations performed by prep_data, restoring
    the original shape and dimension ordering of data that has been processed
    by VectorWind or spherical harmonic methods.

    Parameters
    ----------
    prepared_data : ndarray
        Data array with 2 or 3 dimensions where first two are (lat, lon)
    recovery_info : dict
        Recovery information dictionary from prep_data

    Returns
    -------
    ndarray
        Data restored to original shape and dimension order

    See Also
    --------
    prep_data : Prepare data for processing
    get_recovery : Create recovery function for multiple arrays

    Examples
    --------
    >>> # Recover data processed with prep_data
    >>> original_data = recover_data(processed_data, info)
    >>>
    >>> # Use with VectorWind output
    >>> pdata, info = prep_data(wind_data, 'tzyx')
    >>> vw = VectorWind(u_prepared, v_prepared)
    >>> vorticity = vw.vorticity()
    >>> vort_original = recover_data(vorticity, info)
    """
    if not isinstance(recovery_info, dict):
        raise TypeError("recovery_info must be a dictionary")

    required_keys = {"intermediate_shape", "intermediate_order", "original_order"}
    missing_keys = required_keys - set(recovery_info.keys())
    if missing_keys:
        raise ValueError(f"Missing keys in recovery_info: {missing_keys}")

    # Restore intermediate shape (full dimensionality, spherical harmonic order)
    data = prepared_data.reshape(recovery_info["intermediate_shape"])

    # Determine dimension reordering
    intermediate_order = recovery_info["intermediate_order"]
    original_order = recovery_info["original_order"]

    # Calculate how to reorder dimensions back to original
    roll_dims = np.array(
        [intermediate_order.index(dim) for dim in original_order[::-1]]
    )

    # Apply dimension reordering
    for i in range(len(roll_dims)):
        # Roll axis to the front
        data = np.rollaxis(data, roll_dims[i])
        # Update remaining roll positions
        roll_dims = np.where(roll_dims < roll_dims[i], roll_dims + 1, roll_dims)

    return data


_RECOVERY_DOCSTRING_TEMPLATE = """Shape and dimension recovery function.

Recovers variable shape and dimension order according to:

{recovery_info}

Parameters
----------
*args : ndarray
    Variable number of arrays to recover

Returns
-------
list of ndarray
    Recovered arrays with original shape and dimension order
"""


def get_recovery(recovery_info: RecoveryInfo) -> Callable[..., List[np.ndarray]]:
    """
    Create a recovery function for multiple arrays.

    Returns a function that can recover the original shape and dimension
    order of multiple arrays using a single recovery information dictionary.
    This is useful when processing multiple related arrays.

    Parameters
    ----------
    recovery_info : dict
        Recovery information dictionary from prep_data

    Returns
    -------
    recover_func : callable
        Function that takes variable number of arrays and returns them
        with original shape and dimension order restored

    See Also
    --------
    prep_data : Prepare data for processing
    recover_data : Recover single array

    Examples
    --------
    >>> # Prepare multiple wind components
    >>> u_prep, info = prep_data(u, 'tzyx')
    >>> v_prep, _ = prep_data(v, 'tzyx')
    >>>
    >>> # Process with VectorWind
    >>> vw = VectorWind(u_prep, v_prep)
    >>> sf, vp = vw.sfvp()
    >>>
    >>> # Create recovery function and use it
    >>> recover = get_recovery(info)
    >>> u_orig, v_orig, sf_orig, vp_orig = recover(u_prep, v_prep, sf, vp)
    """

    def _recover_multiple(*arrays: np.ndarray) -> List[np.ndarray]:
        """Recover multiple arrays using the same recovery info."""
        return [recover_data(array, recovery_info) for array in arrays]

    # Create nice documentation for the returned function
    info_str = "\n".join(f"  '{key}': {value}" for key, value in recovery_info.items())

    # Set function metadata
    _recover_multiple.__name__ = "recover"
    _recover_multiple.__doc__ = _RECOVERY_DOCSTRING_TEMPLATE.format(
        recovery_info=info_str
    )

    return _recover_multiple


def reverse_latdim(
    u: ArrayLike, v: ArrayLike, axis: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reverse the latitude dimension of wind components.

    Creates copies of the input wind components with the latitude dimension
    reversed. This is useful for converting between different latitude
    ordering conventions (e.g., north-to-south vs south-to-north).

    Parameters
    ----------
    u, v : array_like
        Zonal and meridional wind components
    axis : int, default 0
        Index of the latitude dimension to reverse

    Returns
    -------
    u_reversed : ndarray
        Zonal wind component with latitude dimension reversed (copy)
    v_reversed : ndarray
        Meridional wind component with latitude dimension reversed (copy)

    See Also
    --------
    order_latdim : Ensure north-to-south latitude ordering

    Examples
    --------
    >>> # Reverse first dimension (default)
    >>> u_rev, v_rev = reverse_latdim(u, v)
    >>>
    >>> # Reverse third dimension
    >>> u_rev, v_rev = reverse_latdim(u, v, axis=2)
    """
    u = np.asarray(u)
    v = np.asarray(v)

    if u.shape != v.shape:
        raise ValueError(
            f"Wind components must have same shape. " f"Got u: {u.shape}, v: {v.shape}"
        )

    if axis >= u.ndim or axis < -u.ndim:
        raise ValueError(
            f"axis {axis} is out of bounds for array of dimension {u.ndim}"
        )

    # Create slice to reverse the specified axis
    slice_list = [slice(None)] * u.ndim
    slice_list[axis] = slice(None, None, -1)

    # Apply reversal and make copies
    u_reversed = u.copy()[tuple(slice_list)]
    v_reversed = v.copy()[tuple(slice_list)]

    return u_reversed, v_reversed


def order_latdim(
    lat_coords: ArrayLike, u: ArrayLike, v: ArrayLike, axis: int = 0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Ensure latitude dimension is ordered north-to-south.

    Checks the latitude coordinate array and ensures it runs from north to south.
    If necessary, reverses the latitude coordinates and the corresponding dimension
    in the wind components. Always returns copies of the input arrays.

    Parameters
    ----------
    lat_coords : array_like
        Array of latitude coordinate values
    u, v : array_like
        Zonal and meridional wind components
    axis : int, default 0
        Index of the latitude dimension in the wind components

    Returns
    -------
    lat_ordered : ndarray
        Latitude coordinates ordered north-to-south (copy)
    u_ordered : ndarray
        Zonal wind component with latitude ordered north-to-south (copy)
    v_ordered : ndarray
        Meridional wind component with latitude ordered north-to-south (copy)

    See Also
    --------
    reverse_latdim : Reverse latitude dimension

    Examples
    --------
    >>> # Order latitude when it's the first dimension
    >>> lat_ord, u_ord, v_ord = order_latdim(lat, u, v)
    >>>
    >>> # Order latitude when it's the third dimension
    >>> lat_ord, u_ord, v_ord = order_latdim(lat, u, v, axis=2)
    """
    lat_coords = np.asarray(lat_coords).copy()
    u = np.asarray(u)
    v = np.asarray(v)

    # Validate inputs
    if u.shape != v.shape:
        raise ValueError(
            f"Wind components must have same shape. " f"Got u: {u.shape}, v: {v.shape}"
        )

    if axis >= u.ndim or axis < -u.ndim:
        raise ValueError(
            f"axis {axis} is out of bounds for array of dimension {u.ndim}"
        )

    if len(lat_coords) != u.shape[axis]:
        raise ValueError(
            f"Length of lat_coords ({len(lat_coords)}) must match "
            f"size of axis {axis} ({u.shape[axis]}) in wind components"
        )

    # Check if latitude needs to be reversed (south-to-north -> north-to-south)
    if lat_coords[0] < lat_coords[-1]:
        # Reverse latitude coordinates
        lat_coords = lat_coords[::-1]
        # Reverse wind components along the latitude axis
        u, v = reverse_latdim(u, v, axis=axis)
    else:
        # Already north-to-south, just make copies
        u, v = u.copy(), v.copy()

    return lat_coords, u, v
