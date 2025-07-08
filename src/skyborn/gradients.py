import numpy as np
from typing import Union, Optional

EARTH_RADIUS = 6371e3  # Earth's mean radius (m)


def calculate_gradient(
    field: np.ndarray,
    coordinates: np.ndarray,
    axis: int = -1,
    radius: float = 6371000.0,
) -> np.ndarray:
    """Calculate gradient of an arbitrary dimensional array along specified coordinates

    Args:
        field: Data field for gradient calculation, can be any dimensional array, e.g., (time, level, lat, lon)
        coordinates: Coordinate array along which to calculate the gradient, such as latitude or longitude values
        axis: Specifies the dimensional axis for gradient calculation, defaults to the last dimension (-1)
        radius: Earth radius, default is 6371000.0 meters, used for latitude gradient calculation

    Returns:
        Gradient field with the same shape as the input field
    """
    # Check if input data dimensions match
    if coordinates.size != field.shape[axis]:
        raise ValueError(
            f"Coordinate array size ({coordinates.size}) does not match field size ({field.shape[axis]}) on the specified axis"
        )

    # Determine if coordinates are latitude or longitude
    is_latitude = False
    if np.min(coordinates) >= -90 and np.max(coordinates) <= 90:
        is_latitude = True
        # For latitude, calculate actual distance (in meters)
        if is_latitude:
            # Convert latitude to actual distance
            distances = coordinates * np.pi / 180.0 * radius
        else:
            # For longitude, we need to consider the latitude effect, but this requires additional latitude information
            # Here we simply use longitude differences directly
            distances = coordinates

    # Create output array with the same shape as input
    gradient = np.zeros_like(field, dtype=float)

    # To use numpy's advanced indexing, we need to create index arrays
    ndim = field.ndim
    idx_ranges = [slice(None)] * ndim

    # Use central difference for interior points
    inner_range = slice(1, field.shape[axis] - 1)
    idx_forward = idx_ranges.copy()
    idx_forward[axis] = slice(2, field.shape[axis])

    idx_center = idx_ranges.copy()
    idx_center[axis] = inner_range

    idx_backward = idx_ranges.copy()
    idx_backward[axis] = slice(0, field.shape[axis] - 2)

    # Use vectorized operations to calculate gradient for interior points
    forward_dists = np.diff(distances[1:])
    backward_dists = np.diff(distances[:-1])
    total_dists = distances[2:] - distances[:-2]

    # Create coefficient arrays with shape suitable for broadcasting
    shape = [1] * ndim
    shape[axis] = len(forward_dists)

    a0 = forward_dists.reshape(shape)
    b0 = backward_dists.reshape(shape)
    c0 = total_dists.reshape(shape)

    # Calculate gradient using weighted difference formula
    gradient[tuple(idx_center)] = (
        b0 / a0 / c0 * field[tuple(idx_forward)]
        - a0 / b0 / c0 * field[tuple(idx_backward)]
        + (a0 - b0) / a0 / b0 * field[tuple(idx_center)]
    )

    # Handle boundary points (forward and backward differences)
    # Left boundary
    left_idx = idx_ranges.copy()
    left_idx[axis] = 0
    left_idx_plus = idx_ranges.copy()
    left_idx_plus[axis] = 1
    gradient[tuple(left_idx)] = (
        field[tuple(left_idx_plus)] - field[tuple(left_idx)]
    ) / (distances[1] - distances[0])

    # Right boundary
    right_idx = idx_ranges.copy()
    right_idx[axis] = -1
    right_idx_minus = idx_ranges.copy()
    right_idx_minus[axis] = -2
    gradient[tuple(right_idx)] = (
        field[tuple(right_idx)] - field[tuple(right_idx_minus)]
    ) / (distances[-1] - distances[-2])

    return gradient


def calculate_meridional_gradient(
    field: np.ndarray,
    latitudes: np.ndarray,
    lat_axis: int = -1,
    radius: float = 6371000.0,
) -> np.ndarray:
    """Calculate meridional gradient (gradient along latitude direction)

    Args:
        field: Data field for gradient calculation, can be any dimensional array
        latitudes: Latitude array (degrees)
        lat_axis: Specifies the axis for latitude, defaults to the last dimension (-1)
        radius: Earth radius, default is 6371000.0 meters

    Returns:
        Meridional gradient field
    """
    return calculate_gradient(field, latitudes, axis=lat_axis, radius=radius)


def calculate_vertical_gradient(
    field: np.ndarray, pressure: np.ndarray, pressure_axis: int = -3
) -> np.ndarray:
    """Calculate vertical gradient (gradient along pressure direction)

    Args:
        field: Data field for gradient calculation
        pressure: Pressure array (Pa), must be monotonically decreasing
        pressure_axis: Specifies the axis for pressure, defaults to the third-to-last dimension (-3)

    Returns:
        Vertical gradient field
    """
    return calculate_gradient(field, pressure, axis=pressure_axis, radius=None)


def calculate_zonal_gradient(
    field: np.ndarray,
    longitudes: np.ndarray,
    latitudes: np.ndarray,
    lon_axis: int = -1,
    lat_axis: int = -2,
    radius: float = 6371000.0,
) -> np.ndarray:
    """Calculate zonal gradient (gradient along longitude direction)

    Args:
        field: Data field for gradient calculation, can be any dimensional array
        longitudes: Longitude array (degrees)
        latitudes: Latitude array (degrees), used to calculate actual distance between longitudes at different latitudes
        lon_axis: Specifies the axis for longitude, defaults to the last dimension (-1)
        lat_axis: Specifies the axis for latitude, defaults to the second-to-last dimension (-2)
        radius: Earth radius, default is 6371000.0 meters

    Returns:
        Zonal gradient field
    """
    # Get latitude factor to adjust actual distance between longitudes at different latitudes
    cos_lat = np.cos(np.radians(latitudes))

    # If field is 4D (time, level, lat, lon)
    if field.ndim == 4 and lon_axis == -1 and lat_axis == -2:
        # Create a latitude factor array with shape suitable for broadcasting
        cos_lat_expanded = cos_lat.reshape(1, 1, -1, 1)

        # Convert longitudes to actual distances considering latitude
        effective_distances = np.radians(longitudes) * radius * cos_lat_expanded

        # Calculate gradient
        return calculate_gradient(field, effective_distances, axis=lon_axis, radius=1.0)

    # If field is 3D (time, lat, lon)
    elif field.ndim == 3 and lon_axis == -1 and lat_axis == -2:
        cos_lat_expanded = cos_lat.reshape(1, -1, 1)
        effective_distances = np.radians(longitudes) * radius * cos_lat_expanded
        return calculate_gradient(field, effective_distances, axis=lon_axis, radius=1.0)

    else:
        # For other dimension combinations, create appropriate broadcasting shape
        broadcast_shape = [1] * field.ndim
        broadcast_shape[lat_axis] = len(latitudes)
        cos_lat_expanded = cos_lat.reshape(broadcast_shape)

        # Create effective distance array
        effective_longitudes = np.radians(longitudes) * radius

        # Calculate gradient for each latitude
        result = np.zeros_like(field)

        # Loop through each latitude (implementation depends on specific data structure, may need adjustment)
        for i in range(len(latitudes)):
            idx = [slice(None)] * field.ndim
            idx[lat_axis] = i

            # Adjust longitude distance for current latitude
            current_effective_dist = effective_longitudes * cos_lat[i]

            # Calculate gradient for current latitude
            result[tuple(idx)] = calculate_gradient(
                field[tuple(idx)], current_effective_dist, axis=lon_axis, radius=1.0
            )

        return result
