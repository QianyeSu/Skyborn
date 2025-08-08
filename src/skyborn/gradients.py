import numpy as np
from typing import Union, Optional

__all__ = [
    "calculate_gradient",
    "calculate_meridional_gradient",
    "calculate_zonal_gradient",
    "calculate_vertical_gradient",
    "EARTH_RADIUS",
]

EARTH_RADIUS = 6371e3  # Earth's mean radius (m)


def calculate_gradient(
    field: np.ndarray,
    coordinates: np.ndarray,
    axis: int = -1,
    radius: Optional[float] = None,
) -> np.ndarray:
    """Calculate gradient of an arbitrary dimensional array along specified coordinates

    Args:
        field: Data field for gradient calculation, can be any dimensional array, e.g., (time, level, lat, lon)
        coordinates: Coordinate array along which to calculate the gradient, such as latitude or longitude values
        axis: Specifies the dimensional axis for gradient calculation, defaults to the last dimension (-1)
        radius: Earth radius, default is None. If provided, coordinates are treated as latitude in degrees
                and converted to distance in meters. If None, coordinates are used directly.

    Returns:
        Gradient field with the same shape as the input field
    """
    # Check if input data dimensions match
    if coordinates.size != field.shape[axis]:
        raise ValueError(
            f"Coordinate array size ({coordinates.size}) does not match field size ({field.shape[axis]}) on the specified axis"
        )

    # If radius is provided, treat coordinates as latitude and convert to distance
    if radius is not None:
        # Convert latitude to actual distance
        distances = coordinates * np.pi / 180.0 * radius
    else:
        # Use coordinates directly
        distances = coordinates

    # Create output array with the same shape as input
    gradient = np.zeros_like(field, dtype=float)

    # Handle edge case: single point (no gradient can be computed)
    if field.shape[axis] < 2:
        return gradient  # Return zero gradient

    # To use numpy's advanced indexing, we need to create index arrays
    ndim = field.ndim
    idx_ranges = [slice(None)] * ndim

    # Handle the simple case where distances are uniform
    dx = np.diff(distances)
    if len(dx) > 0 and np.allclose(dx, dx[0]):
        # Uniform spacing - use simple central differences
        h = dx[0]

        # Central differences for interior points
        for i in range(1, field.shape[axis] - 1):
            idx = idx_ranges.copy()
            idx[axis] = i
            idx_forward = idx_ranges.copy()
            idx_forward[axis] = i + 1
            idx_backward = idx_ranges.copy()
            idx_backward[axis] = i - 1

            gradient[tuple(idx)] = (
                field[tuple(idx_forward)] - field[tuple(idx_backward)]
            ) / (2.0 * h)

        # Forward difference for left boundary
        idx_left = idx_ranges.copy()
        idx_left[axis] = 0
        idx_left_plus = idx_ranges.copy()
        idx_left_plus[axis] = 1
        gradient[tuple(idx_left)] = (
            field[tuple(idx_left_plus)] - field[tuple(idx_left)]
        ) / h

        # Backward difference for right boundary
        idx_right = idx_ranges.copy()
        idx_right[axis] = -1
        idx_right_minus = idx_ranges.copy()
        idx_right_minus[axis] = -2
        gradient[tuple(idx_right)] = (
            field[tuple(idx_right)] - field[tuple(idx_right_minus)]
        ) / h

    else:
        # Non-uniform spacing - use variable spacing formulas
        for i in range(field.shape[axis]):
            idx = idx_ranges.copy()
            idx[axis] = i

            if i == 0:
                # Forward difference for first point
                idx_plus = idx_ranges.copy()
                idx_plus[axis] = 1
                gradient[tuple(idx)] = (field[tuple(idx_plus)] - field[tuple(idx)]) / (
                    distances[1] - distances[0]
                )
            elif i == field.shape[axis] - 1:
                # Backward difference for last point
                idx_minus = idx_ranges.copy()
                idx_minus[axis] = i - 1
                gradient[tuple(idx)] = (field[tuple(idx)] - field[tuple(idx_minus)]) / (
                    distances[i] - distances[i - 1]
                )
            else:
                # Central difference for interior points
                idx_forward = idx_ranges.copy()
                idx_forward[axis] = i + 1
                idx_backward = idx_ranges.copy()
                idx_backward[axis] = i - 1

                h1 = distances[i] - distances[i - 1]
                h2 = distances[i + 1] - distances[i]

                # Use weighted central difference for non-uniform spacing
                gradient[tuple(idx)] = (
                    -h2 / (h1 * (h1 + h2)) * field[tuple(idx_backward)]
                    + (h2 - h1) / (h1 * h2) * field[tuple(idx)]
                    + h1 / (h2 * (h1 + h2)) * field[tuple(idx_forward)]
                )

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
