"""
This scripts contains functions that performs nearest, bilinear, and conservative interpolation
on xarray.Datasets. The original version of this script is available at WeatherBench2.

Qianye Su
suqianye2000@gmail.com

Reference
 - WeatherBench2 regridding:
     https://github.com/google-research/weatherbench2/blob/main/weatherbench2/regridding.py
"""

from __future__ import annotations

import xarray
import numpy as np
from sklearn import neighbors

import functools
import dataclasses
from typing import Union, Tuple, Optional

Array = Union[np.ndarray]


def _detect_coordinate_names(dataset: xarray.Dataset) -> Tuple[str, str]:
    """
    Detect latitude and longitude coordinate names in the dataset.

    Args:
        dataset: xarray Dataset

    Returns:
        Tuple of (longitude_name, latitude_name)

    Raises:
        ValueError: If coordinate names cannot be detected
    """
    # Common variations of coordinate names
    lon_names = ["longitude", "lon", "long", "x"]
    lat_names = ["latitude", "lat", "y"]

    # Find longitude coordinate
    lon_coord = None
    for name in lon_names:
        if name in dataset.dims:
            lon_coord = name
            break

    # Find latitude coordinate
    lat_coord = None
    for name in lat_names:
        if name in dataset.dims:
            lat_coord = name
            break

    if lon_coord is None or lat_coord is None:
        available_dims = list(dataset.dims.keys())
        raise ValueError(
            f"Could not detect longitude/latitude coordinates. "
            f"Available dimensions: {available_dims}. "
            f"Expected one of: lon={lon_names}, lat={lat_names}"
        )

    return lon_coord, lat_coord


@dataclasses.dataclass(frozen=True)
class Grid:
    """Representation of a rectilinear grid."""

    lon: np.ndarray
    lat: np.ndarray

    @classmethod
    def from_degrees(cls, lon: np.ndarray, lat: np.ndarray) -> Grid:
        return cls(np.deg2rad(lon), np.deg2rad(lat))

    @classmethod
    def from_dataset(cls, dataset: xarray.Dataset) -> Grid:
        """Create a Grid from an xarray Dataset by auto-detecting coordinates."""
        lon_name, lat_name = _detect_coordinate_names(dataset)
        lon_values = dataset[lon_name].values
        lat_values = dataset[lat_name].values
        return cls.from_degrees(lon_values, lat_values)

    @property
    def shape(self) -> tuple[int, int]:
        return (len(self.lon), len(self.lat))

    def _to_tuple(self) -> tuple[tuple[float, ...], tuple[float, ...]]:
        return tuple(self.lon.tolist()), tuple(self.lat.tolist())

    def __eq__(self, other):  # needed for hashability
        return isinstance(other, Grid) and self._to_tuple() == other._to_tuple()

    def __hash__(self):
        return hash(self._to_tuple())


@dataclasses.dataclass(frozen=True)
class Regridder:
    """Base class for regridding."""

    source: Grid
    target: Grid

    def regrid_array(self, field: Array) -> np.ndarray:
        """Regrid an array with dimensions (..., lon, lat) from source to target."""
        raise NotImplementedError

    def regrid_dataset(
        self,
        dataset: xarray.Dataset,
        lon_dim: Optional[str] = None,
        lat_dim: Optional[str] = None,
    ) -> xarray.Dataset:
        """
        Regrid an xarray.Dataset from source to target.

        Args:
            dataset: Input xarray Dataset
            lon_dim: Name of longitude dimension (auto-detected if None)
            lat_dim: Name of latitude dimension (auto-detected if None)

        Returns:
            Regridded xarray Dataset with preserved dimension order
        """
        # Auto-detect coordinate names if not provided
        if lon_dim is None or lat_dim is None:
            detected_lon, detected_lat = _detect_coordinate_names(dataset)
            lon_dim = lon_dim or detected_lon
            lat_dim = lat_dim or detected_lat

        # Store original dimension order for each variable
        original_dims = {}
        for var_name in dataset.data_vars:
            original_dims[var_name] = list(dataset[var_name].dims)

        # Ensure latitude is in ascending order
        if not (dataset[lat_dim].diff(lat_dim) > 0).all():
            dataset = dataset.isel({lat_dim: slice(None, None, -1)})  # Reverse
        assert (dataset[lat_dim].diff(lat_dim) > 0).all()

        # Create target grid coordinates
        target_lon_deg = np.rad2deg(self.target.lon)
        target_lat_deg = np.rad2deg(self.target.lat)

        # Process each variable separately to maintain dimension order
        regridded_vars = {}
        for var_name, var in dataset.data_vars.items():
            if lon_dim in var.dims and lat_dim in var.dims:
                # Apply regridding with proper dimension handling
                regridded_var = xarray.apply_ufunc(
                    self.regrid_array,
                    var,
                    input_core_dims=[[lon_dim, lat_dim]],
                    output_core_dims=[[lon_dim, lat_dim]],
                    exclude_dims={lon_dim, lat_dim},
                    vectorize=True,
                    dask="allowed",
                    output_dtypes=[var.dtype],
                    keep_attrs=True,
                )

                # Update coordinates while preserving dimension order
                regridded_var = regridded_var.assign_coords(
                    {lon_dim: target_lon_deg, lat_dim: target_lat_deg}
                )

                # Ensure original dimension order is maintained
                current_dims = list(regridded_var.dims)
                target_dims = original_dims[var_name].copy()

                # Update spatial dimensions in target_dims
                for i, dim in enumerate(target_dims):
                    if dim == lon_dim:
                        target_dims[i] = lon_dim
                    elif dim == lat_dim:
                        target_dims[i] = lat_dim

                # Transpose to match original order if needed
                if current_dims != target_dims:
                    regridded_var = regridded_var.transpose(*target_dims)

                regridded_vars[var_name] = regridded_var
            else:
                # Variables without spatial dimensions remain unchanged
                regridded_vars[var_name] = var

        # Create new dataset with regridded variables
        regridded_dataset = xarray.Dataset(
            regridded_vars,
            coords={
                **{
                    k: v
                    for k, v in dataset.coords.items()
                    if k not in [lon_dim, lat_dim]
                },
                lon_dim: target_lon_deg,
                lat_dim: target_lat_deg,
            },
            attrs=dataset.attrs,
        )

        return regridded_dataset


def nearest_neighbor_indices(source_grid: Grid, target_grid: Grid) -> np.ndarray:
    """Returns Haversine nearest neighbor indices from source_grid to target_grid."""
    # Construct a BallTree to find nearest neighbors on the sphere
    source_mesh = np.meshgrid(source_grid.lat, source_grid.lon, indexing="ij")
    target_mesh = np.meshgrid(target_grid.lat, target_grid.lon, indexing="ij")
    index_coords = np.stack([x.ravel() for x in source_mesh], axis=-1)
    query_coords = np.stack([x.ravel() for x in target_mesh], axis=-1)
    tree = neighbors.BallTree(index_coords, metric="haversine")
    indices = tree.query(query_coords, return_distance=False).squeeze(axis=-1)
    return indices


class NearestRegridder(Regridder):
    """Regrid with nearest neighbor interpolation."""

    def __init__(self, source: Grid, target: Grid):
        super().__init__(source, target)
        self._indices = None

    @property
    def indices(self):
        """The interpolation indices associated with source_grid."""
        if self._indices is None:
            self._indices = nearest_neighbor_indices(self.source, self.target)
        return self._indices

    def _nearest_neighbor_2d(self, array: Array) -> np.ndarray:
        """2D nearest neighbor interpolation using BallTree with optimized indexing."""
        if array.shape != self.source.shape:
            raise ValueError(
                f"Expected array.shape={array.shape} to match source.shape={self.source.shape}"
            )
        # Use advanced indexing for better performance
        array_flat = array.ravel()
        interpolated = array_flat[self.indices]
        return interpolated.reshape(self.target.shape)

    def regrid_array(self, field: Array) -> np.ndarray:
        # Use direct vectorization for better performance
        interp = np.vectorize(self._nearest_neighbor_2d, signature="(a,b)->(c,d)")
        return interp(field)


class BilinearRegridder(Regridder):
    """Regrid with bilinear interpolation."""

    def regrid_array(self, field: Array) -> np.ndarray:
        lat_source = self.source.lat
        lat_target = self.target.lat
        lon_source = self.source.lon
        lon_target = self.target.lon

        # Ensure the field has the correct shape (lon, lat)
        if field.shape != (len(lon_source), len(lat_source)):
            raise ValueError(
                f"Expected field shape {(len(lon_source), len(lat_source))}, "
                f"got {field.shape}"
            )

        # Interpolate over latitude first (for each longitude)
        lat_interp = np.zeros((len(lon_source), len(lat_target)))
        for i, lon_slice in enumerate(field):
            lat_interp[i, :] = np.interp(lat_target, lat_source, lon_slice)

        # Interpolate over longitude (for each target latitude)
        result = np.zeros((len(lon_target), len(lat_target)))
        for j in range(len(lat_target)):
            result[:, j] = np.interp(lon_target, lon_source, lat_interp[:, j])

        return result


def _assert_increasing(x: np.ndarray) -> None:
    if not (np.diff(x) > 0).all():
        raise ValueError(f"Array is not increasing: {x}")


def _latitude_cell_bounds(x: Array) -> np.ndarray:
    pi_over_2 = np.array([np.pi / 2], dtype=x.dtype)
    return np.concatenate((-pi_over_2, (x[:-1] + x[1:]) / 2, pi_over_2))


def _latitude_overlap(
    source_points: Array,
    target_points: Array,
) -> np.ndarray:
    """Calculate the area overlap as a function of latitude."""
    source_bounds = _latitude_cell_bounds(source_points)
    target_bounds = _latitude_cell_bounds(target_points)
    upper = np.minimum(target_bounds[1:, np.newaxis], source_bounds[np.newaxis, 1:])
    lower = np.maximum(target_bounds[:-1, np.newaxis], source_bounds[np.newaxis, :-1])
    # Normalized cell area: integral from lower to upper of cos(latitude)
    overlap = (upper > lower) * (np.sin(upper) - np.sin(lower))
    return overlap


def _conservative_latitude_weights(
    source_points: Array, target_points: Array
) -> np.ndarray:
    """Create a weight matrix for conservative regridding along latitude.

    Args:
        source_points: 1D latitude coordinates in radians for centers of source cells.
        target_points: 1D latitude coordinates in radians for centers of target cells.

    Returns:
        NumPy array with shape (target_size, source_size). Rows sum to 1.
    """
    _assert_increasing(source_points)
    _assert_increasing(target_points)
    weights = _latitude_overlap(source_points, target_points)

    # Handle zero-sum rows to avoid division by zero
    row_sums = np.sum(weights, axis=1, keepdims=True)

    # Avoid in-place division which causes broadcasting issues
    result = np.copy(weights)
    for i in range(result.shape[0]):
        if row_sums[i, 0] > 1e-15:
            result[i, :] /= row_sums[i, 0]
        else:
            # For zero-sum rows, distribute weight equally
            result[i, :] = 1.0 / result.shape[1]

    return result


def _align_phase_with(x, target, period):
    """Align the phase of a periodic number to match another."""
    shift_down = x > target + period / 2
    shift_up = x < target - period / 2
    return x + period * shift_up - period * shift_down


def _periodic_upper_bounds(x, period):
    x_plus = _align_phase_with(np.roll(x, -1), x, period)
    return (x + x_plus) / 2


def _periodic_lower_bounds(x, period):
    x_minus = _align_phase_with(np.roll(x, +1), x, period)
    return (x_minus + x) / 2


def _periodic_overlap(x0, x1, y0, y1, period):
    """Calculate the overlap between two intervals considering periodicity."""
    y0 = _align_phase_with(y0, x0, period)
    y1 = _align_phase_with(y1, x0, period)
    upper = np.minimum(x1, y1)
    lower = np.maximum(x0, y0)
    return np.maximum(upper - lower, 0)


def _longitude_overlap(
    first_points: Array,
    second_points: Array,
    period: float = 2 * np.pi,
) -> np.ndarray:
    """Calculate the area overlap as a function of longitude."""
    first_points = first_points % period
    first_upper = _periodic_upper_bounds(first_points, period)
    first_lower = _periodic_lower_bounds(first_points, period)

    second_points = second_points % period
    second_upper = _periodic_upper_bounds(second_points, period)
    second_lower = _periodic_lower_bounds(second_points, period)

    x0 = first_lower[:, np.newaxis]
    x1 = first_upper[:, np.newaxis]
    y0 = second_lower[np.newaxis, :]
    y1 = second_upper[np.newaxis, :]

    overlap_func = np.vectorize(_periodic_overlap, excluded=["period"])
    overlap = overlap_func(x0, x1, y0, y1, period=period)
    return overlap


def _conservative_longitude_weights(
    source_points: np.ndarray, target_points: np.ndarray
) -> np.ndarray:
    """Create a weight matrix for conservative regridding along longitude.

    Args:
        source_points: 1D longitude coordinates in radians for centers of source cells.
        target_points: 1D longitude coordinates in radians for centers of target cells.

    Returns:
        NumPy array with shape (target_size, source_size). Rows sum to 1.
    """
    _assert_increasing(source_points)
    _assert_increasing(target_points)
    weights = _longitude_overlap(target_points, source_points)

    # Handle zero-sum rows to avoid division by zero
    row_sums = np.sum(weights, axis=1, keepdims=True)
    nonzero_mask = row_sums > 1e-15

    # Avoid in-place division which causes broadcasting issues
    result = np.copy(weights)
    for i in range(result.shape[0]):
        if nonzero_mask[i, 0]:
            result[i, :] /= row_sums[i, 0]
        else:
            # For zero-sum rows, distribute weight equally
            result[i, :] = 1.0 / result.shape[1]

    return result


class ConservativeRegridder(Regridder):
    """Regrid with linear conservative regridding."""

    def __init__(self, source: Grid, target: Grid):
        super().__init__(source, target)
        # Pre-compute weights for better performance
        self._lon_weights = None
        self._lat_weights = None

    @property
    def lon_weights(self):
        """Cached longitude weights for performance."""
        if self._lon_weights is None:
            self._lon_weights = _conservative_longitude_weights(
                self.source.lon, self.target.lon
            )
        return self._lon_weights

    @property
    def lat_weights(self):
        """Cached latitude weights for performance."""
        if self._lat_weights is None:
            self._lat_weights = _conservative_latitude_weights(
                self.source.lat, self.target.lat
            )
        return self._lat_weights

    def _mean(self, field: Array) -> np.ndarray:
        """Computes cell-averages of field on the target grid with optimized einsum."""
        # Use cached weights for better performance
        result = np.einsum(
            "ac,bd,...cd->...ab",
            self.lon_weights,
            self.lat_weights,
            field,
            optimize=True,
        )
        return result

    def _nanmean(self, field: Array) -> np.ndarray:
        """Compute cell-averages skipping NaNs like np.nanmean."""
        nulls = np.isnan(field)
        total = self._mean(np.where(nulls, 0, field))
        count = self._mean(~nulls)
        with np.errstate(divide="ignore", invalid="ignore"):
            result = np.true_divide(total, count)
            result[count == 0] = np.nan  # Set divisions by zero to NaN
        return result

    def regrid_array(self, field: Array) -> np.ndarray:
        return self._nanmean(field)


# Convenience function for easy regridding
def regrid_dataset(
    dataset: xarray.Dataset,
    target_grid: Grid,
    method: str = "bilinear",
    lon_dim: Optional[str] = None,
    lat_dim: Optional[str] = None,
) -> xarray.Dataset:
    """
    Convenience function to regrid a dataset with optimized performance.

    Args:
        dataset: Input xarray Dataset
        target_grid: Target grid for regridding
        method: Interpolation method ('nearest', 'bilinear', 'conservative')
        lon_dim: Name of longitude dimension (auto-detected if None)
        lat_dim: Name of latitude dimension (auto-detected if None)

    Returns:
        Regridded xarray Dataset with preserved dimension order
    """
    # Create source grid from dataset
    source_grid = Grid.from_dataset(dataset)

    # Select regridder based on method with performance optimizations
    if method == "nearest":
        regridder = NearestRegridder(source_grid, target_grid)
    elif method == "bilinear":
        regridder = BilinearRegridder(source_grid, target_grid)
    elif method == "conservative":
        regridder = ConservativeRegridder(source_grid, target_grid)
    else:
        raise ValueError(
            f"Unknown method: {method}. Choose from 'nearest', 'bilinear', 'conservative'"
        )

    return regridder.regrid_dataset(dataset, lon_dim=lon_dim, lat_dim=lat_dim)
