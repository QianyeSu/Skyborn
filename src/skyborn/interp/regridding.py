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

import dataclasses
from typing import Optional, Tuple, Union

import numpy as np
import xarray
from sklearn import neighbors

from . import regrid as _native_regrid

_native_nearest_neighbor_indices = _native_regrid.nearest_neighbor_indices
_native_nearest_regrid_apply = _native_regrid.nearest_regrid_apply
_native_bilinear_regrid = _native_regrid.bilinear_regrid
_native_bilinear_regrid_nd = _native_regrid.bilinear_regrid_nd
_native_conservative_latitude_weights = _native_regrid.conservative_latitude_weights
_native_conservative_longitude_weights = _native_regrid.conservative_longitude_weights
_native_conservative_regrid = _native_regrid.conservative_regrid


# Keep BallTree as the default nearest-neighbor backend until the compiled helper
# matches sklearn's haversine tie-breaking on midpoint-heavy regular grids.
_ENABLE_NATIVE_NEAREST = False

__all__ = [
    "Grid",
    "Regridder",
    "NearestRegridder",
    "BilinearRegridder",
    "ConservativeRegridder",
    "nearest_neighbor_indices",
    "regrid_dataset",
]

Array = Union[np.ndarray]


def _is_strictly_increasing_1d(values: np.ndarray) -> bool:
    """Return True when a 1D coordinate array is strictly increasing."""
    return bool(np.all(np.diff(values) > 0))


def _supports_native_regular_grid(source: Grid, target: Grid) -> bool:
    """Check whether both grids satisfy the monotone 1D requirements of the native helpers."""
    return (
        len(source.lon) > 1
        and len(source.lat) > 1
        and len(target.lon) > 0
        and len(target.lat) > 0
        and _is_strictly_increasing_1d(source.lon)
        and _is_strictly_increasing_1d(source.lat)
        and _is_strictly_increasing_1d(target.lon)
        and _is_strictly_increasing_1d(target.lat)
    )


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
        available_dims = list(dataset.sizes.keys())
        raise ValueError(
            f"Could not detect longitude/latitude coordinates. "
            f"Available dimensions: {available_dims}. "
            f"Expected one of: lon={lon_names}, lat={lat_names}"
        )

    return lon_coord, lat_coord


def _grid_allclose(left: Grid, right: Grid) -> bool:
    """Return True when two grids describe the same 1D lon/lat coordinates."""
    return bool(
        left.lon.shape == right.lon.shape
        and left.lat.shape == right.lat.shape
        and np.allclose(left.lon, right.lon)
        and np.allclose(left.lat, right.lat)
    )


def _has_partial_spatial_dims(
    variable: xarray.DataArray, lon_dim: str, lat_dim: str
) -> bool:
    """Return True when a variable depends on only one horizontal dimension."""
    has_lon = lon_dim in variable.dims
    has_lat = lat_dim in variable.dims
    return has_lon != has_lat


def _validate_no_partial_spatial_variables(
    dataset: xarray.Dataset, lon_dim: str, lat_dim: str
) -> None:
    """Reject variables that would be misaligned after replacing the target grid."""
    for name, variable in dataset.data_vars.items():
        if _has_partial_spatial_dims(variable, lon_dim, lat_dim):
            raise ValueError(
                f"Variable {name!r} has only one horizontal dimension. "
                f"Variables must contain both {lon_dim!r} and {lat_dim!r}, or neither."
            )

    for name, coordinate in dataset.coords.items():
        if name in {lon_dim, lat_dim}:
            continue
        if _has_partial_spatial_dims(coordinate, lon_dim, lat_dim):
            raise ValueError(
                f"Coordinate {name!r} has only one horizontal dimension. "
                f"Coordinates must contain both {lon_dim!r} and {lat_dim!r}, or neither."
            )


def _select_regridder(
    source_grid: Grid,
    target_grid: Grid,
    method: str,
) -> Regridder:
    """Construct the regridder for a supported method."""
    if method == "nearest":
        return NearestRegridder(source_grid, target_grid)
    if method == "bilinear":
        return BilinearRegridder(source_grid, target_grid)
    if method == "conservative":
        return ConservativeRegridder(source_grid, target_grid)
    raise ValueError(
        f"Unknown method: {method}. Choose from 'nearest', 'bilinear', 'conservative'"
    )


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
        lat_diff = dataset[lat_dim].diff(lat_dim)
        if not (lat_diff > 0).all():
            if not (lat_diff < 0).all():
                raise ValueError(
                    f"Latitude coordinate {lat_dim!r} must be strictly monotonic"
                )
            dataset = dataset.isel({lat_dim: slice(None, None, -1)})  # Reverse
        assert (dataset[lat_dim].diff(lat_dim) > 0).all()
        _validate_no_partial_spatial_variables(dataset, lon_dim, lat_dim)

        dataset_source_grid = Grid.from_degrees(
            dataset[lon_dim].values,
            dataset[lat_dim].values,
        )
        active_regridder = self
        if not _grid_allclose(self.source, dataset_source_grid):
            reversed_source_grid = Grid(self.source.lon, self.source.lat[::-1])
            if not _grid_allclose(reversed_source_grid, dataset_source_grid):
                raise ValueError(
                    "Regridder source grid does not match dataset coordinates"
                )
            active_regridder = self.__class__(dataset_source_grid, self.target)

        # Create target grid coordinates
        target_lon_deg = np.rad2deg(self.target.lon)
        target_lat_deg = np.rad2deg(self.target.lat)

        # Process each variable separately to maintain dimension order
        regridded_vars = {}
        for var_name, var in dataset.data_vars.items():
            if lon_dim in var.dims and lat_dim in var.dims:
                # Apply regridding with proper dimension handling
                regridded_var = xarray.apply_ufunc(
                    active_regridder.regrid_array,
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
    if _ENABLE_NATIVE_NEAREST and _native_nearest_neighbor_indices is not None:
        return _native_nearest_neighbor_indices(
            source_grid.lon,
            source_grid.lat,
            target_grid.lon,
            target_grid.lat,
        )

    # Construct a BallTree to find nearest neighbors on the sphere
    source_mesh = np.meshgrid(source_grid.lon, source_grid.lat, indexing="ij")
    target_mesh = np.meshgrid(target_grid.lon, target_grid.lat, indexing="ij")
    index_coords = np.stack([source_mesh[1].ravel(), source_mesh[0].ravel()], axis=-1)
    query_coords = np.stack([target_mesh[1].ravel(), target_mesh[0].ravel()], axis=-1)
    tree = neighbors.BallTree(index_coords, metric="haversine")
    indices = tree.query(query_coords, return_distance=False).squeeze(axis=-1)
    return indices


def _gather_flat_spatial(
    field_arr: np.ndarray,
    indices: np.ndarray,
    source_shape: tuple[int, int],
    target_shape: tuple[int, int],
) -> np.ndarray:
    """Apply flat spatial indices across arbitrary leading dimensions."""
    src_size = source_shape[0] * source_shape[1]
    flat = field_arr.reshape(-1, src_size)
    gathered = np.take(flat, indices, axis=1)
    return gathered.reshape(field_arr.shape[:-2] + target_shape)


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

    def _nearest_neighbor_nd(self, field: Array) -> np.ndarray:
        """Apply the cached flat indices across leading dimensions."""
        field_arr = np.asarray(field)
        if field_arr.shape[-2:] != self.source.shape:
            raise ValueError(
                f"Expected field shape {self.source.shape}, got {field_arr.shape[-2:]}"
            )

        if _native_nearest_regrid_apply is not None:
            return _native_nearest_regrid_apply(
                field_arr,
                self.indices,
                self.target.shape[0],
                self.target.shape[1],
            )

        return _gather_flat_spatial(
            field_arr,
            self.indices,
            self.source.shape,
            self.target.shape,
        )

    def regrid_array(self, field: Array) -> np.ndarray:
        return self._nearest_neighbor_nd(field)


class BilinearRegridder(Regridder):
    """Regrid with bilinear interpolation."""

    def _bilinear_2d(self, field: Array) -> np.ndarray:
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

        if _native_bilinear_regrid is not None and _supports_native_regular_grid(
            self.source, self.target
        ):
            return _native_bilinear_regrid(
                np.asarray(field, dtype=np.float64),
                lon_source,
                lat_source,
                lon_target,
                lat_target,
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

    def _bilinear_nd(self, field: Array) -> np.ndarray:
        """Apply bilinear interpolation across arbitrary leading dimensions."""
        field_arr = np.asarray(field)
        if field_arr.shape[-2:] != self.source.shape:
            raise ValueError(
                f"Expected field shape {self.source.shape}, got {field_arr.shape[-2:]}"
            )

        if field_arr.ndim == 2:
            return self._bilinear_2d(field_arr)

        if _native_bilinear_regrid_nd is not None and _supports_native_regular_grid(
            self.source, self.target
        ):
            return _native_bilinear_regrid_nd(
                np.asarray(field_arr, dtype=np.float64),
                self.source.lon,
                self.source.lat,
                self.target.lon,
                self.target.lat,
            )

        leading_shape = field_arr.shape[:-2]
        result = np.empty(leading_shape + self.target.shape, dtype=np.float64)
        for index in np.ndindex(leading_shape):
            result[index] = self._bilinear_2d(field_arr[index])
        return result

    def regrid_array(self, field: Array) -> np.ndarray:
        return self._bilinear_nd(field)


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
    if _native_conservative_latitude_weights is not None:
        return _native_conservative_latitude_weights(source_points, target_points)

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
    if _native_conservative_longitude_weights is not None:
        return _native_conservative_longitude_weights(source_points, target_points)

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

    def _python_nanmean(self, field: Array) -> np.ndarray:
        """Compute cell-averages skipping NaNs using the Python fallback path."""
        nulls = np.isnan(field)
        total = self._mean(np.where(nulls, 0, field))
        count = self._mean(~nulls)
        with np.errstate(divide="ignore", invalid="ignore"):
            result = np.true_divide(total, count)
            result[count == 0] = np.nan  # Set divisions by zero to NaN
        return result

    def _conservative_2d(self, field: Array) -> np.ndarray:
        """Apply conservative regridding to a single 2D lon/lat field."""
        field_arr = np.asarray(field)
        if field_arr.shape != self.source.shape:
            raise ValueError(
                f"Expected field shape {self.source.shape}, got {field_arr.shape}"
            )

        if _native_conservative_regrid is not None:
            return _native_conservative_regrid(
                np.asarray(field_arr, dtype=np.float64),
                self.lon_weights,
                self.lat_weights,
            )

        return self._python_nanmean(field_arr)

    def _conservative_nd(self, field: Array) -> np.ndarray:
        """Apply conservative regridding across arbitrary leading dimensions."""
        field_arr = np.asarray(field)
        if field_arr.shape[-2:] != self.source.shape:
            raise ValueError(
                f"Expected field shape {self.source.shape}, got {field_arr.shape[-2:]}"
            )

        if field_arr.ndim == 2:
            return self._conservative_2d(field_arr)

        if _native_conservative_regrid is not None:
            return _native_conservative_regrid(
                np.asarray(field_arr, dtype=np.float64),
                self.lon_weights,
                self.lat_weights,
            )

        return self._python_nanmean(field_arr)

    def regrid_array(self, field: Array) -> np.ndarray:
        return self._conservative_nd(field)


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

    return _select_regridder(source_grid, target_grid, method).regrid_dataset(
        dataset,
        lon_dim=lon_dim,
        lat_dim=lat_dim,
    )
