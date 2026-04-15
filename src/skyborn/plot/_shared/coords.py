"""Shared coordinate and field-normalization helpers for Skyborn plots."""

from __future__ import annotations

from typing import Any

import numpy as np
import xarray as xr

from .types import CoordinateShape, TargetDims


def _filled_float_array(values: Any) -> np.ndarray:
    """Return a plain float array with masked values replaced by NaN."""
    array = np.ma.asarray(values, dtype=float)
    if np.ma.isMaskedArray(array):
        return np.asarray(array.filled(np.nan), dtype=float)
    return np.asarray(array, dtype=float)


_filled_scalar_field_array = _filled_float_array


def _subset_ready_array(value: Any) -> np.ndarray:
    """Normalize one array so it can be subset safely by retained points."""
    array = np.ma.asarray(value)
    if np.ma.isMaskedArray(array):
        fill_value = np.nan if np.issubdtype(array.dtype, np.number) else None
        return np.asarray(array.filled(fill_value))
    return np.asarray(array)


def _maybe_squeezed_dataarray(value: Any) -> xr.DataArray | None:
    """Return a squeezed DataArray or ``None`` for non-xarray inputs."""
    if not isinstance(value, xr.DataArray):
        return None
    return value.squeeze(drop=True)


def _transpose_dataarray_if_possible(value: Any, target_dims: TargetDims) -> Any:
    """Align a DataArray to the target plotting dimension order when possible."""
    da = _maybe_squeezed_dataarray(value)
    if da is None or target_dims is None or da.ndim != len(target_dims):
        return da if da is not None else value
    if da.dims == target_dims or set(da.dims) != set(target_dims):
        return da
    return da.transpose(*target_dims)


def _normalized_field_array(value: Any, target_dims: TargetDims = None) -> np.ndarray:
    """Return one style or mask field as an ndarray in plotting dimension order."""
    value = _transpose_dataarray_if_possible(value, target_dims)
    if isinstance(value, xr.DataArray):
        value = value.data
    return _subset_ready_array(value)


def _infer_grid_from_reference(
    value: Any,
    grid_shape: CoordinateShape,
    target_dims: TargetDims,
) -> bool:
    """Infer whether 1D x/y should be expanded to a grid from companion fields."""
    if value is None:
        return False
    array = _normalized_field_array(value, target_dims=target_dims)
    return array.ndim == 2 and array.shape == grid_shape


def _normalize_coordinates(
    x: Any,
    y: Any,
    *,
    reference_fields: tuple[Any, ...],
) -> tuple[np.ndarray, np.ndarray, CoordinateShape, bool, TargetDims]:
    """Normalize paired points or grid coordinates into a plotting form."""
    x_obj = _maybe_squeezed_dataarray(x)
    y_obj = _maybe_squeezed_dataarray(y)

    target_dims: TargetDims = None
    if x_obj is not None and y_obj is not None and x_obj.ndim == y_obj.ndim == 2:
        y_obj = _transpose_dataarray_if_possible(y_obj, tuple(x_obj.dims))
        target_dims = tuple(x_obj.dims)
    elif x_obj is not None and y_obj is not None and x_obj.ndim == y_obj.ndim == 1:
        target_dims = (y_obj.dims[0], x_obj.dims[0])

    x_values = _filled_float_array(x_obj.data if x_obj is not None else x)
    y_values = _filled_float_array(y_obj.data if y_obj is not None else y)

    if x_values.ndim == 2 or y_values.ndim == 2:
        if x_values.ndim != 2 or y_values.ndim != 2:
            raise ValueError("x and y must both be 1D or both be 2D")
        if x_values.shape != y_values.shape:
            raise ValueError("2D x and y coordinates must have the same shape")
        return x_values, y_values, x_values.shape, True, target_dims

    if x_values.ndim != 1 or y_values.ndim != 1:
        raise ValueError("x and y must each be 1D or 2D arrays")

    grid_shape = (int(y_values.size), int(x_values.size))
    grid_like = x_values.size != y_values.size
    if not grid_like:
        grid_like = any(
            _infer_grid_from_reference(field, grid_shape, target_dims)
            for field in reference_fields
        )

    if grid_like:
        x_grid, y_grid = np.meshgrid(x_values, y_values, indexing="xy")
        return x_grid, y_grid, grid_shape, True, target_dims

    if x_values.shape != y_values.shape:
        raise ValueError(
            "Paired scatter points require x and y to have the same 1D length"
        )
    return x_values, y_values, x_values.shape, False, target_dims


def _normalize_selection_mask(
    *,
    where: Any,
    mask: Any,
    candidate_shape: CoordinateShape,
    target_dims: TargetDims,
) -> np.ndarray:
    """Return the boolean candidate-selection mask in plotting coordinate order."""
    if where is not None and mask is not None:
        raise ValueError("Use only one of 'where' or 'mask'")

    selector = where if where is not None else mask
    if selector is None:
        return np.ones(candidate_shape, dtype=bool)

    selector_array = _normalized_field_array(selector, target_dims=target_dims)
    if selector_array.shape != candidate_shape:
        raise ValueError(
            f"Selection mask shape {selector_array.shape} must match candidate shape "
            f"{candidate_shape}"
        )

    if np.issubdtype(selector_array.dtype, np.bool_):
        return np.asarray(selector_array, dtype=bool)

    finite = np.isfinite(selector_array)
    return finite & (selector_array != 0)


def _axis_coordinate_1d(values: Any, axis_name: str) -> np.ndarray | None:
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        return values
    if values.ndim != 2:
        return None
    if axis_name == "x":
        return values[0, :]
    if axis_name == "y":
        return values[:, 0]
    raise ValueError(f"Unsupported axis_name {axis_name!r}")


def _axis_is_uniform(values: Any, rtol: float = 1e-6, atol: float = 1e-10) -> bool:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size < 3:
        return True
    diffs = np.diff(values)
    if not np.all(np.isfinite(diffs)):
        return False
    reference = float(np.nanmedian(diffs))
    tolerance = max(float(atol), abs(reference) * float(rtol))
    return np.all(np.abs(diffs - reference) <= tolerance)


def _cell_edges_from_axis(values: Any) -> np.ndarray:
    """Return inferred cell edges for one center-coordinate axis."""
    axis = np.asarray(values, dtype=float)
    if axis.ndim != 1 or axis.size == 0:
        raise ValueError("Cell edges require a non-empty 1D axis")

    edges = np.empty(axis.size + 1, dtype=float)
    if axis.size == 1:
        center = float(axis[0])
        edges[0] = center - 0.5
        edges[1] = center + 0.5
        return edges

    midpoints = 0.5 * (axis[:-1] + axis[1:])
    edges[1:-1] = midpoints
    edges[0] = axis[0] - 0.5 * (axis[1] - axis[0])
    edges[-1] = axis[-1] + 0.5 * (axis[-1] - axis[-2])
    return edges


def _rectilinear_cell_edges(
    x: Any,
    y: Any,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Return cell edges for 1D or meshgrid-like rectilinear coordinates."""
    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    x_axis = _axis_coordinate_1d(x_values, "x")
    y_axis = _axis_coordinate_1d(y_values, "y")
    if x_axis is None or y_axis is None:
        return None

    if x_values.ndim == 2 and not np.allclose(
        x_values, x_axis[np.newaxis, :], equal_nan=True
    ):
        return None
    if y_values.ndim == 2 and not np.allclose(
        y_values, y_axis[:, np.newaxis], equal_nan=True
    ):
        return None

    return _cell_edges_from_axis(x_axis), _cell_edges_from_axis(y_axis)


def _center_grid_to_corner_grid(values: Any) -> np.ndarray:
    """Infer cell-corner coordinates from a 2D grid of cell-center coordinates."""
    array = np.asarray(values, dtype=float)
    if array.ndim != 2 or array.size == 0:
        raise ValueError("Corner inference requires a non-empty 2D center grid")

    ny, nx = array.shape
    padded = np.empty((ny + 2, nx + 2), dtype=float)
    padded[1:-1, 1:-1] = array

    if nx > 1:
        padded[1:-1, 0] = 2.0 * array[:, 0] - array[:, 1]
        padded[1:-1, -1] = 2.0 * array[:, -1] - array[:, -2]
    else:
        padded[1:-1, 0] = array[:, 0]
        padded[1:-1, -1] = array[:, 0]

    if ny > 1:
        padded[0, 1:-1] = 2.0 * array[0, :] - array[1, :]
        padded[-1, 1:-1] = 2.0 * array[-1, :] - array[-2, :]
    else:
        padded[0, 1:-1] = array[0, :]
        padded[-1, 1:-1] = array[0, :]

    if ny > 1 and nx > 1:
        padded[0, 0] = 3.0 * array[0, 0] - array[0, 1] - array[1, 0]
        padded[0, -1] = 3.0 * array[0, -1] - array[0, -2] - array[1, -1]
        padded[-1, 0] = 3.0 * array[-1, 0] - array[-1, 1] - array[-2, 0]
        padded[-1, -1] = 3.0 * array[-1, -1] - array[-1, -2] - array[-2, -1]
    elif nx > 1:
        padded[0, 0] = 2.0 * array[0, 0] - array[0, 1]
        padded[0, -1] = 2.0 * array[0, -1] - array[0, -2]
        padded[-1, 0] = padded[0, 0]
        padded[-1, -1] = padded[0, -1]
    elif ny > 1:
        padded[0, 0] = 2.0 * array[0, 0] - array[1, 0]
        padded[-1, 0] = 2.0 * array[-1, 0] - array[-2, 0]
        padded[0, -1] = padded[0, 0]
        padded[-1, -1] = padded[-1, 0]
    else:
        padded[0, 0] = array[0, 0]
        padded[0, -1] = array[0, 0]
        padded[-1, 0] = array[0, 0]
        padded[-1, -1] = array[0, 0]

    return 0.25 * (
        padded[:-1, :-1] + padded[1:, :-1] + padded[:-1, 1:] + padded[1:, 1:]
    )


def _scatter_cell_geometry(x: Any, y: Any) -> dict[str, np.ndarray] | None:
    """Return cell geometry metadata for rectilinear or curvilinear grids."""
    rectilinear_edges = _rectilinear_cell_edges(x, y)
    if rectilinear_edges is not None:
        x_edges, y_edges = rectilinear_edges
        return {
            "kind": np.asarray("rectilinear"),
            "x_edges": x_edges,
            "y_edges": y_edges,
        }

    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    if x_values.ndim != 2 or y_values.ndim != 2 or x_values.shape != y_values.shape:
        return None

    return {
        "kind": np.asarray("curvilinear"),
        "x_corners": _center_grid_to_corner_grid(x_values),
        "y_corners": _center_grid_to_corner_grid(y_values),
    }


def _coerce_matching_plot_field(
    values: Any,
    expected_shape: CoordinateShape,
) -> tuple[np.ndarray | None, bool]:
    """Return a 2D scalar field when the input matches the vector-grid shape."""
    if values is None or isinstance(values, str) or np.isscalar(values):
        return None, False

    array = np.asarray(values)
    if array.shape == expected_shape:
        return _filled_float_array(array), True
    return None, array.ndim >= 2


def _extract_meshgrid_axes(x: Any, y: Any) -> tuple[np.ndarray, np.ndarray]:
    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)

    if x_values.ndim == 1 and y_values.ndim == 1:
        return x_values, y_values

    if x_values.ndim != 2 or y_values.ndim != 2:
        raise ValueError(
            "allow_non_uniform_grid requires 1D axes or meshgrid-like 2D x/y coordinates"
        )
    if x_values.shape != y_values.shape:
        raise ValueError("2D x and y coordinates must have the same shape")

    x_axis = np.asarray(x_values[0, :], dtype=float)
    y_axis = np.asarray(y_values[:, 0], dtype=float)
    if not np.allclose(x_values, x_axis[np.newaxis, :], equal_nan=True):
        raise ValueError(
            "2D x coordinates must be meshgrid-like when allow_non_uniform_grid=True"
        )
    if not np.allclose(y_values, y_axis[:, np.newaxis], equal_nan=True):
        raise ValueError(
            "2D y coordinates must be meshgrid-like when allow_non_uniform_grid=True"
        )
    return x_axis, y_axis


def _normalize_regular_grid_orientation(
    x: Any,
    y: Any,
    u: Any,
    v: Any,
    color: Any = None,
    linewidth: Any = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, Any, Any]:
    """Normalize descending rectilinear axes to ascending order."""
    try:
        x_axis, y_axis = _extract_meshgrid_axes(x, y)
    except ValueError:
        return x, y, u, v, color, linewidth

    u_values = _filled_float_array(u)
    v_values = _filled_float_array(v)
    expected_shape = (y_axis.size, x_axis.size)
    if u_values.shape != expected_shape or v_values.shape != expected_shape:
        return x, y, u, v, color, linewidth

    color_field, _ = _coerce_matching_plot_field(color, expected_shape)
    linewidth_field, _ = _coerce_matching_plot_field(linewidth, expected_shape)

    x_norm = np.asarray(x_axis, dtype=float)
    y_norm = np.asarray(y_axis, dtype=float)
    u_norm = u_values
    v_norm = v_values

    if x_norm.size > 1 and x_norm[0] > x_norm[-1]:
        x_norm = x_norm[::-1]
        u_norm = u_norm[:, ::-1]
        v_norm = v_norm[:, ::-1]
        if color_field is not None:
            color_field = color_field[:, ::-1]
        if linewidth_field is not None:
            linewidth_field = linewidth_field[:, ::-1]

    if y_norm.size > 1 and y_norm[0] > y_norm[-1]:
        y_norm = y_norm[::-1]
        u_norm = u_norm[::-1, :]
        v_norm = v_norm[::-1, :]
        if color_field is not None:
            color_field = color_field[::-1, :]
        if linewidth_field is not None:
            linewidth_field = linewidth_field[::-1, :]

    if color_field is not None:
        color = color_field
    if linewidth_field is not None:
        linewidth = linewidth_field

    return x_norm, y_norm, u_norm, v_norm, color, linewidth
