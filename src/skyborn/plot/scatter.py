"""Scatter plotting and public API for display-space-thinned stippling.

Author: Qianye Su <suqianye2000@gmail.com>
Copyright (c) 2025-2026 Qianye Su
Created: 2026-04-03
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from .ncl_vector import _is_cartopy_crs_like, _looks_like_axes
from .vector_plot import (
    _map_ncl_display_points_to_viewport,
    _resolve_ncl_min_distance_fraction,
    _thin_ncl_mapped_candidates,
)

__all__ = ["scatter"]


def _filled_array(value: Any, *, dtype=float) -> np.ndarray:
    """Return a plain ndarray with masked values replaced by NaN."""
    array = np.ma.asarray(value, dtype=dtype)
    if np.ma.isMaskedArray(array):
        return np.asarray(array.filled(np.nan), dtype=dtype)
    return np.asarray(array, dtype=dtype)


def _subset_ready_array(value: Any) -> np.ndarray:
    """Normalize one style array so it can be subset safely by retained points."""
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


def _transpose_dataarray_if_possible(
    value: Any, target_dims: tuple[Any, ...] | None
) -> Any:
    """Align a DataArray to the target plotting dimension order when possible."""
    da = _maybe_squeezed_dataarray(value)
    if da is None or target_dims is None or da.ndim != len(target_dims):
        return da if da is not None else value
    if da.dims == target_dims or set(da.dims) != set(target_dims):
        return da
    return da.transpose(*target_dims)


def _normalized_field_array(
    value: Any, target_dims: tuple[Any, ...] | None = None
) -> np.ndarray:
    """Return one style or mask field as an ndarray in plotting dimension order."""
    value = _transpose_dataarray_if_possible(value, target_dims)
    if isinstance(value, xr.DataArray):
        value = value.data
    return _subset_ready_array(value)


def _infer_grid_from_reference(
    value: Any,
    grid_shape: tuple[int, int],
    target_dims: tuple[Any, ...] | None,
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
) -> tuple[np.ndarray, np.ndarray, tuple[int, ...], bool, tuple[Any, ...] | None]:
    """Normalize paired points or grid coordinates into a consistent plotting form.

    This helper supports three common inputs:

    1. Paired 1D points: ``x[i], y[i]`` describe one point each.
    2. 1D grid axes: ``x`` and ``y`` describe a rectilinear grid and must be
       expanded with ``meshgrid`` before masking and thinning.
    3. 2D meshgrid-like coordinates: already expanded and ready to use.

    The ambiguous case is two 1D arrays. When their lengths differ, they are
    necessarily grid axes. When their lengths match, we inspect companion
    fields such as ``where`` / ``mask`` / ``s`` / ``c`` to see whether the user
    is actually passing grid-shaped metadata, in which case the 1D coordinates
    should also be interpreted as grid axes.
    """
    x_obj = _maybe_squeezed_dataarray(x)
    y_obj = _maybe_squeezed_dataarray(y)

    target_dims = None
    if x_obj is not None and y_obj is not None and x_obj.ndim == y_obj.ndim == 2:
        y_obj = _transpose_dataarray_if_possible(y_obj, tuple(x_obj.dims))
        target_dims = tuple(x_obj.dims)
    elif x_obj is not None and y_obj is not None and x_obj.ndim == y_obj.ndim == 1:
        target_dims = (y_obj.dims[0], x_obj.dims[0])

    x_values = _filled_array(x_obj.data if x_obj is not None else x, dtype=float)
    y_values = _filled_array(y_obj.data if y_obj is not None else y, dtype=float)

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
        # Two equally long 1D arrays could be paired points or grid axes. Use
        # the shapes of companion fields to disambiguate that case.
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
    candidate_shape: tuple[int, ...],
    target_dims: tuple[Any, ...] | None,
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

    # Numeric masks follow NumPy truthiness, but NaN should never create a
    # candidate point implicitly.
    finite = np.isfinite(selector_array)
    return finite & (selector_array != 0)


def _resolve_spacing_fraction(
    *,
    density: Any,
    distance: float | None,
    grid_like: bool,
) -> float | None:
    """Resolve the viewport-space spacing threshold used for thinning."""
    if distance is not None:
        return max(float(distance), 1e-6)
    if density is None:
        return _resolve_ncl_min_distance_fraction(1.0, None) if grid_like else None
    return _resolve_ncl_min_distance_fraction(density, None)


def _subset_scatter_value(
    value: Any,
    *,
    candidate_shape: tuple[int, ...],
    candidate_indices: np.ndarray,
    retained_indices: np.ndarray,
    target_dims: tuple[Any, ...] | None,
) -> Any:
    """Subset one scatter style argument to the points that survive thinning.

    Matplotlib accepts many array-like scatter kwargs such as ``s``, ``c``,
    ``linewidths``, and similar per-point arrays. This helper keeps those
    arrays aligned with the retained coordinates regardless of whether the
    original data was supplied:

    - on the full 2D grid,
    - already flattened to the full candidate count,
    - or already pre-filtered to the masked candidate count.
    """
    if value is None or isinstance(value, str) or np.isscalar(value):
        return value

    array = _normalized_field_array(value, target_dims=target_dims)
    full_count = int(np.prod(candidate_shape))
    candidate_count = int(candidate_indices.size)

    if array.shape == candidate_shape:
        return np.ravel(array)[candidate_indices][retained_indices]

    if array.ndim >= 1 and len(array) == full_count:
        return np.asarray(array)[candidate_indices][retained_indices]

    if (
        candidate_count != full_count
        and array.ndim >= 1
        and len(array) == candidate_count
    ):
        return np.asarray(array)[retained_indices]

    return value


def _subset_scatter_kwargs(
    kwargs: dict[str, Any],
    *,
    candidate_shape: tuple[int, ...],
    candidate_indices: np.ndarray,
    retained_indices: np.ndarray,
    target_dims: tuple[Any, ...] | None,
) -> dict[str, Any]:
    """Subset all array-like scatter kwargs together with the retained points."""
    return {
        key: _subset_scatter_value(
            value,
            candidate_shape=candidate_shape,
            candidate_indices=candidate_indices,
            retained_indices=retained_indices,
            target_dims=target_dims,
        )
        for key, value in kwargs.items()
    }


def _scatter_impl(
    ax,
    x: Any,
    y: Any,
    s: Any = None,
    c: Any = None,
    *,
    where: Any = None,
    mask: Any = None,
    density: Any = None,
    distance: float | None = None,
    min_distance: float | None = None,
    transform: Any = None,
    zorder: float | None = None,
    **kwargs: Any,
):
    """Scatter gridded or paired points with optional display-space thinning.

    The rendering pipeline is:

    1. Normalize coordinates into either paired points or an expanded grid.
    2. Build the candidate mask from ``where`` / ``mask`` plus finite-point
       checks on the coordinates themselves.
    3. Transform candidates into display space so thinning is controlled by
       what the viewer actually sees on the canvas, not by raw grid stride.
    4. Reuse the existing NCL-style viewport thinning helper to keep only one
       point inside each local spacing neighborhood.
    5. Subset all array-like scatter style arguments so they still match the
       retained coordinates exactly.
    """

    if transform is None:
        transform = ax.transData
    if distance is not None and min_distance is not None:
        raise TypeError("Use only one of 'distance' or 'min_distance'")
    if distance is None:
        distance = min_distance

    x_values, y_values, candidate_shape, grid_like, target_dims = (
        _normalize_coordinates(
            x,
            y,
            reference_fields=(where, mask, s, c, *tuple(kwargs.values())),
        )
    )
    selection = _normalize_selection_mask(
        where=where,
        mask=mask,
        candidate_shape=candidate_shape,
        target_dims=target_dims,
    )

    # Candidate extraction always happens on flattened arrays so the same code
    # path works for paired points, rectilinear grids, and 2D coordinates.
    x_flat = np.ravel(np.asarray(x_values, dtype=float))
    y_flat = np.ravel(np.asarray(y_values, dtype=float))
    valid = selection.ravel() & np.isfinite(x_flat) & np.isfinite(y_flat)
    candidate_indices = np.flatnonzero(valid)

    if candidate_indices.size == 0:
        empty_indices = np.empty(0, dtype=int)
        empty_s = _subset_scatter_value(
            s,
            candidate_shape=candidate_shape,
            candidate_indices=candidate_indices,
            retained_indices=empty_indices,
            target_dims=target_dims,
        )
        empty_c = _subset_scatter_value(
            c,
            candidate_shape=candidate_shape,
            candidate_indices=candidate_indices,
            retained_indices=empty_indices,
            target_dims=target_dims,
        )
        plot_kwargs = _subset_scatter_kwargs(
            kwargs,
            candidate_shape=candidate_shape,
            candidate_indices=candidate_indices,
            retained_indices=empty_indices,
            target_dims=target_dims,
        )
        if zorder is not None:
            plot_kwargs["zorder"] = zorder
        return ax.scatter(
            np.empty(0, dtype=float),
            np.empty(0, dtype=float),
            s=empty_s,
            c=empty_c,
            transform=transform,
            **plot_kwargs,
        )

    candidate_points = np.column_stack(
        [x_flat[candidate_indices], y_flat[candidate_indices]]
    )
    if transform is ax.transData:
        ax.update_datalim(candidate_points)
        ax.autoscale_view()

    # Thinning happens in display coordinates so the retained density follows
    # the current axes geometry and projection rather than array index spacing.
    display_points = np.asarray(transform.transform(candidate_points), dtype=float)
    display_valid = np.isfinite(display_points).all(axis=1)
    candidate_indices = candidate_indices[display_valid]
    candidate_points = candidate_points[display_valid]
    display_points = display_points[display_valid]

    spacing_fraction = _resolve_spacing_fraction(
        density=density,
        distance=distance,
        grid_like=grid_like,
    )
    if spacing_fraction is None or candidate_indices.size <= 1:
        retained_indices = np.arange(candidate_indices.size, dtype=int)
    else:
        mapped_points = _map_ncl_display_points_to_viewport(display_points, ax.bbox)
        retained_indices = np.asarray(
            _thin_ncl_mapped_candidates(mapped_points, spacing_fraction),
            dtype=int,
        )

    s_selected = _subset_scatter_value(
        s,
        candidate_shape=candidate_shape,
        candidate_indices=candidate_indices,
        retained_indices=retained_indices,
        target_dims=target_dims,
    )
    c_selected = _subset_scatter_value(
        c,
        candidate_shape=candidate_shape,
        candidate_indices=candidate_indices,
        retained_indices=retained_indices,
        target_dims=target_dims,
    )
    plot_kwargs = _subset_scatter_kwargs(
        kwargs,
        candidate_shape=candidate_shape,
        candidate_indices=candidate_indices,
        retained_indices=retained_indices,
        target_dims=target_dims,
    )
    if zorder is not None:
        plot_kwargs["zorder"] = zorder

    retained_points = candidate_points[retained_indices]
    return ax.scatter(
        retained_points[:, 0],
        retained_points[:, 1],
        s=s_selected,
        c=c_selected,
        transform=transform,
        **plot_kwargs,
    )


_array_scatter = _scatter_impl


def scatter(*args: Any, **kwargs: Any):
    """Scatter points with optional NCL-style display-space thinning.

    This is the public plotting entry point exposed by ``skyborn.plot``.
    It keeps a Matplotlib-compatible calling convention while adding
    display-space thinning for gridded stippling use cases such as
    significance masks on maps or vertical cross-sections.

    Supported call styles include ``scatter(ax, x, y, ...)``,
    ``scatter(x, y, ..., ax=ax)``, ``scatter(x, y, ...)``, and
    ``scatter(x, y, s, c, ...)``.

    Key arguments:

    - ``x`` and ``y`` may be paired 1D points, 1D rectilinear grid axes, or
      2D meshgrid-like coordinates.
    - ``where`` or ``mask`` can select candidate stipple points from a gridded
      field. Only one of them may be supplied.
    - ``density`` controls display-space thinning. Higher values retain more
      points. If it is omitted, gridded inputs use the default NCL-style
      spacing rule while paired 1D points keep all points by default.
    - ``distance`` is an explicit viewport-space thinning threshold.
      ``min_distance`` remains available as a backward-compatible alias.
    - ``transform`` accepts ordinary Matplotlib transforms and Cartopy CRS-like
      objects, which are converted to the matching Matplotlib transform
      automatically.
    """

    if not args:
        raise TypeError("scatter() expects at least x and y positional arguments")

    ax = kwargs.pop("ax", None)
    remaining_args = list(args)

    if remaining_args and _looks_like_axes(remaining_args[0]):
        if ax is not None:
            raise TypeError("scatter() received Axes both positionally and via ax=")
        ax = remaining_args.pop(0)

    if ax is None:
        ax = plt.gca()

    if len(remaining_args) < 2:
        raise TypeError("scatter() requires x and y arguments")
    if len(remaining_args) > 4:
        raise TypeError("scatter() received too many positional arguments")

    x = remaining_args.pop(0)
    y = remaining_args.pop(0)
    s = remaining_args.pop(0) if remaining_args else kwargs.pop("s", None)
    c = remaining_args.pop(0) if remaining_args else kwargs.pop("c", None)

    transform = kwargs.pop("transform", None)
    if _is_cartopy_crs_like(transform):
        # Keep the public API ergonomic for Cartopy users: they can pass a CRS
        # object just like in ordinary Matplotlib/Cartopy plotting code.
        transform = transform._as_mpl_transform(ax)

    return _array_scatter(
        ax,
        x,
        y,
        s=s,
        c=c,
        transform=transform,
        **kwargs,
    )
