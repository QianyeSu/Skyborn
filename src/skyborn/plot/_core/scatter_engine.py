"""Internal scatter rendering helpers for display-space-thinned stippling."""

from __future__ import annotations

from typing import Any

import numpy as np

from .._core.native import (
    _call_native_generate_cell_candidates,
    _call_native_thin_ncl_display_candidates,
)
from .._core.thinning import (
    _map_ncl_display_points_to_viewport,
    _resolve_ncl_min_distance_fraction,
)
from .._shared import coords as _coord_helpers
from .._shared.style import _resolve_scatter_aliases
from ..nclcurly_core import generate_cell_candidates as _generate_cell_candidates_native
from ..nclcurly_core import thin_display_candidates as _thin_display_candidates


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


def _resolve_placement(
    placement: Any,
    *,
    cell_geometry: dict[str, np.ndarray] | None,
) -> str:
    """Resolve how gridded stipple candidates should be positioned."""
    placement_name = _normalize_placement_name(placement)
    if placement_name == "cells" and cell_geometry is None:
        raise ValueError(
            "placement='cells' requires gridded coordinates with inferable "
            "cell geometry"
        )

    return placement_name


def _normalize_placement_name(placement: Any) -> str:
    """Normalize the requested candidate-placement mode."""
    if placement is None:
        placement_name = "auto"
    else:
        placement_name = str(placement).strip().lower()

    if placement_name not in {"auto", "points", "cells"}:
        raise ValueError("placement must be one of 'auto', 'points', or 'cells'")

    if placement_name == "auto":
        # NCL polymarker/scatter primitives draw the supplied marker
        # coordinates directly. Cell-interior stippling is a separate fill-like
        # mode in Skyborn, so callers must request it explicitly.
        return "points"

    return placement_name


def _as_1d_axis(values: Any) -> tuple[np.ndarray, str | None] | None:
    """Return a 1D numeric axis and optional DataArray dimension name."""
    data_array = _coord_helpers._maybe_squeezed_dataarray(values)
    if data_array is not None:
        if data_array.ndim != 1:
            return None
        return _coord_helpers._filled_float_array(data_array.data), data_array.dims[0]

    axis = _coord_helpers._filled_float_array(values)
    if axis.ndim != 1:
        return None
    return axis, None


def _field_implies_grid(
    value: Any,
    *,
    grid_shape: tuple[int, int],
    target_dims: tuple[Any, ...] | None,
) -> bool:
    """Return whether a companion field requests 1D axes as a rectilinear grid."""
    if value is None:
        return False
    array = _coord_helpers._normalized_field_array(value, target_dims=target_dims)
    return array.ndim == 2 and array.shape == grid_shape


def _prepare_1d_grid_candidates(
    x: Any,
    y: Any,
    *,
    reference_fields: tuple[Any, ...],
    where: Any,
    mask: Any,
    resolved_placement: str,
) -> dict[str, Any] | None:
    """Prepare selected points from 1D rectilinear axes without a full meshgrid."""
    x_axis_info = _as_1d_axis(x)
    y_axis_info = _as_1d_axis(y)
    if x_axis_info is None or y_axis_info is None:
        return None

    x_axis, x_dim = x_axis_info
    y_axis, y_dim = y_axis_info
    target_dims = (y_dim, x_dim) if x_dim is not None and y_dim is not None else None
    candidate_shape = (int(y_axis.size), int(x_axis.size))
    grid_like = x_axis.size != y_axis.size or any(
        _field_implies_grid(
            field,
            grid_shape=candidate_shape,
            target_dims=target_dims,
        )
        for field in reference_fields
    )
    if not grid_like:
        return None

    selection = _coord_helpers._normalize_selection_mask(
        where=where,
        mask=mask,
        candidate_shape=candidate_shape,
        target_dims=target_dims,
    )
    source_indices = np.flatnonzero(selection.ravel())
    if source_indices.size == 0:
        source_points = np.empty((0, 2), dtype=float)
    else:
        iy, ix = np.unravel_index(source_indices, candidate_shape)
        finite = np.isfinite(x_axis[ix]) & np.isfinite(y_axis[iy])
        if not np.all(finite):
            source_indices = source_indices[finite]
            iy = iy[finite]
            ix = ix[finite]
        source_points = np.column_stack([x_axis[ix], y_axis[iy]])

    cell_geometry = (
        _coord_helpers._scatter_cell_geometry(x_axis, y_axis)
        if resolved_placement == "cells"
        else None
    )
    return {
        "candidate_shape": candidate_shape,
        "grid_like": True,
        "target_dims": target_dims,
        "cell_geometry": cell_geometry,
        "source_indices": source_indices,
        "source_points": source_points,
    }


def _subset_scatter_value(
    value: Any,
    *,
    candidate_shape: tuple[int, ...],
    source_indices: np.ndarray,
    candidate_source_positions: np.ndarray,
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

    array = _coord_helpers._normalized_field_array(value, target_dims=target_dims)
    full_count = int(np.prod(candidate_shape))
    source_count = int(source_indices.size)
    generated_count = int(candidate_source_positions.size)

    if array.shape == candidate_shape:
        source_values = np.ravel(array)[source_indices]
        return source_values[candidate_source_positions][retained_indices]

    if array.ndim >= 1 and len(array) == full_count:
        source_values = np.asarray(array)[source_indices]
        return source_values[candidate_source_positions][retained_indices]

    if source_count != full_count and array.ndim >= 1 and len(array) == source_count:
        return np.asarray(array)[candidate_source_positions][retained_indices]

    if (
        generated_count != source_count
        and array.ndim >= 1
        and len(array) == generated_count
    ):
        return np.asarray(array)[retained_indices]

    return value


def _subset_scatter_kwargs(
    kwargs: dict[str, Any],
    *,
    candidate_shape: tuple[int, ...],
    source_indices: np.ndarray,
    candidate_source_positions: np.ndarray,
    retained_indices: np.ndarray,
    target_dims: tuple[Any, ...] | None,
) -> dict[str, Any]:
    """Subset all array-like scatter kwargs together with the retained points."""
    return {
        key: _subset_scatter_value(
            value,
            candidate_shape=candidate_shape,
            source_indices=source_indices,
            candidate_source_positions=candidate_source_positions,
            retained_indices=retained_indices,
            target_dims=target_dims,
        )
        for key, value in kwargs.items()
    }


def _generate_cell_candidates(
    *,
    ax,
    transform: Any,
    source_points: np.ndarray,
    source_indices: np.ndarray,
    candidate_shape: tuple[int, ...],
    cell_geometry: dict[str, np.ndarray],
    spacing_fraction: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Expand selected cells into interior stipple candidates.

    NCL contour stippling behaves more like a fill pattern than a simple
    polymarker-on-node pass. To approximate that behavior while staying inside
    Matplotlib ``Axes.scatter(...)``, we first populate each selected cell with
    a small regular interior lattice and then run the usual display-space
    thinning pass on those candidates.
    """
    if source_indices.size == 0:
        empty = np.empty((0, 2), dtype=float)
        return empty, np.empty(0, dtype=int)

    iy, ix = np.unravel_index(source_indices, candidate_shape)
    kind = str(np.asarray(cell_geometry["kind"]).item())
    if kind == "rectilinear":
        x_edges = np.asarray(cell_geometry["x_edges"], dtype=float)
        y_edges = np.asarray(cell_geometry["y_edges"], dtype=float)
        corners = np.stack(
            [
                np.column_stack([x_edges[ix], y_edges[iy]]),
                np.column_stack([x_edges[ix + 1], y_edges[iy]]),
                np.column_stack([x_edges[ix + 1], y_edges[iy + 1]]),
                np.column_stack([x_edges[ix], y_edges[iy + 1]]),
            ],
            axis=1,
        )
    elif kind == "curvilinear":
        x_corners = np.asarray(cell_geometry["x_corners"], dtype=float)
        y_corners = np.asarray(cell_geometry["y_corners"], dtype=float)
        corners = np.stack(
            [
                np.column_stack([x_corners[iy, ix], y_corners[iy, ix]]),
                np.column_stack([x_corners[iy, ix + 1], y_corners[iy, ix + 1]]),
                np.column_stack([x_corners[iy + 1, ix + 1], y_corners[iy + 1, ix + 1]]),
                np.column_stack([x_corners[iy + 1, ix], y_corners[iy + 1, ix]]),
            ],
            axis=1,
        )
    else:
        raise ValueError(f"Unsupported cell geometry kind {kind!r}")

    corner_points = corners.reshape(-1, 2)
    try:
        display_corners = np.asarray(transform.transform(corner_points), dtype=float)
    except Exception:
        return source_points.copy(), np.arange(source_indices.size, dtype=int)

    display_corners = display_corners.reshape(4, source_indices.size, 2).transpose(
        1, 0, 2
    )
    mapped_corners = _map_ncl_display_points_to_viewport(
        display_corners.reshape(-1, 2),
        ax.bbox,
    ).reshape(source_indices.size, 4, 2)

    spacing = max(
        float(spacing_fraction if spacing_fraction is not None else 0.03), 1e-6
    )
    return _call_native_generate_cell_candidates(
        _generate_cell_candidates_native,
        corners,
        mapped_corners,
        source_points,
        spacing,
    )


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
    placement: Any = "auto",
    transform: Any = None,
    zorder: float | None = None,
    **kwargs: Any,
):
    """Render scatter points after optional NCL-style display-space thinning."""

    c, kwargs = _resolve_scatter_aliases(c, kwargs)

    if transform is None:
        transform = ax.transData
    if distance is not None and min_distance is not None:
        raise TypeError("Use only one of 'distance' or 'min_distance'")
    if distance is None:
        distance = min_distance

    resolved_placement = _normalize_placement_name(placement)
    reference_fields = (where, mask, s, c, *tuple(kwargs.values()))
    grid_state = _prepare_1d_grid_candidates(
        x,
        y,
        reference_fields=reference_fields,
        where=where,
        mask=mask,
        resolved_placement=resolved_placement,
    )
    if grid_state is None:
        x_values, y_values, candidate_shape, grid_like, target_dims = (
            _coord_helpers._normalize_coordinates(
                x,
                y,
                reference_fields=reference_fields,
            )
        )
        selection = _coord_helpers._normalize_selection_mask(
            where=where,
            mask=mask,
            candidate_shape=candidate_shape,
            target_dims=target_dims,
        )
        cell_geometry = (
            _coord_helpers._scatter_cell_geometry(x_values, y_values)
            if resolved_placement == "cells"
            else None
        )

        x_flat = np.ravel(np.asarray(x_values, dtype=float))
        y_flat = np.ravel(np.asarray(y_values, dtype=float))
        source_valid = selection.ravel() & np.isfinite(x_flat) & np.isfinite(y_flat)
        source_indices = np.flatnonzero(source_valid)
        source_points = np.column_stack(
            [x_flat[source_indices], y_flat[source_indices]]
        )
    else:
        candidate_shape = grid_state["candidate_shape"]
        grid_like = grid_state["grid_like"]
        target_dims = grid_state["target_dims"]
        cell_geometry = grid_state["cell_geometry"]
        source_indices = grid_state["source_indices"]
        source_points = grid_state["source_points"]

    resolved_placement = _resolve_placement(
        resolved_placement,
        cell_geometry=cell_geometry,
    )

    if source_indices.size == 0:
        empty_indices = np.empty(0, dtype=int)
        empty_s = _subset_scatter_value(
            s,
            candidate_shape=candidate_shape,
            source_indices=empty_indices,
            candidate_source_positions=empty_indices,
            retained_indices=empty_indices,
            target_dims=target_dims,
        )
        empty_c = _subset_scatter_value(
            c,
            candidate_shape=candidate_shape,
            source_indices=empty_indices,
            candidate_source_positions=empty_indices,
            retained_indices=empty_indices,
            target_dims=target_dims,
        )
        plot_kwargs = _subset_scatter_kwargs(
            kwargs,
            candidate_shape=candidate_shape,
            source_indices=empty_indices,
            candidate_source_positions=empty_indices,
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

    spacing_fraction = _resolve_spacing_fraction(
        density=density,
        distance=distance,
        grid_like=grid_like,
    )
    if transform is ax.transData:
        if resolved_placement == "cells" and cell_geometry is not None:
            kind = str(np.asarray(cell_geometry["kind"]).item())
            if kind == "rectilinear":
                extent_points = np.array(
                    [
                        [
                            float(np.nanmin(cell_geometry["x_edges"])),
                            float(np.nanmin(cell_geometry["y_edges"])),
                        ],
                        [
                            float(np.nanmax(cell_geometry["x_edges"])),
                            float(np.nanmax(cell_geometry["y_edges"])),
                        ],
                    ],
                    dtype=float,
                )
            else:
                extent_points = np.column_stack(
                    [
                        np.ravel(np.asarray(cell_geometry["x_corners"], dtype=float)),
                        np.ravel(np.asarray(cell_geometry["y_corners"], dtype=float)),
                    ]
                )
            extent_points = extent_points[np.isfinite(extent_points).all(axis=1)]
            if len(extent_points) > 0:
                ax.update_datalim(extent_points)
                ax.autoscale_view()
        else:
            ax.update_datalim(source_points)
            ax.autoscale_view()

    if resolved_placement == "cells" and cell_geometry is not None:
        candidate_points, candidate_source_positions = _generate_cell_candidates(
            ax=ax,
            transform=transform,
            source_points=source_points,
            source_indices=source_indices,
            candidate_shape=candidate_shape,
            cell_geometry=cell_geometry,
            spacing_fraction=spacing_fraction,
        )
    else:
        candidate_points = source_points
        candidate_source_positions = np.arange(source_indices.size, dtype=int)

    display_points = np.asarray(transform.transform(candidate_points), dtype=float)
    display_valid = np.isfinite(display_points).all(axis=1)

    if not np.all(display_valid):
        candidate_points = candidate_points[display_valid]
        candidate_source_positions = candidate_source_positions[display_valid]
        display_points = display_points[display_valid]

    if spacing_fraction is None or candidate_points.shape[0] <= 1:
        retained_indices = np.arange(candidate_points.shape[0], dtype=int)
    else:
        retained_indices = np.asarray(
            _call_native_thin_ncl_display_candidates(
                _thin_display_candidates,
                display_points,
                ax.bbox,
                spacing_fraction,
            ),
            dtype=int,
        )

    s_selected = _subset_scatter_value(
        s,
        candidate_shape=candidate_shape,
        source_indices=source_indices,
        candidate_source_positions=candidate_source_positions,
        retained_indices=retained_indices,
        target_dims=target_dims,
    )
    c_selected = _subset_scatter_value(
        c,
        candidate_shape=candidate_shape,
        source_indices=source_indices,
        candidate_source_positions=candidate_source_positions,
        retained_indices=retained_indices,
        target_dims=target_dims,
    )
    plot_kwargs = _subset_scatter_kwargs(
        kwargs,
        candidate_shape=candidate_shape,
        source_indices=source_indices,
        candidate_source_positions=candidate_source_positions,
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
