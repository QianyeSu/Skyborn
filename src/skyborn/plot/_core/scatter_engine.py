"""Internal scatter rendering helpers for display-space-thinned stippling."""

from __future__ import annotations

from typing import Any

import numpy as np

from .._core.thinning import (
    _map_ncl_display_points_to_viewport,
    _resolve_ncl_min_distance_fraction,
    _thin_ncl_mapped_candidates,
)
from .._shared import coords as _coord_helpers
from .._shared.style import _resolve_scatter_aliases


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
    grid_like: bool,
    where: Any,
    mask: Any,
    distance: float | None,
    cell_edges: tuple[np.ndarray, np.ndarray] | None,
) -> str:
    """Resolve how gridded stipple candidates should be positioned."""
    if placement is None:
        placement_name = "auto"
    else:
        placement_name = str(placement).strip().lower()

    if placement_name not in {"auto", "points", "cells"}:
        raise ValueError("placement must be one of 'auto', 'points', or 'cells'")

    supports_cells = cell_edges is not None
    if placement_name == "auto":
        if (
            grid_like
            and supports_cells
            and (where is not None or mask is not None)
            and distance is None
        ):
            return "cells"
        return "points"

    if placement_name == "cells" and not supports_cells:
        raise ValueError(
            "placement='cells' requires 1D axes or meshgrid-like rectilinear 2D x/y coordinates"
        )

    return placement_name


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
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    spacing_fraction: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Expand selected rectilinear cells into interior stipple candidates.

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
    x0 = np.asarray(x_edges[ix], dtype=float)
    x1 = np.asarray(x_edges[ix + 1], dtype=float)
    y0 = np.asarray(y_edges[iy], dtype=float)
    y1 = np.asarray(y_edges[iy + 1], dtype=float)

    corner_points = np.column_stack(
        [
            np.concatenate([x0, x1, x0, x1]),
            np.concatenate([y0, y0, y1, y1]),
        ]
    )
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

    corner_valid = np.isfinite(mapped_corners).all(axis=(1, 2))
    widths = np.zeros(source_indices.size, dtype=float)
    heights = np.zeros(source_indices.size, dtype=float)
    if np.any(corner_valid):
        widths[corner_valid] = np.ptp(mapped_corners[corner_valid, :, 0], axis=1)
        heights[corner_valid] = np.ptp(mapped_corners[corner_valid, :, 1], axis=1)

    spacing = max(
        float(spacing_fraction if spacing_fraction is not None else 0.03), 1e-6
    )
    sub_x = np.clip(np.ceil(widths / spacing).astype(int), 1, 12)
    sub_y = np.clip(np.ceil(heights / spacing).astype(int), 1, 12)

    point_chunks: list[np.ndarray] = []
    source_pos_chunks: list[np.ndarray] = []
    for source_pos in range(source_indices.size):
        if not corner_valid[source_pos]:
            point_chunks.append(source_points[source_pos][np.newaxis, :])
            source_pos_chunks.append(np.array([source_pos], dtype=int))
            continue

        x_edges_sub = np.linspace(
            float(x0[source_pos]),
            float(x1[source_pos]),
            int(sub_x[source_pos]) + 1,
            dtype=float,
        )
        y_edges_sub = np.linspace(
            float(y0[source_pos]),
            float(y1[source_pos]),
            int(sub_y[source_pos]) + 1,
            dtype=float,
        )
        xs = 0.5 * (x_edges_sub[:-1] + x_edges_sub[1:])
        ys = 0.5 * (y_edges_sub[:-1] + y_edges_sub[1:])
        grid_x, grid_y = np.meshgrid(xs, ys, indexing="xy")
        cell_points = np.column_stack([grid_x.ravel(), grid_y.ravel()])
        point_chunks.append(cell_points)
        source_pos_chunks.append(np.full(cell_points.shape[0], source_pos, dtype=int))

    return np.concatenate(point_chunks, axis=0), np.concatenate(source_pos_chunks)


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

    x_values, y_values, candidate_shape, grid_like, target_dims = (
        _coord_helpers._normalize_coordinates(
            x,
            y,
            reference_fields=(where, mask, s, c, *tuple(kwargs.values())),
        )
    )
    selection = _coord_helpers._normalize_selection_mask(
        where=where,
        mask=mask,
        candidate_shape=candidate_shape,
        target_dims=target_dims,
    )
    cell_edges = _coord_helpers._rectilinear_cell_edges(x_values, y_values)

    x_flat = np.ravel(np.asarray(x_values, dtype=float))
    y_flat = np.ravel(np.asarray(y_values, dtype=float))
    source_valid = selection.ravel() & np.isfinite(x_flat) & np.isfinite(y_flat)
    source_indices = np.flatnonzero(source_valid)

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
    resolved_placement = _resolve_placement(
        placement,
        grid_like=grid_like,
        where=where,
        mask=mask,
        distance=distance,
        cell_edges=cell_edges,
    )

    source_points = np.column_stack([x_flat[source_indices], y_flat[source_indices]])
    if transform is ax.transData:
        if resolved_placement == "cells" and cell_edges is not None:
            x_edges, y_edges = cell_edges
            extent_points = np.array(
                [
                    [float(np.nanmin(x_edges)), float(np.nanmin(y_edges))],
                    [float(np.nanmax(x_edges)), float(np.nanmax(y_edges))],
                ],
                dtype=float,
            )
            extent_points = extent_points[np.isfinite(extent_points).all(axis=1)]
            if len(extent_points) > 0:
                ax.update_datalim(extent_points)
                ax.autoscale_view()
        else:
            ax.update_datalim(source_points)
            ax.autoscale_view()

    if resolved_placement == "cells" and cell_edges is not None:
        x_edges, y_edges = cell_edges
        candidate_points, candidate_source_positions = _generate_cell_candidates(
            ax=ax,
            transform=transform,
            source_points=source_points,
            source_indices=source_indices,
            candidate_shape=candidate_shape,
            x_edges=x_edges,
            y_edges=y_edges,
            spacing_fraction=spacing_fraction,
        )
    else:
        candidate_points = source_points
        candidate_source_positions = np.arange(source_indices.size, dtype=int)

    display_points = np.asarray(transform.transform(candidate_points), dtype=float)
    display_valid = np.isfinite(display_points).all(axis=1)

    candidate_points = candidate_points[display_valid]
    candidate_source_positions = candidate_source_positions[display_valid]
    display_points = display_points[display_valid]

    if spacing_fraction is None or candidate_points.shape[0] <= 1:
        retained_indices = np.arange(candidate_points.shape[0], dtype=int)
    else:
        mapped_points = _map_ncl_display_points_to_viewport(display_points, ax.bbox)
        retained_indices = np.asarray(
            _thin_ncl_mapped_candidates(mapped_points, spacing_fraction),
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
