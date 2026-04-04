"""Internal scatter rendering helpers for display-space-thinned stippling."""

from __future__ import annotations

from typing import Any

import numpy as np

from .._core.thinning import (
    _map_ncl_display_points_to_viewport,
    _resolve_ncl_min_distance_fraction,
)
from .._shared import coords as _coord_helpers
from .._shared.style import _resolve_scatter_aliases
from ..vector import _thin_ncl_mapped_candidates


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

    array = _coord_helpers._normalized_field_array(value, target_dims=target_dims)
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
