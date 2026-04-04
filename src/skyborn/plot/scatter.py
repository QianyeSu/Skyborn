"""Scatter plotting and public API for display-space-thinned stippling.

Author: Qianye Su <suqianye2000@gmail.com>
Copyright (c) 2025-2026 Qianye Su
Created: 2026-04-03
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from ._shared import coords as _coord_helpers
from ._shared.axes import _is_cartopy_crs_like, _looks_like_axes
from ._shared.style import _resolve_scatter_aliases
from .vector_plot import (
    _map_ncl_display_points_to_viewport,
    _resolve_ncl_min_distance_fraction,
    _thin_ncl_mapped_candidates,
)

__all__ = ["scatter"]

_filled_array = _coord_helpers._filled_float_array
_subset_ready_array = _coord_helpers._subset_ready_array
_maybe_squeezed_dataarray = _coord_helpers._maybe_squeezed_dataarray
_transpose_dataarray_if_possible = _coord_helpers._transpose_dataarray_if_possible
_normalized_field_array = _coord_helpers._normalized_field_array
_infer_grid_from_reference = _coord_helpers._infer_grid_from_reference
_normalize_coordinates = _coord_helpers._normalize_coordinates
_normalize_selection_mask = _coord_helpers._normalize_selection_mask


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
    """Render scatter points after optional NCL-style display-space thinning.

    This is the core implementation behind the public :func:`scatter` wrapper.
    Its job is not just to call ``Axes.scatter``. It also has to reconcile the
    looser input conventions used by Skyborn plotting helpers:

    - ``x`` / ``y`` may describe paired 1D points,
    - or they may describe a rectilinear grid that should be expanded first,
    - while ``where`` / ``mask`` / ``s`` / ``c`` may arrive as NumPy arrays or
      ``xarray.DataArray`` objects that still follow the original grid layout.

    Once those inputs are normalized, the thinning itself is deliberately done
    in display space rather than data space. That is the key design choice:
    stipple density should follow what the viewer sees on the final axes or map
    projection, not the raw array stride of the source field. This matches the
    same idea used by the curly-vector implementation and avoids repeated manual
    ``[::step]`` tuning when the projection, extent, or panel aspect changes.
    """

    c, kwargs = _resolve_scatter_aliases(c, kwargs)

    if transform is None:
        transform = ax.transData
    if distance is not None and min_distance is not None:
        raise TypeError("Use only one of 'distance' or 'min_distance'")
    if distance is None:
        distance = min_distance

    # Normalize all supported coordinate styles into one common representation.
    # After this step, paired points and grid-like inputs both flow through the
    # same candidate-selection and thinning code path.
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

    # Candidate extraction always happens on flattened arrays. This keeps the
    # downstream logic identical for 1D paired points, 1D rectilinear axes that
    # were expanded to a grid, and already-meshed 2D coordinates.
    x_flat = np.ravel(np.asarray(x_values, dtype=float))
    y_flat = np.ravel(np.asarray(y_values, dtype=float))
    valid = selection.ravel() & np.isfinite(x_flat) & np.isfinite(y_flat)
    candidate_indices = np.flatnonzero(valid)

    if candidate_indices.size == 0:
        # Preserve Matplotlib behavior even when nothing survives selection:
        # return a real PathCollection and keep array-like style arguments in
        # the shape that an empty scatter call expects.
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
        # Ordinary Matplotlib scatter updates data limits automatically. Because
        # Skyborn may thin points before the final scatter call, we update the
        # limits from the full candidate cloud first so autoscaling still
        # reflects the user's actual input domain rather than the retained subset.
        ax.update_datalim(candidate_points)
        ax.autoscale_view()

    # Thinning happens in display coordinates so the retained density follows
    # the current axes geometry and projection rather than array index spacing.
    # This is what makes the result stable across different map projections,
    # panel sizes, and vertical cross-section aspect ratios.
    display_points = np.asarray(transform.transform(candidate_points), dtype=float)
    display_valid = np.isfinite(display_points).all(axis=1)

    # Projection transforms can legitimately map some candidate points outside
    # the drawable region or to non-finite coordinates. Those points cannot
    # participate in viewport thinning and would also fail in the final scatter.
    candidate_indices = candidate_indices[display_valid]
    candidate_points = candidate_points[display_valid]
    display_points = display_points[display_valid]

    spacing_fraction = _resolve_spacing_fraction(
        density=density,
        distance=distance,
        grid_like=grid_like,
    )
    if spacing_fraction is None or candidate_indices.size <= 1:
        # Paired point input keeps all points by default, and any explicit
        # thinning request is irrelevant when there is at most one candidate.
        retained_indices = np.arange(candidate_indices.size, dtype=int)
    else:
        # Convert from absolute display pixels into the normalized viewport
        # coordinate system expected by the shared NCL-style thinning helper.
        # The helper returns indices into the current candidate list, not the
        # original flattened grid.
        mapped_points = _map_ncl_display_points_to_viewport(display_points, ax.bbox)
        retained_indices = np.asarray(
            _thin_ncl_mapped_candidates(mapped_points, spacing_fraction),
            dtype=int,
        )

    # Any per-point style arrays must be subset with the same retained indices.
    # Otherwise the final ``ax.scatter`` call would receive coordinates and
    # style metadata with inconsistent lengths.
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

    # The final drawing step is just plain Matplotlib scatter. All of the
    # Skyborn-specific behavior lives in the normalization, thinning, and
    # style-subsetting stages above.
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
    It keeps a Matplotlib-compatible calling convention while adding the
    display-space thinning workflow needed for gridded stippling use cases such
    as significance masks on map projections and vertical cross-sections.

    Unlike a simple ``x[::step], y[::step]`` subsampling strategy, Skyborn
    first transforms candidate points into the current display geometry and then
    thins them there. That means the visible stipple density responds to the
    actual projection, axes aspect ratio, and panel extent seen by the user.

    Supported call styles
    ---------------------
    Matplotlib-style
        ``scatter(ax, x, y, ...)``
        ``scatter(x, y, ..., ax=ax)``
        ``scatter(x, y, ...)``
        ``scatter(x, y, s, c, ...)``

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Target axes. If omitted, ``matplotlib.pyplot.gca()`` is used.
    x, y : array-like
        Coordinate specification for the points to draw. Supported forms are:

        - paired 1D point coordinates, where ``x[i]`` and ``y[i]`` already
          describe one point each,
        - 1D rectilinear grid axes, which are expanded internally with
          ``numpy.meshgrid`` when a gridded mask or style field is supplied,
        - or 2D meshgrid-like coordinate arrays.
    s : scalar or array-like, optional
        Marker size passed through to ``matplotlib.axes.Axes.scatter``. If an
        array is supplied, it may be defined on the full grid, on the masked
        candidate set, or on the already flattened point list.
    c : color-like or array-like, optional
        Marker color argument passed through to ``Axes.scatter``. Array-like
        color fields follow the same subsetting rules as ``s``.
    color : color-like or array-like, optional
        Compatibility alias for ``c``. This is useful when matching plotting
        helpers that prefer the ``color=`` spelling while still preserving
        Matplotlib scatter semantics for numeric color fields.
    linewidth : scalar or array-like, optional
        Compatibility alias for ``linewidths``.
    facecolor : color-like or array-like, optional
        Compatibility alias for ``facecolors``.
    edgecolor : color-like or array-like, optional
        Compatibility alias for ``edgecolors``.
    linewidths, facecolors, edgecolors
        Native Matplotlib scatter styling arguments. Array-like values follow
        the same retained-point subsetting rules as ``s`` and ``c``.
    vmin, vmax : float, optional
        Colormap range controls forwarded to ``Axes.scatter``.
    where, mask : array-like, optional
        Candidate-selection mask. Use only one of them. For gridded stippling,
        these are typically 2D boolean or numeric fields aligned with the input
        grid. Numeric masks follow NumPy truthiness, while NaN values are
        treated as invalid rather than truthy.
    density : float, optional
        Relative stipple density. Larger values retain more points. When
        omitted, gridded inputs use the default NCL-style spacing rule and
        paired 1D points keep all points.
    distance : float, optional
        Explicit display-space thinning distance in normalized viewport units.
        Use this when an exact spacing threshold is preferred over the relative
        ``density`` control.
    min_distance : float, optional
        Backward-compatible alias for ``distance``.
    transform : optional
        Source coordinate transform. Standard Matplotlib transforms are passed
        through directly. Cartopy CRS-like objects are converted to the
        matching Matplotlib transform automatically.
    zorder : float, optional
        Matplotlib z-order of the generated scatter collection.
    **kwargs
        Additional keyword arguments forwarded to ``Axes.scatter`` after any
        array-like values have been subset to the retained points.

    Returns
    -------
    matplotlib.collections.PathCollection
        The scatter collection returned by ``Axes.scatter``.
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
