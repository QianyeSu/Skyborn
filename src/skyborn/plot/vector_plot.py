"""Core NCL-like curly-vector rendering for 2D vector fields.

Author: Qianye Su <suqianye2000@gmail.com>
Copyright (c) 2025-2026 Qianye Su
Created: 2026-03-01 14:58:56
"""

from __future__ import annotations

import warnings
from typing import Any

import matplotlib as mpl
import matplotlib.collections as mcollections
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import numpy as np
from matplotlib import cm

from ._artists import vector_artists as _artist_helpers
from ._artists.vector_artists import _ncl_arrow_edge_size_px, _resolve_open_arrow_size
from ._core import geometry as _geometry
from ._core import legacy_stream as _legacy_stream
from ._core import native as _native_helpers
from ._core import sampling as _sampling
from ._core import thinning as _thinning
from ._core import vector_engine as _vector_engine
from ._core.result import CurlyVectorPlotSet
from ._core.vector_engine import (
    _curve_length_from_magnitude,
    _default_ncl_max_length_px,
    _finite_plot_field_values,
    _resolve_artist_coordinate_context,
    _resolve_ncl_length_scale,
    _resolve_ncl_reference_length_px,
)
from ._shared.coords import (
    _axis_coordinate_1d,
    _axis_is_uniform,
    _coerce_matching_plot_field,
    _extract_meshgrid_axes,
    _filled_float_array,
    _normalize_regular_grid_orientation,
)
from ._shared.style import (
    _collect_named_kwargs,
    _normalize_artist_alpha,
    _normalize_supported_arrowstyle,
    _resolve_curly_anchor_alias,
    _resolve_curly_style_aliases,
)

__all__ = ["curly_vector", "CurlyVectorPlotSet"]

_candidate_data_from_display_step = _geometry._candidate_data_from_display_step
_curve_shape_is_acceptable = _geometry._curve_shape_is_acceptable
_display_points_to_data = _geometry._display_points_to_data
_display_step_to_data = _geometry._display_step_to_data
_default_ncl_box_center_candidates = _vector_engine._default_ncl_box_center_candidates
_default_ncl_candidate_shape = _vector_engine._default_ncl_candidate_shape
_density_scalar = _vector_engine._density_scalar
_density_xy = _vector_engine._density_xy
_evaluate_ncl_display_curve = _geometry._evaluate_ncl_display_curve
_finite_difference_step = _geometry._finite_difference_step
_fit_single_bend_display_curve = _geometry._fit_single_bend_display_curve
DomainMap = _legacy_stream.DomainMap
interpgrid = _legacy_stream.interpgrid
InvalidIndexError = _legacy_stream.InvalidIndexError
_local_display_jacobian = _geometry._local_display_jacobian
_map_ncl_display_points_to_viewport = _thinning._map_ncl_display_points_to_viewport
_NCLDisplaySampler = _thinning._NCLDisplaySampler
_NCLNativeTraceContext = _thinning._NCLNativeTraceContext
_ncl_step_length_px = _vector_engine._ncl_step_length_px
OutOfBounds = _legacy_stream.OutOfBounds
_point_at_arc_distance_from_end = _geometry._point_at_arc_distance_from_end
_point_within_grid_data = _geometry._point_within_grid_data
_prepare_ncl_display_sampler = _thinning._prepare_ncl_display_sampler
_resolve_ncl_min_distance_fraction = _thinning._resolve_ncl_min_distance_fraction
_resolve_curly_anchor = _vector_engine._resolve_curly_anchor
_valid_ncl_center_candidates = _vector_engine._valid_ncl_center_candidates
_euler_step = _legacy_stream._euler_step
_gen_starting_points = _legacy_stream._gen_starting_points
_get_integrator = _legacy_stream._get_integrator
_integrate_rk12 = _legacy_stream._integrate_rk12
Grid = _vector_engine.Grid
StreamMask = _legacy_stream.StreamMask
TerminateTrajectory = _legacy_stream.TerminateTrajectory
_tip_display_geometry_from_display_curve = (
    _geometry._tip_display_geometry_from_display_curve
)
_trim_display_curve_from_end = _geometry._trim_display_curve_from_end
_corrected_ncl_display_origin = _vector_engine._corrected_ncl_display_origin
_clip_display_step_to_viewport = _vector_engine._clip_display_step_to_viewport

_ISSUED_NATIVE_WARNINGS: set[str] = set()
_NATIVE_IMPORT_ERROR: Exception | None = None
_CURLY_VECTOR_NCL_KWARG_NAMES = (
    "density",
    "linewidth",
    "color",
    "vmin",
    "vmax",
    "cmap",
    "norm",
    "alpha",
    "facecolor",
    "edgecolor",
    "rasterized",
    "arrowsize",
    "arrowstyle",
    "transform",
    "zorder",
    "start_points",
    "integration_direction",
    "grains",
    "broken_streamlines",
    "anchor",
    "ref_magnitude",
    "ref_length",
    "min_frac_length",
    "min_distance",
    "allow_non_uniform_grid",
    "ncl_preset",
)

try:
    from .nclcurly_native import sample_grid_field as _sample_grid_field_native
    from .nclcurly_native import (
        sample_grid_field_array as _sample_grid_field_array_native,
    )
    from .nclcurly_native import (
        thin_mapped_candidates as _thin_ncl_mapped_candidates_native,
    )
    from .nclcurly_native import trace_ncl_direction as _trace_ncl_direction_native
except Exception as err:
    _NATIVE_IMPORT_ERROR = err
    _sample_grid_field_native = None
    _sample_grid_field_array_native = None
    _thin_ncl_mapped_candidates_native = None
    _trace_ncl_direction_native = None


def _warn_native_once(key: str, message: str) -> None:
    if key in _ISSUED_NATIVE_WARNINGS:
        return
    _ISSUED_NATIVE_WARNINGS.add(key)
    warnings.warn(message, RuntimeWarning, stacklevel=3)


def _warn_if_native_backend_unavailable() -> None:
    if _NATIVE_IMPORT_ERROR is None:
        return
    detail = f"{type(_NATIVE_IMPORT_ERROR).__name__}: {_NATIVE_IMPORT_ERROR}"
    _warn_native_once(
        "native_backend_unavailable",
        "skyborn.plot curly_vector native backend is unavailable "
        f"({detail}); falling back to the Python implementation.",
    )


def _disable_native_helper(global_name: str, helper_label: str, err: Exception) -> None:
    globals()[global_name] = None
    detail = f"{type(err).__name__}: {err}"
    _warn_native_once(
        f"{global_name}_failed",
        "skyborn.plot curly_vector native "
        f"{helper_label} failed ({detail}); disabling that accelerated path "
        "for the rest of this session and falling back to Python.",
    )


def _normalize_ncl_preset(ncl_preset):
    """Normalize supported NCL-like preset names."""
    if ncl_preset is None:
        return None

    preset = str(ncl_preset).strip().lower()
    if preset in {"profile", "vertical_profile", "vertical-profile", "lat_pressure"}:
        return "profile"

    raise ValueError(f"Unsupported ncl_preset {ncl_preset!r}")


def _resolve_default_ncl_preset(x, y, allow_non_uniform_grid, ncl_preset):
    preset = _normalize_ncl_preset(ncl_preset)
    if preset is not None:
        return True if preset == "profile" else allow_non_uniform_grid, preset

    x_axis = _axis_coordinate_1d(x, "x")
    y_axis = _axis_coordinate_1d(y, "y")
    if x_axis is None or y_axis is None:
        return allow_non_uniform_grid, preset

    if not _axis_is_uniform(x_axis) or not _axis_is_uniform(y_axis):
        return True, "profile"
    return allow_non_uniform_grid, preset


def _infer_profile_ncl_ref_magnitude(u, v, percentile=97.0):
    """Estimate a conservative reference magnitude for lat-pressure sections."""
    magnitude = np.hypot(np.asarray(u, dtype=float), np.asarray(v, dtype=float))
    valid = magnitude[np.isfinite(magnitude)]
    valid = valid[valid > 0.0]

    if valid.size == 0:
        return None

    ref_magnitude = float(np.nanpercentile(valid, float(percentile)))
    if not np.isfinite(ref_magnitude) or ref_magnitude <= 0.0:
        ref_magnitude = float(np.nanmax(valid))

    return ref_magnitude if np.isfinite(ref_magnitude) and ref_magnitude > 0.0 else None


def _apply_ncl_preset_defaults(
    ncl_preset,
    allow_non_uniform_grid,
    ref_magnitude,
    ref_length,
    min_distance,
    u,
    v,
):
    """Apply optional NCL-style presets without changing the global defaults."""
    preset = _normalize_ncl_preset(ncl_preset)
    if preset != "profile":
        return allow_non_uniform_grid, ref_magnitude, ref_length, min_distance, preset

    allow_non_uniform_grid = True
    if ref_magnitude is None:
        ref_magnitude = _infer_profile_ncl_ref_magnitude(u, v)
    if ref_length is None:
        ref_length = 0.06

    return allow_non_uniform_grid, ref_magnitude, ref_length, min_distance, preset


def _regrid_non_uniform_vectors_to_uniform(
    x: Any, y: Any, u: Any, v: Any, *scalars: Any
) -> tuple[np.ndarray, ...]:
    try:
        from scipy.interpolate import RegularGridInterpolator
    except ImportError as err:
        raise ImportError(
            "scipy is required for non-uniform grid support. Please install scipy."
        ) from err

    x_axis, y_axis = _extract_meshgrid_axes(x, y)
    u_values = _filled_float_array(u)
    v_values = _filled_float_array(v)
    scalar_values = [_filled_float_array(field) for field in scalars]
    expected_shape = (y_axis.size, x_axis.size)
    if u_values.shape != expected_shape or v_values.shape != expected_shape:
        raise ValueError(
            f"u and v must match the non-uniform grid shape {expected_shape}"
        )
    for scalar_field in scalar_values:
        if scalar_field.shape != expected_shape:
            raise ValueError(
                "Non-uniform scalar style fields must match the source vector-grid "
                f"shape {expected_shape}"
            )

    x_sorted = x_axis.copy()
    y_sorted = y_axis.copy()
    u_sorted = u_values.copy()
    v_sorted = v_values.copy()
    scalar_sorted = [field.copy() for field in scalar_values]

    if x_sorted.size > 1 and np.any(np.diff(x_sorted) <= 0):
        x_idx = np.argsort(x_sorted)
        x_sorted = x_sorted[x_idx]
        u_sorted = u_sorted[:, x_idx]
        v_sorted = v_sorted[:, x_idx]
        scalar_sorted = [field[:, x_idx] for field in scalar_sorted]
    if y_sorted.size > 1 and np.any(np.diff(y_sorted) <= 0):
        y_idx = np.argsort(y_sorted)
        y_sorted = y_sorted[y_idx]
        u_sorted = u_sorted[y_idx, :]
        v_sorted = v_sorted[y_idx, :]
        scalar_sorted = [field[y_idx, :] for field in scalar_sorted]

    if x_sorted.size > 1 and np.any(np.diff(x_sorted) <= 0):
        raise ValueError("x coordinates must be strictly monotonic after sorting")
    if y_sorted.size > 1 and np.any(np.diff(y_sorted) <= 0):
        raise ValueError("y coordinates must be strictly monotonic after sorting")

    x_uniform = np.linspace(
        float(np.nanmin(x_sorted)), float(np.nanmax(x_sorted)), x_sorted.size
    )
    y_uniform = np.linspace(
        float(np.nanmin(y_sorted)), float(np.nanmax(y_sorted)), y_sorted.size
    )
    X_uniform, Y_uniform = np.meshgrid(x_uniform, y_uniform, indexing="xy")
    points = np.column_stack([Y_uniform.ravel(), X_uniform.ravel()])

    u_interp = RegularGridInterpolator(
        (y_sorted, x_sorted),
        u_sorted,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    v_interp = RegularGridInterpolator(
        (y_sorted, x_sorted),
        v_sorted,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    scalar_interps = [
        RegularGridInterpolator(
            (y_sorted, x_sorted),
            field,
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        for field in scalar_sorted
    ]

    u_uniform = u_interp(points).reshape(Y_uniform.shape)
    v_uniform = v_interp(points).reshape(Y_uniform.shape)
    scalar_uniform = [
        interp(points).reshape(Y_uniform.shape) for interp in scalar_interps
    ]
    return (x_uniform, y_uniform, u_uniform, v_uniform, *scalar_uniform)


def curly_vector(
    axes: Any,
    x: Any,
    y: Any,
    u: Any,
    v: Any,
    density: Any = 1,
    linewidth: Any = None,
    linewidths: Any = None,
    color: Any = None,
    c: Any = None,
    cmap: Any = None,
    norm: Any = None,
    vmin: float | None = None,
    vmax: float | None = None,
    alpha: float | None = None,
    facecolor: Any = None,
    facecolors: Any = None,
    edgecolor: Any = None,
    edgecolors: Any = None,
    rasterized: bool | None = None,
    arrowsize: float = 1,
    arrowstyle: str = "->",
    transform: Any = None,
    zorder: float | None = None,
    start_points: Any = None,
    integration_direction: str = "both",
    grains: Any = 15,
    broken_streamlines: bool = True,
    allow_non_uniform_grid: bool = False,
    anchor: str | None = None,
    pivot: str | None = None,
    ref_magnitude: float | None = None,
    ref_length: float | None = None,
    min_frac_length: float = 0.0,
    min_distance: float | None = None,
    ncl_preset: str | None = None,
) -> CurlyVectorPlotSet:
    """
    Draw NCL-like curved vector glyphs for a 2D vector flow.

    Parameters
    ----------
    x, y : 1D/2D arrays
        Evenly spaced strictly increasing arrays to make a grid.  If 2D, all
        rows of *x* must be equal and all columns of *y* must be equal; i.e.,
        they must be as if generated by ``np.meshgrid(x_1d, y_1d)``.
        For non-uniform grids (e.g., vertical profiles), set allow_non_uniform_grid=True.
    u, v : 2D arrays
        *x* and *y*-velocities. The number of rows and columns must match
        the length of *y* and *x*, respectively.
    density : float or (float, float)
        Controls the closeness of streamlines. When ``density = 1``, the domain
        is divided into a 30x30 grid. *density* linearly scales this grid.
        Each cell in the grid can have, at most, one traversing streamline.
        For different densities in each direction, use a tuple
        (density_x, density_y).
    linewidth : float or 2D array
        The width of the streamlines. With a 2D array the line width can be
        varied across the grid. The array must have the same shape as *u*
        and *v*.
    linewidths : float or 2D array, optional
        Matplotlib ``quiver``-style alias for ``linewidth``.
    color : color or 2D array
        The streamline color. If given an array, its values are converted to
        colors using *cmap* and *norm*.  The array must have the same shape
        as *u* and *v*.
    c : color or 2D array, optional
        Matplotlib-style alias for ``color``.
    cmap, norm
        Data normalization and colormapping parameters for *color*; only used
        if *color* is an array of floats. See `~.Axes.imshow` for a detailed
        description.
    vmin, vmax : float, optional
        Lower and upper normalization bounds used when ``color``/``c`` is a
        scalar field and ``norm`` is omitted.
    alpha : float, optional
        Matplotlib artist alpha applied to both the curved shafts and the
        arrow heads.
    facecolor, edgecolor : color-like, optional
        Explicit arrow-head fill and edge colors, similar to
        ``matplotlib.pyplot.quiver``. These mainly affect the filled
        ``arrowstyle="-|>"`` head. When omitted, the resolved shaft color is
        reused. Open ``"->"`` heads remain line-based and therefore ignore
        ``facecolor``.
    facecolors, edgecolors : color-like, optional
        Matplotlib-style aliases for ``facecolor`` and ``edgecolor``.
    rasterized : bool, optional
        Whether to rasterize the generated curly-vector artists when exporting
        to vector formats such as PDF or SVG. This changes output rendering,
        not the underlying curly-vector algorithm.
    arrowsize : float
        Scaling factor for the arrow size.
    arrowstyle : str
        Supported arrow-head style. Use ``"->"`` for the open NCL-like line
        head or ``"-|>"`` for a filled triangular head.
    transform : Transform, optional
        Coordinate transformation for the plot. Defaults to axes.transData.
    zorder : float
        The zorder of the streamlines and arrows.
        Artists with lower zorder values are drawn first.
    start_points : (N, 2) array
        Coordinates of starting points for the streamlines in data coordinates
        (the same coordinates as the *x* and *y* arrays).
    integration_direction : {'forward', 'backward', 'both'}, default: 'both'
        Integrate the streamline in forward, backward or both directions.
    grains : int, default: 15
        Number of grains used in streamline integration.
    broken_streamlines : boolean, default: True
        If False, forces streamlines to continue until they
        leave the plot domain.  If True, they may be terminated if they
        come too close to another streamline.
    allow_non_uniform_grid : boolean, default: False
        If True, allows non-uniform grids like vertical profiles. The function
        will attempt to create a uniform interpolation grid for streamline calculation.
    anchor : {'tail', 'center', 'head'} or None, default: None
        Anchor point for the NCL-like curved-glyph renderer. If omitted, the anchor is
        inferred from ``integration_direction``: ``'forward'`` -> ``'tail'``,
        ``'backward'`` -> ``'head'``, ``'both'`` -> ``'center'``.
    pivot : {'tail', 'mid', 'middle', 'tip'} or None, default: None
        Matplotlib ``quiver``-style alias for ``anchor``. ``'mid'`` and
        ``'middle'`` map to ``'center'`` and ``'tip'`` maps to ``'head'``.
    ref_magnitude : float or None, default: None
        Reference magnitude used when mapping a
        physical vector magnitude to a display-space glyph length. If omitted,
        the maximum field magnitude is used.
    ref_length : float or None, default: None
        Reference glyph length as a fraction of the axes width for
        the NCL-like curved-glyph renderer. If omitted, a NCL-like default scaled by
        ``arrowsize`` is used.
    min_frac_length : float, default: 0.0
        Minimum glyph length as a fraction of the reference length for
        the NCL-like curved-glyph renderer.
    min_distance : float or None, default: None
        Minimum glyph-center spacing as a fraction of the axes width for
        the NCL-like curved-glyph renderer. If omitted, it is inferred from ``density``.
    ncl_preset : {None, 'profile'}, default: None
        Optional preset override for NCL-like glyph tuning. In most cases you
        can leave this as ``None``: regular lat-lon map grids keep the default
        map-style tuning, while non-uniform/profile-like grids are
        automatically promoted to the conservative ``'profile'`` preset.

    Returns
    -------
    CurlyVectorPlotSet
        Container object with attributes

        - ``lines``: `.LineCollection` of streamlines

        - ``arrows``: tuple of the actual filled arrow-head patches added to
          the axes. Open arrow styles use line segments only and therefore
          return an empty tuple.
    """
    color, linewidth, facecolor, edgecolor, vmin, vmax = _resolve_curly_style_aliases(
        color=color,
        c=c,
        linewidth=linewidth,
        linewidths=linewidths,
        facecolor=facecolor,
        facecolors=facecolors,
        edgecolor=edgecolor,
        edgecolors=edgecolors,
        norm=norm,
        vmin=vmin,
        vmax=vmax,
    )

    allow_non_uniform_grid, ncl_preset = _resolve_default_ncl_preset(
        x=x,
        y=y,
        allow_non_uniform_grid=allow_non_uniform_grid,
        ncl_preset=ncl_preset,
    )

    allow_non_uniform_grid, ref_magnitude, ref_length, min_distance, ncl_preset = (
        _apply_ncl_preset_defaults(
            ncl_preset=ncl_preset,
            allow_non_uniform_grid=allow_non_uniform_grid,
            ref_magnitude=ref_magnitude,
            ref_length=ref_length,
            min_distance=min_distance,
            u=u,
            v=v,
        )
    )
    arrowstyle = _normalize_supported_arrowstyle(arrowstyle)
    alpha = _normalize_artist_alpha(alpha)
    anchor = _resolve_curly_anchor_alias(anchor, pivot)

    if not allow_non_uniform_grid:
        x, y, u, v, color, linewidth = _normalize_regular_grid_orientation(
            x,
            y,
            u,
            v,
            color=color,
            linewidth=linewidth,
        )

    # Handle non-uniform grids by creating a uniform interpolation grid
    if allow_non_uniform_grid:
        expected_shape = np.shape(u)
        color_field, color_is_field = _coerce_matching_plot_field(color, expected_shape)
        linewidth_field, linewidth_is_field = _coerce_matching_plot_field(
            linewidth, expected_shape
        )
        if color_field is None and color_is_field:
            raise ValueError(
                "If 'color' is given, it must match the shape of the (x, y) grid"
            )
        if linewidth_field is None and linewidth_is_field:
            raise ValueError(
                "If 'linewidth' is given, it must match the shape of the (x, y) grid"
            )

        regridded = _regrid_non_uniform_vectors_to_uniform(
            x,
            y,
            u,
            v,
            *([field for field in (color_field, linewidth_field) if field is not None]),
        )
        x, y, u, v, *scalar_fields = regridded
        scalar_iter = iter(scalar_fields)
        if color_field is not None:
            color = next(scalar_iter)
        if linewidth_field is not None:
            linewidth = next(scalar_iter)

    return _curly_vector_ncl(
        axes,
        x,
        y,
        u,
        v,
        **_collect_named_kwargs(locals(), _CURLY_VECTOR_NCL_KWARG_NAMES),
    )


def _curly_vector_ncl(
    axes,
    x,
    y,
    u,
    v,
    density=1,
    linewidth=None,
    color=None,
    vmin=None,
    vmax=None,
    cmap=None,
    norm=None,
    alpha=None,
    facecolor=None,
    edgecolor=None,
    rasterized=None,
    arrowsize=1,
    arrowstyle="->",
    transform=None,
    zorder=None,
    start_points=None,
    integration_direction="both",
    grains=15,
    broken_streamlines=True,
    anchor=None,
    ref_magnitude=None,
    ref_length=None,
    min_frac_length=0.0,
    min_distance=None,
    allow_non_uniform_grid=False,
    ncl_preset=None,
):
    return _vector_engine._curly_vector_ncl_impl(
        axes=axes,
        x=x,
        y=y,
        u=u,
        v=v,
        density=density,
        linewidth=linewidth,
        color=color,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        norm=norm,
        alpha=alpha,
        facecolor=facecolor,
        edgecolor=edgecolor,
        rasterized=rasterized,
        arrowsize=arrowsize,
        arrowstyle=arrowstyle,
        transform=transform,
        zorder=zorder,
        start_points=start_points,
        integration_direction=integration_direction,
        grains=grains,
        broken_streamlines=broken_streamlines,
        anchor=anchor,
        ref_magnitude=ref_magnitude,
        ref_length=ref_length,
        min_frac_length=min_frac_length,
        min_distance=min_distance,
        allow_non_uniform_grid=allow_non_uniform_grid,
        ncl_preset=ncl_preset,
        warn_if_native_backend_unavailable_fn=_warn_if_native_backend_unavailable,
        grid_cls=Grid,
        prepare_ncl_display_sampler_fn=_prepare_ncl_display_sampler,
        prepare_ncl_native_trace_context_fn=_prepare_ncl_native_trace_context,
        select_ncl_centers_fn=_select_ncl_centers,
        build_ncl_curve_fn=_build_ncl_curve,
        sample_grid_field_fn=_sample_grid_field,
        build_ncl_arrow_artists_fn=_build_ncl_arrow_artists,
        display_points_to_data_fn=_display_points_to_data,
    )


def _prepare_ncl_native_trace_context(grid, u, v, viewport, display_sampler):
    if _trace_ncl_direction_native is None or display_sampler is None:
        return None
    try:
        return _NCLNativeTraceContext(
            grid=grid,
            u=u,
            v=v,
            viewport=viewport,
            display_sampler=display_sampler,
        )
    except Exception as err:
        _warn_native_once(
            "native_trace_context_failed",
            "skyborn.plot curly_vector native trace context setup failed "
            f"({type(err).__name__}: {err}); falling back to Python tracing.",
        )
        return None


def _select_ncl_centers(
    grid,
    magnitude,
    transform,
    axes,
    density,
    start_points,
    min_distance,
    display_sampler=None,
    ncl_preset=None,
):
    return _vector_engine._select_ncl_centers(
        grid=grid,
        magnitude=magnitude,
        transform=transform,
        axes=axes,
        density=density,
        start_points=start_points,
        min_distance=min_distance,
        display_sampler=display_sampler,
        ncl_preset=ncl_preset,
        sample_grid_field_array=_sample_grid_field_array,
        thin_ncl_mapped_candidates=_thin_ncl_mapped_candidates,
    )


def _prepare_ncl_center_candidates(grid, magnitude, density, start_points, ncl_preset):
    return _vector_engine._prepare_ncl_center_candidates(
        grid=grid,
        magnitude=magnitude,
        density=density,
        start_points=start_points,
        ncl_preset=ncl_preset,
        sample_grid_field_array=_sample_grid_field_array,
    )


def _thin_ncl_mapped_candidates(mapped_points, spacing_frac):
    mapped_points = np.asarray(mapped_points, dtype=float)
    if len(mapped_points) == 0:
        return []

    spacing_frac = max(float(spacing_frac), 1e-6)
    selected = _native_helpers._try_native_thin_ncl_mapped_candidates(
        native_thinner=_thin_ncl_mapped_candidates_native,
        mapped_points=mapped_points,
        spacing_frac=spacing_frac,
        on_error=_disable_native_helper,
    )
    if selected is not None:
        return selected

    spacing_sq = spacing_frac * spacing_frac
    bucket_scale = 1.0 / spacing_frac
    bucket_map = {}

    for idx, mapped_point in enumerate(mapped_points):
        bucket = tuple(np.floor(mapped_point * bucket_scale).astype(int))
        bucket_map.setdefault(bucket, []).append(idx)

    culled = np.zeros(len(mapped_points), dtype=bool)
    selected = []

    for idx, mapped_point in enumerate(mapped_points):
        if culled[idx]:
            continue

        selected.append(idx)
        bucket = tuple(np.floor(mapped_point * bucket_scale).astype(int))

        for ix in range(bucket[0] - 1, bucket[0] + 2):
            for iy in range(bucket[1] - 1, bucket[1] + 2):
                for other_idx in bucket_map.get((ix, iy), ()):
                    if other_idx <= idx or culled[other_idx]:
                        continue
                    offset = mapped_point - mapped_points[other_idx]
                    if float(np.dot(offset, offset)) < spacing_sq:
                        culled[other_idx] = True

    return selected


def _trace_ncl_curve(
    start_point,
    total_length_px,
    anchor,
    grid,
    u,
    v,
    transform,
    step_px,
    speed_scale,
    viewport,
    display_sampler=None,
    native_trace_context=None,
):
    return _vector_engine._trace_ncl_curve(
        start_point=start_point,
        total_length_px=total_length_px,
        anchor=anchor,
        grid=grid,
        u=u,
        v=v,
        transform=transform,
        step_px=step_px,
        speed_scale=speed_scale,
        viewport=viewport,
        display_sampler=display_sampler,
        native_trace_context=native_trace_context,
        trace_ncl_direction_fn=_trace_ncl_direction,
    )


def _build_ncl_curve(
    start_point,
    total_length_px,
    anchor,
    grid,
    u,
    v,
    transform,
    step_px,
    speed_scale,
    viewport,
    display_sampler=None,
    native_trace_context=None,
):
    return _vector_engine._build_ncl_curve(
        start_point=start_point,
        total_length_px=total_length_px,
        anchor=anchor,
        grid=grid,
        u=u,
        v=v,
        transform=transform,
        step_px=step_px,
        speed_scale=speed_scale,
        viewport=viewport,
        display_sampler=display_sampler,
        native_trace_context=native_trace_context,
        trace_ncl_curve_fn=_trace_ncl_curve,
        evaluate_ncl_display_curve_fn=_evaluate_ncl_display_curve,
    )


def _trace_ncl_direction_via_native(
    start_point,
    max_length_px,
    direction_sign,
    step_px,
    speed_scale,
    native_trace_context=None,
):
    return _native_helpers._try_native_trace_ncl_direction(
        native_tracer=_trace_ncl_direction_native,
        native_trace_context=native_trace_context,
        start_point=start_point,
        max_length_px=max_length_px,
        direction_sign=direction_sign,
        step_px=step_px,
        speed_scale=speed_scale,
        on_error=_disable_native_helper,
    )


def _trace_ncl_direction(
    start_point,
    max_length_px,
    direction_sign,
    grid,
    u,
    v,
    transform,
    step_px,
    speed_scale,
    viewport,
    display_sampler=None,
    native_trace_context=None,
):
    start_point = np.asarray(start_point, dtype=float)
    if native_trace_context is None and display_sampler is not None:
        native_trace_context = _prepare_ncl_native_trace_context(
            grid=grid,
            u=u,
            v=v,
            viewport=viewport,
            display_sampler=display_sampler,
        )
    native_curve = _trace_ncl_direction_via_native(
        start_point=start_point,
        max_length_px=max_length_px,
        direction_sign=direction_sign,
        step_px=step_px,
        speed_scale=speed_scale,
        native_trace_context=native_trace_context,
    )
    if native_curve is not None:
        return native_curve

    return _trace_ncl_direction_python(
        start_point=start_point,
        max_length_px=max_length_px,
        direction_sign=direction_sign,
        grid=grid,
        u=u,
        v=v,
        transform=transform,
        step_px=step_px,
        speed_scale=speed_scale,
        viewport=viewport,
        display_sampler=display_sampler,
    )


def _trace_ncl_direction_python(
    start_point,
    max_length_px,
    direction_sign,
    grid,
    u,
    v,
    transform,
    step_px,
    speed_scale,
    viewport,
    display_sampler=None,
):
    return _vector_engine._trace_ncl_direction_python(
        start_point=start_point,
        max_length_px=max_length_px,
        direction_sign=direction_sign,
        grid=grid,
        u=u,
        v=v,
        transform=transform,
        step_px=step_px,
        speed_scale=speed_scale,
        viewport=viewport,
        display_sampler=display_sampler,
        sample_local_vector_state_fn=_sample_local_vector_state,
        ncl_step_length_px_fn=_ncl_step_length_px,
        corrected_ncl_display_origin_fn=_corrected_ncl_display_origin,
        clip_display_step_to_viewport_fn=_clip_display_step_to_viewport,
        candidate_data_from_display_step_fn=_candidate_data_from_display_step,
        point_within_grid_data_fn=_point_within_grid_data,
    )


def _sample_local_vector_state(grid, u, v, transform, point, display_sampler=None):
    return _vector_engine._sample_local_vector_state(
        grid=grid,
        u=u,
        v=v,
        transform=transform,
        point=point,
        display_sampler=display_sampler,
        sample_grid_field_fn=_sample_grid_field,
        local_display_jacobian_fn=_local_display_jacobian,
    )


def _sample_grid_field(grid, field, xd, yd):
    value = _native_helpers._try_native_sample_grid_field(
        native_sampler=_sample_grid_field_native,
        grid=grid,
        field=field,
        xd=xd,
        yd=yd,
        on_error=_disable_native_helper,
    )
    if value is not None:
        return value

    return _sampling._sample_grid_field_python(
        grid=grid,
        field=field,
        xd=xd,
        yd=yd,
        interpgrid_fn=interpgrid,
        terminate_trajectory_exc=TerminateTrajectory,
    )


def _sample_grid_field_array(grid, field, points):
    points = np.asarray(points, dtype=float)
    if points.ndim == 1:
        points = points[np.newaxis, :]

    sampled = np.full(len(points), np.nan, dtype=float)
    if len(points) == 0:
        return sampled

    sampled_native = _native_helpers._try_native_sample_grid_field_array(
        native_sampler=_sample_grid_field_array_native,
        grid=grid,
        field=field,
        points=points,
        expected_shape=sampled.shape,
        on_error=_disable_native_helper,
    )
    if sampled_native is not None:
        return sampled_native

    return _sampling._sample_grid_field_array_python(
        grid=grid,
        field=field,
        points=points,
        interpgrid_fn=interpgrid,
    )


def _build_arrow_polygon(
    curve,
    grid,
    transform,
    head_length_px,
    head_width_px,
    facecolor,
    edgecolor,
    linewidth,
    alpha,
    zorder,
    display_curve=None,
    display_sampler=None,
    inverse_transform=None,
):
    return _artist_helpers._build_arrow_polygon(
        curve=curve,
        grid=grid,
        transform=transform,
        head_length_px=head_length_px,
        head_width_px=head_width_px,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        alpha=alpha,
        zorder=zorder,
        display_curve=display_curve,
        display_sampler=display_sampler,
        inverse_transform=inverse_transform,
        tip_display_geometry_fn=_tip_display_geometry,
        display_points_to_data_fn=_display_points_to_data,
    )


def _build_ncl_arrow_artists(
    curve,
    grid,
    transform,
    arrowstyle,
    head_length_px,
    head_width_px,
    facecolor,
    edgecolor,
    linewidth,
    alpha,
    zorder,
    display_curve=None,
    display_sampler=None,
    inverse_transform=None,
):
    return _artist_helpers._build_ncl_arrow_artists(
        curve=curve,
        grid=grid,
        transform=transform,
        arrowstyle=arrowstyle,
        head_length_px=head_length_px,
        head_width_px=head_width_px,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        alpha=alpha,
        zorder=zorder,
        display_curve=display_curve,
        display_sampler=display_sampler,
        inverse_transform=inverse_transform,
        uses_open_arrow_head_fn=_uses_open_arrow_head,
        build_open_arrow_segments_fn=_build_open_arrow_segments,
        build_arrow_polygon_fn=_build_arrow_polygon,
    )


def _uses_open_arrow_head(arrowstyle):
    return _artist_helpers._uses_open_arrow_head(arrowstyle)


def _trim_curve_for_open_head(curve, transform, head_length_px, display_sampler=None):
    curve = np.asarray(curve, dtype=float)
    if len(curve) < 2 or head_length_px <= 1e-6:
        return curve

    geometry = _open_arrow_geometry(
        curve,
        transform,
        head_length_px,
        display_sampler=display_sampler,
    )
    if geometry is None:
        return curve
    display_curve = geometry["display_curve"]

    trimmed_display = _trim_display_curve_from_end(display_curve, head_length_px)
    if trimmed_display is None or len(trimmed_display) < 2:
        return curve
    trimmed_display = trimmed_display.copy()
    trimmed_display[-1] = geometry["base_center_display"]

    try:
        trimmed_curve = transform.inverted().transform(trimmed_display)
    except Exception:
        return curve
    if not np.all(np.isfinite(trimmed_curve)):
        return curve
    return trimmed_curve


def _build_open_arrow_segments(
    curve,
    grid,
    transform,
    head_length_px,
    head_width_px,
    display_curve=None,
    display_sampler=None,
    inverse_transform=None,
):
    geometry = _open_arrow_geometry(
        curve,
        transform,
        head_length_px,
        head_width_px,
        display_curve=display_curve,
        display_sampler=display_sampler,
    )
    if geometry is None:
        return []

    display_vertices = np.vstack(
        [
            geometry["left_display"],
            geometry["tip_display"],
            geometry["right_display"],
        ]
    )

    data_vertices = _display_points_to_data(
        transform,
        display_vertices,
        inverse_transform=inverse_transform,
    )
    if data_vertices is None:
        return []

    left, tip_data, right = data_vertices
    return [
        np.vstack([left, tip_data]),
        np.vstack([right, tip_data]),
    ]


def _open_arrow_geometry(
    curve,
    transform,
    head_length_px,
    head_width_px=None,
    display_curve=None,
    display_sampler=None,
):
    curve = np.asarray(curve, dtype=float)
    if len(curve) < 2:
        return None

    if display_curve is None:
        try:
            display_curve = transform.transform(curve)
        except Exception:
            return None
        if not np.all(np.isfinite(display_curve)):
            return None
    else:
        display_curve = np.asarray(display_curve, dtype=float)
        if not np.all(np.isfinite(display_curve)):
            return None

    tip_geometry = _tip_display_geometry_from_display_curve(
        display_curve, head_length_px * 1.35
    )
    if tip_geometry is None:
        return None
    tip_display, unit = tip_geometry
    base_center = tip_display - unit * head_length_px

    geometry = {
        "display_curve": display_curve,
        "tip_display": tip_display,
        "base_center_display": base_center,
        "unit": unit,
    }
    if head_width_px is not None:
        normal = np.array([-unit[1], unit[0]])
        geometry["left_display"] = base_center + normal * head_width_px / 2.0
        geometry["right_display"] = base_center - normal * head_width_px / 2.0
    return geometry


def _tip_display_geometry(
    curve,
    transform,
    backoff_px,
    display_curve=None,
    display_sampler=None,
):
    curve = np.asarray(curve, dtype=float)
    if len(curve) < 2:
        return None

    if display_curve is None:
        try:
            display_curve = transform.transform(curve)
        except Exception:
            return None
        if not np.all(np.isfinite(display_curve)):
            return None
    else:
        display_curve = np.asarray(display_curve, dtype=float)
        if not np.all(np.isfinite(display_curve)):
            return None

    return _tip_display_geometry_from_display_curve(display_curve, backoff_px)
