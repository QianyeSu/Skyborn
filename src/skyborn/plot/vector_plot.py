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
from matplotlib import cm, patches

from ._artists.vector_artists import _ncl_arrow_edge_size_px, _resolve_open_arrow_size
from ._core.result import CurlyVectorPlotSet
from ._core.thinning import (
    _display_jump_threshold,
    _NCLDisplaySampler,
    _NCLNativeTraceContext,
    _prepare_ncl_display_sampler,
)
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
    _warn_if_native_backend_unavailable()
    grid = Grid(x, y, allow_non_uniform=allow_non_uniform_grid)

    if zorder is None:
        zorder = mlines.Line2D.zorder
    if transform is None:
        transform = axes.transData
    if color is None:
        color = axes._get_lines.get_next_color()
    if linewidth is None:
        linewidth = mpl.rcParams["lines.linewidth"]

    artist_transform, artist_inverse_transform, bake_display_geometry = (
        _resolve_artist_coordinate_context(axes, transform)
    )

    # The NCL-like branch sizes glyphs in display space. Ensure the data limits
    # are initialized before sampling ``transData``; otherwise the default
    # 0..1 axes limits can shrink global lon/lat glyphs into near-invisible
    # dots before autoscaling happens at the end.
    axes.update_datalim(
        np.array(
            [
                [grid.x_origin, grid.y_origin],
                [grid.x_origin + grid.width, grid.y_origin + grid.height],
            ]
        )
    )
    axes.autoscale_view()

    if u.shape != grid.shape or v.shape != grid.shape:
        raise ValueError("'u' and 'v' must match the shape of the (x, y) grid")

    u = np.asarray(np.ma.masked_invalid(u).filled(np.nan), dtype=float)
    v = np.asarray(np.ma.masked_invalid(v).filled(np.nan), dtype=float)
    magnitude = np.hypot(u, v)
    valid_magnitude = magnitude[np.isfinite(magnitude)]

    resolved_anchor = _resolve_curly_anchor(anchor, integration_direction)

    if valid_magnitude.size == 0:
        lc = mcollections.LineCollection(
            [],
            transform=artist_transform,
            zorder=zorder,
            alpha=alpha,
        )
        if rasterized is not None:
            lc.set_rasterized(bool(rasterized))
        axes.add_collection(lc, autolim=False)
        return CurlyVectorPlotSet(
            lc,
            (),
            0.0,
            magnitude,
            zorder,
            artist_transform,
            axes,
            linewidth,
            color,
            cm._ensure_cmap(cmap) if cmap is not None else None,
            arrowsize,
            arrowstyle,
            start_points,
            integration_direction,
            grains,
            broken_streamlines,
            allow_non_uniform_grid,
            density=density,
            anchor=resolved_anchor,
            length_scale=None,
            rasterized=rasterized,
        )

    color_field, color_is_field = _coerce_matching_plot_field(color, grid.shape)
    if color_field is None and color_is_field:
        raise ValueError(
            "If 'color' is given, it must match the shape of the (x, y) grid"
        )
    line_width_field, linewidth_is_field = _coerce_matching_plot_field(
        linewidth, grid.shape
    )
    if line_width_field is None and linewidth_is_field:
        raise ValueError(
            "If 'linewidth' is given, it must match the shape of the (x, y) grid"
        )

    use_multicolor_lines = color_field is not None
    color_default = None
    line_width_default = linewidth
    if use_multicolor_lines:
        finite_color_values = _finite_plot_field_values(color_field, "color")
        color_default = float(np.mean(finite_color_values))
        if norm is None:
            norm = mcolors.Normalize(
                float(np.min(finite_color_values)) if vmin is None else float(vmin),
                float(np.max(finite_color_values)) if vmax is None else float(vmax),
            )
        cmap = cm._ensure_cmap(cmap)
    if line_width_field is not None:
        finite_line_width_values = _finite_plot_field_values(
            line_width_field, "linewidth"
        )
        line_width_default = float(np.mean(finite_line_width_values))

    default_max_length_px = _default_ncl_max_length_px(axes.bbox, density)
    requested_ref_length_px = (
        0.0 if ref_length is None else max(axes.bbox.width * float(ref_length), 1.0)
    )
    ref_mag = 0.0 if ref_magnitude is None else float(ref_magnitude)
    min_frac_length = float(np.clip(min_frac_length, 0.0, 1.0))
    min_mag = float(np.min(valid_magnitude))
    max_mag = float(np.max(valid_magnitude))
    ref_length_px = _resolve_ncl_reference_length_px(
        min_mag=min_mag,
        max_mag=max_mag,
        ref_mag=ref_mag,
        requested_ref_length_px=requested_ref_length_px,
        min_frac_length=min_frac_length,
        default_max_length_px=default_max_length_px,
    )
    length_scale = _resolve_ncl_length_scale(
        min_mag=min_mag,
        max_mag=max_mag,
        ref_mag=ref_mag,
        requested_ref_length_px=requested_ref_length_px,
        min_frac_length=min_frac_length,
        default_max_length_px=default_max_length_px,
    )
    ref_length_frac = ref_length_px / max(float(axes.bbox.width), 1.0)
    step_px = max(1.5, axes.bbox.width * 0.0045)
    arrow_min_edge_px = max(axes.bbox.width * 0.003 * max(float(arrowsize), 0.1), 1.2)
    arrow_max_edge_px = max(
        axes.bbox.width * 0.012 * max(float(arrowsize), 0.1), arrow_min_edge_px
    )
    display_sampler = _prepare_ncl_display_sampler(grid, transform)
    native_trace_context = _prepare_ncl_native_trace_context(
        grid=grid,
        u=u,
        v=v,
        viewport=axes.bbox,
        display_sampler=display_sampler,
    )

    selected_centers = _select_ncl_centers(
        grid=grid,
        magnitude=magnitude,
        transform=transform,
        axes=axes,
        density=density,
        start_points=start_points,
        min_distance=min_distance,
        display_sampler=display_sampler,
        ncl_preset=ncl_preset,
    )

    streamlines = []
    line_colors = []
    line_widths = []
    arrows = []

    for center, center_mag in selected_centers:
        target_length_px = _curve_length_from_magnitude(center_mag, length_scale)
        curve_result = _build_ncl_curve(
            start_point=np.asarray(center, dtype=float),
            total_length_px=target_length_px,
            anchor=resolved_anchor,
            grid=grid,
            u=u,
            v=v,
            transform=transform,
            step_px=step_px,
            speed_scale=max_mag,
            viewport=axes.bbox,
            display_sampler=display_sampler,
            native_trace_context=native_trace_context,
        )
        if curve_result is None:
            continue
        curve, display_curve = curve_result
        if len(curve) < 2:
            continue

        if line_width_field is not None:
            sampled_width = _sample_grid_field(
                grid, line_width_field, center[0], center[1]
            )
            current_linewidth = (
                float(sampled_width)
                if sampled_width is not None
                else line_width_default
            )
        else:
            current_linewidth = linewidth

        if use_multicolor_lines:
            sampled_color = _sample_grid_field(grid, color_field, center[0], center[1])
            if sampled_color is None:
                sampled_color = color_default
            curve_color = cmap(norm(sampled_color))
        else:
            curve_color = color

        head_facecolor = curve_color if facecolor is None else facecolor
        if edgecolor is None:
            head_edgecolor = curve_color
        elif isinstance(edgecolor, str) and edgecolor.strip().lower() == "face":
            head_edgecolor = head_facecolor
        else:
            head_edgecolor = edgecolor

        head_length_px, head_width_px = _resolve_open_arrow_size(
            _ncl_arrow_edge_size_px(
                center_mag,
                max_mag=max_mag,
                min_edge_px=arrow_min_edge_px,
                max_edge_px=arrow_max_edge_px,
            )
        )
        artist_curve = curve
        if bake_display_geometry:
            baked_curve = _display_points_to_data(
                artist_transform,
                display_curve,
                inverse_transform=artist_inverse_transform,
            )
            if baked_curve is not None and len(baked_curve) >= 2:
                artist_curve = baked_curve

        streamlines.append(artist_curve)
        if use_multicolor_lines:
            line_colors.append(curve_color)
        if line_width_field is not None:
            line_widths.append(current_linewidth)

        head_segments, arrow = _build_ncl_arrow_artists(
            curve=artist_curve,
            grid=grid,
            transform=artist_transform,
            arrowstyle=arrowstyle,
            head_length_px=head_length_px,
            head_width_px=head_width_px,
            facecolor=head_facecolor,
            edgecolor=head_edgecolor,
            linewidth=current_linewidth,
            alpha=alpha,
            zorder=zorder,
            inverse_transform=artist_inverse_transform,
            display_curve=display_curve,
        )
        if head_segments:
            streamlines.extend(head_segments)
            if use_multicolor_lines:
                line_colors.extend([curve_color] * len(head_segments))
            if line_width_field is not None:
                line_widths.extend([current_linewidth] * len(head_segments))
        if arrow is not None:
            arrows.append(arrow)

    line_kw = {"zorder": zorder}
    if use_multicolor_lines:
        line_kw["colors"] = line_colors
    else:
        line_kw["color"] = color

    if line_width_field is not None:
        line_kw["linewidths"] = line_widths
    else:
        line_kw["linewidth"] = linewidth
    if alpha is not None:
        line_kw["alpha"] = alpha

    lc = mcollections.LineCollection(streamlines, transform=artist_transform, **line_kw)
    if rasterized is not None:
        lc.set_rasterized(bool(rasterized))
    # The axes limits were already seeded from the full grid extent above, so
    # avoid Cartopy reprojecting every segment again just to recompute datalim.
    axes.add_collection(lc, autolim=False)

    for patch in arrows:
        if rasterized is not None:
            patch.set_rasterized(bool(rasterized))
        axes.add_patch(patch)

    axes.autoscale_view()
    return CurlyVectorPlotSet(
        lc,
        tuple(arrows),
        ref_length_frac,
        magnitude,
        zorder,
        artist_transform,
        axes,
        linewidth,
        color,
        cmap,
        arrowsize,
        arrowstyle,
        start_points,
        integration_direction,
        grains,
        broken_streamlines,
        allow_non_uniform_grid,
        density=density,
        anchor=resolved_anchor,
        ncl_preset=ncl_preset,
        length_scale=length_scale,
        rasterized=rasterized,
    )


def _resolve_curly_anchor(anchor, integration_direction):
    mpl._api.check_in_list(
        ["forward", "backward", "both"], integration_direction=integration_direction
    )
    if anchor is None:
        anchor = {
            "forward": "tail",
            "backward": "head",
            "both": "center",
        }[integration_direction]
    mpl._api.check_in_list(["tail", "center", "head"], anchor=anchor)
    return anchor


def _density_xy(density):
    density_xy = np.broadcast_to(np.asarray(density, dtype=float), 2).astype(float)
    return np.maximum(density_xy, 0.1)


def _density_scalar(density):
    if np.isscalar(density):
        return max(float(density), 0.1)
    density_xy = _density_xy(density)
    return float(np.mean(density_xy))


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
    candidates, candidate_magnitudes = _prepare_ncl_center_candidates(
        grid=grid,
        magnitude=magnitude,
        density=density,
        start_points=start_points,
        ncl_preset=ncl_preset,
    )
    if display_sampler is not None:
        display_points = display_sampler.sample_display_points(candidates)
    else:
        display_points = transform.transform(candidates)
    valid = _valid_ncl_center_candidates(
        grid=grid,
        candidates=candidates,
        candidate_magnitudes=candidate_magnitudes,
        display_points=display_points,
        start_points=start_points,
    )

    candidates = candidates[valid]
    display_points = display_points[valid]
    candidate_magnitudes = candidate_magnitudes[valid]
    if candidates.size == 0:
        return []

    mapped_points = _map_ncl_display_points_to_viewport(display_points, axes.bbox)
    spacing_frac = _resolve_ncl_min_distance_fraction(
        density=density,
        min_distance=min_distance,
        ncl_preset=ncl_preset,
    )
    selected_indices = _thin_ncl_mapped_candidates(mapped_points, spacing_frac)
    return [(candidates[idx], candidate_magnitudes[idx]) for idx in selected_indices]


def _prepare_ncl_center_candidates(grid, magnitude, density, start_points, ncl_preset):
    if start_points is None:
        candidates = _default_ncl_box_center_candidates(
            grid,
            density=density,
            ncl_preset=ncl_preset,
        )
        candidate_magnitudes = _sample_grid_field_array(grid, magnitude, candidates)
        return candidates, candidate_magnitudes

    candidates = np.asanyarray(start_points, dtype=float)
    if candidates.ndim != 2 or candidates.shape[1] != 2:
        raise ValueError("'start_points' must be an (N, 2) array")
    candidate_magnitudes = _sample_grid_field_array(grid, magnitude, candidates)
    return candidates, candidate_magnitudes


def _default_ncl_candidate_shape(grid, density):
    density_xy = _density_xy(density)
    candidate_nx = min(grid.nx - 1, max(int(np.ceil(30.0 * density_xy[0])), 1))
    candidate_ny = min(grid.ny - 1, max(int(np.ceil(30.0 * density_xy[1])), 1))
    return candidate_nx, candidate_ny


def _default_ncl_box_center_candidates(grid, density=1, ncl_preset=None):
    if grid.nx < 2 or grid.ny < 2:
        return np.empty((0, 2), dtype=float)

    candidate_nx, candidate_ny = _default_ncl_candidate_shape(grid, density)
    if candidate_nx == grid.nx - 1 and candidate_ny == grid.ny - 1:
        xs = grid.x_origin + (np.arange(grid.nx - 1) + 0.5) * grid.dx
        ys = grid.y_origin + (np.arange(grid.ny - 1) + 0.5) * grid.dy
    else:
        x_edges = np.linspace(
            grid.x_origin, grid.x_origin + grid.width, candidate_nx + 1
        )
        y_edges = np.linspace(
            grid.y_origin, grid.y_origin + grid.height, candidate_ny + 1
        )
        xs = 0.5 * (x_edges[:-1] + x_edges[1:])
        ys = 0.5 * (y_edges[:-1] + y_edges[1:])
    grid_x, grid_y = np.meshgrid(xs, ys, indexing="xy")
    return np.column_stack([grid_x.ravel(), grid_y.ravel()])


def _valid_ncl_center_candidates(
    grid, candidates, candidate_magnitudes, display_points, start_points
):
    valid = np.isfinite(display_points).all(axis=1) & np.isfinite(candidate_magnitudes)

    for idx, (xd, yd) in enumerate(candidates):
        if _point_within_grid_data(grid, np.array([xd, yd], dtype=float)):
            continue
        if start_points is not None:
            raise ValueError(f"Starting point ({xd}, {yd}) outside of data boundaries")
        valid[idx] = False

    return valid


def _map_ncl_display_points_to_viewport(display_points, viewport):
    width = max(float(viewport.width), 1.0)
    height = max(float(viewport.height), 1.0)
    mapped = np.empty_like(display_points, dtype=float)
    mapped[:, 0] = (display_points[:, 0] - float(viewport.x0)) / width
    mapped[:, 1] = (display_points[:, 1] - float(viewport.y0)) / height
    return mapped


def _resolve_ncl_min_distance_fraction(density, min_distance, ncl_preset=None):
    if min_distance is not None:
        return max(float(min_distance), 1e-6)
    spacing_frac = 0.9 / (30.0 * _density_scalar(density))
    if _normalize_ncl_preset(ncl_preset) == "profile":
        spacing_frac *= 0.6
    return spacing_frac


def _thin_ncl_mapped_candidates(mapped_points, spacing_frac):
    mapped_points = np.asarray(mapped_points, dtype=float)
    if len(mapped_points) == 0:
        return []

    spacing_frac = max(float(spacing_frac), 1e-6)
    if _thin_ncl_mapped_candidates_native is not None:
        try:
            selected = _thin_ncl_mapped_candidates_native(
                mapped_points=mapped_points,
                spacing_frac=spacing_frac,
            )
        except Exception as err:
            _disable_native_helper(
                "_thin_ncl_mapped_candidates_native",
                "candidate thinning",
                err,
            )
            selected = None
        if selected is not None:
            return np.asarray(selected, dtype=int).tolist()

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
    if total_length_px <= 0:
        return None

    if anchor == "center":
        backward = _trace_ncl_direction(
            start_point,
            total_length_px / 2.0,
            -1.0,
            grid,
            u,
            v,
            transform,
            step_px,
            speed_scale,
            viewport,
            display_sampler=display_sampler,
            native_trace_context=native_trace_context,
        )
        forward = _trace_ncl_direction(
            start_point,
            total_length_px / 2.0,
            1.0,
            grid,
            u,
            v,
            transform,
            step_px,
            speed_scale,
            viewport,
            display_sampler=display_sampler,
            native_trace_context=native_trace_context,
        )
        if backward is None and forward is None:
            return None
        if backward is None:
            return forward
        if forward is None:
            return backward[::-1]
        return np.vstack([backward[::-1], forward[1:]])

    if anchor == "tail":
        return _trace_ncl_direction(
            start_point,
            total_length_px,
            1.0,
            grid,
            u,
            v,
            transform,
            step_px,
            speed_scale,
            viewport,
            display_sampler=display_sampler,
            native_trace_context=native_trace_context,
        )

    backward = _trace_ncl_direction(
        start_point,
        total_length_px,
        -1.0,
        grid,
        u,
        v,
        transform,
        step_px,
        speed_scale,
        viewport,
        display_sampler=display_sampler,
        native_trace_context=native_trace_context,
    )
    if backward is None:
        return None
    return backward[::-1]


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
    current_length_px = float(total_length_px)

    for _ in range(4):
        curve = _trace_ncl_curve(
            start_point=start_point,
            total_length_px=current_length_px,
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
        )
        if curve is not None and len(curve) >= 2:
            display_curve, transform_failed = _evaluate_ncl_display_curve(
                curve,
                transform,
                viewport=viewport,
            )
            if display_curve is not None:
                return curve, display_curve
            if transform_failed:
                return curve, None

        current_length_px *= 0.78
        if current_length_px <= step_px:
            break

    return None


def _trace_ncl_direction_via_native(
    start_point,
    max_length_px,
    direction_sign,
    step_px,
    speed_scale,
    native_trace_context=None,
):
    if native_trace_context is None:
        return None

    try:
        curve = _trace_ncl_direction_native(
            u=native_trace_context.u,
            v=native_trace_context.v,
            display_grid=native_trace_context.display_grid,
            cell_valid=native_trace_context.cell_valid,
            x_origin=native_trace_context.x_origin,
            y_origin=native_trace_context.y_origin,
            dx=native_trace_context.dx,
            dy=native_trace_context.dy,
            start_x=float(start_point[0]),
            start_y=float(start_point[1]),
            max_length_px=float(max_length_px),
            direction_sign=float(direction_sign),
            step_px=float(step_px),
            speed_scale=float(speed_scale),
            viewport_x0=native_trace_context.viewport_x0,
            viewport_y0=native_trace_context.viewport_y0,
            viewport_x1=native_trace_context.viewport_x1,
            viewport_y1=native_trace_context.viewport_y1,
            max_steps=512,
        )
    except Exception as err:
        _disable_native_helper("_trace_ncl_direction_native", "trace", err)
        return None

    if curve is None:
        return None
    curve = np.asarray(curve, dtype=float)
    return (
        curve
        if curve.ndim == 2 and curve.shape[0] >= 2 and curve.shape[1] == 2
        else None
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
    start_point = np.asarray(start_point, dtype=float)
    initial_state = _sample_local_vector_state(
        grid,
        u,
        v,
        transform,
        start_point,
        display_sampler=display_sampler,
    )
    if initial_state is None:
        return None

    points = [start_point]
    current_data = start_point
    current_display = initial_state[0]
    previous_display = None
    travelled = 0.0

    for step_index in range(512):
        remaining = max_length_px - travelled
        if remaining <= 1e-6:
            break

        if step_index == 0:
            state = initial_state
        else:
            state = _sample_local_vector_state(
                grid,
                u,
                v,
                transform,
                current_data,
                display_sampler=display_sampler,
            )
        if state is None:
            break
        current_display, current_jacobian, current_direction, current_speed = state

        step_length_px = min(
            remaining,
            _ncl_step_length_px(step_px, current_speed, speed_scale),
        )
        if step_length_px <= 1e-6:
            break

        corrected_display = _corrected_ncl_display_origin(
            current_display, previous_display
        )
        candidate_display = (
            corrected_display + direction_sign * current_direction * step_length_px
        )
        candidate_display, clipped = _clip_display_step_to_viewport(
            corrected_display, candidate_display, viewport
        )
        candidate = _candidate_data_from_display_step(
            current_data=current_data,
            current_display=current_display,
            candidate_display=candidate_display,
            jacobian=current_jacobian,
            transform=transform,
        )
        if candidate is None or not _point_within_grid_data(grid, candidate):
            break

        next_state = _sample_local_vector_state(
            grid,
            u,
            v,
            transform,
            candidate,
            display_sampler=display_sampler,
        )
        if next_state is not None:
            _, next_jacobian, next_direction, next_speed = next_state
            average_direction = direction_sign * (current_direction + next_direction)
            average_norm = np.hypot(*average_direction)
            if average_norm > 1e-12:
                average_speed = 0.5 * (current_speed + next_speed)
                step_length_px = min(
                    remaining,
                    _ncl_step_length_px(step_px, average_speed, speed_scale),
                )
                candidate_display = (
                    corrected_display
                    + average_direction / average_norm * step_length_px
                )
                candidate_display, clipped = _clip_display_step_to_viewport(
                    corrected_display,
                    candidate_display,
                    viewport,
                )
                candidate = _candidate_data_from_display_step(
                    current_data=current_data,
                    current_display=current_display,
                    candidate_display=candidate_display,
                    jacobian=0.5 * (current_jacobian + next_jacobian),
                    transform=transform,
                )
                if candidate is None or not _point_within_grid_data(grid, candidate):
                    break

        actual_step = np.hypot(*(candidate_display - current_display))
        if actual_step <= 0.2:
            break

        points.append(candidate)
        previous_display = current_display
        current_display = candidate_display
        current_data = candidate
        travelled += actual_step

        if clipped:
            break

    if len(points) < 2:
        return None
    return np.asarray(points)


def _ncl_step_length_px(base_step_px, local_speed, speed_scale):
    speed_scale = max(float(speed_scale), 1e-12)
    speed_fraction = np.clip(float(local_speed) / speed_scale, 0.0, 1.0)
    return max(0.35, float(base_step_px) * speed_fraction * speed_fraction)


def _corrected_ncl_display_origin(current_display, previous_display):
    if previous_display is None:
        return np.asarray(current_display, dtype=float)
    current_display = np.asarray(current_display, dtype=float)
    previous_display = np.asarray(previous_display, dtype=float)
    return current_display - (current_display - previous_display) / 3.0


def _clip_display_step_to_viewport(start_display, end_display, viewport):
    start_display = np.asarray(start_display, dtype=float)
    end_display = np.asarray(end_display, dtype=float)
    if np.all(np.isfinite(end_display)) and viewport.contains(*end_display):
        return end_display, False

    delta = end_display - start_display
    if not np.all(np.isfinite(delta)):
        return start_display, True

    factors = [1.0]
    if delta[0] < -1e-12:
        factors.append((viewport.x0 - start_display[0]) / delta[0])
    elif delta[0] > 1e-12:
        factors.append((viewport.x1 - start_display[0]) / delta[0])

    if delta[1] < -1e-12:
        factors.append((viewport.y0 - start_display[1]) / delta[1])
    elif delta[1] > 1e-12:
        factors.append((viewport.y1 - start_display[1]) / delta[1])

    factor = float(np.clip(min(factors), 0.0, 1.0))
    return start_display + delta * factor, True


def _curve_to_display(curve, transform, display_sampler=None):
    curve = np.asarray(curve, dtype=float)
    if curve.ndim != 2 or curve.shape[1] != 2:
        return None

    if display_sampler is not None:
        try:
            display_curve = display_sampler.sample_display_points(curve)
        except Exception:
            display_curve = None
        if display_curve is not None:
            display_curve = np.asarray(display_curve, dtype=float)
            if np.all(np.isfinite(display_curve)):
                return display_curve

    try:
        display_curve = transform.transform(curve)
    except Exception:
        return None
    if not np.all(np.isfinite(display_curve)):
        return None
    return np.asarray(display_curve, dtype=float)


def _display_points_to_data(transform, display_points, inverse_transform=None):
    display_points = np.asarray(display_points, dtype=float)
    if display_points.ndim == 1:
        display_points = display_points[np.newaxis, :]
    if display_points.ndim != 2 or display_points.shape[1] != 2:
        return None
    if not np.all(np.isfinite(display_points)):
        return None

    try:
        inverse = (
            transform.inverted() if inverse_transform is None else inverse_transform
        )
        data_points = inverse.transform(display_points)
    except Exception:
        return None
    if not np.all(np.isfinite(data_points)):
        return None
    return np.asarray(data_points, dtype=float)


def _display_to_data(transform, display_point, inverse_transform=None):
    data_points = _display_points_to_data(
        transform,
        np.asarray([display_point], dtype=float),
        inverse_transform=inverse_transform,
    )
    if data_points is None:
        return None
    return data_points[0]


def _candidate_data_from_display_step(
    current_data, current_display, candidate_display, jacobian, transform
):
    display_step = np.asarray(candidate_display, dtype=float) - np.asarray(
        current_display,
        dtype=float,
    )
    data_step = _display_step_to_data(jacobian, display_step)
    if data_step is not None:
        candidate = np.asarray(current_data, dtype=float) + data_step
        if np.all(np.isfinite(candidate)):
            return candidate
    return _display_to_data(transform, candidate_display)


def _postprocess_ncl_curve(curve, transform, display_sampler=None):
    display_curve = _curve_to_display(curve, transform, display_sampler=display_sampler)
    if display_curve is None:
        return None

    smoothed_display = _chaikin_smooth_display_curve(display_curve, refinements=1)
    smoothed_display = _fit_single_bend_display_curve(smoothed_display)
    try:
        return transform.inverted().transform(smoothed_display)
    except Exception:
        return curve


def _chaikin_smooth_display_curve(display_curve, refinements=1):
    curve = np.asarray(display_curve, dtype=float)
    if len(curve) < 3:
        return curve

    for _ in range(refinements):
        refined = [curve[0]]
        for left, right in zip(curve[:-1], curve[1:]):
            refined.append(0.75 * left + 0.25 * right)
            refined.append(0.25 * left + 0.75 * right)
        refined.append(curve[-1])
        curve = np.asarray(refined, dtype=float)
    return curve


def _fit_single_bend_display_curve(display_curve, samples=11, max_offset_frac=0.22):
    curve = np.asarray(display_curve, dtype=float)
    if len(curve) < 3:
        return curve

    start = curve[0]
    end = curve[-1]
    chord = end - start
    chord_length = np.hypot(*chord)
    if chord_length <= 1e-6:
        return curve

    unit = chord / chord_length
    normal = np.array([-unit[1], unit[0]])
    midpoint = _point_at_arc_fraction(curve, 0.5)
    chord_midpoint = 0.5 * (start + end)
    midpoint_offset = midpoint - chord_midpoint

    perpendicular_offset = float(np.dot(midpoint_offset, normal))
    max_offset = chord_length * max_offset_frac
    perpendicular_offset = float(np.clip(perpendicular_offset, -max_offset, max_offset))

    target_midpoint = chord_midpoint + normal * perpendicular_offset
    control = 2.0 * target_midpoint - 0.5 * (start + end)

    t = np.linspace(0.0, 1.0, max(int(samples), 3))
    one_minus_t = 1.0 - t
    return (
        (one_minus_t**2)[:, None] * start
        + (2.0 * one_minus_t * t)[:, None] * control
        + (t**2)[:, None] * end
    )


def _point_at_arc_fraction(curve, fraction):
    curve = np.asarray(curve, dtype=float)
    if len(curve) == 0:
        raise ValueError("curve must contain at least one point")
    if len(curve) == 1:
        return curve[0]

    segment_vectors = np.diff(curve, axis=0)
    segment_lengths = np.hypot(segment_vectors[:, 0], segment_vectors[:, 1])
    cumulative = np.concatenate([[0.0], np.cumsum(segment_lengths)])
    total_length = cumulative[-1]
    if total_length <= 1e-12:
        return curve[0]

    target = float(np.clip(fraction, 0.0, 1.0)) * total_length
    idx = np.searchsorted(cumulative, target, side="right") - 1
    idx = int(np.clip(idx, 0, len(segment_lengths) - 1))

    seg_length = segment_lengths[idx]
    if seg_length <= 1e-12:
        return curve[idx]

    local_frac = (target - cumulative[idx]) / seg_length
    return curve[idx] + local_frac * segment_vectors[idx]


def _evaluate_ncl_display_curve(curve, transform, viewport=None):
    try:
        display_curve = transform.transform(np.asarray(curve, dtype=float))
    except Exception:
        return None, True

    if not np.all(np.isfinite(display_curve)) or len(display_curve) < 2:
        return None, False

    segments = np.diff(display_curve, axis=0)
    seg_lengths = np.hypot(segments[:, 0], segments[:, 1])
    valid = seg_lengths > 1e-6
    if np.count_nonzero(valid) < 1:
        return None, False

    if viewport is not None:
        viewport_diag = float(np.hypot(viewport.width, viewport.height))
        if np.isfinite(viewport_diag) and viewport_diag > 0.0:
            jump_limit = max(viewport_diag * 0.35, 1e-6)
            if np.any(seg_lengths > jump_limit):
                return None, False

    segments = segments[valid]
    seg_lengths = seg_lengths[valid]
    arc_length = float(np.sum(seg_lengths))
    chord_length = float(np.hypot(*(display_curve[-1] - display_curve[0])))
    if arc_length <= 1e-6 or chord_length <= 1e-6:
        return None, False

    directions = segments / seg_lengths[:, None]
    if len(directions) < 2:
        return display_curve, False

    turn_cos = np.clip(np.sum(directions[:-1] * directions[1:], axis=1), -1.0, 1.0)
    turn_angles = np.degrees(np.arccos(turn_cos))
    max_turn = float(np.max(turn_angles, initial=0.0))
    total_turn = float(np.sum(turn_angles))
    tortuosity = arc_length / chord_length

    cross = (
        directions[:-1, 0] * directions[1:, 1] - directions[:-1, 1] * directions[1:, 0]
    )
    turn_sign = np.sign(cross[np.abs(cross) > 1e-6])
    sign_changes = (
        int(np.count_nonzero(turn_sign[1:] * turn_sign[:-1] < 0))
        if len(turn_sign) > 1
        else 0
    )

    if tortuosity > 1.28:
        return None, False
    if max_turn > 52.0:
        return None, False
    if total_turn > 120.0:
        return None, False
    if sign_changes > 0 and total_turn > 70.0:
        return None, False
    return display_curve, False


def _acceptable_ncl_display_curve(curve, transform):
    display_curve, _ = _evaluate_ncl_display_curve(curve, transform)
    return display_curve


def _curve_shape_is_acceptable(curve, transform, display_sampler=None):
    display_curve, transform_failed = _evaluate_ncl_display_curve(curve, transform)
    if display_curve is not None:
        return True
    return transform_failed


def _sample_local_vector_state(grid, u, v, transform, point, display_sampler=None):
    point = np.asarray(point, dtype=float)
    u_value = _sample_grid_field(grid, u, point[0], point[1])
    v_value = _sample_grid_field(grid, v, point[0], point[1])
    if u_value is None or v_value is None:
        return None

    if display_sampler is not None:
        sampled_mapping = display_sampler.sample(point)
        if sampled_mapping is None:
            return None
        origin_display, jacobian = sampled_mapping
    else:
        jacobian = _local_display_jacobian(transform, point, grid)
        if jacobian is None:
            return None

        origin_display = transform.transform(np.asarray([point]))[0]
        if not np.all(np.isfinite(origin_display)):
            return None

    display_vector = jacobian @ np.array([u_value, v_value], dtype=float)
    display_norm = np.hypot(*display_vector)
    if not np.isfinite(display_norm) or display_norm <= 1e-12:
        return None

    return (
        origin_display,
        jacobian,
        display_vector / display_norm,
        np.hypot(u_value, v_value),
    )


def _sample_grid_field(grid, field, xd, yd):
    if _sample_grid_field_native is not None and not np.ma.isMaskedArray(field):
        try:
            value = _sample_grid_field_native(
                field=np.asarray(field, dtype=float),
                x_origin=float(grid.x_origin),
                y_origin=float(grid.y_origin),
                dx=float(grid.dx),
                dy=float(grid.dy),
                x=float(xd),
                y=float(yd),
            )
        except Exception as err:
            _disable_native_helper("_sample_grid_field_native", "scalar sampling", err)
            value = None
        if value is not None:
            value = float(value)
            return value if np.isfinite(value) else None

    xi = (xd - grid.x_origin) * grid.inv_dx
    yi = (yd - grid.y_origin) * grid.inv_dy
    if not grid.within_grid(xi, yi):
        return None
    if np.ma.isMaskedArray(field):
        try:
            value = interpgrid(field, xi, yi)
        except TerminateTrajectory:
            return None
        if np.ma.is_masked(value) or not np.isfinite(value):
            return None
        return float(value)

    x = int(xi)
    y = int(yi)
    xn = x if x == grid.nx - 1 else x + 1
    yn = y if y == grid.ny - 1 else y + 1

    xt = float(xi - x)
    yt = float(yi - y)
    a00 = field[y, x]
    a01 = field[y, xn]
    a10 = field[yn, x]
    a11 = field[yn, xn]
    value = (a00 * (1.0 - xt) + a01 * xt) * (1.0 - yt) + (
        a10 * (1.0 - xt) + a11 * xt
    ) * yt
    if not np.isfinite(value):
        return None
    return float(value)


def _sample_grid_field_array(grid, field, points):
    points = np.asarray(points, dtype=float)
    if points.ndim == 1:
        points = points[np.newaxis, :]

    sampled = np.full(len(points), np.nan, dtype=float)
    if len(points) == 0:
        return sampled

    if _sample_grid_field_array_native is not None and not np.ma.isMaskedArray(field):
        try:
            sampled_native = _sample_grid_field_array_native(
                field=np.asarray(field, dtype=float),
                x_origin=float(grid.x_origin),
                y_origin=float(grid.y_origin),
                dx=float(grid.dx),
                dy=float(grid.dy),
                points=points,
            )
        except Exception as err:
            _disable_native_helper(
                "_sample_grid_field_array_native",
                "vectorized sampling",
                err,
            )
            sampled_native = None
        if sampled_native is not None:
            sampled_native = np.asarray(sampled_native, dtype=float)
            if sampled_native.shape == sampled.shape:
                sampled_native[~np.isfinite(sampled_native)] = np.nan
                return sampled_native

    xi = (points[:, 0] - grid.x_origin) * grid.inv_dx
    yi = (points[:, 1] - grid.y_origin) * grid.inv_dy
    valid = (0.0 <= xi) & (xi <= grid.nx - 1) & (0.0 <= yi) & (yi <= grid.ny - 1)
    if not np.any(valid):
        return sampled

    values = interpgrid(field, xi[valid], yi[valid])
    values = np.asarray(np.ma.filled(values, np.nan), dtype=float)
    sampled[valid] = values
    sampled[~np.isfinite(sampled)] = np.nan
    return sampled


def _local_display_jacobian(transform, point, grid):
    point = np.asarray(point, dtype=float)
    try:
        origin_display = transform.transform(np.asarray([point]))[0]
    except Exception:
        return None

    x_step = _finite_difference_step(
        point[0], grid.x_origin, grid.x_origin + grid.width, grid.dx / 2.0
    )
    y_step = _finite_difference_step(
        point[1], grid.y_origin, grid.y_origin + grid.height, grid.dy / 2.0
    )
    if abs(x_step) <= 1e-12 or abs(y_step) <= 1e-12:
        return None

    try:
        x_display = transform.transform(np.asarray([[point[0] + x_step, point[1]]]))[0]
        y_display = transform.transform(np.asarray([[point[0], point[1] + y_step]]))[0]
    except Exception:
        return None

    if not (np.all(np.isfinite(x_display)) and np.all(np.isfinite(y_display))):
        return None

    return np.column_stack(
        (
            (x_display - origin_display) / x_step,
            (y_display - origin_display) / y_step,
        )
    )


def _display_step_to_data(jacobian, display_step):
    jacobian = np.asarray(jacobian, dtype=float)
    display_step = np.asarray(display_step, dtype=float)
    determinant = jacobian[0, 0] * jacobian[1, 1] - jacobian[0, 1] * jacobian[1, 0]
    if not np.isfinite(determinant) or abs(determinant) <= 1e-12:
        return None
    inverse = (
        np.array(
            [
                [jacobian[1, 1], -jacobian[0, 1]],
                [-jacobian[1, 0], jacobian[0, 0]],
            ],
            dtype=float,
        )
        / determinant
    )
    data_step = inverse @ display_step
    if not np.all(np.isfinite(data_step)):
        return None
    return data_step


def _finite_difference_step(value, lower, upper, base_step):
    base_step = abs(float(base_step))
    if base_step <= 1e-12:
        return 0.0
    if value + base_step <= upper:
        return base_step
    if value - base_step >= lower:
        return -base_step
    span = upper - lower
    return span / 4.0 if span > 0 else 0.0


def _point_within_grid_data(grid, point):
    return (
        grid.x_origin <= point[0] <= grid.x_origin + grid.width
        and grid.y_origin <= point[1] <= grid.y_origin + grid.height
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
    geometry = _tip_display_geometry(
        curve,
        transform,
        head_length_px * 1.25,
        display_curve=display_curve,
        display_sampler=display_sampler,
    )
    if geometry is None:
        return None

    tip_display, unit = geometry
    normal = np.array([-unit[1], unit[0]])
    base_center = tip_display - unit * head_length_px
    display_vertices = np.vstack(
        [
            tip_display,
            base_center + normal * head_width_px / 2.0,
            base_center - normal * head_width_px / 2.0,
        ]
    )

    data_vertices = _display_points_to_data(
        transform,
        display_vertices,
        inverse_transform=inverse_transform,
    )
    if data_vertices is None:
        return None

    return patches.Polygon(
        data_vertices,
        closed=True,
        transform=transform,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=max(float(linewidth) * 0.5, 0.5),
        alpha=alpha,
        zorder=zorder,
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
    if _uses_open_arrow_head(arrowstyle):
        return (
            _build_open_arrow_segments(
                curve=curve,
                grid=grid,
                transform=transform,
                head_length_px=head_length_px,
                head_width_px=head_width_px,
                display_curve=display_curve,
                display_sampler=display_sampler,
                inverse_transform=inverse_transform,
            ),
            None,
        )

    return (
        [],
        _build_arrow_polygon(
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
        ),
    )


def _uses_open_arrow_head(arrowstyle):
    return str(arrowstyle).strip() == "->"


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


def _trim_display_curve_from_end(curve, distance):
    curve = np.asarray(curve, dtype=float)
    if len(curve) < 2:
        return curve

    segment_vectors = np.diff(curve, axis=0)
    segment_lengths = np.hypot(segment_vectors[:, 0], segment_vectors[:, 1])
    total_length = float(np.sum(segment_lengths))
    if total_length <= 1e-12:
        return curve

    distance = float(np.clip(distance, 0.0, total_length * 0.55))
    if distance <= 1e-6:
        return curve

    remaining = distance
    for idx in range(len(curve) - 1, 0, -1):
        right = curve[idx]
        left = curve[idx - 1]
        segment = right - left
        segment_length = float(np.hypot(*segment))
        if segment_length <= 1e-12:
            continue
        if remaining < segment_length:
            new_end = right - segment * (remaining / segment_length)
            return np.vstack([curve[:idx], new_end])
        remaining -= segment_length

    midpoint = 0.5 * (curve[0] + curve[1])
    return np.vstack([curve[0], midpoint])


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


def _tip_display_geometry_from_display_curve(display_curve, backoff_px):
    tip_display = display_curve[-1]
    tail_display = _point_at_arc_distance_from_end(display_curve, backoff_px)
    direction = tip_display - tail_display
    direction_norm = np.hypot(*direction)
    if direction_norm <= 1e-12:
        return None
    return tip_display, direction / direction_norm


def _point_at_arc_distance_from_end(curve, distance):
    curve = np.asarray(curve, dtype=float)
    if len(curve) == 0:
        raise ValueError("curve must contain at least one point")
    if len(curve) == 1:
        return curve[0]

    target = max(float(distance), 0.0)
    remaining = target
    for idx in range(len(curve) - 1, 0, -1):
        right = curve[idx]
        left = curve[idx - 1]
        segment = right - left
        segment_length = np.hypot(*segment)
        if segment_length <= 1e-12:
            continue
        if remaining <= segment_length:
            return right - segment * (remaining / segment_length)
        remaining -= segment_length
    return curve[0]


# Coordinate definitions
# ========================


class DomainMap:
    """
    Map representing different coordinate systems.

    Coordinate definitions:

    * axes-coordinates goes from 0 to 1 in the domain.
    * data-coordinates are specified by the input x-y coordinates.
    * grid-coordinates goes from 0 to N and 0 to M for an N x M grid,
      where N and M match the shape of the input data.
    * mask-coordinates goes from 0 to N and 0 to M for an N x M mask,
      where N and M are user-specified to control the density of streamlines.

    This class also has methods for adding trajectories to the StreamMask.
    Before adding a trajectory, run `start_trajectory` to keep track of regions
    crossed by a given trajectory. Later, if you decide the trajectory is bad
    (e.g., if the trajectory is very short) just call `undo_trajectory`.
    """

    def __init__(self, grid, mask):
        self.grid = grid
        self.mask = mask
        # Constants for conversion between grid- and mask-coordinates
        self.x_grid2mask = (mask.nx - 1) / (grid.nx - 1)
        self.y_grid2mask = (mask.ny - 1) / (grid.ny - 1)

        self.x_mask2grid = 1.0 / self.x_grid2mask
        self.y_mask2grid = 1.0 / self.y_grid2mask

        self.x_data2grid = 1.0 / grid.dx
        self.y_data2grid = 1.0 / grid.dy

    def grid2mask(self, xi, yi):
        """Return nearest space in mask-coords from given grid-coords."""
        return round(xi * self.x_grid2mask), round(yi * self.y_grid2mask)

    def mask2grid(self, xm, ym):
        return xm * self.x_mask2grid, ym * self.y_mask2grid

    def data2grid(self, xd, yd):
        return xd * self.x_data2grid, yd * self.y_data2grid

    def grid2data(self, xg, yg):
        return xg / self.x_data2grid, yg / self.y_data2grid

    def start_trajectory(self, xg, yg, broken_streamlines=True):
        xm, ym = self.grid2mask(xg, yg)
        self.mask._start_trajectory(xm, ym, broken_streamlines)

    def reset_start_point(self, xg, yg):
        xm, ym = self.grid2mask(xg, yg)
        self.mask._current_xy = (xm, ym)

    def update_trajectory(self, xg, yg, broken_streamlines=True):
        if not self.grid.within_grid(xg, yg):
            raise InvalidIndexError
        xm, ym = self.grid2mask(xg, yg)
        self.mask._update_trajectory(xm, ym, broken_streamlines)

    def undo_trajectory(self):
        self.mask._undo_trajectory()


class Grid:
    """Grid of data."""

    def __init__(self, x, y, allow_non_uniform=False):

        if np.ndim(x) == 1:
            pass
        elif np.ndim(x) == 2:
            x_row = x[0]
            if not np.allclose(x_row, x):
                raise ValueError("The rows of 'x' must be equal")
            x = x_row
        else:
            raise ValueError("'x' can have at maximum 2 dimensions")

        if np.ndim(y) == 1:
            pass
        elif np.ndim(y) == 2:
            yt = np.transpose(y)  # Also works for nested lists.
            y_col = yt[0]
            if not np.allclose(y_col, yt):
                raise ValueError("The columns of 'y' must be equal")
            y = y_col
        else:
            raise ValueError("'y' can have at maximum 2 dimensions")

        if not (np.diff(x) > 0).all():
            raise ValueError("'x' must be strictly increasing")
        if not (np.diff(y) > 0).all():
            raise ValueError("'y' must be strictly increasing")

        self.nx = len(x)
        self.ny = len(y)

        if self.nx < 2 or self.ny < 2:
            raise ValueError("'x' and 'y' must each contain at least 2 points")

        self.dx = x[1] - x[0]
        self.dy = y[1] - y[0]

        self.x_origin = x[0]
        self.y_origin = y[0]

        self.width = x[-1] - x[0]
        self.height = y[-1] - y[0]

        # Only check for equal spacing if not allowing non-uniform grids
        if not allow_non_uniform:
            if not np.allclose(np.diff(x), self.width / (self.nx - 1)):
                raise ValueError("'x' values must be equally spaced")
            if not np.allclose(np.diff(y), self.height / (self.ny - 1)):
                raise ValueError("'y' values must be equally spaced")
        else:
            # For non-uniform grids, use average spacing
            self.dx = self.width / (self.nx - 1)
            self.dy = self.height / (self.ny - 1)

        self.inv_dx = 1.0 / max(float(self.dx), 1e-12)
        self.inv_dy = 1.0 / max(float(self.dy), 1e-12)

    @property
    def shape(self):
        return self.ny, self.nx

    def within_grid(self, xi, yi):
        """Return whether (*xi*, *yi*) is a valid index of the grid."""
        # Note that xi/yi can be floats; so, for example, we can't simply check
        # `xi < self.nx` since *xi* can be `self.nx - 1 < xi < self.nx`
        return 0 <= xi <= self.nx - 1 and 0 <= yi <= self.ny - 1


class StreamMask:
    """
    Mask to keep track of discrete regions crossed by streamlines.

    The resolution of this grid determines the approximate spacing between
    trajectories. Streamlines are only allowed to pass through zeroed cells:
    When a streamline enters a cell, that cell is set to 1, and no new
    streamlines are allowed to enter.
    """

    def __init__(self, density):
        try:
            self.nx, self.ny = (30 * np.broadcast_to(density, 2)).astype(int)
        except ValueError as err:
            raise ValueError("'density' must be a scalar or be of length " "2") from err
        if self.nx < 0 or self.ny < 0:
            raise ValueError("'density' must be positive")
        self._mask = np.zeros((self.ny, self.nx))
        self.shape = self._mask.shape

        self._current_xy = None

    def __getitem__(self, args):
        return self._mask[args]

    def _start_trajectory(self, xm, ym, broken_streamlines=True):
        """Start recording streamline trajectory"""
        self._traj = []
        self._update_trajectory(xm, ym, broken_streamlines)

    def _undo_trajectory(self):
        """Remove current trajectory from mask"""
        for t in self._traj:
            self._mask[t] = 0

    def _update_trajectory(self, xm, ym, broken_streamlines=True):
        """
        Update current trajectory position in mask.

        If the new position has already been filled, raise `InvalidIndexError`.
        """
        if self._current_xy != (xm, ym):
            if self[ym, xm] == 0:
                self._traj.append((ym, xm))
                self._mask[ym, xm] = 1
                self._current_xy = (xm, ym)
            else:
                if broken_streamlines:
                    raise InvalidIndexError
                else:
                    pass


class InvalidIndexError(Exception):
    pass


class TerminateTrajectory(Exception):
    pass


# Integrator definitions
# =======================


def _get_integrator(u, v, dmap, resolution, magnitude, integration_direction):

    # rescale velocity onto grid-coordinates for integrations.
    u, v = dmap.data2grid(u, v)

    # speed (path length) will be in axes-coordinates
    u_ax = u / (dmap.grid.nx - 1)
    v_ax = v / (dmap.grid.ny - 1)
    speed = np.ma.sqrt(u_ax**2 + v_ax**2)

    def forward_time(xi, yi):
        if not dmap.grid.within_grid(xi, yi):
            raise OutOfBounds
        ds_dt = interpgrid(speed, xi, yi)
        if ds_dt == 0:
            raise TerminateTrajectory()
        dt_ds = 1.0 / ds_dt
        ui = interpgrid(u, xi, yi)
        vi = interpgrid(v, xi, yi)
        return ui * dt_ds, vi * dt_ds

    def backward_time(xi, yi):
        dxi, dyi = forward_time(xi, yi)
        return -dxi, -dyi

    def integrate(x0, y0, broken_streamlines=True):
        """
        Return x, y grid-coordinates of trajectory based on starting point.

        Integrate both forward and backward in time from starting point in
        grid coordinates.

        Integration is terminated when a trajectory reaches a domain boundary
        or when it crosses into an already occupied cell in the StreamMask. The
        resulting trajectory is None if it is shorter than `minlength`.
        """

        stotal, xy_traj = 0.0, []

        try:
            dmap.start_trajectory(x0, y0, broken_streamlines)
        except InvalidIndexError:
            return None
        if integration_direction in ["both", "backward"]:
            s, xyt = _integrate_rk12(
                x0, y0, dmap, backward_time, resolution, magnitude, broken_streamlines
            )
            stotal += s
            xy_traj += xyt[::-1]

        if integration_direction in ["both", "forward"]:
            dmap.reset_start_point(x0, y0)
            s, xyt = _integrate_rk12(
                x0, y0, dmap, forward_time, resolution, magnitude, broken_streamlines
            )
            stotal += s
            xy_traj += xyt[1:]

        if len(xy_traj) > 1:
            return np.broadcast_arrays(xy_traj, np.empty((1, 2)))[0]
        else:  # reject short trajectories
            dmap.undo_trajectory()
            return None

    return integrate


class OutOfBounds(IndexError):
    pass


def _integrate_rk12(x0, y0, dmap, f, resolution, magnitude, broken_streamlines=True):
    """
    2nd-order Runge-Kutta algorithm with adaptive step size.

    This method is also referred to as the improved Euler's method, or Heun's
    method. This method is favored over higher-order methods because:

    1. To get decent looking trajectories and to sample every mask cell
       on the trajectory we need a small timestep, so a lower order
       solver doesn't hurt us unless the data is *very* high resolution.
       In fact, for cases where the user inputs
       data smaller or of similar grid size to the mask grid, the higher
       order corrections are negligible because of the very fast linear
       interpolation used in `interpgrid`.

    2. For high resolution input data (i.e. beyond the mask
       resolution), we must reduce the timestep. Therefore, an adaptive
       timestep is more suited to the problem as this would be very hard
       to judge automatically otherwise.

    This integrator is about 1.5 - 2x as fast as RK4 and RK45 solvers (using
    similar Python implementations) in most setups.
    """
    # This error is below that needed to match the RK4 integrator. It
    # is set for visual reasons -- too low and corners start
    # appearing ugly and jagged. Can be tuned.
    maxerror = 0.003

    # This limit is important (for all integrators) to avoid the
    # trajectory skipping some mask cells. We could relax this
    # condition if we use the code which is commented out below to
    # increment the location gradually. However, due to the efficient
    # nature of the interpolation, this doesn't boost speed by much
    # for quite a bit of complexity.
    maxds = min(1.0 / dmap.mask.nx, 1.0 / dmap.mask.ny, 0.1)

    ds = maxds
    stotal = 0
    xi = x0
    yi = y0
    xyf_traj = []
    m_total = []

    while True:
        try:
            if dmap.grid.within_grid(xi, yi):
                xyf_traj.append((xi, yi))
                m_total.append(interpgrid(magnitude, xi, yi))
                maxlength = resolution * np.mean(m_total)
            else:
                raise OutOfBounds

            # Compute the two intermediate gradients.
            # f should raise OutOfBounds if the locations given are
            # outside the grid.
            k1x, k1y = f(xi, yi)
            k2x, k2y = f(xi + ds * k1x, yi + ds * k1y)

        except OutOfBounds:
            # Out of the domain during this step.
            # Take an Euler step to the boundary to improve neatness
            # unless the trajectory is currently empty.
            if xyf_traj:
                ds, xyf_traj = _euler_step(xyf_traj, dmap, f)
                stotal += ds
            break
        except TerminateTrajectory:
            break

        dx1 = ds * k1x
        dy1 = ds * k1y
        dx2 = ds * 0.5 * (k1x + k2x)
        dy2 = ds * 0.5 * (k1y + k2y)

        ny, nx = dmap.grid.shape
        # Error is normalized to the axes coordinates
        error = np.hypot((dx2 - dx1) / (nx - 1), (dy2 - dy1) / (ny - 1))

        # Only save step if within error tolerance
        if error < maxerror:
            xi += dx2
            yi += dy2
            try:
                dmap.update_trajectory(xi, yi, broken_streamlines)
            except InvalidIndexError:
                break
            if stotal + ds > maxlength:
                break
            stotal += ds

        # recalculate stepsize based on step error
        if error == 0:
            ds = maxds
        else:
            ds = min(maxds, 0.85 * ds * (maxerror / error) ** 0.5)

    return stotal, xyf_traj


def _euler_step(xyf_traj, dmap, f):
    """Simple Euler integration step that extends streamline to boundary."""
    ny, nx = dmap.grid.shape
    xi, yi = xyf_traj[-1]
    cx, cy = f(xi, yi)
    if cx == 0:
        dsx = np.inf
    elif cx < 0:
        dsx = xi / -cx
    else:
        dsx = (nx - 1 - xi) / cx
    if cy == 0:
        dsy = np.inf
    elif cy < 0:
        dsy = yi / -cy
    else:
        dsy = (ny - 1 - yi) / cy
    ds = min(dsx, dsy)
    xyf_traj.append((xi + cx * ds, yi + cy * ds))
    return ds, xyf_traj


# Utility functions
# ========================


def interpgrid(a, xi, yi):
    """Fast 2D, linear interpolation on an integer grid"""

    Ny, Nx = np.shape(a)
    if isinstance(xi, np.ndarray):
        x = np.clip(np.asarray(xi, dtype=int), 0, Nx - 1)
        y = np.clip(np.asarray(yi, dtype=int), 0, Ny - 1)
        # Check that xn, yn don't exceed max index
        xn = np.clip(x + 1, 0, Nx - 1)
        yn = np.clip(y + 1, 0, Ny - 1)
    else:
        x = int(np.clip(int(xi), 0, Nx - 1))
        y = int(np.clip(int(yi), 0, Ny - 1))
        # conditional is faster than clipping for integers
        if x == (Nx - 1):
            xn = x
        else:
            xn = x + 1
        if y == (Ny - 1):
            yn = y
        else:
            yn = y + 1

    a00 = a[y, x]
    a01 = a[y, xn]
    a10 = a[yn, x]
    a11 = a[yn, xn]
    xt = xi - x
    yt = yi - y
    a0 = a00 * (1 - xt) + a01 * xt
    a1 = a10 * (1 - xt) + a11 * xt
    ai = a0 * (1 - yt) + a1 * yt

    if not isinstance(xi, np.ndarray):
        if np.ma.is_masked(ai):
            raise TerminateTrajectory

    return ai


def _gen_starting_points(x, y, grains):
    if isinstance(grains, tuple):
        nx, ny = grains
    elif isinstance(grains, int):
        nx = ny = grains

    eps = np.finfo(np.float32).eps

    tmp_x = np.linspace(x.min() + eps, x.max() - eps, nx)
    tmp_y = np.linspace(y.min() + eps, y.max() - eps, ny)

    xs = np.tile(tmp_x, ny)
    ys = np.repeat(tmp_y, nx)

    seed_points = np.array([xs, ys])

    return seed_points.T
