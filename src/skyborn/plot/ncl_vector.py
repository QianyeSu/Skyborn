"""Dataset and Cartopy adapters for Skyborn curly-vector plots.

Author: Qianye Su <suqianye2000@gmail.com>
Copyright (c) 2025-2026 Qianye Su
Created: 2026-03-01 14:58:56
"""

from __future__ import annotations

from collections.abc import Hashable
from typing import TYPE_CHECKING, Any

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from . import _ncl_vector_dataset as _dataset_helpers
from . import _ncl_vector_regrid as _regrid_helpers
from .vector_key import CurlyVectorKey, curly_vector_key
from .vector_plot import CurlyVectorPlotSet, _resolve_curly_style_aliases
from .vector_plot import curly_vector as _array_curly_vector

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = ["curly_vector", "CurlyVectorKey", "curly_vector_key"]
_ARRAY_CURLY_VECTOR_KWARG_NAMES = (
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
    "pivot",
    "ref_magnitude",
    "ref_length",
    "min_frac_length",
    "min_distance",
    "ncl_preset",
)


def _is_cartopy_crs_like(transform) -> bool:
    return transform is not None and hasattr(transform, "_as_mpl_transform")


def _looks_like_axes(value: Any) -> bool:
    return hasattr(value, "add_collection") and hasattr(value, "transData")


def _collect_named_kwargs(
    scope: dict[str, Any], names: tuple[str, ...]
) -> dict[str, Any]:
    """Collect a stable subset of keyword arguments from a local scope."""
    return {name: scope[name] for name in names}


_ISSUED_PLOT_WARNINGS = _dataset_helpers._ISSUED_PLOT_WARNINGS
_warn_plot_once = _dataset_helpers._warn_plot_once
_apply_dataset_isel = _dataset_helpers._apply_dataset_isel
_get_plot_dataarray = _dataset_helpers._get_plot_dataarray
_transpose_2d_dataarray_to_dims = _dataset_helpers._transpose_2d_dataarray_to_dims
_filled_scalar_field_array = _dataset_helpers._filled_scalar_field_array
_extract_curly_vector_dataset_source = (
    _dataset_helpers._extract_curly_vector_dataset_source
)
_prepare_dataset_style_field = _dataset_helpers._prepare_dataset_style_field


_default_cartopy_target_extent = _regrid_helpers._default_cartopy_target_extent
_normalize_regrid_shape = _regrid_helpers._normalize_regrid_shape
_normalize_density_pair = _regrid_helpers._normalize_density_pair
_is_curvilinear_grid = _regrid_helpers._is_curvilinear_grid
_default_curvilinear_regrid_shape = _regrid_helpers._default_curvilinear_regrid_shape
_build_curvilinear_target_grid = _regrid_helpers._build_curvilinear_target_grid
_rcm2rgrid_2d = _regrid_helpers._rcm2rgrid_2d
_rcm2rgrid_fields = _regrid_helpers._rcm2rgrid_fields
_maybe_as_scalar_field = _regrid_helpers._maybe_as_scalar_field
_regrid_curvilinear_vectors = _regrid_helpers._regrid_curvilinear_vectors
_build_projection_target_grid = _regrid_helpers._build_projection_target_grid
_prepare_source_vector_grid = _regrid_helpers._prepare_source_vector_grid
_grid_spans_full_longitude = _regrid_helpers._grid_spans_full_longitude
_has_cyclic_longitude_endpoint = _regrid_helpers._has_cyclic_longitude_endpoint
_append_cyclic_column = _regrid_helpers._append_cyclic_column
_wrap_periodic_grid_queries = _regrid_helpers._wrap_periodic_grid_queries
_extract_regular_grid_from_regridded_vectors = (
    _regrid_helpers._extract_regular_grid_from_regridded_vectors
)
_regrid_cartopy_vectors = _regrid_helpers._regrid_cartopy_vectors


def _prepare_curly_vector_dataset_inputs(
    ax,
    x,
    y,
    u,
    v,
    transform,
    regrid_shape,
    curvilinear_regrid_shape,
    target_extent,
    color,
    linewidth,
    density,
):
    src_crs = transform if _is_cartopy_crs_like(transform) else None
    color_field = _maybe_as_scalar_field(color, np.shape(u))
    linewidth_field = _maybe_as_scalar_field(linewidth, np.shape(u))
    if color_field is not None:
        color = color_field
    if linewidth_field is not None:
        linewidth = linewidth_field

    if _is_curvilinear_grid(x, y):
        scalar_fields = []
        scalar_keys = []
        if color_field is not None:
            scalar_fields.append(color_field)
            scalar_keys.append("color")
        if linewidth_field is not None:
            scalar_fields.append(linewidth_field)
            scalar_keys.append("linewidth")

        target_shape = (
            _normalize_regrid_shape(curvilinear_regrid_shape)
            if curvilinear_regrid_shape is not None
            else _default_curvilinear_regrid_shape(x, density)
        )
        regrid_result = _regrid_curvilinear_vectors(
            x,
            y,
            u,
            v,
            *scalar_fields,
            target_shape=target_shape,
        )
        x, y, u, v, *scalar_grids = regrid_result
        for key, scalar_grid in zip(scalar_keys, scalar_grids):
            if key == "color":
                color = scalar_grid
            elif key == "linewidth":
                linewidth = scalar_grid

    if regrid_shape is None:
        if src_crs is not None:
            transform = src_crs._as_mpl_transform(ax)
        return x, y, u, v, color, linewidth, transform

    if src_crs is None:
        raise ValueError("regrid_shape requires a Cartopy CRS passed via transform")
    if not hasattr(ax, "projection"):
        raise ValueError("regrid_shape requires a Cartopy GeoAxes with a projection")

    scalar_fields = []
    scalar_keys = []
    color_field = _maybe_as_scalar_field(color, np.shape(u))
    linewidth_field = _maybe_as_scalar_field(linewidth, np.shape(u))
    if color_field is not None:
        scalar_fields.append(color_field)
        scalar_keys.append("color")
    if linewidth_field is not None:
        scalar_fields.append(linewidth_field)
        scalar_keys.append("linewidth")

    regrid_result = _regrid_cartopy_vectors(
        src_crs=src_crs,
        target_proj=ax.projection,
        regrid_shape=regrid_shape,
        x=x,
        y=y,
        u=u,
        v=v,
        *scalar_fields,
        target_extent=_default_cartopy_target_extent(ax, target_extent),
    )

    x_grid, y_grid, u_grid, v_grid, *scalar_grids = regrid_result
    x, y, u, v, scalar_grids = _extract_regular_grid_from_regridded_vectors(
        x_grid,
        y_grid,
        u_grid,
        v_grid,
        scalar_grids,
    )

    for key, scalar_grid in zip(scalar_keys, scalar_grids):
        if key == "color":
            color = scalar_grid
        elif key == "linewidth":
            linewidth = scalar_grid

    return x, y, u, v, color, linewidth, None


def _curly_vector_from_dataset(
    ds: xr.Dataset,
    x: Hashable,
    y: Hashable,
    u: Hashable,
    v: Hashable,
    ax: Axes | None = None,
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
    arrowsize=1,
    arrowstyle="->",
    transform: Any = None,
    zorder: float | None = None,
    start_points: Any = None,
    integration_direction="both",
    grains=15,
    broken_streamlines=True,
    anchor: str | None = None,
    pivot: str | None = None,
    ref_magnitude: float | None = None,
    ref_length: float | None = None,
    min_frac_length=0.0,
    min_distance: float | None = None,
    ncl_preset: str | None = None,
    regrid_shape: Any = None,
    curvilinear_regrid_shape: Any = None,
    target_extent: Any = None,
    isel: Any = None,
) -> CurlyVectorPlotSet:
    """
    Plot NCL-like curved vector glyphs for a 2D vector flow.

    .. warning::

        This function is experimental and the API is subject to change. Please use with caution.

    Parameters
    ----------
    ds : :py:class:`xarray.Dataset`.
        Wind dataset.
    x : Hashable or None, optional.
        Variable name for x-axis.
    y : Hashable or None, optional.
        Variable name for y-axis.
    u : Hashable or None, optional.
        Variable name for the u velocity (in `x` direction).
    v : Hashable or None, optional.
        Variable name for the v velocity (in `y` direction).
    ax : :py:class:`matplotlib.axes.Axes`, optional.
        Axes on which to plot. By default, use the current axes.
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
        Matplotlib artist alpha applied to both shafts and arrow heads.
    facecolor, edgecolor : color-like, optional
        Explicit arrow-head fill and edge colors. These mainly affect the
        filled ``arrowstyle="-|>"`` head. When omitted, the resolved shaft
        color is reused.
    facecolors, edgecolors : color-like, optional
        Matplotlib-style aliases for ``facecolor`` and ``edgecolor``.
    rasterized : bool, optional
        Whether to rasterize the generated curly-vector artists when exporting
        to vector formats such as PDF or SVG.
    arrowsize : float
        Scaling factor for the arrow size.
    arrowstyle : str
        Supported arrow-head style. Use ``"->"`` for the open NCL-like line
        head or ``"-|>"`` for a filled triangular head.
    start_points : (N, 2) array
        Coordinates of starting points for the streamlines in data coordinates
        (the same coordinates as the *x* and *y* arrays).
    zorder : float
        The zorder of the streamlines and arrows.
        Artists with lower zorder values are drawn first.
    integration_direction : {'forward', 'backward', 'both'}, default: 'both'
        Integrate the streamline in forward, backward or both directions.
    broken_streamlines : boolean, default: True
        If False, forces streamlines to continue until they
        leave the plot domain.  If True, they may be terminated if they
        come too close to another streamline.
    anchor : {'tail', 'center', 'head'} or None, default: None
        Anchor point used by the NCL-like curved-glyph renderer. If omitted, the value
        is inferred from ``integration_direction``.
    pivot : {'tail', 'mid', 'middle', 'tip'} or None, default: None
        Matplotlib ``quiver``-style alias for ``anchor``.
    ref_magnitude : float or None, default: None
        Reference vector magnitude for mapping data magnitude to glyph length.
    ref_length : float or None, default: None
        Reference glyph length as a fraction of the axes width.
    min_frac_length : float, default: 0.0
        Minimum glyph length as a fraction of ``ref_length``.
    min_distance : float or None, default: None
        Minimum center-to-center spacing as a fraction of the axes width.
    ncl_preset : {None, 'profile'}, default: None
        Optional preset for NCL-like glyph tuning. ``'profile'`` applies a
        conservative lat-pressure-style reference scaling while preserving the
        default lat-lon behavior when omitted.
    regrid_shape : int or (int, int), optional
        If provided together with a Cartopy CRS in ``transform``, first regrid
        the vector field onto a regular grid in the target map projection using
        Cartopy's projection-aware vector interpolation. This is especially
        helpful for polar projections where native lat-lon vectors can become
        visually over-dense near the pole.
    curvilinear_regrid_shape : int or (int, int), optional
        If ``x`` and ``y`` are 2D curvilinear coordinates, first regularize the
        vector field onto a rectilinear lon/lat grid of this shape using
        ``skyborn.interp.rcm2rgrid``. If omitted, a conservative shape is
        derived from ``density`` to keep rendering cost under control.
    target_extent : tuple, optional
        Optional target-projection extent forwarded to Cartopy regridding as
        ``(x0, x1, y0, y1)``. If omitted and the axes can report an extent,
        the current axes extent is used.
    isel : mapping, optional
        Optional positional indexers applied to ``x``, ``y``, ``u``, and ``v``
        before plotting. This is useful when the dataset stores extra
        dimensions and you want to select a 2D slice without pre-slicing the
        dataset yourself. Indexer keys that do not match the selected
        DataArray dimensions are ignored with a warning.

    Returns
    -------
    CurlyVectorPlotSet
        Container object with attributes

        - ``lines``: `.LineCollection` of streamlines

        - ``arrows``: tuple of the actual filled arrow-head patches added to
          the axes. Open arrow styles use line segments only and therefore
          return an empty tuple.

    .. seealso::
        - https://github.com/matplotlib/matplotlib/issues/20038
        - https://github.com/kieranmrhunt/curved-quivers
        - https://github.com/Deltares/dfm_tools/issues/483
        - https://github.com/NCAR/geocat-viz/issues/4
        - https://docs.xarray.dev/en/stable/generated/xarray.Dataset.plot.streamplot.html#xarray.Dataset.plot.streamplot
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

    x, y, u, v, source_metadata = _extract_curly_vector_dataset_source(
        ds,
        x,
        y,
        u,
        v,
        isel=isel,
        return_metadata=True,
    )
    color = _prepare_dataset_style_field(
        color,
        isel=isel,
        expected_shape=np.shape(u),
        vector_dims=source_metadata["vector_dims"],
        x_descending=bool(source_metadata["x_descending"]),
        y_descending=bool(source_metadata["y_descending"]),
        role="color",
    )
    linewidth = _prepare_dataset_style_field(
        linewidth,
        isel=isel,
        expected_shape=np.shape(u),
        vector_dims=source_metadata["vector_dims"],
        x_descending=bool(source_metadata["x_descending"]),
        y_descending=bool(source_metadata["y_descending"]),
        role="linewidth",
    )

    if ax is None:
        ax = plt.gca()
    x, y, u, v, color, linewidth, transform = _prepare_curly_vector_dataset_inputs(
        ax=ax,
        x=x,
        y=y,
        u=u,
        v=v,
        transform=transform,
        regrid_shape=regrid_shape,
        curvilinear_regrid_shape=curvilinear_regrid_shape,
        target_extent=target_extent,
        color=color,
        linewidth=linewidth,
        density=density,
    )

    obj = _array_curly_vector(
        ax,
        x,
        y,
        u,
        v,
        **_collect_named_kwargs(locals(), _ARRAY_CURLY_VECTOR_KWARG_NAMES),
    )
    return obj


def _curly_vector_from_arrays(
    ax: Axes,
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
    arrowsize=1,
    arrowstyle="->",
    transform: Any = None,
    zorder: float | None = None,
    start_points: Any = None,
    integration_direction="both",
    grains=15,
    broken_streamlines=True,
    anchor: str | None = None,
    pivot: str | None = None,
    ref_magnitude: float | None = None,
    ref_length: float | None = None,
    min_frac_length=0.0,
    min_distance: float | None = None,
    ncl_preset: str | None = None,
    regrid_shape: Any = None,
    curvilinear_regrid_shape: Any = None,
    target_extent: Any = None,
) -> CurlyVectorPlotSet:
    """Array-input adapter that preserves Cartopy and curvilinear support."""
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

    x, y, u, v, color, linewidth, transform = _prepare_curly_vector_dataset_inputs(
        ax=ax,
        x=np.asarray(x),
        y=np.asarray(y),
        u=np.asarray(u),
        v=np.asarray(v),
        transform=transform,
        regrid_shape=regrid_shape,
        curvilinear_regrid_shape=curvilinear_regrid_shape,
        target_extent=target_extent,
        color=color,
        linewidth=linewidth,
        density=density,
    )

    return _array_curly_vector(
        ax,
        x,
        y,
        u,
        v,
        **_collect_named_kwargs(locals(), _ARRAY_CURLY_VECTOR_KWARG_NAMES),
    )


def curly_vector(*args: Any, **kwargs: Any) -> CurlyVectorPlotSet:
    """Plot NCL-like curly vectors from arrays or an xarray dataset.

    This is the public plotting entry point exposed by ``skyborn.plot``.
    It accepts either raw array inputs or an ``xarray.Dataset`` and routes
    them to the same Matplotlib-compatible curly-vector renderer.

    Supported call styles
    ---------------------
    Array/Matplotlib style
        ``curly_vector(ax, x, y, u, v, ...)``
        ``curly_vector(x, y, u, v, ..., ax=ax)``
        ``curly_vector(x, y, u, v, ...)``

    xarray dataset style
        ``curly_vector(ds, x='lon', y='lat', u='u', v='v', ax=ax, ...)``

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Target axes. If omitted, ``matplotlib.pyplot.gca()`` is used.
    x, y : array-like or hashable
        Coordinate specification. For array-style calls, pass the x/y
        coordinates directly. For dataset-style calls, pass the variable or
        coordinate names to read from ``ds``.
    u, v : array-like or hashable
        Vector components on the plotting grid. For dataset-style calls, pass
        the variable names.
    ds : xarray.Dataset, optional
        Dataset source used by the dataset-style call form.
    density : float or tuple of float, default: 1
        Controls glyph density. Larger values produce more curly vectors.
    linewidth : float or 2D array, optional
        Line width for the vector shafts. A 2D field can be used for varying
        line width on the source grid.
    color : color-like or 2D array, optional
        Shaft color. A 2D scalar field may also be supplied together with
        ``cmap`` and ``norm``.
    cmap, norm : optional
        Matplotlib colormap and normalization used when ``color`` is a scalar
        field.
    alpha : float, optional
        Matplotlib artist alpha applied to both shafts and arrow heads.
    facecolor, edgecolor : color-like, optional
        Arrow-head fill and edge colors, similar to ``plt.quiver``. These are
        mainly relevant for the filled ``arrowstyle="-|>"`` head.
    rasterized : bool, optional
        Whether to rasterize the generated curly-vector artists when exporting
        to vector formats such as PDF or SVG.
    arrowsize : float, default: 1
        Scales the arrow-head size.
    arrowstyle : str, default: ``"->"``
        Supported arrow-head style. Use ``"->"`` for the open NCL-like line
        head or ``"-|>"`` for a filled triangular head.
    transform : optional
        Source coordinate transform. For Cartopy usage, this is typically
        ``ccrs.PlateCarree()``.
    zorder : float, optional
        Matplotlib z-order of the generated artists.
    integration_direction : {"forward", "backward", "both"}, default: ``"both"``
        Controls whether each glyph grows forward, backward, or symmetrically
        about its seed point.
    broken_streamlines : bool, default: True
        Whether glyph tracing may terminate early when the local geometry
        becomes unsuitable.
    anchor : {"tail", "center", "head"}, optional
        Explicit glyph anchor override. If omitted, it is inferred from
        ``integration_direction``.
    pivot : {"tail", "mid", "middle", "tip"}, optional
        Matplotlib ``quiver``-style alias for ``anchor``.
    ref_magnitude : float, optional
        Reference vector magnitude used for NCL-like length scaling.
    ref_length : float, optional
        Reference glyph length as a fraction of axes width.
    min_frac_length : float, default: 0.0
        Minimum glyph length as a fraction of the reference length.
    min_distance : float, optional
        Minimum accepted spacing between glyph centers in axes-fraction units.
    ncl_preset : {None, "profile"}, optional
        Optional NCL-style preset. ``"profile"`` is intended for vertical
        cross-sections such as latitude-pressure plots.
    regrid_shape : int or tuple of int, optional
        For Cartopy GeoAxes, first regrid vectors in the target projection
        before rendering. This is especially useful for polar projections.
    curvilinear_regrid_shape : int or tuple of int, optional
        Optional rectification target shape for 2D curvilinear source grids.
    target_extent : tuple, optional
        Explicit Cartopy target-projection extent used by ``regrid_shape``.
    isel : mapping, optional
        Dataset-style positional indexers applied before plotting. Unused
        indexer keys are ignored with a warning.

    Returns
    -------
    CurlyVectorPlotSet
        Container holding the generated line collection, arrow artists, and
        the metadata needed by ``curly_vector_key``.
    """
    if not args:
        raise TypeError(
            "curly_vector() expects either (ax, x, y, u, v, ...) or "
            "(ds, x='...', y='...', u='...', v='...', ...)"
        )

    first = args[0]
    if _looks_like_axes(first):
        return _curly_vector_from_arrays(*args, **kwargs)

    if isinstance(first, xr.Dataset):
        return _curly_vector_from_dataset(*args, **kwargs)

    if len(args) >= 4:
        ax = kwargs.pop("ax", None)
        if ax is None:
            ax = plt.gca()
        return _curly_vector_from_arrays(ax, *args, **kwargs)

    raise TypeError(
        "Unsupported arguments for curly_vector(). Expected either "
        "(ax, x, y, u, v, ...), (x, y, u, v, ...), or "
        "(ds, x='...', y='...', u='...', v='...', ...)."
    )
