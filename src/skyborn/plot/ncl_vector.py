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

from .vector_key import CurlyVectorKey, curly_vector_key
from .vector_plot import CurlyVectorPlotSet
from .vector_plot import curly_vector as _array_curly_vector

if TYPE_CHECKING:
    from matplotlib.axes import Axes

__all__ = ["curly_vector", "CurlyVectorKey", "curly_vector_key"]


def _is_cartopy_crs_like(transform) -> bool:
    return transform is not None and hasattr(transform, "_as_mpl_transform")


def _looks_like_axes(value: Any) -> bool:
    return hasattr(value, "add_collection") and hasattr(value, "transData")


def _apply_dataset_isel(da: xr.DataArray, isel):
    if not isel:
        return da
    indexers = {dim: value for dim, value in dict(isel).items() if dim in da.dims}
    return da.isel(indexers) if indexers else da


def _get_plot_dataarray(ds: xr.Dataset, value, *, isel=None, role: str) -> xr.DataArray:
    if isinstance(value, xr.DataArray):
        da = value
    else:
        da = ds[value]
    da = _apply_dataset_isel(da, isel)
    da = da.squeeze(drop=True)
    if da.ndim == 0:
        raise ValueError(f"{role} must resolve to at least one dimension")
    if da.ndim > 2:
        remaining_dims = [dim for dim, size in da.sizes.items() if size > 1]
        raise ValueError(
            f"{role} must resolve to 1D or 2D data after selection; remaining dims: "
            f"{remaining_dims}. Pass a 2D slice or use isel=..."
        )
    return da


def _extract_curly_vector_dataset_source(ds, x, y, u, v, *, isel=None):
    x_da = _get_plot_dataarray(ds, x, isel=isel, role="x")
    y_da = _get_plot_dataarray(ds, y, isel=isel, role="y")
    u_da = _get_plot_dataarray(ds, u, isel=isel, role="u")
    v_da = _get_plot_dataarray(ds, v, isel=isel, role="v")

    if u_da.ndim != 2 or v_da.ndim != 2:
        raise ValueError("u and v must each resolve to 2D arrays before plotting")
    if u_da.shape != v_da.shape:
        raise ValueError(
            "u and v must share the same 2D shape. If your vector components live on "
            "different staggered or unmatched grids, align them onto the same physical "
            "grid before calling curly_vector()."
        )

    if x_da.ndim == 1 and y_da.ndim == 1:
        x_values = np.asarray(x_da.data, dtype=float)
        y_values = np.asarray(y_da.data, dtype=float)
        u_values = np.asarray(u_da.data, dtype=float)
        v_values = np.asarray(v_da.data, dtype=float)

        expected_shape = (y_da.size, x_da.size)
        if u_da.shape != expected_shape:
            raise ValueError(
                f"u/v shape {u_da.shape} does not match the rectilinear x/y grid "
                f"shape {expected_shape}"
            )

        if x_values.size > 1 and x_values[0] > x_values[-1]:
            x_values = x_values[::-1]
            u_values = u_values[:, ::-1]
            v_values = v_values[:, ::-1]
        if y_values.size > 1 and y_values[0] > y_values[-1]:
            y_values = y_values[::-1]
            u_values = u_values[::-1, :]
            v_values = v_values[::-1, :]

        return x_values, y_values, u_values, v_values
    elif x_da.ndim == 2 and y_da.ndim == 2:
        if x_da.shape != y_da.shape:
            raise ValueError("2D x and y coordinates must have the same shape")
        if x_da.shape != u_da.shape:
            raise ValueError(
                f"2D x/y shape {x_da.shape} must match the u/v shape {u_da.shape}"
            )
    else:
        raise ValueError("x and y must both be 1D or both be 2D")

    return (
        np.asarray(x_da.data, dtype=float),
        np.asarray(y_da.data, dtype=float),
        np.asarray(u_da.data, dtype=float),
        np.asarray(v_da.data, dtype=float),
    )


def _default_cartopy_target_extent(ax, target_extent):
    if target_extent is not None or not hasattr(ax, "get_extent"):
        return target_extent
    try:
        return ax.get_extent(ax.projection)
    except Exception:
        return target_extent


def _regrid_cartopy_vectors(
    src_crs,
    target_proj,
    regrid_shape,
    x,
    y,
    u,
    v,
    *scalars,
    target_extent=None,
):
    try:
        from scipy.interpolate import RegularGridInterpolator
    except ImportError as err:
        raise ImportError(
            "scipy is required for Cartopy projection regridding support."
        ) from err

    x_1d, y_1d, u, v, scalar_list = _prepare_source_vector_grid(
        x, y, u, v, scalars=scalars
    )
    x_target, y_target = _build_projection_target_grid(
        target_proj, regrid_shape, target_extent
    )
    source_points = src_crs.transform_points(target_proj, x_target, y_target)
    x_source = source_points[..., 0]
    y_source = source_points[..., 1]
    valid = np.isfinite(x_source) & np.isfinite(y_source)

    x_query = _wrap_periodic_grid_queries(x_source, x_1d)
    sample_points = np.column_stack([y_source[valid], x_query[valid]])

    u_interp = RegularGridInterpolator(
        (y_1d, x_1d),
        u,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    v_interp = RegularGridInterpolator(
        (y_1d, x_1d),
        v,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )

    u_sampled = np.full(x_target.shape, np.nan, dtype=float)
    v_sampled = np.full(y_target.shape, np.nan, dtype=float)
    if np.any(valid):
        u_sampled[valid] = u_interp(sample_points)
        v_sampled[valid] = v_interp(sample_points)

    valid_vectors = valid & np.isfinite(u_sampled) & np.isfinite(v_sampled)
    u_target = np.full_like(u_sampled, np.nan)
    v_target = np.full_like(v_sampled, np.nan)
    if np.any(valid_vectors):
        ut, vt = target_proj.transform_vectors(
            src_crs,
            x_source[valid_vectors],
            y_source[valid_vectors],
            u_sampled[valid_vectors],
            v_sampled[valid_vectors],
        )
        u_target[valid_vectors] = ut
        v_target[valid_vectors] = vt

    scalar_targets = []
    for scalar in scalar_list:
        scalar_interp = RegularGridInterpolator(
            (y_1d, x_1d),
            scalar,
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        scalar_target = np.full(x_target.shape, np.nan, dtype=float)
        if np.any(valid):
            scalar_target[valid] = scalar_interp(sample_points)
        scalar_targets.append(scalar_target)

    return (x_target, y_target, u_target, v_target, *scalar_targets)


def _normalize_regrid_shape(regrid_shape):
    if np.isscalar(regrid_shape):
        size = max(int(regrid_shape), 2)
        return size, size
    nx, ny = regrid_shape
    return max(int(nx), 2), max(int(ny), 2)


def _normalize_density_pair(density):
    try:
        density_pair = np.broadcast_to(np.asarray(density, dtype=float), 2)
    except ValueError as err:
        raise ValueError("density must be a scalar or a length-2 sequence") from err
    return float(density_pair[0]), float(density_pair[1])


def _is_curvilinear_grid(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 2 or y.ndim != 2:
        return False
    if x.shape != y.shape:
        raise ValueError("2D x and y coordinates must have matching shapes")

    x_mesh = np.broadcast_to(x[0, :], x.shape)
    y_mesh = np.broadcast_to(y[:, 0][:, None], y.shape)
    return not (
        np.allclose(x, x_mesh, equal_nan=True)
        and np.allclose(y, y_mesh, equal_nan=True)
    )


def _default_curvilinear_regrid_shape(x, density):
    density_x, density_y = _normalize_density_pair(density)
    ny_src, nx_src = np.asarray(x).shape
    nx = min(nx_src, max(24, int(np.clip(round(72.0 * density_x), 24, 144))))
    ny = min(ny_src, max(18, int(np.clip(round(56.0 * density_y), 18, 112))))
    return nx, ny


def _build_curvilinear_target_grid(x, y, target_shape):
    nx, ny = _normalize_regrid_shape(target_shape)
    lon2d = np.asarray(x, dtype=float)
    lat2d = np.asarray(y, dtype=float)
    finite = np.isfinite(lon2d) & np.isfinite(lat2d)
    if not np.any(finite):
        raise ValueError("Curvilinear x/y coordinates contain no finite points")

    lon_work = np.rad2deg(np.unwrap(np.deg2rad(lon2d), axis=1))
    lon1d = np.linspace(
        float(np.nanmin(lon_work[finite])), float(np.nanmax(lon_work[finite])), nx
    )
    lat1d = np.linspace(
        float(np.nanmin(lat2d[finite])), float(np.nanmax(lat2d[finite])), ny
    )
    return lon1d, lat1d


def _rcm2rgrid_2d(lat2d, lon2d, field, lat1d, lon1d):
    try:
        from skyborn.interp import rcm2rgrid
    except ImportError as err:
        raise ImportError(
            "Curvilinear vector support requires skyborn.interp.rcm2rgrid to be available."
        ) from err

    field = np.asarray(field, dtype=float)
    regridded = rcm2rgrid(
        np.asarray(lat2d, dtype=float),
        np.asarray(lon2d, dtype=float),
        field[np.newaxis, :, :],
        np.asarray(lat1d, dtype=float),
        np.asarray(lon1d, dtype=float),
        msg=np.nan,
        meta=False,
    )
    regridded = np.asarray(regridded, dtype=float)
    if regridded.ndim == 3:
        regridded = regridded[0]
    return regridded


def _maybe_as_scalar_field(value, expected_shape):
    if value is None or isinstance(value, str) or np.isscalar(value):
        return None
    array = np.asarray(value)
    if array.shape != expected_shape:
        return None
    return np.asarray(array, dtype=float)


def _regrid_curvilinear_vectors(x, y, u, v, *scalars, target_shape):
    lon1d, lat1d = _build_curvilinear_target_grid(x, y, target_shape)
    u_reg = _rcm2rgrid_2d(y, x, u, lat1d, lon1d)
    v_reg = _rcm2rgrid_2d(y, x, v, lat1d, lon1d)
    scalar_regs = [_rcm2rgrid_2d(y, x, scalar, lat1d, lon1d) for scalar in scalars]
    return (lon1d, lat1d, u_reg, v_reg, *scalar_regs)


def _build_projection_target_grid(target_proj, regrid_shape, target_extent):
    if target_extent is None:
        target_extent = target_proj.x_limits + target_proj.y_limits
    nx, ny = _normalize_regrid_shape(regrid_shape)
    x_target = np.linspace(float(target_extent[0]), float(target_extent[1]), nx)
    y_target = np.linspace(float(target_extent[2]), float(target_extent[3]), ny)
    return np.meshgrid(x_target, y_target)


def _prepare_source_vector_grid(x, y, u, v, scalars=()):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    scalar_list = [np.asarray(s, dtype=float) for s in scalars]

    if x.ndim == 2:
        x = x[0, :]
    if y.ndim == 2:
        y = y[:, 0]

    if x.ndim != 1 or y.ndim != 1:
        raise ValueError(
            "Cartopy projection regridding requires 1D or meshgrid x/y coordinates"
        )

    if u.shape != (len(y), len(x)) or v.shape != (len(y), len(x)):
        raise ValueError(
            "u and v must match the source y/x grid shape for projection regridding"
        )

    if len(x) > 1 and x[0] > x[-1]:
        x = x[::-1]
        u = u[:, ::-1]
        v = v[:, ::-1]
        scalar_list = [s[:, ::-1] for s in scalar_list]

    if len(y) > 1 and y[0] > y[-1]:
        y = y[::-1]
        u = u[::-1, :]
        v = v[::-1, :]
        scalar_list = [s[::-1, :] for s in scalar_list]

    if _grid_spans_full_longitude(x):
        x_base = x.copy()
        x, u = _append_cyclic_column(x_base, u)
        _, v = _append_cyclic_column(x_base, v)
        scalar_list = [_append_cyclic_column(x_base, s)[1] for s in scalar_list]

    return x, y, u, v, scalar_list


def _grid_spans_full_longitude(x):
    if len(x) < 2:
        return False
    spacing = float(np.nanmedian(np.diff(x)))
    span = float(x[-1] - x[0] + spacing)
    return np.isfinite(span) and 300.0 <= span <= 370.0


def _append_cyclic_column(x, field, period=360.0):
    x = np.asarray(x, dtype=float)
    field = np.asarray(field, dtype=float)
    x_cyclic = np.concatenate([x, [x[0] + float(period)]])
    field_cyclic = np.concatenate([field, field[:, :1]], axis=1)
    return x_cyclic, field_cyclic


def _wrap_periodic_grid_queries(x_query, x_1d, period=360.0):
    x_query = np.asarray(x_query, dtype=float)
    x_1d = np.asarray(x_1d, dtype=float)
    if len(x_1d) < 2 or not _grid_spans_full_longitude(x_1d):
        return x_query

    x_min = float(x_1d[0])
    x_max = float(x_1d[-1])
    wrapped = ((x_query - x_min) % float(period)) + x_min
    wrapped = np.where(wrapped > x_max, wrapped - float(period), wrapped)
    return wrapped


def _extract_regular_grid_from_regridded_vectors(
    x_grid,
    y_grid,
    u_grid,
    v_grid,
    scalar_grids,
):
    x_grid = np.asarray(x_grid)
    y_grid = np.asarray(y_grid)
    u_grid = np.asarray(u_grid)
    v_grid = np.asarray(v_grid)

    x_1d = np.asarray(x_grid[0, :], dtype=float)
    y_1d = np.asarray(y_grid[:, 0], dtype=float)
    scalar_grids = [np.asarray(s) for s in scalar_grids]

    if x_1d.size > 1 and x_1d[0] > x_1d[-1]:
        x_1d = x_1d[::-1]
        u_grid = u_grid[:, ::-1]
        v_grid = v_grid[:, ::-1]
        scalar_grids = [s[:, ::-1] for s in scalar_grids]

    if y_1d.size > 1 and y_1d[0] > y_1d[-1]:
        y_1d = y_1d[::-1]
        u_grid = u_grid[::-1, :]
        v_grid = v_grid[::-1, :]
        scalar_grids = [s[::-1, :] for s in scalar_grids]

    return x_1d, y_1d, u_grid, v_grid, scalar_grids


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

    if _is_curvilinear_grid(x, y):
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
    if isinstance(color, np.ndarray):
        scalar_fields.append(color)
        scalar_keys.append("color")
    if isinstance(linewidth, np.ndarray):
        scalar_fields.append(linewidth)
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
    color: Any = None,
    cmap: Any = None,
    norm: Any = None,
    arrowsize=1,
    arrowstyle="-|>",
    transform: Any = None,
    zorder: float | None = None,
    start_points: Any = None,
    integration_direction="both",
    grains=15,
    broken_streamlines=True,
    anchor: str | None = None,
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
        Axes on which to plot. By default, use the current axes. Mutually exclusive with `size` and `figsize`.
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
    color : color or 2D array
        The streamline color. If given an array, its values are converted to
        colors using *cmap* and *norm*.  The array must have the same shape
        as *u* and *v*.
    cmap, norm
        Data normalization and colormapping parameters for *color*; only used
        if *color* is an array of floats. See `~.Axes.imshow` for a detailed
        description.
    arrowsize : float
        Scaling factor for the arrow size.
    arrowstyle : str
        Arrow style specification.
        See `~matplotlib.patches.FancyArrowPatch`.
        ``'->'`` uses an open line arrowhead while other styles fall back to a
        filled polygon head.
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
        dataset yourself.

    Returns
    -------
    CurlyVectorPlotSet
        Container object with attributes

        - ``lines``: `.LineCollection` of streamlines

        - ``arrows``: `.PatchCollection` containing `.FancyArrowPatch`
          objects representing the arrows half-way along streamlines.

            This container will probably change in the future to allow changes
            to the colormap, alpha, etc. for both lines and arrows, but these
            changes should be backward compatible.

    .. seealso::
        - https://github.com/matplotlib/matplotlib/issues/20038
        - https://github.com/kieranmrhunt/curved-quivers
        - https://github.com/Deltares/dfm_tools/issues/483
        - https://github.com/NCAR/geocat-viz/issues/4
        - https://docs.xarray.dev/en/stable/generated/xarray.Dataset.plot.streamplot.html#xarray.Dataset.plot.streamplot
    """
    x, y, u, v = _extract_curly_vector_dataset_source(ds, x, y, u, v, isel=isel)

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
        density=density,
        linewidth=linewidth,
        color=color,
        cmap=cmap,
        norm=norm,
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
        ncl_preset=ncl_preset,
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
    color: Any = None,
    cmap: Any = None,
    norm: Any = None,
    arrowsize=1,
    arrowstyle="-|>",
    transform: Any = None,
    zorder: float | None = None,
    start_points: Any = None,
    integration_direction="both",
    grains=15,
    broken_streamlines=True,
    anchor: str | None = None,
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
        density=density,
        linewidth=linewidth,
        color=color,
        cmap=cmap,
        norm=norm,
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
        ncl_preset=ncl_preset,
    )


def curly_vector(*args: Any, **kwargs: Any) -> CurlyVectorPlotSet:
    """Plot NCL-like curly vectors from arrays or an xarray dataset.

    Supported call styles
    ---------------------
    Array/Matplotlib style
        ``curly_vector(ax, x, y, u, v, ...)``
        ``curly_vector(x, y, u, v, ..., ax=ax)``
        ``curly_vector(x, y, u, v, ...)``

    xarray dataset style
        ``curly_vector(ds, x='lon', y='lat', u='u', v='v', ax=ax, ...)``
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
