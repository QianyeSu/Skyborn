"""Cartopy transform and projection-regridding helpers for Skyborn plots."""

from __future__ import annotations

import numpy as np

from .grid_prepare import (
    _normalize_regrid_shape,
    _prepare_source_vector_grid,
    _wrap_periodic_grid_queries,
)


def _default_cartopy_target_extent(ax, target_extent):
    if target_extent is not None or not hasattr(ax, "get_extent"):
        return target_extent
    try:
        return ax.get_extent(ax.projection)
    except Exception:
        return target_extent


def _build_projection_target_grid(target_proj, regrid_shape, target_extent):
    if target_extent is None:
        target_extent = target_proj.x_limits + target_proj.y_limits
    nx, ny = _normalize_regrid_shape(regrid_shape)
    x_target = np.linspace(float(target_extent[0]), float(target_extent[1]), nx)
    y_target = np.linspace(float(target_extent[2]), float(target_extent[3]), ny)
    return np.meshgrid(x_target, y_target)


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
