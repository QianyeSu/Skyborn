"""Internal grid and regridding helpers for Skyborn curly-vector plotting."""

from __future__ import annotations

import numpy as np

from ._ncl_vector_dataset import _filled_scalar_field_array


def _default_cartopy_target_extent(ax, target_extent):
    if target_extent is not None or not hasattr(ax, "get_extent"):
        return target_extent
    try:
        return ax.get_extent(ax.projection)
    except Exception:
        return target_extent


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


def _rcm2rgrid_fields(lat2d, lon2d, fields, lat1d, lon1d):
    try:
        from skyborn.interp import rcm2rgrid
    except ImportError as err:
        raise ImportError(
            "Curvilinear vector support requires skyborn.interp.rcm2rgrid to be available."
        ) from err

    field_stack = np.asarray(
        [np.asarray(field, dtype=float) for field in fields],
        dtype=float,
    )
    regridded = rcm2rgrid(
        np.asarray(lat2d, dtype=float),
        np.asarray(lon2d, dtype=float),
        field_stack,
        np.asarray(lat1d, dtype=float),
        np.asarray(lon1d, dtype=float),
        msg=np.nan,
        meta=False,
    )
    regridded = np.asarray(regridded, dtype=float)
    if regridded.ndim == 2:
        regridded = regridded[np.newaxis, :, :]
    return [regridded[idx] for idx in range(regridded.shape[0])]


def _maybe_as_scalar_field(value, expected_shape):
    if value is None or isinstance(value, str) or np.isscalar(value):
        return None
    array = np.asarray(value)
    if array.shape != expected_shape:
        return None
    return _filled_scalar_field_array(array)


def _regrid_curvilinear_vectors(x, y, u, v, *scalars, target_shape):
    lon1d, lat1d = _build_curvilinear_target_grid(x, y, target_shape)
    regridded_fields = _rcm2rgrid_fields(y, x, (u, v, *scalars), lat1d, lon1d)
    u_reg, v_reg, *scalar_regs = regridded_fields
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


def _has_cyclic_longitude_endpoint(x, period=360.0):
    x = np.asarray(x, dtype=float)
    if len(x) < 2:
        return False

    spacing = float(np.nanmedian(np.diff(x)))
    tolerance = max(1e-6, abs(spacing) * 1e-3)
    span = float(x[-1] - x[0])
    return np.isfinite(span) and abs(span - float(period)) <= tolerance


def _append_cyclic_column(x, field, period=360.0):
    x = np.asarray(x, dtype=float)
    field = np.asarray(field, dtype=float)
    x_cyclic = np.concatenate([x, [x[0] + float(period)]])
    field_cyclic = np.concatenate([field, field[:, :1]], axis=1)
    return x_cyclic, field_cyclic


def _wrap_periodic_grid_queries(x_query, x_1d, period=360.0):
    x_query = np.asarray(x_query, dtype=float)
    x_1d = np.asarray(x_1d, dtype=float)
    if len(x_1d) < 2:
        return x_query

    has_cyclic_endpoint = _has_cyclic_longitude_endpoint(x_1d, period=period)
    if not has_cyclic_endpoint and not _grid_spans_full_longitude(x_1d):
        return x_query

    x_min = float(x_1d[0])
    x_max = float(x_1d[-1])
    wrapped = ((x_query - x_min) % float(period)) + x_min
    if not has_cyclic_endpoint:
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
