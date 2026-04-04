"""Sampling helpers extracted from ``vector_plot.py``."""

from __future__ import annotations

import numpy as np


def _sample_grid_field_python(
    grid,
    field,
    xd,
    yd,
    interpgrid_fn,
    terminate_trajectory_exc,
):
    xi = (xd - grid.x_origin) * grid.inv_dx
    yi = (yd - grid.y_origin) * grid.inv_dy
    if not grid.within_grid(xi, yi):
        return None
    if np.ma.isMaskedArray(field):
        try:
            value = interpgrid_fn(field, xi, yi)
        except terminate_trajectory_exc:
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


def _sample_grid_field_array_python(grid, field, points, interpgrid_fn):
    points = np.asarray(points, dtype=float)
    if points.ndim == 1:
        points = points[np.newaxis, :]

    sampled = np.full(len(points), np.nan, dtype=float)
    if len(points) == 0:
        return sampled

    xi = (points[:, 0] - grid.x_origin) * grid.inv_dx
    yi = (points[:, 1] - grid.y_origin) * grid.inv_dy
    valid = (0.0 <= xi) & (xi <= grid.nx - 1) & (0.0 <= yi) & (yi <= grid.ny - 1)
    if not np.any(valid):
        return sampled

    values = interpgrid_fn(field, xi[valid], yi[valid])
    values = np.asarray(np.ma.filled(values, np.nan), dtype=float)
    sampled[valid] = values
    sampled[~np.isfinite(sampled)] = np.nan
    return sampled
