"""Display-space thinning helpers for Skyborn vector plotting."""

from __future__ import annotations

import numpy as np

from ..nclcurly_native import (
    compute_display_cell_valid as _compute_display_cell_valid_native,
)
from ..nclcurly_native import sample_display_grid as _sample_display_grid_native


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

    if np.isscalar(density):
        density_scalar = max(float(density), 0.1)
    else:
        density_xy = np.broadcast_to(np.asarray(density, dtype=float), 2).astype(float)
        density_scalar = float(np.mean(np.maximum(density_xy, 0.1)))

    spacing_frac = 0.9 / (30.0 * density_scalar)
    if ncl_preset is not None and str(ncl_preset).strip().lower() in {
        "profile",
        "vertical_profile",
        "vertical-profile",
        "lat_pressure",
    }:
        spacing_frac *= 0.6
    return spacing_frac


class _NCLDisplaySampler:
    """Sample a precomputed data->display mapping on the plotting grid."""

    def __init__(self, grid, display_grid):
        self.grid = grid
        self.x_origin = float(grid.x_origin)
        self.y_origin = float(grid.y_origin)
        self.inv_dx = float(grid.inv_dx)
        self.inv_dy = float(grid.inv_dy)
        self.max_x_index = float(grid.nx - 1)
        self.max_y_index = float(grid.ny - 1)
        self.max_cell_x = float(grid.nx - 2)
        self.max_cell_y = float(grid.ny - 2)
        self.dx_safe = max(float(grid.dx), 1e-12)
        self.dy_safe = max(float(grid.dy), 1e-12)
        self.display_grid = np.asarray(display_grid, dtype=float)
        self.cell_valid = np.asarray(
            _compute_display_cell_valid_native(display_grid=self.display_grid),
            dtype=bool,
        )

    def sample(self, point):
        point = np.asarray(point, dtype=float)
        xi = (point[0] - self.x_origin) * self.inv_dx
        yi = (point[1] - self.y_origin) * self.inv_dy
        if not (
            0.0 <= xi <= self.max_x_index
            and 0.0 <= yi <= self.max_y_index
            and self.max_cell_x >= 0.0
            and self.max_cell_y >= 0.0
        ):
            return None

        ix = int(np.floor(np.clip(xi, 0.0, self.max_cell_x)))
        iy = int(np.floor(np.clip(yi, 0.0, self.max_cell_y)))
        if not self.cell_valid[iy, ix]:
            return None

        sx = float(xi - ix)
        sy = float(yi - iy)
        p00 = self.display_grid[iy, ix]
        p01 = self.display_grid[iy, ix + 1]
        p10 = self.display_grid[iy + 1, ix]
        p11 = self.display_grid[iy + 1, ix + 1]

        display_point = (
            p00 * (1.0 - sx) * (1.0 - sy)
            + p01 * sx * (1.0 - sy)
            + p10 * (1.0 - sx) * sy
            + p11 * sx * sy
        )
        dfdx = ((p01 - p00) * (1.0 - sy) + (p11 - p10) * sy) / self.dx_safe
        dfdy = ((p10 - p00) * (1.0 - sx) + (p11 - p01) * sx) / self.dy_safe
        jacobian = np.column_stack((dfdx, dfdy))
        determinant = jacobian[0, 0] * jacobian[1, 1] - jacobian[0, 1] * jacobian[1, 0]
        if (
            not np.all(np.isfinite(display_point))
            or not np.all(np.isfinite(jacobian))
            or abs(determinant) <= 1e-12
        ):
            return None
        return display_point, jacobian

    def sample_display_points(self, points):
        display_points, _, valid = self.sample_points(
            points,
            include_jacobian=False,
        )
        display_points[~valid] = np.nan
        return display_points

    def sample_points(self, points, include_jacobian):
        points = np.asarray(points, dtype=float)
        if points.ndim == 1:
            points = points[np.newaxis, :]
        display_points, jacobians, valid = _sample_display_grid_native(
            display_grid=self.display_grid,
            cell_valid=self.cell_valid,
            x_origin=self.x_origin,
            y_origin=self.y_origin,
            dx=self.dx_safe,
            dy=self.dy_safe,
            points=points,
            include_jacobian=bool(include_jacobian),
        )
        return (
            np.asarray(display_points, dtype=float),
            None if jacobians is None else np.asarray(jacobians, dtype=float),
            np.asarray(valid, dtype=bool),
        )


class _NCLNativeTraceContext:
    """Prepacked native tracing inputs reused across many glyph traces."""

    def __init__(self, grid, u, v, viewport, display_sampler):
        self.u = np.asarray(u, dtype=float)
        self.v = np.asarray(v, dtype=float)
        self.display_grid = display_sampler.display_grid
        self.cell_valid = display_sampler.cell_valid
        self.x_origin = float(grid.x_origin)
        self.y_origin = float(grid.y_origin)
        self.dx = float(grid.dx)
        self.dy = float(grid.dy)
        self.viewport_x0 = float(viewport.x0)
        self.viewport_y0 = float(viewport.y0)
        self.viewport_x1 = float(viewport.x1)
        self.viewport_y1 = float(viewport.y1)


def _prepare_ncl_display_sampler(grid, transform):
    if grid.nx < 2 or grid.ny < 2:
        return None

    xs = grid.x_origin + np.arange(grid.nx, dtype=float) * grid.dx
    ys = grid.y_origin + np.arange(grid.ny, dtype=float) * grid.dy
    node_x, node_y = np.meshgrid(xs, ys, indexing="xy")

    try:
        display_points = transform.transform(
            np.column_stack([node_x.ravel(), node_y.ravel()])
        )
    except Exception:
        return None

    display_grid = np.asarray(display_points, dtype=float).reshape(grid.ny, grid.nx, 2)
    if not np.isfinite(display_grid).any():
        return None
    return _NCLDisplaySampler(grid, display_grid)
