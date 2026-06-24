"""Display-space thinning helpers for Skyborn vector plotting."""

from __future__ import annotations

import numpy as np

from ..nclcurly_core import (
    compute_display_cell_valid as _compute_display_cell_valid_native,
)
from ..nclcurly_core import sample_display_grid as _sample_display_grid_native


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
        self.dx_safe = max(float(grid.dx), 1e-12)
        self.dy_safe = max(float(grid.dy), 1e-12)
        self.display_grid = np.asarray(display_grid, dtype=float)
        self.cell_valid = np.asarray(
            _compute_display_cell_valid_native(display_grid=self.display_grid),
            dtype=bool,
        )

    def sample(self, point):
        display_points, jacobians, valid = self.sample_points(
            np.asarray(point, dtype=float),
            include_jacobian=True,
        )
        if not valid[0]:
            return None
        return display_points[0], jacobians[0]

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
