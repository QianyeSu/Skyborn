"""Display-space thinning helpers for Skyborn vector plotting."""

from __future__ import annotations

import numpy as np


def _display_jump_threshold(lengths):
    lengths = np.asarray(lengths, dtype=float)
    finite = lengths[np.isfinite(lengths) & (lengths > 1e-9)]
    if finite.size == 0:
        return None

    sorted_lengths = np.sort(finite)
    lower_half = sorted_lengths[: max(sorted_lengths.size // 2, 1)]
    baseline = float(np.nanmedian(lower_half))
    return max(baseline * 12.0, 1e-6)


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
        finite = np.isfinite(self.display_grid).all(axis=2)
        self.cell_valid = self._compute_cell_valid(finite)

    def _compute_cell_valid(self, finite):
        cell_valid = (
            finite[:-1, :-1] & finite[:-1, 1:] & finite[1:, :-1] & finite[1:, 1:]
        )
        if not np.any(cell_valid):
            return cell_valid

        horizontal_edges = np.hypot(
            np.diff(self.display_grid[:, :, 0], axis=1),
            np.diff(self.display_grid[:, :, 1], axis=1),
        )
        vertical_edges = np.hypot(
            np.diff(self.display_grid[:, :, 0], axis=0),
            np.diff(self.display_grid[:, :, 1], axis=0),
        )
        jump_limit = _display_jump_threshold(
            np.concatenate([horizontal_edges.ravel(), vertical_edges.ravel()])
        )
        if jump_limit is None:
            return cell_valid

        p00 = self.display_grid[:-1, :-1]
        p01 = self.display_grid[:-1, 1:]
        p10 = self.display_grid[1:, :-1]
        p11 = self.display_grid[1:, 1:]

        top_edges = np.hypot(p01[:, :, 0] - p00[:, :, 0], p01[:, :, 1] - p00[:, :, 1])
        bottom_edges = np.hypot(
            p11[:, :, 0] - p10[:, :, 0],
            p11[:, :, 1] - p10[:, :, 1],
        )
        left_edges = np.hypot(p10[:, :, 0] - p00[:, :, 0], p10[:, :, 1] - p00[:, :, 1])
        right_edges = np.hypot(
            p11[:, :, 0] - p01[:, :, 0],
            p11[:, :, 1] - p01[:, :, 1],
        )
        max_edge = np.maximum.reduce([top_edges, bottom_edges, left_edges, right_edges])
        return cell_valid & (max_edge <= jump_limit)

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

        count = len(points)
        display_points = np.full((count, 2), np.nan, dtype=float)
        jacobians = (
            np.full((count, 2, 2), np.nan, dtype=float) if include_jacobian else None
        )
        valid = np.zeros(count, dtype=bool)

        if count == 0 or self.grid.nx < 2 or self.grid.ny < 2:
            return display_points, jacobians, valid

        xi = (points[:, 0] - self.x_origin) * self.inv_dx
        yi = (points[:, 1] - self.y_origin) * self.inv_dy
        within = (
            (0.0 <= xi)
            & (xi <= self.max_x_index)
            & (0.0 <= yi)
            & (yi <= self.max_y_index)
        )
        if not np.any(within):
            return display_points, jacobians, valid

        within_indices = np.nonzero(within)[0]
        xi = np.clip(xi[within], 0.0, self.max_x_index)
        yi = np.clip(yi[within], 0.0, self.max_y_index)
        ix = np.floor(np.clip(xi, 0.0, self.max_cell_x)).astype(int)
        iy = np.floor(np.clip(yi, 0.0, self.max_cell_y)).astype(int)
        cell_valid = self.cell_valid[iy, ix]
        if not np.any(cell_valid):
            return display_points, jacobians, valid

        keep = within_indices[cell_valid]
        xi = xi[cell_valid]
        yi = yi[cell_valid]
        ix = ix[cell_valid]
        iy = iy[cell_valid]
        sx = (xi - ix)[:, None]
        sy = (yi - iy)[:, None]

        p00 = self.display_grid[iy, ix]
        p01 = self.display_grid[iy, ix + 1]
        p10 = self.display_grid[iy + 1, ix]
        p11 = self.display_grid[iy + 1, ix + 1]

        one_minus_sx = 1.0 - sx
        one_minus_sy = 1.0 - sy
        display_points[keep] = (
            p00 * one_minus_sx * one_minus_sy
            + p01 * sx * one_minus_sy
            + p10 * one_minus_sx * sy
            + p11 * sx * sy
        )
        valid[keep] = np.isfinite(display_points[keep]).all(axis=1)

        if not include_jacobian:
            return display_points, jacobians, valid

        dfdx = ((p01 - p00) * one_minus_sy + (p11 - p10) * sy) / self.dx_safe
        dfdy = ((p10 - p00) * one_minus_sx + (p11 - p01) * sx) / self.dy_safe
        jacobians[keep, :, 0] = dfdx
        jacobians[keep, :, 1] = dfdy

        determinants = (
            jacobians[keep, 0, 0] * jacobians[keep, 1, 1]
            - jacobians[keep, 0, 1] * jacobians[keep, 1, 0]
        )
        jacobian_valid = np.isfinite(jacobians[keep]).all(axis=(1, 2)) & (
            np.abs(determinants) > 1e-12
        )
        valid[keep] &= jacobian_valid
        display_points[keep[~jacobian_valid]] = np.nan
        jacobians[keep[~jacobian_valid]] = np.nan

        return display_points, jacobians, valid


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
