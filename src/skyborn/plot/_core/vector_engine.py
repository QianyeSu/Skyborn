"""Pure vector-engine helpers extracted from ``vector_plot.py``."""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.collections as mcollections
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import numpy as np
from matplotlib import cm

from .._shared.coords import _coerce_matching_plot_field
from .._shared.style import _ncl_arrow_edge_size_px, _resolve_open_arrow_size
from .geometry import _point_within_grid_data
from .thinning import (
    _map_ncl_display_points_to_viewport,
    _resolve_ncl_min_distance_fraction,
)


def _finite_plot_field_values(field, field_name):
    """Return finite values from a style field or fail fast with a clear error."""
    finite_values = np.asarray(field, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size == 0:
        raise ValueError(
            f"{field_name} field must contain at least one finite value after masking"
        )
    return finite_values


def _resolve_artist_coordinate_context(axes, transform):
    """Choose the final artist transform for the rendered curly-vector geometry."""
    artist_transform = transform
    artist_inverse_transform = None
    bake_display_geometry = bool(
        hasattr(axes, "projection")
        and transform is not None
        and transform is not axes.transData
    )
    if not bake_display_geometry:
        return artist_transform, artist_inverse_transform, False

    try:
        artist_inverse_transform = axes.transData.inverted()
    except Exception:
        return transform, None, False

    return axes.transData, artist_inverse_transform, True


def _default_ncl_max_length_px(axes_bbox, density):
    density_xy = np.broadcast_to(np.asarray(density, dtype=float), 2).astype(float)
    density_xy = np.maximum(density_xy, 0.1)
    nx = max(30.0 * density_xy[0], 1.0)
    ny = max(30.0 * density_xy[1], 1.0)
    sx = float(axes_bbox.width) / nx
    sy = float(axes_bbox.height) / ny
    return max(np.sqrt((sx * sx + sy * sy) / 2.0), 1.0)


def _resolve_ncl_reference_length_px(
    min_mag,
    max_mag,
    ref_mag,
    requested_ref_length_px,
    min_frac_length,
    default_max_length_px,
):
    min_mag = float(min_mag)
    max_mag = float(max_mag)
    ref_mag = float(ref_mag)
    requested_ref_length_px = float(requested_ref_length_px)
    min_frac_length = float(np.clip(min_frac_length, 0.0, 1.0))
    default_max_length_px = max(float(default_max_length_px), 1.0)

    if requested_ref_length_px > 0.0:
        return requested_ref_length_px
    if max_mag - min_mag <= 1e-12:
        return default_max_length_px
    if ref_mag > 0.0 and min_frac_length > 0.0:
        ratio = np.clip((ref_mag - min_mag) / max(max_mag - min_mag, 1e-12), 0.0, 1.0)
        denominator = 1.0 - min_frac_length + min_frac_length * ratio
        if denominator > 1e-12:
            return default_max_length_px * ratio / denominator
    if ref_mag > 0.0:
        return default_max_length_px * ref_mag / max(max_mag, 1e-12)
    return default_max_length_px


def _resolve_ncl_length_scale(
    min_mag,
    max_mag,
    ref_mag,
    requested_ref_length_px,
    min_frac_length,
    default_max_length_px,
):
    uvmn = float(min_mag)
    uvmx = max(float(max_mag), 1e-12)
    rvrm = float(ref_mag)
    rvrl = float(requested_ref_length_px)
    vfr = float(np.clip(min_frac_length, 0.0, 1.0))
    rdmx = max(float(default_max_length_px), 1.0)
    rdmn = 0.0
    vrl = rdmx
    vfl = rdmn
    iav = 0

    if uvmx - uvmn <= 1e-12:
        if rvrl > 0.0 and rvrm > 0.0:
            vrl = rvrl
            rdmx = vrl * uvmx / max(rvrm, 1e-12)
        elif rvrm > 0.0:
            vrl = rdmx * rvrm / uvmx
        elif rvrl > 0.0:
            rdmx = rvrl
            vrl = rdmx
        vfl = vrl
        rdmn = vfl
    elif rvrm <= 0.0:
        if rvrl > 0.0:
            rdmx = rvrl
        vrl = rdmx
        if vfr > 0.0:
            iav = 1
            vfl = vfr * rdmx
            rdmn = vfl
        else:
            rdmn = rdmx * (uvmn / uvmx)
            vfl = rdmn
    elif rvrm <= uvmn:
        iav = 1
        if rvrl > 0.0:
            vrl = rvrl
            rdmx = vrl * uvmx / max(rvrm, 1e-12)
        elif vfr > 0.0:
            vrl = rdmx * vfr
        else:
            vrl = rdmx * rvrm / uvmx
        rdmn = vrl * uvmn / max(rvrm, 1e-12)
        vfl = rdmn
    elif vfr > 0.0:
        iav = 1
        if rvrl > 0.0:
            vrl = rvrl
            vfl = vfr * vrl
            rdmn = vfl
            rdmx = rdmn + (vrl - rdmn) * (uvmx - uvmn) / max(rvrm - uvmn, 1e-12)
        else:
            ratio = (rvrm - uvmn) / max(uvmx - uvmn, 1e-12)
            denominator = 1.0 - vfr + vfr * ratio
            if denominator > 1e-12:
                vrl = rdmx * ratio / denominator
            vfl = vfr * vrl
            rdmn = vfl
    else:
        if rvrl > 0.0:
            vrl = rvrl
            rdmx = vrl * uvmx / max(rvrm, 1e-12)
            vfl = vrl * uvmn / max(rvrm, 1e-12)
            rdmn = vfl
        else:
            vrl = rdmx * rvrm / uvmx
            vfl = rdmx * uvmn / uvmx
            rdmn = vfl

    return {
        "min_mag": uvmn,
        "max_mag": uvmx,
        "ref_mag": rvrm,
        "ref_length_px": vrl,
        "min_length_px": rdmn,
        "max_length_px": rdmx,
        "adjust_min": bool(iav),
        "min_frac_length": vfr,
    }


def _curve_length_from_magnitude(magnitude_value, length_scale):
    magnitude_value = float(magnitude_value)
    min_mag = length_scale["min_mag"]
    max_mag = length_scale["max_mag"]
    max_length_px = length_scale["max_length_px"]
    min_length_px = length_scale["min_length_px"]

    if max_mag - min_mag <= 1e-12:
        return max_length_px

    if not length_scale["adjust_min"]:
        return np.clip(magnitude_value / max_mag, 0.0, 1.0) * max_length_px

    scale = np.clip(
        (magnitude_value - min_mag) / max(max_mag - min_mag, 1e-12), 0.0, 1.0
    )
    return min_length_px + (max_length_px - min_length_px) * scale


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
            yt = np.transpose(y)
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

        if not allow_non_uniform:
            if not np.allclose(np.diff(x), self.width / (self.nx - 1)):
                raise ValueError("'x' values must be equally spaced")
            if not np.allclose(np.diff(y), self.height / (self.ny - 1)):
                raise ValueError("'y' values must be equally spaced")
        else:
            self.dx = self.width / (self.nx - 1)
            self.dy = self.height / (self.ny - 1)

        self.inv_dx = 1.0 / max(float(self.dx), 1e-12)
        self.inv_dy = 1.0 / max(float(self.dy), 1e-12)

    @property
    def shape(self):
        return self.ny, self.nx

    def within_grid(self, xi, yi):
        return 0 <= xi <= self.nx - 1 and 0 <= yi <= self.ny - 1


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


def _sample_local_vector_state(
    grid,
    u,
    v,
    transform,
    point,
    display_sampler=None,
    sample_grid_field_fn=None,
    local_display_jacobian_fn=None,
):
    point = np.asarray(point, dtype=float)
    u_value = sample_grid_field_fn(grid, u, point[0], point[1])
    v_value = sample_grid_field_fn(grid, v, point[0], point[1])
    if u_value is None or v_value is None:
        return None

    if display_sampler is not None:
        sampled_mapping = display_sampler.sample(point)
        if sampled_mapping is None:
            return None
        origin_display, jacobian = sampled_mapping
    else:
        jacobian = local_display_jacobian_fn(transform, point, grid)
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


def _density_xy(density):
    density_xy = np.broadcast_to(np.asarray(density, dtype=float), 2).astype(float)
    return np.maximum(density_xy, 0.1)


def _density_scalar(density):
    if np.isscalar(density):
        return max(float(density), 0.1)
    return float(np.mean(_density_xy(density)))


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
    sample_grid_field_array=None,
    thin_ncl_mapped_candidates=None,
):
    candidates, candidate_magnitudes = _prepare_ncl_center_candidates(
        grid=grid,
        magnitude=magnitude,
        density=density,
        start_points=start_points,
        ncl_preset=ncl_preset,
        sample_grid_field_array=sample_grid_field_array,
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
    selected_indices = thin_ncl_mapped_candidates(mapped_points, spacing_frac)
    return [(candidates[idx], candidate_magnitudes[idx]) for idx in selected_indices]


def _prepare_ncl_center_candidates(
    grid,
    magnitude,
    density,
    start_points,
    ncl_preset,
    sample_grid_field_array,
):
    if start_points is None:
        candidates = _default_ncl_box_center_candidates(
            grid,
            density=density,
            ncl_preset=ncl_preset,
        )
        candidate_magnitudes = sample_grid_field_array(grid, magnitude, candidates)
        return candidates, candidate_magnitudes

    candidates = np.asanyarray(start_points, dtype=float)
    if candidates.ndim != 2 or candidates.shape[1] != 2:
        raise ValueError("'start_points' must be an (N, 2) array")
    candidate_magnitudes = sample_grid_field_array(grid, magnitude, candidates)
    return candidates, candidate_magnitudes


def _default_ncl_candidate_shape(grid, density):
    density_xy = _density_xy(density)
    candidate_nx = min(grid.nx - 1, max(int(np.ceil(30.0 * density_xy[0])), 1))
    candidate_ny = min(grid.ny - 1, max(int(np.ceil(30.0 * density_xy[1])), 1))
    return candidate_nx, candidate_ny


def _default_ncl_box_center_candidates(grid, density=1, ncl_preset=None):
    del ncl_preset
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
    trace_ncl_direction_fn=None,
):
    if total_length_px <= 0:
        return None

    if anchor == "center":
        backward = trace_ncl_direction_fn(
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
        forward = trace_ncl_direction_fn(
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
        return trace_ncl_direction_fn(
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

    backward = trace_ncl_direction_fn(
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
    trace_ncl_curve_fn=None,
    evaluate_ncl_display_curve_fn=None,
):
    current_length_px = float(total_length_px)

    for _ in range(4):
        curve = trace_ncl_curve_fn(
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
            display_curve, transform_failed = evaluate_ncl_display_curve_fn(
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
    sample_local_vector_state_fn=None,
    ncl_step_length_px_fn=None,
    corrected_ncl_display_origin_fn=None,
    clip_display_step_to_viewport_fn=None,
    candidate_data_from_display_step_fn=None,
    point_within_grid_data_fn=None,
):
    start_point = np.asarray(start_point, dtype=float)
    initial_state = sample_local_vector_state_fn(
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
            state = sample_local_vector_state_fn(
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
            ncl_step_length_px_fn(step_px, current_speed, speed_scale),
        )
        if step_length_px <= 1e-6:
            break

        corrected_display = corrected_ncl_display_origin_fn(
            current_display, previous_display
        )
        candidate_display = (
            corrected_display + direction_sign * current_direction * step_length_px
        )
        candidate_display, clipped = clip_display_step_to_viewport_fn(
            corrected_display, candidate_display, viewport
        )
        candidate = candidate_data_from_display_step_fn(
            current_data=current_data,
            current_display=current_display,
            candidate_display=candidate_display,
            jacobian=current_jacobian,
            transform=transform,
        )
        if candidate is None or not point_within_grid_data_fn(grid, candidate):
            break

        next_state = sample_local_vector_state_fn(
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
                    ncl_step_length_px_fn(step_px, average_speed, speed_scale),
                )
                candidate_display = (
                    corrected_display
                    + average_direction / average_norm * step_length_px
                )
                candidate_display, clipped = clip_display_step_to_viewport_fn(
                    corrected_display,
                    candidate_display,
                    viewport,
                )
                candidate = candidate_data_from_display_step_fn(
                    current_data=current_data,
                    current_display=current_display,
                    candidate_display=candidate_display,
                    jacobian=0.5 * (current_jacobian + next_jacobian),
                    transform=transform,
                )
                if candidate is None or not point_within_grid_data_fn(grid, candidate):
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


def _curly_vector_ncl_impl(
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
    *,
    warn_if_native_backend_unavailable_fn=None,
    grid_cls=Grid,
    prepare_ncl_display_sampler_fn=None,
    prepare_ncl_native_trace_context_fn=None,
    select_ncl_centers_fn=None,
    build_ncl_curve_fn=None,
    sample_grid_field_fn=None,
    build_ncl_arrow_artists_fn=None,
    display_points_to_data_fn=None,
    result_cls=None,
):
    warn_if_native_backend_unavailable_fn()
    grid = grid_cls(x, y, allow_non_uniform=allow_non_uniform_grid)

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
        return result_cls(
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
    display_sampler = prepare_ncl_display_sampler_fn(grid, transform)
    native_trace_context = prepare_ncl_native_trace_context_fn(
        grid=grid,
        u=u,
        v=v,
        viewport=axes.bbox,
        display_sampler=display_sampler,
    )

    selected_centers = select_ncl_centers_fn(
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
        curve_result = build_ncl_curve_fn(
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
            sampled_width = sample_grid_field_fn(
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
            sampled_color = sample_grid_field_fn(
                grid, color_field, center[0], center[1]
            )
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
            baked_curve = display_points_to_data_fn(
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

        head_segments, arrow = build_ncl_arrow_artists_fn(
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
    axes.add_collection(lc, autolim=False)

    for patch in arrows:
        if rasterized is not None:
            patch.set_rasterized(bool(rasterized))
        axes.add_patch(patch)

    axes.autoscale_view()
    return result_cls(
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
