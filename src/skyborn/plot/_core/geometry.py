"""Display-space geometry helpers extracted from ``vector_plot.py``."""

from __future__ import annotations

import numpy as np


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


def _point_within_grid_data(grid, point):
    return (
        grid.x_origin <= point[0] <= grid.x_origin + grid.width
        and grid.y_origin <= point[1] <= grid.y_origin + grid.height
    )


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


def _tip_display_geometry_from_display_curve(display_curve, backoff_px):
    tip_display = display_curve[-1]
    tail_display = _point_at_arc_distance_from_end(display_curve, backoff_px)
    direction = tip_display - tail_display
    direction_norm = np.hypot(*direction)
    if direction_norm <= 1e-12:
        return None
    return tip_display, direction / direction_norm


def _tip_display_geometry(
    curve,
    transform,
    backoff_px,
    display_curve=None,
    display_sampler=None,
):
    del display_sampler
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
