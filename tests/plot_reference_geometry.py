"""Independent Python reference geometry for curly-vector head tests.

This module is intentionally test-side only. It provides a stable reference
implementation for the open-line and filled-triangle head geometry so native
and runtime-path changes can still be checked against an algorithmically
independent baseline after runtime batch fallbacks are removed.
"""

from __future__ import annotations

import numpy as np


def point_at_arc_distance_from_end_reference(
    curve: np.ndarray, distance: float
) -> np.ndarray:
    """Return the point located ``distance`` back from the curve tip."""

    array = np.asarray(curve, dtype=float)
    if array.ndim != 2 or array.shape[1] != 2 or len(array) == 0:
        raise ValueError("curve must be a non-empty (N, 2) array")
    if len(array) == 1:
        return array[0]

    remaining = max(float(distance), 0.0)
    for idx in range(len(array) - 1, 0, -1):
        right = array[idx]
        left = array[idx - 1]
        segment = right - left
        segment_length = float(np.hypot(segment[0], segment[1]))
        if segment_length <= 1e-12:
            continue
        if remaining <= segment_length:
            return right - segment * (remaining / segment_length)
        remaining -= segment_length
    return array[0]


def open_arrow_geometry_reference(
    display_curve: np.ndarray,
    head_length_px: float,
    head_width_px: float,
) -> dict[str, np.ndarray] | None:
    """Independent display-space reference for open-arrow head geometry."""

    array = np.asarray(display_curve, dtype=float)
    if array.ndim != 2 or array.shape[1] != 2 or len(array) < 2:
        return None
    if not np.all(np.isfinite(array)):
        return None

    tip_display = array[-1]
    tail_display = point_at_arc_distance_from_end_reference(
        array, head_length_px * 1.35
    )
    direction = tip_display - tail_display
    norm = float(np.hypot(direction[0], direction[1]))
    if norm <= 1e-12:
        return None

    unit = direction / norm
    base_center = tip_display - unit * float(head_length_px)
    normal = np.array([-unit[1], unit[0]], dtype=float)
    return {
        "tip_display": tip_display,
        "base_center_display": base_center,
        "left_display": base_center + normal * float(head_width_px) / 2.0,
        "right_display": base_center - normal * float(head_width_px) / 2.0,
        "unit": unit,
    }


def filled_arrow_geometry_reference(
    display_curve: np.ndarray,
    head_length_px: float,
    head_width_px: float,
) -> dict[str, np.ndarray] | None:
    """Independent display-space reference for filled-arrow triangle geometry."""

    array = np.asarray(display_curve, dtype=float)
    if array.ndim != 2 or array.shape[1] != 2 or len(array) < 2:
        return None
    if not np.all(np.isfinite(array)):
        return None

    tip_display = array[-1]
    tail_display = point_at_arc_distance_from_end_reference(
        array, head_length_px * 1.25
    )
    direction = tip_display - tail_display
    norm = float(np.hypot(direction[0], direction[1]))
    if norm <= 1e-12:
        return None

    unit = direction / norm
    normal = np.array([-unit[1], unit[0]], dtype=float)
    base_center = tip_display - unit * float(head_length_px)
    return {
        "tip_display": tip_display,
        "left_display": base_center + normal * float(head_width_px) / 2.0,
        "right_display": base_center - normal * float(head_width_px) / 2.0,
    }
