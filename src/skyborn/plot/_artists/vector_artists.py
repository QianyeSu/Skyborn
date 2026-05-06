"""Low-level vector artist sizing helpers for Skyborn plots."""

from __future__ import annotations

import numpy as np


def _pack_display_curves(display_curves):
    """Pack variable-length display curves into one flat point array plus offsets."""

    count = len(display_curves)
    offsets = np.zeros(count + 1, dtype=np.intp)
    normalized_curves = []
    total_points = 0
    for idx, display_curve in enumerate(display_curves):
        curve = np.asarray(display_curve, dtype=float)
        if curve.ndim != 2 or curve.shape[1] != 2:
            curve = np.empty((0, 2), dtype=float)
        normalized_curves.append(curve)
        total_points += len(curve)
        offsets[idx + 1] = total_points

    packed_points = np.empty((total_points, 2), dtype=float)
    cursor = 0
    for curve in normalized_curves:
        next_cursor = cursor + len(curve)
        packed_points[cursor:next_cursor] = curve
        cursor = next_cursor
    return packed_points, offsets


def _build_open_arrow_segments_batch(
    *,
    transform,
    display_curves,
    head_lengths_px,
    head_widths_px,
    inverse_transform=None,
    build_open_arrow_segments_batch_fn=None,
):
    """Build many open-arrow head line segments in one batch.

    The runtime batch implementation is native-only. The Python reference
    implementation lives under ``tests/`` and is used for correctness checks
    rather than as a production fallback path.
    """

    if build_open_arrow_segments_batch_fn is None:
        return np.empty((0, 2, 2), dtype=float), np.empty(0, dtype=int)
    packed_points, curve_offsets = _pack_display_curves(display_curves)
    native_result = build_open_arrow_segments_batch_fn(
        packed_points,
        curve_offsets,
        np.asarray(head_lengths_px, dtype=float),
        np.asarray(head_widths_px, dtype=float),
    )
    if native_result is None:
        return np.empty((0, 2, 2), dtype=float), np.empty(0, dtype=int)

    display_segments, source_positions = native_result
    if len(display_segments) == 0:
        return np.empty((0, 2, 2), dtype=float), np.empty(0, dtype=int)
    try:
        inverse = (
            transform.inverted() if inverse_transform is None else inverse_transform
        )
        data_points = inverse.transform(
            np.asarray(display_segments, dtype=float).reshape(-1, 2)
        )
    except Exception:
        return np.empty((0, 2, 2), dtype=float), np.empty(0, dtype=int)

    data_segments = np.asarray(data_points, dtype=float).reshape(-1, 2, 2)
    finite_valid = np.isfinite(data_segments).all(axis=(1, 2))
    if not np.any(finite_valid):
        return np.empty((0, 2, 2), dtype=float), np.empty(0, dtype=int)
    return (
        data_segments[finite_valid],
        np.asarray(source_positions, dtype=int)[finite_valid],
    )


def _build_filled_arrow_polygons_batch(
    *,
    transform,
    display_curves,
    head_lengths_px,
    head_widths_px,
    inverse_transform=None,
    build_filled_arrow_polygons_batch_fn=None,
    display_points_to_data_fn=None,
):
    """Build many filled-arrow polygons in one batch.

    The runtime batch implementation is native-only. The Python reference
    implementation lives under ``tests/`` and is used for correctness checks
    rather than as a production fallback path.
    """

    count = len(display_curves)
    head_lengths_px = np.asarray(head_lengths_px, dtype=float)
    head_widths_px = np.asarray(head_widths_px, dtype=float)
    if head_lengths_px.shape != (count,) or head_widths_px.shape != (count,):
        raise ValueError("head_lengths_px and head_widths_px must match display_curves")
    if count == 0:
        return np.empty((0, 3, 2), dtype=float), np.empty(0, dtype=int)

    if build_filled_arrow_polygons_batch_fn is None:
        return np.empty((0, 3, 2), dtype=float), np.empty(0, dtype=int)
    packed_points, curve_offsets = _pack_display_curves(display_curves)
    native_result = build_filled_arrow_polygons_batch_fn(
        packed_points,
        curve_offsets,
        head_lengths_px,
        head_widths_px,
    )
    if native_result is None:
        return np.empty((0, 3, 2), dtype=float), np.empty(0, dtype=int)

    display_polygons, source_positions = native_result
    data_points = display_points_to_data_fn(
        transform,
        np.asarray(display_polygons, dtype=float).reshape(-1, 2),
        inverse_transform=inverse_transform,
    )
    if data_points is None:
        return np.empty((0, 3, 2), dtype=float), np.empty(0, dtype=int)

    data_polygons = np.asarray(data_points, dtype=float).reshape(-1, 3, 2)
    finite_valid = np.isfinite(data_polygons).all(axis=(1, 2))
    if not np.any(finite_valid):
        return np.empty((0, 3, 2), dtype=float), np.empty(0, dtype=int)
    return (
        data_polygons[finite_valid],
        np.asarray(source_positions, dtype=int)[finite_valid],
    )
