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


def _assemble_open_head_streamlines(
    *,
    shafts,
    shaft_colors,
    shaft_linewidths,
    head_segments,
    source_positions,
    use_multicolor_lines,
    use_linewidth_field,
):
    """Merge shaft curves and batched open-head segments into line-collection inputs.

    Parameters
    ----------
    shafts : sequence of ndarray
        Data-space shaft polylines, one per glyph.
    shaft_colors : sequence
        Per-glyph shaft colors aligned with ``shafts``.
    shaft_linewidths : sequence of float
        Per-glyph shaft linewidths aligned with ``shafts``.
    head_segments : ndarray
        Batched head segments with shape ``(N, 2, 2)`` in data coordinates.
    source_positions : array-like of int
        Source glyph index for each segment.
    use_multicolor_lines : bool
        Whether a per-segment color array should be emitted.
    use_linewidth_field : bool
        Whether a per-segment linewidth array should be emitted.

    Returns
    -------
    tuple
        ``(streamlines, line_colors, line_widths)`` where ``streamlines`` is a
        list of line segments/polyline arrays ready for ``LineCollection`` and
        the style arrays are either lists or ``None`` depending on the flags.
    """

    shaft_count = len(shafts)
    head_segments = np.asarray(head_segments, dtype=float)
    source_positions = np.asarray(source_positions, dtype=int)

    if source_positions.shape != (len(head_segments),):
        raise ValueError("source_positions must match head_segments")

    if len(source_positions) > 1 and np.any(
        source_positions[1:] < source_positions[:-1]
    ):
        order = np.argsort(source_positions, kind="stable")
        source_positions = source_positions[order]
        head_segments = head_segments[order]

    if len(source_positions) == 0:
        counts = np.zeros(shaft_count, dtype=np.intp)
    else:
        counts = np.bincount(source_positions, minlength=shaft_count)
    offsets = np.zeros(shaft_count + 1, dtype=np.intp)
    np.cumsum(counts, out=offsets[1:])

    total_streamlines = shaft_count + len(head_segments)
    streamlines = [None] * total_streamlines
    line_colors = [None] * total_streamlines if use_multicolor_lines else None
    line_widths = [None] * total_streamlines if use_linewidth_field else None

    cursor = 0
    for idx in range(shaft_count):
        seg_start = int(offsets[idx])
        seg_end = int(offsets[idx + 1])
        seg_count = seg_end - seg_start

        streamlines[cursor] = shafts[idx]
        if seg_count > 0:
            streamlines[cursor + 1 : cursor + 1 + seg_count] = list(
                head_segments[seg_start:seg_end]
            )

        style_span = 1 + seg_count
        if line_colors is not None:
            line_colors[cursor : cursor + style_span] = [shaft_colors[idx]] * style_span
        if line_widths is not None:
            line_widths[cursor : cursor + style_span] = [
                shaft_linewidths[idx]
            ] * style_span

        cursor += style_span

    return streamlines, line_colors, line_widths


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


def _assemble_filled_head_artists(
    *,
    shafts,
    shaft_colors,
    shaft_linewidths,
    filled_polygons,
    polygon_sources,
    facecolors,
    edgecolors,
    use_multicolor_lines,
    use_linewidth_field,
):
    """Build shaft and polygon style payloads for filled-arrow batch rendering.

    Parameters
    ----------
    shafts : sequence of ndarray
        Data-space shaft polylines, one per glyph.
    shaft_colors : sequence
        Per-glyph shaft colors aligned with ``shafts``.
    shaft_linewidths : sequence of float
        Per-glyph shaft linewidths aligned with ``shafts``.
    filled_polygons : ndarray
        Batched filled head polygons with shape ``(N, 3, 2)``.
    polygon_sources : array-like of int
        Source glyph index for each polygon.
    facecolors, edgecolors : sequence
        Per-glyph head colors aligned with ``shafts``.
    use_multicolor_lines : bool
        Whether a per-shaft line color list should be emitted.
    use_linewidth_field : bool
        Whether a per-shaft line width list should be emitted.

    Returns
    -------
    tuple
        ``(streamlines, line_colors, line_widths, polygon_facecolors,
        polygon_edgecolors, polygon_linewidths)``.
    """

    streamlines = list(shafts)
    line_colors = list(shaft_colors) if use_multicolor_lines else None
    line_widths = list(shaft_linewidths) if use_linewidth_field else None

    polygon_sources = np.asarray(polygon_sources, dtype=int)
    if polygon_sources.shape != (len(filled_polygons),):
        raise ValueError("polygon_sources must match filled_polygons")
    if len(polygon_sources) == 0:
        return streamlines, line_colors, line_widths, [], [], []

    polygon_facecolors = [facecolors[int(source_idx)] for source_idx in polygon_sources]
    polygon_edgecolors = [edgecolors[int(source_idx)] for source_idx in polygon_sources]
    polygon_linewidths = [
        max(float(shaft_linewidths[int(source_idx)]) * 0.5, 0.5)
        for source_idx in polygon_sources
    ]
    return (
        streamlines,
        line_colors,
        line_widths,
        polygon_facecolors,
        polygon_edgecolors,
        polygon_linewidths,
    )
