"""Native-helper adapters extracted from ``vector_plot.py``."""

from __future__ import annotations

import numpy as np


def _native_field_array(field):
    if np.ma.isMaskedArray(field):
        return np.asarray(np.ma.filled(field, np.nan), dtype=float)
    return np.asarray(field, dtype=float)


def _call_native_sample_grid_field(native_sampler, grid, field, xd, yd):
    value = native_sampler(
        field=_native_field_array(field),
        x_origin=float(grid.x_origin),
        y_origin=float(grid.y_origin),
        dx=float(grid.dx),
        dy=float(grid.dy),
        x=float(xd),
        y=float(yd),
    )
    if value is None:
        return None

    value = float(value)
    return value if np.isfinite(value) else None


def _call_native_sample_grid_field_array(
    native_sampler,
    grid,
    field,
    points,
    expected_shape,
):
    sampled = native_sampler(
        field=_native_field_array(field),
        x_origin=float(grid.x_origin),
        y_origin=float(grid.y_origin),
        dx=float(grid.dx),
        dy=float(grid.dy),
        points=np.asarray(points, dtype=float),
    )
    if sampled is None:
        return None

    sampled = np.asarray(sampled, dtype=float)
    if tuple(sampled.shape) != tuple(expected_shape):
        return None
    sampled[~np.isfinite(sampled)] = np.nan
    return sampled


def _call_native_thin_ncl_mapped_candidates(
    native_thinner,
    mapped_points,
    spacing_frac,
):
    selected = native_thinner(
        mapped_points=np.asarray(mapped_points, dtype=float),
        spacing_frac=float(spacing_frac),
    )
    if selected is None:
        return None
    return np.asarray(selected, dtype=int).tolist()


def _call_native_thin_ncl_display_candidates(
    native_thinner,
    display_points,
    viewport,
    spacing_frac,
):
    selected = native_thinner(
        display_points=np.asarray(display_points, dtype=float),
        viewport_x0=float(viewport.x0),
        viewport_y0=float(viewport.y0),
        viewport_width=float(viewport.width),
        viewport_height=float(viewport.height),
        spacing_frac=float(spacing_frac),
    )
    if selected is None:
        return None
    return np.asarray(selected, dtype=int).tolist()


def _call_native_generate_cell_candidates(
    native_generator,
    corners,
    mapped_corners,
    source_points,
    spacing_fraction,
):
    generated = native_generator(
        corners=np.asarray(corners, dtype=float),
        mapped_corners=np.asarray(mapped_corners, dtype=float),
        source_points=np.asarray(source_points, dtype=float),
        spacing_frac=float(spacing_fraction),
    )
    candidate_points, source_positions = generated
    return (
        np.asarray(candidate_points, dtype=float),
        np.asarray(source_positions, dtype=int),
    )


def _call_native_build_open_arrow_segments(
    native_builder,
    display_points,
    curve_offsets,
    head_lengths_px,
    head_widths_px,
):
    """Build compact batched open-arrow display segments through the native kernel."""

    segments, source_positions = native_builder(
        display_points=np.asarray(display_points, dtype=float),
        curve_offsets=np.asarray(curve_offsets, dtype=np.intp),
        head_lengths_px=np.asarray(head_lengths_px, dtype=float),
        head_widths_px=np.asarray(head_widths_px, dtype=float),
    )
    segments = np.asarray(segments, dtype=float)
    source_positions = np.asarray(source_positions, dtype=int)
    if segments.ndim != 3 or segments.shape[1:] != (2, 2):
        return None
    if source_positions.shape != (segments.shape[0],):
        return None
    return segments, source_positions


def _call_native_build_filled_arrow_polygons(
    native_builder,
    display_points,
    curve_offsets,
    head_lengths_px,
    head_widths_px,
):
    """Build compact batched filled-arrow polygons through the native kernel."""

    polygons, source_positions = native_builder(
        display_points=np.asarray(display_points, dtype=float),
        curve_offsets=np.asarray(curve_offsets, dtype=np.intp),
        head_lengths_px=np.asarray(head_lengths_px, dtype=float),
        head_widths_px=np.asarray(head_widths_px, dtype=float),
    )
    polygons = np.asarray(polygons, dtype=float)
    source_positions = np.asarray(source_positions, dtype=int)
    if polygons.ndim != 3 or polygons.shape[1:] != (3, 2):
        return None
    if source_positions.shape != (polygons.shape[0],):
        return None
    return polygons, source_positions


def _call_native_trace_ncl_direction(
    native_tracer,
    native_trace_context,
    start_point,
    max_length_px,
    direction_sign,
    step_px,
    speed_scale,
):
    if native_trace_context is None:
        return None

    curve = native_tracer(
        **_trace_kwargs(
            native_trace_context,
            start_point,
            max_length_px,
            direction_sign,
            step_px,
            speed_scale,
        )
    )
    if curve is None:
        return None
    curve = np.asarray(curve, dtype=float)
    return (
        curve
        if curve.ndim == 2 and curve.shape[0] >= 2 and curve.shape[1] == 2
        else None
    )


def _trace_kwargs(
    native_trace_context,
    start_point,
    max_length_px,
    direction_sign,
    step_px,
    speed_scale,
):
    return {
        "u": native_trace_context.u,
        "v": native_trace_context.v,
        "display_grid": native_trace_context.display_grid,
        "cell_valid": native_trace_context.cell_valid,
        "x_origin": native_trace_context.x_origin,
        "y_origin": native_trace_context.y_origin,
        "dx": native_trace_context.dx,
        "dy": native_trace_context.dy,
        "start_x": float(start_point[0]),
        "start_y": float(start_point[1]),
        "max_length_px": float(max_length_px),
        "direction_sign": float(direction_sign),
        "step_px": float(step_px),
        "speed_scale": float(speed_scale),
        "viewport_x0": native_trace_context.viewport_x0,
        "viewport_y0": native_trace_context.viewport_y0,
        "viewport_x1": native_trace_context.viewport_x1,
        "viewport_y1": native_trace_context.viewport_y1,
        "max_steps": 512,
    }


def _call_native_trace_ncl_direction_with_display(
    native_tracer,
    native_trace_context,
    start_point,
    max_length_px,
    direction_sign,
    step_px,
    speed_scale,
):
    if native_trace_context is None:
        return None

    traced = native_tracer(
        **_trace_kwargs(
            native_trace_context,
            start_point,
            max_length_px,
            direction_sign,
            step_px,
            speed_scale,
        ),
        return_display=True,
    )
    if traced is None:
        return None

    curve, display_curve = traced
    curve = np.asarray(curve, dtype=float)
    display_curve = np.asarray(display_curve, dtype=float)
    if (
        curve.ndim != 2
        or display_curve.ndim != 2
        or curve.shape != display_curve.shape
        or curve.shape[0] < 2
        or curve.shape[1] != 2
    ):
        return None
    return curve, display_curve


def _call_native_validate_display_curve(native_validator, display_curve, viewport):
    return bool(
        native_validator(
            display_curve=np.asarray(display_curve, dtype=float),
            viewport_width=float(viewport.width) if viewport is not None else 0.0,
            viewport_height=float(viewport.height) if viewport is not None else 0.0,
        )
    )
