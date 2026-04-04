"""Native-helper adapters extracted from ``vector_plot.py``."""

from __future__ import annotations

import numpy as np


def _try_native_sample_grid_field(native_sampler, grid, field, xd, yd, on_error):
    if native_sampler is None or np.ma.isMaskedArray(field):
        return None

    try:
        value = native_sampler(
            field=np.asarray(field, dtype=float),
            x_origin=float(grid.x_origin),
            y_origin=float(grid.y_origin),
            dx=float(grid.dx),
            dy=float(grid.dy),
            x=float(xd),
            y=float(yd),
        )
    except Exception as err:
        on_error("_sample_grid_field_native", "scalar sampling", err)
        return None

    if value is None:
        return None

    value = float(value)
    return value if np.isfinite(value) else None


def _try_native_sample_grid_field_array(
    native_sampler,
    grid,
    field,
    points,
    expected_shape,
    on_error,
):
    if native_sampler is None or np.ma.isMaskedArray(field):
        return None

    try:
        sampled = native_sampler(
            field=np.asarray(field, dtype=float),
            x_origin=float(grid.x_origin),
            y_origin=float(grid.y_origin),
            dx=float(grid.dx),
            dy=float(grid.dy),
            points=np.asarray(points, dtype=float),
        )
    except Exception as err:
        on_error("_sample_grid_field_array_native", "vectorized sampling", err)
        return None

    if sampled is None:
        return None

    sampled = np.asarray(sampled, dtype=float)
    if tuple(sampled.shape) != tuple(expected_shape):
        return None
    sampled[~np.isfinite(sampled)] = np.nan
    return sampled


def _try_native_thin_ncl_mapped_candidates(
    native_thinner,
    mapped_points,
    spacing_frac,
    on_error,
):
    if native_thinner is None:
        return None

    try:
        selected = native_thinner(
            mapped_points=np.asarray(mapped_points, dtype=float),
            spacing_frac=float(spacing_frac),
        )
    except Exception as err:
        on_error("_thin_ncl_mapped_candidates_native", "candidate thinning", err)
        return None

    if selected is None:
        return None
    return np.asarray(selected, dtype=int).tolist()


def _try_native_trace_ncl_direction(
    native_tracer,
    native_trace_context,
    start_point,
    max_length_px,
    direction_sign,
    step_px,
    speed_scale,
    on_error,
):
    if native_tracer is None or native_trace_context is None:
        return None

    try:
        curve = native_tracer(
            u=native_trace_context.u,
            v=native_trace_context.v,
            display_grid=native_trace_context.display_grid,
            cell_valid=native_trace_context.cell_valid,
            x_origin=native_trace_context.x_origin,
            y_origin=native_trace_context.y_origin,
            dx=native_trace_context.dx,
            dy=native_trace_context.dy,
            start_x=float(start_point[0]),
            start_y=float(start_point[1]),
            max_length_px=float(max_length_px),
            direction_sign=float(direction_sign),
            step_px=float(step_px),
            speed_scale=float(speed_scale),
            viewport_x0=native_trace_context.viewport_x0,
            viewport_y0=native_trace_context.viewport_y0,
            viewport_x1=native_trace_context.viewport_x1,
            viewport_y1=native_trace_context.viewport_y1,
            max_steps=512,
        )
    except Exception as err:
        on_error("_trace_ncl_direction_native", "trace", err)
        return None

    if curve is None:
        return None
    curve = np.asarray(curve, dtype=float)
    return (
        curve
        if curve.ndim == 2 and curve.shape[0] >= 2 and curve.shape[1] == 2
        else None
    )
