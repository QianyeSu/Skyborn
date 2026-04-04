"""Low-level vector artist sizing helpers for Skyborn plots."""

from __future__ import annotations

import numpy as np
from matplotlib import patches


def _uses_open_arrow_head(arrowstyle):
    return str(arrowstyle).strip() == "->"


def _open_arrow_geometry(
    curve,
    transform,
    head_length_px,
    head_width_px=None,
    display_curve=None,
    display_sampler=None,
    tip_display_geometry_from_display_curve_fn=None,
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

    tip_geometry = tip_display_geometry_from_display_curve_fn(
        display_curve,
        head_length_px * 1.35,
    )
    if tip_geometry is None:
        return None
    tip_display, unit = tip_geometry
    base_center = tip_display - unit * head_length_px

    geometry = {
        "display_curve": display_curve,
        "tip_display": tip_display,
        "base_center_display": base_center,
        "unit": unit,
    }
    if head_width_px is not None:
        normal = np.array([-unit[1], unit[0]])
        geometry["left_display"] = base_center + normal * head_width_px / 2.0
        geometry["right_display"] = base_center - normal * head_width_px / 2.0
    return geometry


def _trim_curve_for_open_head(
    curve,
    transform,
    head_length_px,
    display_sampler=None,
    open_arrow_geometry_fn=None,
    trim_display_curve_from_end_fn=None,
):
    curve = np.asarray(curve, dtype=float)
    if len(curve) < 2 or head_length_px <= 1e-6:
        return curve

    geometry = open_arrow_geometry_fn(
        curve,
        transform,
        head_length_px,
        display_sampler=display_sampler,
    )
    if geometry is None:
        return curve
    display_curve = geometry["display_curve"]

    trimmed_display = trim_display_curve_from_end_fn(display_curve, head_length_px)
    if trimmed_display is None or len(trimmed_display) < 2:
        return curve
    trimmed_display = trimmed_display.copy()
    trimmed_display[-1] = geometry["base_center_display"]

    try:
        trimmed_curve = transform.inverted().transform(trimmed_display)
    except Exception:
        return curve
    if not np.all(np.isfinite(trimmed_curve)):
        return curve
    return trimmed_curve


def _build_open_arrow_segments(
    curve,
    grid,
    transform,
    head_length_px,
    head_width_px,
    display_curve=None,
    display_sampler=None,
    inverse_transform=None,
    open_arrow_geometry_fn=None,
    display_points_to_data_fn=None,
):
    del grid
    geometry = open_arrow_geometry_fn(
        curve,
        transform,
        head_length_px,
        head_width_px,
        display_curve=display_curve,
        display_sampler=display_sampler,
    )
    if geometry is None:
        return []

    display_vertices = np.vstack(
        [
            geometry["left_display"],
            geometry["tip_display"],
            geometry["right_display"],
        ]
    )

    data_vertices = display_points_to_data_fn(
        transform,
        display_vertices,
        inverse_transform=inverse_transform,
    )
    if data_vertices is None:
        return []

    left, tip_data, right = data_vertices
    return [
        np.vstack([left, tip_data]),
        np.vstack([right, tip_data]),
    ]


def _build_arrow_polygon(
    curve,
    grid,
    transform,
    head_length_px,
    head_width_px,
    facecolor,
    edgecolor,
    linewidth,
    alpha,
    zorder,
    display_curve=None,
    display_sampler=None,
    inverse_transform=None,
    tip_display_geometry_fn=None,
    display_points_to_data_fn=None,
):
    del grid
    geometry = tip_display_geometry_fn(
        curve,
        transform,
        head_length_px * 1.25,
        display_curve=display_curve,
        display_sampler=display_sampler,
    )
    if geometry is None:
        return None

    tip_display, unit = geometry
    normal = np.array([-unit[1], unit[0]])
    base_center = tip_display - unit * head_length_px
    display_vertices = np.vstack(
        [
            tip_display,
            base_center + normal * head_width_px / 2.0,
            base_center - normal * head_width_px / 2.0,
        ]
    )

    data_vertices = display_points_to_data_fn(
        transform,
        display_vertices,
        inverse_transform=inverse_transform,
    )
    if data_vertices is None:
        return None

    return patches.Polygon(
        data_vertices,
        closed=True,
        transform=transform,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=max(float(linewidth) * 0.5, 0.5),
        alpha=alpha,
        zorder=zorder,
    )


def _build_ncl_arrow_artists(
    curve,
    grid,
    transform,
    arrowstyle,
    head_length_px,
    head_width_px,
    facecolor,
    edgecolor,
    linewidth,
    alpha,
    zorder,
    display_curve=None,
    display_sampler=None,
    inverse_transform=None,
    uses_open_arrow_head_fn=None,
    build_open_arrow_segments_fn=None,
    build_arrow_polygon_fn=None,
):
    if uses_open_arrow_head_fn(arrowstyle):
        return (
            build_open_arrow_segments_fn(
                curve=curve,
                grid=grid,
                transform=transform,
                head_length_px=head_length_px,
                head_width_px=head_width_px,
                display_curve=display_curve,
                display_sampler=display_sampler,
                inverse_transform=inverse_transform,
            ),
            None,
        )

    return (
        [],
        build_arrow_polygon_fn(
            curve=curve,
            grid=grid,
            transform=transform,
            head_length_px=head_length_px,
            head_width_px=head_width_px,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            alpha=alpha,
            zorder=zorder,
            display_curve=display_curve,
            display_sampler=display_sampler,
            inverse_transform=inverse_transform,
        ),
    )
