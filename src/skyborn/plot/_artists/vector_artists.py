"""Low-level vector artist sizing helpers for Skyborn plots."""

from __future__ import annotations

import numpy as np
from matplotlib import patches


def _ncl_arrow_edge_size_px(magnitude_value, max_mag, min_edge_px, max_edge_px):
    max_mag = max(float(max_mag), 1e-12)
    vmf = np.clip(float(magnitude_value) / max_mag, 0.0, 1.0)
    return float(min_edge_px) + (float(max_edge_px) - float(min_edge_px)) * vmf


def _resolve_open_arrow_size(edge_size_px):
    edge_size_px = max(float(edge_size_px), 1.0)
    half_angle = 0.5
    shaft_length_px = edge_size_px * np.cos(half_angle)
    head_width_px = 2.0 * edge_size_px * np.sin(half_angle)
    return shaft_length_px, head_width_px


def _uses_open_arrow_head(arrowstyle):
    return str(arrowstyle).strip() == "->"


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
