"""Result containers for Skyborn vector plotting."""

from __future__ import annotations

from typing import Any

import numpy as np

from .._artists.vector_artists import _ncl_arrow_edge_size_px, _resolve_open_arrow_size
from .vector_engine import _curve_length_from_magnitude


class CurlyVectorPlotSet:
    """Container returned by :func:`skyborn.plot.curly_vector`."""

    def __init__(
        self,
        lines: Any,
        arrows: Any,
        resolution: float,
        magnitude: Any,
        zorder: float | None,
        transform: Any,
        axes: Any,
        linewidth: Any,
        color: Any,
        cmap: Any,
        arrowsize: float,
        arrowstyle: str,
        start_points: Any,
        integration_direction: str,
        grains: Any,
        broken_streamlines: bool,
        allow_non_uniform_grid: bool = False,
        density: Any = None,
        anchor: str | None = None,
        ncl_preset: str | None = None,
        length_scale: Any = None,
        rasterized: bool | None = None,
    ) -> None:
        self.lines = lines
        self.arrows = tuple(arrows) if arrows is not None else ()
        self.resolution = resolution
        self.magnitude = magnitude
        self.zorder = zorder
        self.transform = transform
        self.axes = axes
        self.linewidth = linewidth
        self.color = color
        self.cmap = cmap
        self.arrowsize = arrowsize
        self.arrowstyle = arrowstyle
        self.start_points = start_points
        self.integration_direction = integration_direction
        self.grains = grains
        self.broken_streamlines = broken_streamlines
        self.allow_non_uniform_grid = allow_non_uniform_grid
        self.density = density
        self.render_mode = "ncl_curly"
        self.anchor = anchor
        self.ncl_preset = ncl_preset
        self.length_scale = length_scale
        self.ref_length_fraction = resolution
        self.rasterized = rasterized

        self.max_magnitude = (
            float(np.nanmax(magnitude))
            if magnitude is not None and np.isfinite(magnitude).any()
            else None
        )
        self.scale_factor = 2.0 if integration_direction == "both" else 1.0

    def get_scale_factor(self) -> float:
        return self.scale_factor

    def scale_value(self, value: float) -> float:
        return value / self.scale_factor

    def unscale_value(self, value: float) -> float:
        return value * self.scale_factor

    def glyph_length_axes_fraction(self, magnitude_value: float) -> float:
        try:
            magnitude_value = float(magnitude_value)
        except (TypeError, ValueError):
            return 0.0

        if not np.isfinite(magnitude_value) or magnitude_value <= 0.0:
            return 0.0

        if self.length_scale is not None and self.axes is not None:
            axes_width_px = max(float(self.axes.bbox.width), 1.0)
            length_px = _curve_length_from_magnitude(magnitude_value, self.length_scale)
            return max(float(length_px) / axes_width_px, 0.0)

        if self.max_magnitude is not None and float(self.max_magnitude) > 0.0:
            return max(
                float(self.resolution)
                * np.clip(magnitude_value / float(self.max_magnitude), 0.0, 1.0),
                0.0,
            )

        return 0.0

    def glyph_head_axes_dimensions(self, magnitude_value: float) -> tuple[float, float]:
        try:
            magnitude_value = float(magnitude_value)
        except (TypeError, ValueError):
            return 0.0, 0.0

        if (
            self.axes is None
            or not np.isfinite(magnitude_value)
            or magnitude_value <= 0.0
        ):
            return 0.0, 0.0

        axes_width_px = max(float(self.axes.bbox.width), 1.0)
        axes_height_px = max(float(self.axes.bbox.height), 1.0)
        max_mag = (
            float(self.max_magnitude)
            if self.max_magnitude is not None
            else magnitude_value
        )
        max_mag = max(max_mag, 1e-12)
        arrowsize = max(float(getattr(self, "arrowsize", 1.0)), 0.1)

        min_edge_px = max(axes_width_px * 0.003 * arrowsize, 1.2)
        max_edge_px = max(axes_width_px * 0.012 * arrowsize, min_edge_px)
        head_length_px, head_width_px = _resolve_open_arrow_size(
            _ncl_arrow_edge_size_px(
                magnitude_value,
                max_mag=max_mag,
                min_edge_px=min_edge_px,
                max_edge_px=max_edge_px,
            )
        )
        return (
            max(float(head_length_px) / axes_width_px, 0.0),
            max(float(head_width_px) / axes_height_px, 0.0),
        )
