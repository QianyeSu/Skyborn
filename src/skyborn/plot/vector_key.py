"""Reference-vector annotation artists for Skyborn curly-vector plots.

Author: Qianye Su <suqianye2000@gmail.com>
Copyright (c) 2025-2026 Qianye Su
Created: 2026-03-01 14:58:56
"""

from __future__ import annotations

from typing import Any, Literal

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.artist import Artist, allow_rasterization
from matplotlib.axes import Axes
from matplotlib.backend_bases import RendererBase
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle
from matplotlib.text import Text

from .vector_plot import CurlyVectorPlotSet, _normalize_supported_arrowstyle

__all__ = ["CurlyVectorKey", "curly_vector_key"]


def _looks_like_axes(value: Any) -> bool:
    return hasattr(value, "add_artist") and hasattr(value, "transAxes")


class CurlyVectorKey(Artist):
    """NCL-like boxed reference-vector annotation for curly-vector plots.

    The annotation prefers the active :class:`CurlyVectorPlotSet` glyph scaling.
    ``reference_speed`` and ``max_arrow_length`` are only used as fallback
    controls when the plot set cannot provide a usable reference-length mapping.
    """

    def __init__(
        self,
        ax: Axes,
        curly_vector_set: CurlyVectorPlotSet,
        U: float,
        units: str = "m/s",
        label: str | None = None,
        description: str | None = None,
        width: float = 0.15,
        height: float = 0.10,
        loc: Literal[
            "lower left", "lower right", "upper left", "upper right"
        ] = "lower right",
        x: float | None = None,
        y: float | None = None,
        labelpos: Literal["N", "S", "E", "W"] = "N",
        max_arrow_length: float = 0.08,
        arrow_props: dict[str, Any] | None = None,
        patch_props: dict[str, Any] | None = None,
        text_props: dict[str, Any] | None = None,
        padding: float = 0.012,
        margin: float = 0.02,
        reference_speed: float = 2.0,
        center_label: bool | None = None,
        frameon: bool = True,
        show_description: bool = True,
    ) -> None:
        super().__init__()

        if loc not in ["lower left", "lower right", "upper left", "upper right"]:
            raise ValueError(
                "loc must be one of ['lower left', 'lower right', "
                f"'upper left', 'upper right'], got {loc}"
            )
        if labelpos not in ["N", "S", "E", "W"]:
            raise ValueError("labelpos must be one of ['N', 'S', 'E', 'W']")
        if (x is None) != (y is None):
            raise ValueError("x and y must both be provided together")

        self.ax = ax
        self.curly_vector_set = curly_vector_set
        self.U = float(U)
        if not np.isfinite(self.U) or self.U <= 0.0:
            raise ValueError("U must be a finite positive reference magnitude")
        self.units = units
        self.label = label
        self.description = description
        self.labelpos = labelpos
        self.reference_speed = reference_speed
        self.margin = margin
        self.loc = loc
        self.width = float(width)
        self.height = float(height)
        self.x = None if x is None else float(x)
        self.y = None if y is None else float(y)
        self.max_arrow_length = float(max_arrow_length)
        self.padding = float(padding)
        self.frameon = bool(frameon)
        self.show_description = bool(show_description)
        self.show_units = units != ""
        if (
            self.x is not None
            and self.y is not None
            and (not np.isfinite(self.x) or not np.isfinite(self.y))
        ):
            raise ValueError("x and y must be finite axes coordinates")

        if center_label is None:
            self.center_label = bool(label is not None and description is None)
        else:
            self.center_label = bool(center_label)

        self.arrow_props = self._setup_arrow_props(arrow_props)
        self.patch_props = self._setup_patch_props(patch_props)
        self.text_props = self._setup_text_props(text_props)
        self.arrowstyle = _normalize_supported_arrowstyle(
            self.arrow_props.get("arrowstyle", "->")
        )
        self.arrow_props["arrowstyle"] = self.arrowstyle

        self.arrow_length = self._calculate_arrow_length()
        self.patch = Rectangle(
            xy=(0.0, 0.0),
            width=0.0,
            height=0.0,
            transform=ax.transAxes,
            clip_on=False,
            **self.patch_props,
        )
        self.patch.set_visible(self.frameon)
        line_props = self._setup_line_props()
        self.arrow = Line2D([], [], transform=ax.transAxes, clip_on=False, **line_props)
        self.head_left = Line2D(
            [], [], transform=ax.transAxes, clip_on=False, **line_props
        )
        self.head_right = Line2D(
            [], [], transform=ax.transAxes, clip_on=False, **line_props
        )
        self.head_fill = Polygon(
            np.zeros((3, 2), dtype=float),
            closed=True,
            transform=ax.transAxes,
            clip_on=False,
            visible=False,
            **self._setup_head_fill_props(),
        )
        self.text = Text(0.0, 0.0, "", transform=ax.transAxes, clip_on=False)
        self.text2 = Text(0.0, 0.0, "", transform=ax.transAxes, clip_on=False)

        self.set_zorder(10)
        self.ax.add_artist(self)

    def _format_reference_value(self) -> str:
        rounded = round(self.U)
        if np.isfinite(self.U) and abs(self.U - rounded) <= 1e-8:
            value_text = str(int(rounded))
        else:
            value_text = f"{self.U:g}"
        return value_text if not self.show_units else f"{value_text} {self.units}"

    def _resolve_text_blocks(self) -> tuple[str | None, str | None, str | None]:
        description_text: str | None
        if self.label is not None:
            main_text = str(self.label)
            description_text = None if self.center_label else self.description
        else:
            main_text = self._format_reference_value()
            description_text = None if self.center_label else self.description
            if (
                description_text is None
                and not self.center_label
                and self.show_description
            ):
                description_text = "Reference Vector"

        if not self.show_description:
            description_text = None

        if self.labelpos in ("E", "W"):
            parts = [main_text]
            if description_text:
                parts.append(str(description_text))
            return None, None, "\n".join(parts)

        if self.labelpos == "S":
            return description_text, main_text, None
        return main_text, description_text, None

    def _measure_text(
        self, text: str | None, renderer: RendererBase
    ) -> tuple[float, float]:
        if not text:
            return 0.0, 0.0

        temp = Text(0.0, 0.0, str(text), **self.text_props)
        temp.set_figure(self.figure or self.ax.figure)
        bbox = temp.get_window_extent(renderer=renderer)
        return (
            float(bbox.width) / max(float(self.ax.bbox.width), 1.0),
            float(bbox.height) / max(float(self.ax.bbox.height), 1.0),
        )

    def _calculate_head_geometry(self) -> tuple[float, float]:
        head_dims = getattr(self.curly_vector_set, "glyph_head_axes_dimensions", None)
        if callable(head_dims):
            try:
                head_length, head_width = head_dims(self.U)
            except (TypeError, ValueError):
                head_length, head_width = np.nan, np.nan
            if np.isfinite(head_length) and np.isfinite(head_width):
                head_length = max(float(head_length), 0.0)
                head_width = max(float(head_width), 0.0)
                if head_length > 0.0 and head_width > 0.0:
                    return (
                        min(head_length, self.arrow_length * 0.5),
                        head_width,
                    )

        fallback_length = min(
            max(self.arrow_length * 0.18, 0.010), self.arrow_length * 0.5
        )
        fallback_width = max(fallback_length * 1.1, 0.018)
        return fallback_length, fallback_width

    def _set_text_artist(
        self,
        artist: Text,
        text: str | None,
        x: float,
        y: float,
        ha: str,
        va: str,
    ) -> None:
        artist.set_text("" if text is None else str(text))
        artist.update(self.text_props)
        artist.set_position((x, y))
        artist.set_horizontalalignment(ha)
        artist.set_verticalalignment(va)
        artist.set_visible(bool(text))
        artist.set_clip_on(False)

    def _update_arrow_geometry(
        self,
        start_x: float,
        tip_x: float,
        center_y: float,
        head_length: float,
        head_width: float,
    ) -> None:
        tip_x = max(float(tip_x), float(start_x))
        shaft_end_x = max(start_x, tip_x - max(head_length, 0.0))
        if self.arrowstyle == "->":
            # Keep the shaft running beneath the open head so the key arrow
            # reads as one connected glyph instead of a detached line plus head.
            shaft_end_x = tip_x
        self.arrow.set_data([start_x, shaft_end_x], [center_y, center_y])
        self.arrow.set_visible(True)

        if self.arrowstyle == "->":
            base_x = tip_x - head_length
            self.head_left.set_data(
                [base_x, tip_x], [center_y + head_width / 2.0, center_y]
            )
            self.head_right.set_data(
                [base_x, tip_x], [center_y - head_width / 2.0, center_y]
            )
            self.head_left.set_visible(True)
            self.head_right.set_visible(True)
            self.head_fill.set_visible(False)
            return

        self.head_left.set_visible(False)
        self.head_right.set_visible(False)
        self.head_fill.set_xy(
            np.asarray(
                [
                    [tip_x, center_y],
                    [tip_x - head_length, center_y + head_width / 2.0],
                    [tip_x - head_length, center_y - head_width / 2.0],
                ],
                dtype=float,
            )
        )
        self.head_fill.set_visible(True)

    def _layout_vertical(
        self,
        renderer: RendererBase,
        top_text: str | None,
        bottom_text: str | None,
    ) -> None:
        top_width, top_height = self._measure_text(top_text, renderer)
        bottom_width, bottom_height = self._measure_text(bottom_text, renderer)
        head_length, head_width = self._calculate_head_geometry()
        arrow_height = max(head_width, 0.010)
        content_gap = max(self.padding * 0.8, 0.008)
        top_block = top_height + (content_gap if top_text else 0.0)
        bottom_block = bottom_height + (content_gap if bottom_text else 0.0)
        content_width = max(self.arrow_length, top_width, bottom_width)
        content_height = top_block + arrow_height + bottom_block

        width = max(self.width, content_width + 2.0 * self.padding)
        height = max(self.height, content_height + 2.0 * self.padding)
        x, y = self._calculate_position(width, height)
        center_x = x + width / 2.0
        center_y = y + self.padding + bottom_block + arrow_height / 2.0

        self.patch.set_bounds(x, y, width, height)
        self._update_arrow_geometry(
            center_x - self.arrow_length / 2.0,
            center_x + self.arrow_length / 2.0,
            center_y,
            head_length,
            head_width,
        )

        if top_text:
            self._set_text_artist(
                self.text,
                top_text,
                center_x,
                y + height - self.padding,
                "center",
                "top",
            )
        else:
            self._set_text_artist(
                self.text, None, center_x, center_y, "center", "center"
            )

        if bottom_text:
            self._set_text_artist(
                self.text2,
                bottom_text,
                center_x,
                y + self.padding,
                "center",
                "bottom",
            )
        else:
            self._set_text_artist(
                self.text2, None, center_x, center_y, "center", "center"
            )

    def _layout_horizontal(self, renderer: RendererBase, side_text: str | None) -> None:
        side_width, side_height = self._measure_text(side_text, renderer)
        head_length, head_width = self._calculate_head_geometry()
        arrow_height = max(head_width, 0.010)
        content_gap = max(self.padding * 0.9, 0.010)
        has_side = bool(side_text)
        content_width = self.arrow_length + (
            content_gap + side_width if has_side else 0.0
        )
        content_height = max(arrow_height, side_height)

        width = max(self.width, content_width + 2.0 * self.padding)
        height = max(self.height, content_height + 2.0 * self.padding)
        x, y = self._calculate_position(width, height)
        center_y = y + height / 2.0

        self.patch.set_bounds(x, y, width, height)
        if self.labelpos == "W" and has_side:
            arrow_start = x + self.padding + side_width + content_gap
            text_x = x + self.padding + side_width
            text_ha = "right"
        else:
            arrow_start = x + self.padding
            text_x = arrow_start + self.arrow_length + content_gap
            text_ha = "left"

        self._update_arrow_geometry(
            arrow_start,
            arrow_start + self.arrow_length,
            center_y,
            head_length,
            head_width,
        )
        self._set_text_artist(
            self.text,
            side_text,
            text_x,
            center_y,
            text_ha,
            "center",
        )
        self._set_text_artist(self.text2, None, text_x, center_y, text_ha, "center")

    def _update_layout(self, renderer: RendererBase) -> None:
        top_text, bottom_text, side_text = self._resolve_text_blocks()
        if self.labelpos in ("E", "W"):
            self._layout_horizontal(renderer, side_text)
        else:
            self._layout_vertical(renderer, top_text, bottom_text)

        self.patch.set_visible(self.frameon)
        self.patch.set_zorder(self.get_zorder() - 1)
        for artist in (
            self.arrow,
            self.head_left,
            self.head_right,
            self.head_fill,
            self.text,
            self.text2,
        ):
            artist.set_zorder(self.get_zorder())

    def _calculate_arrow_length(self) -> float:
        """Calculate annotation arrow length from the active glyph scaling."""
        try:
            glyph_length = getattr(
                self.curly_vector_set, "glyph_length_axes_fraction", None
            )
            if callable(glyph_length):
                try:
                    arrow_length = float(glyph_length(self.U))
                except (TypeError, ValueError):
                    arrow_length = np.nan
                if np.isfinite(arrow_length) and arrow_length > 0.0:
                    return max(arrow_length, 0.012)

            max_magnitude = getattr(self.curly_vector_set, "max_magnitude", None)
            ref_length_fraction = getattr(
                self.curly_vector_set, "ref_length_fraction", None
            )
            if max_magnitude is not None and ref_length_fraction is not None:
                try:
                    arrow_length = np.clip(
                        float(self.U) / max(float(max_magnitude), 1e-12), 0.0, 1.0
                    ) * float(ref_length_fraction)
                except (TypeError, ValueError):
                    arrow_length = np.nan
                if np.isfinite(arrow_length) and arrow_length > 0.0:
                    return max(arrow_length, 0.012)

            reference_speed = max(float(getattr(self, "reference_speed", 2.0)), 1e-12)
            scale_factor = float(self.U) / reference_speed
            arrow_length = min(
                scale_factor * self.max_arrow_length, self.max_arrow_length * 4.0
            )
            return max(arrow_length, self.max_arrow_length * 0.2)

        except Exception:
            return self.max_arrow_length * 0.6

    def _calculate_position(
        self,
        width: float | None = None,
        height: float | None = None,
    ) -> tuple[float, float]:
        width = self.width if width is None else float(width)
        height = self.height if height is None else float(height)
        if self.x is not None and self.y is not None:
            return self.x, self.y
        margin = getattr(self, "margin", 0.02)
        if self.loc == "lower left":
            return margin, margin
        if self.loc == "lower right":
            return 1.0 - width - margin, margin
        if self.loc == "upper left":
            return margin, 1.0 - height - margin
        return 1.0 - width - margin, 1.0 - height - margin

    def _setup_props(
        self,
        user_props: dict[str, Any] | None,
        defaults: dict[str, Any],
    ) -> dict[str, Any]:
        merged = defaults.copy()
        if user_props:
            merged.update(user_props)
        return merged

    def _setup_arrow_props(self, arrow_props: dict[str, Any] | None) -> dict[str, Any]:
        defaults = {
            "arrowstyle": getattr(self.curly_vector_set, "arrowstyle", "->"),
            "linewidth": 1.5,
            "color": "black",
        }
        return self._setup_props(arrow_props, defaults)

    def _setup_patch_props(self, patch_props: dict[str, Any] | None) -> dict[str, Any]:
        defaults = {
            "linewidth": 1,
            "edgecolor": "black",
            "facecolor": "white",
            "alpha": 0.95,
        }
        return self._setup_props(patch_props, defaults)

    def _setup_text_props(self, text_props: dict[str, Any] | None) -> dict[str, Any]:
        defaults = {
            "fontsize": 10,
            "color": "black",
        }
        return self._setup_props(text_props, defaults)

    def _setup_line_props(self) -> dict[str, Any]:
        return {
            "linewidth": float(self.arrow_props.get("linewidth", 1.5)),
            "color": self.arrow_props.get("color", "black"),
            "alpha": self.arrow_props.get("alpha", 1.0),
            "solid_capstyle": "round",
            "solid_joinstyle": "round",
        }

    def _setup_head_fill_props(self) -> dict[str, Any]:
        color = self.arrow_props.get("color", "black")
        return {
            "facecolor": self.arrow_props.get("facecolor", color),
            "edgecolor": self.arrow_props.get("edgecolor", color),
            "linewidth": max(float(self.arrow_props.get("linewidth", 1.5)) * 0.75, 0.5),
            "alpha": self.arrow_props.get("alpha", 1.0),
        }

    def set_figure(self, fig: Figure) -> None:
        super().set_figure(fig)
        for artist in (
            self.patch,
            self.arrow,
            self.head_left,
            self.head_right,
            self.head_fill,
            self.text,
            self.text2,
        ):
            artist.set_figure(fig)

    @allow_rasterization
    def draw(self, renderer: RendererBase) -> None:
        if not self.get_visible():
            return

        self._update_layout(renderer)
        for artist in (
            self.patch,
            self.arrow,
            self.head_left,
            self.head_right,
            self.head_fill,
            self.text,
            self.text2,
        ):
            if artist.get_visible():
                artist.draw(renderer)
        self.stale = False


def curly_vector_key(
    *args: Any,
    **kwargs: Any,
) -> CurlyVectorKey:
    """Add an NCL-like reference-vector annotation to axes.

    This is the public companion to :func:`skyborn.plot.curly_vector`. It
    creates a boxed reference-vector annotation that reuses the active
    curly-vector length and head-size scaling.

    Supported call styles
    ---------------------
    - ``curly_vector_key(ax, curly_vector_set, U=..., ...)``
    - ``curly_vector_key(curly_vector_set, U=..., ax=ax, ...)``
    - ``curly_vector_key(curly_vector_set, U=..., ...)``

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Target axes. If omitted, ``matplotlib.pyplot.gca()`` is used.
    curly_vector_set : CurlyVectorPlotSet
        The object returned by :func:`skyborn.plot.curly_vector`.
    U : float, default: 2.0
        Finite positive reference magnitude represented by the annotation.
    units : str, default: ``"m/s"``
        Unit label appended to ``U``.
    label : str, optional
        Optional explicit main label. If omitted, the label is derived from
        ``U`` and ``units``.
    description : str, optional
        Secondary descriptive text shown with the reference vector.
    loc : {"lower left", "lower right", "upper left", "upper right"}, default: ``"lower right"``
        Corner placement used when ``x``/``y`` are not given.
    x, y : float, optional
        Explicit axes-fraction location of the annotation anchor. These behave
        similarly to ``quiverkey`` positional arguments.
    labelpos : {"N", "S", "E", "W"}, default: ``"N"``
        Label layout relative to the vector symbol.
    width, height : float, optional
        Box size in axes-fraction units.
    reference_speed, max_arrow_length : float, optional
        Fallback scaling controls used only when ``curly_vector_set`` cannot
        provide a valid glyph-length mapping for the requested reference
        magnitude.
    frameon : bool, default: True
        Whether to draw the surrounding box.
    show_description : bool, default: True
        Whether to render the secondary description text.

    Returns
    -------
    CurlyVectorKey
        The annotation artist added to the axes.
    """
    if not args and "curly_vector_set" not in kwargs:
        raise TypeError(
            "curly_vector_key() expects either (ax, curly_vector_set, ...) or "
            "(curly_vector_set, ...)"
        )

    ax_kwarg = kwargs.pop("ax", None)
    units = kwargs.pop("units", "m/s")
    label = kwargs.pop("label", None)
    loc = kwargs.pop("loc", "lower right")
    labelpos = kwargs.pop("labelpos", "N")

    remaining_args = list(args)
    if remaining_args and _looks_like_axes(remaining_args[0]):
        ax = remaining_args.pop(0)
    else:
        ax = ax_kwarg if ax_kwarg is not None else plt.gca()

    if remaining_args:
        curly_vector_set = remaining_args.pop(0)
    elif "curly_vector_set" in kwargs:
        curly_vector_set = kwargs.pop("curly_vector_set")
    else:
        raise TypeError("curly_vector_key() missing required curly_vector_set argument")

    if remaining_args:
        U = float(remaining_args.pop(0))
    else:
        U = float(kwargs.pop("U", 2.0))

    if remaining_args:
        units = str(remaining_args.pop(0))

    if remaining_args:
        raise TypeError("curly_vector_key() received too many positional arguments")

    return CurlyVectorKey(
        ax=ax,
        curly_vector_set=curly_vector_set,
        U=U,
        units=units,
        label=label,
        loc=loc,
        labelpos=labelpos,
        **kwargs,
    )
