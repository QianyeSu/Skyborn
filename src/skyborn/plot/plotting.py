from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Literal

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
from matplotlib.axes import Axes
from matplotlib.image import AxesImage
from matplotlib.patches import PathPatch
from matplotlib.path import Path

__all__ = [
    "GradientFillBetween",
    "add_equal_axes",
    "createFigure",
    "gradient_fill_between",
]


@dataclass
class GradientFillBetween:
    """Artists created by :func:`gradient_fill_between`."""

    images: list[AxesImage]
    clip_paths: list[PathPatch]

    def remove(self) -> None:
        """Remove the gradient images and their clipping patches from the axes."""
        for artist in [*self.images, *self.clip_paths]:
            try:
                artist.remove()
            except ValueError:
                pass


def add_equal_axes(ax, loc, pad, width):
    """
    Add a new Axes with equal height or width next to the original Axes and return that object.

    Parameters
    ----------
    ax : Axes or array_like of Axes
        The original Axes, or can be an array of Axes.

    loc : {'left', 'right', 'bottom', 'top'}
        Position of the new Axes relative to the old Axes.

    pad : float
        Spacing between the new Axes and the old Axes.

    width: float
        When loc='left' or 'right', width represents the width of the new Axes.
        When loc='bottom' or 'top', width represents the height of the new Axes.

    Returns
    -------
    ax_new : Axes
        New Axes object.
    """
    # Whether ax is a single Axes or a group of Axes, get the size and position of ax.
    axes = np.atleast_1d(ax).ravel()
    bbox = mtransforms.Bbox.union([ax.get_position() for ax in axes])

    # Determine the size and position of the new Axes.
    if loc == "left":
        x0_new = bbox.x0 - pad - width
        x1_new = x0_new + width
        y0_new = bbox.y0
        y1_new = bbox.y1
    elif loc == "right":
        x0_new = bbox.x1 + pad
        x1_new = x0_new + width
        y0_new = bbox.y0
        y1_new = bbox.y1
    elif loc == "bottom":
        x0_new = bbox.x0
        x1_new = bbox.x1
        y0_new = bbox.y0 - pad - width
        y1_new = y0_new + width
    elif loc == "top":
        x0_new = bbox.x0
        x1_new = bbox.x1
        y0_new = bbox.y1 + pad
        y1_new = y0_new + width
    else:
        raise ValueError(
            f"Invalid location '{loc}'. Must be one of: 'left', 'right', 'bottom', 'top'"
        )

    # Create new Axes.
    fig = axes[0].get_figure()
    bbox_new = mtransforms.Bbox.from_extents(x0_new, y0_new, x1_new, y1_new)
    ax_new = fig.add_axes(bbox_new)

    return ax_new


def createFigure(figsize=(12, 8), dpi=300, subplotAdj=None, **kwargs):
    figsize = figsize
    figure = plt.figure(figsize=figsize, dpi=dpi, **kwargs)
    if subplotAdj is not None:
        plt.subplots_adjust(**subplotAdj)
    return figure


def _convert_x_to_float(x: Any) -> np.ndarray:
    values = np.asarray(x)
    if np.issubdtype(values.dtype, np.datetime64):
        return mpl.dates.date2num(values.astype("datetime64[us]").astype(object))

    if values.dtype == object:
        try:
            return mpl.dates.date2num(values)
        except (TypeError, ValueError):
            pass

    return values.astype(np.float64)


def _iter_valid_runs(valid: np.ndarray) -> Iterable[slice]:
    indices = np.flatnonzero(valid)
    if indices.size == 0:
        return

    breaks = np.flatnonzero(np.diff(indices) > 1)
    starts = np.r_[indices[0], indices[breaks + 1]]
    stops = np.r_[indices[breaks] + 1, indices[-1] + 1]
    for start, stop in zip(starts, stops):
        if stop - start >= 2:
            yield slice(int(start), int(stop))


def _auto_gradient_resolution(ax: Axes, resolution: Any) -> tuple[int, int]:
    if resolution == "auto":
        bbox = ax.get_window_extent()
        nx = max(128, int(np.ceil(bbox.width)))
        ny = max(128, int(np.ceil(bbox.height)))
        return nx, ny

    if np.isscalar(resolution):
        value = int(resolution)
        if value < 2:
            raise ValueError("resolution must be at least 2")
        return value, value

    try:
        nx, ny = resolution
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "resolution must be 'auto', an integer, or an (nx, ny) pair"
        ) from exc

    nx = int(nx)
    ny = int(ny)
    if nx < 2 or ny < 2:
        raise ValueError("resolution values must be at least 2")
    return nx, ny


def _gradient_segment(
    ax: Axes,
    x: np.ndarray,
    y1: np.ndarray,
    y2: np.ndarray,
    *,
    cmap: Any,
    alpha: float,
    nx: int,
    ny: int,
    zorder: float | None,
    interpolation: str,
    edge_mode: Literal["auto", "y1", "y2"],
) -> tuple[AxesImage, PathPatch]:
    order = np.argsort(x)
    x = x[order]
    y1 = y1[order]
    y2 = y2[order]

    xmin, xmax = float(x[0]), float(x[-1])
    ymin = float(min(np.min(y1), np.min(y2)))
    ymax = float(max(np.max(y1), np.max(y2)))
    if xmin == xmax or ymin == ymax:
        raise ValueError("gradient_fill_between requires non-zero x and y ranges")

    xi = np.linspace(xmin, xmax, nx)
    yi = np.linspace(ymin, ymax, ny)
    y1i = np.interp(xi, x, y1)
    y2i = np.interp(xi, x, y2)

    if edge_mode == "y1":
        lower = y1i
        upper = y2i
    elif edge_mode == "y2":
        lower = y2i
        upper = y1i
    else:
        lower = np.minimum(y1i, y2i)
        upper = np.maximum(y1i, y2i)

    bottom = np.minimum(y1i, y2i)
    top = np.maximum(y1i, y2i)
    thickness = np.maximum(upper - lower, np.finfo(float).eps)
    values = (yi[:, None] - lower[None, :]) / thickness[None, :]
    values = np.ma.masked_where(
        (yi[:, None] < bottom[None, :]) | (yi[:, None] > top[None, :]),
        values,
    )
    values = np.ma.clip(values, 0.0, 1.0)

    image = ax.imshow(
        values,
        extent=(xmin, xmax, ymin, ymax),
        origin="lower",
        aspect="auto",
        cmap=cmap,
        alpha=alpha,
        interpolation=interpolation,
        zorder=zorder,
    )

    vertices = np.column_stack([x, y1]).tolist()
    vertices.extend(np.column_stack([x[::-1], y2[::-1]]).tolist())
    vertices.append(vertices[0])
    codes = [Path.MOVETO] + [Path.LINETO] * (len(vertices) - 2) + [Path.CLOSEPOLY]
    patch = PathPatch(
        Path(vertices, codes),
        facecolor="none",
        edgecolor="none",
        transform=ax.transData,
    )
    ax.add_patch(patch)
    image.set_clip_path(patch)
    return image, patch


def gradient_fill_between(
    ax: Axes,
    x: Any,
    y1: Any,
    y2: Any = 0,
    *,
    cmap: Any = "YlOrBr",
    alpha: float = 1.0,
    resolution: Any = "auto",
    interpolation: str = "bicubic",
    zorder: float | None = None,
    edge_mode: Literal["auto", "y1", "y2"] = "auto",
) -> GradientFillBetween:
    """Fill between two curves with a smooth relative gradient.

    Unlike stacking many ``fill_between`` strips, this draws a single gradient
    image per continuous valid segment and clips it to the polygon bounded by
    ``y1`` and ``y2``. The gradient is relative to the two curves at each x
    position: one boundary is mapped to the low end of the colormap and the
    other boundary to the high end.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    x, y1, y2 : array-like
        Coordinates of the two bounding curves. Datetime-like x values are
        supported.
    cmap : Colormap or str, default: "YlOrBr"
        Colormap used from the lower relative boundary to the upper boundary.
    alpha : float, default: 1.0
        Overall image alpha.
    resolution : "auto", int, or tuple(int, int), default: "auto"
        Pixel grid used for the clipped gradient. ``"auto"`` follows the
        current axes size in display pixels.
    interpolation : str, default: "bicubic"
        Matplotlib image interpolation mode.
    zorder : float, optional
        Image z-order.
    edge_mode : {"auto", "y1", "y2"}, default: "auto"
        Which curve receives the low end of the colormap. ``"auto"`` uses the
        lower y-value at each x. Use ``"y1"`` or ``"y2"`` when the semantic
        lower-color boundary should follow a named input curve even if curves
        cross.

    Returns
    -------
    GradientFillBetween
        Container with the created image and clipping-patch artists.
    """
    if not isinstance(ax, Axes):
        raise TypeError("ax must be a matplotlib.axes.Axes instance")
    if edge_mode not in {"auto", "y1", "y2"}:
        raise ValueError("edge_mode must be one of: 'auto', 'y1', or 'y2'")

    x_values = _convert_x_to_float(x).ravel()
    y1_values = np.asarray(y1, dtype=np.float64).ravel()
    y2_values = np.broadcast_to(
        np.asarray(y2, dtype=np.float64), y1_values.shape
    ).ravel()

    if x_values.shape != y1_values.shape:
        raise ValueError("x and y1 must have the same shape")
    if y2_values.shape != y1_values.shape:
        raise ValueError("y2 must be scalar or have the same shape as y1")

    nx, ny = _auto_gradient_resolution(ax, resolution)
    cmap = mpl.colormaps.get_cmap(cmap) if isinstance(cmap, str) else cmap
    valid = np.isfinite(x_values) & np.isfinite(y1_values) & np.isfinite(y2_values)

    images: list[AxesImage] = []
    clip_paths: list[PathPatch] = []
    for run in _iter_valid_runs(valid):
        image, patch = _gradient_segment(
            ax,
            x_values[run],
            y1_values[run],
            y2_values[run],
            cmap=cmap,
            alpha=float(alpha),
            nx=nx,
            ny=ny,
            zorder=zorder,
            interpolation=interpolation,
            edge_mode=edge_mode,
        )
        images.append(image)
        clip_paths.append(patch)

    if not images:
        raise ValueError("gradient_fill_between requires at least two finite points")

    return GradientFillBetween(images=images, clip_paths=clip_paths)
