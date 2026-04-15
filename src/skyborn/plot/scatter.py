"""Scatter plotting and public API for display-space-thinned stippling.

Author: Qianye Su <suqianye2000@gmail.com>
Copyright (c) 2025-2026 Qianye Su
Created: 2026-04-03
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt

from ._core import scatter_engine as _scatter_engine
from ._shared.axes import _is_cartopy_crs_like, _looks_like_axes

__all__ = ["scatter"]


def _array_scatter(
    ax,
    x: Any,
    y: Any,
    s: Any = None,
    c: Any = None,
    **kwargs: Any,
):
    """Compatibility wrapper around the internal scatter rendering engine."""
    return _scatter_engine._scatter_impl(
        ax,
        x,
        y,
        s=s,
        c=c,
        **kwargs,
    )


def scatter(*args: Any, **kwargs: Any):
    """Scatter points with optional NCL-style display-space thinning.

    This is the public plotting entry point exposed by ``skyborn.plot``.
    It keeps a Matplotlib-compatible calling convention while adding the
    display-space thinning workflow needed for gridded stippling use cases such
    as significance masks on map projections and vertical cross-sections.

    Unlike a simple ``x[::step], y[::step]`` subsampling strategy, Skyborn
    first transforms candidate points into the current display geometry and then
    thins them there. That means the visible stipple density responds to the
    actual projection, axes aspect ratio, and panel extent seen by the user.

    Supported call styles
    ---------------------
    Matplotlib-style
        ``scatter(ax, x, y, ...)``
        ``scatter(x, y, ..., ax=ax)``
        ``scatter(x, y, ...)``
        ``scatter(x, y, s, c, ...)``

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Target axes. If omitted, ``matplotlib.pyplot.gca()`` is used.
    x, y : array-like
        Coordinate specification for the points to draw. Supported forms are:

        - paired 1D point coordinates, where ``x[i]`` and ``y[i]`` already
          describe one point each,
        - 1D rectilinear grid axes, which are expanded internally with
          ``numpy.meshgrid`` when a gridded mask or style field is supplied,
        - or 2D meshgrid-like coordinate arrays.
    s : scalar or array-like, optional
        Marker size passed through to ``matplotlib.axes.Axes.scatter``. If an
        array is supplied, it may be defined on the full grid, on the masked
        candidate set, or on the already flattened point list.
    c : color-like or array-like, optional
        Marker color argument passed through to ``Axes.scatter``. Array-like
        color fields follow the same subsetting rules as ``s``.
    color : color-like or array-like, optional
        Compatibility alias for ``c``. This is useful when matching plotting
        helpers that prefer the ``color=`` spelling while still preserving
        Matplotlib scatter semantics for numeric color fields.
    linewidth : scalar or array-like, optional
        Compatibility alias for ``linewidths``.
    facecolor : color-like or array-like, optional
        Compatibility alias for ``facecolors``.
    edgecolor : color-like or array-like, optional
        Compatibility alias for ``edgecolors``.
    linewidths, facecolors, edgecolors
        Native Matplotlib scatter styling arguments. Array-like values follow
        the same retained-point subsetting rules as ``s`` and ``c``.
    vmin, vmax : float, optional
        Colormap range controls forwarded to ``Axes.scatter``.
    where, mask : array-like, optional
        Candidate-selection mask. Use only one of them. For gridded stippling,
        these are typically 2D boolean or numeric fields aligned with the input
        grid. Numeric masks follow NumPy truthiness, while NaN values are
        treated as invalid rather than truthy.
    density : float, optional
        Relative stipple density. Larger values retain more points. When
        omitted, gridded inputs use the default NCL-style spacing rule and
        paired 1D points keep all points.
    distance : float, optional
        Explicit display-space thinning distance in normalized viewport units.
        Use this when an exact spacing threshold is preferred over the relative
        ``density`` control.
    min_distance : float, optional
        Backward-compatible alias for ``distance``.
    placement : {"auto", "points", "cells"}, optional
        Controls where gridded stipple candidates are generated before
        thinning. ``"points"`` keeps the original node-based behavior.
        ``"cells"`` fills the selected grid cells with interior candidates so
        dots can appear between coordinate centers, which is closer to NCL
        stipple fill behavior. The default ``"auto"`` enables ``"cells"`` for
        masked grids with inferable cell geometry when spacing is controlled by
        ``density`` or the default NCL-style rule, and otherwise falls back to
        ``"points"``.
    transform : optional
        Source coordinate transform. Standard Matplotlib transforms are passed
        through directly. Cartopy CRS-like objects are converted to the
        matching Matplotlib transform automatically.
    zorder : float, optional
        Matplotlib z-order of the generated scatter collection.
    **kwargs
        Additional keyword arguments forwarded to ``Axes.scatter`` after any
        array-like values have been subset to the retained points.

    Returns
    -------
    matplotlib.collections.PathCollection
        The scatter collection returned by ``Axes.scatter``.
    """

    if not args:
        raise TypeError("scatter() expects at least x and y positional arguments")

    ax = kwargs.pop("ax", None)
    remaining_args = list(args)

    if remaining_args and _looks_like_axes(remaining_args[0]):
        if ax is not None:
            raise TypeError("scatter() received Axes both positionally and via ax=")
        ax = remaining_args.pop(0)

    if ax is None:
        ax = plt.gca()

    if len(remaining_args) < 2:
        raise TypeError("scatter() requires x and y arguments")
    if len(remaining_args) > 4:
        raise TypeError("scatter() received too many positional arguments")

    x = remaining_args.pop(0)
    y = remaining_args.pop(0)
    s = remaining_args.pop(0) if remaining_args else kwargs.pop("s", None)
    c = remaining_args.pop(0) if remaining_args else kwargs.pop("c", None)

    transform = kwargs.pop("transform", None)
    if _is_cartopy_crs_like(transform):
        # Keep the public API ergonomic for Cartopy users: they can pass a CRS
        # object just like in ordinary Matplotlib/Cartopy plotting code.
        transform = transform._as_mpl_transform(ax)

    return _array_scatter(
        ax,
        x,
        y,
        s=s,
        c=c,
        transform=transform,
        **kwargs,
    )
