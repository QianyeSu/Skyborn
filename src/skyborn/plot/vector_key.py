"""Public reference-key facade for Skyborn curly-vector plots.

Author: Qianye Su <suqianye2000@gmail.com>
Copyright (c) 2025-2026 Qianye Su
Created: 2026-03-01 14:58:56
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt

from ._artists.vector_key_artist import CurlyVectorKey
from ._shared.axes import _looks_like_axes

__all__ = ["CurlyVectorKey", "curly_vector_key"]


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
