"""Contour plotting helpers."""

from __future__ import annotations

import types
from dataclasses import dataclass
from numbers import Integral
from typing import Any, List, Optional, Tuple, Union

import contourpy
import matplotlib as mpl
import matplotlib.patheffects as mpatheffects
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import numpy.ma as ma
from contourpy import FillType
from matplotlib.collections import PathCollection
from matplotlib.contour import ContourSet, QuadContourSet

from ._core.contour_arrows import (
    _add_contour_arrows,
    _arrow_contour_label_positions,
    _install_arrow_remove_hook,
    _validate_contour_direction,
    _validate_positive_float,
    _validate_positive_int,
)
from ._shared.axes import _looks_like_axes

__all__ = ["arrow_contour", "arrow_contour_clabel", "shadow_contourf"]


class _PrecomputedQuadContourSet(QuadContourSet):
    """QuadContourSet variant backed by already-computed contour polygons."""

    def _process_args(self, *args: Any, **kwargs: Any) -> dict[str, Any]:
        return ContourSet._process_args(self, *args, **kwargs)


@dataclass(frozen=True)
class _ContourpyCall:
    x: np.ndarray
    y: np.ndarray
    z: ma.MaskedArray
    levels_arg: Any = None


@dataclass(frozen=True)
class _ContourpyInput:
    x: np.ndarray
    y: np.ndarray
    z: ma.MaskedArray
    levels: np.ndarray


_CONTOURPY_FALLBACK_KWARGS = frozenset({"data", "filled", "locator"})


def _validate_shadow_offset(
    offset: Union[Tuple[float, float], List[float]],
) -> Tuple[float, float]:
    try:
        dx, dy = offset
    except (TypeError, ValueError) as exc:
        raise ValueError("shadow_offset must be a two-item (x, y) sequence") from exc
    return float(dx), float(dy)


def _normalize_shadow_backend(value: Any, name: str) -> str:
    backend = str(value).strip().lower()
    aliases = {
        "standard": "standard",
        "matplotlib": "standard",
        "fast": "fast",
        "contourpy": "fast",
        "auto": "auto",
    }
    if backend not in aliases:
        raise ValueError(f"{name} must be one of: 'standard', 'fast', or 'auto'")
    return aliases[backend]


def _blur_filter_factory(sigma_points: float):
    if sigma_points <= 0:
        return None

    def _blur_filter(image, dpi):
        from scipy.ndimage import gaussian_filter

        sigma_pixels = sigma_points * dpi / 72.0
        return gaussian_filter(image, sigma=(sigma_pixels, sigma_pixels, 0)), 0, 0

    return _blur_filter


def _apply_artist_filter(artist: Any, agg_filter: Any) -> None:
    if agg_filter is not None and hasattr(artist, "set_agg_filter"):
        artist.set_agg_filter(agg_filter)


def _hide_contour_artists(contour_set: Any) -> None:
    if hasattr(contour_set, "set_visible"):
        contour_set.set_visible(False)

    for collection in getattr(contour_set, "collections", ()):
        collection.set_visible(False)


def _iter_contour_layers(contour_set: Any):
    hatches = list(getattr(contour_set, "hatches", ()) or ())
    if hasattr(contour_set, "get_paths") and hasattr(contour_set, "get_facecolors"):
        paths = list(contour_set.get_paths())
        facecolors = list(contour_set.get_facecolors())
        for index, (path, facecolor) in enumerate(zip(paths, facecolors)):
            hatch = hatches[index % len(hatches)] if hatches else None
            yield path, facecolor, hatch
        return

    for index, collection in enumerate(getattr(contour_set, "collections", ())):
        facecolors = collection.get_facecolors()
        facecolor = facecolors[0] if len(facecolors) else None
        hatch = hatches[index % len(hatches)] if hatches else collection.get_hatch()
        for path in collection.get_paths():
            yield path, facecolor, hatch


def _path_touches_view_boundary(path: Any, ax: Any, margin: float = 0.1) -> bool:
    """Check if a path touches or extends beyond the axes view limits.

    Parameters
    ----------
    path : matplotlib.path.Path
        The path to check.
    ax : matplotlib.axes.Axes
        The axes containing the view limits.
    margin : float, default: 0.1
        Margin in data coordinates. Paths within this margin of the boundary
        are considered to touch it.

    Returns
    -------
    bool
        True if the path touches or extends beyond the view boundary.
    """
    vertices = path.vertices
    if len(vertices) == 0:
        return False

    x_coords = vertices[:, 0]
    y_coords = vertices[:, 1]

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()

    # Check if any vertex is near or beyond the boundaries
    touches = (
        np.any(x_coords <= x_min + margin)
        or np.any(x_coords >= x_max - margin)
        or np.any(y_coords <= y_min + margin)
        or np.any(y_coords >= y_max - margin)
    )

    return touches


def _add_layered_shadow_artists(
    contour_set: Any,
    ax: Any,
    offset: tuple[float, float],
    color: Any,
    alpha: float,
    blur: float,
) -> List[Any]:
    transform = contour_set.get_transform()
    base_zorder = contour_set.get_zorder()
    dx, dy = offset
    offset_transform = mtransforms.ScaledTranslation(
        dx / 72.0,
        dy / 72.0,
        ax.figure.dpi_scale_trans,
    )
    agg_filter = _blur_filter_factory(float(blur))

    artists: List[Any] = []
    fill_alpha = contour_set.get_alpha()
    for index, (path, facecolor, hatch) in enumerate(_iter_contour_layers(contour_set)):
        layer_zorder = base_zorder + index * 0.02

        # Only add shadow for paths that don't touch the view boundary
        # to avoid artifacts from truncated contours
        if not _path_touches_view_boundary(path, ax, margin=0.1):
            shadow_collection = PathCollection(
                [path],
                facecolors=color,
                edgecolors="none",
                alpha=alpha,
                transform=transform + offset_transform,
                zorder=layer_zorder,
            )
            _apply_artist_filter(shadow_collection, agg_filter)
            ax.add_collection(shadow_collection)
            artists.append(shadow_collection)

        # Use matching edge color to eliminate white gaps between adjacent layers
        # caused by antialiasing artifacts. Disable antialiasing on fills to
        # eliminate gaps, but keep edges smooth with matching color and width.
        edge_color = mpl.rcParams["hatch.color"] if hatch else facecolor

        # Compute adaptive edge width based on figure DPI for consistent appearance
        # across different display densities. Use 0 for hatched contours.
        if hatch:
            edge_width = 0
        else:
            # Get DPI-aware edge width: 1pt at 72dpi, scales with DPI
            # Clamp to [0.25, 1.0] points to avoid too thick/thin edges
            fig_dpi = ax.figure.dpi if hasattr(ax, "figure") else 72.0
            edge_width = np.clip(72.0 / fig_dpi, 0.25, 1.0)

        fill_collection = PathCollection(
            [path],
            facecolors=[facecolor],
            edgecolors=[edge_color],
            linewidths=edge_width,
            hatch=hatch,
            alpha=fill_alpha,
            transform=transform,
            zorder=layer_zorder + 0.01,
            antialiaseds=False,  # Disable AA to eliminate white line gaps
        )
        ax.add_collection(fill_collection)
        artists.append(fill_collection)

    _hide_contour_artists(contour_set)
    return artists


def _install_layered_remove_hook(contour_set: Any, artists: List[Any]) -> None:
    original_remove = contour_set.remove

    def _remove(self):
        for artist in list(artists):
            try:
                artist.remove()
            except ValueError:
                pass
        artists.clear()
        return original_remove()

    contour_set.remove = types.MethodType(_remove, contour_set)


def _initialize_contour_xy(
    z: ma.MaskedArray,
    origin: Any,
    extent: Any,
) -> tuple[np.ndarray, np.ndarray]:
    if z.ndim != 2:
        raise TypeError(f"Input z must be 2D, not {z.ndim}D")
    if z.shape[0] < 2 or z.shape[1] < 2:
        raise TypeError(
            f"Input z must be at least a (2, 2) shaped array, but has shape {z.shape}"
        )

    ny, nx = z.shape
    if origin is None:
        if extent is None:
            return np.meshgrid(np.arange(nx), np.arange(ny))
        x0, x1, y0, y1 = extent
        return np.meshgrid(np.linspace(x0, x1, nx), np.linspace(y0, y1, ny))

    if extent is None:
        x0, x1, y0, y1 = (0, nx, 0, ny)
    else:
        x0, x1, y0, y1 = extent

    dx = (x1 - x0) / nx
    dy = (y1 - y0) / ny
    x = x0 + (np.arange(nx) + 0.5) * dx
    y = y0 + (np.arange(ny) + 0.5) * dy
    if origin == "upper":
        y = y[::-1]
    return np.meshgrid(x, y)


def _read_z_contourpy_call(
    ax: Any,
    args: list[Any],
    kwargs: dict[str, Any],
) -> _ContourpyCall:
    del ax
    z = ma.asarray(args[0])
    x, y = _initialize_contour_xy(z, kwargs.get("origin"), kwargs.get("extent"))
    return _ContourpyCall(x=x, y=y, z=z)


def _read_z_levels_contourpy_call(
    ax: Any,
    args: list[Any],
    kwargs: dict[str, Any],
) -> _ContourpyCall:
    call = _read_z_contourpy_call(ax, args, kwargs)
    return _ContourpyCall(x=call.x, y=call.y, z=call.z, levels_arg=args[1])


def _check_contour_xyz(
    ax: Any,
    x: Any,
    y: Any,
    z: Any,
    kwargs: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray, ma.MaskedArray]:
    x, y = ax._process_unit_info([("x", x), ("y", y)], kwargs)
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    z = ma.asarray(z)

    if z.ndim != 2:
        raise TypeError(f"Input z must be 2D, not {z.ndim}D")
    if z.shape[0] < 2 or z.shape[1] < 2:
        raise TypeError(
            f"Input z must be at least a (2, 2) shaped array, but has shape {z.shape}"
        )

    ny, nx = z.shape
    if x.ndim != y.ndim:
        raise TypeError(
            f"Number of dimensions of x ({x.ndim}) and y ({y.ndim}) do not match"
        )
    if x.ndim == 1:
        if x.shape[0] != nx:
            raise TypeError(
                f"Length of x ({x.shape[0]}) must match number of columns in z ({nx})"
            )
        if y.shape[0] != ny:
            raise TypeError(
                f"Length of y ({y.shape[0]}) must match number of rows in z ({ny})"
            )
        x, y = np.meshgrid(x, y)
    elif x.ndim == 2:
        if x.shape != z.shape:
            raise TypeError(f"Shapes of x {x.shape} and z {z.shape} do not match")
        if y.shape != z.shape:
            raise TypeError(f"Shapes of y {y.shape} and z {z.shape} do not match")
    else:
        raise TypeError(f"Inputs x and y must be 1D or 2D, not {x.ndim}D")

    return x, y, z


def _read_xyz_contourpy_call(
    ax: Any,
    args: list[Any],
    kwargs: dict[str, Any],
) -> _ContourpyCall:
    x, y, z = _check_contour_xyz(ax, args[0], args[1], args[2], kwargs)
    return _ContourpyCall(x=x, y=y, z=z)


def _read_xyz_levels_contourpy_call(
    ax: Any,
    args: list[Any],
    kwargs: dict[str, Any],
) -> _ContourpyCall:
    call = _read_xyz_contourpy_call(ax, args, kwargs)
    return _ContourpyCall(x=call.x, y=call.y, z=call.z, levels_arg=args[3])


_CONTOURPY_CALL_READERS = {
    1: _read_z_contourpy_call,
    2: _read_z_levels_contourpy_call,
    3: _read_xyz_contourpy_call,
    4: _read_xyz_levels_contourpy_call,
}


def _auto_contour_levels(
    zmin: float,
    zmax: float,
    levels: int,
) -> np.ndarray:
    locator = mpl.ticker.MaxNLocator(levels + 1, min_n_ticks=1)
    values = np.asarray(locator.tick_values(zmin, zmax), dtype=np.float64)
    under = np.nonzero(values < zmin)[0]
    start = under[-1] if len(under) else 0
    over = np.nonzero(values > zmax)[0]
    stop = over[0] + 1 if len(over) else len(values)
    if stop - start < 3:
        start, stop = 0, len(values)
    return values[start:stop]


def _prepare_contourpy_z(
    z: ma.MaskedArray,
    norm: Any,
) -> tuple[ma.MaskedArray, bool]:
    z = ma.masked_invalid(z, copy=False)
    logscale = norm is not None and isinstance(norm, mpl.colors.LogNorm)
    if logscale:
        z = ma.masked_where(z <= 0, z)
    return z, logscale


def _resolve_contourpy_levels(
    z: ma.MaskedArray,
    levels_arg: Any,
    logscale: bool,
) -> Optional[np.ndarray]:
    zmin = float(z.min())
    zmax = float(z.max())

    if levels_arg is None:
        if np.issubdtype(z.dtype, bool):
            return np.asarray([0.0, 0.5, 1.0], dtype=np.float64)
        return None if logscale else _auto_contour_levels(zmin, zmax, 7)

    if isinstance(levels_arg, Integral):
        return None if logscale else _auto_contour_levels(zmin, zmax, int(levels_arg))

    return np.asarray(levels_arg, dtype=np.float64)


def _validate_contourpy_levels(levels: np.ndarray) -> None:
    if len(levels) < 2:
        raise ValueError("Filled contours require at least 2 levels.")
    if len(levels) > 1 and np.min(np.diff(levels)) <= 0.0:
        raise ValueError("Contour levels must be increasing")


def _resolve_contourpy_input(
    ax: Any,
    args: list[Any],
    kwargs: dict[str, Any],
) -> Optional[_ContourpyInput]:
    reader = _CONTOURPY_CALL_READERS.get(len(args))
    if reader is None:
        return None

    call = reader(ax, args, kwargs)
    keyword_levels = kwargs.get("levels", None)
    if call.levels_arg is not None and keyword_levels is not None:
        return None

    z, logscale = _prepare_contourpy_z(call.z, kwargs.get("norm"))
    levels_arg = keyword_levels if keyword_levels is not None else call.levels_arg
    levels = _resolve_contourpy_levels(z, levels_arg, logscale)
    if levels is None:
        return None

    _validate_contourpy_levels(levels)
    return _ContourpyInput(x=call.x, y=call.y, z=z, levels=levels)


def _contourpy_supported(kwargs: dict[str, Any]) -> bool:
    return bool(
        not _CONTOURPY_FALLBACK_KWARGS.intersection(kwargs)
        and kwargs.get("extend", "neither") == "neither"
    )


def _contourpy_generator_kwargs(kwargs: dict[str, Any]) -> dict[str, Any]:
    algorithm = kwargs.get("algorithm", mpl.rcParams["contour.algorithm"])
    mpl.rcParams.validate["contour.algorithm"](algorithm)
    corner_mask = kwargs.get("corner_mask", None)
    if corner_mask is None:
        corner_mask = (
            False if algorithm == "mpl2005" else mpl.rcParams["contour.corner_mask"]
        )

    generator_kwargs = {
        "name": algorithm,
        "corner_mask": corner_mask,
        "fill_type": FillType.OuterCode,
    }
    nchunk = int(kwargs.get("nchunk", 0) or 0)
    if nchunk > 0:
        generator_kwargs["chunk_size"] = nchunk
    return generator_kwargs


def _precomputed_contour_kwargs(kwargs: dict[str, Any]) -> dict[str, Any]:
    contour_kwargs = dict(kwargs)
    for key in ("levels", "algorithm", "corner_mask", "nchunk"):
        contour_kwargs.pop(key, None)
    return contour_kwargs


def _contourpy_contourf(
    ax: Any,
    args: list[Any],
    kwargs: dict[str, Any],
) -> Optional[QuadContourSet]:
    if not _contourpy_supported(kwargs):
        return None

    contour_input = _resolve_contourpy_input(ax, args, kwargs)
    if contour_input is None:
        return None

    contour_generator = contourpy.contour_generator(
        contour_input.x,
        contour_input.y,
        contour_input.z,
        **_contourpy_generator_kwargs(kwargs),
    )
    allsegs = []
    allkinds = []
    for lower, upper in zip(contour_input.levels[:-1], contour_input.levels[1:]):
        segs, kinds = contour_generator.filled(float(lower), float(upper))
        allsegs.append(segs)
        allkinds.append(kinds)
    if not any(segs for segs in allsegs):
        return None

    return _PrecomputedQuadContourSet(
        ax,
        contour_input.levels,
        allsegs,
        allkinds,
        filled=True,
        **_precomputed_contour_kwargs(kwargs),
    )


def arrow_contour(*args: Any, **kwargs: Any):
    """Draw contour lines with directional arrowheads along each line.

    Ordinary positional and keyword arguments are forwarded to ``Axes.contour``.
    For closed contour paths, positive levels are oriented clockwise and
    negative levels are oriented counterclockwise in the current displayed
    axes coordinates before arrowheads are placed.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Target axes. If omitted, ``matplotlib.pyplot.gca()`` is used. The axes
        may also be supplied as the first positional argument.
    arrows : bool, default: True
        Whether to add arrowheads.
    arrow_count : int, default: 1
        Number of arrowheads to place on each contour segment.
    arrow_size : float, default: 0.45
        Arrowhead width relative to the local arrowhead length.
    arrow_length_fraction : float, default: 0.035
        Fraction of the contour-segment length used for each arrow body.
        Mutually exclusive with ``arrow_length_points``.
    arrow_length_points : float, optional
        Absolute arrowhead length in points for cross-figure comparison. When
        specified, all arrows will have this exact length regardless of contour
        segment length, enabling consistent visual comparison across different
        plots. Takes priority over ``arrow_length_fraction`` if both are provided.
    arrow_max_length : float, default: 10.0
        Maximum arrowhead side length in points. Only used when
        ``arrow_length_fraction`` is active (ignored if ``arrow_length_points``
        is specified).
    positive_direction : {"clockwise", "counterclockwise"}, default: "clockwise"
        Direction used for positive closed contours. Negative closed contours
        use the opposite direction. Direction is evaluated in the current
        display coordinates, so set final axis limits or inverted pressure axes
        before calling ``arrow_contour``.
    arrow_color : color-like, optional
        Arrow color. Defaults to the matching contour line color.
    arrow_linewidth : float, optional
        Arrow line width. Defaults to the matching contour line width.
    arrow_zorder : float, optional
        Arrow z-order. Defaults to slightly above the contour line.
    arrow_style : {"swept", "line", "filled"}, default: "swept"
        Arrow rendering style.
        "swept" uses NCL/NCAR Graphics-style swept-back barbs (see
        ``drwvec.f`` in NCAR Graphics, default style): both barbs originate
        at the arrow tip and sweep back at ``arrow_angle`` degrees from the
        reversed shaft direction, giving a narrow, sharp look.
        "line" uses a flat-backed V-shaped line segment (original style),
        with both barbs fanning out perpendicular to the shaft from its
        base.
        "filled" uses solid triangular arrowheads.
    arrow_angle : float, default: 22.5
        Sweep-back angle in degrees for each barb, measured from the
        reversed shaft direction. Only used when ``arrow_style="swept"``.
        NCAR Graphics uses 22.5 degrees by default.
    arrow_min_spacing : float, optional
        Minimum spacing in display pixels between arrows, both within the same
        contour line and across nearby contour lines. When contour lines are
        densely packed (e.g., tight gradients), this prevents arrow crowding by
        skipping arrows that would be too close to existing ones. Similar to NCL's
        ``vcMinDistanceF``. If not specified, uses an automatic spacing based on
        arrow length (typically 2.5× arrow length).

    Returns
    -------
    matplotlib.contour.QuadContourSet
        The contour set returned by ``Axes.contour``. Added arrow artists are
        available as ``result._skyborn_arrow_contour_artists``.
    """
    ax = kwargs.pop("ax", None)
    arrows_enabled = bool(kwargs.pop("arrows", True))
    arrow_count = _validate_positive_int(kwargs.pop("arrow_count", 1), "arrow_count")
    arrow_size = _validate_positive_float(kwargs.pop("arrow_size", 0.45), "arrow_size")
    arrow_length_fraction = _validate_positive_float(
        kwargs.pop("arrow_length_fraction", 0.035),
        "arrow_length_fraction",
    )
    arrow_length_points = kwargs.pop("arrow_length_points", None)
    if arrow_length_points is not None:
        arrow_length_points = _validate_positive_float(
            arrow_length_points, "arrow_length_points"
        )
    arrow_max_length = _validate_positive_float(
        kwargs.pop("arrow_max_length", 10.0),
        "arrow_max_length",
    )
    positive_direction = _validate_contour_direction(
        kwargs.pop("positive_direction", "clockwise"),
        "positive_direction",
    )
    arrow_color = kwargs.pop("arrow_color", None)
    arrow_linewidth = kwargs.pop("arrow_linewidth", None)
    arrow_zorder = kwargs.pop("arrow_zorder", None)
    arrow_style = kwargs.pop("arrow_style", "swept")
    if arrow_style not in {"filled", "line", "swept"}:
        raise ValueError("arrow_style must be 'line', 'swept', or 'filled'")
    arrow_angle = _validate_positive_float(
        kwargs.pop("arrow_angle", 22.5), "arrow_angle"
    )
    arrow_min_spacing = kwargs.pop("arrow_min_spacing", None)
    if arrow_min_spacing is not None:
        arrow_min_spacing = _validate_positive_float(
            arrow_min_spacing, "arrow_min_spacing"
        )

    remaining_args = list(args)
    if remaining_args and _looks_like_axes(remaining_args[0]):
        if ax is not None:
            raise TypeError(
                "arrow_contour() received Axes both positionally and via ax="
            )
        ax = remaining_args.pop(0)

    if ax is None:
        ax = plt.gca()

    contour_set = ax.contour(*remaining_args, **kwargs)
    arrow_artists: List[Any] = []
    if arrows_enabled:
        arrow_artists = _add_contour_arrows(
            contour_set,
            ax,
            arrow_count=arrow_count,
            arrow_size=arrow_size,
            arrow_length_fraction=arrow_length_fraction,
            arrow_length_points=arrow_length_points,
            arrow_max_length=arrow_max_length,
            positive_direction=positive_direction,
            arrow_color=arrow_color,
            arrow_linewidth=arrow_linewidth,
            zorder=arrow_zorder,
            arrow_style=arrow_style,
            arrow_angle=arrow_angle,
            arrow_min_spacing=arrow_min_spacing,
        )
        _install_arrow_remove_hook(contour_set, arrow_artists)
        _hide_contour_artists(contour_set)

    contour_set._skyborn_arrow_contour_artists = arrow_artists
    contour_set._skyborn_contour_arrows = arrow_artists
    return contour_set


def arrow_contour_clabel(contour_set: Any, *args: Any, **kwargs: Any):
    """Label ``arrow_contour`` lines, placing labels clear of the arrowheads.

    This is a thin wrapper around ``Axes.clabel``. For each contour line
    produced by :func:`arrow_contour`, it picks a label anchor point (by arc
    length along the line) that is as far as possible from that line's
    arrowheads, then passes those anchors to ``clabel`` via ``manual=``. This
    mirrors NCL's behavior of routing contour labels around vector/arrow
    glyphs so labels and arrowheads do not overlap.

    Parameters
    ----------
    contour_set : matplotlib.contour.QuadContourSet
        The contour set returned by :func:`arrow_contour`.
    *args, **kwargs :
        Forwarded to ``Axes.clabel``. If ``manual`` is already supplied, it is
        used as-is and no automatic label placement is performed.

    Returns
    -------
    dict
        The mapping of label text artists returned by ``Axes.clabel``.
    """
    if "manual" not in kwargs:
        label_positions = _arrow_contour_label_positions(contour_set)
        if label_positions:
            kwargs["manual"] = label_positions
    return contour_set.axes.clabel(contour_set, *args, **kwargs)


def shadow_contourf(*args: Any, **kwargs: Any):
    """Draw a filled contour plot with an efficient drop-shadow effect.

    The positional and ordinary keyword arguments are forwarded to
    ``Axes.contourf``. The returned object is the main Matplotlib contour set,
    so it can still be passed to ``colorbar`` or adjusted like a normal
    ``contourf`` result.

    Parameters
    ----------
    ax : matplotlib.axes.Axes, optional
        Target axes. If omitted, ``matplotlib.pyplot.gca()`` is used. The axes
        may also be supplied as the first positional argument.
    shadow : bool, default: True
        Whether to draw the shadow.
    shadow_offset : tuple of float, default: (2.0, -2.0)
        Shadow offset in points. Positive x moves right and positive y moves up.
    shadow_alpha : float, default: 0.35
        Shadow opacity.
    shadow_color : color-like, default: "black"
        Shadow color.
    shadow_blur : float, default: 1.2
        Blur radius in points. Use ``0`` for a sharper and faster stepped shadow.
    shadow_backend : {"standard", "fast", "auto"}, default: "standard"
        Geometry backend used for the main filled contour. ``"standard"`` follows
        Matplotlib's ordinary ``Axes.contourf`` path. ``"fast"`` precomputes
        supported filled-contour geometry before building the visible layers.
        ``"auto"`` uses the fast path when supported and otherwise falls back to
        the standard path.

    Returns
    -------
    matplotlib.contour.QuadContourSet
        The main filled contour set.
    """
    ax = kwargs.pop("ax", None)
    shadow = bool(kwargs.pop("shadow", True))
    shadow_offset = _validate_shadow_offset(kwargs.pop("shadow_offset", (2.0, -2.0)))
    shadow_alpha = float(kwargs.pop("shadow_alpha", 0.35))
    shadow_color = kwargs.pop("shadow_color", "black")
    shadow_blur = float(kwargs.pop("shadow_blur", 1.2))
    if "shadow_backend" in kwargs and "shadow_engine" in kwargs:
        raise TypeError(
            "shadow_contourf() received both shadow_backend and shadow_engine"
        )
    if "shadow_engine" in kwargs:
        import warnings

        warnings.warn(
            "shadow_engine is deprecated, use shadow_backend instead",
            DeprecationWarning,
            stacklevel=2,
        )
        shadow_backend_name = "shadow_engine"
    else:
        shadow_backend_name = "shadow_backend"
    shadow_backend = _normalize_shadow_backend(
        kwargs.pop(shadow_backend_name, "standard"),
        shadow_backend_name,
    )

    remaining_args = list(args)
    if remaining_args and _looks_like_axes(remaining_args[0]):
        if ax is not None:
            raise TypeError(
                "shadow_contourf() received Axes both positionally and via ax="
            )
        ax = remaining_args.pop(0)

    if ax is None:
        ax = plt.gca()

    contour_kwargs = dict(kwargs)
    contour_set = None
    if shadow_backend in {"fast", "auto"}:
        contour_set = _contourpy_contourf(ax, remaining_args, contour_kwargs)

    if contour_set is None:
        contour_set = ax.contourf(*remaining_args, **contour_kwargs)

    shadow_artists: List[Any] = []
    if shadow:
        shadow_artists = _add_layered_shadow_artists(
            contour_set,
            ax,
            shadow_offset,
            shadow_color,
            shadow_alpha,
            shadow_blur,
        )
        _install_layered_remove_hook(contour_set, shadow_artists)

    contour_set._skyborn_shadow_artists = shadow_artists
    contour_set._skyborn_shadow_backend = (
        "fast" if isinstance(contour_set, _PrecomputedQuadContourSet) else "standard"
    )
    contour_set._skyborn_shadow_engine = (
        "contourpy"
        if isinstance(contour_set, _PrecomputedQuadContourSet)
        else "matplotlib"
    )
    return contour_set
