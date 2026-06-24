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
    _install_arrow_remove_hook,
    _validate_contour_direction,
    _validate_positive_float,
    _validate_positive_int,
)
from ._shared.axes import _looks_like_axes

__all__ = ["arrow_contour", "shadow_contourf"]


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


def _set_path_effects(contour_set: Any, effects: List[Any]) -> None:
    if hasattr(contour_set, "set_path_effects"):
        contour_set.set_path_effects(effects)

    for collection in getattr(contour_set, "collections", ()):
        if hasattr(collection, "set_path_effects"):
            collection.set_path_effects(effects)


def _hide_contour_artists(contour_set: Any) -> None:
    if hasattr(contour_set, "set_visible"):
        contour_set.set_visible(False)

    for collection in getattr(contour_set, "collections", ()):
        collection.set_visible(False)


def _overlay_shadow_kwargs(
    kwargs: dict[str, Any], color: Any, alpha: float
) -> dict[str, Any]:
    shadow_kwargs = dict(kwargs)
    shadow_kwargs.pop("cmap", None)
    shadow_kwargs.pop("norm", None)
    shadow_kwargs.pop("vmin", None)
    shadow_kwargs.pop("vmax", None)
    shadow_kwargs["colors"] = color
    shadow_kwargs["alpha"] = alpha

    zorder = shadow_kwargs.get("zorder", None)
    if isinstance(zorder, (int, float)):
        shadow_kwargs["zorder"] = zorder - 0.1

    return shadow_kwargs


def _apply_overlay_shadow_style(
    contour_set: Any,
    ax: Any,
    offset: tuple[float, float],
    blur: float,
) -> None:
    dx, dy = offset
    transform = contour_set.get_transform()
    offset_transform = mtransforms.ScaledTranslation(
        dx / 72.0,
        dy / 72.0,
        ax.figure.dpi_scale_trans,
    )
    contour_set.set_transform(transform + offset_transform)

    agg_filter = _blur_filter_factory(float(blur))
    if agg_filter is not None and hasattr(contour_set, "set_agg_filter"):
        contour_set.set_agg_filter(agg_filter)

    for collection in getattr(contour_set, "collections", ()):
        collection.set_transform(collection.get_transform() + offset_transform)
        if agg_filter is not None and hasattr(collection, "set_agg_filter"):
            collection.set_agg_filter(agg_filter)


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

        fill_collection = PathCollection(
            [path],
            facecolors=[facecolor],
            edgecolors=mpl.rcParams["hatch.color"] if hatch else "none",
            hatch=hatch,
            alpha=fill_alpha,
            transform=transform,
            zorder=layer_zorder + 0.01,
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
    arrow_max_length : float, default: 10.0
        Maximum arrowhead side length in points.
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
            arrow_max_length=arrow_max_length,
            positive_direction=positive_direction,
            arrow_color=arrow_color,
            arrow_linewidth=arrow_linewidth,
            zorder=arrow_zorder,
        )
        _install_arrow_remove_hook(contour_set, arrow_artists)
        _hide_contour_artists(contour_set)

    contour_set._skyborn_arrow_contour_artists = arrow_artists
    contour_set._skyborn_contour_arrows = arrow_artists
    return contour_set


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
    shadow_method : {"layered", "path_effect", "overlay"}, default: "layered"
        ``"layered"`` computes contour geometry once, then interleaves one
        shadow collection and one filled collection per contour layer. This
        preserves the internal stepped-shadow effect without recomputing
        ``contourf`` for every level. ``"path_effect"`` uses one ``contourf``
        call and draws a simpler outer shadow through Matplotlib's path-effect
        renderer.
        ``"overlay"`` draws one additional contour set underneath the main plot
        and can apply a soft Gaussian blur through an Agg filter.
    shadow_offset : tuple of float, default: (2.0, -2.0)
        Shadow offset in points. Positive x moves right and positive y moves up.
    shadow_alpha : float, default: 0.35
        Shadow opacity.
    shadow_color : color-like, default: "black"
        Shadow color.
    shadow_rho : float, default: 0.3
        Matplotlib ``SimplePatchShadow`` rho value used by ``shadow_method`` =
        ``"path_effect"``.
    shadow_blur : float, default: 1.2
        Blur radius in points used by ``shadow_method="layered"`` and
        ``shadow_method="overlay"``. Use ``0`` for a sharper and faster
        stepped shadow.
    shadow_backend : {"standard", "fast", "auto"}, default: "standard"
        Geometry backend used for the main filled contour when
        ``shadow_method="layered"``. ``"standard"`` follows Matplotlib's
        ordinary ``Axes.contourf`` path. ``"fast"`` precomputes supported
        filled-contour geometry before building the visible layers. ``"auto"``
        uses the fast path when supported and otherwise falls back to the
        standard path.

    Returns
    -------
    matplotlib.contour.QuadContourSet
        The main filled contour set.
    """
    ax = kwargs.pop("ax", None)
    shadow = bool(kwargs.pop("shadow", True))
    shadow_method = kwargs.pop("shadow_method", "layered")
    shadow_offset = _validate_shadow_offset(kwargs.pop("shadow_offset", (2.0, -2.0)))
    shadow_alpha = float(kwargs.pop("shadow_alpha", 0.35))
    shadow_color = kwargs.pop("shadow_color", "black")
    shadow_rho = float(kwargs.pop("shadow_rho", 0.3))
    shadow_blur = float(kwargs.pop("shadow_blur", 1.2))
    if "shadow_backend" in kwargs and "shadow_engine" in kwargs:
        raise TypeError(
            "shadow_contourf() received both shadow_backend and shadow_engine"
        )
    shadow_backend_name = (
        "shadow_engine" if "shadow_engine" in kwargs else "shadow_backend"
    )
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

    if shadow_method not in {"layered", "path_effect", "overlay"}:
        raise ValueError(
            "shadow_method must be one of: 'layered', 'path_effect', 'overlay'"
        )
    shadow_set = None
    contour_kwargs = dict(kwargs)
    if shadow and shadow_method == "overlay":
        shadow_set = ax.contourf(
            *remaining_args,
            **_overlay_shadow_kwargs(contour_kwargs, shadow_color, shadow_alpha),
        )
        _apply_overlay_shadow_style(shadow_set, ax, shadow_offset, shadow_blur)

    contour_set = None
    if shadow_method == "layered" and shadow_backend in {"fast", "auto"}:
        contour_set = _contourpy_contourf(ax, remaining_args, contour_kwargs)

    if contour_set is None:
        contour_set = ax.contourf(*remaining_args, **contour_kwargs)

    shadow_artists: List[Any] = []
    if shadow and shadow_method == "layered":
        shadow_artists = _add_layered_shadow_artists(
            contour_set,
            ax,
            shadow_offset,
            shadow_color,
            shadow_alpha,
            shadow_blur,
        )
        _install_layered_remove_hook(contour_set, shadow_artists)
    elif shadow and shadow_method == "path_effect":
        effects = [
            mpatheffects.SimplePatchShadow(
                offset=shadow_offset,
                shadow_rgbFace=shadow_color,
                alpha=shadow_alpha,
                rho=shadow_rho,
            ),
            mpatheffects.Normal(),
        ]
        _set_path_effects(contour_set, effects)

    contour_set._skyborn_shadow_contour_set = shadow_set
    contour_set._skyborn_shadow_method = shadow_method if shadow else None
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
