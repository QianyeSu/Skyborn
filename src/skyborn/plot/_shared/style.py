"""Shared style-alias and kwargs helpers for Skyborn plots."""

from __future__ import annotations

from typing import Any

import numpy as np

_SUPPORTED_CURLY_ARROWSTYLES = ("->", "-|>")


def _collect_named_kwargs(
    scope: dict[str, Any], names: tuple[str, ...]
) -> dict[str, Any]:
    """Collect a stable subset of keyword arguments from a local scope."""
    return {name: scope[name] for name in names}


def _normalize_supported_arrowstyle(arrowstyle: Any) -> str:
    """Validate and normalize the supported curly-vector arrow styles."""
    normalized = str(arrowstyle).strip()
    if normalized not in _SUPPORTED_CURLY_ARROWSTYLES:
        supported = ", ".join(repr(value) for value in _SUPPORTED_CURLY_ARROWSTYLES)
        raise ValueError(f"arrowstyle must be one of {supported}; got {arrowstyle!r}")
    return normalized


def _normalize_artist_alpha(alpha: Any) -> float | None:
    """Validate a Matplotlib-style artist alpha value."""
    if alpha is None:
        return None

    alpha_value = float(alpha)
    if not np.isfinite(alpha_value) or not 0.0 <= alpha_value <= 1.0:
        raise ValueError("alpha must be a finite number between 0 and 1")
    return alpha_value


def _ncl_arrow_edge_size_px(
    magnitude_value: Any,
    max_mag: Any,
    min_edge_px: Any,
    max_edge_px: Any,
) -> float:
    """Scale the NCL-like arrow edge size between configured pixel bounds."""
    max_mag = max(float(max_mag), 1e-12)
    vmf = np.clip(float(magnitude_value) / max_mag, 0.0, 1.0)
    return float(min_edge_px) + (float(max_edge_px) - float(min_edge_px)) * vmf


def _resolve_open_arrow_size(edge_size_px: Any) -> tuple[float, float]:
    """Convert an NCL-like edge size into open-head length/width pixels."""
    edge_size_px = max(float(edge_size_px), 1.0)
    half_angle = 0.5
    shaft_length_px = edge_size_px * np.cos(half_angle)
    head_width_px = 2.0 * edge_size_px * np.sin(half_angle)
    return shaft_length_px, head_width_px


def _curly_head_axes_dimensions(
    axes: Any,
    magnitude_value: Any,
    *,
    max_magnitude: Any,
    arrowsize: Any,
) -> tuple[float, float]:
    """Resolve reference-key head dimensions in axes fractions."""
    try:
        magnitude_value = float(magnitude_value)
    except (TypeError, ValueError):
        return 0.0, 0.0

    if axes is None or not np.isfinite(magnitude_value) or magnitude_value <= 0.0:
        return 0.0, 0.0

    try:
        axes_width_px = max(float(axes.bbox.width), 1.0)
        axes_height_px = max(float(axes.bbox.height), 1.0)
    except Exception:
        return 0.0, 0.0

    try:
        max_mag = float(max_magnitude) if max_magnitude is not None else magnitude_value
    except (TypeError, ValueError):
        max_mag = magnitude_value
    if not np.isfinite(max_mag) or max_mag <= 0.0:
        max_mag = magnitude_value

    try:
        arrowsize = max(float(arrowsize), 0.1)
    except (TypeError, ValueError):
        arrowsize = 1.0

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


def _normalize_curly_pivot(pivot: Any) -> str | None:
    """Map Matplotlib quiver-style pivots to curly-vector anchors."""
    if pivot is None:
        return None

    normalized = str(pivot).strip().lower()
    mapping = {
        "tail": "tail",
        "mid": "center",
        "middle": "center",
        "tip": "head",
    }
    if normalized not in mapping:
        supported = ", ".join(repr(value) for value in mapping)
        raise ValueError(f"pivot must be one of {supported}; got {pivot!r}")
    return mapping[normalized]


def _resolve_curly_anchor_alias(anchor: Any, pivot: Any) -> Any:
    """Resolve Matplotlib-style ``pivot`` into the native ``anchor`` setting."""
    pivot_anchor = _normalize_curly_pivot(pivot)
    if anchor is None:
        return pivot_anchor
    if pivot_anchor is None:
        return anchor
    if str(anchor).strip().lower() != pivot_anchor:
        raise ValueError(
            "anchor and pivot refer to different glyph anchors; "
            f"got anchor={anchor!r} and pivot={pivot!r}"
        )
    return anchor


def _resolve_curly_style_aliases(
    *,
    color: Any,
    c: Any,
    linewidth: Any,
    linewidths: Any,
    facecolor: Any,
    facecolors: Any,
    edgecolor: Any,
    edgecolors: Any,
    norm: Any,
    vmin: Any,
    vmax: Any,
) -> tuple[Any, Any, Any, Any, float | None, float | None]:
    """Resolve Matplotlib-style curly-vector aliases onto canonical names."""
    if color is not None and c is not None:
        raise ValueError("Use only one of 'c' or 'color'")
    if linewidth is not None and linewidths is not None:
        raise ValueError("Use only one of 'linewidth' or 'linewidths'")
    if facecolor is not None and facecolors is not None:
        raise ValueError("Use only one of 'facecolor' or 'facecolors'")
    if edgecolor is not None and edgecolors is not None:
        raise ValueError("Use only one of 'edgecolor' or 'edgecolors'")
    if norm is not None and (vmin is not None or vmax is not None):
        raise ValueError("Use only one of 'norm' or 'vmin'/'vmax'")

    if color is None:
        color = c
    if linewidth is None:
        linewidth = linewidths
    if facecolor is None:
        facecolor = facecolors
    if edgecolor is None:
        edgecolor = edgecolors

    vmin_value = None if vmin is None else float(vmin)
    vmax_value = None if vmax is None else float(vmax)
    if vmin_value is not None and not np.isfinite(vmin_value):
        raise ValueError("vmin must be finite")
    if vmax_value is not None and not np.isfinite(vmax_value):
        raise ValueError("vmax must be finite")

    return color, linewidth, facecolor, edgecolor, vmin_value, vmax_value


def _resolve_scatter_aliases(
    c: Any,
    kwargs: dict[str, Any],
) -> tuple[Any, dict[str, Any]]:
    """Resolve singular compatibility aliases onto Matplotlib scatter kwargs."""
    resolved = dict(kwargs)

    if "color" in resolved:
        if c is not None:
            raise ValueError("Use only one of 'c' or 'color'")
        c = resolved.pop("color")

    if "linewidth" in resolved:
        if "linewidths" in resolved:
            raise ValueError("Use only one of 'linewidth' or 'linewidths'")
        resolved["linewidths"] = resolved.pop("linewidth")

    if "facecolor" in resolved:
        if "facecolors" in resolved:
            raise ValueError("Use only one of 'facecolor' or 'facecolors'")
        resolved["facecolors"] = resolved.pop("facecolor")

    if "edgecolor" in resolved:
        if "edgecolors" in resolved:
            raise ValueError("Use only one of 'edgecolor' or 'edgecolors'")
        resolved["edgecolors"] = resolved.pop("edgecolor")

    return c, resolved
