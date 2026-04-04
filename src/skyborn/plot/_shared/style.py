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
