"""Shared axes and CRS-like detection helpers for Skyborn plots."""

from __future__ import annotations

from typing import Any


def _is_cartopy_crs_like(transform: Any) -> bool:
    """Return whether *transform* looks like a Cartopy CRS object."""
    return transform is not None and hasattr(transform, "_as_mpl_transform")


def _looks_like_axes(value: Any) -> bool:
    """Return whether *value* looks like a Matplotlib axes object."""
    return (
        hasattr(value, "transData")
        and hasattr(value, "transAxes")
        and (hasattr(value, "add_collection") or hasattr(value, "add_artist"))
    )
