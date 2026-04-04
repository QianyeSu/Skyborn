"""Small shared type aliases for Skyborn plotting internals."""

from __future__ import annotations

from collections.abc import Hashable

CoordinateShape = tuple[int, ...]
TargetDims = tuple[Hashable, ...] | None
