"""Small shared type aliases for Skyborn plotting internals."""

from __future__ import annotations

from collections.abc import Hashable
from typing import Optional, Tuple

CoordinateShape = tuple[int, ...]
TargetDims = Optional[Tuple[Hashable, ...]]
