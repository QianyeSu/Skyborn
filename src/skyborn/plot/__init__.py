"""Compatibility namespace re-exporting :mod:`skybornplot` (``skyborn-plot``)."""

from __future__ import annotations

import importlib
import sys
from typing import Any

import skybornplot as _skybornplot

from .._package_alias import register_subpackage_alias

register_subpackage_alias("skyborn.plot", "skybornplot")

for _submod in ("contour", "plotting", "scatter", "vector"):
    sys.modules[f"{__name__}.{_submod}"] = importlib.import_module(
        f"skybornplot.{_submod}"
    )

__all__ = [name for name in _skybornplot.__all__ if name != "__version__"]

for _name in __all__:
    globals()[_name] = getattr(_skybornplot, _name)


def __getattr__(name: str) -> Any:
    try:
        return getattr(_skybornplot, name)
    except AttributeError:
        try:
            mod = importlib.import_module(f"skyborn.plot.{name}")
            globals()[name] = mod
            return mod
        except ImportError as exc:
            raise AttributeError(
                f"module {__name__!r} has no attribute {name!r}"
            ) from exc


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
