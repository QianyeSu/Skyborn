"""Compatibility namespace re-exporting :mod:`skybornspharm` (``skyborn-spharm``)."""

from __future__ import annotations

import importlib
from typing import Any

import skybornspharm as _skybornspharm

from .._package_alias import register_subpackage_alias

register_subpackage_alias("skyborn.spharm", "skybornspharm")

__author__ = getattr(_skybornspharm, "__author__", "Qianye Su")
__license__ = getattr(_skybornspharm, "__license__", "BSD-3-Clause")
__all__ = [name for name in _skybornspharm.__all__ if name != "__version__"]

for _name in __all__:
    globals()[_name] = getattr(_skybornspharm, _name)


def __getattr__(name: str) -> Any:
    try:
        return getattr(_skybornspharm, name)
    except AttributeError:
        try:
            mod = importlib.import_module(f"skyborn.spharm.{name}")
            globals()[name] = mod
            return mod
        except ImportError as exc:
            raise AttributeError(
                f"module {__name__!r} has no attribute {name!r}"
            ) from exc


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
