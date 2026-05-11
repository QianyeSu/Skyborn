"""Spherical harmonic vector wind analysis for Skyborn.

This module is based on the ajdawson/windspharm project (https://github.com/ajdawson/windspharm),
originally authored by Andrew Dawson. The current version is maintained by Qianye Su and is licensed under the BSD-3-Clause.

Main Classes:
    VectorWind: Enhanced interface for wind analysis on regular grids
    ReducedVectorWind: Packed reduced-Gaussian wind analysis

Example:
    >>> from skyborn.windspharm import VectorWind
    >>> import numpy as np
    >>>
    >>> # Create sample wind data
    >>> nlat, nlon = 73, 144
    >>> u = np.random.randn(nlat, nlon)
    >>> v = np.random.randn(nlat, nlon)
    >>>
    >>> # Initialize VectorWind with type hints and modern interface
    >>> vw = VectorWind(u, v, gridtype='gaussian')
    >>>
    >>> # Calculate various fields with improved documentation
    >>> vorticity = vw.vorticity()
    >>> divergence = vw.divergence()
    >>> psi, chi = vw.sfvp()

Notes:
    The NumPy-facing ``VectorWind`` interface expects latitude to run
    north-to-south. If your array is south-to-north, reorder it first with
    ``skyborn.windspharm.tools.order_latdim`` or use the xarray interface,
    which handles latitude reordering automatically.
"""

from __future__ import absolute_import

import importlib

from . import standard, tools, xarray

# Import main class for easier access
from .standard import VectorWind

__all__ = [
    "VectorWind",
    "ReducedVectorWind",
    "standard",
    "reduced",
    "tools",
    "xarray",
]


def __getattr__(name):
    if name == "ReducedVectorWind":
        from .reduced import ReducedVectorWind

        return ReducedVectorWind
    if name == "reduced":
        return importlib.import_module(f"{__name__}.reduced")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(__all__)


__author__ = "Qianye Su"
__license__ = "BSD-3-Clause"
