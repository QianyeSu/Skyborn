"""
spharm - Spherical harmonic transforms for atmospheric/oceanic modeling

This module provides Python interfaces to NCAR's SPHEREPACK Fortran library
for spherical harmonic transforms on the sphere. It is particularly useful
for atmospheric and oceanic modeling applications.

Main classes:
    Spharmt: Main interface for spherical harmonic transforms

Example:
    >>> from skyborn.spharm import Spharmt
    >>> sht = Spharmt(nlon=144, nlat=73)
    >>> spec = sht.grdtospec(data)  # Grid to spectral transform
    >>> data_back = sht.spectogrd(spec)  # Spectral to grid transform
"""

from .spherical_harmonics import *
from .spherical_harmonics import __all__ as _spherical_harmonics_all

__all__ = [*_spherical_harmonics_all, "ReducedGaussianSpharmt"]


def __getattr__(name):
    if name == "ReducedGaussianSpharmt":
        from .reduced_gaussian import ReducedGaussianSpharmt

        return ReducedGaussianSpharmt
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(__all__)


__author__ = "Qianye Su"
__license__ = "BSD-3-Clause"
