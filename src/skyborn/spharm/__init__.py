"""
spharm - Spherical harmonic transforms for atmospheric/oceanic modeling

This module provides Python interfaces to NCAR's SPHEREPACK Fortran library
for spherical harmonic transforms on the sphere. It is particularly useful
for atmospheric and oceanic modeling applications.

Main classes:
    Spharmt: Main interface for spherical harmonic transforms
    ReducedGaussianSpharmt: Interface for packed reduced Gaussian transforms

Example:
    >>> from skyborn.spharm import Spharmt
    >>> sht = Spharmt(nlon=144, nlat=73)
    >>> spec = sht.grdtospec(data)  # Grid to spectral transform
    >>> data_back = sht.spectogrd(spec)  # Spectral to grid transform
"""

from .reduced_gaussian import ReducedGaussianSpharmt
from .spherical_harmonics import (
    Spharmt,
    SpheremackError,
    ValidationError,
    gaussian_lats_wts,
    getgeodesicpts,
    getspecindx,
    legendre,
    regrid,
    specintrp,
)

__all__ = [
    "Spharmt",
    "SpheremackError",
    "ValidationError",
    "regrid",
    "gaussian_lats_wts",
    "getspecindx",
    "getgeodesicpts",
    "legendre",
    "specintrp",
    "ReducedGaussianSpharmt",
]


__author__ = "Qianye Su"
__license__ = "BSD-3-Clause"
