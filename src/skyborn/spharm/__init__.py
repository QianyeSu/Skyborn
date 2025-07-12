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

# Import from the main spherical harmonics module (renamed to avoid conflicts)
from .spherical_harmonics import *

__version__ = "1.0.9"
__author__ = "Jeff Whitaker"
__license__ = "BSD-3-Clause"
