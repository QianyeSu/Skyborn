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

# Try to import the main module, with graceful fallback for environments
# where the Fortran extension cannot be compiled (e.g., Read the Docs)
try:
    from .spherical_harmonics import *

    _spharm_available = True
except ImportError as e:
    # Create placeholder classes/functions for documentation purposes
    _spharm_available = False

    class Spharmt:
        """
        Placeholder Spharmt class for environments without Fortran compiler.

        This is a documentation placeholder. The actual implementation requires
        a compiled Fortran extension (_spherepack) which is not available in
        this environment.
        """

        def __init__(self, *args, **kwargs):
            raise ImportError(
                "spharm module requires compiled Fortran extensions. "
                "Please ensure gfortran and proper build environment are available."
            )

    def regrid(*args, **kwargs):
        """Placeholder function for regrid."""
        raise ImportError("spharm module not available")

    def gaussian_lats_wts(*args, **kwargs):
        """Placeholder function for gaussian_lats_wts."""
        raise ImportError("spharm module not available")


__version__ = "1.0.9"
__author__ = "Jeff Whitaker"
__license__ = "BSD-3-Clause"
