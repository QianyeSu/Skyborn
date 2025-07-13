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

    import warnings

    warnings.warn(
        "spharm Fortran extensions not available. "
        "To build extensions, ensure you have meson, ninja, and gfortran installed, "
        "then reinstall skyborn.",
        ImportWarning,
    )

    class Spharmt:
        """
        Placeholder Spharmt class for environments without Fortran compiler.

        This is a documentation placeholder. The actual implementation requires
        a compiled Fortran extension (_spherepack) which is not available in
        this environment.

        To build the Fortran extensions:
        1. Install build dependencies: pip install meson ninja
        2. Install Fortran compiler: apt-get install gfortran (Ubuntu/Debian)
        3. Reinstall skyborn: pip install --force-reinstall skyborn
        """

        def __init__(self, *args, **kwargs):
            raise ImportError(
                "spharm module requires compiled Fortran extensions.\n\n"
                "To install required dependencies:\n"
                "  Linux/macOS: pip install meson ninja && apt-get install gfortran\n"
                "  Windows: Install MSYS2 and MINGW64 toolchain\n\n"
                "Then reinstall skyborn:\n"
                "  pip install --force-reinstall --no-binary=skyborn skyborn"
            )

    def regrid(*args, **kwargs):
        """Placeholder function for regrid."""
        raise ImportError("spharm module not available - Fortran extensions required")

    def gaussian_lats_wts(*args, **kwargs):
        """Placeholder function for gaussian_lats_wts."""
        raise ImportError("spharm module not available - Fortran extensions required")


__author__ = "Qianye Su"
__license__ = "BSD-3-Clause"
