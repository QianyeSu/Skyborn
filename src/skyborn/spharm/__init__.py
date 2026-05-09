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
        """
        Regrid data using spectral interpolation with optional smoothing and truncation.

        .. note::
           This function requires compiled Fortran extensions which are not available
           in this environment.

        Args:
            grdin: Input grid Spharmt instance
            grdout: Output grid Spharmt instance
            datagrid: Data on input grid
            ntrunc: Optional spectral truncation limit
            smooth: Optional smoothing factors

        Returns:
            Interpolated data on output grid

        Raises:
            ImportError: If Fortran extensions are not available
        """
        raise ImportError("spharm module not available - Fortran extensions required")

    def gaussian_lats_wts(*args, **kwargs):
        """
        Compute Gaussian latitudes and quadrature weights.

        .. note::
           This function requires compiled Fortran extensions which are not available
           in this environment.

        Args:
            nlat: Number of Gaussian latitudes desired

        Returns:
            Tuple of (latitudes_in_degrees, quadrature_weights)

        Raises:
            ImportError: If Fortran extensions are not available
        """
        raise ImportError("spharm module not available - Fortran extensions required")

    def getspecindx(*args, **kwargs):
        """
        Compute indices of zonal wavenumber and degree for spherical harmonic coefficients.

        .. note::
           This function requires compiled Fortran extensions which are not available
           in this environment.

        Args:
            ntrunc: Spherical harmonic triangular truncation limit

        Returns:
            Tuple of (zonal_wavenumber_indices, degree_indices)

        Raises:
            ImportError: If Fortran extensions are not available
        """
        raise ImportError("spharm module not available - Fortran extensions required")

    def getgeodesicpts(*args, **kwargs):
        """
        Compute lat/lon values for icosahedral geodesic points.

        .. note::
           This function requires compiled Fortran extensions which are not available
           in this environment.

        Args:
            m: Number of points on edge of single geodesic triangle

        Returns:
            Tuple of (latitudes, longitudes) in degrees

        Raises:
            ImportError: If Fortran extensions are not available
        """
        raise ImportError("spharm module not available - Fortran extensions required")

    def legendre(*args, **kwargs):
        """
        Calculate associated Legendre functions for triangular truncation.

        .. note::
           This function requires compiled Fortran extensions which are not available
           in this environment.

        Args:
            lat: Latitude in degrees
            ntrunc: Triangular truncation limit

        Returns:
            Associated Legendre functions array

        Raises:
            ImportError: If Fortran extensions are not available
        """
        raise ImportError("spharm module not available - Fortran extensions required")

    def specintrp(*args, **kwargs):
        """
        Spectral interpolation to arbitrary point on sphere.

        .. note::
           This function requires compiled Fortran extensions which are not available
           in this environment.

        Args:
            lon: Longitude in degrees
            dataspec: Spectral coefficients
            legfuncs: Associated Legendre functions

        Returns:
            Interpolated value

        Raises:
            ImportError: If Fortran extensions are not available
        """
        raise ImportError("spharm module not available - Fortran extensions required")


__author__ = "Qianye Su"
__license__ = "BSD-3-Clause"

"""
Note:
    This __init__.py is designed to enable API cross-references to work correctly in Read the Docs.
"""
