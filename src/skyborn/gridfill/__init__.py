"""
GridFill Module for Skyborn
============================

Fill missing values in grids using iterative relaxation.

This module provides functions to fill missing values in gridded data using
iterative relaxation methods to solve Poisson's equation. It is particularly
useful for meteorological and oceanographic data where spatial interpolation
of missing values is needed while preserving physical constraints.

This module is based on the gridfill package by Andrew Dawson:
https://github.com/ajdawson/gridfill

Mathematical Background:
    The algorithm solves the 2D Poisson equation:
    ∇²φ = 0
    where φ represents the field to be filled. The iterative relaxation
    scheme converges to a solution that smoothly interpolates missing values
    while preserving the boundary conditions from observed data.

Main Functions:
    fill: Fill missing values in numpy arrays or masked arrays
    fill_cube: Fill missing values in iris cubes (if iris is available)

Interfaces:
    - Standard numpy/masked array interface via fill() function
    - Iris cube interface via fill_cube() function
    - Modern xarray interface via skyborn.gridfill.xarray submodule

Examples:
    Basic usage with numpy masked arrays:

    >>> import numpy as np
    >>> import numpy.ma as ma
    >>> from skyborn.gridfill import fill
    >>>
    >>> # Create test data with missing values
    >>> data = np.random.rand(50, 100)
    >>> mask = np.zeros_like(data, dtype=bool)
    >>> mask[20:30, 40:60] = True  # Create a gap
    >>> masked_data = ma.array(data, mask=mask)
    >>>
    >>> # Fill the missing values
    >>> filled_data, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-4)
    >>> print(f"Convergence: {converged[0]}")

    Using with iris cubes (if iris is installed):

    >>> import iris
    >>> from skyborn.gridfill import fill_cube
    >>>
    >>> # Load a cube with missing data
    >>> cube = iris.load_cube('temperature_data.nc')
    >>>
    >>> # Fill missing values
    >>> filled_cube = fill_cube(cube, eps=1e-4, verbose=True)

    Modern xarray interface with automatic coordinate detection:

    >>> import xarray as xr
    >>> from skyborn.gridfill.xarray import fill
    >>>
    >>> # Load data with missing values
    >>> data = xr.open_dataarray('sst_with_gaps.nc')
    >>>
    >>> # Fill missing values preserving metadata
    >>> filled = fill(data, eps=1e-4)
    >>> print(f"Attributes preserved: {filled.attrs == data.attrs}")

Notes:
    This implementation requires compiled Cython extensions for optimal
    performance. The core computational routines are implemented in
    `_gridfill.pyx` and compiled during package installation.
"""

from __future__ import absolute_import

from .gridfill import fill, fill_cube

# Try to import xarray interface if xarray is available
try:
    from . import xarray

    _has_xarray = True
except ImportError:
    _has_xarray = False

try:
    from skyborn import __version__
except ImportError:
    __version__ = "1.0.0"

# Define the objects to be imported by imports of the form:
#   from skyborn.gridfill import *
__all__ = ["fill", "fill_cube"]

# Add xarray to __all__ if available
if _has_xarray:
    __all__.append("xarray")

# Note: xarray interface is available as a submodule:
#   from skyborn.gridfill.xarray import fill
