"""
Grid filling procedures for missing value interpolation.

This module provides functions to fill missing values in gridded data using
iterative relaxation methods to solve Poisson's equation. It is particularly
useful for meteorological and oceanographic data where spatial interpolation
of missing values is needed while preserving physical constraints.

This module is based on the gridfill package by Andrew Dawson:
https://github.com/ajdawson/gridfill

Key Features:
    - Poisson equation solver for gap filling
    - Support for cyclic and non-cyclic boundaries
    - Configurable initialization (zeros or zonal mean)
    - Multi-dimensional array support
    - Integration with iris cubes

Mathematical Background:
    The algorithm solves the 2D Poisson equation:
    ∇²φ = 0
    where φ represents the field to be filled. The iterative relaxation
    scheme converges to a solution that smoothly interpolates missing values
    while preserving the boundary conditions from observed data.

Examples:
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
    >>> print(f"Convergence: {converged}")

Note:
    This implementation requires compiled Cython extensions for optimal
    performance. The core computational routines are implemented in
    `_gridfill.pyx` and compiled during package installation.
"""

from __future__ import absolute_import, print_function

import warnings
from typing import Tuple, Optional, Union, Any, Dict, TYPE_CHECKING

import numpy as np
import numpy.ma as ma

from ._gridfill import poisson_fill_grids as _poisson_fill_grids


def _order_dims(grid: np.ndarray, xpos: int, ypos: int) -> Tuple[np.ndarray, list]:
    outorder = list(range(grid.ndim))
    try:
        outorder.remove(xpos)
        outorder.remove(ypos)
    except ValueError:
        raise ValueError(
            "xdim and ydim must be the numbers of "
            "the array dimensions corresponding to the "
            "x-coordinate and y-coordinate respectively"
        )
    outorder = [ypos, xpos] + outorder
    grid = np.rollaxis(grid, xpos)
    if ypos < xpos:
        ypos += 1
    grid = np.rollaxis(grid, ypos)
    return grid, outorder


def _prep_data(grid: np.ndarray, xdim: int, ydim: int) -> Tuple[np.ndarray, dict]:
    origndim = grid.ndim
    grid, intorder = _order_dims(grid, xdim, ydim)
    intshape = grid.shape
    grid = grid.reshape(grid.shape[:2] + (int(np.prod(grid.shape[2:])),))
    info = dict(intshape=intshape, intorder=intorder, origndim=origndim)
    grid = grid.astype(np.float64)
    return grid, info


def _recover_data(grid: np.ndarray, info: dict) -> np.ndarray:
    grid = grid.reshape(info["intshape"])
    rolldims = np.array(
        [info["intorder"].index(dim) for dim in range(info["origndim"] - 1, -1, -1)]
    )
    for i in range(rolldims.size):
        grid = np.rollaxis(grid, rolldims[i])
        rolldims = np.where(rolldims < rolldims[i], rolldims + 1, rolldims)
    return grid


def fill(
    grids: np.ma.MaskedArray,
    xdim: int,
    ydim: int,
    eps: float,
    relax: float = 0.6,
    itermax: int = 100,
    initzonal: bool = False,
    cyclic: bool = False,
    verbose: bool = False,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fill missing values in grids using iterative Poisson equation solver.

    This function fills missing values in gridded data by solving Poisson's
    equation (∇²φ = 0) using a successive over-relaxation (SOR) iterative
    scheme. The algorithm converges to a solution that smoothly interpolates
    missing values while preserving boundary conditions from observed data.

    Parameters
    ----------
    grids : numpy.ma.MaskedArray
        A masked array containing the gridded data with missing values to fill.
        The mask indicates which values are missing (True) or observed (False).
        Must be at least 2-dimensional.
    xdim : int
        The axis number corresponding to the x-coordinate (typically longitude).
        For arrays with shape (lat, lon), this would be 1.
    ydim : int
        The axis number corresponding to the y-coordinate (typically latitude).
        For arrays with shape (lat, lon), this would be 0.
    eps : float
        Convergence tolerance. The iteration stops when the maximum change
        between successive iterations falls below this threshold.
    relax : float, optional
        Relaxation constant for the SOR iteration scheme. Must be in the range
        (0, 2) for convergence, with values between 0.45 and 0.6 typically
        providing good performance. Default is 0.6.
    itermax : int, optional
        Maximum number of iterations allowed before stopping. Default is 100.
    initzonal : bool, optional
        Initialization method for missing values:
        - False: Initialize missing values to zero
        - True: Initialize missing values to the zonal (x-direction) mean
        Default is False.
    cyclic : bool, optional
        Boundary condition for the x-coordinate:
        - False: Non-cyclic boundary (e.g., for regional grids)
        - True: Cyclic boundary (e.g., for global longitude grids)
        Default is False.
    verbose : bool, optional
        If True, print convergence information including iteration count
        and residual values. Default is False.

    Returns
    -------
    filled_grids : numpy.ndarray
        The input grids with missing values filled by interpolation.
        Same shape as input but as a regular numpy array.
    converged : numpy.ndarray
        Boolean array indicating convergence status for each grid.
        True if the algorithm converged within the specified tolerance,
        False otherwise.

    Raises
    ------
    TypeError
        If the input is not a masked array.
    ValueError
        If xdim or ydim are invalid dimension numbers.

    Examples
    --------
    >>> import numpy as np
    >>> import numpy.ma as ma
    >>> from skyborn.gridfill import fill
    >>>
    >>> # Create test data with missing values
    >>> data = np.random.rand(10, 20)
    >>> mask = np.zeros_like(data, dtype=bool)
    >>> mask[4:6, 8:12] = True  # Create a rectangular gap
    >>> masked_data = ma.array(data, mask=mask)
    >>>
    >>> # Fill missing values
    >>> filled_data, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-4)
    >>> print(f"Convergence: {converged[0]}")

    Notes
    -----
    The algorithm uses the finite difference approximation to Poisson's equation:

    .. math::
        \\frac{\\partial^2 \\phi}{\\partial x^2} + \\frac{\\partial^2 \\phi}{\\partial y^2} = 0

    The successive over-relaxation scheme updates each grid point using:

    .. math::
        \\phi_{i,j}^{(n+1)} = \\phi_{i,j}^{(n)} + \\omega \\cdot r_{i,j}^{(n)}

    where ω is the relaxation parameter and r is the residual.

    Performance Tips:
    - For meteorological data, consider initzonal=True if there are strong
      zonal patterns
    - Use cyclic=True for global longitude grids spanning 0-360°
    - Adjust eps based on your precision requirements vs. computation time
    - The relaxation parameter affects convergence speed; values around 0.6
      typically work well

    References
    ----------
    This implementation is based on the algorithm described in numerical
    analysis textbooks for solving elliptic PDEs, particularly the treatment
    in Press et al., "Numerical Recipes".
    """
    # re-shape to 3-D leaving the grid dimensions at the front:
    grids, info = _prep_data(grids, xdim, ydim)
    # fill missing values:
    fill_value = 1.0e20
    try:
        masks = grids.mask.astype(np.int32)
        grids = grids.filled(fill_value=fill_value)
    except AttributeError:
        raise TypeError("grids must be a masked array")
    # Call the computation subroutine:
    niter, resmax = _poisson_fill_grids(
        grids, masks, relax, eps, itermax, 1 if cyclic else 0, 1 if initzonal else 0
    )
    grids = _recover_data(grids, info)
    converged = np.logical_not(resmax > eps)
    # optional performance information:
    if verbose:
        for i, c in enumerate(converged):
            if c:
                converged_string = "converged"
            else:
                converged_string = "did not converge"
            print(
                "[{:d}] relaxation {:s} ({:d} iterations "
                "with maximum residual {:.3e})".format(
                    i, converged_string, int(niter[i]), resmax[i]
                )
            )
    return grids, converged


def fill_cube(
    cube: Any,
    xdim: Optional[int] = None,
    ydim: Optional[int] = None,
    eps: float = 1e-4,
    relax: float = 0.6,
    itermax: int = 100,
    initzonal: bool = False,
    cyclic: Optional[bool] = None,
    verbose: bool = False,
    inplace: bool = False,
    full_output: bool = False,
) -> Union[Any, Tuple[Any, np.ndarray]]:
    """
    Fill missing values in an iris cube using Poisson equation solver.

    This function fills missing values in an iris.cube.Cube by solving
    Poisson's equation in the horizontal (x-y) plane. It automatically
    detects coordinate dimensions and can handle cyclic boundaries.

    Parameters
    ----------
    cube : iris.cube.Cube
        The iris cube containing gridded data with missing values to fill.
        Must have at least 2 spatial dimensions.
    xdim : int, optional
        The axis number corresponding to the x-coordinate (longitude).
        If None, automatically detected from cube coordinates with axis='x'.
    ydim : int, optional
        The axis number corresponding to the y-coordinate (latitude).
        If None, automatically detected from cube coordinates with axis='y'.
    eps : float, optional
        Convergence tolerance for the iterative solver. Default is 1e-4.
    relax : float, optional
        Relaxation constant for the SOR iteration. Should be between
        0.45 and 0.6 for optimal performance. Default is 0.6.
    itermax : int, optional
        Maximum number of iterations before stopping. Default is 100.
    initzonal : bool, optional
        If True, initialize missing values with zonal (x-direction) mean.
        If False, initialize with zeros. Default is False.
    cyclic : bool, optional
        If True, treat x-coordinate as cyclic (e.g., longitude wrapping).
        If None, automatically detected from cube coordinate properties.
        Default is None (auto-detect).
    verbose : bool, optional
        If True, print convergence information. Default is False.
    inplace : bool, optional
        If True, modify the cube's data in-place. If False, return a copy.
        Default is False.
    full_output : bool, optional
        If True, return tuple (filled_cube, convergence_info).
        If False, return only the filled cube. Default is False.

    Returns
    -------
    filled_cube : iris.cube.Cube
        The cube with missing values filled by interpolation.
    convergence_info : numpy.ndarray, optional
        Only returned if full_output=True. Boolean array indicating
        which slices failed to converge (True = did not converge).

    Raises
    ------
    ValueError
        If the cube lacks required spatial coordinates.
    UserWarning
        If the algorithm fails to converge on some slices.

    Examples
    --------
    >>> import iris
    >>> from skyborn.gridfill import fill_cube
    >>>
    >>> # Load a cube with missing data
    >>> cube = iris.load_cube('temperature_data.nc')
    >>>
    >>> # Fill missing values
    >>> filled_cube = fill_cube(cube, eps=1e-4, verbose=True)
    >>>
    >>> # Get convergence information
    >>> filled_cube, not_converged = fill_cube(
    ...     cube, full_output=True, eps=1e-5
    ... )
    >>> if not_converged.any():
    ...     print(f"Failed to converge on {not_converged.sum()} slices")

    Notes
    -----
    This function is a high-level wrapper around the `fill()` function
    that integrates with the iris data model. It automatically:

    - Detects spatial coordinate dimensions
    - Determines if longitude is cyclic
    - Handles multi-dimensional cubes by applying filling to each 2D slice
    - Preserves cube metadata and coordinate information

    For cubes with more than 2 spatial dimensions, the algorithm is applied
    independently to each horizontal slice (x-y plane).

    The automatic coordinate detection looks for coordinates with
    axis='x' and axis='y' attributes. For cubes without these attributes,
    you must specify xdim and ydim explicitly.
    """
    # Auto-detect coordinate dimensions if not specified
    if xdim is None:
        try:
            xdim = cube.coord_dims(cube.coord(axis="x"))[0]
        except Exception:
            raise ValueError(
                "Could not auto-detect x-dimension. Please specify xdim parameter."
            )

    if ydim is None:
        try:
            ydim = cube.coord_dims(cube.coord(axis="y"))[0]
        except Exception:
            raise ValueError(
                "Could not auto-detect y-dimension. Please specify ydim parameter."
            )

    # Auto-detect cyclic boundary if not specified
    if cyclic is None:
        try:
            cyclic = cube.coord(axis="x").circular
        except Exception:
            cyclic = False

    filled_data, cnv = fill(
        cube.data,
        xdim,
        ydim,
        eps=eps,
        relax=relax,
        itermax=itermax,
        initzonal=initzonal,
        cyclic=cyclic,
        verbose=verbose,
    )

    not_converged = np.logical_not(cnv)
    if np.any(not_converged):
        warnings.warn(
            "gridfill did not converge on {} out of {} "
            "slices".format(not_converged.sum(), not_converged.size)
        )
    if inplace:
        cube.data = filled_data
        retcube = cube
    else:
        retcube = cube.copy(data=filled_data)
    if full_output:
        return retcube, not_converged
    else:
        return retcube
