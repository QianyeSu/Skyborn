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
from typing import Tuple, Optional, Union, Any, Dict

import numpy as np
import numpy.ma as ma

from ._gridfill import poisson_fill_grids as _poisson_fill_grids

# Type alias for iris cube (may not be available)
try:
    import iris

    IrisCube = iris.cube.Cube
except (ImportError, RuntimeError):
    IrisCube = Any


def _order_dims(grid: np.ndarray, xpos: int, ypos: int) -> Tuple[np.ndarray, list]:
    """
    Reorder array dimensions to put spatial dimensions first.

    This internal function reorganizes the input array so that the y-dimension
    comes first, followed by the x-dimension, then all other dimensions. This
    standardized ordering is required by the Cython computational routines.

    Parameters
    ----------
    grid : np.ndarray
        Input array to reorder
    xpos : int
        Position of x-coordinate dimension in original array
    ypos : int
        Position of y-coordinate dimension in original array

    Returns
    -------
    grid : np.ndarray
        Reordered array with shape (ny, nx, ...)
    outorder : list
        List tracking the dimension reordering for later restoration

    Raises
    ------
    ValueError
        If xpos or ypos are invalid dimension indices

    Examples
    --------
    >>> data = np.random.rand(10, 20, 30)  # (time, lat, lon)
    >>> reordered, order = _order_dims(data, xpos=2, ypos=1)
    >>> print(reordered.shape)  # (20, 30, 10) -> (lat, lon, time)
    """
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


def _prep_data(grid: np.ndarray, xdim: int, ydim: int) -> Tuple[np.ndarray, Dict]:
    """
    Prepare input data for the Cython filling algorithm.

    This function performs several preprocessing steps:
    1. Reorders dimensions to put spatial coordinates first
    2. Reshapes to 3D format (ny, nx, nother) for efficient processing
    3. Converts to float64 for numerical precision
    4. Stores metadata needed for result reconstruction

    Parameters
    ----------
    grid : np.ndarray
        Input array to prepare
    xdim : int
        X-coordinate dimension index
    ydim : int
        Y-coordinate dimension index

    Returns
    -------
    grid : np.ndarray
        Preprocessed array with shape (ny, nx, nother) as float64
    info : dict
        Metadata dictionary containing:
        - 'intshape': intermediate shape after reordering
        - 'intorder': dimension reordering information
        - 'origndim': original number of dimensions

    Notes
    -----
    The reshaping flattens all non-spatial dimensions into a single dimension
    to enable efficient batch processing of multiple 2D grids.
    """
    origndim = grid.ndim
    grid, intorder = _order_dims(grid, xdim, ydim)
    intshape = grid.shape
    grid = grid.reshape(grid.shape[:2] + (int(np.prod(grid.shape[2:])),))
    info = dict(intshape=intshape, intorder=intorder, origndim=origndim)
    grid = grid.astype(np.float64)
    return grid, info


def _recover_data(grid: np.ndarray, info: Dict) -> np.ndarray:
    """
    Restore original array structure after processing.

    This function reverses the preprocessing steps applied by `_prep_data`:
    1. Reshapes from 3D back to original dimensional structure
    2. Reorders dimensions back to original layout

    Parameters
    ----------
    grid : np.ndarray
        Processed array with shape (ny, nx, nother)
    info : dict
        Metadata dictionary from `_prep_data` containing restoration info

    Returns
    -------
    np.ndarray
        Array restored to original dimensional structure and ordering

    Notes
    -----
    This function must be called with the same info dictionary that was
    generated by the corresponding `_prep_data` call to ensure correct
    reconstruction.
    """
    grid = grid.reshape(info["intshape"])
    rolldims = np.array(
        [info["intorder"].index(dim) for dim in range(info["origndim"] - 1, -1, -1)]
    )
    for i in range(rolldims.size):
        grid = np.rollaxis(grid, rolldims[i])
        rolldims = np.where(rolldims < rolldims[i], rolldims + 1, rolldims)
    return grid


def fill(
    grids: ma.MaskedArray,
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
    Fill missing values in grids using Poisson equation solver.

    This function fills missing values in masked arrays by solving Poisson's
    equation (∇²φ = 0) using an iterative relaxation scheme. The method provides
    smooth interpolation of gaps while preserving boundary conditions from
    observed data.

    Parameters
    ----------
    grids : numpy.ma.MaskedArray
        Masked array containing data with missing values to fill. Missing
        values are indicated by the mask.
    xdim : int
        Index of the dimension representing the x-coordinate (longitude)
    ydim : int
        Index of the dimension representing the y-coordinate (latitude)
    eps : float
        Convergence tolerance. Iteration stops when the maximum residual
        falls below this threshold.
    relax : float, optional
        Relaxation parameter for the iterative scheme. Must be in range
        (0, 1). Values between 0.45-0.6 typically work well. Default: 0.6
    itermax : int, optional
        Maximum number of iterations. Default: 100
    initzonal : bool, optional
        Initialization method for missing values:
        - False: Initialize with zeros
        - True: Initialize with zonal (x-direction) mean
        Default: False
    cyclic : bool, optional
        Whether the x-coordinate is cyclic (e.g., longitude wrapping).
        Default: False
    verbose : bool, optional
        Print convergence information for each grid. Default: False

    Returns
    -------
    filled_grids : numpy.ndarray
        Array with missing values filled, same shape as input
    converged : numpy.ndarray
        Boolean array indicating convergence status for each 2D grid slice

    Raises
    ------
    TypeError
        If input is not a masked array
    ValueError
        If xdim or ydim are invalid dimension indices

    Notes
    -----
    The algorithm solves:
    ∇²φ = (∂²φ/∂x²) + (∂²φ/∂y²) = 0

    using a finite difference relaxation scheme:
    φ[i,j]^(k+1) = φ[i,j]^k + relax * (residual[i,j] / 4)

    where the residual is computed using a 5-point stencil.

    For best performance, use `fill_cube()` when working with iris cubes
    as it automatically handles coordinate detection and cyclic boundaries.

    Examples
    --------
    Fill a simple 2D grid:

    >>> import numpy as np
    >>> import numpy.ma as ma
    >>> from skyborn.gridfill import fill
    >>>
    >>> # Create test data with gaps
    >>> data = np.random.rand(50, 100)
    >>> mask = np.zeros_like(data, dtype=bool)
    >>> mask[20:30, 40:60] = True
    >>> masked_data = ma.array(data, mask=mask)
    >>>
    >>> # Fill missing values
    >>> filled, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-4)
    >>> print(f"Converged: {converged[0]}")

    Fill multiple time steps:

    >>> # 3D data: (time, lat, lon)
    >>> data_3d = np.random.rand(12, 50, 100)
    >>> mask_3d = np.zeros_like(data_3d, dtype=bool)
    >>> mask_3d[:, 20:30, 40:60] = True
    >>> masked_3d = ma.array(data_3d, mask=mask_3d)
    >>>
    >>> filled_3d, converged_3d = fill(masked_3d, xdim=2, ydim=1, eps=1e-4)
    >>> print(f"Convergence rate: {converged_3d.mean():.1%}")

    See Also
    --------
    fill_cube : Fill iris cubes with automatic coordinate detection
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
    cube: IrisCube,
    eps: float,
    relax: float = 0.6,
    itermax: int = 100,
    initzonal: bool = False,
    full_output: bool = False,
    verbose: bool = False,
    inplace: bool = False,
) -> Union[IrisCube, Tuple[IrisCube, np.ndarray]]:
    """
    Fill missing values in an iris cube using Poisson equation solver.

    This is a convenience function that wraps `fill()` for iris cubes,
    automatically detecting spatial coordinates and handling cyclic boundaries.
    It preserves all cube metadata while filling gaps in the data.

    Parameters
    ----------
    cube : iris.cube.Cube
        Input cube containing masked data to fill
    eps : float
        Convergence tolerance for the iterative solver
    relax : float, optional
        Relaxation parameter (0 < relax < 1). Default: 0.6
    itermax : int, optional
        Maximum iterations. Default: 100
    initzonal : bool, optional
        Initialize missing values with zonal mean instead of zeros. Default: False
    full_output : bool, optional
        Return convergence information along with filled cube. Default: False
    verbose : bool, optional
        Print convergence diagnostics. Default: False
    inplace : bool, optional
        Modify input cube in-place instead of creating copy. Default: False

    Returns
    -------
    filled_cube : iris.cube.Cube
        Cube with missing values filled
    not_converged : numpy.ndarray, optional
        Boolean array indicating non-converged slices (only if full_output=True)

    Raises
    ------
    CoordinateNotFoundError
        If cube lacks recognizable x or y coordinates

    Warnings
    --------
    Issues warning if algorithm fails to converge on any slices

    Notes
    -----
    - Automatically detects longitude/latitude coordinates using iris conventions
    - Handles cyclic longitude boundaries automatically
    - Preserves all cube attributes, coordinates, and metadata
    - For global grids, cyclic boundary conditions are applied automatically

    Examples
    --------
    Fill missing values in a temperature cube:

    >>> import iris
    >>> from skyborn.gridfill import fill_cube
    >>>
    >>> # Load cube with missing data
    >>> cube = iris.load_cube('temperature_with_gaps.nc')
    >>>
    >>> # Fill gaps with tight convergence
    >>> filled_cube = fill_cube(cube, eps=1e-5, verbose=True)
    >>>
    >>> # Get convergence information
    >>> filled_cube, not_converged = fill_cube(
    ...     cube, eps=1e-4, full_output=True
    ... )
    >>> print(f"Failed slices: {not_converged.sum()}")

    See Also
    --------
    fill : Lower-level function for numpy arrays
    """
    # Detect if longitude coordinate is cyclic
    cyclic = cube.coord(axis="x").circular

    # Get coordinate dimension indices
    x_dim = cube.coord_dims(cube.coord(axis="x"))[0]
    y_dim = cube.coord_dims(cube.coord(axis="y"))[0]

    # Fill the data
    filled_data, cnv = fill(
        cube.data,
        x_dim,
        y_dim,
        eps=eps,
        relax=relax,
        itermax=itermax,
        initzonal=initzonal,
        cyclic=cyclic,
        verbose=verbose,
    )

    # Check convergence and issue warnings
    not_converged = np.logical_not(cnv)
    if np.any(not_converged):
        warnings.warn(
            "gridfill did not converge on {} out of {} "
            "slices".format(not_converged.sum(), not_converged.size)
        )

    # Create output cube
    if inplace:
        cube.data = filled_data
        retcube = cube
    else:
        retcube = cube.copy(data=filled_data)

    if full_output:
        return retcube, not_converged
    else:
        return retcube
