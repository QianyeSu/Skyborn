"""Define a procedure for filling missing values."""

# Copyright (c) 2012-2016 Andrew Dawson
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
from __future__ import absolute_import, print_function

import warnings

import numpy as np

from ._gridfill_upgrade import poisson_fill_grids as _poisson_fill_grids


def _order_dims(grid, xpos, ypos):
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


def _prep_data(grid, xdim, ydim):
    origndim = grid.ndim
    grid, intorder = _order_dims(grid, xdim, ydim)
    intshape = grid.shape
    grid = grid.reshape(grid.shape[:2] + (int(np.prod(grid.shape[2:])),))
    info = dict(intshape=intshape, intorder=intorder, origndim=origndim)
    grid = grid.astype(np.float64)
    return grid, info


def _recover_data(grid, info):
    grid = grid.reshape(info["intshape"])
    rolldims = np.array(
        [info["intorder"].index(dim) for dim in range(info["origndim"] - 1, -1, -1)]
    )
    for i in range(rolldims.size):
        grid = np.rollaxis(grid, rolldims[i])
        rolldims = np.where(rolldims < rolldims[i], rolldims + 1, rolldims)
    return grid


def fill(
    grids,
    xdim,
    ydim,
    eps,
    relax=0.6,
    itermax=100,
    initzonal=False,
    initzonal_linear=False,
    cyclic=False,
    initial_value=0.0,
    verbose=False,
):
    """
    Fill missing values in grids with values derived by solving
    Poisson's equation using a relaxation scheme.

    **Arguments:**

    *grid*
        A masked array with missing values to fill.

    *xdim*, *ydim*
        The numbers of the dimensions in *grid* that represent the
        x-coordinate and y-coordinate respectively.

    *eps*
        Tolerance for determining the solution complete.

    **Keyword arguments:**

    *relax*
        Relaxation constant. Usually 0.45 <= *relax* <= 0.6. Defaults to
        0.6.

    *itermax*
        Maximum number of iterations of the relaxation scheme. Defaults
        to 100 iterations.

    *initzonal*
        If *False* missing values will be initialized to zero, if *True*
        missing values will be initialized to the zonal mean. Defaults
        to *False*.

    *initzonal_linear*
        If *False* missing values will be initialized to zero, if *True*
        missing values will be initialized to the zonal mean. Defaults
        to *False*.

    *cyclic*
        Set to *False* if the x-coordinate of the grid is not cyclic,
        set to *True* if it is cyclic. Defaults to *False*.

    *verbose*
        If *True* information about algorithm performance will be
        printed to stdout, if *False* nothing is printed. Defaults to
        *False*.

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
        grids,
        masks,
        relax,
        eps,
        itermax,
        1 if cyclic else 0,
        1 if initzonal else 0,
        1 if initzonal_linear else 0,
        initial_value,
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
    cube,
    eps,
    relax=0.6,
    itermax=100,
    initzonal=False,
    full_output=False,
    verbose=False,
    inplace=False,
):
    """Fill missing values in a cube with values derived by solving
    Poisson's equation using a relaxation scheme in the horizontal
    (i.e., x-y-) plane.

    **Arguments:**

    *cube*
        The iris.cube.Cube to be filled

    *eps*
        Tolerance for determining the solution complete.

    **Keyword arguments:**

    *relax*
        Relaxation constant. Usually 0.45 <= *relax* <= 0.6. Defaults to
        0.6.

    *itermax*
        Maximum number of iterations of the relaxation scheme. Defaults
        to 100 iterations.

    *initzonal*
        If *False* missing values will be initialized to zero, if *True*
        missing values will be initialized to the zonal mean. Defaults
        to *False*.

    *full_output*
        If *True* this function returns a tuple (filled_cube,
        not_converged), where not_converged is a boolean array
        indicating slices where the algorithm did not converge.
        Defaults to *False*.

    *verbose*
        If *True* information about algorithm performance will be
        printed to stdout, if *False* nothing is printed. Defaults to
        *False*.

    *inplace*
        If *True*, modify cube.data in-place; if *False*, return a
        copy. Defaults to *False*.

    """
    cyclic = cube.coord(axis="x").circular
    filled_data, cnv = fill(
        cube.data,
        cube.coord_dims(cube.coord(axis="x"))[0],
        cube.coord_dims(cube.coord(axis="y"))[0],
        eps=eps,
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
