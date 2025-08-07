#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False
#cython: cdivision=True, nonecheck=False, overflowcheck=False
#cython: profile=False, linetrace=False, binding=False
"""
Optimized Cython implementation for grid filling using iterative relaxation.

This module is based on the gridfill package by Andrew Dawson:
https://github.com/ajdawson/gridfill

Performance optimizations include:
- Disabled bounds checking and wraparound
- Enabled C-style division
- Disabled overflow checking for maximum speed
- Used nogil contexts where possible
"""

from libc.math cimport fabs, fmax
import numpy as np
cimport numpy as np


# Type definitions for numpy array types.
ctypedef np.float64_t FLOAT64_t
ctypedef np.uint32_t UINT32_t


cdef int int_sum(int[:, :] a, unsigned int ny, unsigned int nx) noexcept nogil:
    """
    Compute the sum of the elements in a 2-dimensional integer array.

    """
    cdef unsigned int i, j
    cdef int s
    s = 0
    for i in range(ny):
        for j in range(nx):
            s += a[i, j]
    return s


cdef void latitude_indices(unsigned int index, unsigned int nlat,
                           unsigned int *im1, unsigned int *ip1) noexcept nogil:
    """
    The indices of the neighbouring points for a given index in the
    latitude direction, taking into account grid edges.

    **Arguments:**

    * index [unsigned int]
        An index into the latitude dimension.

    * nlat [unsigned int]
        The size of the latitude dimension.

    * im1 [unsigned int *]
        The resulting index - 1 value.

    * ip1 [unsigned int *]
        The resulting index + 1 value.

    """
    im1[0] = index - 1
    ip1[0] = index + 1
    if index == 0:
        im1[0] = 1
    if index == nlat - 1:
        ip1[0] = nlat - 2


cdef void longitude_indices(unsigned int index, unsigned int nlon, int cyclic,
                            unsigned int *jm1, unsigned int *jp1) noexcept nogil:
    """
    The indices of the neighbouring points for a given index in the
    longitude direction, taking into account grid edges and cyclicity.

    **Arguments:**

    * index [unsigned int]
        An index into the longitude dimension.

    * nlon [unsigned int]
        The size of the longitude dimension.

    * cyclic [int]
        If `0` the input grid is assumed not to be cyclic. If non-zero
        the input grid is assumed to be cyclic in its second dimension.

    * jm1 [unsigned int *]
        The resulting index - 1 value.

    * jp1 [unsigned int *]
        The resulting index + 1 value.

    """
    jm1[0] = index - 1
    jp1[0] = index + 1
    if index == 0:
        if not cyclic:
            jm1[0] = 1
        else:
            jm1[0] = nlon - 1
    if index == nlon - 1:
        if not cyclic:
            jp1[0] = nlon - 2
        else:
            jp1[0] = 0


cdef void initialize_missing(double[:, :] grid,
                             int[:, :] mask,
                             unsigned int nlat,
                             unsigned int nlon,
                             int initialize_zonal) noexcept nogil:
    """
    Initialize the missing values in a grid in-place.

    **Arguments:**

    * grid [double[:, :]]
        A 2-dimensional array` containing a grid of values. The grid
        dimensions are [y, x]. The missing values will be modified
        in-place.

    * mask [int[:. :]]
        A 2-dimensional array the same shape as *grid* that contains the
        grid mask. Valid grid points are indicated with a zero value,
        and missing values are indicated by a non-zero value.

    * nlat [unsigned int]
        The size of the latitude (first) grid dimension.

    * nlon [unsigned int]
        The size of the longitude (second) grid dimension.

    * initialize_zonal [int]
        If non-zero, take the zonal mean as the initial guess for missing
        values. If zero, use the value `0` as the initial guess for all
        missing values.

    """
    cdef unsigned int i, j, n
    cdef double zonal_mean, initial_value
    for i in range(nlat):
        if initialize_zonal:
            n = 0
            zonal_mean = 0
            for j in range(nlon):
                if not mask[i, j]:
                    n += 1
                    zonal_mean += grid[i, j]
            if n > 0:
                zonal_mean /= n
            initial_value = zonal_mean
        else:
            initial_value = 0
        for j in range(nlon):
            if mask[i, j]:
                grid[i, j] = initial_value


cdef void poisson_fill(double[:, :] grid,
                       int[:, :] mask,
                       unsigned int nlat,
                       unsigned int nlon,
                       double relaxc,
                       double tolerance,
                       unsigned int itermax,
                       int cyclic,
                       int initialize_zonal,
                       unsigned int *numiter,
                       double *resmax) noexcept nogil:
    """
    Fill missing values in a grid by iterative relaxation.

    **Arguments:**

    * grid [double[:, :]]
        A 2-dimensional array containing a grid of values. The grid
        dimensions are [y, x].

    * mask [int[:, :]]
        A 2-dimensional array the same shape as *grid* that contains the
        grid mask. Valid grid points are indicated with a zero value,
        and missing values are indicated by a non-zero value.

    * nlat [unsigned int]
        The size of the latitude (first) grid dimension.

    * nlon [unsigned int]
        The size of the longitude (second) grid dimension.

    * relaxc [double]
        The relaxation constant, typically 0.45 <= *relaxc* <= 0.6.

    * tolerance [double]
        Numerical tolerance for convergence of the relaxation scheme.
        This value is data dependent.

    * itermax [unsigned int]
        The maximum number of iterations allowed for the relaxation
        scheme.

    * cyclic [int]
        If `0` the input grid is assumed not to be cyclic. If non-zero
        the input grid is assumed to be cyclic in its second dimension.

    * initialize_zonal [int]
        If non-zero, take the zonal mean as the initial guess for missing
        values. If zero, use the value `0` as the initial guess for all
        missing values.

    * numiter [unsigned int *]
        The number of iterations used to fill the grid.

    * resmax [double *]
        The maximum residual value at the end of the iteration. If this
        value is larger than the specified *tolerance* then the iteration
        did not converge.

    """
    cdef unsigned int _numiter
    cdef unsigned int i, j, im1, ip1, jm1, jp1
    cdef double _resmax, dp25, residual
    # Exit early if there are no missing values in the grid.
    if int_sum(mask, nlat, nlon) == 0:
        numiter[0] = 0
        resmax[0] = 0
        return
    # Set initial values for all missing values.
    initialize_missing(grid, mask, nlat, nlon, initialize_zonal)
    dp25 = 0.25
    _numiter = 0
    _resmax = 0
    while _numiter < itermax:
        _resmax = 0
        _numiter += 1
        for i in range(nlat):
            latitude_indices(i, nlat, &im1, &ip1)
            for j in range(nlon):
                if mask[i, j]:
                    longitude_indices(j, nlon, cyclic, &jm1, &jp1)
                    residual = dp25 * (grid[im1, j] + grid[ip1, j] +
                                       grid[i, jm1] + grid[i, jp1]) - grid[i, j]
                    residual *= relaxc
                    grid[i, j] += residual
                    _resmax = fmax(fabs(residual), _resmax)
        if _resmax <= tolerance:
            break
    numiter[0] = _numiter
    resmax[0] = _resmax


def poisson_fill_grids(double[:, :, :] grids,
                       int[:, :, :] masks,
                       double relaxc,
                       double tolerance,
                       unsigned int itermax,
                       int cyclic,
                       int initialize_zonal):
    """
    Fill missing values in grids by iterative relaxation.

    **Arguments:**

    * grids [double[:, :]]
       A 3-dimensional array containing grids of values. The grid
       dimensions are [y, x, n] where n is a non-grid dimension.

    * masks [int[:, :]]
       A 3-dimensional array the same shape as *grids* that contains the
       grid mask. Valid grid points are indicated with a zero value,
       and missing values are indicated by a non-zero value.

    * relaxc [double]
       The relaxation constant, typically 0.45 <= *relaxc* <= 0.6.

    * tolerance [double]
       Numerical tolerance for convergence of the relaxation scheme.
       This value is data dependent.

    * itermax [unsigned int]
       The maximum number of iterations allowed for the relaxation
       scheme.

    * cyclic [int]
       If `0` the input grid is assumed not to be cyclic. If non-zero
       the input grid is assumed to be cyclic in its second dimension.

    * initialize_zonal [int]
       If non-zero, take the zonal mean as the initial guess for missing
       values. If zero, use the value `0` as the initial guess for all
       missing values.

    **Returns:**

    * numiter [numpy.ndarray[unsigned int, ndim=1]
       The number of iterations used to fill each grid.

    * resmax [numpy.ndarray[double, ndim=1]]
       The maximum residual value at the end of the iteration for each
       grid. If this value is larger than the specified *tolerance* then
       the iteration did not converge.

    """
    cdef unsigned int grid_num
    cdef unsigned int nlat = grids.shape[0]
    cdef unsigned int nlon = grids.shape[1]
    cdef unsigned int ngrid = grids.shape[2]
    cdef np.ndarray[UINT32_t, ndim=1] numiter = np.empty([ngrid],
                                                         dtype=np.uint32)
    cdef np.ndarray[FLOAT64_t, ndim=1] resmax = np.empty([ngrid],
                                                         dtype=np.float64)
    if nlat < 3 or nlon < 3:
        raise ValueError('The x and y directions must have at least 3 points, '
                         'got x: {} and y: {}'.format(nlon, nlat))
    if masks.shape[0] != nlat or masks.shape[1] != nlon or \
            masks.shape[2] != ngrid:
        raise ValueError('The dimensions of the grids and the masks must '
                         'match, but found ({}, {}, {}) != ({}, {}, {})'
                         ''.format(nlat, nlon, ngrid,
                                   masks.shape[0],
                                   masks.shape[1],
                                   masks.shape[2]))
    for grid_num in range(ngrid):
        poisson_fill(
            grids[:, :, grid_num],
            masks[:, :, grid_num],
            nlat,
            nlon,
            relaxc,
            tolerance,
            itermax,
            cyclic,
            initialize_zonal,
            &numiter[grid_num],
            &resmax[grid_num])
    return (numiter, resmax)
