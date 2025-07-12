Spherical Harmonics (spharm)
=============================

.. automodule:: skyborn.spharm
   :members:
   :undoc-members:
   :show-inheritance:

Introduction
------------

The ``skyborn.spharm`` module provides a Python interface to the NCAR SPHEREPACK library
for spherical harmonic transforms. This module is adapted from the original
`pyspharm project <https://github.com/jswhit/pyspharm>`_.

It provides a simple interface oriented toward working with atmospheric general circulation model (GCM) data.

Requirements
------------

* numpy, and a Fortran compiler supported by numpy.f2py
* NCAR SPHEREPACK library (included)

Basic Usage
-----------

.. code-block:: python

    import skyborn.spharm as spharm

    # Create instance for 144x72 gaussian grid
    x = spharm.Spharmt(144, 72, rsphere=8e6, gridtype='gaussian', legfunc='computed')

This creates a class instance for spherical harmonic calculations on a 144x72 gaussian grid
on a sphere with radius 8000 km. The associated Legendre functions are recomputed on the fly
(instead of pre-computed and stored).

Default values:
* ``rsphere``: 6.3712e6 (Earth's radius in meters)
* ``gridtype``: 'regular'
* ``legfunc``: 'stored'

Main Class
----------

Spharmt Class
~~~~~~~~~~~~~

.. class:: skyborn.spharm.Spharmt(nlon, nlat, rsphere=6371200.0, gridtype='regular', legfunc='stored')

   Spherical harmonic transform class.

   :param int nlon: Number of longitudes. The grid must be oriented from east to west, with the first point at the Greenwich meridian and the last point at 360-delta degrees east (where delta = 360/nlon degrees). Must be >= 4. Transforms will be faster when nlon is the product of small primes.
   :param int nlat: Number of latitudes. The grid must be oriented from north to south. If nlat is odd the equator is included. If nlat is even the equator will lie half way between points nlat/2 and (nlat/2)+1. Must be >=3.
   :param float rsphere: The radius of the sphere in meters. Default 6371200 (the value for Earth).
   :param str gridtype: 'regular' (default) or 'gaussian'. Regular grids will include the poles and equator if nlat is odd. Gaussian grids never include the poles, but will include the equator if nlat is odd.
   :param str legfunc: 'stored' (default) or 'computed'. If 'stored', associated legendre functions are precomputed and stored when the class instance is created. This uses O(nlat**3) memory, but speeds up the spectral transforms. If 'computed', associated legendre functions are computed on the fly when transforms are requested. This uses O(nlat**2) memory, but slows down the spectral transforms a bit.

   **Instance Variables:**

   .. attribute:: nlon

      Number of longitudes (read-only, set at creation)

   .. attribute:: nlat

      Number of latitudes (read-only, set at creation)

   .. attribute:: rsphere

      The radius of the sphere in meters (read-only, set at creation)

   .. attribute:: gridtype

      'regular' or 'gaussian' (read-only, set at creation)

   .. attribute:: legfunc

      'stored' or 'computed' (read-only, set at creation)

   **Methods:**

   .. method:: grdtospec(datagrid, ntrunc=None)

      Grid to spectral transform (spherical harmonic analysis).

      :param numpy.ndarray datagrid: Rank 2 or 3 numpy float32 array with shape (nlat,nlon) or (nlat,nlon,nt), where nt is the number of grids to be transformed. If datagrid is rank 2, nt is assumed to be 1.
      :param int ntrunc: Optional spectral truncation limit. Default is self.nlat-1.
      :returns: Rank 1 or 2 numpy complex array with shape (ntrunc+1)*(ntrunc+2)/2 or ((ntrunc+1)*(ntrunc+2)/2,nt) containing complex spherical harmonic coefficients resulting from the spherical harmonic analysis of datagrid.
      :rtype: numpy.ndarray

   .. method:: spectogrd(dataspec)

      Spectral to grid transform (spherical harmonic synthesis).

      :param numpy.ndarray dataspec: Rank 1 or 2 numpy complex array with shape (ntrunc+1)*(ntrunc+2)/2 or ((ntrunc+1)*(ntrunc+2)/2,nt) containing complex spherical harmonic coefficients (where ntrunc is the triangular truncation limit and nt is the number of spectral arrays to be transformed). If dataspec is rank 1, nt is assumed to be 1.
      :returns: Rank 2 or 3 numpy float32 array with shape (nlat,nlon) or (nlat,nlon,nt) containing the gridded data resulting from the spherical harmonic synthesis of dataspec.
      :rtype: numpy.ndarray

   .. method:: getuv(vrtspec, divspec)

      Compute vector wind on grid given complex spectral coefficients of vorticity and divergence.

      :param numpy.ndarray vrtspec: Rank 1 or 2 numpy complex array of vorticity spectral coefficients, with shape (ntrunc+1)*(ntrunc+2)/2 or ((ntrunc+1)*(ntrunc+2)/2,nt) (where ntrunc is the triangular truncation and nt is the number of spectral arrays to be transformed). If vrtspec is rank 1, nt is assumed to be 1.
      :param numpy.ndarray divspec: Rank 1 or 2 numpy complex array of divergence spectral coefficients, with shape (ntrunc+1)*(ntrunc+2)/2 or ((ntrunc+1)*(ntrunc+2)/2,nt) (where ntrunc is the triangular truncation and nt is the number of spectral arrays to be transformed). Both vrtspec and divspec must have the same shape.
      :returns: Tuple of rank 2 or 3 numpy float32 arrays containing gridded zonal and meridional winds. Shapes are either (nlat,nlon) or (nlat,nlon,nt).
      :rtype: tuple[numpy.ndarray, numpy.ndarray]

   .. method:: getvrtdivspec(ugrid, vgrid, ntrunc=None)

      Compute spectral coefficients of vorticity and divergence given vector wind.

      :param numpy.ndarray ugrid: Rank 2 or 3 numpy float32 array containing grid of zonal winds. Must have shape (nlat,nlon) or (nlat,nlon,nt), where nt is the number of grids to be transformed. If ugrid is rank 2, nt is assumed to be 1.
      :param numpy.ndarray vgrid: Rank 2 or 3 numpy float32 array containing grid of meridional winds. Must have shape (nlat,nlon) or (nlat,nlon,nt), where nt is the number of grids to be transformed. Both ugrid and vgrid must have the same shape.
      :param int ntrunc: Optional spectral truncation limit. Default is self.nlat-1.
      :returns: Tuple of rank 1 or 2 numpy complex arrays of vorticity and divergence spherical harmonic coefficients with shape (ntrunc+1)*(ntrunc+2)/2 or ((ntrunc+1)*(ntrunc+2)/2,nt).
      :rtype: tuple[numpy.ndarray, numpy.ndarray]

   .. method:: getgrad(chispec)

      Compute vector gradient on grid given complex spectral coefficients.

      :param numpy.ndarray chispec: Rank 1 or 2 numpy complex array with shape (ntrunc+1)*(ntrunc+2)/2 or ((ntrunc+1)*(ntrunc+2)/2,nt) containing complex spherical harmonic coefficients (where ntrunc is the triangular truncation limit and nt is the number of spectral arrays to be transformed). If chispec is rank 1, nt is assumed to be 1.
      :returns: Tuple of rank 2 or 3 numpy float32 arrays containing gridded zonal and meridional components of the vector gradient. Shapes are either (nlat,nlon) or (nlat,nlon,nt).
      :rtype: tuple[numpy.ndarray, numpy.ndarray]

   .. method:: getpsichi(ugrid, vgrid, ntrunc=None)

      Compute streamfunction and velocity potential on grid given vector wind.

      :param numpy.ndarray ugrid: Rank 2 or 3 numpy float32 array containing grid of zonal winds. Must have shape (nlat,nlon) or (nlat,nlon,nt), where nt is the number of grids to be transformed. If ugrid is rank 2, nt is assumed to be 1.
      :param numpy.ndarray vgrid: Rank 2 or 3 numpy float32 array containing grid of meridional winds. Must have shape (nlat,nlon) or (nlat,nlon,nt), where nt is the number of grids to be transformed. Both ugrid and vgrid must have the same shape.
      :param int ntrunc: Optional spectral truncation limit. Default is self.nlat-1.
      :returns: Tuple of rank 2 or 3 numpy float32 arrays of gridded streamfunction and velocity potential. Shapes are either (nlat,nlon) or (nlat,nlon,nt).
      :rtype: tuple[numpy.ndarray, numpy.ndarray]

   .. method:: specsmooth(datagrid, smooth)

      Isotropic spectral smoothing on a sphere.

      :param numpy.ndarray datagrid: Rank 2 or 3 numpy float32 array with shape (nlat,nlon) or (nlat,nlon,nt), where nt is the number of grids to be smoothed. If datagrid is rank 2, nt is assumed to be 1.
      :param numpy.ndarray smooth: Rank 1 array of length nlat containing smoothing factors as a function of total wavenumber.
      :returns: Rank 2 or 3 numpy float32 array with shape (nlat,nlon) or (nlat,nlon,nt) containing the smoothed grids.
      :rtype: numpy.ndarray

Standalone Functions
--------------------

.. function:: skyborn.spharm.gaussian_lats_wts(nlat)

   Compute the gaussian latitudes (in degrees) and quadrature weights.

   :param int nlat: Number of gaussian latitudes desired.
   :returns: Tuple of rank 1 numpy float64 arrays containing gaussian latitudes (in degrees north) and gaussian quadrature weights.
   :rtype: tuple[numpy.ndarray, numpy.ndarray]

.. function:: skyborn.spharm.getspecindx(ntrunc)

   Compute indices of zonal wavenumber and degree for complex spherical harmonic coefficients.

   :param int ntrunc: Spherical harmonic triangular truncation limit.
   :returns: Tuple of rank 1 numpy Int32 arrays containing zonal wavenumber (indxm) and degree (indxn) of spherical harmonic coefficients.
   :rtype: tuple[numpy.ndarray, numpy.ndarray]

.. function:: skyborn.spharm.legendre(lat, ntrunc)

   Calculate associated legendre functions for triangular truncation T(ntrunc), at a given latitude.

   :param float lat: The latitude (in degrees) to compute the associate legendre functions.
   :param int ntrunc: The triangular truncation limit.
   :returns: Rank 1 numpy float32 array containing the (ntrunc+1)*(ntrunc+2)/2 associated legendre functions at latitude lat.
   :rtype: numpy.ndarray

.. function:: skyborn.spharm.getgeodesicpts(m)

   Computes the lat/lon values of the points on the surface of the sphere corresponding to a twenty-sided (icosahedral) geodesic.

   :param int m: The number of points on the edge of a single geodesic triangle. There are 10*(m-1)**2+2 total geodesic points, including the poles.
   :returns: Tuple of rank 1 numpy float32 arrays containing the latitudes and longitudes of the geodesic points (in degrees). These points are nearly evenly distributed on the surface of the sphere.
   :rtype: tuple[numpy.ndarray, numpy.ndarray]

.. function:: skyborn.spharm.regrid(grdin, grdout, datagrid, ntrunc=None, smooth=None)

   Regrid data using spectral interpolation, while performing optional spectral smoothing and/or truncation.

   :param Spharmt grdin: Spharmt class instance describing input grid.
   :param Spharmt grdout: Spharmt class instance describing output grid.
   :param numpy.ndarray datagrid: Data on input grid (grdin.nlat x grdin.nlon). If datagrid is rank 3, last dimension is the number of grids to interpolate.
   :param int ntrunc: Optional spectral truncation limit for datagrid. Default is min(grdin.nlat-1,grdout.nlat-1).
   :param numpy.ndarray smooth: Rank 1 array of length grdout.nlat containing smoothing factors as a function of total wavenumber. Default is no smoothing.
   :returns: Interpolated (and optionally smoothed) array(s) on grdout.nlon x grdout.nlat grid.
   :rtype: numpy.ndarray

.. function:: skyborn.spharm.specintrp(lon, dataspec, legfuncs)

   Spectral interpolation given spherical harmonic coefficients.

   :param float lon: Longitude (in degrees) of point on a sphere to interpolate to.
   :param numpy.ndarray dataspec: Spectral coefficients of function to interpolate.
   :param numpy.ndarray legfuncs: Associated legendre functions with same triangular truncation as dataspec (computed using legendre), computed at latitude of interpolation point.
   :returns: Interpolated value.
   :rtype: float

Grid Conventions
----------------

The gridded data is assumed to be oriented such that:

* ``i=1`` is the Greenwich meridian
* ``j=1`` is the northernmost point
* Grid indices increase eastward and southward
* If ``nlat`` is odd, the equator is included
* If ``nlat`` is even, the equator lies halfway between points ``nlat/2`` and ``(nlat/2)+1``
* ``nlat`` must be at least 3

For regular grids (``gridtype='regular'``), the poles are included when ``nlat`` is odd.

The grid increment in longitude is ``2*pi/nlon`` radians. For example, ``nlon = 72`` for a five degree grid.
``nlon`` must be greater than or equal to 4. The efficiency of the computation is improved when ``nlon``
is a product of small prime numbers.

Spectral Data Conventions
--------------------------

The spectral data is assumed to be in a complex array of dimension ``(ntrunc+1)*(ntrunc+2)/2``.

* ``ntrunc`` is the triangular truncation limit (e.g., ``ntrunc = 42`` for T42)
* ``ntrunc`` must be <= ``nlat-1``
* Coefficients are ordered so that:

  * First (``nm=0``) is ``m=0, n=0``
  * Second is ``m=0, n=1``
  * ``nm=ntrunc`` is ``m=0, n=ntrunc``
  * ``nm=ntrunc+1`` is ``m=1, n=1``
  * etc.

The values of ``m`` (degree) and ``n`` (order) as a function of the index ``nm`` are given by the arrays
``indxm``, ``indxn`` returned by ``getspecindx``.

Performance Notes
-----------------

The associated Legendre polynomials are normalized so that the integral
``(pbar(n,m,theta)**2)*sin(theta)`` on the interval ``theta=0`` to ``pi`` is 1.

Quantities needed to compute spherical harmonics are:

* **Pre-computed and stored** when ``legfunc='stored'`` (default) - faster for repeated use but uses more memory
* **Recomputed on the fly** when ``legfunc='computed'`` - slower but uses less memory

Storage requirements:
* ``legfunc='stored'``: increases like ``nlat**2``
* ``legfunc='computed'``: increases like ``nlat**3``

For repeated method invocations on a single class instance, ``legfunc='stored'`` will always be faster.

Usage Examples
--------------

Basic Transform Operations
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import numpy as np
    import skyborn.spharm as spharm

    # Create a spherical harmonic transform instance
    sht = spharm.Spharmt(144, 72, gridtype='gaussian')

    # Generate some sample data on the grid
    lats, lons = np.meshgrid(sht.lats, sht.lons, indexing='ij')
    data = np.sin(2 * np.pi * lons / 360) * np.cos(np.pi * lats / 180)

    # Transform to spectral space
    spec_coeffs = sht.grdtospec(data)

    # Transform back to grid space
    reconstructed_data = sht.spectogrd(spec_coeffs)

Wind Field Analysis
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Compute vorticity and divergence from wind fields
    ugrid = ...  # zonal wind component
    vgrid = ...  # meridional wind component

    vrt_spec, div_spec = sht.getvrtdivspec(ugrid, vgrid)

    # Compute streamfunction and velocity potential
    psi_grid, chi_grid = sht.getpsichi(ugrid, vgrid)

    # Reconstruct winds from vorticity and divergence
    u_reconstructed, v_reconstructed = sht.getuv(vrt_spec, div_spec)

Spectral Regridding
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Create input and output grids
    grid_in = spharm.Spharmt(72, 36)    # T36 grid
    grid_out = spharm.Spharmt(144, 72)  # T72 grid

    # Regrid data with spectral smoothing
    smooth_factors = np.ones(72)
    smooth_factors[30:] = 0.5  # Smooth high wavenumbers

    regridded_data = spharm.regrid(grid_in, grid_out, input_data, smooth=smooth_factors)

Attribution
-----------

This code is adapted from the original `pyspharm project <https://github.com/jswhit/pyspharm>`_
by Jeff Whitaker. The original SPHEREPACK library was developed by NCAR.

API Reference
-------------
