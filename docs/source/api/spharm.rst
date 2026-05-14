Spherical Harmonics (spharm)
=============================

.. automodule:: skyborn.spharm
   :members:
   :show-inheritance:

Overview
--------

The ``skyborn.spharm`` module provides a Python interface to the NCAR SPHEREPACK library for spherical harmonic transforms and analysis. This module is particularly useful for atmospheric and oceanic data analysis, offering efficient spectral transforms and vector field decomposition capabilities.

Origins And Lineage
-------------------

Skyborn's ``spharm`` module is derived from the original
`pyspharm <https://github.com/jswhit/pyspharm>`_ interface written by
Jeff Whitaker. The original package exposed a high-level Python API over
NCAR SPHEREPACK and documented it with an ``epydoc`` HTML bundle.

Skyborn keeps the same core user-facing concepts:

- ``Spharmt`` as the main transform object
- scalar transforms between grid space and spectral space
- vector diagnostics such as vorticity, divergence, streamfunction, and
  velocity potential
- spectral utilities such as regridding, Gaussian latitudes, geodesic points,
  Legendre functions, and point interpolation

Attribution
~~~~~~~~~~~

- Original Python interface: Jeff Whitaker
- Original spectral backend: NCAR SPHEREPACK
- Skyborn integration and maintenance: Qianye Su

Documentation Note
~~~~~~~~~~~~~~~~~~

The historical ``epydoc`` pages summarized:

- module introduction and usage conventions
- method and function inventories
- grid and spectral storage conventions
- Legendre normalization details

Those topics are now maintained here in the Sphinx docs so that the published
documentation, API reference, and packaging stay in one place.

Key Features
------------

- **Spherical Harmonic Transforms**: Forward and inverse transforms between grid and spectral space
- **Vector Field Analysis**: Decomposition of vector fields into rotational and divergent components
- **Spectral Filtering**: Isotropic smoothing and truncation in spectral space
- **Grid Support**: Both regular and Gaussian grids with flexible resolution
- **Optimized Performance**: Modern Fortran kernels with OpenMP parallelization
- **Type Safety**: Comprehensive input validation and error handling
- **Reduced Gaussian Backend**: Experimental packed reduced-Gaussian transforms and regridding without first forcing data onto a rectangular full-Gaussian user workflow

Installation Requirements
--------------------------

The spharm module requires:

- NumPy with F2PY support
- A Fortran compiler (gfortran, ifort, etc.)
- SPHEREPACK library (automatically handled during installation)

Compared With Legacy pyspharm Packaging
---------------------------------------

The original ``pyspharm`` installation flow relied on ``setup.py``-driven
fetch/build steps around SPHEREPACK. In Skyborn this module is built as part of
the package's normal extension build, so users should install ``skyborn``
instead of following the old standalone ``pyspharm`` instructions from the
historical HTML pages.

Legacy Interface Notes
----------------------

The historical ``pyspharm`` HTML pages were not just API stubs; they also
described what the module was for, how the original ``Spharmt`` object was
intended to be used, and which data conventions users were expected to follow.
Those interface notes are kept here so the current documentation preserves the
same context.

Original Interface Scope
~~~~~~~~~~~~~~~~~~~~~~~~

The original ``pyspharm`` docs described the module as a high-level Python
interface to NCAR SPHEREPACK, rather than a one-to-one wrapper around every
Fortran routine. That design is still true in Skyborn:

- the public API is organized around atmospheric and oceanic analysis tasks
- the ``Spharmt`` class is the main working object
- scalar and vector transforms are exposed as direct analysis/synthesis methods
- helper functions handle Gaussian grids, spectral indexing, interpolation, and
  spectral regridding

Original Usage Pattern
~~~~~~~~~~~~~~~~~~~~~~

The legacy HTML documentation used a compact example like this:

.. code-block:: python

    from skyborn.spharm import Spharmt

    spharm = Spharmt(
        144,
        72,
        rsphere=8.0e6,
        gridtype="gaussian",
        legfunc="computed",
    )

This example still captures the intended interface:

- ``nlon`` and ``nlat`` define the transform grid
- ``rsphere`` sets the sphere radius
- ``gridtype`` selects ``"regular"`` or ``"gaussian"``
- ``legfunc`` selects ``"stored"`` or ``"computed"``

The historical note that accompanied this example is still relevant:

- default values are tuned for Earth-scale applications
- ``legfunc="computed"`` reduces persistent memory usage
- ``legfunc="stored"`` is usually faster when transforms are repeated many times

Legacy Method Inventory
~~~~~~~~~~~~~~~~~~~~~~~

The original class page grouped the core ``Spharmt`` methods as the standard
analysis workflow:

- ``grdtospec``: grid-to-spectral transform
- ``spectogrd``: spectral-to-grid transform
- ``getvrtdivspec``: vorticity/divergence spectra from ``u``/``v`` winds
- ``getuv``: wind components from vorticity/divergence spectra
- ``getgrad``: vector gradient from scalar spectral coefficients
- ``getpsichi``: streamfunction and velocity potential from winds
- ``specsmooth``: isotropic spectral smoothing

Legacy Utility Inventory
~~~~~~~~~~~~~~~~~~~~~~~~

The original module page also highlighted the same utility functions now
documented below:

- ``regrid``: spectral regridding with optional smoothing and truncation
- ``gaussian_lats_wts``: Gaussian latitudes and quadrature weights
- ``getspecindx``: spectral coefficient index mapping
- ``legendre``: associated Legendre functions
- ``getgeodesicpts``: icosahedral geodesic points on the sphere
- ``specintrp``: interpolation to an arbitrary point from spectral coefficients

Legacy Conventions Carried Forward
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several conventions in the old HTML docs remain important and are preserved in
Skyborn:

- gridded data is treated as north-to-south in latitude and eastward in longitude
- the Greenwich meridian is the first longitude
- the northernmost latitude is the first row
- odd ``nlat`` includes the equator; even ``nlat`` places it between two rows
- regular grids include poles only when the latitude count is odd
- spectral coefficients follow triangular truncation ordering

These conventions are restated in the formal sections below so users do not
need the historical HTML bundle to understand the interface.

Quick Start
-----------

.. code-block:: python

    import numpy as np
    from skyborn.spharm import Spharmt

    # Create a spherical harmonic transform instance
    nlon, nlat = 144, 72
    spharm = Spharmt(nlon, nlat, gridtype='gaussian', legfunc='stored')

    # Generate some sample data
    data = np.random.rand(nlat, nlon)

    # Transform to spectral space
    spec_coeffs = spharm.grdtospec(data)

    # Transform back to grid space
    data_reconstructed = spharm.spectogrd(spec_coeffs)

    # Analyze vector fields
    u_wind = np.random.rand(nlat, nlon)
    v_wind = np.random.rand(nlat, nlon)

    # Compute streamfunction and velocity potential
    streamfunction, velocity_potential = spharm.getpsichi(u_wind, v_wind)

Classes
-------

Spharmt
~~~~~~~

.. autoclass:: skyborn.spharm.Spharmt
   :members:
   :undoc-members:
   :show-inheritance:

   The main class for spherical harmonic transforms and analysis.

   **Constructor Parameters:**

   - ``nlon`` (int): Number of longitude points
   - ``nlat`` (int): Number of latitude points
   - ``rsphere`` (float): Sphere radius in meters (default: Earth radius)
   - ``gridtype`` (str): Grid type - 'regular' or 'gaussian' (default: 'regular')
   - ``legfunc`` (str): Legendre function handling - 'stored' or 'computed' (default: 'stored')

   **Key Methods:**

   - :meth:`grdtospec`: Grid to spectral transform (spherical harmonic analysis)
   - :meth:`spectogrd`: Spectral to grid transform (spherical harmonic synthesis)
   - :meth:`getvrtdivspec`: Compute vorticity and divergence spectra from wind components
   - :meth:`getuv`: Compute wind components from vorticity and divergence spectra
   - :meth:`getpsichi`: Compute streamfunction and velocity potential from winds
   - :meth:`getgrad`: Compute gradient from scalar spectral coefficients
   - :meth:`specsmooth`: Apply spectral smoothing to grid data

ReducedGaussianSpharmt
~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: skyborn.spharm.ReducedGaussianSpharmt
   :members:
   :undoc-members:
   :show-inheritance:

   Experimental transform interface for packed reduced-Gaussian grids.

   Unlike :class:`skyborn.spharm.Spharmt`, this backend works on packed arrays
   with leading dimension ``ngptot`` instead of rectangular
   ``(nlat, nlon, ...)`` arrays.

   **Constructor Parameters:**

   - ``nloen`` (array-like): number of longitude points for each Gaussian latitude circle, ordered north-to-south
   - ``rsphere`` (float): sphere radius in meters
   - ``backend`` (str): currently ``"ectrans"``

Functions
---------

Grid and Coordinate Utilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.spharm.gaussian_lats_wts

.. autofunction:: skyborn.spharm.getgeodesicpts

.. autofunction:: skyborn.spharm.getspecindx

Mathematical Functions
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.spharm.legendre

.. autofunction:: skyborn.spharm.specintrp

Data Processing
~~~~~~~~~~~~~~~

.. autofunction:: skyborn.spharm.regrid

Unified Workflow Examples
-------------------------

Rectangular-grid Example
~~~~~~~~~~~~~~~~~~~~~~~~

The top-level workflow helpers can now be used directly with the standard
rectangular-grid :class:`skyborn.spharm.Spharmt` interface.

.. code-block:: python

    import numpy as np
    from skyborn.spharm import Spharmt, regrid

    src = Spharmt(72, 37, gridtype="gaussian")
    dst = Spharmt(36, 19, gridtype="gaussian")

    scalar_in = np.random.randn(src.nlat, src.nlon)
    scalar_out = regrid(src, dst, scalar_in, ntrunc=18)

Reduced-grid Workflow
---------------------

Skyborn now also provides an experimental reduced-Gaussian backend for cases
where the source data already lives on a packed reduced Gaussian grid and the
workflow should stay on that native layout as long as possible.

Packed Reduced-grid Convention
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For :class:`skyborn.spharm.ReducedGaussianSpharmt`, gridpoint data is packed as
``(ngptot,)`` or ``(ngptot, *extra_dims)`` where:

- ``nloen[i]`` is the longitude count on latitude circle ``i``
- latitude circles are ordered north-to-south
- each latitude circle is stored as one contiguous longitude block
- the total packed size is ``ngptot = sum(nloen)``

Minimal Example
~~~~~~~~~~~~~~~

.. code-block:: python

    import numpy as np
    from skyborn.spharm import ReducedGaussianSpharmt, regrid

    src = ReducedGaussianSpharmt(
        np.array([20, 24, 28, 24, 20], dtype=np.int32)
    )
    dst = ReducedGaussianSpharmt(
        np.array([12, 16, 20, 24, 20, 16, 12], dtype=np.int32)
    )

    scalar_field = np.zeros(src.ngptot, dtype=np.float64)
    scalar_out = regrid(src, dst, scalar_field, ntrunc=4)

The same top-level ``regrid(...)`` helper is used for both rectangular-grid
and reduced-grid scalar workflows. The backend is selected from the transform
objects passed as ``grdin`` and ``grdout``.

Current Scope
~~~~~~~~~~~~~

The reduced-grid backend currently covers the core packed-grid transform and
diagnostic workflow:

- scalar analysis/synthesis
- vorticity/divergence spectra from ``u``/``v``
- wind synthesis from ``vrt/div`` spectra
- scalar gradient synthesis
- streamfunction / velocity-potential diagnostics
- packed-grid ``regrid(...)``

It is still narrower in scope than the main rectangular-grid
:class:`skyborn.spharm.Spharmt` interface:

- it is reduced-Gaussian-only
- it is still marked experimental
- it is not yet a full replacement for all rectangular-grid ``spharm`` workflows

ECTRANS Status
--------------

The current reduced-grid backend should be understood as:

- a **Skyborn-native** reduced-Gaussian backend
- shaped by ECMWF/OpenIFS-style reduced-grid goals
- **not** yet a direct vendoring of the real OpenIFS ``trans`` runtime

Feature Status Matrix
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 34 18 48

   * - Capability
     - Status
     - Notes
   * - ``grdtospec(...)``
     - Implemented
     - Native reduced-grid scalar analysis path is wired
   * - ``spectogrd(...)``
     - Implemented
     - Native reduced-grid scalar synthesis path is wired
   * - ``getvrtdivspec(...)``
     - Implemented
     - Native reduced-grid vector analysis path is wired
   * - ``getvrtspec(...)`` / ``getdivspec(...)``
     - Implemented
     - Derived from the native ``getvrtdivspec(...)`` path
   * - ``getuv(...)``
     - Implemented
     - Native reduced-grid wind synthesis path is wired
   * - ``getgrad(...)``
     - Implemented
     - Native reduced-grid scalar-gradient synthesis path is wired
   * - ``getpsispec(...)`` / ``getchispec(...)``
     - Implemented
     - Uses local reduced-grid ``vrt/div`` spectra plus inverse-Laplacian path
   * - ``getpsichispec(...)``
     - Implemented
     - Native reduced-grid main path
   * - ``getpsi(...)`` / ``getchi(...)`` / ``getpsichi(...)``
     - Implemented
     - Built from native reduced-grid spectra plus reduced-grid synthesis
   * - ``specsmooth(...)``
     - Implemented
     - Runs through the local reduced-grid spectral chain
   * - ``regrid(...)``
     - Implemented
     - Packed reduced-grid scalar regridding through spectral space
   * - Regular lat-lon backend
     - Not implemented
     - ``ectrans`` currently targets reduced Gaussian only
   * - Full-Gaussian rectangular backend
     - Not implemented
     - Main rectangular-grid workflows still belong to ``Spharmt``
   * - FULLPOS-style change-resolution workflow
     - Not implemented
     - Only the first reduced-grid regrid helpers are present so far
   * - Full public utility surface parity
     - Partial
     - The backend does not yet mirror every ``spharm``-side utility as a reduced-grid-native public workflow

OpenIFS Down-sinking
~~~~~~~~~~~~~~~~~~~~

The current backend does **not** directly vendor the real OpenIFS
``trans`` / ``transi`` runtime. These areas are still outside the current
implementation:

- most of the real OpenIFS ``ifs-source/trans`` runtime sources
- FULLPOS / FIELDS business workflows
- OpenIFS handle/setup/cache runtime model
- MPI / distributed transform infrastructure
- full OpenIFS change-resolution and spectral interpolation workflow

So the present state is better described as a Skyborn-native reduced-grid
backend with OpenIFS-style direction, not as a completed OpenIFS runtime import.

Data Conventions
----------------

Grid Data Format
~~~~~~~~~~~~~~~~~

The gridded data is assumed to be oriented such that:

- ``i=0`` corresponds to the Greenwich meridian (0° longitude)
- ``j=0`` corresponds to the northernmost point (90° latitude for regular grids)
- Grid indices increase eastward and southward
- For odd ``nlat``, the equator is included
- For even ``nlat``, the equator lies halfway between points ``nlat/2`` and ``nlat/2+1``

**Grid Requirements:**

- ``nlat`` must be at least 3
- ``nlon`` must be at least 4
- For optimal performance, ``nlon`` should be a product of small prime numbers

Spectral Data Format
~~~~~~~~~~~~~~~~~~~~

The spectral data uses triangular truncation with coefficients stored in a complex array of size ``(ntrunc+1)*(ntrunc+2)/2``, where ``ntrunc`` is the triangular truncation limit.

**Coefficient Ordering:**

- First coefficient (``nm=0``): ``m=0, n=0``
- Second coefficient (``nm=1``): ``m=0, n=1``
- At ``nm=ntrunc``: ``m=0, n=ntrunc``
- At ``nm=ntrunc+1``: ``m=1, n=1``
- And so on...

The arrays ``indxm`` and ``indxn`` returned by :func:`getspecindx` provide the mapping between linear index ``nm`` and spherical harmonic degrees ``m`` and ``n``.

Legendre Function Normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The associated Legendre polynomials are normalized such that:

.. math::

   \int_0^\pi P_{nm}(\theta)^2 \sin(\theta) \, d\theta = 1

where:

.. math::

   P_{nm}(\theta) = \sqrt{\frac{2n+1}{2} \frac{(n-m)!}{(n+m)!}} \frac{\sin^m(\theta)}{2^n n!} \frac{d^{n+m}}{d(\cos\theta)^{n+m}}(\cos^2\theta - 1)^n

**Coordinate Convention:**

- ``θ = π/2 - φ`` where ``φ`` is latitude and ``θ`` is colatitude
- ``cos(θ) = sin(φ)`` and ``sin(θ) = cos(φ)``
- Special values: ``P_{00}(θ) = √2/2`` and ``P_{10}(θ) = (√6/2)sin(φ)``

Grid Types
----------

Regular Grid
~~~~~~~~~~~~

Equally spaced latitude points from pole to pole. When ``nlat`` is odd, the poles are included.

**Characteristics:**

- Simple and uniform spacing
- Suitable for most applications
- May have convergence issues near poles for high-resolution data

Gaussian Grid
~~~~~~~~~~~~~

Latitudes correspond to the roots of Legendre polynomials, providing optimal quadrature properties.

**Characteristics:**

- Optimal for spherical harmonic transforms
- Better numerical properties than regular grids
- Recommended for high-accuracy applications
- Latitudes computed using :func:`gaussian_lats_wts`

Legendre Function Handling
--------------------------

Stored Mode (``legfunc='stored'``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Legendre functions are precomputed and stored during initialization
- Faster execution for repeated transforms
- Higher memory usage (scales as ``nlat²``)
- Recommended for applications with many repeated transforms

Computed Mode (``legfunc='computed'``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Legendre functions are computed on-the-fly for each transform
- Lower memory usage
- Slower execution for repeated transforms
- Suitable for memory-constrained applications or one-time transforms

Performance Considerations
--------------------------

Memory Usage
~~~~~~~~~~~~

- **Stored mode**: Memory scales as ``O(nlat²)`` for Legendre functions
- **Computed mode**: Memory scales as ``O(nlat)`` for temporary arrays
- **Work arrays**: Additional memory for SPHEREPACK internal arrays

Computational Complexity
~~~~~~~~~~~~~~~~~~~~~~~~~

- **Forward/inverse transforms**: ``O(nlat² × nlon × log(nlon))``
- **Vector transforms**: Approximately 2× the cost of scalar transforms
- **Spectral operations**: ``O(ntrunc²)`` where ``ntrunc ≤ nlat-1``

Optimization Tips
~~~~~~~~~~~~~~~~~

1. **Choose appropriate grid type**: Use Gaussian grids for high-accuracy applications
2. **Optimize nlon**: Use values that are products of small primes (2, 3, 5, 7)
3. **Select legfunc wisely**: Use 'stored' for repeated transforms, 'computed' for one-time use
4. **Minimize truncation**: Higher truncation limits increase computational cost
5. **Batch operations**: Process multiple fields together when possible

Error Handling
--------------

The module provides comprehensive error handling with custom exception types:

- :exc:`ValidationError`: Input parameter validation failures
- :exc:`SpheremackError`: SPHEREPACK library errors
- Standard Python exceptions for other error conditions

Common error scenarios include:

- Invalid grid dimensions (``nlat < 3``, ``nlon < 4``)
- Mismatched array shapes between input arrays
- Spectral truncation limits exceeding ``nlat-1``
- Insufficient memory for large transforms

Examples
--------

Basic Transform Operations
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import numpy as np
    from skyborn.spharm import Spharmt

    # Create transform instance
    nlon, nlat = 128, 64
    spharm = Spharmt(nlon, nlat, gridtype='regular')

    # Create sample data
    lons = np.linspace(0, 360, nlon, endpoint=False)
    lats = np.linspace(90, -90, nlat)
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    # Simple sinusoidal pattern
    data = np.sin(2 * np.pi * lon_grid / 360) * np.cos(np.pi * lat_grid / 180)

    # Transform to spectral space
    spec_coeffs = spharm.grdtospec(data)
    print(f"Spectral coefficients shape: {spec_coeffs.shape}")

    # Transform back to grid space
    data_reconstructed = spharm.spectogrd(spec_coeffs)

    # Check reconstruction accuracy
    error = np.max(np.abs(data - data_reconstructed))
    print(f"Reconstruction error: {error:.2e}")

Vector Field Analysis
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import numpy as np
    from skyborn.spharm import Spharmt

    # Create transform instance
    nlon, nlat = 144, 72
    spharm = Spharmt(nlon, nlat, gridtype='gaussian')

    # Create coordinate grids
    lons = np.linspace(0, 360, nlon, endpoint=False)
    lats, weights = spharm.gaussian_lats_wts(nlat)  # Get Gaussian latitudes
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    # Convert to radians
    lon_rad = np.deg2rad(lon_grid)
    lat_rad = np.deg2rad(lat_grid)

    # Create solid body rotation (example atmospheric flow)
    omega = 7.27e-5  # Earth's rotation rate (rad/s)
    u_wind = -omega * spharm.rsphere * np.sin(lon_rad) * np.cos(lat_rad)
    v_wind = omega * spharm.rsphere * np.cos(lon_rad) * np.cos(lat_rad)

    # Compute streamfunction and velocity potential
    streamfunction, velocity_potential = spharm.getpsichi(u_wind, v_wind)

    # Analyze the flow
    print(f"Streamfunction range: [{np.min(streamfunction):.2e}, {np.max(streamfunction):.2e}]")
    print(f"Velocity potential range: [{np.min(velocity_potential):.2e}, {np.max(velocity_potential):.2e}]")

    # For solid body rotation, velocity potential should be nearly zero
    chi_psi_ratio = np.std(velocity_potential) / np.std(streamfunction)
    print(f"Chi/Psi ratio: {chi_psi_ratio:.2e} (should be << 1 for rotational flow)")

Spectral Filtering
~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import numpy as np
    from skyborn.spharm import Spharmt

    # Create transform instance
    nlon, nlat = 256, 128
    spharm = Spharmt(nlon, nlat, gridtype='gaussian')

    # Create noisy data
    lons = np.linspace(0, 360, nlon, endpoint=False)
    lats, _ = spharm.gaussian_lats_wts(nlat)
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    # Signal: large-scale pattern
    signal = np.sin(2 * np.pi * lon_grid / 360) * np.cos(np.pi * lat_grid / 180)

    # Add high-frequency noise
    noise = 0.1 * np.random.randn(nlat, nlon)
    noisy_data = signal + noise

    # Create smoothing filter (exponential decay)
    wavenumbers = np.arange(nlat)
    smooth_factors = np.exp(-wavenumbers / 20.0)  # Aggressive smoothing

    # Apply spectral smoothing
    filtered_data = spharm.specsmooth(noisy_data, smooth_factors)

    # Compare results
    original_std = np.std(noisy_data)
    filtered_std = np.std(filtered_data)
    print(f"Original std: {original_std:.3f}")
    print(f"Filtered std: {filtered_std:.3f}")
    print(f"Noise reduction: {((original_std - filtered_std) / original_std * 100):.1f}%")

Regridding Between Different Grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import numpy as np
    from skyborn.spharm import Spharmt, regrid

    # Create input grid (coarse regular)
    nlon_in, nlat_in = 72, 36
    spharm_in = Spharmt(nlon_in, nlat_in, gridtype='regular')

    # Create output grid (fine Gaussian)
    nlon_out, nlat_out = 144, 72
    spharm_out = Spharmt(nlon_out, nlat_out, gridtype='gaussian')

    # Create sample data on input grid
    lons_in = np.linspace(0, 360, nlon_in, endpoint=False)
    lats_in = np.linspace(90, -90, nlat_in)
    lon_grid_in, lat_grid_in = np.meshgrid(lons_in, lats_in)

    data_in = np.sin(4 * np.pi * lon_grid_in / 360) * np.cos(2 * np.pi * lat_grid_in / 180)

    # Regrid to output grid
    data_out = regrid(spharm_in, spharm_out, data_in)

    print(f"Input grid shape: {data_in.shape}")
    print(f"Output grid shape: {data_out.shape}")
    print(f"Data range preserved: [{np.min(data_out):.3f}, {np.max(data_out):.3f}]")

Related Modules
---------------

- :mod:`skyborn.interp`: Interpolation functions for atmospheric data
- :mod:`skyborn.gradients`: Gradient calculations using finite differences
- :mod:`skyborn.plot`: Visualization tools for spherical data

References
----------

1. Swarztrauber, P. N. (2003), **"On computing the points and weights for Gauss–Legendre quadrature"**, SIAM Journal on Scientific Computing, 24(3), 945–954.

2. Swarztrauber, P. N., and W. F. Spotz (2000), **"Generalized discrete spherical harmonic transforms"**, Journal of Computational Physics, 159(2), 213–230.

3. Adams, J. C., and P. N. Swarztrauber (1999), **"SPHEREPACK 3.0: A model development facility"**, Monthly Weather Review, 127(8), 1872–1878.

4. Swarztrauber, P. N. (1996), **"Spectral transform methods for solving the shallow-water equations on the sphere"**, Monthly Weather Review, 124(4), 730–744.

5. Williamson, D. L., et al. (1992), **"A standard test set for numerical approximations to the shallow water equations in spherical geometry"**, Journal of Computational Physics, 102(1), 211–224.

.. note::
   This module is based on the NCAR SPHEREPACK library and includes optimizations for modern computing environments. The interface has been enhanced with comprehensive type checking and error handling while maintaining compatibility with the original SPHEREPACK functionality.
