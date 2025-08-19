Spherical Harmonics (spharm)
=============================

.. automodule:: skyborn.spharm
   :members:
   :show-inheritance:

Overview
--------

The ``skyborn.spharm`` module provides a Python interface to the NCAR SPHEREPACK library for spherical harmonic transforms and analysis. This module is particularly useful for atmospheric and oceanic data analysis, offering efficient spectral transforms and vector field decomposition capabilities.

Key Features
------------

- **Spherical Harmonic Transforms**: Forward and inverse transforms between grid and spectral space
- **Vector Field Analysis**: Decomposition of vector fields into rotational and divergent components
- **Spectral Filtering**: Isotropic smoothing and truncation in spectral space
- **Grid Support**: Both regular and Gaussian grids with flexible resolution
- **Optimized Performance**: Modern Fortran kernels with OpenMP parallelization
- **Type Safety**: Comprehensive input validation and error handling

Installation Requirements
--------------------------

The spharm module requires:

- NumPy with F2PY support
- A Fortran compiler (gfortran, ifort, etc.)
- SPHEREPACK library (automatically handled during installation)

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

Functions
---------

Grid and Coordinate Utilities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.spharm.gaussian_lats_wts

.. autofunction:: skyborn.spharm.getgeodesicpts
    :no-index:

.. autofunction:: skyborn.spharm.getspecindx
    :no-index:

Mathematical Functions
~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.spharm.legendre
    :no-index:

.. autofunction:: skyborn.spharm.specintrp
    :no-index:

Data Processing
~~~~~~~~~~~~~~~

.. autofunction:: skyborn.spharm.regrid

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
