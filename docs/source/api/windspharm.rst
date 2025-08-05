windspharm
==========

.. automodule:: skyborn.windspharm
   :members:
   :undoc-members:
   :show-inheritance:

Overview
--------

The ``skyborn.windspharm`` module provides spherical harmonic vector wind analysis capabilities
for atmospheric and oceanic sciences. This module is based on the ajdawson/windspharm project
and has been enhanced with modern Python features including type hints, comprehensive error
handling, and improved documentation.

Main Features
-------------

* **Vector Wind Analysis**: Compute vorticity, divergence, streamfunction, and velocity potential
* **Helmholtz Decomposition**: Separate wind fields into rotational and divergent components
* **Multiple Grid Support**: Works with regular and Gaussian latitude grids
* **Efficient Computation**: Uses optimized spherical harmonic transforms
* **Comprehensive Validation**: Thorough input validation and error handling

Core Classes
------------

VectorWind
~~~~~~~~~~

.. autoclass:: skyborn.windspharm.VectorWind
   :members:
   :inherited-members:
   :show-inheritance:

   .. rubric:: Methods

   **Field Computation**

   .. autosummary::

      ~VectorWind.vorticity
      ~VectorWind.divergence
      ~VectorWind.streamfunction
      ~VectorWind.velocitypotential
      ~VectorWind.sfvp

   **Vector Field Operations**

   .. autosummary::

      ~VectorWind.helmholtz
      ~VectorWind.nondivergentcomponent
      ~VectorWind.irrotationalcomponent
      ~VectorWind.gradient

   **Utility Methods**

   .. autosummary::

      ~VectorWind.truncate
      ~VectorWind.nlat
      ~VectorWind.nlon

Standard Interface
------------------

.. automodule:: skyborn.windspharm.standard
   :members:
   :undoc-members:
   :show-inheritance:

Tools and Utilities
-------------------

.. automodule:: skyborn.windspharm.tools
   :members:
   :undoc-members:
   :show-inheritance:

Common Functions
~~~~~~~~~~~~~~~~

.. automodule:: skyborn.windspharm._common
   :members:
   :undoc-members:
   :show-inheritance:

xarray Interface
----------------

.. automodule:: skyborn.windspharm.xarray
   :members:
   :undoc-members:
   :show-inheritance:

Examples
--------

Basic Wind Analysis
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   from skyborn.windspharm import VectorWind

   # Create sample wind data on a Gaussian grid
   nlat, nlon = 73, 144
   u = np.random.randn(nlat, nlon)  # Zonal wind
   v = np.random.randn(nlat, nlon)  # Meridional wind

   # Initialize VectorWind
   vw = VectorWind(u, v, gridtype='gaussian')

   # Calculate dynamical quantities
   vorticity = vw.vorticity()
   divergence = vw.divergence()
   streamfunction = vw.streamfunction()
   velocity_potential = vw.velocitypotential()

Helmholtz Decomposition
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Decompose wind into rotational and divergent components
   u_rot, v_rot, u_div, v_div = vw.helmholtz()

   # Or get components separately
   u_nondiv, v_nondiv = vw.nondivergentcomponent()
   u_irrot, v_irrot = vw.irrotationalcomponent()

Time Series Analysis
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # For time series data
   nt = 100  # number of time steps
   u_time = np.random.randn(nlat, nlon, nt)
   v_time = np.random.randn(nlat, nlon, nt)

   vw_time = VectorWind(u_time, v_time, gridtype='gaussian')
   vorticity_time = vw_time.vorticity()  # Shape: (nlat, nlon, nt)

Error Handling
--------------

The module provides comprehensive error handling for common issues:

.. code-block:: python

   try:
       # This will raise an error if shapes don't match
       vw = VectorWind(u[:, :50], v, gridtype='gaussian')
   except ValueError as e:
       print(f"Shape mismatch error: {e}")

   try:
       # This will raise an error for invalid grid type
       vw = VectorWind(u, v, gridtype='invalid')
   except ValueError as e:
       print(f"Invalid grid type: {e}")

Performance Considerations
--------------------------

* **Grid Type**: Gaussian grids are generally more efficient for spherical harmonic transforms
* **Legendre Functions**: Use 'stored' mode for repeated computations, 'computed' for memory-constrained environments
* **Array Layout**: Ensure latitude dimension decreases from north to south
* **Data Type**: Use float64 for best numerical accuracy

Mathematical Background
-----------------------

The module implements spherical harmonic analysis for vector fields on the sphere.
Key mathematical concepts include:

**Vorticity** (ζ):
   .. math:: \zeta = \frac{1}{a\cos\phi}\left(\frac{\partial v}{\partial \lambda} - \frac{\partial (u\cos\phi)}{\partial \phi}\right)

**Divergence** (δ):
   .. math:: \delta = \frac{1}{a\cos\phi}\left(\frac{\partial u}{\partial \lambda} + \frac{\partial (v\cos\phi)}{\partial \phi}\right)

**Streamfunction** (ψ) and **Velocity Potential** (χ):
   .. math:: u = -\frac{1}{a}\frac{\partial \psi}{\partial \phi} + \frac{1}{a\cos\phi}\frac{\partial \chi}{\partial \lambda}

   .. math:: v = \frac{1}{a\cos\phi}\frac{\partial \psi}{\partial \lambda} + \frac{1}{a}\frac{\partial \chi}{\partial \phi}

where:
- a is Earth's radius
- φ is latitude
- λ is longitude
- u, v are zonal and meridional wind components

See Also
--------

* :doc:`../modules/spharm` : Spherical harmonic transform utilities
* :doc:`interpolation` : Interpolation and regridding tools
* :doc:`calculations` : Additional atmospheric calculations
