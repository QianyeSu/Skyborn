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

**Key Use Case**: Divergent Wind Component Analysis
The module is particularly useful for calculating divergent wind components from 4D atmospheric
data (time, level, lat, lon), which is essential for understanding atmospheric circulation patterns
and energy transport mechanisms.

Main Features
-------------

* **Vector Wind Analysis**: Compute vorticity, divergence, streamfunction, and velocity potential
* **Helmholtz Decomposition**: Separate wind fields into rotational and divergent components
* **Data Preparation Tools**: Robust ``prep_data`` and ``recover_data`` utilities for handling multi-dimensional arrays
* **Parallel Processing Support**: Efficient computation for large time-series datasets
* **Multiple Grid Support**: Works with regular and Gaussian latitude grids
* **Efficient Computation**: Uses optimized spherical harmonic transforms
* **Comprehensive Validation**: Thorough input validation and error handling

Real-World Application Examples
------------------------------

**Computing Divergent Wind Components from 4D Data**

A common use case involves processing atmospheric wind data with shape ``(time, level, lat, lon)``
to extract divergent wind components for circulation analysis. This requires careful data preparation
and can benefit from parallel processing for large datasets.

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

**Key Functions for Multi-dimensional Data Processing**

* **prep_data**: Prepares wind data for spherical harmonic analysis by handling dimension reordering
* **recover_data**: Recovers original data structure after spherical harmonic processing

These tools are essential when working with 4D atmospheric data (time, level, lat, lon) and need
to process individual time slices while preserving the original data structure.

Real-World Usage Pattern
~~~~~~~~~~~~~~~~~~~~~~~~

For processing large atmospheric datasets with parallel computation:

.. code-block:: python

    from skyborn.windspharm.tools import prep_data, recover_data
    from skyborn.windspharm.standard import VectorWind
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from tqdm import tqdm
    import numpy as np

    def calculate_divergent_wind_component(time_idx: int,
                                         zonal_wind: np.ndarray,
                                         meridional_wind: np.ndarray) -> tuple:
        """
        Calculate divergent wind components for a single time step.

        Args:
            time_idx: Time index for processing
            zonal_wind: Zonal wind component array (time, level, lat, lon)
            meridional_wind: Meridional wind component array (time, level, lat, lon)

        Returns:
            tuple: (time_idx, divergent_zonal_wind, divergent_meridional_wind)
        """
        try:
            # Prepare data for spherical harmonic analysis
            # Convert from (level, lat, lon) to analysis format
            zonal_prepared, wind_info = prep_data(zonal_wind[time_idx], 'zyx')
            meridional_prepared, _ = prep_data(meridional_wind[time_idx], 'zyx')

            # Create vector wind object and compute irrotational component
            vector_wind = VectorWind(zonal_prepared, meridional_prepared)
            divergent_u, divergent_v = vector_wind.irrotationalcomponent()

            # Recover original data structure
            divergent_u_recovered = recover_data(divergent_u, wind_info)
            divergent_v_recovered = recover_data(divergent_v, wind_info)

            return time_idx, divergent_u_recovered, divergent_v_recovered

        except Exception as e:
            print(f"Error processing time step {time_idx}: {e}")
            raise

    def compute_divergent_wind_components(zonal_wind: np.ndarray,
                                        meridional_wind: np.ndarray,
                                        max_workers: int = 4) -> tuple:
        """
        Compute divergent wind components using parallel processing.

        Args:
            zonal_wind: 4D array of zonal wind (time, level, lat, lon)
            meridional_wind: 4D array of meridional wind (time, level, lat, lon)
            max_workers: Maximum number of worker threads

        Returns:
            tuple: (divergent_zonal_wind, divergent_meridional_wind)
        """
        input_shape = zonal_wind.shape
        n_timesteps = input_shape[0]

        # Initialize output arrays with same shape as input
        divergent_zonal = np.zeros_like(zonal_wind)
        divergent_meridional = np.zeros_like(meridional_wind)

        # Process all time steps in parallel
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all tasks
            future_to_time = {
                executor.submit(calculate_divergent_wind_component, t,
                              zonal_wind, meridional_wind): t
                for t in range(n_timesteps)
            }

            # Collect results with progress bar
            for future in tqdm(as_completed(future_to_time),
                             total=n_timesteps,
                             desc='Computing divergent wind components'):
                try:
                    time_idx, div_u, div_v = future.result()
                    divergent_zonal[time_idx] = div_u
                    divergent_meridional[time_idx] = div_v

                except Exception as e:
                    print(f"Error in future result: {e}")
                    raise

        return divergent_zonal, divergent_meridional

**Usage Example**:

.. code-block:: python

    # Load your 4D wind data (time, level, lat, lon)
    # zonal_wind.shape = (365, 37, 181, 360)  # Daily data, 37 levels
    # meridional_wind.shape = (365, 37, 181, 360)

    # Import skyborn windspharm components
    from skyborn.windspharm.tools import prep_data, recover_data
    from skyborn.windspharm.standard import VectorWind

    # Compute divergent components efficiently
    div_u, div_v = compute_divergent_wind_components(
        zonal_wind, meridional_wind, max_workers=8
    )

    print(f"Processed {zonal_wind.shape[0]} time steps")
    print(f"Output shape: {div_u.shape}")

This approach is particularly effective for:

* Long time series atmospheric data analysis
* Climate model output processing
* Reanalysis data divergent circulation studies
* Large-scale atmospheric dynamics research

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
