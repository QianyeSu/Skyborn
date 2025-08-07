GridFill
========

The :mod:`skyborn.gridfill` module provides functions to fill missing values in gridded data using iterative relaxation methods to solve Poisson's equation.

This module is particularly useful for meteorological and oceanographic data where spatial interpolation of missing values is needed while preserving physical constraints.

Overview
--------

The gridfill algorithm solves the 2D Poisson equation:

.. math::
   \nabla^2 \phi = 0

where :math:`\phi` represents the field to be filled. The iterative relaxation scheme converges to a solution that smoothly interpolates missing values while preserving the boundary conditions from observed data.

Key Features
~~~~~~~~~~~~

* Poisson equation solver for gap filling
* Support for cyclic and non-cyclic boundaries
* Configurable initialization (zeros or zonal mean)
* Multi-dimensional array support
* Integration with iris cubes and xarray DataArrays

Interfaces
~~~~~~~~~~

The gridfill module provides three main interfaces:

1. **Standard Interface**: Works with numpy arrays and masked arrays
2. **Iris Interface**: Preserves Iris cube metadata and coordinates
3. **xarray Interface**: Modern interface with automatic coordinate detection

Standard Interface
------------------

.. automodule:: skyborn.gridfill.gridfill
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Main Functions
~~~~~~~~~~~~~~

.. autofunction:: skyborn.gridfill.fill
.. autofunction:: skyborn.gridfill.fill_cube

xarray Interface
----------------

The xarray interface provides modern, metadata-aware gap filling with automatic coordinate detection.

.. automodule:: skyborn.gridfill.xarray
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Main Functions
~~~~~~~~~~~~~~

.. autofunction:: skyborn.gridfill.xarray.fill
.. autofunction:: skyborn.gridfill.xarray.fill_multiple
.. autofunction:: skyborn.gridfill.xarray.validate_grid_coverage

Examples
--------

Basic Usage with Numpy Arrays
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   import numpy.ma as ma
   from skyborn.gridfill import fill

   # Create test data with missing values
   data = np.random.rand(50, 100)
   mask = np.zeros_like(data, dtype=bool)
   mask[20:30, 40:60] = True  # Create a gap
   masked_data = ma.array(data, mask=mask)

   # Fill the missing values
   filled_data, converged = fill(masked_data, xdim=1, ydim=0, eps=1e-4)
   print(f"Convergence: {converged[0]}")

Modern xarray Interface
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import xarray as xr
   import numpy as np
   from skyborn.gridfill.xarray import fill

   # Load data with missing values
   data = xr.open_dataarray('sst_with_gaps.nc')

   # Fill missing values preserving metadata
   filled = fill(data, eps=1e-4)
   print(f"Attributes preserved: {bool(filled.attrs)}")

   # Advanced usage with validation
   from skyborn.gridfill.xarray import validate_grid_coverage

   validation = validate_grid_coverage(data, min_coverage=0.2)
   if validation['valid']:
       filled = fill(data, eps=1e-5, verbose=True)
   else:
       print("Data quality issues:", validation['messages'])

Working with Time Series
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   import pandas as pd
   import xarray as xr
   from skyborn.gridfill.xarray import fill

   # Create time series data with gaps
   lons = np.linspace(0, 360, 72, endpoint=False)
   lats = np.linspace(-90, 90, 36)
   time = pd.date_range('2020-01-01', periods=12, freq='M')

   # Create DataArray with metadata
   data = xr.DataArray(
       np.random.rand(12, 36, 72),
       coords={'time': time, 'lat': lats, 'lon': lons},
       dims=['time', 'lat', 'lon'],
       attrs={'units': 'K', 'long_name': 'temperature'}
   )

   # Add some missing values
   data = data.where(np.random.rand(*data.shape) > 0.1)

   # Fill with custom settings
   filled = fill(
       data,
       eps=1e-5,
       relax=0.55,
       initzonal=True,
       verbose=True
   )

Multiple Dataset Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from skyborn.gridfill.xarray import fill_multiple

   # Fill multiple related variables consistently
   temp_filled, humid_filled = fill_multiple(
       [temperature_data, humidity_data],
       eps=1e-4,
       verbose=True
   )

Cyclic Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # For global data, cyclic boundaries are detected automatically
   global_sst = xr.open_dataarray('global_sst.nc')
   filled_sst = fill(global_sst, eps=1e-4)  # Automatic cyclic detection

   # Or specify explicitly
   filled_sst = fill(global_sst, eps=1e-4, cyclic=True)

Performance Considerations
--------------------------

Algorithm Parameters
~~~~~~~~~~~~~~~~~~~~

* **eps**: Convergence tolerance. Smaller values give higher accuracy but require more iterations
* **relax**: Relaxation parameter (0 < relax < 1). Values 0.45-0.6 typically work well
* **itermax**: Maximum iterations. Increase for difficult convergence cases
* **initzonal**: Use zonal mean initialization for data with strong zonal patterns

Best Practices
~~~~~~~~~~~~~~

1. **Validate data quality** before filling using :func:`validate_grid_coverage`
2. **Use appropriate convergence tolerance** - 1e-4 is often sufficient
3. **Enable cyclic boundaries** for global longitude data
4. **Consider zonal initialization** for atmospheric data with strong zonal patterns
5. **Monitor convergence** with verbose output for difficult cases

Mathematical Details
--------------------

Algorithm Description
~~~~~~~~~~~~~~~~~~~~~

The gridfill algorithm uses a finite difference relaxation scheme to solve:

.. math::
   \nabla^2 \phi = \frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2} = 0

The iterative update formula is:

.. math::
   \phi_{i,j}^{(k+1)} = \phi_{i,j}^k + \text{relax} \times \frac{\text{residual}_{i,j}}{4}
png
where the residual is computed using a 5-point stencil:

.. math::
   \text{residual}_{i,j} = \phi_{i+1,j} + \phi_{i-1,j} + \phi_{i,j+1} + \phi_{i,j-1} - 4\phi_{i,j}

Boundary Conditions
~~~~~~~~~~~~~~~~~~~

* **Non-cyclic**: Natural boundary conditions at edges
* **Cyclic**: Periodic boundary conditions for longitude (wrapping)
* **Missing data**: Dirichlet boundary conditions from observed values

Convergence Criteria
~~~~~~~~~~~~~~~~~~~~

The algorithm converges when:

.. math::
   \max(|\text{residual}_{i,j}|) < \epsilon

for all missing data points.

References
----------

This implementation is based on the gridfill package by Andrew Dawson:
https://github.com/ajdawson/gridfill

The algorithm is commonly used in meteorological and oceanographic applications for:

* Sea surface temperature gap filling
* Satellite data reconstruction
* Climate model data processing
* Observational data quality control

See Also
--------

* :mod:`skyborn.interp` - For other interpolation methods
* :mod:`skyborn.windspharm` - For spherical harmonic analysis
* :mod:`skyborn.calc` - For atmospheric calculations
