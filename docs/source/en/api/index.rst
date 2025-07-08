API Reference
=============

This section contains the complete API documentation for all Skyborn modules.

.. toctree::
   :maxdepth: 2

   conversion
   calculations
   gradients
   causality
   interpolation
   plotting
   rof

Core Functions
--------------

.. currentmodule:: skyborn

Main Functions
~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   convert_grib_to_nc
   convert_grib_to_nc_simple
   batch_convert_grib_to_nc
   grib2nc
   grib_to_netcdf

Data Processing
~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   convert_longitude_range
   linear_regression

Gradient Calculations
~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   calculate_gradient
   calculate_meridional_gradient
   calculate_zonal_gradient
   calculate_vertical_gradient

Causality Analysis
~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   liang_causality
   granger_causality
