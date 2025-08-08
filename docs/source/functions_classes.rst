Functions and Classes Reference
================================

.. rst-class:: functions-classes-page

This page provides a comprehensive list of all functions and classes available in Skyborn, organized by module for easy navigation.

.. contents:: Quick Navigation
   :local:
   :depth: 2

Core Calculations
-----------------

**Statistical Functions:**

.. currentmodule:: skyborn.calc

.. autosummary::
   :nosignatures:

   linear_regression
   pearson_correlation
   spearman_correlation
   kendall_correlation
   calculate_potential_temperature
   convert_longitude_range

**Emergent Constraint Functions:**

.. autosummary::
   :nosignatures:

   gaussian_pdf
   emergent_constraint_posterior
   emergent_constraint_prior
   calc_GAUSSIAN_PDF
   calc_PDF_EC
   find_std_from_PDF
   calc_PDF_EC_PRIOR

**Mann-Kendall Trend Analysis:**

.. autosummary::
   :nosignatures:

   mann_kendall_test
   mann_kendall_multidim
   mann_kendall_xarray
   trend_analysis
   mk_test
   mk_multidim

Data Conversion
---------------

**Format Conversion Functions:**

.. currentmodule:: skyborn.conversion

.. autosummary::
   :nosignatures:

   convert_grib_to_nc
   convert_grib_to_nc_simple
   batch_convert_grib_to_nc
   grib2nc
   grib_to_netcdf

**Classes:**

.. autosummary::
   :nosignatures:

   GribToNetCDFError

GridFill - Atmospheric Data Interpolation
------------------------------------------

**Core Functions:**

.. currentmodule:: skyborn.gridfill

.. autosummary::
   :nosignatures:

   fill
   fill_cube

**XArray Interface:**

.. currentmodule:: skyborn.gridfill.xarray

.. autosummary::
   :nosignatures:

   fill
   fill_multiple
   validate_grid_coverage

Interpolation and Regridding
----------------------------

**Regridding Classes:**

.. currentmodule:: skyborn.interp

.. autosummary::
   :nosignatures:

   Grid
   Regridder
   NearestRegridder
   BilinearRegridder
   ConservativeRegridder

**Interpolation Functions:**

.. autosummary::
   :nosignatures:

   interp_hybrid_to_pressure
   interp_sigma_to_hybrid
   interp_multidim
   nearest_neighbor_indices
   regrid_dataset

Gradient Calculations
---------------------

**Spatial Gradient Functions:**

.. currentmodule:: skyborn.gradients

.. autosummary::
   :nosignatures:

   calculate_gradient
   calculate_meridional_gradient
   calculate_zonal_gradient
   calculate_vertical_gradient

Causality Analysis
------------------

**Causality Functions:**

.. currentmodule:: skyborn.causality

.. autosummary::
   :nosignatures:

   liang_causality
   granger_causality

Plotting and Visualization
---------------------------

**Core Plotting Functions:**

.. currentmodule:: skyborn.plot

.. autosummary::
   :nosignatures:

   add_equal_axes
   createFigure
   curved_quiver
   add_curved_quiverkey

Spherical Harmonics Analysis
----------------------------

**Core Classes:**

.. currentmodule:: skyborn.spharm

.. autosummary::
   :nosignatures:

   Spharmt

**Spherical Harmonics Functions:**

.. autosummary::
   :nosignatures:

   regrid
   gaussian_lats_wts
   getspecindx
   getgeodesicpts
   legendre
   specintrp

WindSpharm - Vector Wind Analysis
---------------------------------

**Standard Interface:**

.. currentmodule:: skyborn.windspharm.standard

.. autosummary::
   :nosignatures:

   VectorWind

**XArray Interface:**

.. currentmodule:: skyborn.windspharm.xarray

.. autosummary::
   :nosignatures:

   VectorWind

ROF - (Development)
---------------------------------------

.. note::
   The ROF module is currently under development. Documentation will be available in future releases.

Quick Function Lookup
---------------------

**Most Commonly Used Functions:**

**Statistical Analysis:**
- :func:`~skyborn.calc.pearson_correlation` - Pearson correlation coefficient
- :func:`~skyborn.calc.spearman_correlation` - Spearman rank correlation
- :func:`~skyborn.calc.linear_regression` - Linear regression analysis
- :func:`~skyborn.calc.gaussian_pdf` - Gaussian probability density function

**Trend Analysis:**
- :func:`~skyborn.calc.mann_kendall_test` - Mann-Kendall trend test
- :func:`~skyborn.calc.mann_kendall_multidim` - Multidimensional Mann-Kendall test
- :func:`~skyborn.calc.trend_analysis` - Comprehensive trend analysis

**Data Processing:**
- :func:`~skyborn.conversion.grib_to_netcdf` - Convert GRIB to NetCDF format
- :func:`~skyborn.gridfill.fill` - Fill missing data using Poisson solver
- :func:`~skyborn.gridfill.xarray.fill` - XArray interface for grid filling

**Spatial Analysis:**
- :func:`~skyborn.gradients.calculate_gradient` - Calculate spatial gradients
- :func:`~skyborn.gradients.calculate_meridional_gradient` - Meridional gradients
- :func:`~skyborn.gradients.calculate_zonal_gradient` - Zonal gradients

**Causality Analysis:**
- :func:`~skyborn.causality.liang_causality` - Liang causality analysis
- :func:`~skyborn.causality.granger_causality` - Granger causality test

**Vector Wind Analysis:**
- :func:`~skyborn.windspharm.standard.VectorWind` - Standard vector wind analysis
- :func:`~skyborn.windspharm.xarray.VectorWind` - XArray-based vector wind analysis

**Emergent Constraints:**
- :func:`~skyborn.calc.emergent_constraint_posterior` - Calculate posterior distributions
- :func:`~skyborn.calc.emergent_constraint_prior` - Calculate prior distributions

Function Categories
--------------------

**By Data Type:**

*Time Series Analysis:*
- Mann-Kendall functions
- Trend analysis
- Correlation functions

*Spatial Data:*
- Grid filling functions
- Gradient calculations
- Vector wind analysis
- Spherical harmonics

*Data Conversion:*
- GRIB to NetCDF conversion
- Format conversion utilities

*Statistical Methods:*
- Correlation analysis
- Regression methods
- PDF calculations
- Emergent constraints

**By Performance Level:**

*High-Performance (Optimized):*
- :func:`~skyborn.calc.mann_kendall_multidim` - Vectorized trend analysis
- :func:`~skyborn.gridfill.fill` - Cython-accelerated interpolation
- Vector wind analysis functions

*Standard Performance:*
- Statistical correlation functions
- Data conversion utilities
- Plotting functions

.. note::
   This list reflects the actual functions available in Skyborn v0.3.11.
   Some functions may appear as plain text rather than clickable links because:

   - They are internal/private functions not intended for public use
   - They are not properly exported in the module's ``__all__`` list
   - They are part of development branches and not yet fully documented

   All clickable functions have complete documentation and are fully supported.

.. tip::
   Use Ctrl+F (Cmd+F on Mac) to quickly search for specific functions on this page.
   For detailed documentation of each function, click on the function name or visit the specific module documentation.
