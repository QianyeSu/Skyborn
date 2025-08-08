Functions and Classes Reference
================================

.. rst-class:: functions-classes-page

This page provides a comprehensive list of all functions and classes available in Skyborn, organized by module for easy navigation.

.. contents:: Quick Navigation
   :local:
   :depth: 2

Core Calculations
-----------------

Statistical Functions
~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: skyborn.calc

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function
     - Description
   * - :func:`skyborn.calc.linear_regression`
     - Linear regression analysis
   * - :func:`skyborn.calc.pearson_correlation`
     - Pearson correlation coefficient
   * - :func:`skyborn.calc.spearman_correlation`
     - Spearman rank correlation
   * - :func:`skyborn.calc.kendall_correlation`
     - Kendall's tau correlation
   * - :func:`skyborn.calc.calculate_potential_temperature`
     - Potential temperature calculation
   * - :func:`skyborn.calc.convert_longitude_range`
     - Longitude coordinate conversion

Emergent Constraints
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function
     - Description
   * - :func:`skyborn.calc.gaussian_pdf`
     - Gaussian probability density function
   * - :func:`skyborn.calc.emergent_constraint_posterior`
     - Posterior probability calculation
   * - :func:`skyborn.calc.emergent_constraint_prior`
     - Prior probability calculation
   * - :func:`skyborn.calc.calc_GAUSSIAN_PDF`
     - Gaussian PDF calculation
   * - :func:`skyborn.calc.calc_PDF_EC`
     - Emergent constraint PDF
   * - :func:`skyborn.calc.find_std_from_PDF`
     - Standard deviation from PDF
   * - :func:`skyborn.calc.calc_PDF_EC_PRIOR`
     - Prior PDF for emergent constraints

Trend Analysis
~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function
     - Description
   * - :func:`skyborn.calc.mann_kendall_test`
     - Mann-Kendall trend test
   * - :func:`skyborn.calc.mann_kendall_multidim`
     - Multidimensional trend analysis
   * - :func:`skyborn.calc.mann_kendall_xarray`
     - XArray Mann-Kendall implementation
   * - :func:`skyborn.calc.trend_analysis`
     - Comprehensive trend analysis

Data Conversion
---------------

Conversion Functions
~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: skyborn.conversion

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function
     - Description
   * - :func:`skyborn.conversion.convert_grib_to_nc`
     - Convert GRIB to NetCDF
   * - :func:`skyborn.conversion.convert_grib_to_nc_simple`
     - Simple GRIB to NetCDF conversion
   * - :func:`skyborn.conversion.batch_convert_grib_to_nc`
     - Batch conversion utility
   * - :func:`skyborn.conversion.grib2nc`
     - GRIB to NetCDF converter
   * - :func:`skyborn.conversion.grib_to_netcdf`
     - GRIB to NetCDF transformation

Conversion Classes
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Class
     - Description
   * - :class:`skyborn.conversion.GribToNetCDFError`
     - Conversion error handling

GridFill - Atmospheric Data Interpolation
------------------------------------------

Core GridFill Functions
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: skyborn.gridfill

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function
     - Description
   * - :func:`skyborn.gridfill.fill`
     - Fill missing data using Poisson grid filling
   * - :func:`skyborn.gridfill.fill_multiple`
     - Fill multiple arrays simultaneously

XArray GridFill Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function
     - Description
   * - :func:`skyborn.gridfill.xarray.fill`
     - XArray interface for data filling
   * - :func:`skyborn.gridfill.xarray.fill_multiple`
     - Fill multiple variables
   * - :func:`skyborn.gridfill.xarray.validate_grid_coverage`
     - Validate grid coverage

Interpolation and Regridding
-----------------------------

Regridding Classes
~~~~~~~~~~~~~~~~~~

.. currentmodule:: skyborn.interp

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Class
     - Description
   * - :class:`skyborn.interp.Grid`
     - Grid definition class
   * - :class:`skyborn.interp.Regridder`
     - Base regridding class
   * - :class:`skyborn.interp.NearestRegridder`
     - Nearest neighbor regridding
   * - :class:`skyborn.interp.BilinearRegridder`
     - Bilinear interpolation regridding
   * - :class:`skyborn.interp.ConservativeRegridder`
     - Conservative regridding

Interpolation Functions
~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function
     - Description
   * - :func:`skyborn.interp.interp_hybrid_to_pressure`
     - Hybrid to pressure interpolation
   * - :func:`skyborn.interp.interp_sigma_to_hybrid`
     - Sigma to hybrid interpolation
   * - :func:`skyborn.interp.interp_multidim`
     - Multidimensional interpolation
   * - :func:`skyborn.interp.nearest_neighbor_indices`
     - Find nearest neighbor indices
   * - :func:`skyborn.interp.regrid_dataset`
     - Regrid entire dataset

Spatial Gradients
-----------------

Gradient Functions
~~~~~~~~~~~~~~~~~~

.. currentmodule:: skyborn.gradients

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function
     - Description
   * - :func:`skyborn.gradients.calculate_gradient`
     - Calculate spatial gradients
   * - :func:`skyborn.gradients.calculate_meridional_gradient`
     - Meridional gradient calculation
   * - :func:`skyborn.gradients.calculate_zonal_gradient`
     - Zonal gradient calculation
   * - :func:`skyborn.gradients.calculate_vertical_gradient`
     - Vertical gradient calculation

Plotting and Visualization
---------------------------

Core Plotting Functions
~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: skyborn.plot

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function
     - Description
   * - :func:`skyborn.plot.plot_field`
     - Plot 2D scalar fields
   * - :func:`skyborn.plot.plot_vector_field`
     - Plot vector fields
   * - :func:`skyborn.plot.plot_streamlines`
     - Plot streamlines
   * - :func:`skyborn.plot.plot_contour`
     - Contour plotting

Specialized Plotting
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function
     - Description
   * - :func:`skyborn.plot.curved_quiver_plot`
     - Curved quiver plots
   * - :func:`skyborn.plot.modplot`
     - Module-specific plotting utilities

Causality Analysis
------------------

Causality Functions
~~~~~~~~~~~~~~~~~~~

.. currentmodule:: skyborn.causality

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function
     - Description
   * - :func:`skyborn.causality.granger_causality`
     - Granger causality analysis
   * - :func:`skyborn.causality.liang_causality`
     - Liang information flow analysis

Spherical Harmonics
-------------------

Spherical Harmonics Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: skyborn.spharm

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function
     - Description
   * - :func:`skyborn.spharm.Spharmt`
     - Main class for spherical harmonic transforms with grid-spectral conversions
   * - :func:`skyborn.spharm.gaussian_lats_wts`
     - Compute gaussian latitudes (in degrees) and quadrature weights for spherical grids
   * - :func:`skyborn.spharm.regrid`
     - Spectral re-gridding with optional smoothing and/or truncation for data interpolation
   * - :func:`skyborn.spharm.getspecindx`
     - Compute indices of zonal wavenumber and degree for complex spherical harmonic coefficients
   * - :func:`skyborn.spharm.getgeodesicpts`
     - Compute points on sphere surface corresponding to icosahedral geodesic
   * - :func:`skyborn.spharm.legendre`
     - Compute associated Legendre functions for spherical harmonic calculations
   * - :func:`skyborn.spharm.specintrp`
     - Spectral interpolation to arbitrary point on sphere given harmonic coefficients

Vector Wind Analysis
--------------------

Standard Interface
~~~~~~~~~~~~~~~~~~

.. currentmodule:: skyborn.windspharm.standard

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Class
     - Description
   * - :class:`skyborn.windspharm.standard.VectorWind`
     - Standard vector wind analysis

XArray Interface
~~~~~~~~~~~~~~~~

.. currentmodule:: skyborn.windspharm.xarray

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Class
     - Description
   * - :class:`skyborn.windspharm.xarray.VectorWind`
     - XArray vector wind analysis

Documentation Structure Guide
-----------------------------

Complete Page Reference
~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60
   :class: skyborn-function-table

   * - Page/Section
     - Content & Purpose
   * - **Main Documentation Pages**
     -
   * - :doc:`index`
     - Main documentation homepage
   * - :doc:`installation`
     - Installation instructions
   * - :doc:`quickstart`
     - Quick start guide
   * - :doc:`functions_classes`
     - This page - complete function reference
   * - :doc:`contributing`
     - Contribution guidelines
   * - :doc:`changelog`
     - Version history and changes
   * - **API Reference**
     -
   * - :doc:`api/index`
     - API documentation overview
   * - :doc:`api/calculations`
     - Statistical and calculation functions
   * - :doc:`api/mann_kendall`
     - Mann-Kendall trend analysis
   * - :doc:`api/causality`
     - Causality analysis methods
   * - :doc:`api/conversion`
     - Data format conversion
   * - :doc:`api/gridfill`
     - Grid interpolation functions
   * - :doc:`api/interpolation`
     - Regridding and interpolation
   * - :doc:`api/gradients`
     - Spatial gradient calculations
   * - :doc:`api/plotting`
     - Visualization functions
   * - :doc:`api/spharm`
     - Spherical harmonics analysis
   * - :doc:`api/windspharm`
     - Vector wind analysis
   * - **Tutorials & Examples**
     -
   * - :doc:`notebooks/index`
     - Jupyter notebook tutorials overview
   * - :doc:`notebooks/gridfill_tutorial`
     - GridFill interpolation tutorial
   * - :doc:`notebooks/mann_kendall_tutorial`
     - Trend analysis tutorial
   * - :doc:`notebooks/windspharm_tutorial`
     - Vector wind analysis tutorial
   * - :doc:`notebooks/ecs_emergent_constraints_analysis`
     - Emergent constraints tutorial
   * - **Module Documentation**
     -
   * - :doc:`modules/index`
     - Auto-generated module index
   * - :doc:`modules/skyborn`
     - Complete module structureFunction Usage Categories
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 70
   :class: skyborn-function-table

   * - Function/Class
     - Usage Category & Description
   * - **Statistical Analysis**
     -
   * - :func:`skyborn.calc.pearson_correlation`
     - Correlation analysis
   * - :func:`skyborn.calc.spearman_correlation`
     - Non-parametric correlation
   * - :func:`skyborn.calc.linear_regression`
     - Regression analysis
   * - :func:`skyborn.calc.gaussian_pdf`
     - Probability distributions
   * - **Trend Analysis**
     -
   * - :func:`skyborn.calc.mann_kendall_test`
     - Single series trend detection
   * - :func:`skyborn.calc.mann_kendall_multidim`
     - Multidimensional trend analysis
   * - :func:`skyborn.calc.trend_analysis`
     - Comprehensive trend analysis
   * - **Data Processing**
     -
   * - :func:`skyborn.conversion.grib_to_netcdf`
     - Format conversion
   * - :func:`skyborn.gridfill.fill`
     - Missing data interpolation
   * - :func:`skyborn.gridfill.xarray.fill`
     - XArray data filling
   * - **Spatial Analysis**
     -
   * - :func:`skyborn.gradients.calculate_gradient`
     - Spatial derivatives
   * - :func:`skyborn.gradients.calculate_meridional_gradient`
     - Meridional gradients
   * - :func:`skyborn.gradients.calculate_zonal_gradient`
     - Zonal gradients
   * - **Causality Analysis**
     -
   * - :func:`skyborn.causality.liang_causality`
     - Information flow analysis
   * - :func:`skyborn.causality.granger_causality`
     - Granger causality testing
   * - **Vector Analysis**
     -
   * - :class:`skyborn.windspharm.standard.VectorWind`
     - Standard vector wind analysis
   * - :class:`skyborn.windspharm.xarray.VectorWind`
     - XArray-based vector analysis
   * - **Emergent Constraints**
     -
   * - :func:`skyborn.calc.emergent_constraint_posterior`
     - Posterior distributions
   * - :func:`skyborn.calc.emergent_constraint_prior`
     - Prior distributions

.. note::

   This page provides a comprehensive list of all functions and classes available in Skyborn. Each entry links to detailed API documentation with parameter descriptions, examples, and usage guidelines.

.. tip::

   **Quick Navigation**: Use Ctrl+F (Cmd+F on Mac) to quickly search for specific functions or keywords within this page. You can also use the Quick Navigation menu above to jump to different sections.
