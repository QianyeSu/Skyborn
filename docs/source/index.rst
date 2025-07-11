.. Skyborn documentation master file

.. image:: _static/SkyBornLogo.svg
   :align: center
   :width: 600px
   :alt: Skyborn Logo

Skyborn is a comprehensive Python package for climate data analysis, featuring advanced statistical methods, emergent constraint techniques, and data conversion utilities.

.. toctree::
   :maxdepth: 2

   quickstart
   installation
   notebooks/index
   gallery
   api/index
   modules/index
   changelog
   contributing

Features
--------

* **Emergent Constraint Methods**: Advanced statistical techniques for reducing uncertainty in climate projections
* **GRIB to NetCDF Conversion**: Robust data format conversion utilities
* **Statistical Analysis**: Comprehensive statistical and correlation functions
* **Gradient Calculations**: Spatial and temporal gradient analysis tools
* **Causality Analysis**: Granger and Liang causality methods
* **Visualization**: High-quality plotting and visualization capabilities

Quick Start
-----------

Install Skyborn and start analyzing climate data:

.. code-block:: bash

   pip install skyborn

.. code-block:: python

   import skyborn as skb

   # Use emergent constraint methods
   pdf = skb.gaussian_pdf(mu=0, sigma=1, x=x_values)

   # Convert GRIB to NetCDF
   skb.grib_to_netcdf('input.grib', 'output.nc')

   # Statistical analysis
   correlation = skb.pearson_correlation(x_data, y_data)

Reference Implementation
========================

Our emergent constraint methods are adapted from the work by blackcata:
https://github.com/blackcata/Emergent_Constraints/tree/master

Based on Cox, P. M., et al. (2013). Nature, 494(7437), 341-344.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
