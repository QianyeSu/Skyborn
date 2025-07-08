.. Skyborn documentation master file

Welcome to Skyborn's Documentation!
====================================

Skyborn is a comprehensive Python library for atmospheric data processing, analysis, and visualization. 
It provides tools for working with meteorological and climate datasets, including GRIB to NetCDF conversion, 
data interpolation, gradient calculations, and specialized plotting functions.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :numbered:

   installation
   quickstart
   api/index
   modules/index
   examples/index
   changelog
   contributing

.. toctree::
   :maxdepth: 1
   :caption: Language:

   English <self>
   中文 <../zh_CN/index>

Features
========

* **Data Conversion**: GRIB to NetCDF conversion using eccodes
* **Interpolation**: Advanced interpolation and regridding capabilities  
* **Gradient Calculations**: Meridional, zonal, and vertical gradients
* **Causality Analysis**: Liang and Granger causality analysis
* **Visualization**: Specialized atmospheric data plotting tools
* **ROF Analysis**: River flow and hydrological analysis tools

Quick Example
=============

.. code-block:: python

   import skyborn
   
   # Convert GRIB to NetCDF
   skyborn.grib2nc('input.grib', 'output.nc')
   
   # Calculate gradients
   import xarray as xr
   ds = xr.open_dataset('output.nc')
   grad = skyborn.calculate_gradient(ds['temperature'], 'latitude')

Getting Help
============

* **GitHub Issues**: Report bugs and request features at `GitHub Issues <https://github.com/yourusername/skyborn/issues>`_
* **Documentation**: Full API documentation and examples
* **Email**: Contact the author at suqianye2000@gmail.com

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
