Conversion Module
================

.. currentmodule:: skyborn.conversion

The conversion module provides tools for converting between different atmospheric data formats.

GRIB to NetCDF Conversion
-------------------------

.. automodule:: skyborn.conversion.grib_to_netcdf
   :members:
   :undoc-members:
   :show-inheritance:

Main Functions
~~~~~~~~~~~~~~

.. autofunction:: convert_grib_to_nc

.. autofunction:: convert_grib_to_nc_simple

.. autofunction:: batch_convert_grib_to_nc

Convenience Aliases
~~~~~~~~~~~~~~~~~~~

.. autodata:: grib2nc
   :annotation: = convert_grib_to_nc_simple

.. autodata:: grib_to_netcdf
   :annotation: = convert_grib_to_nc

Exceptions
~~~~~~~~~~

.. autoexception:: GribToNetCDFError

Examples
--------

Basic Usage
~~~~~~~~~~~

Convert a single GRIB file to NetCDF:

.. code-block:: python

   from skyborn.conversion import convert_grib_to_nc_simple
   
   convert_grib_to_nc_simple('input.grib', 'output.nc')

Advanced Usage
~~~~~~~~~~~~~~

Convert with custom options:

.. code-block:: python

   from skyborn.conversion import convert_grib_to_nc
   
   convert_grib_to_nc(
       grib_files=['file1.grib', 'file2.grib'],
       output_file='merged.nc',
       data_type='NC_FLOAT',
       ignore_keys=['type', 'step'],
       unlimited_dimension='time',
       file_kind=4,
       deflate_level=6
   )

Batch Processing
~~~~~~~~~~~~~~~~

Convert multiple files:

.. code-block:: python

   from skyborn.conversion import batch_convert_grib_to_nc
   
   converted = batch_convert_grib_to_nc(
       input_directory='./grib_data/',
       output_directory='./netcdf_data/',
       pattern='*.grb',
       compress=True
   )
