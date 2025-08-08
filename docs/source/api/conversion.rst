Conversion
====================

The Skyborn conversion module provides utilities for converting between different data formats,
particularly GRIB to NetCDF conversion.

GRIB to NetCDF Conversion
-------------------------

.. autofunction:: skyborn.conversion.grib_to_netcdf

.. autofunction:: skyborn.conversion.convert_grib_to_nc

.. autofunction:: skyborn.conversion.convert_grib_to_nc_simple

.. autofunction:: skyborn.conversion.batch_convert_grib_to_nc

.. autofunction:: skyborn.conversion.grib2nc


Example Usage
-------------

.. code-block:: python

   import skyborn as skb

   # Convert single GRIB file to NetCDF
   skb.grib_to_netcdf('input.grib', 'output.nc')

   # Batch convert multiple files
   input_files = ['file1.grib', 'file2.grib', 'file3.grib']
   output_files = ['file1.nc', 'file2.nc', 'file3.nc']
   skb.conversion.batch_convert_grib_to_nc(input_files, output_files)

   # Convert longitude range
   data = skb.conversion.convert_longitude_range(dataset, '180')

Basic GRIB to NetCDF conversion:

.. code-block:: python

   import skyborn as skb

   # Convert GRIB to NetCDF
   success = skb.conversion.grib_to_netcdf('input.grib', 'output.nc')

   # Convert longitude range
   data_converted = skb.conversion.convert_longitude_range(data, '180')

For more examples, see the :doc:`../notebooks/index` section.
