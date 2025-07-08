Quick Start Guide
=================

This guide will help you get started with Skyborn in just a few minutes.

First Steps
-----------

Import Skyborn
~~~~~~~~~~~~~~

.. code-block:: python

   import skyborn
   import xarray as xr
   import numpy as np

Basic Data Conversion
~~~~~~~~~~~~~~~~~~~~~

Convert GRIB files to NetCDF format:

.. code-block:: python

   # Simple conversion
   skyborn.grib2nc('input.grib', 'output.nc')

   # With options
   skyborn.convert_grib_to_nc_simple(
       'input.grib',
       'output.nc',
       high_precision=True,
       compress=True
   )

Working with Data
-----------------

Load and Analyze Data
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Load NetCDF data
   ds = xr.open_dataset('output.nc')

   # Basic data info
   print(ds)
   print(ds.data_vars)

Calculate Gradients
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Meridional gradient (∂/∂lat)
   temp_grad_lat = skyborn.calculate_meridional_gradient(
       ds['temperature'],
       'latitude'
   )

   # Zonal gradient (∂/∂lon)
   temp_grad_lon = skyborn.calculate_zonal_gradient(
       ds['temperature'],
       'longitude'
   )

   # General gradient
   pressure_grad = skyborn.calculate_gradient(
       ds['pressure'],
       'level'
   )

Data Manipulation
~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Convert longitude range
   ds_converted = skyborn.convert_longitude_range(ds, target_range=(-180, 180))

   # Linear regression
   slope, intercept, r_value = skyborn.linear_regression(
       ds['temperature'].values.flatten(),
       ds['pressure'].values.flatten()
   )

Advanced Features
-----------------

Interpolation and Regridding
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from skyborn.interp import interpolation, regridding

   # Interpolate data
   interpolated = interpolation.interpolate_data(
       ds['temperature'],
       target_coords={'latitude': np.arange(-90, 91, 1)}
   )

   # Regrid to new grid
   regridded = regridding.regrid_data(
       ds['temperature'],
       target_grid=(180, 360)  # new resolution
   )

Causality Analysis
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Liang causality analysis
   causality_result = skyborn.liang_causality(
       ds['temperature'].values,
       ds['precipitation'].values
   )

   # Granger causality
   granger_result = skyborn.granger_causality(
       ds['temperature'].values,
       ds['pressure'].values,
       max_lag=5
   )

Visualization
~~~~~~~~~~~~~

.. code-block:: python

   from skyborn.plot import plotting, modplot

   # Basic plotting
   fig, ax = plotting.plot_contour(
       ds['temperature'].isel(time=0),
       levels=20,
       title='Temperature Distribution'
   )

   # Specialized atmospheric plots
   fig = modplot.plot_wind_field(
       ds['u_wind'].isel(time=0),
       ds['v_wind'].isel(time=0)
   )

Batch Processing
----------------

Process Multiple Files
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Batch convert GRIB files
   converted_files = skyborn.batch_convert_grib_to_nc(
       input_directory='./grib_data/',
       output_directory='./netcdf_data/',
       pattern='*.grb',
       high_precision=True,
       compress=True
   )

   print(f"Converted {len(converted_files)} files")

Common Workflows
----------------

ERA5 Data Processing
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Optimized for ERA5 data
   skyborn.convert_grib_to_nc(
       grib_files='era5_data.grib',
       output_file='era5_processed.nc',
       ignore_keys=['method', 'type', 'stream'],
       split_keys=['param', 'levtype'],
       data_type='NC_FLOAT',
       unlimited_dimension='time',
       file_kind=4,
       deflate_level=4
   )

Climate Analysis Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Complete analysis workflow
   def analyze_climate_data(grib_file, output_dir):
       # 1. Convert data
       nc_file = f"{output_dir}/converted_data.nc"
       skyborn.grib2nc(grib_file, nc_file, compress=True)

       # 2. Load and process
       ds = xr.open_dataset(nc_file)

       # 3. Calculate gradients
       temp_grad = skyborn.calculate_meridional_gradient(
           ds['temperature'], 'latitude'
       )

       # 4. Causality analysis
       causality = skyborn.liang_causality(
           ds['temperature'].values,
           ds['pressure'].values
       )

       # 5. Save results
       results = xr.Dataset({
           'temperature_gradient': temp_grad,
           'causality_strength': (['time'], causality)
       })
       results.to_netcdf(f"{output_dir}/analysis_results.nc")

       return results

Best Practices
--------------

Performance Tips
~~~~~~~~~~~~~~~~

1. **Use compression** for large datasets
2. **Specify data types** appropriately (NC_SHORT vs NC_FLOAT)
3. **Process in chunks** for very large files
4. **Use appropriate ignore_keys** for your data type

Memory Management
~~~~~~~~~~~~~~~~~

.. code-block:: python

   # For large datasets, use dask
   import dask.array as da

   # Load data lazily
   ds = xr.open_dataset('large_file.nc', chunks={'time': 100})

   # Process with dask
   result = skyborn.calculate_gradient(ds['temperature'], 'latitude')
   result = result.compute()  # Execute computation

Error Handling
~~~~~~~~~~~~~~

.. code-block:: python

   from skyborn.conversion import GribToNetCDFError

   try:
       skyborn.grib2nc('input.grib', 'output.nc')
   except GribToNetCDFError as e:
       print(f"Conversion failed: {e}")
   except FileNotFoundError as e:
       print(f"File not found: {e}")

Next Steps
----------

* Explore the :doc:`api/index` for detailed function documentation
* Check out :doc:`examples/index` for more comprehensive examples
* Read about specific modules in :doc:`modules/index`
