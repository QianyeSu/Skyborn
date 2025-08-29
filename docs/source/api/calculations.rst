Calculation
=====================

The Skyborn calculation module provides statistical, atmospheric, and mathematical functions for climate data analysis.

Atmospheric Physics Functions
-----------------------------

Tropopause Calculations
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.calc.troposphere.tropopause.trop_wmo

.. autofunction:: skyborn.calc.troposphere.tropopause.trop_wmo_profile

.. autofunction:: skyborn.calc.troposphere.xarray.trop_wmo

Geostrophic Wind Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.calc.geostrophic.interface.geostrophic_wind

.. autoclass:: skyborn.calc.geostrophic.interface.GeostrophicWind

.. autofunction:: skyborn.calc.geostrophic.xarray.geostrophic_wind

.. autoclass:: skyborn.calc.geostrophic.xarray.GeostrophicWind

Statistical Functions
---------------------

.. autofunction:: skyborn.calc.linear_regression

.. autofunction:: skyborn.calc.spatial_correlation

.. autofunction:: skyborn.calc.pearson_correlation

.. autofunction:: skyborn.calc.spearman_correlation

.. autofunction:: skyborn.calc.kendall_correlation


Emergent Constraint Functions
-----------------------------

.. autofunction:: skyborn.calc.gaussian_pdf

.. autofunction:: skyborn.calc.emergent_constraint_posterior

.. autofunction:: skyborn.calc.emergent_constraint_prior


Legacy Functions (for backward compatibility)
---------------------------------------------

.. autofunction:: skyborn.calc.calc_GAUSSIAN_PDF

.. autofunction:: skyborn.calc.calc_PDF_EC

.. autofunction:: skyborn.calc.find_std_from_PDF

.. autofunction:: skyborn.calc.calc_PDF_EC_PRIOR


Utility Functions
-----------------

.. autofunction:: skyborn.calc.calculate_potential_temperature

.. autofunction:: skyborn.calc.convert_longitude_range


Example Usage
-------------

**Tropopause Calculation Examples**

.. code-block:: python

   import skyborn as skb
   import numpy as np
   import xarray as xr

   # === 1D Profile Analysis ===
   # Single atmospheric profile
   pressure = np.array([50, 100, 200, 300, 500, 700, 850, 1000])  # hPa
   temperature = np.array([200, 210, 220, 230, 250, 270, 280, 288])  # K

   result = skb.calc.trop_wmo_profile(temperature, pressure)
   print(f"Tropopause: {result['pressure']:.1f} hPa at {result['height']:.0f} m")

   # === xarray Interface - Simplified ===
   # Load atmospheric data with level coordinate
   ds = xr.open_dataset('era5_data.nc')  # Temperature with level coordinate in hPa

   # Auto-generate pressure from level coordinate - no pressure array needed!
   result = skb.calc.troposphere.xarray.trop_wmo(ds.temperature)
   print(f"Global tropopause calculated for {result.pressure.shape}")

   # === 2D Cross-section Analysis ===
   # Meridional cross-section (level, lat)
   temp_meridional = ds.temperature.isel(time=0, lon=0)  # (level, lat)
   result_2d = skb.calc.troposphere.xarray.trop_wmo(temp_meridional)
   # Result shape: (lat,) - tropopause height at each latitude

   # === 4D Climate Analysis ===
   # Multi-year dataset (time, level, lat, lon)
   result_4d = skb.calc.troposphere.xarray.trop_wmo(ds.temperature)
   # Result shape: (time, lat, lon) - preserves time and spatial dimensions

   # Seasonal analysis
   seasonal_mean = result_4d.height.groupby('time.season').mean()

   # === Advanced Usage ===
   # Custom pressure field and WMO criterion
   result = skb.calc.troposphere.xarray.trop_wmo(
       temperature=ds.temperature,
       pressure=ds.pressure,  # Custom pressure field
       lapse_criterion=2.5,   # Custom WMO threshold (K/km)
       pressure_unit='Pa'     # If pressure is in Pascals
   )

**Geostrophic Wind Examples**

.. code-block:: python

   import skyborn as skb
   import numpy as np
   import xarray as xr

   # === NumPy Interface ===
   # Traditional interface with manual parameter specification
   nlat, nlon = 73, 144
   z_data = np.random.randn(nlat, nlon) * 100 + 5500  # Geopotential height [gpm]
   lat = np.linspace(-90, 90, nlat)  # Degrees north
   lon = np.linspace(0, 360, nlon)[:-1]  # Degrees east (0-357.5)

   # Calculate geostrophic wind components
   ug, vg = skb.calc.geostrophic_wind(z_data, lon, lat, 'yx')
   print(f"Wind components: ug{ug.shape}, vg{vg.shape}")

   # Class interface for derived quantities
   gw = skb.calc.GeostrophicWind(z_data, lon, lat, 'yx')
   speed = gw.speed()
   print(f"Max wind speed: {speed.max():.1f} m/s")

   # === xarray Interface - Simplified ===
   # Load geopotential height data
   z = xr.DataArray(
       z_data,
       dims=['lat', 'lon'],
       coords={
           'lat': (['lat'], lat, {'units': 'degrees_north'}),
           'lon': (['lon'], lon, {'units': 'degrees_east'})
       },
       attrs={'long_name': '500 hPa geopotential height', 'units': 'gpm'}
   )

   # Automatic coordinate detection and parameter inference
   from skyborn.calc.geostrophic.xarray import geostrophic_wind
   result = geostrophic_wind(z)  # That's it!

   print(f"Automatic features:")
   print(f"  Longitude cyclic: {result.attrs['longitude_cyclic']}")
   print(f"  Latitude ordering: {result.attrs['latitude_ordering']}")
   print(f"  Output: ug{result.ug.shape}, vg{result.vg.shape}")

   # === Multi-dimensional Examples ===
   # 3D time series (time, lat, lon)
   nt = 12
   z_3d_data = np.random.randn(nt, nlat, nlon) * 100 + 5500
   z_3d = xr.DataArray(
       z_3d_data,
       dims=['time', 'lat', 'lon'],
       coords={
           'time': pd.date_range('2023-01-01', periods=nt, freq='MS'),
           'lat': lat, 'lon': lon
       }
   )

   result_3d = geostrophic_wind(z_3d)
   print(f"3D result: ug{result_3d.ug.shape}, vg{result_3d.vg.shape}")

   # Seasonal analysis
   seasonal_winds = result_3d.ug.groupby('time.season').mean()

   # 4D multi-level data (time, level, lat, lon)
   nz = 17
   levels = [50, 100, 150, 200, 250, 300, 400, 500, 600, 700, 775, 850, 925, 1000]
   z_4d_data = np.random.randn(nt, len(levels), nlat, nlon) * 100 + 5500
   z_4d = xr.DataArray(
       z_4d_data,
       dims=['time', 'level', 'lat', 'lon'],
       coords={
           'time': pd.date_range('2023-01-01', periods=nt, freq='MS'),
           'level': (['level'], levels, {'units': 'hPa'}),
           'lat': lat, 'lon': lon
       }
   )

   result_4d = geostrophic_wind(z_4d)
   print(f"4D result: ug{result_4d.ug.shape}, vg{result_4d.vg.shape}")

   # Class interface with xarray
   from skyborn.calc.geostrophic.xarray import GeostrophicWind
   gw_xr = GeostrophicWind(z)
   speed_xr = gw_xr.speed()
   print(f"xarray class speed: {speed_xr.attrs['standard_name']}")

**Statistical Analysis Examples**

.. code-block:: python

   # Statistical analysis
   x = np.random.randn(100)
   y = 2 * x + np.random.randn(100) * 0.5

   slope, intercept, r_value = skb.linear_regression(x, y)
   correlation = skb.pearson_correlation(x, y)

   # Spatial correlation analysis
   # Create sample spatial data (time, lat, lon)
   n_time, n_lat, n_lon = 120, 36, 72
   spatial_data = np.random.randn(n_time, n_lat, n_lon)
   time_series = np.random.randn(n_time)

   # Calculate spatial correlations efficiently
   corr_map, p_values = skb.spatial_correlation(spatial_data, time_series)

   # Works with xarray too
   data_xr = xr.DataArray(spatial_data, dims=['time', 'lat', 'lon'])
   predictor_xr = xr.DataArray(time_series, dims=['time'])
   corr_xr, p_xr = skb.spatial_correlation(data_xr, predictor_xr)

   # Emergent constraint analysis
   x_values = np.linspace(-3, 3, 100)
   pdf = skb.gaussian_pdf(mu=0, sigma=1, x=x_values)

   # Apply emergent constraint
   posterior_mean, posterior_std = skb.emergent_constraint_posterior(
       prior_mean=3.0, prior_std=1.5,
       obs_mean=0.5, obs_std=0.2,
       relationship_slope=2.0, relationship_intercept=0.1
   )

Emergent Constraints
====================

The emergent constraint module implements methods for reducing uncertainty in climate projections
by leveraging relationships between observable present-day quantities and uncertain future projections.

For a complete interactive example, see :doc:`../notebooks/ecs_emergent_constraints_analysis`.
