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

Genesis Potential Index (GPI) / Tropical Cyclone Potential Intensity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.calc.GPI.interface.potential_intensity

.. autoclass:: skyborn.calc.GPI.interface.PotentialIntensityCalculator

.. autofunction:: skyborn.calc.GPI.xarray.potential_intensity

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

**Tropical Cyclone Potential Intensity (GPI) Examples**

.. code-block:: python

   import skyborn as skb
   import numpy as np
   import xarray as xr

   # === Single Profile Calculation ===
   # Atmospheric profile data
   pressure_levels = np.array([1000, 850, 700, 500, 300, 200, 100])  # mb
   temperature = np.array([300, 290, 280, 270, 250, 230, 210])  # K
   mixing_ratio = np.array([0.015, 0.010, 0.005, 0.002, 0.0005, 0.0001, 0.00005])  # kg/kg

   # Surface conditions
   sst = 302.0  # K (29°C - warm tropical ocean)
   psl = 101325.0  # Pa (sea level pressure)

   # NumPy interface
   from skyborn.calc.GPI.interface import potential_intensity
   min_p, pi, err = potential_intensity(
       sst, psl, pressure_levels, temperature, mixing_ratio
   )
   print(f"Potential Intensity: {pi:.1f} m/s, Min Pressure: {min_p:.1f} mb")
   print(f"Error flag: {err} (1 = success)")

   # === xarray Interface with Automatic Unit Conversion ===
   from skyborn.calc.GPI.xarray import potential_intensity

   # Create xarray data with various units (will be auto-converted)
   temp_xr = xr.DataArray(
       temperature - 273.15,  # Convert to Celsius for demo
       dims=['level'],
       attrs={'units': '°C'}  # Will be auto-converted to K
   )

   mixr_xr = xr.DataArray(
       mixing_ratio * 1000,  # Convert to g/kg for demo
       dims=['level'],
       attrs={'units': 'g/kg'}  # Will be auto-converted to kg/kg
   )

   pres_xr = xr.DataArray(
       pressure_levels,
       dims=['level'],
       attrs={'units': 'mb'}
   )

   # Automatic unit conversion and dimension detection
   result = potential_intensity(sst, psl, pres_xr, temp_xr, mixr_xr)
   print(f"xarray PI: {result.pi.values:.1f} m/s")
   print(f"Minimum pressure: {result.min_pressure.values:.1f} mb")

   # === 3D Gridded Data (Spatial Analysis) ===
   nlat, nlon = 20, 30
   nlevs = len(pressure_levels)

   # Create 3D atmospheric data
   temp_3d = xr.DataArray(
       np.random.randn(nlevs, nlat, nlon) * 5 + 280,
       dims=['level', 'lat', 'lon'],
       coords={
           'level': pressure_levels,
           'lat': np.linspace(-30, 30, nlat),
           'lon': np.linspace(0, 360, nlon, endpoint=False)
       },
       attrs={'units': 'K', 'long_name': 'Temperature'}
   )

   mixr_3d = xr.DataArray(
       np.random.rand(nlevs, nlat, nlon) * 0.01 + 0.005,
       dims=['level', 'lat', 'lon'],
       coords=temp_3d.coords,
       attrs={'units': 'kg/kg', 'long_name': 'Water vapor mixing ratio'}
   )

   sst_3d = xr.DataArray(
       np.random.randn(nlat, nlon) * 2 + 300,
       dims=['lat', 'lon'],
       coords={'lat': temp_3d.lat, 'lon': temp_3d.lon},
       attrs={'units': 'K', 'long_name': 'Sea surface temperature'}
   )

   psl_3d = xr.DataArray(
       np.random.randn(nlat, nlon) * 500 + 101325,
       dims=['lat', 'lon'],
       coords={'lat': temp_3d.lat, 'lon': temp_3d.lon},
       attrs={'units': 'Pa', 'long_name': 'Sea level pressure'}
   )

   # Calculate potential intensity for entire grid
   result_3d = potential_intensity(sst_3d, psl_3d, pres_xr, temp_3d, mixr_3d)
   print(f"3D PI shape: {result_3d.pi.shape}")
   print(f"Max PI: {result_3d.pi.max().values:.1f} m/s")
   print(f"Min central pressure: {result_3d.min_pressure.min().values:.1f} mb")

   # Analyze tropical regions only
   tropical_mask = np.abs(result_3d.lat) <= 25
   tropical_pi = result_3d.pi.where(tropical_mask)
   print(f"Mean tropical PI: {tropical_pi.mean().values:.1f} m/s")

   # === 4D Time Series Analysis ===
   ntimes = 12  # Monthly data

   # Create 4D data (time, level, lat, lon)
   temp_4d = xr.DataArray(
       np.random.randn(ntimes, nlevs, nlat, nlon) * 5 + 280,
       dims=['time', 'level', 'lat', 'lon'],
       coords={
           'time': pd.date_range('2023-01', periods=ntimes, freq='MS'),
           'level': pressure_levels,
           'lat': np.linspace(-30, 30, nlat),
           'lon': np.linspace(0, 360, nlon, endpoint=False)
       },
       attrs={'units': 'K'}
   )

   mixr_4d = xr.DataArray(
       np.random.rand(ntimes, nlevs, nlat, nlon) * 0.01 + 0.005,
       dims=['time', 'level', 'lat', 'lon'],
       coords=temp_4d.coords,
       attrs={'units': 'kg/kg'}
   )

   sst_4d = xr.DataArray(
       np.random.randn(ntimes, nlat, nlon) * 2 + 300,
       dims=['time', 'lat', 'lon'],
       coords={'time': temp_4d.time, 'lat': temp_4d.lat, 'lon': temp_4d.lon},
       attrs={'units': 'K'}
   )

   psl_4d = xr.DataArray(
       np.random.randn(ntimes, nlat, nlon) * 500 + 101325,
       dims=['time', 'lat', 'lon'],
       coords={'time': temp_4d.time, 'lat': temp_4d.lat, 'lon': temp_4d.lon},
       attrs={'units': 'Pa'}
   )

   # Calculate monthly potential intensity
   result_4d = potential_intensity(sst_4d, psl_4d, pres_xr, temp_4d, mixr_4d)
   print(f"4D result shape: {result_4d.pi.shape}")  # (time, lat, lon)

   # Seasonal analysis
   seasonal_pi = result_4d.pi.groupby('time.season').mean()
   print(f"Summer mean PI: {seasonal_pi.sel(season='JJA').mean().values:.1f} m/s")
   print(f"Winter mean PI: {seasonal_pi.sel(season='DJF').mean().values:.1f} m/s")

   # === Class Interface for Batch Processing ===
   from skyborn.calc.GPI.interface import PotentialIntensityCalculator

   # Process multiple profiles efficiently
   calculator = PotentialIntensityCalculator(
       sst_3d.values, psl_3d.values,
       pressure_levels, temp_3d.values, mixr_3d.values
   )
   calculator.calculate()
   results = calculator.results

   print(f"\nClass interface results:")
   print(f"  Mean PI: {results['pi'].mean():.1f} m/s")
   print(f"  Std PI: {results['pi'].std():.1f} m/s")
   print(f"  Success rate: {(results['error_flag'] == 1).mean() * 100:.1f}%")


**Statistical Analysis Examples**

.. code-block:: python

   # Statistical analysis
   x = np.random.randn(100)
   y = 2 * x + np.random.randn(100) * 0.5

   slope, p_value = skb.linear_regression(x, y)
   print(f"Linear regression: slope={slope:.4f}, p_value={p_value:.6f}")

   correlation = skb.pearson_correlation(x, y)
   print(f"Pearson correlation: {correlation:.4f}")

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
