Calculation
=====================

The Skyborn calculation module provides statistical, atmospheric, and mathematical functions for climate data analysis.

Atmospheric Physics Functions
-----------------------------

.. automodule:: skyborn.calc.tropopause
   :members:

.. automodule:: skyborn.calc.tropopause_xarray
   :members:

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
   pressure = np.array([1000, 850, 700, 500, 300, 200, 100, 50])  # hPa
   temperature = np.array([288, 280, 270, 250, 230, 220, 210, 200])  # K

   result = skb.calc.trop_wmo_profile(temperature, pressure)
   print(f"Tropopause: {result['pressure']:.1f} hPa at {result['height']:.0f} m")

   # === xarray Interface - Simplified ===
   # Load atmospheric data with level coordinate
   ds = xr.open_dataset('era5_data.nc')  # Temperature with level coordinate in hPa

   # Auto-generate pressure from level coordinate - no pressure array needed!
   result = skb.calc.tropopause_xarray.trop_wmo(ds.temperature)
   print(f"Global tropopause calculated for {result.pressure.shape}")

   # === 2D Cross-section Analysis ===
   # Meridional cross-section (level, lat)
   temp_meridional = ds.temperature.isel(time=0, lon=0)  # (level, lat)
   result_2d = skb.calc.tropopause_xarray.trop_wmo(temp_meridional)
   # Result shape: (lat,) - tropopause height at each latitude

   # === 4D Climate Analysis ===
   # Multi-year dataset (time, level, lat, lon)
   result_4d = skb.calc.tropopause_xarray.trop_wmo(ds.temperature)
   # Result shape: (time, lat, lon) - preserves time and spatial dimensions

   # Seasonal analysis
   seasonal_mean = result_4d.height.groupby('time.season').mean()

   # === Advanced Usage ===
   # Custom pressure field and WMO criterion
   result = skb.calc.tropopause_xarray.trop_wmo(
       temperature=ds.temperature,
       pressure=ds.pressure,  # Custom pressure field
       lapse_criterion=2.5,   # Custom WMO threshold (K/km)
       pressure_unit='Pa'     # If pressure is in Pascals
   )

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
