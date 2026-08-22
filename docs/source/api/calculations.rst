Calculation
=====================

The Skyborn calculation module provides statistical, atmospheric, and mathematical functions for climate data analysis.

Atmospheric Physics Functions
-----------------------------

Tropopause Calculations
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.calc.troposphere.core.trop_wmo

.. autofunction:: skyborn.calc.troposphere.core.trop_wmo_profile

.. autofunction:: skyborn.calc.troposphere.xarray.trop_wmo

Growth-rate Diagnostics
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.calc.baroc_growth_rate

.. autofunction:: skyborn.calc.barot_growth_rate

.. note::

   ``skyborn.calc.baroc_growth_rate`` uses ``method="log"`` by
   default. Here ``"log"`` means interpolation that is linear in
   log-pressure, not a logarithmic transform of the field values.

   The public growth-rate wrappers return ``s^-1`` directly. Missing or
   unusable profiles return ``NaN`` instead of a user-supplied
   missing-value marker.

   By default, the compiled solver uses the fixed high-resolution
   Chemke/MATLAB-style zonal-wavenumber grid

   .. math::

      k = 0, 1\times10^{-7}, \dots, 1.99\times10^{-5}\ \mathrm{m}^{-1}.

   If you want the lower-resolution Chemke Python-share behavior instead,
   pass ``wavenumber_mode="low"``. When ``lon`` is omitted, Skyborn uses
   the original Chemke Python-share longitude grid
   ``np.arange(0.0, 361.0, 1.5)`` by default; you may still pass a custom
   one-dimensional longitude coordinate ``lon`` if needed. In that mode the
   solver uses

   .. math::

      k = \frac{2\pi w}{\left|\lambda_0 - \lambda_1\right|\pi a \cos\phi / 180},

   with :math:`w = 0, 1, \dots, \mathrm{round}(N_{\lambda}/3)-1`.

   When a latitude band is more appropriate than a single representative
   latitude, pass the band endpoints and the compiled solver uses the same
   Chemke-style averages as the original reference scripts:

   .. math::

      \sin_{\mathrm{avg}} = \frac{\sin(\phi_1) + \sin(\phi_2)}{2},
      \qquad
      \cos_{\mathrm{avg}} = \frac{\cos(\phi_1) + \cos(\phi_2)}{2}

   .. math::

      f = 2\Omega \sin_{\mathrm{avg}},
      \qquad
      \beta = \frac{2\Omega\cos_{\mathrm{avg}}}{a}

   If the input fields still carry a latitude axis, ``lat_bounds`` can also
   drive the latitude-band reduction step directly. Xarray inputs such as
   ``(time, level, lat)`` are selected and cosine-weighted over that band
   automatically before the compiled profile solver runs. For raw NumPy
   arrays with an explicit latitude axis, provide the one-dimensional
   latitude coordinate through ``lat=...`` together with ``lat_bounds`` so the
   requested band can be selected unambiguously. The recommended NumPy layout
   is ``(time, level, lat)``. More generally, the wrapper expects exactly one
   axis matching the pressure-coordinate length and one non-pressure axis
   matching the supplied latitude-coordinate length.

   The function diagnoses the WMO tropopause pressure by default, or uses an
   explicit ``tropopause_pressure`` when given, and then builds a fixed solver
   grid with
   ``solver_levels=45`` by default between that tropopause and the lower
   troposphere. This is a practical resolution choice for the compiled
   eigenvalue problem, not a universal optimum: increasing ``solver_levels``
   can improve vertical-grid convergence up to a point, but it also raises
   the per-profile cost and is not automatically better once the solution is
   converged.

   In the original Chemke-style analysis scripts, a parameter named
   ``window_size`` is often used for a centered running mean applied to the
   growth-rate spectrum over zonal wavenumber. In Skyborn this role is
   represented by ``smooth_window`` in ``baroc_growth_rate``. That smoothing
   width is not the same thing as ``solver_levels``. Even for a single
   atmospheric column, the compiled solver still sweeps many zonal
   wavenumbers and solves a generalized eigenvalue problem for each one, so
   the dominant cost is the repeated linear-algebra solve rather than
   tropopause diagnosis or pressure interpolation. The public default is
   ``smooth_window=1``; values greater than 1 apply the centered running mean
   inside the Fortran backend before the final Chemke-style maximum-growth
   diagnostic is taken.

.. rubric:: References

The current growth-rate implementation is documented against the following
Chemke publications that are directly related to baroclinic wave growth and
storm-track changes:

* Chemke, R., and Ming, Y. (2020): *Large Atmospheric Waves Will Get
  Stronger, While Small Waves Will Get Weaker by the End of the 21st
  Century*. Geophysical Research Letters, 47, ``e2020GL090441``.
  https://doi.org/10.1029/2020GL090441
* Chemke, R., Zanna, L., Orbe, C., Sentman, L. T., and Polvani, L. M.
  (2022): *The Future Intensification of the North Atlantic Winter Storm
  Track: The Key Role of Dynamic Ocean Coupling*. Journal of Climate,
  35(8), 2407-2421. https://doi.org/10.1175/JCLI-D-21-0407.1
* Chemke, R. (2022): *The future poleward shift of Southern Hemisphere
  summer mid-latitude storm tracks stems from ocean coupling*. Nature
  Communications, 13, 2531. https://doi.org/10.1038/s41467-022-29392-4
* Chemke, R., Ming, Y., and Yuval, J. (2022): *The intensification of
  winter mid-latitude storm tracks in the Southern Hemisphere*. Nature
  Climate Change, 12, 553-557. https://doi.org/10.1038/s41558-022-01368-8

Geostrophic Wind Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.calc.geostrophic.core.geostrophic_wind

.. autoclass:: skyborn.calc.geostrophic.core.GeostrophicWind

.. autofunction:: skyborn.calc.geostrophic.xarray.geostrophic_wind

.. autoclass:: skyborn.calc.geostrophic.xarray.GeostrophicWind

Genesis Potential Index (GPI) / Tropical Cyclone Potential Intensity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.calc.GPI.core.potential_intensity

.. autofunction:: skyborn.calc.GPI.core.log_decompose_pi

.. autofunction:: skyborn.calc.GPI.core.pi_log_decomposition

.. autofunction:: skyborn.calc.GPI.xarray.potential_intensity

.. autofunction:: skyborn.calc.GPI.xarray.pi_log_decomposition

.. note::

   ``error_flag == 1`` indicates a successful PI solve.

   ``potential_intensity(...)`` remains the pure PI entry point.

   Use ``pi_log_decomposition(...)`` for supported 1D, 3D, and 4D inputs
   when you need ``tcpyPI``-style outflow diagnostics
   (``outflow_source={"cape_star", "cape_env"}``, ``t0``, ``otl``) and the
   Wing et al. (2015) logarithmic decomposition where ``lnpi = ln(V^2)``.

Convective Diagnostics
~~~~~~~~~~~~~~~~~~~~~~~

CAPE and CIN Calculations
**************************

.. autofunction:: skyborn.calc.cape.calculate_cape_cin

.. autofunction:: skyborn.calc.cape.calculate_most_unstable_parcel

.. autofunction:: skyborn.calc.cape.calculate_most_unstable_cape_cin

.. autofunction:: skyborn.calc.cape.calculate_parcel_profile

.. autofunction:: skyborn.calc.cape.cape_grid

.. autofunction:: skyborn.calc.cape.most_unstable_parcel_grid

.. autofunction:: skyborn.calc.cape.most_unstable_cape_cin_grid

.. autofunction:: skyborn.calc.cape.parcel_profile_grid

.. note::

   The CAPE/CIN module provides both 1D profile and 3D gridded interfaces:

   - **1D functions** (``calculate_*``): Accept 1D NumPy arrays and return scalar values.
     Optimized for single atmospheric profile analysis.

   - **3D functions** (``*_grid``): Accept 3D NumPy arrays ``(nz, ny, nx)`` and return 2D fields ``(ny, nx)``.
     Optimized for gridded data analysis with vectorized operations.

   - **xarray wrappers**: Automatically detect input dimensions and apply appropriate backend.
     Handle coordinate preservation and metadata.

   Performance: ~10,000 profiles/second, 195× faster than MetPy's per-profile loops.

   Validation: CAPE within +1.1% of MetPy 1.7.1 (difference due to Bolton vs Romps LCL
   and RK4 vs LSODA integration methods).

Storm Relative Helicity
************************

.. autofunction:: skyborn.calc.srh.calculate_storm_relative_helicity

.. autofunction:: skyborn.calc.srh.srh_grid

.. note::

   The SRH module calculates storm-relative helicity for tornado and severe storm forecasting:

   - Supports custom layer depths (default 0-3 km AGL)
   - Returns positive, negative, and total helicity components
   - Handles AGL height conversion and layer interpolation
   - Matches MetPy's ``interpolate_1d`` NaN semantics

   Performance: ~2,000,000 profiles/second, 5258× faster than MetPy.

   Validation: Bit-for-bit identical to MetPy 1.7.1 reference soundings.

Ventilation Index for Tropical Cyclones
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: skyborn.calc.ventilation.vertical_wind_shear

.. autofunction:: skyborn.calc.ventilation.entropy_deficit

.. autofunction:: skyborn.calc.ventilation.air_sea_disequilibrium

.. autofunction:: skyborn.calc.ventilation.ventilation_index

.. autofunction:: skyborn.calc.ventilation.ventilated_pi

.. autofunction:: skyborn.calc.ventilation.absolute_vorticity_850

.. autofunction:: skyborn.calc.ventilation.genesis_potential_index

.. autofunction:: skyborn.calc.ventilation.ventilation_components

.. note::

   The ventilation index module implements the ventilated potential intensity (vPI)
   framework from Chavas, Camargo & Tippett (2025, J. Climate):

   - **Ventilation Index**: VI = VWS × Chi / PI following Tang & Emanuel (2012)
   - **Ventilated PI**: Analytic cubic solution using Cardano formula
   - **Genesis Potential Index**: GPIv incorporating ventilation effects

   All functions accept and return xarray.DataArray for coordinate-aware workflows.

   Performance: 16× speedup through NumPy vectorization with complex number support.

   Use ``ventilation_components()`` as a high-level pipeline function that calculates
   all ventilation-related diagnostics from a single ERA5-style dataset.

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
   from skyborn.calc.GPI.core import potential_intensity
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

   # === PI log decomposition diagnostics ===
   sst_profile = xr.DataArray(sst, attrs={'units': 'K'})
   psl_profile = xr.DataArray(psl, attrs={'units': 'Pa'})
   temp_profile_xr = xr.DataArray(temperature, dims=['level'], coords={'level': pressure_levels}, attrs={'units': 'K'})
   mixr_profile_xr = xr.DataArray(mixing_ratio, dims=['level'], coords={'level': pressure_levels}, attrs={'units': 'kg/kg'})
   pres_profile_xr = xr.DataArray(pressure_levels, dims=['level'], attrs={'units': 'mb'})

   from skyborn.calc.GPI.xarray import pi_log_decomposition

   profile_full = pi_log_decomposition(
       sst_profile,
       psl_profile,
       pres_profile_xr,
       temp_profile_xr,
       mixr_profile_xr,
       outflow_source='cape_env'
   )
   print(profile_full[['pi', 't0', 'otl', 'lnpi', 'lneff', 'lndiseq', 'lnCKCD']])


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

**CAPE/CIN Calculation Examples**

.. code-block:: python

   import skyborn as skb
   import numpy as np
   import xarray as xr

   # === 1D Profile Calculation ===
   # Single atmospheric sounding
   height = np.array([0, 500, 1000, 1500, 2000, 3000, 5000, 8000, 10000])  # meters AGL
   temperature = np.array([25, 20, 15, 10, 5, -5, -20, -40, -50]) + 273.15  # K
   dewpoint = np.array([20, 15, 10, 5, 0, -10, -30, -50, -60]) + 273.15  # K

   # Calculate surface-based CAPE and CIN
   from skyborn.calc.cape import calculate_cape_cin
   cape, cin = calculate_cape_cin(height, temperature, dewpoint)
   print(f"CAPE: {cape:.1f} J/kg, CIN: {cin:.1f} J/kg")

   # Calculate most unstable parcel
   from skyborn.calc.cape import calculate_most_unstable_parcel
   mu_height, mu_temp, mu_dewpoint = calculate_most_unstable_parcel(
       height, temperature, dewpoint
   )
   print(f"Most unstable parcel at {mu_height:.0f} m AGL")

   # Calculate most unstable CAPE/CIN
   from skyborn.calc.cape import calculate_most_unstable_cape_cin
   mu_cape, mu_cin = calculate_most_unstable_cape_cin(height, temperature, dewpoint)
   print(f"MU-CAPE: {mu_cape:.1f} J/kg, MU-CIN: {mu_cin:.1f} J/kg")

   # === 3D Gridded Data ===
   # ERA5-style data (level, lat, lon)
   nz, ny, nx = 25, 50, 100
   height_3d = np.random.randn(nz, ny, nx) * 500 + np.arange(nz)[:, None, None] * 500
   temp_3d = np.random.randn(nz, ny, nx) * 3 + 280 - np.arange(nz)[:, None, None] * 6
   dewpoint_3d = temp_3d - np.random.rand(nz, ny, nx) * 10 - 5

   from skyborn.calc.cape import cape_grid, most_unstable_cape_cin_grid

   # Calculate CAPE/CIN for entire grid
   cape_field, cin_field = cape_grid(height_3d, temp_3d, dewpoint_3d)
   print(f"CAPE field shape: {cape_field.shape}")  # (ny, nx)
   print(f"Maximum CAPE: {cape_field.max():.1f} J/kg")

   # Most unstable CAPE/CIN
   mu_cape_field, mu_cin_field = most_unstable_cape_cin_grid(
       height_3d, temp_3d, dewpoint_3d
   )

   # === xarray Interface ===
   # Automatic dimension detection and coordinate preservation
   ds = xr.open_dataset('era5_data.nc')

   # Compute CAPE/CIN with coordinate awareness
   cape_xr = skb.calc.cape.cape_grid(
       ds.geopotential_height,
       ds.temperature,
       ds.dewpoint
   )
   # Result preserves lat/lon coordinates

**Storm Relative Helicity Examples**

.. code-block:: python

   import skyborn as skb
   import numpy as np

   # === 1D Profile Calculation ===
   # Hodograph data
   height = np.array([0, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000])  # m AGL
   u_wind = np.array([5, 8, 12, 15, 18, 20, 22, 25, 28])  # m/s
   v_wind = np.array([0, 2, 5, 8, 10, 12, 13, 15, 16])  # m/s

   # Calculate 0-3 km storm relative helicity
   from skyborn.calc.srh import calculate_storm_relative_helicity

   # Bunkers right-moving storm motion
   storm_u, storm_v = 15.0, 8.0

   srh_positive, srh_negative, srh_total = calculate_storm_relative_helicity(
       height, u_wind, v_wind,
       storm_u=storm_u, storm_v=storm_v,
       depth=3000.0,  # 0-3 km layer
       bottom=0.0
   )

   print(f"0-3 km SRH: {srh_total:.1f} m²/s²")
   print(f"Positive: {srh_positive:.1f} m²/s², Negative: {srh_negative:.1f} m²/s²")

   # Calculate 0-1 km SRH (for tornado potential)
   srh_pos_01, srh_neg_01, srh_tot_01 = calculate_storm_relative_helicity(
       height, u_wind, v_wind,
       storm_u=storm_u, storm_v=storm_v,
       depth=1000.0,
       bottom=0.0
   )
   print(f"0-1 km SRH: {srh_tot_01:.1f} m²/s²")

   # === 3D Gridded Data ===
   from skyborn.calc.srh import srh_grid

   nz, ny, nx = 30, 50, 100
   height_3d = np.arange(nz)[:, None, None] * 300.0  # 300m intervals
   u_3d = np.random.randn(nz, ny, nx) * 5 + 10
   v_3d = np.random.randn(nz, ny, nx) * 3 + 5

   # Storm motion fields (2D)
   storm_u_2d = np.ones((ny, nx)) * 15.0
   storm_v_2d = np.ones((ny, nx)) * 8.0

   # Calculate SRH for entire grid
   srh_pos_grid, srh_neg_grid, srh_tot_grid = srh_grid(
       height_3d, u_3d, v_3d,
       storm_u=storm_u_2d, storm_v=storm_v_2d,
       depth=3000.0, bottom=0.0
   )

   print(f"SRH grid shape: {srh_tot_grid.shape}")  # (ny, nx)

**Ventilation Index Examples**

.. code-block:: python

   import skyborn as skb
   import xarray as xr

   # === Load ERA5 Data ===
   ds = xr.open_dataset('era5_tropical.nc')
   # Expected variables: U, V, T, Q, SSTK, SP
   # Expected levels: 200, 600, 850 hPa (and others)

   # === High-Level Pipeline (Recommended) ===
   from skyborn.calc.ventilation import ventilation_components

   # Calculate all ventilation diagnostics in one call
   result = ventilation_components(ds)

   # Result contains:
   print(result.data_vars)
   # Dict: PI, vPI, VWS, Chi, ventilation_index, eta_c, GPIv, air_sea_disequilibrium

   # Visualize ventilated genesis potential index
   result.GPIv.plot()

   # === Step-by-Step Calculation ===
   from skyborn.calc.ventilation import (
       vertical_wind_shear,
       entropy_deficit,
       air_sea_disequilibrium,
       ventilation_index,
       ventilated_pi,
       absolute_vorticity_850,
       genesis_potential_index
   )

   # Step 1: Calculate vertical wind shear (200-850 hPa)
   vws = vertical_wind_shear(ds)
   print(f"VWS: {vws.mean().values:.2f} m/s")

   # Step 2: Calculate air-sea disequilibrium
   chi_star = air_sea_disequilibrium(ds)

   # Step 3: Calculate entropy deficit (600 hPa)
   chi = entropy_deficit(ds, chi_star)
   print(f"Entropy deficit: {chi.mean().values:.3f}")

   # Step 4: Calculate potential intensity (from GPI module)
   from skyborn.calc.GPI.xarray import potential_intensity
   pi_result = potential_intensity(
       ds.SSTK, ds.SP,
       ds.level, ds.T, ds.Q
   )
   pi = pi_result.pi

   # Step 5: Calculate ventilation index
   vi = ventilation_index(vws, chi, pi)
   print(f"VI: {vi.mean().values:.4f}")

   # Step 6: Calculate ventilated potential intensity
   vpi = ventilated_pi(pi, vi)
   print(f"Mean PI: {pi.mean().values:.1f} m/s")
   print(f"Mean vPI: {vpi.mean().values:.1f} m/s")
   print(f"Mean reduction: {(1 - vpi/pi).mean().values * 100:.1f}%")

   # Step 7: Calculate absolute vorticity at 850 hPa
   eta_c = absolute_vorticity_850(ds)

   # Step 8: Calculate ventilated genesis potential index
   gpiv = genesis_potential_index(vpi, eta_c)
   print(f"Mean GPIv: {gpiv.mean().values:.2e}")

   # === Sensitivity Analysis ===
   # Test different VI_MAX thresholds
   result_sensitive = ventilation_components(
       ds,
       vi_max=0.20,  # Default is 0.145
       vorticity_cap=5.0e-5  # Default is 3.7e-5
   )

Emergent Constraints
====================

The emergent constraint module implements methods for reducing uncertainty in climate projections
by leveraging relationships between observable present-day quantities and uncertain future projections.

For a complete interactive example, see :doc:`../notebooks/ecs_emergent_constraints_analysis`.
