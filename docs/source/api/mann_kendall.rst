Mann-Kendall
============================

.. currentmodule:: skyborn.calc

The Mann-Kendall test is a non-parametric statistical test used to detect trends in time series data. This module provides optimized implementations for both single time series and multidimensional arrays, making it ideal for climate data analysis.

Key Features
------------

* **High Performance**: Vectorized implementation with 15-30x speed improvements
* **Climate Data Optimized**: Efficient processing of typical climate grids (40×192×288)
* **Memory Efficient**: Smart chunking with minimal memory overhead (~25MB for full climate datasets)
* **Robust Statistics**: Handles missing data and provides comprehensive trend statistics
* **Multiple Interfaces**: Support for NumPy arrays, xarray DataArrays, and Dask arrays

Supported Test Families
-----------------------

The public API uses ``test=...`` to select the Mann-Kendall test family:

* ``test="original"``: original Mann-Kendall test
* ``test="yue_wang"``: Yue-Wang (2004) modified variance correction
* ``test="seasonal"``: Hirsch-Slack (1984) seasonal Mann-Kendall test
* ``test="correlated_seasonal"``: Hipel (1994) correlated seasonal Mann-Kendall test
* ``test="correlated_multivariate"``: Libiseller-Grimvall correlated multivariate Mann-Kendall test
* ``test="multivariate"``: grouped multivariate Mann-Kendall test
* ``test="regional"``: grouped regional Mann-Kendall test
* ``test="hamed_rao"``: Hamed-Rao (1998) variance correction
* ``test="pre_whitening"``: Yue-Wang (2002) pre-whitening modification
* ``test="trend_free_pre_whitening"``: Yue-Wang (2002) trend-free pre-whitening

All interfaces default to ``test="original"``.

The optional ``lag=...`` argument is only used by
``test="yue_wang"`` and ``test="hamed_rao"``. Other test families ignore it.

Partial Mann-Kendall
--------------------

Partial Mann-Kendall is exposed through a separate API because it needs two
inputs: a response series and a covariate series. It is therefore not selected
through ``test=...`` in the single-series MK interfaces.

Available functions:

* ``partial_mann_kendall_test``: one response series plus one covariate series
* ``partial_mann_kendall_multidim``: multidimensional NumPy arrays
* ``partial_mann_kendall_xarray``: xarray ``DataArray`` inputs

Quick Start
-----------

Basic Trend Detection
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import numpy as np
    from skyborn.calc import mann_kendall_test

    # Create sample time series with trend
    time = np.arange(50)
    data = 0.02 * time + np.random.randn(50) * 0.5

    # Perform Mann-Kendall test
    result = mann_kendall_test(data)

    print(f"Trend: {result['trend']:.4f} units/year")
    print(f"Significant: {result['h']} (p={result['p']:.3f})")
    print(f"Mann-Kendall tau: {result['tau']:.3f}")

Multidimensional Climate Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import numpy as np
    from skyborn.calc import mann_kendall_multidim

    # Climate data: (time, lat, lon)
    climate_data = np.random.randn(40, 192, 288)

    # Add realistic warming trend
    years = np.arange(40)
    warming_pattern = np.random.randn(192, 288) * 0.02
    for t, year in enumerate(years):
        climate_data[t] += warming_pattern * year

    # Analyze trends across all grid points
    results = mann_kendall_multidim(climate_data, axis=0)

    print(f"Grid shape: {results['trend'].shape}")  # (192, 288)
    print(f"Significant trends: {np.sum(results['h'])} grid points")
    print(f"Mean warming: {np.nanmean(results['trend']):.4f} ± {np.nanstd(results['trend']):.4f}")

Using xarray Interface
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import xarray as xr
    import pandas as pd
    from skyborn.calc import mann_kendall_xarray

    # Create xarray DataArray
    time = pd.date_range('1980', periods=40, freq='AS')
    data = xr.DataArray(
        np.random.randn(40, 96, 144),
        coords={'time': time, 'lat': np.arange(96), 'lon': np.arange(144)},
        dims=['time', 'lat', 'lon']
    )

    # Perform trend analysis
    trends = mann_kendall_xarray(data, dim='time')

    # Results are returned as xarray Dataset
    print(trends.trend.attrs)  # Includes metadata
    trends.trend.plot()  # Easy visualization

Grouped Multivariate / Regional MK
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Grouped MK families collapse both the analyzed time dimension and one grouping
dimension, while preserving the remaining spatial dimensions.

.. code-block:: python

    grouped = xr.DataArray(
        data,
        dims=["time", "member", "lat", "lon"],
    )

    regional = mann_kendall_xarray(
        grouped,
        dim="time",
        group_dim="member",
        test="regional",
    )

    print(regional.trend.dims)  # ('lat', 'lon')

Partial MK With a Covariate
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Partial MK is useful when the trend in a response variable should be tested
after accounting for a second covariate, for example precipitation after
controlling for an ENSO index or water quality after controlling for streamflow.

.. code-block:: python

    import numpy as np
    from skyborn.calc import partial_mann_kendall_multidim

    # Response field: (time, lat, lon)
    response = np.random.randn(40, 96, 144)

    # One global covariate shared by all grid points
    covariate = np.linspace(-1.0, 1.0, 40)

    partial = partial_mann_kendall_multidim(
        response,
        covariate,
        axis=0,
    )

    print(partial["trend"].shape)  # (96, 144)

Performance Optimization
------------------------

For Large Datasets
~~~~~~~~~~~~~~~~~~

The implementation automatically optimizes performance for large climate datasets:

.. code-block:: python

    # For very large datasets, control memory usage
    results = mann_kendall_multidim(
        large_climate_data,
        axis=0,
        chunk_size=2000  # Process 2000 grid points at once
    )

Expected Performance
~~~~~~~~~~~~~~~~~~~~

Typical processing speeds for climate data:

* **Small grids** (50×20×30): ~1,500 grid points/second
* **Climate grids** (40×192×288): ~1,800 grid points/second (~30 seconds total)
* **Large grids** (100×360×720): ~600 grid points/second

Memory usage is minimal (~25MB) regardless of grid size due to intelligent chunking.

API Reference
-------------

Single Time Series Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mann_kendall_test

Multidimensional Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mann_kendall_multidim

XArray Interface
~~~~~~~~~~~~~~~~

.. autofunction:: mann_kendall_xarray

Partial Single-Series Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: partial_mann_kendall_test

Partial Multidimensional Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: partial_mann_kendall_multidim

Partial XArray Interface
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: partial_mann_kendall_xarray

Unified Interface
~~~~~~~~~~~~~~~~~

.. autofunction:: trend_analysis

Statistical Background
----------------------

The Mann-Kendall Test
~~~~~~~~~~~~~~~~~~~~~

The Mann-Kendall test is based on the Mann-Kendall statistic S:

.. math::

    S = \sum_{i=1}^{n-1} \sum_{j=i+1}^{n} \text{sign}(x_j - x_i)

where :math:`\text{sign}(x)` is the sign function. Under the null hypothesis of no trend,
S follows approximately a normal distribution for large n.

**Test Statistics:**

* **S**: Mann-Kendall statistic
* **tau**: Kendall's tau coefficient (normalized correlation)
* **z**: Standardized test statistic
* **p**: Two-tailed p-value
* **h**: Boolean significance test result

**Slope Estimation:**

The trend magnitude is estimated using either:

* **Theil-Sen estimator** (default): Robust, non-parametric slope estimation
* **Linear regression**: Ordinary least squares for comparison

Serial-Correlation-Aware Tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For autocorrelated data, explicit test families can be selected:

.. code-block:: python

    result_yw = mann_kendall_test(data, test="yue_wang")
    result_seasonal = mann_kendall_test(data, test="seasonal", period=12)
    result_corr_seasonal = mann_kendall_test(data, test="correlated_seasonal", period=12)
    result_corr_multi = mann_kendall_test(grouped_matrix, test="correlated_multivariate")
    result_multi = mann_kendall_test(grouped_matrix, test="multivariate")
    result_regional = mann_kendall_xarray(
        grouped_da, dim="time", group_dim="member", test="regional"
    )
    result_hr = mann_kendall_test(data, test="hamed_rao")
    result_pw = mann_kendall_test(data, test="pre_whitening")
    result_tfpw = mann_kendall_test(data, test="trend_free_pre_whitening")

These variants account for serial correlation using different correction strategies.

Advantages
----------

* **Non-parametric**: No assumptions about data distribution
* **Robust**: Resistant to outliers and non-linear trends
* **Missing data**: Handles gaps in time series naturally
* **Significance testing**: Built-in statistical inference
* **Climate appropriate**: Designed for geophysical time series

Use Cases
---------

Climate Science Applications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Temperature trends**: Global and regional warming patterns
* **Precipitation changes**: Long-term rainfall trend detection
* **Sea level rise**: Coastal monitoring and analysis
* **Extreme events**: Frequency and intensity trend analysis
* **Model validation**: Comparing observed vs. simulated trends

Quality Control
~~~~~~~~~~~~~~~

* **Data homogeneity**: Detecting artificial trends in observations
* **Station records**: Long-term consistency checking
* **Reanalysis validation**: Trend comparison across datasets

Best Practices
--------------

Data Preparation
~~~~~~~~~~~~~~~~

1. **Quality control**: Remove obviously erroneous values
2. **Homogenization**: Ensure data consistency over time
3. **Gap analysis**: Understand missing data patterns
4. **Deseasonalization**: Remove seasonal cycles if needed

Interpretation
~~~~~~~~~~~~~~

1. **Physical significance**: Consider whether trends are physically meaningful
2. **Spatial coherence**: Look for consistent patterns across neighboring regions
3. **Multiple variables**: Cross-validate trends across different measurements
4. **Uncertainty**: Report confidence intervals and significance levels

References
----------

* Bari, S. H., Rahman, M. T. U., Hoque, M. A., & Hussain, M. M. (2016). Analysis of seasonal and annual rainfall trends in the northern region of Bangladesh. *Atmospheric Research*, 176, 148-158. DOI: `<https://doi.org/10.1016/j.atmosres.2016.02.008>`_
* Conover, W. J. (1980). Some methods based on ranks (Chapter 5). In *Practical nonparametric statistics* (2nd ed.). John Wiley and Sons.
* Cox, D. R., & Stuart, A. (1955). Some quick sign tests for trend in location and dispersion. *Biometrika*, 42(1/2), 80-95. DOI: `<https://doi.org/10.2307/2333424>`_
* Dietz, E. J. (1987). A comparison of robust estimators in simple linear regression. *Communications in Statistics - Simulation and Computation*, 16(4), 1209-1227. DOI: `<https://doi.org/10.1080/03610918708812645>`_
* Hamed, K. H., & Rao, A. R. (1998). A modified Mann-Kendall trend test for autocorrelated data. *Journal of Hydrology*, 204(1-4), 182-196. DOI: `<https://doi.org/10.1016/S0022-1694(97)00125-X>`_
* Helsel, D. R., & Frans, L. M. (2006). Regional Kendall test for trend. *Environmental Science & Technology*, 40(13), 4066-4073. DOI: `<https://doi.org/10.1021/es051650b>`_
* Hipel, K. W., & McLeod, A. I. (1994). *Time series modelling of water resources and environmental systems* (Vol. 45). Elsevier.
* Hirsch, R. M., Slack, J. R., & Smith, R. A. (1982). Techniques of trend analysis for monthly water quality data. *Water Resources Research*, 18(1), 107-121. DOI: `<https://doi.org/10.1029/WR018i001p00107>`_
* Kendall, M. (1975). *Rank correlation measures*. Charles Griffin, London.
* Libiseller, C., & Grimvall, A. (2002). Performance of partial Mann-Kendall tests for trend detection in the presence of covariates. *Environmetrics*, 13(1), 71-84. DOI: `<https://doi.org/10.1002/env.507>`_
* Mann, H. B. (1945). Nonparametric tests against trend. *Econometrica*, 13(3), 245-259. DOI: `<https://doi.org/10.2307/1907187>`_
* Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association*, 63(324), 1379-1389. DOI: `<https://doi.org/10.1080/01621459.1968.10480934>`_
* Theil, H. (1950). A rank-invariant method of linear and polynomial regression analysis (Parts 1-3). In *Ned. Akad. Wetensch. Proc. Ser. A* (Vol. 53, pp. 1397-1412).
* Yue, S., & Wang, C. (2004). The Mann-Kendall test modified by effective sample size to detect trend in serially correlated hydrological series. *Water Resources Management*, 18(3), 201-218. DOI: `<https://doi.org/10.1023/B:WARM.0000043140.61082.60>`_
* Yue, S., & Wang, C. Y. (2002). Applicability of prewhitening to eliminate the influence of serial correlation on the Mann-Kendall test. *Water Resources Research*, 38(6), 4-1. DOI: `<https://doi.org/10.1029/2001WR000861>`_
* Yue, S., Pilon, P., Phinney, B., & Cavadias, G. (2002). The influence of autocorrelation on the ability to detect trend in hydrological series. *Hydrological Processes*, 16(9), 1807-1829. DOI: `<https://doi.org/10.1002/hyp.1095>`_
