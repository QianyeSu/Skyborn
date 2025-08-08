Mann-Kendall Trend Analysis
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
    results = mann_kendall_multidim(climate_data, time_axis=0)

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

Performance Optimization
------------------------

For Large Datasets
~~~~~~~~~~~~~~~~~~

The implementation automatically optimizes performance for large climate datasets:

.. code-block:: python

    # For very large datasets, control memory usage
    results = mann_kendall_multidim(
        large_climate_data,
        time_axis=0,
        chunk_size=2000  # Process 2000 grid points at once
    )

Expected Performance
~~~~~~~~~~~~~~~~~~~

Typical processing speeds for climate data:

* **Small grids** (50×20×30): ~1,500 grid points/second
* **Climate grids** (40×192×288): ~1,800 grid points/second (~30 seconds total)
* **Large grids** (100×360×720): ~600 grid points/second

Memory usage is minimal (~25MB) regardless of grid size due to intelligent chunking.

API Reference
-------------

Single Time Series Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mann_kendall_test

Multidimensional Analysis
~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: mann_kendall_multidim

XArray Interface
~~~~~~~~~~~~~~~

.. autofunction:: mann_kendall_xarray

Unified Interface
~~~~~~~~~~~~~~~~

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

Modified Mann-Kendall Test
~~~~~~~~~~~~~~~~~~~~~~~~~~

For autocorrelated data, the modified Mann-Kendall test (Yue & Wang, 2004) can be used:

.. code-block:: python

    result = mann_kendall_test(data, modified=True)

This accounts for serial correlation by adjusting the variance calculation.

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Examples Gallery
----------------

The following examples demonstrate various applications:

.. toctree::
   :maxdepth: 1

   examples/climate_warming_trends
   examples/precipitation_variability
   examples/sea_level_analysis
   examples/model_comparison

Best Practices
--------------

Data Preparation
~~~~~~~~~~~~~~~

1. **Quality control**: Remove obviously erroneous values
2. **Homogenization**: Ensure data consistency over time
3. **Gap analysis**: Understand missing data patterns
4. **Deseasonalization**: Remove seasonal cycles if needed

Interpretation
~~~~~~~~~~~~~

1. **Physical significance**: Consider whether trends are physically meaningful
2. **Spatial coherence**: Look for consistent patterns across neighboring regions
3. **Multiple variables**: Cross-validate trends across different measurements
4. **Uncertainty**: Report confidence intervals and significance levels

References
----------

* Mann, H. B. (1945). Nonparametric tests against trend. *Econometrica*, 13(3), 245-259.
* Kendall, M. G. (1948). *Rank correlation methods*. Griffin, London.
* Theil, H. (1950). A rank-invariant method of linear and polynomial regression analysis. *Indagationes Mathematicae*, 12(85), 173.
* Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association*, 63(324), 1379-1389.
* Yue, S., & Wang, C. (2004). The Mann-Kendall test modified by effective sample size to detect trend in serially correlated hydrological series. *Water Resources Management*, 18(3), 201-218.
